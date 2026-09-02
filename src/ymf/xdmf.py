"""Convert between a YMF data tree and an XDMF 2.0 XML document, round-trip.

Per ``docs/ymf-schema.md`` §5 ("XDMF Preservation") and the design
principle that *XDMF is a consumer, not a dependency* of YMF: this module
reads a plain Python dict describing a time-collection of grids (mesh
topology + geometry, optionally HDF5-backed) and emits an ``.xmf`` file
that ParaView and other XDMF-aware tools can open directly, using only the
parts of the YMF document that XDMF actually understands.

Round-trip strategy (write -> read -> get the same YMF document back)
-----------------------------------------------------------------------
XDMF's grammar (the ``Xdmf.dtd``/``Xdmf.xsd``, referenced in
``examples/notebooks/YMF-original.ipynb`` cell 1) defines an
``<Information Name="..." Value="...">`` element that is legal as a child
of essentially any XDMF element (``Domain``, ``Grid``, ``DataItem``, ...)
specifically as an extension point: generic XDMF consumers (ParaView,
VisIt) parse and then ignore ``Information`` elements they don't
recognize, but the element itself survives round-trips through any
XDMF-compliant XML tool, unlike XML comments (which many parsers,
including a plain ``xml.etree.ElementTree`` reader with default settings,
silently drop) or non-standard custom elements (which a strict
schema-validating XDMF reader could legitimately reject).

So: everything in a YMF document that *does* map onto XDMF concepts
(``Domain``/``Grid``/``Topology``/``Geometry``/``DataItem`` — i.e. the mesh
archive) is written as real XDMF elements. Everything else (``Problem``,
``solution_paths``, ``vvuq``, and any ``archive``-level metadata XDMF has
no concept of) is serialized to JSON, base64-encoded (safe inside an XML
attribute value with no escaping headaches), and stashed in a single
``<Information Name="YMF" Value="...">`` child of ``<Domain>``.

JSON rather than YAML is used for this *specific* embedded payload:
strictyaml's ``as_document()`` cannot serialize empty lists/dicts without
an explicit schema (raises ``YAMLSerializationError`` on
``solution_paths.analytical: []``, which every current example uses), and
since this payload is base64-encoded anyway — not meant to be read
directly out of the XML — YAML's human-authoring niceties (comments, flow
style, block scalars) buy nothing here. The rest of the package still uses
YAML/strictyaml everywhere a human or LLM actually edits a YMF document;
only this internal round-trip encoding uses JSON.

**Verified**: checked directly against the live ``Xdmf.dtd`` fetched
from gitlab.kitware.com/xdmf/xdmf (2026-08-12, outside this sandbox's own
no-network environment). The grammar confirms
``<!ELEMENT Domain (Information*, Grid+)>`` — ``Information`` is
explicitly legal as a direct child of ``Domain``, ahead of ``Grid+`` in
its content model. No fallback needed (``Grid`` also permits
``Information*``, so the originally-documented ``<Grid>`` fallback would
have worked too, it just isn't required).

Lossiness
---------
The round-trip is **not** claimed to be fully lossless in general:

- Mesh/archive data (``Domain``) round-trips exactly (it *is* the native
  XDMF representation already).
- Everything under ``Problem``/``solution_paths``/``vvuq`` round-trips
  exactly *as long as it survives a YAML dump/load cycle* — which is true
  for everything currently in the schema (plain strings, numbers, nested
  maps/sequences). :func:`round_trip_equal` is provided to verify this for
  any given document rather than asserting it unconditionally.
- If a caller writes a ``domain`` dict that *itself* isn't representable
  in XDMF (e.g. an unsupported ``Topology`` type), that part is lossy by
  construction — this module only writes the ``Topology``/``Geometry``
  shapes it knows about (see :func:`build_xdmf_tree`).
"""

from __future__ import annotations

import base64
import json
from pathlib import Path
from typing import Any, Dict, Optional, Tuple
from xml.etree.ElementTree import Element, ElementTree, SubElement, parse as et_parse

from ymf.archive import canonicalize_domain

XDMF_HEADER = b'<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'

# Name used for the <Information Name="YMF" Value="..."> extension element
# that carries the non-mesh part of a YMF document through an XDMF file.
YMF_INFORMATION_NAME = "YMF"


def _indent_xml(elem: Element, level: int = 0) -> None:
    """In-place pretty-indent an ElementTree ``Element`` (no external deps)."""
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        child = None
        for child in elem:
            _indent_xml(child, level + 1)
            if not child.tail or not child.tail.strip():
                child.tail = i + "  "
        if child is not None and (not child.tail or not child.tail.strip()):
            child.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


XI_NAMESPACE = "http://www.w3.org/2001/XInclude"


def _data_item_attrs(data_item: Dict[str, Any]) -> Dict[str, str]:
    return {
        "Format": data_item["Format"],
        "DataType": data_item["DataType"],
        "Precision": str(data_item["Precision"]),
        "Dimensions": " ".join(str(d) for d in data_item["Dimensions"]),
    }


def _add_data_item(parent: Element, data_item: Dict[str, Any]) -> Element:
    """Write one DataItem, HDF5-referencing or text-including.

    ``Format="HDF"`` puts the dataset reference in the element's text.
    ``Format="XML"`` instead gets an ``<xi:include parse="text">`` child
    naming a sidecar file -- the no-HDF5 fallback.
    """
    elem = SubElement(parent, "DataItem", _data_item_attrs(data_item))
    if "Include" in data_item:
        SubElement(
            elem,
            "xi:include",
            {"parse": "text", "href": str(data_item["Include"])},
        )
    else:
        elem.text = str(data_item["Data"])
    return elem


def _parse_data_item(elem: Element) -> Dict[str, Any]:
    precision = elem.attrib.get("Precision", "4")
    item: Dict[str, Any] = {
        "Format": elem.attrib.get("Format", "HDF"),
        "DataType": elem.attrib.get("DataType", "Float"),
        "Precision": int(precision) if precision.isdigit() else precision,
        "Dimensions": [int(d) for d in elem.attrib["Dimensions"].split()],
    }
    # An xi:include child means the values live in a sidecar file. Match on
    # the local name so the lookup works whether the document used the
    # literal "xi:include" this module writes or a properly namespaced form
    # that a round-trip through another XML tool may have produced.
    include = None
    for child in elem:
        tag = child.tag
        if tag == "xi:include" or tag == f"{{{XI_NAMESPACE}}}include":
            include = child.attrib.get("href")
            break
    if include is not None:
        item["Include"] = include
    else:
        item["Data"] = (elem.text or "").strip()
    return item


def _add_grid_body(grid_elem: Element, grid: Dict[str, Any]) -> None:
    """Write a grid's Topology, Geometry and Attributes into ``grid_elem``."""
    topo = grid["Topology"]
    topo_attrs = {
        "Type": topo["Type"],
        "NumberOfElements": str(topo["NumberOfElements"]),
    }
    if "NodesPerElement" in topo:
        topo_attrs["NodesPerElement"] = str(topo["NodesPerElement"])
    _add_data_item(SubElement(grid_elem, "Topology", topo_attrs), topo["DataItem"])

    geom = grid["Geometry"]
    _add_data_item(SubElement(grid_elem, "Geometry", {"Type": geom["Type"]}), geom["DataItem"])

    for attr in grid.get("Attributes", []):
        attr_elem = SubElement(
            grid_elem,
            "Attribute",
            {
                "Name": attr["Name"],
                "AttributeType": attr["AttributeType"],
                "Center": attr["Center"],
            },
        )
        _add_data_item(attr_elem, attr["DataItem"])


def _parse_grid_body(grid_elem: Element) -> Dict[str, Any]:
    """Inverse of :func:`_add_grid_body`."""
    grid: Dict[str, Any] = {}

    topology_elem = grid_elem.find("Topology")
    if topology_elem is not None:
        topo: Dict[str, Any] = {
            "Type": topology_elem.attrib["Type"],
            "NumberOfElements": int(topology_elem.attrib["NumberOfElements"]),
            "DataItem": _parse_data_item(topology_elem.find("DataItem")),
        }
        if "NodesPerElement" in topology_elem.attrib:
            topo["NodesPerElement"] = int(topology_elem.attrib["NodesPerElement"])
        grid["Topology"] = topo

    geometry_elem = grid_elem.find("Geometry")
    if geometry_elem is not None:
        grid["Geometry"] = {
            "Type": geometry_elem.attrib["Type"],
            "DataItem": _parse_data_item(geometry_elem.find("DataItem")),
        }

    grid["Attributes"] = [
        {
            "Name": attr_elem.attrib["Name"],
            "AttributeType": attr_elem.attrib.get("AttributeType", "Scalar"),
            "Center": attr_elem.attrib.get("Center", "Node"),
            "DataItem": _parse_data_item(attr_elem.find("DataItem")),
        }
        for attr_elem in grid_elem.findall("Attribute")
    ]

    return grid


def parse_grid_element(grid_elem: Element) -> Dict[str, Any]:
    """Convert one XDMF ``<Grid>`` element into a grid dict.

    The inverse of what :func:`build_xdmf_tree` writes for a single uniform
    grid: ``Name`` (when present), ``Topology``, ``Geometry`` and
    ``Attributes``. A ``<Time>`` child is *not* included -- time belongs to
    the step that holds the grid, not to the grid itself.

    This exists for producers that still build XDMF elements directly and
    want to hand the result to the archive core as data. It is a bridge:
    a producer that builds grid dicts in the first place has no use for it.
    """
    grid: Dict[str, Any] = {}
    if "Name" in grid_elem.attrib:
        grid["Name"] = grid_elem.attrib["Name"]
    grid.update(_parse_grid_body(grid_elem))
    return grid


def grid_element_time(grid_elem: Element) -> Optional[float]:
    """Return the ``<Time>`` value of a grid element, or ``None``."""
    time_elem = grid_elem.find("Time")
    if time_elem is None:
        return None
    return float(time_elem.attrib["Value"])


def _encode_ymf_extra(ymf_extra: Any) -> str:
    """Serialize a YMF (sub-)document to base64-encoded JSON text.

    See the module docstring ("Round-trip strategy") for why JSON is used
    here specifically, rather than YAML/strictyaml as used everywhere else
    in the package.
    """
    json_text = json.dumps(ymf_extra, ensure_ascii=False)
    return base64.b64encode(json_text.encode("utf-8")).decode("ascii")


def _decode_ymf_extra(value: str) -> Any:
    """Inverse of :func:`_encode_ymf_extra`."""
    json_text = base64.b64decode(value.encode("ascii")).decode("utf-8")
    return json.loads(json_text)


def build_xdmf_tree(
    domain: Dict[str, Any], ymf_extra: Optional[Any] = None
) -> ElementTree:
    """Build an :class:`xml.etree.ElementTree.ElementTree` for an XDMF document.

    Parameters
    ----------
    domain:
        A plain dict shaped like the ``Domain`` block described in
        ``docs/ymf-schema.md`` §5, e.g.::

            {
                "TimeCollection": {
                    "Name": "Mesh Spatial_Domain",
                    "Data": [
                        {
                            "Time": 0.0,
                            "Topology": {"Type": "Tetrahedron",
                                         "NumberOfElements": 687737,
                                         "DataItem": {...}},
                            "Geometry": {"Type": "XYZ",
                                         "DataItem": {...}},
                        },
                        ...
                    ],
                }
            }
    ymf_extra:
        Optional additional YMF content with no XDMF equivalent (e.g. the
        ``Problem``/``solution_paths``/``vvuq`` sections of a full YMF
        document). If given, it's serialized to YAML, base64-encoded, and
        embedded as an ``<Information Name="YMF" Value="...">`` child of
        ``<Domain>`` so :func:`read_xdmf` can recover it later. See the
        module docstring for why this mechanism (rather than XML comments
        or custom elements) was chosen.

    Returns
    -------
    ElementTree
        Ready to ``.write(fileobj, encoding="utf-8")``. Callers typically
        prepend :data:`XDMF_HEADER` to the file before writing the tree,
        since ``ElementTree`` doesn't emit the XDMF DOCTYPE on its own.
    """
    root = Element(
        "Xdmf", {"Version": "2.0", "xmlns:xi": "http://www.w3.org/2001/XInclude"}
    )
    tree = ElementTree(root)
    domain_elem = SubElement(root, "Domain")

    if ymf_extra is not None:
        SubElement(
            domain_elem,
            "Information",
            {"Name": YMF_INFORMATION_NAME, "Value": _encode_ymf_extra(ymf_extra)},
        )

    # Canonicalize up front so this function reads one well-defined shape
    # rather than re-deriving defaults inline at every use.
    canonical = canonicalize_domain(domain)
    for time_collection in canonical.get("TimeCollections", []):
        collection_elem = SubElement(
            domain_elem,
            "Grid",
            {
                "Name": time_collection["Name"],
                "GridType": "Collection",
                "CollectionType": "Temporal",
            },
        )
        for i, step in enumerate(time_collection["Data"]):
            if "SpatialCollection" in step:
                # One grid per subdomain. Time hangs off the spatial
                # collection, not off each subdomain grid, so a viewer sees
                # a single instant made of several pieces.
                spatial_elem = SubElement(
                    collection_elem,
                    "Grid",
                    {"GridType": "Collection", "CollectionType": "Spatial"},
                )
                SubElement(
                    spatial_elem, "Time", {"Value": str(step["Time"]), "Name": str(i)}
                )
                for subdomain in step["SpatialCollection"]:
                    attrs = {"GridType": "Uniform"}
                    if "Name" in subdomain:
                        attrs["Name"] = subdomain["Name"]
                    _add_grid_body(SubElement(spatial_elem, "Grid", attrs), subdomain)
            else:
                grid_elem = SubElement(collection_elem, "Grid", {"GridType": "Uniform"})
                SubElement(
                    grid_elem, "Time", {"Value": str(step["Time"]), "Name": str(i)}
                )
                _add_grid_body(grid_elem, step)

    return tree


def write_xdmf(
    domain: Dict[str, Any], path: str | Path, ymf_extra: Optional[Any] = None
) -> None:
    """Write an XDMF 2.0 ``.xmf`` file for the given ``domain`` data tree.

    Convenience wrapper around :func:`build_xdmf_tree` that also writes the
    ``XDMF_HEADER`` (DOCTYPE) and pretty-indents the output. See
    :func:`build_xdmf_tree` for ``ymf_extra``.
    """
    tree = build_xdmf_tree(domain, ymf_extra=ymf_extra)
    _indent_xml(tree.getroot())
    path = Path(path)
    with open(path, "wb") as xml_file:
        xml_file.write(XDMF_HEADER)
        tree.write(xml_file, encoding="utf-8")


def parse_xdmf_domain(root: Element) -> Dict[str, Any]:
    """Parse an ``<Xdmf><Domain>...`` tree back into a ``domain`` dict.

    Inverse of the mesh-archive half of :func:`build_xdmf_tree` (i.e. of
    ``domain``, not ``ymf_extra`` — see :func:`read_xdmf` for the combined
    round-trip). Only understands the ``TimeCollection`` shape this module
    writes; other valid XDMF structures (non-temporal collections, multiple
    top-level grids, ``Attribute`` fields, etc.) are not yet parsed back —
    this is a round-trip for *what YMF itself writes*, not a general XDMF
    reader.
    """
    domain_elem = root.find("Domain")
    if domain_elem is None:
        return {}

    result: Dict[str, Any] = {}
    collections = []
    # A Domain holds one temporal collection per finite-element space -- a
    # Proteus archive commonly has both the linear base mesh and a
    # quadratic space -- so every child Grid is parsed, not just the first.
    for collection_elem in domain_elem.findall("Grid"):
        if collection_elem.attrib.get("CollectionType") != "Temporal":
            continue
        steps = []
        for step_elem in collection_elem.findall("Grid"):
            time_elem = step_elem.find("Time")
            step: Dict[str, Any] = {
                "Time": float(time_elem.attrib["Value"]) if time_elem is not None else 0.0,
            }
            if step_elem.attrib.get("CollectionType") == "Spatial":
                step["SpatialCollection"] = [
                    _parse_grid_body(sub_elem) for sub_elem in step_elem.findall("Grid")
                ]
                for sub_elem, parsed in zip(
                    step_elem.findall("Grid"), step["SpatialCollection"]
                ):
                    if "Name" in sub_elem.attrib:
                        # Name goes first, matching the key order the core's
                        # canonicalize_domain() produces, so a parsed
                        # document compares equal to a canonicalized one.
                        parsed_with_name = {"Name": sub_elem.attrib["Name"]}
                        parsed_with_name.update(parsed)
                        parsed.clear()
                        parsed.update(parsed_with_name)
            else:
                step.update(_parse_grid_body(step_elem))
            steps.append(step)
        collections.append(
            {
                "Name": collection_elem.attrib.get("Name", "TimeCollection"),
                "Data": steps,
            }
        )
    if collections:
        result["TimeCollections"] = collections

    return result


def parse_xdmf_extra(root: Element) -> Optional[Any]:
    """Extract and decode the ``ymf_extra`` payload from an XDMF tree, if present."""
    domain_elem = root.find("Domain")
    if domain_elem is None:
        return None
    for info_elem in domain_elem.findall("Information"):
        if info_elem.attrib.get("Name") == YMF_INFORMATION_NAME:
            return _decode_ymf_extra(info_elem.attrib["Value"])
    return None


def read_xdmf(path: str | Path) -> Tuple[Dict[str, Any], Optional[Any]]:
    """Read an ``.xmf`` file written by :func:`write_xdmf` back into Python.

    Returns a ``(domain, ymf_extra)`` tuple:

    - ``domain``: the mesh archive dict, in the same shape
      :func:`build_xdmf_tree` accepts (round-trips exactly for anything
      *this module* wrote — see :func:`parse_xdmf_domain`).
    - ``ymf_extra``: whatever was passed as ``ymf_extra`` to
      :func:`write_xdmf`/:func:`build_xdmf_tree`, decoded back from the
      embedded ``<Information Name="YMF">`` element, or ``None`` if the
      file has no such element (e.g. it's a plain XDMF file not written by
      this module, or was written with ``ymf_extra=None``).
    """
    tree = et_parse(path)
    root = tree.getroot()
    domain = parse_xdmf_domain(root)
    ymf_extra = parse_xdmf_extra(root)
    return domain, ymf_extra


def round_trip_equal(domain: Dict[str, Any], ymf_extra: Optional[Any] = None) -> bool:
    """Write ``(domain, ymf_extra)`` to a temp file, read it back, and compare.

    Useful for verifying a *specific* document round-trips exactly, rather
    than assuming it does (see the "Lossiness" note in the module
    docstring). Returns ``True`` iff both the parsed ``domain`` and
    ``ymf_extra`` are equal (``==``) to the originals after a full
    write/read cycle.
    """
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".xmf", delete=False) as f:
        tmp_path = f.name
    try:
        write_xdmf(domain, tmp_path, ymf_extra=ymf_extra)
        read_domain, read_extra = read_xdmf(tmp_path)
        return read_domain == canonicalize_domain(domain) and read_extra == ymf_extra
    finally:
        Path(tmp_path).unlink(missing_ok=True)
