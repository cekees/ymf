"""Convert between a YDMF data tree and an XDMF 2.0 XML document, round-trip.

Per ``docs/ydmf-schema.md`` §5 ("XDMF Preservation") and the design
principle that *XDMF is a consumer, not a dependency* of YDMF: this module
reads a plain Python dict describing a time-collection of grids (mesh
topology + geometry, optionally HDF5-backed) and emits an ``.xmf`` file
that ParaView and other XDMF-aware tools can open directly, using only the
parts of the YDMF document that XDMF actually understands.

Round-trip strategy (write -> read -> get the same YDMF document back)
-----------------------------------------------------------------------
XDMF's grammar (the ``Xdmf.dtd``/``Xdmf.xsd``, referenced in
``examples/notebooks/YDMF-original.ipynb`` cell 1) defines an
``<Information Name="..." Value="...">`` element that is legal as a child
of essentially any XDMF element (``Domain``, ``Grid``, ``DataItem``, ...)
specifically as an extension point: generic XDMF consumers (ParaView,
VisIt) parse and then ignore ``Information`` elements they don't
recognize, but the element itself survives round-trips through any
XDMF-compliant XML tool, unlike XML comments (which many parsers,
including a plain ``xml.etree.ElementTree`` reader with default settings,
silently drop) or non-standard custom elements (which a strict
schema-validating XDMF reader could legitimately reject).

So: everything in a YDMF document that *does* map onto XDMF concepts
(``Domain``/``Grid``/``Topology``/``Geometry``/``DataItem`` — i.e. the mesh
archive) is written as real XDMF elements. Everything else (``Problem``,
``solution_paths``, ``vvuq``, and any ``archive``-level metadata XDMF has
no concept of) is serialized to JSON, base64-encoded (safe inside an XML
attribute value with no escaping headaches), and stashed in a single
``<Information Name="YDMF" Value="...">`` child of ``<Domain>``.

JSON rather than YAML is used for this *specific* embedded payload:
strictyaml's ``as_document()`` cannot serialize empty lists/dicts without
an explicit schema (raises ``YAMLSerializationError`` on
``solution_paths.analytical: []``, which every current example uses), and
since this payload is base64-encoded anyway — not meant to be read
directly out of the XML — YAML's human-authoring niceties (comments, flow
style, block scalars) buy nothing here. The rest of the package still uses
YAML/strictyaml everywhere a human or LLM actually edits a YDMF document;
only this internal round-trip encoding uses JSON.

**Honesty note**: this is a reasonable, standards-consistent design given
what's known about XDMF's grammar, but it has *not* been verified against
the live ``Xdmf.dtd``/XSD from gitlab.kitware.com — this sandbox has no
network access to fetch it. If ``Information`` turns out to have stricter
placement rules than assumed here, the fallback is to move the element
under a ``<Grid>`` instead of directly under ``<Domain>`` (also permitted
by every XDMF version this design is aware of).

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

XDMF_HEADER = b'<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'

# Name used for the <Information Name="YDMF" Value="..."> extension element
# that carries the non-mesh part of a YDMF document through an XDMF file.
YDMF_INFORMATION_NAME = "YDMF"


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


def _canonicalize_data_item(data_item: Dict[str, Any]) -> Dict[str, Any]:
    """Fill in XDMF's documented DataItem defaults for missing fields.

    Once a DataItem is written to XML, every attribute has *some* value --
    there's no XML-level distinction between "attribute omitted, consumer
    should assume the default" and "attribute explicitly set to the
    default". So a sparse input dict (e.g. no ``Precision`` given) and the
    dict parsed back after writing will only be equal if the sparse dict is
    first canonicalized to the same defaults the writer used. This helper
    does that canonicalization so callers/tests can compare against a
    well-defined expected form rather than being surprised by it.
    """
    return {
        "Format": data_item.get("Format", "HDF"),
        "DataType": data_item.get("DataType", "Float"),
        "Precision": int(data_item.get("Precision", 4)),
        "Dimensions": list(data_item["Dimensions"]),
        "Data": str(data_item["Data"]),
    }


def canonicalize_domain(domain: Dict[str, Any]) -> Dict[str, Any]:
    """Return ``domain`` with every DataItem's defaultable fields filled in.

    Use this when comparing a hand-written sparse ``domain`` dict against
    the result of :func:`read_xdmf` on a file written from it -- the
    written/read-back form always has ``Format``/``DataType``/``Precision``
    explicit (see :func:`_canonicalize_data_item`), so a direct ``==``
    against the original sparse dict can spuriously fail even though
    nothing was actually lost.
    """
    result: Dict[str, Any] = {}
    time_collection = domain.get("TimeCollection")
    if time_collection is not None:
        result["TimeCollection"] = {
            "Name": time_collection.get("Name", "TimeCollection"),
            "Data": [
                {
                    "Time": float(grid["Time"]),
                    "Topology": {
                        "Type": grid["Topology"]["Type"],
                        "NumberOfElements": int(grid["Topology"]["NumberOfElements"]),
                        "DataItem": _canonicalize_data_item(grid["Topology"]["DataItem"]),
                    },
                    "Geometry": {
                        "Type": grid["Geometry"]["Type"],
                        "DataItem": _canonicalize_data_item(grid["Geometry"]["DataItem"]),
                    },
                }
                for grid in time_collection["Data"]
            ],
        }
    return result


def _data_item_attrs(data_item: Dict[str, Any]) -> Dict[str, str]:
    canonical = _canonicalize_data_item(data_item)
    return {
        "Format": canonical["Format"],
        "DataType": canonical["DataType"],
        "Precision": str(canonical["Precision"]),
        "Dimensions": " ".join(str(d) for d in canonical["Dimensions"]),
    }


def _parse_data_item(elem: Element) -> Dict[str, Any]:
    dims = [int(d) for d in elem.attrib["Dimensions"].split()]
    precision = elem.attrib.get("Precision", "4")
    return {
        "Format": elem.attrib.get("Format", "HDF"),
        "DataType": elem.attrib.get("DataType", "Float"),
        "Precision": int(precision) if precision.isdigit() else precision,
        "Dimensions": dims,
        "Data": (elem.text or "").strip(),
    }


def _add_topology_and_geometry(grid_elem: Element, grid: Dict[str, Any]) -> None:
    topo = grid["Topology"]
    topology_elem = SubElement(
        grid_elem,
        "Topology",
        {"Type": topo["Type"], "NumberOfElements": str(topo["NumberOfElements"])},
    )
    d = topo["DataItem"]
    data_item = SubElement(topology_elem, "DataItem", _data_item_attrs(d))
    data_item.text = str(d["Data"])

    geom = grid["Geometry"]
    geometry_elem = SubElement(grid_elem, "Geometry", {"Type": geom["Type"]})
    d = geom["DataItem"]
    data_item = SubElement(geometry_elem, "DataItem", _data_item_attrs(d))
    data_item.text = str(d["Data"])


def _parse_topology_and_geometry(grid_elem: Element) -> Dict[str, Any]:
    grid: Dict[str, Any] = {}

    topology_elem = grid_elem.find("Topology")
    if topology_elem is not None:
        data_item_elem = topology_elem.find("DataItem")
        grid["Topology"] = {
            "Type": topology_elem.attrib["Type"],
            "NumberOfElements": int(topology_elem.attrib["NumberOfElements"]),
            "DataItem": _parse_data_item(data_item_elem),
        }

    geometry_elem = grid_elem.find("Geometry")
    if geometry_elem is not None:
        data_item_elem = geometry_elem.find("DataItem")
        grid["Geometry"] = {
            "Type": geometry_elem.attrib["Type"],
            "DataItem": _parse_data_item(data_item_elem),
        }

    return grid


def _encode_ydmf_extra(ydmf_extra: Any) -> str:
    """Serialize a YDMF (sub-)document to base64-encoded JSON text.

    See the module docstring ("Round-trip strategy") for why JSON is used
    here specifically, rather than YAML/strictyaml as used everywhere else
    in the package.
    """
    json_text = json.dumps(ydmf_extra, ensure_ascii=False)
    return base64.b64encode(json_text.encode("utf-8")).decode("ascii")


def _decode_ydmf_extra(value: str) -> Any:
    """Inverse of :func:`_encode_ydmf_extra`."""
    json_text = base64.b64decode(value.encode("ascii")).decode("utf-8")
    return json.loads(json_text)


def build_xdmf_tree(
    domain: Dict[str, Any], ydmf_extra: Optional[Any] = None
) -> ElementTree:
    """Build an :class:`xml.etree.ElementTree.ElementTree` for an XDMF document.

    Parameters
    ----------
    domain:
        A plain dict shaped like the ``Domain`` block described in
        ``docs/ydmf-schema.md`` §5, e.g.::

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
    ydmf_extra:
        Optional additional YDMF content with no XDMF equivalent (e.g. the
        ``Problem``/``solution_paths``/``vvuq`` sections of a full YDMF
        document). If given, it's serialized to YAML, base64-encoded, and
        embedded as an ``<Information Name="YDMF" Value="...">`` child of
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

    if ydmf_extra is not None:
        SubElement(
            domain_elem,
            "Information",
            {"Name": YDMF_INFORMATION_NAME, "Value": _encode_ydmf_extra(ydmf_extra)},
        )

    time_collection = domain.get("TimeCollection")
    if time_collection is not None:
        collection_elem = SubElement(
            domain_elem,
            "Grid",
            {
                "Name": time_collection.get("Name", "TimeCollection"),
                "GridType": "Collection",
                "CollectionType": "Temporal",
            },
        )
        for i, grid in enumerate(time_collection["Data"]):
            grid_elem = SubElement(collection_elem, "Grid", {"GridType": "Uniform"})
            SubElement(
                grid_elem,
                "Time",
                {"Value": str(grid["Time"]), "Name": str(i)},
            )
            _add_topology_and_geometry(grid_elem, grid)

    return tree


def write_xdmf(
    domain: Dict[str, Any], path: str | Path, ydmf_extra: Optional[Any] = None
) -> None:
    """Write an XDMF 2.0 ``.xmf`` file for the given ``domain`` data tree.

    Convenience wrapper around :func:`build_xdmf_tree` that also writes the
    ``XDMF_HEADER`` (DOCTYPE) and pretty-indents the output. See
    :func:`build_xdmf_tree` for ``ydmf_extra``.
    """
    tree = build_xdmf_tree(domain, ydmf_extra=ydmf_extra)
    _indent_xml(tree.getroot())
    path = Path(path)
    with open(path, "wb") as xml_file:
        xml_file.write(XDMF_HEADER)
        tree.write(xml_file, encoding="utf-8")


def parse_xdmf_domain(root: Element) -> Dict[str, Any]:
    """Parse an ``<Xdmf><Domain>...`` tree back into a ``domain`` dict.

    Inverse of the mesh-archive half of :func:`build_xdmf_tree` (i.e. of
    ``domain``, not ``ydmf_extra`` — see :func:`read_xdmf` for the combined
    round-trip). Only understands the ``TimeCollection`` shape this module
    writes; other valid XDMF structures (non-temporal collections, multiple
    top-level grids, ``Attribute`` fields, etc.) are not yet parsed back —
    this is a round-trip for *what YDMF itself writes*, not a general XDMF
    reader.
    """
    domain_elem = root.find("Domain")
    if domain_elem is None:
        return {}

    result: Dict[str, Any] = {}
    collection_elem = domain_elem.find("Grid")
    if collection_elem is not None and collection_elem.attrib.get(
        "CollectionType"
    ) == "Temporal":
        grids = []
        for grid_elem in collection_elem.findall("Grid"):
            time_elem = grid_elem.find("Time")
            grid_data = {
                "Time": float(time_elem.attrib["Value"]) if time_elem is not None else 0.0,
            }
            grid_data.update(_parse_topology_and_geometry(grid_elem))
            grids.append(grid_data)
        result["TimeCollection"] = {
            "Name": collection_elem.attrib.get("Name", "TimeCollection"),
            "Data": grids,
        }

    return result


def parse_xdmf_extra(root: Element) -> Optional[Any]:
    """Extract and decode the ``ydmf_extra`` payload from an XDMF tree, if present."""
    domain_elem = root.find("Domain")
    if domain_elem is None:
        return None
    for info_elem in domain_elem.findall("Information"):
        if info_elem.attrib.get("Name") == YDMF_INFORMATION_NAME:
            return _decode_ydmf_extra(info_elem.attrib["Value"])
    return None


def read_xdmf(path: str | Path) -> Tuple[Dict[str, Any], Optional[Any]]:
    """Read an ``.xmf`` file written by :func:`write_xdmf` back into Python.

    Returns a ``(domain, ydmf_extra)`` tuple:

    - ``domain``: the mesh archive dict, in the same shape
      :func:`build_xdmf_tree` accepts (round-trips exactly for anything
      *this module* wrote — see :func:`parse_xdmf_domain`).
    - ``ydmf_extra``: whatever was passed as ``ydmf_extra`` to
      :func:`write_xdmf`/:func:`build_xdmf_tree`, decoded back from the
      embedded ``<Information Name="YDMF">`` element, or ``None`` if the
      file has no such element (e.g. it's a plain XDMF file not written by
      this module, or was written with ``ydmf_extra=None``).
    """
    tree = et_parse(path)
    root = tree.getroot()
    domain = parse_xdmf_domain(root)
    ydmf_extra = parse_xdmf_extra(root)
    return domain, ydmf_extra


def round_trip_equal(domain: Dict[str, Any], ydmf_extra: Optional[Any] = None) -> bool:
    """Write ``(domain, ydmf_extra)`` to a temp file, read it back, and compare.

    Useful for verifying a *specific* document round-trips exactly, rather
    than assuming it does (see the "Lossiness" note in the module
    docstring). Returns ``True`` iff both the parsed ``domain`` and
    ``ydmf_extra`` are equal (``==``) to the originals after a full
    write/read cycle.
    """
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".xmf", delete=False) as f:
        tmp_path = f.name
    try:
        write_xdmf(domain, tmp_path, ydmf_extra=ydmf_extra)
        read_domain, read_extra = read_xdmf(tmp_path)
        return read_domain == canonicalize_domain(domain) and read_extra == ydmf_extra
    finally:
        Path(tmp_path).unlink(missing_ok=True)
