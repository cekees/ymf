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

JSON rather than YAML is used for this embedded payload because it is
compact, in the standard library, and never read by eye: it is base64
text inside an XML attribute. The ``.ymf`` archive itself holds the same
content as readable YAML.

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
from typing import Any, Dict, List, Optional, Tuple
from xml.etree.ElementTree import Element, ElementTree, ParseError, SubElement, parse as et_parse

from ymf.archive import YmfArchiveError, canonicalize_domain

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


# ---------------------------------------------------------------------------
# reading
#
# The reader accepts what XDMF allows wherever the archive model can hold
# it, and refuses the rest by name. It never drops or guesses: an XDMF
# construct the model cannot represent raises YmfArchiveError naming the
# element, rather than coming back empty or subtly wrong. See
# docs/xdmf-model.yaml for the construct-by-construct correspondence.
# ---------------------------------------------------------------------------

#: XDMF 3 renamed several attributes. The reader accepts either name; the
#: writer emits the XDMF 2 one (first in each pair is the XDMF 3 name).
_TOPOLOGY_TYPE = ("TopologyType", "Type")
_GEOMETRY_TYPE = ("GeometryType", "Type")
_NUMBER_TYPE = ("NumberType", "DataType")

#: Topologies defined by Dimensions alone, with implicit connectivity. A
#: YMF topology always has a connectivity DataItem, so these can't be read.
_STRUCTURED_TOPOLOGIES = frozenset(
    {"2DSMesh", "2DRectMesh", "2DCoRectMesh", "3DSMesh", "3DRectMesh", "3DCoRectMesh"})

_REQUIRED = object()


def _attr(elem: Element, names: Tuple[str, ...], where: str, default: Any = _REQUIRED) -> Any:
    """The first of ``names`` present on ``elem``, else ``default``.

    Raises if neither is available and no default was given.
    """
    for name in names:
        if name in elem.attrib:
            return elem.attrib[name]
    if default is _REQUIRED:
        raise YmfArchiveError("%s: missing %s" % (where, " / ".join(names)))
    return default


def _refuse_reference(elem: Element, where: str) -> None:
    if "Reference" in elem.attrib:
        raise YmfArchiveError(
            "%s: Reference=%r points at another element; a YMF archive has no "
            "references (XPath), so the referenced element has to be written "
            "out in place" % (where, elem.attrib["Reference"]))


def _only_data_item(elem: Element, where: str, what: str) -> Element:
    """The single DataItem child of ``elem``, refusing none or several."""
    items = elem.findall("DataItem")
    if not items:
        raise YmfArchiveError("%s: has no DataItem" % (where,))
    if len(items) > 1:
        raise YmfArchiveError(
            "%s: has %d DataItems; a YMF %s holds exactly one"
            % (where, len(items), what))
    return items[0]


def _data_item_attrs(data_item: Dict[str, Any]) -> Dict[str, str]:
    return {
        "Format": data_item["Format"],
        "DataType": data_item["DataType"],
        "Precision": str(data_item["Precision"]),
        "Dimensions": " ".join(str(d) for d in data_item["Dimensions"]),
    }


def _format_values(values: Any, dimensions: Any) -> str:
    """Inline values as DataItem text: one row (last dimension) per line.

    Floats use ``repr``, the shortest text that reads back to the same
    double, so inline values round-trip exactly.
    """
    row = int(dimensions[-1]) if len(dimensions) > 1 and int(dimensions[-1]) > 0 else len(values) or 1
    text = [repr(v) if isinstance(v, float) else str(v) for v in values]
    lines = [" ".join(text[i:i + row]) for i in range(0, len(text), row)]
    return "\n" + "\n".join(lines) + "\n"


def _add_data_item(parent: Element, data_item: Dict[str, Any]) -> Element:
    """Write one DataItem: an HDF5 reference, a text sidecar, or inline values.

    ``Format="HDF"`` puts the dataset reference in the element's text.
    ``Format="XML"`` either gets an ``<xi:include parse="text">`` child
    naming a sidecar file, or, for inline ``Values``, the values themselves
    as the element's text -- a self-contained ``.xmf``.
    """
    elem = SubElement(parent, "DataItem", _data_item_attrs(data_item))
    if "Values" in data_item:
        elem.text = _format_values(data_item["Values"], data_item["Dimensions"])
    elif "Include" in data_item:
        SubElement(
            elem,
            "xi:include",
            {"parse": "text", "href": str(data_item["Include"])},
        )
    else:
        elem.text = str(data_item["Data"])
    return elem


def _parse_data_item(elem: Element, where: str = "DataItem") -> Dict[str, Any]:
    """Read one DataItem element.

    A uniform item, in any of the forms the archive model holds: an HDF5
    reference (``Format="HDF"``, the reference as element text), a text
    sidecar (``Format="XML"`` with an ``xi:include`` child), values inline
    (``Format="XML"`` -- XDMF's default -- with the numbers as element
    text), or a Binary file name. Inline values are parsed by DataType and
    checked against ``Dimensions``.
    """
    _refuse_reference(elem, where)
    item_type = _attr(elem, ("ItemType", "Type"), where, "Uniform")
    if item_type.lower() != "uniform":
        raise YmfArchiveError(
            "%s: ItemType=%r; a YMF archive holds Uniform DataItems only "
            "(no HyperSlab, Coordinate, Function, Collection or Tree)" % (where, item_type))
    precision = elem.attrib.get("Precision", "4")
    fmt = elem.attrib.get("Format", "XML")  # XDMF's default
    item: Dict[str, Any] = {
        "Format": fmt,
        "DataType": _attr(elem, _NUMBER_TYPE, where, "Float"),
        "Precision": int(precision) if precision.isdigit() else precision,
        "Dimensions": [int(d) for d in _attr(elem, ("Dimensions",), where).split()],
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
    elif fmt == "XML":
        item["Values"] = _parse_values(elem.text or "", item, where)
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


def _parse_values(text: str, item: Dict[str, Any], where: str) -> List[Any]:
    """Inline DataItem text as a flat list of numbers, checked for count."""
    integral = item["DataType"] in ("Int", "UInt", "Char", "UChar")
    tokens = text.split()
    try:
        values = [int(t) if integral else float(t) for t in tokens]
    except ValueError as exc:
        raise YmfArchiveError(
            "%s: inline values are not all %s numbers (%s)"
            % (where, item["DataType"], exc)) from exc
    expected = 1
    for d in item["Dimensions"]:
        expected *= d
    if len(values) != expected:
        raise YmfArchiveError(
            "%s: Dimensions %s declare %d values, but %d are inline"
            % (where, item["Dimensions"], expected, len(values)))
    return values


def _parse_grid_body(grid_elem: Element, where: str = "Grid") -> Dict[str, Any]:
    """Inverse of :func:`_add_grid_body`."""
    _refuse_reference(grid_elem, where)
    grid: Dict[str, Any] = {}

    topology_elem = grid_elem.find("Topology")
    if topology_elem is not None:
        at = where + ".Topology"
        _refuse_reference(topology_elem, at)
        topo_type = _attr(topology_elem, _TOPOLOGY_TYPE, at)
        if topo_type in _STRUCTURED_TOPOLOGIES:
            raise YmfArchiveError(
                "%s: %s is a structured topology, defined by Dimensions with "
                "implicit connectivity; a YMF topology is unstructured, with a "
                "connectivity DataItem" % (at, topo_type))
        # XDMF allows Dimensions in place of NumberOfElements; for an
        # unstructured topology it is the element count.
        count = _attr(topology_elem, ("NumberOfElements", "Dimensions"), at)
        if len(count.split()) != 1:
            raise YmfArchiveError(
                "%s: element count %r is not a single number" % (at, count))
        topo: Dict[str, Any] = {
            "Type": topo_type,
            "NumberOfElements": int(count),
            "DataItem": _parse_data_item(
                _only_data_item(topology_elem, at, "Topology"), at + ".DataItem"),
        }
        if "NodesPerElement" in topology_elem.attrib:
            topo["NodesPerElement"] = int(topology_elem.attrib["NodesPerElement"])
        grid["Topology"] = topo

    geometry_elem = grid_elem.find("Geometry")
    if geometry_elem is not None:
        at = where + ".Geometry"
        _refuse_reference(geometry_elem, at)
        grid["Geometry"] = {
            "Type": _attr(geometry_elem, _GEOMETRY_TYPE, at, "XYZ"),
            "DataItem": _parse_data_item(
                _only_data_item(geometry_elem, at, "Geometry"), at + ".DataItem"),
        }

    attributes = []
    for i, attr_elem in enumerate(grid_elem.findall("Attribute")):
        at = "%s.Attribute[%s]" % (where, attr_elem.attrib.get("Name", i))
        _refuse_reference(attr_elem, at)
        if "ItemType" in attr_elem.attrib:
            raise YmfArchiveError(
                "%s: ItemType=%r; a YMF attribute is one array of values, not a "
                "finite-element function (dofmap plus values)"
                % (at, attr_elem.attrib["ItemType"]))
        attributes.append({
            "Name": _attr(attr_elem, ("Name",), at),
            "AttributeType": attr_elem.attrib.get("AttributeType", "Scalar"),
            "Center": attr_elem.attrib.get("Center", "Node"),
            "DataItem": _parse_data_item(
                _only_data_item(attr_elem, at, "Attribute"), at + ".DataItem"),
        })
    grid["Attributes"] = attributes

    for unsupported in ("Set",):
        if grid_elem.find(unsupported) is not None:
            raise YmfArchiveError(
                "%s: has a %s; a YMF archive has no sets (store region markers "
                "as Cell- or Node-centred attributes)" % (where, unsupported))

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
    grid.update(_parse_grid_body(grid_elem, "Grid[%s]" % grid_elem.attrib.get("Name", "")))
    return grid


def grid_element_time(grid_elem: Element, where: str = "Grid") -> Optional[float]:
    """Return the ``<Time>`` value of a grid element, or ``None``.

    Only ``TimeType="Single"`` (XDMF's default) can be held: a step in a
    YMF archive is one instant. List, HyperSlab and Range are refused.
    """
    time_elem = grid_elem.find("Time")
    if time_elem is None:
        return None
    time_type = time_elem.attrib.get("TimeType", "Single")
    if time_type != "Single":
        raise YmfArchiveError(
            "%s.Time: TimeType=%r; a YMF step holds a single time value"
            % (where, time_type))
    return float(_attr(time_elem, ("Value",), where + ".Time"))


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
        document). If given, it's serialized to JSON, base64-encoded, and
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


def _grid_kind(grid_elem: Element) -> str:
    """``"uniform"``, ``"temporal"`` or ``"spatial"``; refuses the rest."""
    grid_type = grid_elem.attrib.get("GridType", "Uniform")
    if grid_type == "Uniform":
        return "uniform"
    if grid_type == "Collection":
        # XDMF's default CollectionType is Spatial.
        return grid_elem.attrib.get("CollectionType", "Spatial").lower()
    raise YmfArchiveError(
        "Grid %r: GridType=%r; a YMF archive holds Uniform grids and "
        "Temporal or Spatial collections only"
        % (grid_elem.attrib.get("Name", ""), grid_type))


def _parse_step(step_elem: Element, where: str) -> Dict[str, Any]:
    """One instant: a uniform grid, or a spatial collection of them."""
    kind = _grid_kind(step_elem)
    time = grid_element_time(step_elem, where)
    step: Dict[str, Any] = {"Time": time if time is not None else 0.0}
    if kind == "uniform":
        step.update(_parse_grid_body(step_elem, where))
    elif kind == "spatial":
        subgrids = []
        for j, sub_elem in enumerate(step_elem.findall("Grid")):
            at = "%s.Grid[%d]" % (where, j)
            if _grid_kind(sub_elem) != "uniform":
                raise YmfArchiveError(
                    "%s: a spatial collection in a YMF archive holds uniform grids "
                    "only, not further collections" % (at,))
            # Name goes first, matching the key order canonicalize_domain()
            # produces, so a parsed document compares equal to a canonical one.
            sub: Dict[str, Any] = {}
            if "Name" in sub_elem.attrib:
                sub["Name"] = sub_elem.attrib["Name"]
            sub.update(_parse_grid_body(sub_elem, at))
            subgrids.append(sub)
        if not subgrids:
            raise YmfArchiveError("%s: empty spatial collection" % (where,))
        step["SpatialCollection"] = subgrids
    else:
        raise YmfArchiveError(
            "%s: a temporal collection inside a time step; a YMF step is one "
            "instant" % (where,))
    return step


def parse_xdmf_domain(root: Element) -> Dict[str, Any]:
    """Parse an ``<Xdmf><Domain>...`` tree back into a ``domain`` dict.

    The inverse of the mesh-archive half of :func:`build_xdmf_tree`, and a
    reader for XDMF written by other tools, as far as the archive model
    reaches. Each top-level ``<Grid>`` of the Domain becomes one time
    collection:

    - a Temporal collection: its steps, as written;
    - a Uniform grid, or a Spatial collection: a collection of a single
      step, at the grid's ``<Time>`` if it has one and 0.0 otherwise.

    A collection takes the grid's ``Name``; unnamed ones are numbered.
    Anything the archive model cannot hold (Tree or Subset grids, Sets,
    references, non-uniform DataItems, inline values, ...) raises
    :class:`ymf.archive.YmfArchiveError` naming the element, so a file is
    never read as emptier, or different, than it is.
    """
    domain_elem = root.find("Domain")
    if domain_elem is None:
        return {}

    collections = []
    # A Domain holds one temporal collection per finite-element space -- a
    # Proteus archive commonly has both the linear base mesh and a
    # quadratic space -- so every child Grid is parsed, not just the first.
    for ci, collection_elem in enumerate(domain_elem.findall("Grid")):
        name = collection_elem.attrib.get("Name")
        where = "Domain.Grid[%s]" % (name if name is not None else ci)
        kind = _grid_kind(collection_elem)
        if kind == "temporal":
            steps = [
                _parse_step(step_elem, "%s.Grid[%d]" % (where, i))
                for i, step_elem in enumerate(collection_elem.findall("Grid"))
            ]
            default_name = "TimeCollection"
        else:
            steps = [_parse_step(collection_elem, where)]
            default_name = "Grid %d" % (ci,)
        collections.append({
            "Name": name if name is not None else default_name,
            "Data": steps,
        })

    return {"TimeCollections": collections} if collections else {}


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
    """Read an ``.xmf`` file back into Python.

    Files :func:`write_xdmf` wrote round-trip exactly. XDMF from other
    tools is read as far as the archive model reaches, and anything beyond
    it raises :class:`ymf.archive.YmfArchiveError` naming the element; see
    :func:`parse_xdmf_domain`.

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
    try:
        tree = et_parse(path)
    except ParseError as exc:
        raise YmfArchiveError("%s is not an XML document (%s)" % (path, exc)) from exc
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
