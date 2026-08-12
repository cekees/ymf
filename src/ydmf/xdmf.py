"""Convert a YDMF ``Archive``/``Domain`` data tree to an XDMF 2.0 XML document.

Per ``docs/ydmf-schema.md`` §5 ("XDMF Preservation") and the design
principle that *XDMF is a consumer, not a dependency* of YDMF: this module
reads a plain Python dict describing a time-collection of grids (mesh
topology + geometry, optionally HDF5-backed) and emits an ``.xmf`` file
that ParaView and other XDMF-aware tools can open directly.

This generalizes the prototype conversion logic originally written
ad hoc in ``examples/notebooks/YDMF-original.ipynb`` (which used a
bespoke ``Model``/``ModelDomain``/``Domain`` strictyaml schema that
predates the current ``ydmf.schema`` package). The updated notebook,
``examples/notebooks/YDMF-updated.ipynb``, uses this module instead of
its own private XML-writing code.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict
from xml.etree.ElementTree import Element, ElementTree, SubElement

XDMF_HEADER = (
    b'<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'
)


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


def _data_item_attrs(data_item: Dict[str, Any]) -> Dict[str, str]:
    dims = " ".join(str(d) for d in data_item["Dimensions"])
    attrs = {
        "Format": data_item.get("Format", "HDF"),
        "DataType": data_item.get("DataType", "Float"),
        "Precision": str(data_item.get("Precision", "4")),
        "Dimensions": dims,
    }
    return attrs


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


def build_xdmf_tree(domain: Dict[str, Any]) -> ElementTree:
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


def write_xdmf(domain: Dict[str, Any], path: str | Path) -> None:
    """Write an XDMF 2.0 ``.xmf`` file for the given ``domain`` data tree.

    Convenience wrapper around :func:`build_xdmf_tree` that also writes the
    ``XDMF_HEADER`` (DOCTYPE) and pretty-indents the output.
    """
    tree = build_xdmf_tree(domain)
    _indent_xml(tree.getroot())
    path = Path(path)
    with open(path, "wb") as xml_file:
        xml_file.write(XDMF_HEADER)
        tree.write(xml_file, encoding="utf-8")
