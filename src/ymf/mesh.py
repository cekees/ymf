"""Helpers for writing and checking a mesh stored as a ymf grid.

A mesh is a grid (``ymf.archive``): ``topology``, ``geometry`` and
``attributes``, plus two optional keys that only a mesh needs -- its tagged
``facets`` and its ``provenance``. The schema, and the cheap structural
checks, live in ``ymf.archive``, which depends on pyyaml alone; the helpers
here read files and arrays, so they import numpy (and h5py, to read a
referenced array) when called.

Nothing here names an external mesh format. A converter from some code's
input files to ymf states that code's conventions as data -- for instance
the vertex pairs its quadratic cells' edge nodes sit on -- and these
helpers do the rest.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from .archive import (
    QUADRATIC_EDGES,
    YmfArchiveError,
    check_facet_tags,
    load_array,
    provenance,
)

__all__ = [
    "file_sha256",
    "file_provenance",
    "edge_node_permutation",
    "reorder_cells",
    "edge_midpoint_offsets",
    "check_facets",
]

#: The faces of each linear 3D cell and the edges of each linear 2D cell, as
#: tuples of the cell's vertex numbers (VTK's order, as for the node order
#: in :data:`ymf.archive.QUADRATIC_EDGES`). A facet is matched against them
#: as a set of vertices, so the order within a face does not matter here.
_CELL_FACETS = {
    "Polyline": ((0,), (1,)),
    "Triangle": ((0, 1), (1, 2), (2, 0)),
    "Quadrilateral": ((0, 1), (1, 2), (2, 3), (3, 0)),
    "Tetrahedron": ((0, 1, 3), (1, 2, 3), (2, 0, 3), (0, 2, 1)),
    "Hexahedron": ((0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4), (3, 7, 6, 2),
                   (0, 3, 2, 1), (4, 5, 6, 7)),
    "Wedge": ((0, 1, 2), (3, 5, 4), (0, 3, 4, 1), (1, 4, 5, 2), (2, 5, 3, 0)),
    "Pyramid": ((0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)),
}
_LINEAR = {
    "Edge_3": "Polyline", "Tri_6": "Triangle", "Triangle_6": "Triangle",
    "Quad_8": "Quadrilateral", "Quadrilateral_8": "Quadrilateral",
    "Tet_10": "Tetrahedron", "Tetrahedron_10": "Tetrahedron",
    "Pyramid_13": "Pyramid", "Wedge_15": "Wedge",
    "Hex_20": "Hexahedron", "Hexahedron_20": "Hexahedron",
}
#: Vertices of a facet of each type (the leading nodes of a quadratic facet).
_FACET_VERTICES = {"Polyvertex": 1, "Polyline": 2, "Edge_3": 2, "Triangle": 3,
                   "Tri_6": 3, "Triangle_6": 3, "Quadrilateral": 4, "Quad_8": 4,
                   "Quadrilateral_8": 4}


def file_sha256(path: str | Path) -> str:
    """The SHA-256 digest of a file's bytes, as 64 lowercase hex digits."""
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_provenance(
    tool: str,
    paths: Iterable[str | Path],
    *,
    version: Optional[str] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """A :func:`ymf.archive.provenance` block for a mesh made from ``paths``.

    Each file is recorded by its base name and the SHA-256 of its bytes, so
    a reader can tell whether the mesh still matches the files beside it
    without trusting their timestamps.
    """
    files = [{"name": Path(p).name, "sha256": file_sha256(p)} for p in paths]
    return provenance(tool, version=version, files=files, extra=extra)


def edge_node_permutation(topology_type: str,
                          source_edges: Sequence[Tuple[int, int]]) -> List[int]:
    """The column order that takes a source's quadratic cells to ymf's.

    ``source_edges`` lists, for the source's edge nodes in its own order,
    the pair of vertices each sits on. The result ``perm`` gives
    ``ymf_cells = source_cells[:, perm]``: vertices unchanged, edge nodes
    reordered to :data:`ymf.archive.QUADRATIC_EDGES`. An edge is the same
    whichever way round its pair is written.
    """
    if topology_type not in QUADRATIC_EDGES:
        raise YmfArchiveError(
            f"{topology_type!r} is not a quadratic cell type; quadratic types are "
            f"{sorted(QUADRATIC_EDGES)}")
    target = QUADRATIC_EDGES[topology_type]
    n_vertices = max(max(e) for e in target) + 1
    source = [frozenset(e) for e in source_edges]
    if len(source) != len(target) or len(set(source)) != len(source):
        raise YmfArchiveError(
            f"{topology_type} has {len(target)} distinct edges; the source lists "
            f"{list(source_edges)}")
    perm = list(range(n_vertices))
    for edge in target:
        try:
            perm.append(n_vertices + source.index(frozenset(edge)))
        except ValueError:
            raise YmfArchiveError(
                f"{topology_type}: edge {edge} is not among the source's edges "
                f"{list(source_edges)}") from None
    return perm


def reorder_cells(cells, topology_type: str, source_edges: Sequence[Tuple[int, int]]):
    """``cells`` (one row per cell) with the edge nodes put in ymf's order."""
    import numpy
    return numpy.asarray(cells)[:, edge_node_permutation(topology_type, source_edges)]


def edge_midpoint_offsets(topology_type: str, nodes, cells):
    """How far each cell's edge nodes sit from the midpoints of their edges.

    One number per cell: the largest distance from an edge node to the
    midpoint of the edge :data:`ymf.archive.QUADRATIC_EDGES` puts it on,
    relative to that edge's length. For a straight-sided mesh it is zero up
    to rounding when the nodes are in ymf's order, and of order one when
    two edge nodes are swapped -- which makes it the check a converter runs
    on a straight-sided sample of its source. A curved mesh legitimately
    gives nonzero values.
    """
    import numpy
    if topology_type not in QUADRATIC_EDGES:
        raise YmfArchiveError(f"{topology_type!r} is not a quadratic cell type")
    edges = numpy.array(QUADRATIC_EDGES[topology_type])
    nodes = numpy.asarray(nodes, dtype=float)
    cells = numpy.asarray(cells)
    n_vertices = int(edges.max()) + 1
    a = nodes[cells[:, edges[:, 0]]]
    b = nodes[cells[:, edges[:, 1]]]
    mid = nodes[cells[:, n_vertices:n_vertices + len(edges)]]
    length = numpy.linalg.norm(b - a, axis=-1)
    offset = numpy.linalg.norm(mid - 0.5 * (a + b), axis=-1) / numpy.where(length > 0, length, 1.0)
    return offset.max(axis=1)


def check_facets(grid: Dict[str, Any], base_dir: str | Path = ".") -> Dict[str, int]:
    """Read a grid's facet list and check what structure alone cannot.

    The tags are legal for the view (:func:`ymf.archive.check_facet_tags`);
    every node number is a node of the grid; and every facet is a face (an
    edge, in 2D) of some cell, matched on its vertices. Returns
    ``{tag: count}``. ``grid`` is a canonical grid or uniform space, as
    :func:`ymf.archive.validate_approximation` returns it.
    """
    import numpy
    if "facets" not in grid:
        raise YmfArchiveError("grid has no facets")
    f = grid["facets"]
    where = "facets"
    tags = numpy.asarray(load_array(f["tags"]["data_item"], base_dir)).reshape(-1)
    check_facet_tags(tags.tolist(), f["view"], f"{where}.tags")
    facet_type = f["topology"]["type"]
    facet_nodes = numpy.asarray(load_array(f["topology"]["data_item"], base_dir))
    n_facets = f["topology"]["number_of_elements"]
    facet_nodes = facet_nodes.reshape(n_facets, -1)
    if len(tags) != n_facets:
        raise YmfArchiveError(f"{where}: {n_facets} facets but {len(tags)} tags")
    geometry = numpy.asarray(load_array(grid["geometry"]["data_item"], base_dir))
    n_nodes = geometry.shape[0]
    if facet_nodes.size and (facet_nodes.min() < 0 or facet_nodes.max() >= n_nodes):
        raise YmfArchiveError(
            f"{where}.topology: node numbers must lie in [0, {n_nodes}); got "
            f"{int(facet_nodes.min())}..{int(facet_nodes.max())}")
    cell_type = grid["topology"]["type"]
    linear = _LINEAR.get(cell_type, cell_type)
    faces = _CELL_FACETS.get(linear)
    if faces is not None and facet_type != "Mixed":
        cells = numpy.asarray(load_array(grid["topology"]["data_item"], base_dir))
        cells = cells.reshape(grid["topology"]["number_of_elements"], -1)
        known = set()
        for face in faces:
            for row in cells[:, list(face)].tolist():
                known.add(frozenset(row))
        k = _FACET_VERTICES[facet_type]
        for i, row in enumerate(facet_nodes[:, :k].tolist()):
            if frozenset(row) not in known:
                raise YmfArchiveError(
                    f"{where}: facet {i} (nodes {row}) is not a face of any "
                    f"{cell_type} cell of the grid")
    values, counts = numpy.unique(tags, return_counts=True)
    return {int(v): int(c) for v, c in zip(values, counts)}
