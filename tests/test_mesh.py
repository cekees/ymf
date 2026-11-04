"""A mesh as a ymf grid: tagged facets, provenance, and the order carried by the cell type."""

import pytest

from ymf.archive import (
    KNOWN_TOPOLOGY_TYPES,
    QUADRATIC_EDGES,
    YmfArchiveError,
    add_instant,
    attribute,
    data_item_for,
    facets,
    geometry,
    geometry_order,
    inline_approximation,
    new_approximation,
    nodes_per_cell,
    provenance,
    read_ymf,
    select_approximation,
    topology,
    uniform_grid,
    validate_approximation,
    write_ymf,
)

np = pytest.importorskip("numpy")

from ymf.mesh import (  # noqa: E402
    check_facets,
    edge_midpoint_offsets,
    edge_node_permutation,
    file_provenance,
    reorder_cells,
)


# --------------------------------------------------------------------------
# fixtures: a unit square in two triangles, its four sides tagged
# --------------------------------------------------------------------------

NODES = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
CELLS = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int32)
SIDES = np.array([[0, 1], [1, 2], [2, 3], [3, 0]], dtype=np.int32)
SIDE_TAGS = np.array([1, 2, 0, 4], dtype=np.int32)    # top side untagged


def square_space(tags=SIDE_TAGS, view="global", source=None, sides=SIDES):
    return uniform_grid(
        topology("Triangle", len(CELLS), data_item_for(CELLS, inline=True)),
        geometry(data_item_for(NODES, inline=True)),
        [attribute("elementMaterialTypes",
                   data_item_for(np.ones(2, dtype=np.int32), inline=True), center="Cell")],
        facets=facets(topology("Polyline", len(sides), data_item_for(sides, inline=True),
                               nodes_per_element=2),
                      data_item_for(np.asarray(tags, dtype=np.int32), inline=True), view=view),
        provenance=source)


def mesh_domain(space):
    domain = new_approximation()
    add_instant(domain, 0.0, {"Mesh Spatial_Domain": space})
    return domain


SOURCE = provenance("a mesh converter", version="0.1",
                    files=[{"name": "square.msh", "sha256": "ab" * 32}],
                    extra={"units": "m"})


# --------------------------------------------------------------------------
# facets and provenance survive the archive's round trips
# --------------------------------------------------------------------------


def test_a_mesh_with_facets_and_provenance_round_trips(tmp_path):
    domain = validate_approximation(mesh_domain(square_space(source=SOURCE)))
    write_ymf(domain, tmp_path / "square.ymf")
    back, _ = read_ymf(tmp_path / "square.ymf")
    assert validate_approximation(back) == domain
    space = back["times"][0]["spaces"]["Mesh Spatial_Domain"]
    assert space["provenance"] == {"tool": "a mesh converter", "version": "0.1",
                                   "files": [{"name": "square.msh", "sha256": "ab" * 32}],
                                   "extra": {"units": "m"}}
    assert space["facets"]["view"] == "global"
    assert check_facets(space) == {0: 1, 1: 1, 2: 1, 4: 1}


def test_facets_held_in_hdf5_are_read_and_checked(tmp_path):
    h5py = pytest.importorskip("h5py")
    with h5py.File(tmp_path / "square.h5", "w") as f:
        for name, array in (("nodes", NODES), ("elements", CELLS),
                            ("facets", SIDES), ("facet_tags", SIDE_TAGS)):
            f[name] = array
    ref = lambda name, a: data_item_for(a, "square.h5:/" + name)  # noqa: E731
    space = uniform_grid(
        topology("Triangle", 2, ref("elements", CELLS)), geometry(ref("nodes", NODES)),
        facets=facets(topology("Polyline", 4, ref("facets", SIDES), nodes_per_element=2),
                      ref("facet_tags", SIDE_TAGS)))
    domain = validate_approximation(mesh_domain(space))
    write_ymf(domain, tmp_path / "square.ymf")
    back, _ = read_ymf(tmp_path / "square.ymf")
    space = validate_approximation(back)["times"][0]["spaces"]["Mesh Spatial_Domain"]
    assert check_facets(space, tmp_path) == {0: 1, 1: 1, 2: 1, 4: 1}
    # inlining carries the facet arrays along with the rest
    inlined = inline_approximation(back, tmp_path)
    facet_list = inlined["times"][0]["spaces"]["Mesh Spatial_Domain"]["facets"]
    assert facet_list["tags"]["data_item"]["values"] == SIDE_TAGS.tolist()
    assert facet_list["topology"]["data_item"]["values"] == SIDES.reshape(-1).tolist()


def test_selecting_fields_keeps_the_facets():
    domain = mesh_domain(square_space())
    picked = select_approximation(domain, fields=["elementMaterialTypes"])
    assert "facets" in picked["times"][0]["spaces"]["Mesh Spatial_Domain"]


def test_ymf2xmf_writes_a_mesh_with_facets(tmp_path):
    from ymf.xdmf import write_xdmf
    write_xdmf(mesh_domain(square_space(source=SOURCE)), tmp_path / "square.xmf")
    assert "Triangle" in (tmp_path / "square.xmf").read_text()


def test_file_provenance_records_names_and_checksums(tmp_path):
    (tmp_path / "a.node").write_bytes(b"4 2 0 0\n")
    block = file_provenance("a generator", [tmp_path / "a.node"], version="1.6")
    import hashlib
    assert block == {"tool": "a generator", "version": "1.6",
                     "files": [{"name": "a.node",
                                "sha256": hashlib.sha256(b"4 2 0 0\n").hexdigest()}]}
    validate_approximation(mesh_domain(square_space(source=block)))


# --------------------------------------------------------------------------
# what the facet list refuses
# --------------------------------------------------------------------------


def test_a_global_mesh_never_stores_a_rank_interface():
    with pytest.raises(YmfArchiveError, match="-1 marks a rank interface"):
        validate_approximation(mesh_domain(square_space(tags=[1, -1, 0, 4])))


def test_a_subdomain_view_may_mark_rank_interfaces():
    space = validate_approximation(mesh_domain(
        square_space(tags=[1, -1, 0, 4], view="subdomain")))["times"][0]["spaces"]
    assert check_facets(space["Mesh Spatial_Domain"]) == {-1: 1, 0: 1, 1: 1, 4: 1}


@pytest.mark.parametrize("tags, match", [
    ([1, -2, 0, 4], "tag -2 is not a tag"),
    ([1, 2, 0], r"one tag per facet, \[4\]"),
])
def test_bad_tags_are_refused(tags, match):
    with pytest.raises(YmfArchiveError, match=match):
        validate_approximation(mesh_domain(square_space(tags=tags)))


def test_tags_are_integers():
    space = square_space()
    space["facets"]["tags"]["data_item"] = data_item_for(np.array([1.0, 2.0, 0.0, 4.0]),
                                                         inline=True)
    with pytest.raises(YmfArchiveError, match="tags are integers"):
        validate_approximation(mesh_domain(space))


def test_the_facet_type_follows_the_cell_type():
    space = square_space()
    space["facets"]["topology"]["type"] = "Triangle"
    with pytest.raises(YmfArchiveError, match="facets of Triangle cells are Polyline"):
        validate_approximation(mesh_domain(space))


def test_an_unknown_view_is_refused():
    with pytest.raises(YmfArchiveError, match="'partitioned' is not one of"):
        validate_approximation(mesh_domain(square_space(view="partitioned")))


def test_a_facet_must_be_a_face_of_a_cell():
    # (1, 3) joins two corners of the square but is an edge of neither triangle
    sides = np.array([[0, 1], [1, 3]], dtype=np.int32)
    space = validate_approximation(mesh_domain(square_space(sides=sides, tags=[1, 2])))
    with pytest.raises(YmfArchiveError, match=r"facet 1 \(nodes \[1, 3\]\) is not a face"):
        check_facets(space["times"][0]["spaces"]["Mesh Spatial_Domain"])


def test_a_facet_node_must_be_a_node_of_the_grid():
    space = validate_approximation(mesh_domain(square_space(
        sides=np.array([[0, 1], [3, 4]], dtype=np.int32), tags=[1, 2])))
    with pytest.raises(YmfArchiveError, match=r"must lie in \[0, 4\)"):
        check_facets(space["times"][0]["spaces"]["Mesh Spatial_Domain"])


@pytest.mark.parametrize("block, match", [
    ({"version": "1"}, "provenance has no 'tool'"),
    ({"tool": "t", "codes": ["x"]}, "put anything else under extra"),
    ({"tool": "t", "files": [{"name": "a"}]}, r"a file is \{name, sha256\}"),
    ({"tool": "t", "files": [{"name": "a", "sha256": "1234"}]}, "not a SHA-256 digest"),
])
def test_bad_provenance_is_refused(block, match):
    with pytest.raises(YmfArchiveError, match=match):
        validate_approximation(mesh_domain(square_space(source=block)))


# --------------------------------------------------------------------------
# geometry order is the cell type's
# --------------------------------------------------------------------------


def test_every_known_type_has_an_order():
    for t in KNOWN_TOPOLOGY_TYPES:
        order = geometry_order(t)
        assert order in (1, 2) or (t == "Mixed" and order is None), t
    assert geometry_order("Tri_6") == geometry_order("Tet_10") == 2
    assert geometry_order("Tetrahedron") == 1
    with pytest.raises(YmfArchiveError, match="not a known topology type"):
        geometry_order("Tet_11")


@pytest.mark.parametrize("cell_type, count", [
    ("Triangle", 3), ("Tri_6", 6), ("Quad_8", 8), ("Tetrahedron", 4), ("Tet_10", 10),
    ("Tetrahedron_10", 10), ("Pyramid_13", 13), ("Wedge_15", 15), ("Hex_20", 20),
    ("Polyline", None), ("Mixed", None),
])
def test_nodes_per_cell(cell_type, count):
    assert nodes_per_cell(cell_type) == count


def test_quadratic_edges_name_each_edge_once():
    for t, edges in QUADRATIC_EDGES.items():
        assert len({frozenset(e) for e in edges}) == len(edges), t


# --------------------------------------------------------------------------
# the Tet_10 order: a source that numbers edge nodes differently is permuted
# --------------------------------------------------------------------------

#: Two edge-node conventions met in practice, as test data. The first is a
#: mesh generator's: the last two edge nodes the other way round from ymf's.
#: The second is how Proteus numbers its quadratic tetrahedron's nodes
#: (FemTools: increments of 1, then 2, then 3 between vertex numbers).
GENERATOR_TET10 = ((0, 1), (1, 2), (2, 0), (0, 3), (2, 3), (1, 3))
INCREMENT_TET10 = ((0, 1), (1, 2), (2, 3), (0, 2), (1, 3), (0, 3))


def _straight_tet10(edges):
    """One straight-sided Tet_10 with its edge nodes at the midpoints of ``edges``."""
    vertices = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.0, 1, 0], [0.0, 0, 1]])
    nodes = np.vstack([vertices] + [0.5 * (vertices[a] + vertices[b]) for a, b in edges])
    return nodes, np.arange(10)[None, :]


def test_ymfs_own_order_is_the_identity():
    assert edge_node_permutation("Tet_10", QUADRATIC_EDGES["Tet_10"]) == list(range(10))
    nodes, cells = _straight_tet10(QUADRATIC_EDGES["Tet_10"])
    assert edge_midpoint_offsets("Tet_10", nodes, cells).max() < 1e-15


@pytest.mark.parametrize("source, perm", [
    (GENERATOR_TET10, [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]),
    (INCREMENT_TET10, [0, 1, 2, 3, 4, 5, 7, 9, 8, 6]),
])
def test_a_source_tet10_is_put_in_ymf_order(source, perm):
    assert edge_node_permutation("Tet_10", source) == perm
    nodes, cells = _straight_tet10(source)
    # read as ymf's order the edge nodes are off their edges...
    assert edge_midpoint_offsets("Tet_10", nodes, cells).max() > 0.3
    # ...and after the permutation every one is at its edge's midpoint
    fixed = reorder_cells(cells, "Tet_10", source)
    assert edge_midpoint_offsets("Tet_10", nodes, fixed).max() < 1e-15


def test_a_pair_written_either_way_round_is_the_same_edge():
    flipped = tuple((b, a) for a, b in QUADRATIC_EDGES["Tri_6"])
    assert edge_node_permutation("Tri_6", flipped) == list(range(6))


@pytest.mark.parametrize("cell_type, source, match", [
    ("Tetrahedron", QUADRATIC_EDGES["Tet_10"], "not a quadratic cell type"),
    ("Tet_10", QUADRATIC_EDGES["Tet_10"][:5], "has 6 distinct edges"),
    ("Tri_6", ((0, 1), (1, 2), (1, 3)), r"edge \(2, 0\) is not among"),
])
def test_a_wrong_source_convention_is_refused(cell_type, source, match):
    with pytest.raises(YmfArchiveError, match=match):
        edge_node_permutation(cell_type, source)


def test_a_tet10_mesh_carries_tri6_facets(tmp_path):
    nodes, cells = _straight_tet10(QUADRATIC_EDGES["Tet_10"])
    # the face (0, 1, 2) of the tetrahedron, with its edge nodes 4, 5, 6
    face = np.array([[0, 1, 2, 4, 5, 6]], dtype=np.int32)
    space = uniform_grid(
        topology("Tet_10", 1, data_item_for(cells.astype(np.int32), inline=True)),
        geometry(data_item_for(nodes, inline=True)),
        facets=facets(topology("Tri_6", 1, data_item_for(face, inline=True)),
                      data_item_for(np.array([3], dtype=np.int32), inline=True)))
    space = validate_approximation(mesh_domain(space))["times"][0]["spaces"]["Mesh Spatial_Domain"]
    assert check_facets(space) == {3: 1}
    space["facets"]["topology"]["type"] = "Triangle"
    with pytest.raises(YmfArchiveError, match="facets of Tet_10 cells are Tri_6 or Triangle_6"):
        validate_approximation(mesh_domain(space))
