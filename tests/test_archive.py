"""Tests for ymf.archive -- the pyyaml-only archive core."""

import sys
import warnings

import pytest

from ymf.archive import time_collections  # noqa: E402
from ymf.archive import (
    domain_arrays,
    inline_domain,
    ARCHIVE_FORMAT_VERSION,
    HAVE_LIBYAML,
    YmfArchiveError,
    add_collection,
    add_spatial_step,
    add_uniform_step,
    append_step,
    attribute,
    canonicalize_domain,
    check_dimensions,
    data_item,
    data_item_for,
    dump_grid,
    geometry,
    grid,
    load_grid,
    new_domain,
    read_ymf,
    require_libyaml,
    topology,
    validate_domain,
    write_ymf,
)


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------


def make_uniform_domain():
    """A single-grid archive -- what a serial run produces."""
    domain = new_domain("Mesh Spatial_Domain")
    for t in (0.0, 1.0):
        add_uniform_step(
            domain,
            t,
            topology(
                "Tetrahedron",
                687737,
                data_item([687737, 4], f"poisson.h5:/elements_t{int(t)}", data_type="Int"),
            ),
            geometry(data_item([118084, 3], f"poisson.h5:/nodes_t{int(t)}", precision=8)),
            [
                attribute("u", data_item([118084], f"poisson.h5:/u_t{int(t)}", precision=8)),
                attribute(
                    "velocity",
                    data_item([118084, 3], f"poisson.h5:/vel_t{int(t)}", precision=8),
                    attribute_type="Vector",
                ),
            ],
        )
    return domain


def make_spatial_domain(n_ranks=3):
    """A subdomain-per-rank archive -- what a parallel run produces."""
    domain = new_domain("Mesh Spatial_Domain")
    for t in (0.0, 0.5):
        grids = []
        for r in range(n_ranks):
            grids.append(
                grid(
                    topology(
                        "Triangle",
                        1000 + r,
                        data_item([1000 + r, 3], f"out.h5:/elements_p{r}_t{t}", data_type="Int"),
                    ),
                    geometry(data_item([600 + r, 2], f"out.h5:/nodes_p{r}_t{t}"), geometry_type="XY"),
                    [attribute("p", data_item([600 + r], f"out.h5:/p_p{r}_t{t}", precision=8))],
                    name=f"subdomain_{r}",
                )
            )
        add_spatial_step(domain, t, grids)
    return domain


# --------------------------------------------------------------------------
# constructors
# --------------------------------------------------------------------------


def test_data_item_infers_hdf_format_from_data():
    item = data_item([10], "f.h5:/x")
    assert item["format"] == "HDF"
    assert item["data"] == "f.h5:/x"
    assert "include" not in item


def test_data_item_infers_xml_format_from_include():
    item = data_item([10], include="./data/x.txt")
    assert item["format"] == "XML"
    assert item["include"] == "./data/x.txt"
    assert "data" not in item


def test_data_item_precision_defaults_to_xdmf_default_of_four():
    # XDMF's own documented default. Faithful to the format, and usually
    # *not* what a float64 array wants -- hence data_item_for().
    assert data_item([10], "f.h5:/x")["precision"] == 4


@pytest.mark.parametrize(
    "kwargs",
    [
        {},  # neither data nor include
        {"data": "f.h5:/x", "include": "./x.txt"},  # both
    ],
)
def test_data_item_requires_exactly_one_of_data_or_include(kwargs):
    with pytest.raises(YmfArchiveError, match="exactly one of"):
        data_item([10], **kwargs)


def test_topology_omits_nodes_per_element_when_not_given():
    # There is no correct value to invent for a Tetrahedron, so the key
    # stays absent rather than being defaulted to something wrong.
    assert "nodes_per_element" not in topology("Tetrahedron", 5, data_item([5, 4], "f.h5:/e"))


def test_topology_keeps_nodes_per_element_when_given():
    topo = topology("Polyline", 5, data_item([10], "f.h5:/e"), nodes_per_element=2)
    assert topo["nodes_per_element"] == 2


# --------------------------------------------------------------------------
# data_item_for -- dtype inference
# --------------------------------------------------------------------------


def test_data_item_for_infers_dimensions_datatype_and_precision():
    numpy = pytest.importorskip("numpy")
    arr = numpy.zeros((118084, 3), dtype="float64")
    item = data_item_for(arr, "out.h5:/nodes")
    assert item["dimensions"] == [118084, 3]
    assert item["data_type"] == "Float"
    assert item["precision"] == 8  # read off the dtype, not defaulted to 4


def test_data_item_for_distinguishes_float32_from_float64():
    numpy = pytest.importorskip("numpy")
    assert data_item_for(numpy.zeros(4, dtype="float32"), "f.h5:/x")["precision"] == 4
    assert data_item_for(numpy.zeros(4, dtype="float64"), "f.h5:/x")["precision"] == 8


def test_data_item_for_maps_integer_kinds():
    numpy = pytest.importorskip("numpy")
    assert data_item_for(numpy.zeros(4, dtype="int32"), "f.h5:/x")["data_type"] == "Int"
    assert data_item_for(numpy.zeros(4, dtype="uint8"), "f.h5:/x")["data_type"] == "UInt"


def test_data_item_for_allows_deliberate_flattening():
    # Declaring [N*k] for an (N, k) connectivity array is normal for XDMF.
    numpy = pytest.importorskip("numpy")
    arr = numpy.zeros((100, 4), dtype="int32")
    item = data_item_for(arr, "f.h5:/e", dimensions=[400])
    assert item["dimensions"] == [400]


def test_data_item_for_rejects_a_genuine_shape_mismatch():
    numpy = pytest.importorskip("numpy")
    arr = numpy.zeros((100, 4), dtype="int32")
    with pytest.raises(YmfArchiveError, match="399 values"):
        data_item_for(arr, "f.h5:/e", dimensions=[399])


def test_data_item_for_rejects_a_dtype_with_no_xdmf_equivalent():
    numpy = pytest.importorskip("numpy")
    with pytest.raises(YmfArchiveError, match="no XDMF"):
        data_item_for(numpy.zeros(3, dtype="complex128"), "f.h5:/x")


# --------------------------------------------------------------------------
# check_dimensions
# --------------------------------------------------------------------------


def test_check_dimensions_accepts_an_exact_match():
    check_dimensions(data_item([100, 3], "f.h5:/x"), (100, 3))


def test_check_dimensions_accepts_a_flattened_declaration():
    check_dimensions(data_item([300], "f.h5:/x"), (100, 3))


def test_check_dimensions_rejects_a_mismatch_and_names_the_field():
    with pytest.raises(YmfArchiveError, match=r"u\.data_item\.dimensions"):
        check_dimensions(data_item([99], "f.h5:/x"), (100,), where="u.data_item")


# --------------------------------------------------------------------------
# canonicalization
# --------------------------------------------------------------------------


def test_canonicalize_is_idempotent():
    once = canonicalize_domain(make_uniform_domain())
    assert canonicalize_domain(once) == once


def test_canonicalize_is_idempotent_for_spatial_domains():
    once = canonicalize_domain(make_spatial_domain())
    assert canonicalize_domain(once) == once


def test_canonicalize_fills_every_documented_default():
    sparse = {
        "time_collection": {
            "data": [
                {
                    "time": 0.0,
                    "topology": {
                        "type": "Triangle",
                        "number_of_elements": 2,
                        "data_item": {"dimensions": [2, 3], "data": "f.h5:/e"},
                    },
                    "geometry": {"data_item": {"dimensions": [4, 2], "data": "f.h5:/n"}},
                    "attributes": [{"name": "u", "data_item": {"dimensions": [4], "data": "f.h5:/u"}}],
                }
            ]
        }
    }
    c = canonicalize_domain(sparse)
    # the singular format 1 input form is upgraded to format 2
    assert list(c) == ["times"]
    assert len(time_collections(c)) == 1
    tc = time_collections(c)[0]
    assert tc["name"] == "TimeCollection"
    step = tc["data"][0]
    assert step["geometry"]["type"] == "XYZ"
    assert step["topology"]["data_item"]["format"] == "HDF"
    assert step["topology"]["data_item"]["data_type"] == "Float"
    assert step["topology"]["data_item"]["precision"] == 4
    attr = step["attributes"][0]
    assert attr["attribute_type"] == "Scalar"
    assert attr["center"] == "Node"


def test_canonicalize_supplies_an_empty_attributes_list():
    domain = new_domain()
    add_uniform_step(
        domain, 0.0,
        topology("Triangle", 1, data_item([1, 3], "f.h5:/e")),
        geometry(data_item([3, 2], "f.h5:/n")),
    )
    assert time_collections(canonicalize_domain(domain))[0]["data"][0]["attributes"] == []


def test_canonicalize_of_an_empty_domain_is_empty():
    assert canonicalize_domain({}) == {}


def test_step_cannot_be_both_uniform_and_spatial():
    domain = {
        "time_collection": {
            "data": [
                {
                    "time": 0.0,
                    "spatial_collection": [],
                    "topology": {"type": "Triangle", "number_of_elements": 1,
                                 "data_item": {"dimensions": [1], "data": "f.h5:/e"}},
                }
            ]
        }
    }
    with pytest.raises(YmfArchiveError, match="not both"):
        canonicalize_domain(domain)


def test_step_must_be_either_uniform_or_spatial():
    with pytest.raises(YmfArchiveError, match="neither"):
        canonicalize_domain({"time_collection": {"data": [{"time": 0.0}]}})


@pytest.mark.parametrize("missing", ["type", "number_of_elements", "data_item"])
def test_canonicalize_reports_a_missing_topology_key_by_path(missing):
    topo = topology("Triangle", 1, data_item([1, 3], "f.h5:/e"))
    del topo[missing]
    domain = {
        "time_collection": {
            "data": [{"time": 0.0, "topology": topo,
                      "geometry": geometry(data_item([3, 2], "f.h5:/n"))}]
        }
    }
    with pytest.raises(YmfArchiveError, match=rf"time_collections\[0\]\.data\[0\]\.topology.*{missing}"):
        canonicalize_domain(domain)


# --------------------------------------------------------------------------
# validation -- the cheap half of the trust boundary
# --------------------------------------------------------------------------


def test_validate_accepts_the_uniform_and_spatial_fixtures():
    validate_domain(make_uniform_domain())
    validate_domain(make_spatial_domain())


def test_validate_returns_the_canonicalized_document():
    assert validate_domain(make_uniform_domain()) == canonicalize_domain(make_uniform_domain())


def test_validate_rejects_an_unknown_topology_type():
    domain = make_uniform_domain()
    time_collections(domain)[0]["data"][0]["topology"]["type"] = "Dodecahedron"
    with pytest.raises(YmfArchiveError, match="not a known XDMF topology type"):
        validate_domain(domain)


def test_validate_accepts_every_topology_type_proteus_emits():
    # The union of the Archiver.py writer set and the MeshTools.py
    # topologyid2name reader map. If any of these ever stops validating,
    # Proteus loses the ability to write a mesh it can currently write.
    from ymf.archive import KNOWN_TOPOLOGY_TYPES

    proteus_types = {
        "Polyvertex", "Polyline", "Triangle", "Quadrilateral", "Tetrahedron",
        "Hexahedron", "Edge_3", "Tri_6", "Tet_10", "Wedge", "Mixed",
    }
    assert proteus_types <= KNOWN_TOPOLOGY_TYPES


def test_validate_rejects_a_bad_centering():
    domain = make_uniform_domain()
    time_collections(domain)[0]["data"][0]["attributes"][0]["center"] = "Vertex"
    with pytest.raises(YmfArchiveError, match=r"attributes\[0\]\.center"):
        validate_domain(domain)


def test_validate_rejects_a_bad_attribute_type():
    domain = make_uniform_domain()
    time_collections(domain)[0]["data"][0]["attributes"][0]["attribute_type"] = "Scalarr"
    with pytest.raises(YmfArchiveError, match=r"attributes\[0\]\.attribute_type"):
        validate_domain(domain)


def test_validate_rejects_a_bad_data_type():
    domain = make_uniform_domain()
    time_collections(domain)[0]["data"][0]["attributes"][0]["data_item"]["data_type"] = "Double"
    with pytest.raises(YmfArchiveError, match="data_type"):
        validate_domain(domain)


def test_validate_rejects_duplicate_attribute_names_within_one_grid():
    # Viewers key fields by name within a grid, so a duplicate silently
    # shadows rather than erroring at read time -- worth catching on write.
    domain = make_uniform_domain()
    attrs = time_collections(domain)[0]["data"][0]["attributes"]
    attrs[1]["name"] = attrs[0]["name"]
    attrs[1]["center"] = attrs[0]["center"]
    with pytest.raises(YmfArchiveError, match="duplicate .*-centred attribute name"):
        validate_domain(domain)


def test_validate_allows_the_same_name_on_different_centerings():
    # Point data and cell data are separate namespaces in a viewer; an
    # exporter writing T both cell- and node-centred is legitimate.
    domain = make_uniform_domain()
    attrs = time_collections(domain)[0]["data"][0]["attributes"]
    attrs[1]["name"] = attrs[0]["name"]
    attrs[0]["center"], attrs[1]["center"] = "Cell", "Node"
    validate_domain(domain)


def test_validate_allows_the_same_attribute_name_in_different_subdomains():
    # Each rank writes its own 'p'; that is correct, not a duplicate.
    validate_domain(make_spatial_domain(n_ranks=4))


def test_validate_rejects_an_empty_spatial_collection():
    domain = new_domain()
    add_spatial_step(domain, 0.0, [])
    with pytest.raises(YmfArchiveError, match="must not be empty"):
        validate_domain(domain)


def test_validate_rejects_empty_dimensions():
    domain = make_uniform_domain()
    time_collections(domain)[0]["data"][0]["attributes"][0]["data_item"]["dimensions"] = []
    with pytest.raises(YmfArchiveError, match="must not be empty"):
        validate_domain(domain)


def test_validate_rejects_a_negative_dimension():
    domain = make_uniform_domain()
    time_collections(domain)[0]["data"][0]["attributes"][0]["data_item"]["dimensions"] = [-1]
    with pytest.raises(YmfArchiveError, match="non-negative"):
        validate_domain(domain)


def test_validate_of_an_empty_domain_is_empty():
    assert validate_domain({}) == {}


# --------------------------------------------------------------------------
# YAML round-trip
# --------------------------------------------------------------------------


def test_uniform_domain_round_trips_through_yaml(tmp_path):
    domain = make_uniform_domain()
    path = tmp_path / "archive.ymf"
    write_ymf(domain, path)
    read_domain, extra = read_ymf(path)
    assert read_domain == canonicalize_domain(domain)
    assert extra is None


def test_spatial_domain_round_trips_through_yaml(tmp_path):
    domain = make_spatial_domain()
    path = tmp_path / "archive.ymf"
    write_ymf(domain, path)
    read_domain, _ = read_ymf(path)
    assert read_domain == canonicalize_domain(domain)


def test_text_include_data_items_round_trip_through_yaml(tmp_path):
    domain = new_domain()
    add_uniform_step(
        domain, 0.0,
        topology("Triangle", 2, data_item([2, 3], include="./d/elements.txt", data_type="Int")),
        geometry(data_item([4, 2], include="./d/nodes.txt")),
        [attribute("u", data_item([4], include="./d/u.txt"))],
    )
    path = tmp_path / "archive.ymf"
    write_ymf(domain, path)
    read_domain, _ = read_ymf(path)
    assert read_domain == canonicalize_domain(domain)
    assert time_collections(read_domain)[0]["data"][0]["attributes"][0]["data_item"]["include"] == "./d/u.txt"


def test_extra_payload_round_trips_and_is_stored_readably(tmp_path):
    extra = {"model": {"name": "poisson", "provenance": "proteus_derived"}}
    path = tmp_path / "archive.ymf"
    write_ymf(make_uniform_domain(), path, extra=extra)
    _, read_extra = read_ymf(path)
    assert read_extra == extra
    # Unlike the XDMF path, which has to base64 this into an <Information>
    # attribute, the YAML archive holds it in the clear.
    assert "poisson" in path.read_text()


def test_unicode_survives_the_yaml_round_trip(tmp_path):
    extra = {"strong_form_expression": "-∇·(κ∇u) = f on Ω"}
    path = tmp_path / "archive.ymf"
    write_ymf(make_uniform_domain(), path, extra=extra)
    assert read_ymf(path)[1] == extra


def test_write_ymf_stamps_the_format_version(tmp_path):
    path = tmp_path / "archive.ymf"
    write_ymf(make_uniform_domain(), path)
    assert f"ymf_archive_version: {ARCHIVE_FORMAT_VERSION}" in path.read_text()


def test_read_ymf_rejects_a_future_format_version(tmp_path):
    path = tmp_path / "archive.ymf"
    write_ymf(make_uniform_domain(), path)
    path.write_text(
        path.read_text().replace(
            f"ymf_archive_version: {ARCHIVE_FORMAT_VERSION}",
            f"ymf_archive_version: {ARCHIVE_FORMAT_VERSION + 1}",
        )
    )
    with pytest.raises(YmfArchiveError, match="archive format version"):
        read_ymf(path)


def test_read_ymf_rejects_a_non_mapping_document(tmp_path):
    path = tmp_path / "archive.ymf"
    path.write_text("- just\n- a\n- list\n")
    with pytest.raises(YmfArchiveError, match="mapping"):
        read_ymf(path)


def test_write_ymf_validates_by_default(tmp_path):
    domain = make_uniform_domain()
    time_collections(domain)[0]["data"][0]["attributes"][0]["center"] = "Nowhere"
    with pytest.raises(YmfArchiveError, match="center"):
        write_ymf(domain, tmp_path / "archive.ymf")


def test_write_ymf_can_skip_validation_deliberately(tmp_path):
    domain = make_uniform_domain()
    time_collections(domain)[0]["data"][0]["attributes"][0]["center"] = "Nowhere"
    path = tmp_path / "archive.ymf"
    write_ymf(domain, path, validate=False)
    assert time_collections(read_ymf(path)[0])[0]["data"][0]["attributes"][0]["center"] == "Nowhere"


# --------------------------------------------------------------------------
# grid fragments -- the MPI / in-HDF5 metadata path
# --------------------------------------------------------------------------


def test_a_single_grid_fragment_round_trips_as_a_string():
    g = time_collections(make_spatial_domain())[0]["data"][0]["spatial_collection"][0]
    assert load_grid(dump_grid(g)) == g


def test_grid_fragments_reassemble_into_a_spatial_step():
    # This is the shape of the parallel write: each rank dumps its own grid
    # to a string, the strings are collected, and the master reassembles.
    original = make_spatial_domain(n_ranks=4)
    fragments = [
        dump_grid(g) for g in time_collections(original)[0]["data"][0]["spatial_collection"]
    ]
    rebuilt = new_domain(time_collections(original)[0]["name"])
    add_spatial_step(rebuilt, 0.0, [load_grid(f) for f in fragments])
    assert (
        time_collections(canonicalize_domain(rebuilt))[0]["data"][0]
        == time_collections(canonicalize_domain(original))[0]["data"][0]
    )


# --------------------------------------------------------------------------
# libyaml
# --------------------------------------------------------------------------


@pytest.mark.skipif(not HAVE_LIBYAML, reason="pyyaml built without libyaml")
def test_require_libyaml_is_quiet_when_the_c_extension_is_present():
    require_libyaml()


@pytest.mark.skipif(HAVE_LIBYAML, reason="pyyaml built with libyaml")
def test_require_libyaml_raises_when_the_c_extension_is_missing():  # pragma: no cover
    with pytest.raises(YmfArchiveError, match="libyaml"):
        require_libyaml()


# --------------------------------------------------------------------------
# multiple time collections -- one per finite-element space
# --------------------------------------------------------------------------


def _two_space_domain():
    """The shape a real Proteus archive has: linear mesh + quadratic space."""
    domain = new_domain("Mesh Spatial_Domain")
    add_collection(domain, "Mesh_c0p2_Lagrange")
    for t in (0.0, 1.0):
        add_uniform_step(
            domain, t,
            topology("Triangle", 128, data_item([128, 3], f"o.h5:/e_t{int(t)}", data_type="Int")),
            geometry(data_item([81, 3], f"o.h5:/n_t{int(t)}", precision=8)),
            [attribute("u", data_item([81], f"o.h5:/u_t{int(t)}", precision=8))],
        )
        add_uniform_step(
            domain, t,
            topology("Tri_6", 128, data_item([128, 6], f"o.h5:/e6_t{int(t)}", data_type="Int")),
            geometry(data_item([289, 3], f"o.h5:/n6_t{int(t)}", precision=8)),
            [attribute("u", data_item([289], f"o.h5:/u6_t{int(t)}", precision=8))],
            collection="Mesh_c0p2_Lagrange",
        )
    return domain


def test_new_domain_starts_with_one_named_collection():
    domain = new_domain("Mesh Spatial_Domain")
    assert [c["name"] for c in time_collections(domain)] == ["Mesh Spatial_Domain"]


def test_add_collection_appends_a_second_mesh():
    domain = _two_space_domain()
    names = [c["name"] for c in time_collections(domain)]
    assert names == ["Mesh Spatial_Domain", "Mesh_c0p2_Lagrange"]
    assert [len(c["data"]) for c in time_collections(domain)] == [2, 2]


def test_add_collection_rejects_a_duplicate_name():
    domain = new_domain("Mesh Spatial_Domain")
    with pytest.raises(YmfArchiveError, match="already has a time collection"):
        add_collection(domain, "Mesh Spatial_Domain")


def test_steps_can_target_a_collection_by_name_or_index():
    domain = _two_space_domain()
    by_name = time_collections(domain)[1]["data"]
    assert all(step["topology"]["type"] == "Tri_6" for step in by_name)
    # index 0 is the default target
    assert all(
        step["topology"]["type"] == "Triangle"
        for step in time_collections(domain)[0]["data"]
    )


def test_targeting_an_unknown_collection_name_lists_the_real_ones():
    domain = _two_space_domain()
    with pytest.raises(YmfArchiveError, match=r"Mesh Spatial_Domain.*Mesh_c0p2_Lagrange"):
        add_uniform_step(
            domain, 0.0,
            topology("Triangle", 1, data_item([1, 3], "o.h5:/e")),
            geometry(data_item([3, 3], "o.h5:/n")),
            collection="Mesh_nonexistent",
        )


def test_targeting_an_out_of_range_collection_index_reports_the_count():
    domain = _two_space_domain()
    with pytest.raises(YmfArchiveError, match="2 time collection"):
        add_uniform_step(
            domain, 0.0,
            topology("Triangle", 1, data_item([1, 3], "o.h5:/e")),
            geometry(data_item([3, 3], "o.h5:/n")),
            collection=7,
        )


def test_a_two_space_domain_validates_and_round_trips(tmp_path):
    domain = _two_space_domain()
    validate_domain(domain)
    path = tmp_path / "archive.ymf"
    write_ymf(domain, path)
    assert read_ymf(path)[0] == canonicalize_domain(domain)


def test_the_same_attribute_name_may_appear_in_different_collections():
    # 'u' on the linear mesh and 'u' on the quadratic space are different
    # fields on different meshes, not a duplicate.
    validate_domain(_two_space_domain())


def test_duplicate_collection_names_are_rejected_at_canonicalization():
    domain = {
        "time_collections": [
            {"name": "Mesh", "data": []},
            {"name": "Mesh", "data": []},
        ]
    }
    with pytest.raises(YmfArchiveError, match="duplicate collection name"):
        canonicalize_domain(domain)


def test_singular_time_collection_is_accepted_as_input_sugar():
    singular = {"time_collection": {"name": "Mesh", "data": []}}
    plural = {"time_collections": [{"name": "Mesh", "data": []}]}
    assert canonicalize_domain(singular) == canonicalize_domain(plural)


# ---------------------------------------------------------------------------
# inline values: a self-contained archive, no HDF5 or sidecar
# ---------------------------------------------------------------------------


def test_data_item_takes_values_inline_and_flattens_them_row_major():
    item = data_item([2, 3], values=[[0, 1, 3], [0, 3, 2]], data_type="Int")
    assert item == {"format": "XML", "data_type": "Int", "precision": 4,
                    "dimensions": [2, 3], "values": [0, 1, 3, 0, 3, 2]}


def test_data_item_needs_exactly_one_of_data_include_and_values():
    with pytest.raises(YmfArchiveError, match="exactly one"):
        data_item([1], "x.h5:/a", values=[1.0])


def test_data_item_for_inline_copies_an_array_in():
    np = pytest.importorskip("numpy")
    item = data_item_for(np.arange(6, dtype=np.int32).reshape(2, 3), inline=True)
    assert item["values"] == [0, 1, 2, 3, 4, 5]
    assert (item["format"], item["data_type"], item["precision"]) == ("XML", "Int", 4)
    assert all(type(v) is int for v in item["values"])


def test_data_item_for_inline_takes_no_reference():
    np = pytest.importorskip("numpy")
    with pytest.raises(YmfArchiveError, match="no data= or include="):
        data_item_for(np.zeros(2), "x.h5:/a", inline=True)


@pytest.mark.parametrize("values, data_type, message", [
    ([1, 2], "Float", r"declared \[3\] \(3 values\) but the array has shape \[2\]"),
    ([1, 2, 2.5], "Int", r"values\[2\]: 2.5 is not an integer"),
    ([1, "2", 3], "Float", r"values\[1\]: '2' is not a number"),
    ([1, True, 3], "Int", r"values\[1\]: True is not a number"),
])
def test_inline_values_are_checked(values, data_type, message):
    with pytest.raises(YmfArchiveError, match=message):
        data_item([3], values=values, data_type=data_type)


def test_validate_refuses_inline_values_under_another_format():
    domain = new_domain()
    add_uniform_step(domain, 0.0,
                     topology("Polyvertex", 1, data_item([1], values=[0], data_type="Int"),
                              nodes_per_element=1),
                     geometry(data_item([1, 3], values=[0.0, 0.0, 0.0])))
    time_collections(domain)[0]["data"][0]["geometry"]["data_item"]["format"] = "HDF"
    with pytest.raises(YmfArchiveError, match=r"geometry\.data_item: inline 'values' need format='XML'"):
        validate_domain(domain)


def test_an_inline_archive_round_trips_and_writes_each_array_on_few_lines(tmp_path):
    domain = new_domain("mesh")
    add_uniform_step(
        domain, 0.0,
        topology("Triangle", 2, data_item([2, 3], values=[0, 1, 3, 0, 3, 2], data_type="Int")),
        geometry(data_item([4, 3], values=[0.0] * 12, precision=8)),
        [attribute("u", data_item([4], values=[0.1, 1 / 3, 2.0, -1e300], precision=8))])
    path = tmp_path / "self_contained.ymf"
    write_ymf(domain, path)
    text = path.read_text()
    assert "values: [0, 1, 3, 0, 3, 2]" in text
    read_back, _ = read_ymf(path)
    assert read_back == validate_domain(domain)


# ---------------------------------------------------------------------------
# arrays: loading, inlining, comparing across storage
# ---------------------------------------------------------------------------


def test_dimensions_are_written_on_one_line(tmp_path):
    domain = new_domain("mesh")
    add_uniform_step(domain, 0.0,
                     topology("Triangle", 8, data_item([8, 3], "m.h5:/e", data_type="Int")),
                     geometry(data_item([9, 3], "m.h5:/n", precision=8)))
    write_ymf(domain, tmp_path / "a.ymf")
    text = (tmp_path / "a.ymf").read_text()
    assert "dimensions: [8, 3]" in text and "dimensions: [9, 3]" in text


def test_inline_domain_copies_hdf5_arrays_in_exactly(tmp_path):
    np = pytest.importorskip("numpy")
    h5py = pytest.importorskip("h5py")
    cells = np.array([[0, 1, 3], [0, 3, 2]], dtype=np.int32)
    nodes = np.random.default_rng(1).random((4, 3))
    u = np.linspace(0, 1, 4, dtype=np.float32)
    with h5py.File(tmp_path / "m.h5", "w") as f:
        f["e"], f["n"], f["u"] = cells, nodes, u
    domain = new_domain("mesh")
    add_uniform_step(domain, 0.0, topology("Triangle", 2, data_item_for(cells, "m.h5:/e")),
                     geometry(data_item_for(nodes, "m.h5:/n")),
                     [attribute("u", data_item_for(u, "m.h5:/u"))])
    referenced = domain_arrays(domain, tmp_path)
    inline = inline_domain(domain, tmp_path)
    assert all("values" in a["data_item"] for a in time_collections(inline)[0]["data"][0]["attributes"])
    copied = domain_arrays(inline)            # no base directory needed any more
    assert set(copied) == set(referenced)
    for key in referenced:
        assert copied[key].dtype == referenced[key].dtype
        assert np.array_equal(copied[key], referenced[key])


# --------------------------------------------------------------------------
# append_step: a writer should not rewrite the document to record a step
# --------------------------------------------------------------------------


def _one_step(time):
    """A step dict, built the way a writer builds one."""
    domain = new_domain("scratch")
    add_uniform_step(
        domain, time,
        {"type": "Triangle", "number_of_elements": 4,
         "data_item": data_item([4, 3], "o.h5:/e", data_type="Int")},
        {"type": "XYZ", "data_item": data_item([6, 3], "o.h5:/n")},
        [], collection=0)
    return time_collections(domain)[0]["data"][0]


def _archive(tmp_path, name="a.ymf"):
    path = tmp_path / name
    write_ymf(new_domain("Mesh Spatial_Domain"), path)
    return path


def _times(domain):
    return {c["name"]: [s["time"] for s in c["data"]]
            for c in time_collections(domain)}


def test_appended_steps_are_read_back(tmp_path):
    path = _archive(tmp_path)
    for t in (0.0, 0.5, 1.0):
        append_step(path, "Mesh Spatial_Domain", _one_step(t))
    domain, _ = read_ymf(path)
    assert _times(domain) == {"Mesh Spatial_Domain": [0.0, 0.5, 1.0]}


def test_a_collection_may_appear_after_the_header(tmp_path):
    """Collections are created lazily, so the header cannot know them all.

    A model that first archives at step 40 adds its collection then. That is
    why a step record names its collection instead of indexing into the
    header's list -- an index agreed at step 0 would be wrong by step 40.
    """
    path = _archive(tmp_path)
    append_step(path, "Mesh Spatial_Domain", _one_step(0.0))
    append_step(path, "Mesh_c0p2_Lagrange", _one_step(0.0))   # new here
    append_step(path, "Mesh Spatial_Domain", _one_step(0.5))
    domain, _ = read_ymf(path)
    assert _times(domain) == {"Mesh Spatial_Domain": [0.0, 0.5],
                              "Mesh_c0p2_Lagrange": [0.0]}


def test_appending_matches_writing_the_whole_document(tmp_path):
    """The two ways to produce an archive have to agree.

    Otherwise what a reader gets depends on how the writer happened to
    build the file, which is the bug this is meant to prevent rather than
    introduce.
    """
    appended = _archive(tmp_path, "appended.ymf")
    whole = new_domain("Mesh Spatial_Domain")
    for t in (0.0, 0.5):
        step = _one_step(t)
        append_step(appended, "Mesh Spatial_Domain", step)
        time_collections(whole)[0]["data"].append(step)
    written = tmp_path / "written.ymf"
    write_ymf(whole, written)
    assert read_ymf(appended)[0] == read_ymf(written)[0]


def test_only_the_step_being_written_is_lost_to_a_kill(tmp_path):
    """The reason for appending at all.

    Note what makes this detectable: a truncated YAML document usually still
    *parses* -- cutting inside a mapping leaves a shorter mapping -- so the
    reader cannot judge completeness by parsing. It judges by the explicit
    document-end marker, which is written last.
    """
    path = _archive(tmp_path)
    for t in (0.0, 0.5, 1.0):
        append_step(path, "Mesh Spatial_Domain", _one_step(t))
    with open(path, "a", encoding="utf-8") as stream:
        stream.write("---\ncollection: Mesh Spatial_Domain\n"
                     "step:\n  time: 1.5\n  topology:\n    type: Tri")
    with pytest.warns(RuntimeWarning, match="did not finish"):
        domain, _ = read_ymf(path)
    assert _times(domain) == {"Mesh Spatial_Domain": [0.0, 0.5, 1.0]}


@pytest.mark.parametrize("tail", ["---\n", "---\ncollection: X\n", "---"])
def test_a_partial_record_is_dropped_wherever_it_was_cut(tmp_path, tail):
    path = _archive(tmp_path)
    append_step(path, "Mesh Spatial_Domain", _one_step(0.0))
    with open(path, "a", encoding="utf-8") as stream:
        stream.write(tail)
    with pytest.warns(RuntimeWarning):
        domain, _ = read_ymf(path)
    assert _times(domain) == {"Mesh Spatial_Domain": [0.0]}


def test_a_kill_during_the_very_first_step_leaves_the_header(tmp_path):
    """No end marker exists yet, so the whole appended region is partial."""
    path = _archive(tmp_path)
    with open(path, "a", encoding="utf-8") as stream:
        stream.write("---\ncollection: X\nstep:\n  time: 0.0\n  topo")
    with pytest.warns(RuntimeWarning):
        domain, _ = read_ymf(path)
    # a space exists once an instant holds it; the header held none
    assert _times(domain) == {}


def test_a_complete_archive_does_not_warn(tmp_path):
    path = _archive(tmp_path)
    append_step(path, "Mesh Spatial_Domain", _one_step(0.0))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        read_ymf(path)


def test_corruption_before_the_end_is_an_error_not_a_truncation(tmp_path):
    """A bad document in the middle is not a killed writer, and must not be
    silently treated as one -- that would hide real damage."""
    path = _archive(tmp_path)
    for t in (0.0, 0.5):
        append_step(path, "Mesh Spatial_Domain", _one_step(t))
    text = path.read_text(encoding="utf-8")
    path.write_text(text.replace("time: 0.0", "time: [unclosed", 1),
                    encoding="utf-8")
    with pytest.raises(Exception):
        read_ymf(path)


def test_appending_is_constant_work_per_step(tmp_path):
    """What this exists for: the file grows, the write does not.

    Rewriting is O(steps) per step and so O(steps^2) over a run. Appending
    writes one record whatever the archive already holds -- asserted on
    bytes written, which is deterministic, rather than on time.
    """
    path = _archive(tmp_path)
    sizes = []
    for t in range(12):
        before = path.stat().st_size
        append_step(path, "Mesh Spatial_Domain", _one_step(float(t)))
        sizes.append(path.stat().st_size - before)
    # Not byte-identical: the time value's own text is a character wider at
    # 10.0 than at 0.0. What matters is that the write does not grow with
    # what the archive already holds, so bound the spread rather than
    # demanding equality -- a rewrite would show the last write an order of
    # magnitude above the first.
    assert max(sizes) - min(sizes) <= 8, sizes
    assert path.stat().st_size > 10 * max(sizes)   # the file grew; the write did not
