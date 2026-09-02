"""Tests for ymf.archive -- the pyyaml-only archive core."""

import sys

import pytest

from ymf.archive import (
    ARCHIVE_FORMAT_VERSION,
    HAVE_LIBYAML,
    YmfArchiveError,
    add_collection,
    add_spatial_step,
    add_uniform_step,
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
    assert item["Format"] == "HDF"
    assert item["Data"] == "f.h5:/x"
    assert "Include" not in item


def test_data_item_infers_xml_format_from_include():
    item = data_item([10], include="./data/x.txt")
    assert item["Format"] == "XML"
    assert item["Include"] == "./data/x.txt"
    assert "Data" not in item


def test_data_item_precision_defaults_to_xdmf_default_of_four():
    # XDMF's own documented default. Faithful to the format, and usually
    # *not* what a float64 array wants -- hence data_item_for().
    assert data_item([10], "f.h5:/x")["Precision"] == 4


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
    assert "NodesPerElement" not in topology("Tetrahedron", 5, data_item([5, 4], "f.h5:/e"))


def test_topology_keeps_nodes_per_element_when_given():
    topo = topology("Polyline", 5, data_item([10], "f.h5:/e"), nodes_per_element=2)
    assert topo["NodesPerElement"] == 2


# --------------------------------------------------------------------------
# data_item_for -- dtype inference
# --------------------------------------------------------------------------


def test_data_item_for_infers_dimensions_datatype_and_precision():
    numpy = pytest.importorskip("numpy")
    arr = numpy.zeros((118084, 3), dtype="float64")
    item = data_item_for(arr, "out.h5:/nodes")
    assert item["Dimensions"] == [118084, 3]
    assert item["DataType"] == "Float"
    assert item["Precision"] == 8  # read off the dtype, not defaulted to 4


def test_data_item_for_distinguishes_float32_from_float64():
    numpy = pytest.importorskip("numpy")
    assert data_item_for(numpy.zeros(4, dtype="float32"), "f.h5:/x")["Precision"] == 4
    assert data_item_for(numpy.zeros(4, dtype="float64"), "f.h5:/x")["Precision"] == 8


def test_data_item_for_maps_integer_kinds():
    numpy = pytest.importorskip("numpy")
    assert data_item_for(numpy.zeros(4, dtype="int32"), "f.h5:/x")["DataType"] == "Int"
    assert data_item_for(numpy.zeros(4, dtype="uint8"), "f.h5:/x")["DataType"] == "UInt"


def test_data_item_for_allows_deliberate_flattening():
    # Declaring [N*k] for an (N, k) connectivity array is normal for XDMF.
    numpy = pytest.importorskip("numpy")
    arr = numpy.zeros((100, 4), dtype="int32")
    item = data_item_for(arr, "f.h5:/e", dimensions=[400])
    assert item["Dimensions"] == [400]


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
    with pytest.raises(YmfArchiveError, match=r"u\.DataItem\.Dimensions"):
        check_dimensions(data_item([99], "f.h5:/x"), (100,), where="u.DataItem")


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
        "TimeCollection": {
            "Data": [
                {
                    "Time": 0.0,
                    "Topology": {
                        "Type": "Triangle",
                        "NumberOfElements": 2,
                        "DataItem": {"Dimensions": [2, 3], "Data": "f.h5:/e"},
                    },
                    "Geometry": {"DataItem": {"Dimensions": [4, 2], "Data": "f.h5:/n"}},
                    "Attributes": [{"Name": "u", "DataItem": {"Dimensions": [4], "Data": "f.h5:/u"}}],
                }
            ]
        }
    }
    c = canonicalize_domain(sparse)
    # the singular input form normalizes into the canonical list form
    assert list(c) == ["TimeCollections"]
    assert len(c["TimeCollections"]) == 1
    tc = c["TimeCollections"][0]
    assert tc["Name"] == "TimeCollection"
    step = tc["Data"][0]
    assert step["Geometry"]["Type"] == "XYZ"
    assert step["Topology"]["DataItem"]["Format"] == "HDF"
    assert step["Topology"]["DataItem"]["DataType"] == "Float"
    assert step["Topology"]["DataItem"]["Precision"] == 4
    attr = step["Attributes"][0]
    assert attr["AttributeType"] == "Scalar"
    assert attr["Center"] == "Node"


def test_canonicalize_supplies_an_empty_attributes_list():
    domain = new_domain()
    add_uniform_step(
        domain, 0.0,
        topology("Triangle", 1, data_item([1, 3], "f.h5:/e")),
        geometry(data_item([3, 2], "f.h5:/n")),
    )
    assert canonicalize_domain(domain)["TimeCollections"][0]["Data"][0]["Attributes"] == []


def test_canonicalize_of_an_empty_domain_is_empty():
    assert canonicalize_domain({}) == {}


def test_step_cannot_be_both_uniform_and_spatial():
    domain = {
        "TimeCollection": {
            "Data": [
                {
                    "Time": 0.0,
                    "SpatialCollection": [],
                    "Topology": {"Type": "Triangle", "NumberOfElements": 1,
                                 "DataItem": {"Dimensions": [1], "Data": "f.h5:/e"}},
                }
            ]
        }
    }
    with pytest.raises(YmfArchiveError, match="not both"):
        canonicalize_domain(domain)


def test_step_must_be_either_uniform_or_spatial():
    with pytest.raises(YmfArchiveError, match="neither"):
        canonicalize_domain({"TimeCollection": {"Data": [{"Time": 0.0}]}})


@pytest.mark.parametrize("missing", ["Type", "NumberOfElements", "DataItem"])
def test_canonicalize_reports_a_missing_topology_key_by_path(missing):
    topo = topology("Triangle", 1, data_item([1, 3], "f.h5:/e"))
    del topo[missing]
    domain = {
        "TimeCollection": {
            "Data": [{"Time": 0.0, "Topology": topo,
                      "Geometry": geometry(data_item([3, 2], "f.h5:/n"))}]
        }
    }
    with pytest.raises(YmfArchiveError, match=rf"TimeCollections\[0\]\.Data\[0\]\.Topology.*{missing}"):
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
    domain["TimeCollections"][0]["Data"][0]["Topology"]["Type"] = "Dodecahedron"
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
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["Center"] = "Vertex"
    with pytest.raises(YmfArchiveError, match=r"Attributes\[0\]\.Center"):
        validate_domain(domain)


def test_validate_rejects_a_bad_attribute_type():
    domain = make_uniform_domain()
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["AttributeType"] = "Scalarr"
    with pytest.raises(YmfArchiveError, match=r"Attributes\[0\]\.AttributeType"):
        validate_domain(domain)


def test_validate_rejects_a_bad_data_type():
    domain = make_uniform_domain()
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["DataItem"]["DataType"] = "Double"
    with pytest.raises(YmfArchiveError, match="DataType"):
        validate_domain(domain)


def test_validate_rejects_duplicate_attribute_names_within_one_grid():
    # Viewers key fields by name within a grid, so a duplicate silently
    # shadows rather than erroring at read time -- worth catching on write.
    domain = make_uniform_domain()
    attrs = domain["TimeCollections"][0]["Data"][0]["Attributes"]
    attrs[1]["Name"] = attrs[0]["Name"]
    with pytest.raises(YmfArchiveError, match="duplicate attribute name"):
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
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["DataItem"]["Dimensions"] = []
    with pytest.raises(YmfArchiveError, match="must not be empty"):
        validate_domain(domain)


def test_validate_rejects_a_negative_dimension():
    domain = make_uniform_domain()
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["DataItem"]["Dimensions"] = [-1]
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
    assert read_domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["DataItem"]["Include"] == "./d/u.txt"


def test_extra_payload_round_trips_and_is_stored_readably(tmp_path):
    extra = {"Problem": {"name": "poisson", "provenance": "proteus_derived"}}
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
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["Center"] = "Nowhere"
    with pytest.raises(YmfArchiveError, match="Center"):
        write_ymf(domain, tmp_path / "archive.ymf")


def test_write_ymf_can_skip_validation_deliberately(tmp_path):
    domain = make_uniform_domain()
    domain["TimeCollections"][0]["Data"][0]["Attributes"][0]["Center"] = "Nowhere"
    path = tmp_path / "archive.ymf"
    write_ymf(domain, path, validate=False)
    assert read_ymf(path)[0]["TimeCollections"][0]["Data"][0]["Attributes"][0]["Center"] == "Nowhere"


# --------------------------------------------------------------------------
# grid fragments -- the MPI / in-HDF5 metadata path
# --------------------------------------------------------------------------


def test_a_single_grid_fragment_round_trips_as_a_string():
    g = make_spatial_domain()["TimeCollections"][0]["Data"][0]["SpatialCollection"][0]
    assert load_grid(dump_grid(g)) == g


def test_grid_fragments_reassemble_into_a_spatial_step():
    # This is the shape of the parallel write: each rank dumps its own grid
    # to a string, the strings are collected, and the master reassembles.
    original = make_spatial_domain(n_ranks=4)
    fragments = [
        dump_grid(g) for g in original["TimeCollections"][0]["Data"][0]["SpatialCollection"]
    ]
    rebuilt = new_domain(original["TimeCollections"][0]["Name"])
    add_spatial_step(rebuilt, 0.0, [load_grid(f) for f in fragments])
    assert (
        canonicalize_domain(rebuilt)["TimeCollections"][0]["Data"][0]
        == canonicalize_domain(original)["TimeCollections"][0]["Data"][0]
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
    assert [c["Name"] for c in domain["TimeCollections"]] == ["Mesh Spatial_Domain"]


def test_add_collection_appends_a_second_mesh():
    domain = _two_space_domain()
    names = [c["Name"] for c in domain["TimeCollections"]]
    assert names == ["Mesh Spatial_Domain", "Mesh_c0p2_Lagrange"]
    assert [len(c["Data"]) for c in domain["TimeCollections"]] == [2, 2]


def test_add_collection_rejects_a_duplicate_name():
    domain = new_domain("Mesh Spatial_Domain")
    with pytest.raises(YmfArchiveError, match="already has a time collection"):
        add_collection(domain, "Mesh Spatial_Domain")


def test_steps_can_target_a_collection_by_name_or_index():
    domain = _two_space_domain()
    by_name = domain["TimeCollections"][1]["Data"]
    assert all(step["Topology"]["Type"] == "Tri_6" for step in by_name)
    # index 0 is the default target
    assert all(
        step["Topology"]["Type"] == "Triangle"
        for step in domain["TimeCollections"][0]["Data"]
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
        "TimeCollections": [
            {"Name": "Mesh", "Data": []},
            {"Name": "Mesh", "Data": []},
        ]
    }
    with pytest.raises(YmfArchiveError, match="duplicate collection name"):
        canonicalize_domain(domain)


def test_singular_time_collection_is_accepted_as_input_sugar():
    singular = {"TimeCollection": {"Name": "Mesh", "Data": []}}
    plural = {"TimeCollections": [{"Name": "Mesh", "Data": []}]}
    assert canonicalize_domain(singular) == canonicalize_domain(plural)
