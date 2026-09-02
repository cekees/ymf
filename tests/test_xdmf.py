from xml.etree.ElementTree import Element

import pytest

from ymf.xdmf import (
    build_xdmf_tree,
    canonicalize_domain,
    parse_xdmf_domain,
    parse_xdmf_extra,
    read_xdmf,
    round_trip_equal,
    write_xdmf,
)

SAMPLE_DOMAIN = {
    "TimeCollection": {
        "Name": "Mesh Spatial_Domain",
        "Data": [
            {
                "Time": 0.0,
                "Topology": {
                    "Type": "Tetrahedron",
                    "NumberOfElements": 687737,
                    "DataItem": {
                        "Format": "HDF",
                        "DataType": "Int",
                        "Dimensions": [687737, 4],
                        "Data": "poisson_3.h5:/elementsSpatial_Domain0",
                    },
                },
                "Geometry": {
                    "Type": "XYZ",
                    "DataItem": {
                        "Format": "HDF",
                        "DataType": "Float",
                        "Precision": 8,
                        "Dimensions": [118084, 3],
                        "Data": "poisson_3.h5:/nodesSpatial_Domain0",
                    },
                },
            },
            {
                "Time": 1.0,
                "Topology": {
                    "Type": "Tetrahedron",
                    "NumberOfElements": 687737,
                    "DataItem": {
                        "Format": "HDF",
                        "DataType": "Int",
                        "Dimensions": [687737, 4],
                        "Data": "poisson_3.h5:/elementsSpatial_Domain1",
                    },
                },
                "Geometry": {
                    "Type": "XYZ",
                    "DataItem": {
                        "Format": "HDF",
                        "DataType": "Float",
                        "Precision": 8,
                        "Dimensions": [118084, 3],
                        "Data": "poisson_3.h5:/nodesSpatial_Domain1",
                    },
                },
            },
        ],
    }
}


# ---------------------------------------------------------------------------
# domain-only round trip (no ymf_extra)
# ---------------------------------------------------------------------------


def test_build_xdmf_tree_basic_structure():
    tree = build_xdmf_tree(SAMPLE_DOMAIN)
    root = tree.getroot()
    assert root.tag == "Xdmf"
    assert root.attrib["Version"] == "2.0"
    domain_elem = root.find("Domain")
    assert domain_elem is not None
    collection = domain_elem.find("Grid")
    assert collection.attrib["GridType"] == "Collection"
    assert collection.attrib["CollectionType"] == "Temporal"
    assert len(collection.findall("Grid")) == 2


def test_domain_round_trips_via_parse_xdmf_domain():
    tree = build_xdmf_tree(SAMPLE_DOMAIN)
    parsed = parse_xdmf_domain(tree.getroot())
    assert parsed == canonicalize_domain(SAMPLE_DOMAIN)


def test_write_then_read_xdmf_domain_only(tmp_path):
    path = tmp_path / "test.xmf"
    write_xdmf(SAMPLE_DOMAIN, path)
    domain, ymf_extra = read_xdmf(path)
    assert domain == canonicalize_domain(SAMPLE_DOMAIN)
    assert ymf_extra is None


def test_canonicalize_domain_is_idempotent():
    once = canonicalize_domain(SAMPLE_DOMAIN)
    twice = canonicalize_domain(once)
    assert once == twice


def test_canonicalize_domain_fills_precision_default():
    # SAMPLE_DOMAIN's Topology DataItem has no explicit Precision.
    canonical = canonicalize_domain(SAMPLE_DOMAIN)
    topo_data_item = canonical["TimeCollections"][0]["Data"][0]["Topology"]["DataItem"]
    assert topo_data_item["Precision"] == 4  # XDMF default for Int


def test_round_trip_equal_domain_only():
    assert round_trip_equal(SAMPLE_DOMAIN) is True


def test_parse_xdmf_domain_returns_empty_dict_without_domain_element():
    root = Element("Xdmf", {"Version": "2.0"})
    assert parse_xdmf_domain(root) == {}


def test_no_information_element_when_ymf_extra_omitted(tmp_path):
    path = tmp_path / "test.xmf"
    write_xdmf(SAMPLE_DOMAIN, path)
    tree = __import__("xml.etree.ElementTree", fromlist=["parse"]).parse(path)
    domain_elem = tree.getroot().find("Domain")
    assert domain_elem.find("Information") is None


# ---------------------------------------------------------------------------
# domain + ymf_extra round trip
# ---------------------------------------------------------------------------

SAMPLE_YMF_EXTRA = {
    "Problem": {
        "name": "Poisson 3D (unit cube, tetrahedral mesh)",
        "strong_form": {
            "unknowns": [{"name": "u", "units": "K", "std_name": "land_surface_temperature"}],
            "unknown_provenance": "llm_derived",
            "equation_formulation": "Poisson",
            "strong_form_expression": "Delta u = f",
            "domain": "unit cube",
            "boundary_regions": [{"name": "all", "geometry": "boundary"}],
        },
    },
    "solution_paths": {
        "analytical": [],
        "discretizations": [
            {
                "name": "P1_tetrahedral",
                "provenance": "proteus_derived",
                "solver": {"type": "linear"},
            }
        ],
    },
}


def test_information_element_present_when_ymf_extra_given(tmp_path):
    path = tmp_path / "test.xmf"
    write_xdmf(SAMPLE_DOMAIN, path, ymf_extra=SAMPLE_YMF_EXTRA)
    import xml.etree.ElementTree as ET

    tree = ET.parse(path)
    domain_elem = tree.getroot().find("Domain")
    info = domain_elem.find("Information")
    assert info is not None
    assert info.attrib["Name"] == "YMF"
    # value should be base64 (no raw YAML/unicode leaking into the attribute)
    assert info.attrib["Value"].isascii()


def test_write_then_read_xdmf_with_ymf_extra(tmp_path):
    path = tmp_path / "test.xmf"
    write_xdmf(SAMPLE_DOMAIN, path, ymf_extra=SAMPLE_YMF_EXTRA)
    domain, ymf_extra = read_xdmf(path)
    assert domain == canonicalize_domain(SAMPLE_DOMAIN)
    assert ymf_extra == SAMPLE_YMF_EXTRA


def test_round_trip_equal_with_ymf_extra():
    assert round_trip_equal(SAMPLE_DOMAIN, ymf_extra=SAMPLE_YMF_EXTRA) is True


def test_parse_xdmf_extra_returns_none_without_information():
    tree = build_xdmf_tree(SAMPLE_DOMAIN)
    assert parse_xdmf_extra(tree.getroot()) is None


def test_parse_xdmf_extra_returns_none_without_domain_element():
    root = Element("Xdmf", {"Version": "2.0"})
    assert parse_xdmf_extra(root) is None


def test_ymf_extra_with_unicode_survives_round_trip(tmp_path):
    # strong_form_expression fields routinely contain unicode math (∂, ∇, Ω, ...)
    extra = {
        "strong_form_expression": "∂u/∂t + u·∇u = -∇p + νΔu + f\n∇·u = 0",
        "domain": "Ω = [0,10] × [0,1]",
    }
    path = tmp_path / "test.xmf"
    write_xdmf(SAMPLE_DOMAIN, path, ymf_extra=extra)
    _, ymf_extra = read_xdmf(path)
    assert ymf_extra == extra


def test_ymf_extra_with_special_xml_chars_survives_round_trip(tmp_path):
    # <, >, &, and quotes would break a naive (unescaped) XML attribute;
    # base64 encoding sidesteps XML escaping entirely.
    extra = {"notes": "if a < b & c > d then \"quote\" and 'apostrophe'"}
    path = tmp_path / "test.xmf"
    write_xdmf(SAMPLE_DOMAIN, path, ymf_extra=extra)
    _, ymf_extra = read_xdmf(path)
    assert ymf_extra == extra


# ---------------------------------------------------------------------------
# full YMF document round trip (the actual notebook example)
# ---------------------------------------------------------------------------


def test_full_notebook_example_round_trips(tmp_path):
    """End-to-end: the same Problem doc used in YMF-updated.ipynb, plus the
    mesh archive, written and read back through write_xdmf/read_xdmf."""
    from ymf import validate_ymf

    poisson_problem_yaml = """
Problem:
  name: "Poisson 3D (unit cube, tetrahedral mesh)"
  physical_model:
    provenance: llm_derived
    processes: [diffusion]
    assumptions: [steady, constant_properties]
  strong_form:
    unknowns: [u]
    unknown_provenance: llm_derived
    equation_formulation: "Poisson"
    strong_form_expression: |
      Delta u = f
    domain: |
      unit cube [0,1]^3
    boundary_regions:
      - name: "all_boundaries"
        geometry: "boundary"
    coefficients:
      f: "1.0"
solution_paths:
  analytical: []
  discretizations:
    - name: "P1_tetrahedral"
      provenance: proteus_derived
      finite_element:
        fields: {family: CG, order: 1}
      solver:
        type: linear
        linear_solver: petsc
        tolerance: 1e-10
"""
    problem_doc = validate_ymf(poisson_problem_yaml)

    path = tmp_path / "poisson_full.xmf"
    write_xdmf(SAMPLE_DOMAIN, path, ymf_extra=problem_doc)
    domain, ymf_extra = read_xdmf(path)

    assert domain == canonicalize_domain(SAMPLE_DOMAIN)
    assert ymf_extra == problem_doc
    assert ymf_extra["Problem"]["name"] == "Poisson 3D (unit cube, tetrahedral mesh)"


# ==========================================================================
# Widened domain model: attributes, spatial collections, NodesPerElement,
# and text-include DataItems. These are the parts Proteus needs that the
# original TimeCollection/Topology/Geometry model did not cover.
# ==========================================================================

from ymf.archive import (  # noqa: E402  (grouped with the tests that use it)
    add_spatial_step,
    add_uniform_step,
    attribute,
    data_item,
    geometry,
    grid,
    new_domain,
    topology,
)


def _uniform_with_attributes():
    domain = new_domain("Mesh Spatial_Domain")
    add_uniform_step(
        domain,
        0.0,
        topology("Tetrahedron", 12, data_item([12, 4], "out.h5:/elements_t0", data_type="Int")),
        geometry(data_item([9, 3], "out.h5:/nodes_t0", precision=8)),
        [
            attribute("u", data_item([9], "out.h5:/u_t0", precision=8)),
            attribute(
                "velocity",
                data_item([9, 3], "out.h5:/vel_t0", precision=8),
                attribute_type="Vector",
            ),
            attribute(
                "porosity",
                data_item([12], "out.h5:/por_t0", precision=8),
                center="Cell",
            ),
        ],
    )
    return domain


def _spatial(n_ranks=3):
    domain = new_domain("Mesh Spatial_Domain")
    grids = [
        grid(
            topology("Triangle", 10 + r, data_item([10 + r, 3], f"out.h5:/e_p{r}", data_type="Int")),
            geometry(data_item([8 + r, 2], f"out.h5:/n_p{r}", precision=8), geometry_type="XY"),
            [attribute("p", data_item([8 + r], f"out.h5:/p_p{r}", precision=8))],
            name=f"subdomain_{r}",
        )
        for r in range(n_ranks)
    ]
    domain_with_step = domain
    add_spatial_step(domain_with_step, 0.0, grids)
    return domain_with_step


def test_attributes_are_written_as_xdmf_attribute_elements():
    root = build_xdmf_tree(_uniform_with_attributes()).getroot()
    grid_elem = root.find("Domain").find("Grid").find("Grid")
    attrs = grid_elem.findall("Attribute")
    assert [a.attrib["Name"] for a in attrs] == ["u", "velocity", "porosity"]
    assert attrs[1].attrib["AttributeType"] == "Vector"
    assert attrs[2].attrib["Center"] == "Cell"
    # each carries its own DataItem, with the HDF5 reference as its text
    assert attrs[0].find("DataItem").text == "out.h5:/u_t0"


def test_attributes_round_trip():
    domain = _uniform_with_attributes()
    assert parse_xdmf_domain(build_xdmf_tree(domain).getroot()) == canonicalize_domain(domain)


def test_attributes_round_trip_through_a_file():
    assert round_trip_equal(_uniform_with_attributes()) is True


def test_spatial_collection_is_written_as_a_nested_collection_grid():
    root = build_xdmf_tree(_spatial()).getroot()
    temporal = root.find("Domain").find("Grid")
    assert temporal.attrib["CollectionType"] == "Temporal"
    spatial = temporal.find("Grid")
    assert spatial.attrib["GridType"] == "Collection"
    assert spatial.attrib["CollectionType"] == "Spatial"
    # Time hangs off the spatial collection, once, not off each subdomain.
    assert spatial.find("Time").attrib["Value"] == "0.0"
    subdomains = spatial.findall("Grid")
    assert [g.attrib["Name"] for g in subdomains] == [
        "subdomain_0", "subdomain_1", "subdomain_2"
    ]
    assert all(g.find("Time") is None for g in subdomains)


def test_spatial_collection_round_trips():
    domain = _spatial()
    assert parse_xdmf_domain(build_xdmf_tree(domain).getroot()) == canonicalize_domain(domain)


def test_spatial_collection_round_trips_through_a_file():
    assert round_trip_equal(_spatial(n_ranks=4)) is True


def test_nodes_per_element_is_written_and_round_trips():
    domain = new_domain()
    add_uniform_step(
        domain,
        0.0,
        topology("Polyline", 5, data_item([10], "out.h5:/e", data_type="Int"), nodes_per_element=2),
        geometry(data_item([10, 3], "out.h5:/n", precision=8)),
    )
    root = build_xdmf_tree(domain).getroot()
    topo_elem = root.find("Domain").find("Grid").find("Grid").find("Topology")
    assert topo_elem.attrib["NodesPerElement"] == "2"
    assert parse_xdmf_domain(root) == canonicalize_domain(domain)


def test_absent_nodes_per_element_stays_absent():
    domain = new_domain()
    add_uniform_step(
        domain,
        0.0,
        topology("Tetrahedron", 5, data_item([5, 4], "out.h5:/e", data_type="Int")),
        geometry(data_item([4, 3], "out.h5:/n", precision=8)),
    )
    root = build_xdmf_tree(domain).getroot()
    topo_elem = root.find("Domain").find("Grid").find("Grid").find("Topology")
    assert "NodesPerElement" not in topo_elem.attrib
    assert "NodesPerElement" not in parse_xdmf_domain(root)["TimeCollections"][0]["Data"][0]["Topology"]


def _text_include_domain():
    domain = new_domain()
    add_uniform_step(
        domain,
        0.0,
        topology("Triangle", 2, data_item([2, 3], include="./d/elements.txt", data_type="Int")),
        geometry(data_item([4, 2], include="./d/nodes.txt", precision=8)),
        [attribute("u", data_item([4], include="./d/u.txt", precision=8))],
    )
    return domain


def test_text_include_data_items_get_an_xi_include_child():
    root = build_xdmf_tree(_text_include_domain()).getroot()
    grid_elem = root.find("Domain").find("Grid").find("Grid")
    u_item = grid_elem.find("Attribute").find("DataItem")
    assert u_item.attrib["Format"] == "XML"
    include = list(u_item)[0]
    assert include.tag == "xi:include"
    assert include.attrib == {"parse": "text", "href": "./d/u.txt"}
    # the reference lives in the child, not in the element's own text
    assert not (u_item.text or "").strip()


def test_each_text_include_attaches_to_its_own_data_item():
    # Regression guard for the class of bug found in Proteus's Archiver.py,
    # where a sidecar reference was attached to the previous DataItem,
    # leaving one field with no data and another with two references.
    root = build_xdmf_tree(_text_include_domain()).getroot()
    grid_elem = root.find("Domain").find("Grid").find("Grid")
    items = {
        "topology": grid_elem.find("Topology").find("DataItem"),
        "geometry": grid_elem.find("Geometry").find("DataItem"),
        "attribute": grid_elem.find("Attribute").find("DataItem"),
    }
    hrefs = {}
    for role, item in items.items():
        includes = [c for c in item if c.tag == "xi:include"]
        assert len(includes) == 1, f"{role} DataItem has {len(includes)} xi:include children"
        hrefs[role] = includes[0].attrib["href"]
    assert hrefs == {
        "topology": "./d/elements.txt",
        "geometry": "./d/nodes.txt",
        "attribute": "./d/u.txt",
    }


def test_text_include_data_items_round_trip():
    domain = _text_include_domain()
    assert parse_xdmf_domain(build_xdmf_tree(domain).getroot()) == canonicalize_domain(domain)


def test_text_include_round_trips_through_a_file():
    assert round_trip_equal(_text_include_domain()) is True


def test_mixed_hdf_and_text_data_items_in_one_grid_round_trip():
    domain = new_domain()
    add_uniform_step(
        domain,
        0.0,
        topology("Triangle", 2, data_item([2, 3], "out.h5:/e", data_type="Int")),
        geometry(data_item([4, 2], include="./d/nodes.txt", precision=8)),
        [attribute("u", data_item([4], "out.h5:/u", precision=8))],
    )
    assert round_trip_equal(domain) is True
