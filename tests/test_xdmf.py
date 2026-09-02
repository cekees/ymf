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
    topo_data_item = canonical["TimeCollection"]["Data"][0]["Topology"]["DataItem"]
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
