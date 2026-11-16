"""ymf.spaces: what a finite_element entry says beyond family and order (C8's vocabulary), checked on load."""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf.schema import load_ymf  # noqa: E402
from ymf.spaces import SpaceError, implemented_in_proteus, resolve  # noqa: E402
from ymf.vocabulary import VocabularyError, check  # noqa: E402

EXAMPLES = Path(__file__).parent.parent / "examples"


def test_defaults_are_the_ones_femtools_uses():
    assert resolve({"family": "CG", "order": 3}) == {
        "family": "CG", "order": 3, "map": "identity", "basis": "lagrange", "nodes": "gll"}
    assert resolve({"family": "CG", "order": 3, "basis": "bernstein"}) == {
        "family": "CG", "order": 3, "map": "identity", "basis": "bernstein", "bernstein_operator": False}
    assert resolve({"family": "DG", "order": 1})["basis"] == "monomial"


def test_two_ways_to_say_the_same_space_resolve_equal():
    assert resolve({"family": "CG", "order": 2}) == resolve(
        {"family": "CG", "order": 2, "basis": "lagrange", "nodes": "gll"})


def test_the_vector_families_have_a_map_and_no_basis():
    assert resolve({"family": "RT", "order": 0}) == {"family": "RT", "order": 0, "map": "contravariant_piola"}
    assert resolve({"family": "BDM", "order": 1})["map"] == "contravariant_piola"
    assert resolve({"family": "N1curl", "order": 0})["map"] == "covariant_piola"
    with pytest.raises(SpaceError, match="RT has no basis"):
        resolve({"family": "RT", "order": 0, "basis": "lagrange"})
    with pytest.raises(SpaceError, match="BDM has no nodes"):
        resolve({"family": "BDM", "order": 1, "nodes": "gll"})


@pytest.mark.parametrize("space, message", [
    ({"family": "CG", "order": 0}, "the lowest order of CG is 1"),
    ({"family": "BDM", "order": 0}, "the lowest order of BDM is 1"),
    ({"family": "RT", "order": -1}, "the lowest order of RT is 0"),
    ({"family": "CG", "order": 2, "basis": "monomial"}, "CG has no monomial basis"),
    ({"family": "CG", "order": 2, "basis": "chebyshev"}, "basis chebyshev is not one of"),
    ({"family": "CG", "order": 3, "nodes": "chebyshev"}, "nodes chebyshev is not one of"),
    ({"family": "CG", "order": 3, "basis": "bernstein", "nodes": "gll"}, "nodes belong to a lagrange basis"),
    ({"family": "CG", "order": 3, "bernstein_operator": True}, "bernstein_operator belongs to a bernstein basis"),
    ({"family": "CG", "order": 3, "nodez": "gll"}, "no key nodez"),
])
def test_a_space_that_makes_no_sense_is_refused(space, message):
    with pytest.raises(SpaceError, match=message):
        resolve(space)


def test_proteus_status_is_information_not_a_rule():
    assert implemented_in_proteus({"family": "CG", "order": 5, "nodes": "gll"})
    assert implemented_in_proteus({"family": "RT", "order": 0}) and implemented_in_proteus({"family": "BDM", "order": 1})
    assert not implemented_in_proteus({"family": "RT", "order": 1})
    assert not implemented_in_proteus({"family": "N1curl", "order": 0})
    resolve({"family": "RT", "order": 1})          # still a well-formed space


@pytest.fixture(scope="module")
def example():
    return load_ymf(EXAMPLES / "poisson_spaces.ymf")


def test_the_example_states_the_spaces_and_loads(example):
    by = {e["name"]: e["finite_element"] for e in example["solution_paths"]["discretizations"]}
    assert by["P3_gll"]["fields"]["nodes"] == "gll" and by["P3_gll"]["geometry_order"] == 1
    assert by["P4_equispaced"]["fields"]["nodes"] == "equispaced"
    assert by["P3_bernstein_exact"]["fields"]["bernstein_operator"] is False
    for finite_element in by.values():
        resolve(finite_element["fields"])


def test_a_bad_space_in_a_discretization_is_refused_on_check(example):
    doc = copy.deepcopy(example)
    entry = next(e for e in doc["solution_paths"]["discretizations"] if e["name"] == "P3_gll")
    entry["finite_element"]["fields"]["nodes"] = "chebyshev"
    with pytest.raises(VocabularyError, match=r"P3_gll\].finite_element.fields: nodes chebyshev"):
        check(doc)


def test_the_geometry_order_is_at_least_one(example):
    doc = copy.deepcopy(example)
    entry = next(e for e in doc["solution_paths"]["discretizations"] if e["name"] == "P3_gll")
    entry["finite_element"]["geometry_order"] = 0
    with pytest.raises(VocabularyError, match="geometry_order: 0"):
        check(doc)


def test_the_hdiv_pair_loads_through_the_schema(tmp_path, example):
    """RT0/DG0 is a mixed discretization; here only the finite_element block is exercised."""
    doc = copy.deepcopy(example)
    entry = next(e for e in doc["solution_paths"]["discretizations"] if e["name"] == "P3_gll")
    entry["finite_element"] = {"velocity": {"family": "RT", "order": 0}, "pressure": {"family": "DG", "order": 0}}
    check(doc)                                      # no VocabularyError
    entry["finite_element"]["velocity"]["order"] = 1
    check(doc)                                      # RT1 is well formed (proteus lacks it: information only)
