"""Boundary regions named by the mesh's boundary tags (C27 Phase 5, D1).

A region is ``where`` (a predicate on the point) or ``tags`` (ints, or names
the domain's ``boundary_tags`` resolve), exactly one. A tag region reaches a
solver as ``flags``, tested against the boundary flag of the point, and has
no ``where`` key: a solver that does not know ``flags`` fails on it rather
than reading a missing predicate as the whole boundary.
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf.references import YmfReferenceError, resolve  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402

EXAMPLES = Path(__file__).parent.parent / "examples"

L_SHAPE = {
    "kind": "piecewise_linear", "dimension": 2,
    "vertices": [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0], [1.0, 2.0], [0.0, 2.0]],
    "segments": [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]],
    "segment_flags": [1, 2, 3, 3, 4, 5],
    "boundary_tags": {"bottom": 1, "right": 2, "notch": 3, "top": 4, "left": 5},
}


@pytest.fixture
def lshape():
    doc = load_ymf(EXAMPLES / "poisson.ymf")
    strong = doc["model"]["strong_form"]
    strong["domain"] = copy.deepcopy(L_SHAPE)
    strong["boundary_regions"] = [{"name": "dry", "tags": ["bottom", "left"]},
                                  {"name": "wet", "tags": [2, "notch", "top"]}]
    strong["boundary_conditions"] = [
        {"region": "dry", "variable": "u", "type": "dirichlet", "value": 0.0},
        {"region": "wet", "variable": "u", "type": "neumann", "formula": "x"}]
    doc["solution_paths"]["analytical"] = []
    doc.pop("vvuq", None)
    for d in doc["solution_paths"]["discretizations"]:
        d["mesh"] = {"h": 0.25}
    return doc


def test_a_tag_region_reaches_the_solver_as_flags(lshape):
    pytest.importorskip("sympy")
    from ymf.symbolic import adr_problem
    resolve(lshape)
    problem = adr_problem(lshape)
    (dry,) = problem["dirichlet"]
    assert dry["flags"] == [1, 5] and "where" not in dry
    (wet,) = problem["neumann"]
    assert wet["flags"] == [2, 3, 4] and wet["value"] == "x"


def test_where_and_tags_mix(lshape):
    pytest.importorskip("sympy")
    from ymf.symbolic import adr_problem
    lshape["model"]["strong_form"]["boundary_regions"][1] = {"name": "wet", "where": "x = 2"}
    (wet,) = adr_problem(lshape)["neumann"]
    assert wet["where"] == "(abs(x - 2) <= tol)" and "flags" not in wet


def _with_region(doc, region):
    doc = copy.deepcopy(doc)
    doc["model"]["strong_form"]["boundary_regions"][0] = region
    return doc


@pytest.mark.parametrize("region, message", [
    ({"name": "dry"}, "exactly one; this gives neither"),
    ({"name": "dry", "where": "x = 0", "tags": [5]}, "exactly one; this gives where and tags"),
    ({"name": "dry", "tags": ["floor"]}, "tag 'floor' is not one of the domain's boundary_tags"),
])
def test_a_region_says_how_it_is_found(lshape, region, message):
    with pytest.raises(YmfReferenceError, match=message):
        resolve(_with_region(lshape, region))


def test_a_box_has_no_tags():
    doc = load_ymf(EXAMPLES / "poisson.ymf")
    doc["model"]["strong_form"]["boundary_regions"][0] = {"name": "all_boundaries",
                                                          "tags": [1]}
    with pytest.raises(YmfReferenceError, match="a box carries no boundary tags"):
        resolve(doc)


def test_a_periodic_region_is_found_by_position(lshape):
    lshape["model"]["strong_form"]["boundary_conditions"][0]["type"] = "periodic"
    lshape["model"]["strong_form"]["boundary_conditions"][0].pop("value")
    with pytest.raises(YmfReferenceError, match="needs where:, not tags:"):
        resolve(lshape)


def test_a_tag_region_loads_from_a_file(tmp_path):
    text = (EXAMPLES / "poisson.ymf").read_text(encoding="utf-8")
    text = text.replace('''    domain: "Ω = [0, 1] × [0, 1]"''', '''    domain:
      kind: piecewise_linear
      dimension: 2
      vertices: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
      segments: [[0, 1], [1, 2], [2, 3], [3, 0]]
      segment_flags: [1, 1, 2, 1]
      boundary_tags: {sides: 1, lid: 2}''').replace('''        where: "∂Ω"''', '''        tags: [sides, lid]''')
    text = text.replace("mesh: {cells: [4, 8, 16]}", "mesh: {h: 0.25}")
    (tmp_path / "p.ymf").write_text(text, encoding="utf-8")
    region = load_ymf(tmp_path / "p.ymf")["model"]["strong_form"]["boundary_regions"][0]
    assert region == {"name": "all_boundaries", "tags": ["sides", "lid"]}


def test_split_problem_uses_flags_too(lshape):
    pytest.importorskip("sympy")
    from ymf.symbolic import split_problem
    lshape["model"]["strong_form"]["equations"] = [
        {"name": "poisson", "for": "u", "equation": "-Δu = f  in Ω"}]
    entry = lshape["solution_paths"]["discretizations"][0]
    entry["subproblems"] = [{"name": "one", "solves": ["poisson"]}]
    (sub,) = split_problem(lshape, entry["name"])["subproblems"]
    assert [c.get("flags") for c in sub["boundary_conditions"]] == [[1, 5], [2, 3, 4]]
