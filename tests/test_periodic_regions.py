"""Periodic conditions on a piecewise_linear region (C27 Phase 5, D3).

Proteus pairs a node on a periodic face with the node the translation puts
on the opposite face, so a mesh has partners only if the boundary pieces on
the two faces match. ymf checks that at posing: the segments (2-D) or
facets (3-D) on the lower face, moved along the axis, must be exactly those
on the upper face. The mesher adding no points to them, and the partners
themselves, are checked by the solver after meshing.
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")
sympy = pytest.importorskip("sympy")

from ymf.schema import load_ymf  # noqa: E402
from ymf.symbolic import SymbolicError, adr_problem  # noqa: E402

EXAMPLES = Path(__file__).parent.parent / "examples"

#: a channel with a notch in its top, away from the periodic x faces
NOTCHED = {
    "kind": "piecewise_linear", "dimension": 2,
    "vertices": [[0.0, 0.0], [4.0, 0.0], [4.0, 1.0], [3.0, 1.0],
                 [3.0, 0.5], [1.0, 0.5], [1.0, 1.0], [0.0, 1.0]],
    "segments": [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 0]],
}

#: the unit cube as a PLC, one quadrilateral facet per face
CUBE = {
    "kind": "piecewise_linear", "dimension": 3,
    "vertices": [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                 [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]],
    "facets": [[[0, 1, 2, 3]], [[4, 5, 6, 7]], [[0, 1, 5, 4]],
               [[2, 3, 7, 6]], [[1, 2, 6, 5]], [[0, 3, 7, 4]]],
}


def _periodic(domain, ends="x = 0 or x = 4"):
    doc = load_ymf(EXAMPLES / "poisson.ymf")
    strong = doc["model"]["strong_form"]
    strong["domain"] = copy.deepcopy(domain)
    strong["boundary_regions"] = [{"name": "ends", "where": ends},
                                  {"name": "rest", "where": "∂Ω"}]
    strong["boundary_conditions"] = [
        {"region": "ends", "variable": "u", "type": "periodic"},
        {"region": "rest", "variable": "u", "type": "dirichlet", "value": 0.0}]
    doc["solution_paths"]["analytical"] = []
    doc.pop("vvuq", None)
    return doc


def test_a_notched_channel_is_periodic_along_x():
    problem = adr_problem(_periodic(NOTCHED))
    assert problem["periodic"] == [{"component": "u", "region": "ends", "axes": [0]}]
    assert problem["domain"]["region"]["segments"] == NOTCHED["segments"]


def test_a_face_cut_differently_on_the_two_sides_is_refused():
    split = copy.deepcopy(NOTCHED)
    split["vertices"].append([4.0, 0.5])                 # the x = 4 face in two
    split["segments"][1] = [1, 8]
    split["segments"].append([8, 2])
    with pytest.raises(SymbolicError, match="faces of axis x do not match: 3 of 3"):
        adr_problem(_periodic(split))


def test_a_cube_is_periodic_along_x_and_y():
    doc = _periodic(CUBE, ends="x = 0 or x = 1 or y = 0 or y = 1")
    assert adr_problem(doc)["periodic"][0]["axes"] == [0, 1]


def test_a_cube_face_with_its_facet_rotated_still_matches():
    """Facets match as point sets, not by vertex order."""
    cube = copy.deepcopy(CUBE)
    cube["facets"][4] = [[2, 6, 5, 1]]
    assert adr_problem(_periodic(cube, ends="x = 0 or x = 1"))["periodic"][0]["axes"] == [0]


def test_split_problem_checks_the_faces_too():
    from ymf.symbolic import split_problem
    doc = _periodic(NOTCHED)
    doc["model"]["strong_form"]["equations"] = [
        {"name": "poisson", "for": "u", "equation": "-Δu = f  in Ω"}]
    entry = doc["solution_paths"]["discretizations"][0]
    entry["subproblems"] = [{"name": "one", "solves": ["poisson"]}]
    (sub,) = split_problem(doc, entry["name"])["subproblems"]
    assert [c["axes"] for c in sub["boundary_conditions"] if c["type"] == "periodic"] == [[0]]
