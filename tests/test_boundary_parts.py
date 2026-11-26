"""A state imposed on part of a boundary (C42 round 7, the SWFlow decks' partial inflow).

A condition's ``part`` is a predicate, like a region's ``where``, that says on which points of its region the condition
holds; there it replaces the region's conditions on the same variable that have no part and are on an overlapping
component, and those hold on the rest. A ``component`` may be a Cartesian one (x, y, z) as well as the normal or
tangential part. Predicates take comparisons, so a part may be a piece of a face (``x = 0 and 0.4 ≤ y ≤ 0.6``).
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")
pytest.importorskip("sympy")

from _formulation import SHALLOW_WATER, roundtrip  # noqa: E402
from ymf import closure  # noqa: E402
from ymf.boundaries import kinds  # noqa: E402
from ymf.conditions import excluded, overlap, parts  # noqa: E402
from ymf.normalize import boundary_conditions  # noqa: E402
from ymf.references import YmfReferenceError, resolve  # noqa: E402
from ymf.symbolic import SymbolicError, split_problem  # noqa: E402
from ymf.symbolic.language import Space  # noqa: E402
from ymf.symbolic.notation import NotationError, Parser  # noqa: E402

BUMP = SHALLOW_WATER / "transcritical_bump.ymf"
INFLOW = SHALLOW_WATER / "obstacle_inflow.ymf"


def _load(path):
    spec, _, _ = closure.load(path)
    return spec


def _conditions(spec, discretization="convex_limited"):
    (sub,) = split_problem(spec, discretization)["subproblems"]
    return sub["boundary_conditions"]


def _test(code, x, y=0.0, tol=1e-8):
    return eval(code, {"__builtins__": {"abs": abs}}, {"x": x, "y": y, "tol": tol})


# -- the predicate's comparisons ------------------------------------------------------------------------------------

@pytest.mark.parametrize("text, inside, outside", [
    ("x = 0 and 0.4 ≤ y ≤ 0.6", [(0, 0.4), (0, 0.5), (0, 0.6), (1e-9, 0.6 + 1e-9)], [(0, 0.39), (0, 0.61), (0.1, 0.5)]),
    ("x = 0 and 0.4 <= y <= 0.6", [(0, 0.4), (0, 0.6)], [(0, 0.3)]),
    ("x = 0 and 0.4 < y < 0.6", [(0, 0.5)], [(0, 0.4), (0, 0.6), (0, 0.6 - 1e-9)]),
    ("y ≥ 1 or x > 2", [(0, 1), (2.5, 0)], [(2, 0), (0, 0.5)]),
    ("x >= 1 and y = 0", [(1, 0), (3, 0)], [(0.5, 0), (3, 1)]),
])
def test_a_predicate_compares(text, inside, outside):
    code = Parser(Space(2, {})).predicate(text)
    assert all(_test(code, *p) for p in inside), code
    assert not any(_test(code, *p) for p in outside), code


def test_a_predicate_still_says_what_it_expected():
    with pytest.raises(NotationError, match="expected '=' or a comparison"):
        Parser(Space(2, {})).predicate("x + 1")


# -- the decks' boundaries --------------------------------------------------------------------------------------------

def test_transcritical_bump_states_the_decks_three_partial_conditions():
    spec = _load(BUMP)
    resolve(spec)
    y, x0, h25 = _conditions(spec)
    # q's y component is 0 on all of ∂Ω, x = 0 included: the part at x = 0 is on the x component, which does not overlap
    assert (y["components"], y["value"], y["where"]) == (["q_1"], "0.0", None)
    assert "part" not in y and "excluding" not in y
    assert (x0["components"], x0["value"], x0["component"]) == (["q_0"], "0.18", "x")
    assert _test(x0["part"], 0.0) and not _test(x0["part"], 25.0)
    assert h25["components"] == ["h"] and float(h25["value"]) == 0.28205279813802181
    assert _test(h25["part"], 25.0) and not _test(h25["part"], 0.0)


def test_obstacle_inflow_puts_the_state_on_x_0_and_the_outflow_on_the_rest():
    spec = _load(INFLOW)
    resolve(spec)
    h, q, rest = _conditions(spec)
    assert h["components"] == ["h"] and h["value"] == "1.0" and _test(h["part"], 0.0)
    assert q["value"] == {"q_0": "1.5", "q_1": "0"} and q["part"] == h["part"]
    # the outflow set's member: not replaced by the part (it is the rest), and told so
    assert rest["type"] == "neumann" and rest["expanded_from"] == {"type": "open", "set": "outflow"}
    assert rest["excluding"] == [q["part"]]


def test_a_part_replaces_no_condition_set_member():
    strong = _load(INFLOW)["model"]["strong_form"]
    members = [bc for bc in boundary_conditions(strong) if bc.get("set") == "outflow"]
    assert [bc["variable"] for bc in members] == ["q"]
    # a condition without a part, on the same region and variable, still does
    strong = copy.deepcopy(strong)
    strong["boundary_conditions"].append({"region": "all", "variable": "q", "type": "open"})
    assert not [bc for bc in boundary_conditions(strong) if bc.get("set") == "outflow"]


def test_a_wall_with_an_inflow_window_on_a_piece_of_a_face():
    spec = _load(INFLOW)
    strong = spec["model"]["strong_form"]
    strong["boundary_sets"] = [{"region": "all", "set": "wall"}]
    window = "x = 0 and 4 ≤ y ≤ 6"
    for bc in strong["boundary_conditions"]:
        bc["part"] = window
    resolve(spec)
    h, q, normal, tangential = _conditions(spec)
    assert (normal["component"], tangential["component"]) == ("normal", "tangential")
    # the whole-vector state overlaps both parts of the wall, so the wall holds off the window only
    assert normal["excluding"] == tangential["excluding"] == [q["part"]]
    assert _test(q["part"], 0.0, 5.0) and not _test(q["part"], 0.0, 3.0) and not _test(q["part"], 1.0, 5.0)


@pytest.mark.parametrize("a, b, expected", [
    (None, "x", True), ("x", "x", True), ("x", "y", False), ("normal", "tangential", False),
    ("normal", "x", True), ("tangential", None, True), ("y", "normal", True),
])
def test_which_components_overlap(a, b, expected):
    assert overlap(a, b) is expected and overlap(b, a) is expected


def test_excluded_follows_region_variable_and_component():
    conditions = [
        {"region": "all", "variable": "q", "type": "dirichlet", "component": "y", "value": 0.0},
        {"region": "all", "variable": "q", "type": "dirichlet", "component": "x", "value": 1.0, "part": "x = 0"},
        {"region": "all", "variable": "q", "type": "neumann", "value": 0.0},
        {"region": "all", "variable": "h", "type": "neumann", "value": 0.0},
        {"region": "top", "variable": "q", "type": "neumann", "value": 0.0},
    ]
    claimed = parts(conditions)
    assert claimed == {("all", "q"): [("x = 0", "x")]}
    assert [excluded(bc, claimed) for bc in conditions] == [[], [], ["x = 0"], [], []]


def test_a_part_of_a_tag_region_keeps_its_flags():
    spec = _load(INFLOW)
    strong = spec["model"]["strong_form"]
    strong["domain"] = {
        "kind": "piecewise_linear", "dimension": 2,
        "vertices": [[0.0, 0.0], [40.0, 0.0], [40.0, 10.0], [0.0, 10.0]],
        "segments": [[0, 1], [1, 2], [2, 3], [3, 0]], "segment_flags": [1, 2, 3, 4],
        "boundary_tags": {"bottom": 1, "right": 2, "top": 3, "left": 4}}
    strong["boundary_regions"] = [{"name": "all", "tags": ["bottom", "right", "top", "left"]}]
    resolve(spec)
    h, q, rest = _conditions(spec)
    # the region is found by flag, the part by position: both are given, and no where
    assert h["flags"] == q["flags"] == rest["flags"] == [1, 2, 3, 4]
    assert "where" not in h and "where" not in rest
    assert rest["excluding"] == [q["part"]] and _test(q["part"], 0.0) and not _test(q["part"], 40.0)
    assert [k.get("part") for k in kinds(spec)] == ["x = 0", "x = 0", None]


# -- what is refused ---------------------------------------------------------------------------------------------------

def _with_condition(path, condition):
    spec = _load(path)
    spec["model"]["strong_form"]["boundary_conditions"].append(condition)
    return spec


@pytest.mark.parametrize("condition, message", [
    ({"region": "all", "part": "∂Ω", "variable": "h", "type": "dirichlet", "value": 1.0}, "is the whole region"),
    ({"region": "all", "part": " ", "variable": "h", "type": "dirichlet", "value": 1.0}, "is the whole region"),
    ({"region": "all", "part": "x = 0", "variable": "h", "type": "periodic"}, "periodic condition .* has no part"),
    ({"region": "all", "variable": "q", "type": "slip", "component": "x"}, "Cartesian component takes a value"),
    ({"region": "all", "variable": "q", "type": "robin", "component": "y", "alpha": 1, "beta": 1, "value": 0},
     "Cartesian component"),
    ({"region": "all", "variable": "h", "type": "dirichlet", "component": "x", "value": 0.0},
     "a x component is a vector's, and h is a scalar"),
])
def test_what_a_part_or_an_axis_needs_is_checked(condition, message):
    with pytest.raises(YmfReferenceError, match=message):
        resolve(_with_condition(INFLOW, condition))


def test_a_part_that_does_not_parse_says_which_condition():
    spec = _with_condition(INFLOW, {"region": "all", "part": "x == 0", "variable": "h", "type": "dirichlet",
                                    "value": 1.0})
    with pytest.raises(SymbolicError, match="boundary condition dirichlet on all, part 'x == 0'"):
        split_problem(spec, "convex_limited")


def test_a_2d_vector_has_no_z_component():
    spec = _with_condition(INFLOW, {"region": "all", "variable": "q", "type": "dirichlet", "component": "z",
                                    "value": 0.0})
    with pytest.raises(SymbolicError, match="a 2-dimensional vector has no z component"):
        split_problem(spec, "convex_limited")


def test_one_adr_system_runs_no_part():
    from ymf.schema import load_ymf
    from ymf.symbolic import adr_problem
    doc = load_ymf(Path(__file__).parent.parent / "examples" / "poisson.ymf")
    bc = doc["model"]["strong_form"]["boundary_conditions"][0]
    bc["part"] = "x = 0"
    with pytest.raises(SymbolicError, match="on the part x = 0 runs in a split discretization's solver"):
        adr_problem(doc)


# -- an archive keeps them ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("path, realization", [(BUMP, {"cells": "80x8"}), (INFLOW, {"h": 0.5})])
def test_a_part_and_an_axis_round_trip_through_an_archive(tmp_path, path, realization):
    spec, again, digest = roundtrip(tmp_path, path, "convex_limited", realization)
    assert _conditions(again) == _conditions(spec)
    written = [bc for bc in again["model"]["strong_form"]["boundary_conditions"]]
    assert [bc.get("part") for bc in written if bc.get("part")] == [
        bc["part"] for bc in spec["model"]["strong_form"]["boundary_conditions"] if bc.get("part")]
    # and the part is in the key: another part is another input
    other = _load(path)
    next(bc for bc in other["model"]["strong_form"]["boundary_conditions"] if bc.get("part"))["part"] = "x = 1"
    assert closure.input_digest(other, "convex_limited", realization) != digest
