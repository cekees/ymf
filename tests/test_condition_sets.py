"""Condition sets: written once in a model, applied to regions by a problem.

The 25 air-water-vv two-phase cases are written in proteus's boundary-
condition setters (setFreeSlip, setAtmosphere, ...), each a bundle of
conditions on every model's variables. A ymf model names such bundles over
its own unknowns and coefficients; a case says where each holds.
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf.normalize import boundary_conditions  # noqa: E402
from ymf.references import YmfReferenceError, resolve  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402

TWO_PHASE = Path(__file__).parent / "fixtures" / "two_phase"
EXAMPLES = Path(__file__).parent.parent / "examples"


@pytest.fixture
def dambreak():
    return load_ymf(TWO_PHASE / "dambreak.ymf")


def _strong(doc):
    return doc["model"]["strong_form"]


def test_the_dam_break_applies_free_slip_and_atmosphere(dambreak):
    assert [(bc["region"], bc["variable"], bc["type"], bc["set"])
            for bc in boundary_conditions(_strong(dambreak))] == [
        ("walls", "v", "slip", "free_slip"),
        ("top", "p", "dirichlet", "atmosphere"),
        ("top", "v", "open", "atmosphere"),
        ("top", "χ", "inflow", "atmosphere"),
    ]


def test_a_condition_written_out_replaces_the_sets_member(dambreak):
    strong = _strong(dambreak)
    strong["boundary_conditions"] = [
        {"region": "top", "variable": "p", "type": "dirichlet", "formula": "9.8 t"}]
    conditions = boundary_conditions(strong)
    pressure = [bc for bc in conditions if bc["region"] == "top" and bc["variable"] == "p"]
    assert pressure == [{"region": "top", "variable": "p", "type": "dirichlet",
                         "formula": "9.8 t"}]
    assert len(conditions) == 4


def test_an_undefined_set_is_refused(dambreak):
    _strong(dambreak)["boundary_sets"].append({"region": "walls", "set": "no_slip"})
    with pytest.raises(YmfReferenceError, match="no condition set named 'no_slip'"):
        resolve(dambreak)


def test_two_sets_may_not_both_decide_one_regions_variable(dambreak):
    strong = _strong(dambreak)
    strong["condition_sets"]["no_slip"] = [{"variable": "v", "type": "dirichlet", "value": 0.0}]
    strong["boundary_sets"].append({"region": "walls", "set": "no_slip"})
    with pytest.raises(YmfReferenceError, match="free_slip and no_slip both give v"):
        resolve(dambreak)


def test_a_sets_members_are_checked_where_they_are_defined(dambreak):
    strong = _strong(dambreak)
    strong["condition_sets"]["broken"] = [{"variable": "p", "type": "slip"}]
    with pytest.raises(YmfReferenceError,
                       match=r"condition_sets\.broken\[variable='p'\]: slip .* p is a scalar"):
        resolve(dambreak)


def test_a_set_is_what_its_conditions_would_be(tmp_path):
    """The heat equation's wall condition as a set gives the same problem."""
    pytest.importorskip("sympy")
    from ymf.symbolic import adr_problem
    plain = load_ymf(EXAMPLES / "heat_equation.ymf")
    with_set = copy.deepcopy(plain)
    strong = _strong(with_set)
    (wall,) = strong.pop("boundary_conditions")
    strong["condition_sets"] = {"held": [{k: v for k, v in wall.items() if k != "region"}]}
    strong["boundary_sets"] = [{"region": wall["region"], "set": "held"}]
    a, b = adr_problem(plain), adr_problem(with_set)
    assert a["dirichlet"] == b["dirichlet"] and a["adr"] == b["adr"]


def test_the_units_check_says_which_set_a_condition_came_from(dambreak):
    pytest.importorskip("pint")
    from ymf import check_units
    strong = _strong(dambreak)
    strong["condition_sets"]["atmosphere"][0] = {"variable": "p", "type": "dirichlet",
                                                 "value": 0.0, "units": "m"}
    issues = [i for i in check_units(dambreak).issues if i.severity == "error"]
    assert any("(set atmosphere)" in i.location for i in issues)
