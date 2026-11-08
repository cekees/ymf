"""A condition on a phase indicator no sub-problem carries (C42 round 3,
decision 18): with CLSVOF solving the level set alone, the atmosphere's
χ inflow (air enters) is a condition on the sign of φ.
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("sympy")
pytest.importorskip("strictyaml")

from ymf.schema import load_ymf  # noqa: E402
from ymf.symbolic.language import SymbolicError  # noqa: E402
from ymf.symbolic.split import split_problem  # noqa: E402

TWO_PHASE = Path(__file__).parent / "fixtures" / "two_phase"


@pytest.fixture
def doc():
    return load_ymf(TWO_PHASE / "dambreak_clsvof.ymf")


def _conditions(split, name):
    return next(s for s in split["subproblems"] if s["name"] == name)["boundary_conditions"]


def test_clsvof_takes_the_atmospheres_inflow_as_the_sign_of_phi(doc):
    split = split_problem(doc, "level_set_clsvof")
    (bc,) = _conditions(split, "clsvof")
    assert bc["variable"] == "φ" and bc["components"] == ["φ"]
    assert (bc["type"], bc["only"], bc["region"], bc["sign"]) == ("dirichlet", "inflow", "top", 1)
    assert bc["indicator"] == {"variable": "χ", "value": 1.0}
    assert "value" not in bc
    assert bc["expanded_from"] == {"type": "inflow", "set": "atmosphere"}
    # the flow keeps its own: p = 0 and the open velocity, nothing on χ
    assert {c["variable"] for c in _conditions(split, "flow")} == {"p", "v"}


def test_with_a_carrier_nothing_changes(doc):
    """level_set_five: VOF carries χ, so the inflow is its χh's, and the
    level set and its redistancing take none (as the deck)."""
    split = split_problem(doc, "level_set_five")
    (vof,) = _conditions(split, "vof")
    assert (vof["variable"], vof["carried_by"], vof["value"]) == ("χ", "χh", "1.0")
    assert _conditions(split, "ncls") == [] and _conditions(split, "rdls") == []


@pytest.mark.parametrize("value, sign", [(0.0, -1), (0.5, 0), (1.0, 1)])
def test_each_value_of_the_indicator_is_a_sign(doc, value, sign):
    doc = copy.deepcopy(doc)
    atmosphere = doc["model"]["strong_form"]["condition_sets"]["atmosphere"]
    next(c for c in atmosphere if c["variable"] == "χ")["value"] = value
    (bc,) = _conditions(split_problem(doc, "level_set_clsvof"), "clsvof")
    assert bc["sign"] == sign


def test_a_value_no_sign_says_is_refused(doc):
    doc = copy.deepcopy(doc)
    atmosphere = doc["model"]["strong_form"]["condition_sets"]["atmosphere"]
    next(c for c in atmosphere if c["variable"] == "χ")["value"] = 0.3
    with pytest.raises(SymbolicError, match="χ is heaviside.φ., whose values are 0, 1/2 and 1; "
                                            "0.3 says no sign of φ"):
        split_problem(doc, "level_set_clsvof")
