"""slip, open, inflow and characteristic; conditions on a coefficient.

The kinds a split discretization's solver needs (C42): a free-slip wall,
an open boundary to the air, the phase entering through it, a tailwater
stage. They load and are checked; one ADR system runs none of them, and
adr_problem says which solver does.
"""

import copy
from pathlib import Path

import pytest

pytest.importorskip("strictyaml")

from ymf.references import YmfReferenceError, resolve  # noqa: E402
from ymf.schema import load_ymf  # noqa: E402

EXAMPLES = Path(__file__).parent.parent / "examples"


@pytest.fixture
def heat():
    return load_ymf(EXAMPLES / "heat_equation.ymf")


@pytest.fixture
def poiseuille():
    return load_ymf(EXAMPLES / "navier_stokes" / "plane_poiseuille.ymf")


def _with(doc, *conditions):
    doc = copy.deepcopy(doc)
    doc["model"]["strong_form"]["boundary_conditions"] = list(conditions)
    return doc


def test_the_four_kinds_load(tmp_path):
    text = (EXAMPLES / "navier_stokes" / "plane_poiseuille.ymf").read_text(encoding="utf-8")
    start = text.index("    boundary_conditions:")
    end = text.index("    coefficients:", start)
    (tmp_path / "navier_stokes_model.ymf").write_text(
        (EXAMPLES / "navier_stokes" / "navier_stokes_model.ymf").read_text(encoding="utf-8"),
        encoding="utf-8")
    (tmp_path / "p.ymf").write_text(text[:start] + """    boundary_conditions:
      - {region: walls, variable: v, type: slip}
      - {region: ends, variable: v, type: open}
      - {region: ends, variable: p, type: inflow, value: 0.0}
      - {region: pressure_datum, variable: p, type: characteristic, series: [[0, 0], [1, 1]]}
""" + text[end:], encoding="utf-8")
    kinds = [bc["type"] for bc in load_ymf(tmp_path / "p.ymf")["model"]["strong_form"][
        "boundary_conditions"]]
    assert kinds == ["slip", "open", "inflow", "characteristic"]


@pytest.mark.parametrize("condition, message", [
    ({"region": "walls", "variable": "T", "type": "slip"},
     "slip gives the normal component of a vector, and T is a scalar"),
    ({"region": "walls", "variable": "T", "type": "inflow"},
     "an inflow condition needs a value, formula or series"),
    ({"region": "walls", "variable": "T", "type": "characteristic"},
     "a characteristic condition needs a value"),
    ({"region": "walls", "variable": "T", "type": "open", "value": 1.0},
     "an open condition takes no value"),
    ({"region": "walls", "variable": "q", "type": "dirichlet", "value": 1.0},
     "q is neither an unknown nor a coefficient"),
    ({"region": "lid", "variable": "T", "type": "dirichlet", "value": 1.0},
     "no boundary region named 'lid'"),
])
def test_what_a_condition_needs_is_checked(heat, condition, message):
    with pytest.raises(YmfReferenceError, match=message):
        resolve(_with(heat, condition))


def test_a_condition_may_name_a_coefficient(heat):
    resolve(_with(heat, {"region": "walls", "variable": "κ", "type": "inflow", "value": 1.0}))


def test_slip_on_a_vector_needs_no_value(poiseuille):
    resolve(_with(poiseuille, {"region": "walls", "variable": "v", "type": "slip"}))


def test_one_adr_system_runs_none_of_them(heat, poiseuille):
    pytest.importorskip("sympy")
    from ymf.symbolic import SymbolicError, adr_problem
    with pytest.raises(SymbolicError, match="slip conditions are for a split discretization"):
        adr_problem(_with(poiseuille, {"region": "walls", "variable": "v", "type": "slip"}))
    with pytest.raises(SymbolicError, match="a condition on a coefficient is for a split"):
        adr_problem(_with(heat, {"region": "walls", "variable": "κ", "type": "dirichlet",
                                 "value": 1.0}))


def test_the_units_of_an_inflow_value_are_those_of_its_variable(heat):
    pytest.importorskip("pint")
    from ymf import check_units
    wrong = _with(heat, {"region": "walls", "variable": "T", "type": "inflow",
                         "value": 300.0, "units": "m"})
    assert not check_units(wrong)
    right = _with(heat, {"region": "walls", "variable": "T", "type": "inflow",
                         "value": 300.0, "units": "K"})
    assert check_units(right)
