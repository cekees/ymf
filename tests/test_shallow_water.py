"""The shallow-water family (C42's second, C35's model) as ymf documents.

tests/fixtures/shallow_water: the model (depth h, discharge q; porosity,
Manning friction, Darcy-Forchheimer drag as coefficients of the state), the
numerics Proteus's SW2DCV uses (one sub-problem), a lake at rest over a bump
(well-balancing, exact), and the Torino flume (a dam break down a 3 % slope
into a rod canopy, from torino_flume/dswe/torino.py).
"""

import math
from pathlib import Path

import pytest

sympy = pytest.importorskip("sympy")
pytest.importorskip("strictyaml")

from ymf.schema import load_ymf  # noqa: E402
from ymf.state import reached  # noqa: E402
from ymf.symbolic import D, Dt, split_problem, strong_form_equations  # noqa: E402

SW = Path(__file__).parent / "fixtures" / "shallow_water"


@pytest.fixture(scope="module")
def lake():
    return load_ymf(SW / "lake_at_rest.ymf")


@pytest.fixture(scope="module")
def torino():
    return load_ymf(SW / "torino.ymf")


def _swe(doc):
    (sub,) = split_problem(doc, "convex_limited")["subproblems"]
    return sub


def _equation(sub, component):
    return next(e for e in sub["adr"]["equations"] if e["component"] == component)


def _evaluate(code, **names):
    numpy = pytest.importorskip("numpy")       # not in ymf[symbolic] alone
    return eval(code, {"__builtins__": {}, "numpy": numpy}, names)  # noqa: S307


def test_both_cases_are_discretized(lake, torino):
    assert reached(lake)[0] == "discretized" and reached(torino)[0] == "discretized"
    # a case changes a free-form scheme option of the shared numerics
    assert torino["solution_paths"]["discretizations"][0]["discretization"][
        "entropy_viscosity_cE"] == "5.0"


def test_shallow_water_is_one_adr_system(lake):
    swe = _swe(lake)
    assert swe["unknowns"] == {"h": ["h"], "q": ["q_0", "q_1"]}
    momentum = _equation(swe, "q_0")
    assert momentum["advection"]["f"][0] == "4.905*h**2 + 1.0*q_0**2/h"
    assert "q_0*numpy.sqrt(q_0**2 + q_1**2)/h**(7/3)" in momentum["reaction"]["r"]
    assert swe["coefficients"]["τb"]["depends_on"] == {
        "h": "nonlinear", "q_0": "nonlinear", "q_1": "nonlinear"}
    # the wall set's slip, expanded: q.n = 0 and no tangential traction
    assert [(c["type"], c["variable"], c.get("component")) for c in swe["boundary_conditions"]] == [
        ("dirichlet", "q", "normal"), ("neumann", "q", "tangential")]


def test_the_lake_at_rest_solves_the_model_exactly(lake):
    """Well-balancing of the strong form: h = 0.5 - Z, q = 0 makes every
    residual vanish identically, the bed's kink included."""
    space, residuals = strong_form_equations(lake)
    bed = space.coefficients["Z"]                  # the spec's own bed
    state = {sympy.Symbol("h", real=True): sympy.Float(0.5) - bed,
             sympy.Symbol("q_0", real=True): 0, sympy.Symbol("q_1", real=True): 0}

    def substitute(r):
        r = r.replace(D, lambda e, axis: sympy.diff(sympy.sympify(e).subs(state),
                                                     space.x[int(axis)]))
        r = r.replace(Dt, lambda e: sympy.diff(sympy.sympify(e).subs(state), space.t))
        return sympy.simplify(r.subs(state).doit())

    assert [substitute(r) for r in residuals] == [0, 0, 0]


def test_the_bed_kink_is_not_regularized(lake):
    """Only sharp functions as written are smoothed; the bump's max gives a
    bounded step in ∇Z, which is left as it is."""
    swe = _swe(lake)
    assert "regularized" not in swe
    assert "[0,1/2,1]" in _equation(swe, "q_0")["reaction"]["r"]


def test_the_torino_canopy(torino):
    swe = _swe(torino)
    assert swe["regularized"] == {"heaviside": {"by": "smoothed_heaviside", "width": "he"}}
    porosity = 1 - 194.0 * math.pi * 0.006 ** 2 / 4
    mass = _equation(swe, "h")["mass"]["m"]
    assert math.isclose(_evaluate(mass, h=1.0, x=7.0, he=0.0725), porosity, rel_tol=1e-12)
    assert _evaluate(mass, h=1.0, x=2.0, he=0.0725) == 1.0
    assert [(c["type"], c["region"]) for c in swe["boundary_conditions"]] == [
        ("dirichlet", "walls"), ("neumann", "walls"), ("neumann", "open_end")]


def test_the_torino_drag_is_the_decks_froude_law(torino):
    """Fd = (α + β|q|/h) q with β = C_D min(1, 0.9 Fr^-0.8) m D / 2 in the
    canopy, Fr = |q|/(h sqrt(g h)): SW2DCV's -(α + β|u|) h u."""
    swe = _swe(torino)
    h, q0, x = 0.1, 0.05, 7.0
    fr = (q0 / h) / math.sqrt(9.81 * h)
    beta = 1.0 * min(1.0, 0.9 * fr ** -0.8) * 194.0 * 0.006 / 2
    expected = beta * abs(q0) / h * q0
    fd = swe["coefficients"]["Fd"]["code"][0]
    assert math.isclose(_evaluate(fd, h=h, q_0=q0, q_1=0.0, x=x, he=0.0725), expected,
                        rel_tol=1e-12)
    assert _evaluate(fd, h=h, q_0=q0, q_1=0.0, x=2.0, he=0.0725) == 0.0   # outside


def test_both_cases_pass_the_unit_check(lake, torino):
    pytest.importorskip("pint")
    from ymf import check_units
    assert check_units(lake) and check_units(torino)
