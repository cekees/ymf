"""The plain data keeps a coefficient of the state by name (C42 round 2).

adr_problem and split_problem substitute ρ = ρw (1 − χ) + ρa χ into every
slot, which is all a solver of one system needs. A kernel generator (C19)
wants ρ once, as an inline function of the state, and the slots over the
name. So beside the substituted form they give, for each coefficient of the
state, its `inline` code over the state and the other names it `uses`, and
`named`: the order to evaluate them in and the ADR slot values over them.
Evaluating the named slots after the inline coefficients gives the
substituted slots, which these tests check at random points.
"""

from pathlib import Path

import pytest

numpy = pytest.importorskip("numpy")
pytest.importorskip("sympy")
pytest.importorskip("strictyaml")

from ymf.compose import load_composed  # noqa: E402
from ymf.symbolic import adr_problem, split_problem  # noqa: E402

ROOT = Path(__file__).parent.parent
FIXTURES = Path(__file__).parent / "fixtures"
DAMBREAK = FIXTURES / "two_phase" / "dambreak.ymf"
TORINO = FIXTURES / "shallow_water" / "torino.ymf"


class _Random(dict):
    """Names a code string reads, each a random positive array on first use."""

    def __init__(self, rng):
        super().__init__()
        self.rng = rng

    def __missing__(self, name):
        if name == "numpy" or name.startswith("__"):    # a global, or pytest's own
            raise KeyError(name)
        self[name] = self.rng.uniform(0.5, 1.5, 7)
        return self[name]


def _evaluate(code, names):
    return eval(code, {"numpy": numpy}, names)   # noqa: S307 -- our own generated code


def _slots(equation):
    """(where, code) for every slot value of one equation of an ADR form."""
    for slot, key in (("mass", "m"), ("reaction", "r"), ("hamiltonian", "H")):
        if slot in equation:
            yield (equation["component"], slot), equation[slot][key]
    if "advection" in equation:
        for j, f in enumerate(equation["advection"]["f"]):
            yield (equation["component"], "advection", j), f
    for phi, d in (equation.get("diffusion") or {}).items():
        for j, a in enumerate(d["a"]):
            yield (equation["component"], "diffusion", phi, j), a


def _check_named_equals_substituted(problem):
    names = _Random(numpy.random.default_rng(42))
    for coefficient in problem["named"]["order"]:
        inline = problem["coefficients"][coefficient]["inline"]
        if isinstance(inline, list):
            for i, code in enumerate(inline):
                names["%s_%d" % (coefficient, i)] = _evaluate(code, names)
        else:
            names[coefficient] = _evaluate(inline, names)
    checked = 0
    for substituted, named in zip(problem["adr"]["equations"], problem["named"]["equations"]):
        assert dict(_slots(substituted)).keys() == dict(_slots(named)).keys()
        for (where, code), (_, named_code) in zip(_slots(substituted), _slots(named)):
            a = numpy.broadcast_to(_evaluate(code, names), (7,))
            b = numpy.broadcast_to(_evaluate(named_code, names), (7,))
            assert numpy.allclose(a, b, rtol=1e-13, atol=1e-13), (where, code, named_code)
            checked += 1
    return checked


def _sub(problem, name):
    return next(s for s in problem["subproblems"] if s["name"] == name)


def test_the_dam_break_flow_keeps_density_and_viscosity_by_name():
    flow = _sub(split_problem(load_composed(DAMBREAK), "level_set_five"), "flow")
    assert flow["named"]["order"] == ["χ", "ρ", "μ"]
    coefficients = flow["coefficients"]
    assert coefficients["χ"]["inline"] == "χh"          # read through VOF's carrier
    assert coefficients["ρ"]["inline"] == "998.2 - 996.995*χ"
    assert coefficients["ρ"]["uses"] == ["χ"]
    assert "uses" not in coefficients["χ"]
    momentum = flow["named"]["equations"][0]
    assert momentum["component"] == "v_0"
    assert momentum["mass"] == {"m": "v_0*ρ"}
    assert momentum["advection"] == {"f": ["p + v_0**2*ρ", "v_0*v_1*ρ"]}
    assert momentum["diffusion"]["v_0"]["a"] == ["2*μ", "μ"]
    assert flow["named"]["equations"][1]["reaction"] == {"r": "9.8*ρ"}
    # the derivatives stay the substituted form's, which are exact
    assert "dm" not in momentum["mass"] and "depends_on" not in momentum["mass"]
    assert _check_named_equals_substituted(flow) == 15   # 6 + 7 + 2 slot values


def test_sub_problems_without_coefficients_of_the_state_have_no_named_form():
    problem = split_problem(load_composed(DAMBREAK), "level_set_five")
    for name in ("vof", "ncls", "rdls", "mcorr"):
        assert "named" not in _sub(problem, name), name


def test_torino_keeps_its_drag_law_by_name():
    swe = _sub(split_problem(load_composed(TORINO), "convex_limited"), "swe")
    order = swe["named"]["order"]
    assert set(order) == set(swe["coefficients"])
    # each names only earlier ones
    for i, name in enumerate(order):
        assert set(swe["coefficients"][name].get("uses", [])) <= set(order[:i])
    assert isinstance(swe["coefficients"]["Fd"]["inline"], list)   # a vector
    assert _check_named_equals_substituted(swe) > 0


def test_adr_problem_names_its_coefficients_of_the_state():
    problem = adr_problem(load_composed(ROOT / "examples" / "nonlinear_diffusion.ymf"))
    assert problem["named"]["order"] == ["k"]
    assert problem["coefficients"]["k"]["inline"] == "u**2 + 1"
    assert problem["named"]["equations"][0]["diffusion"]["u"]["a"] == ["k", "k"]
    assert _check_named_equals_substituted(problem) > 0


def test_a_problem_without_coefficients_of_the_state_reads_as_before():
    problem = adr_problem(load_composed(ROOT / "examples" / "heat_equation.ymf"))
    assert "named" not in problem and "coefficients" not in problem
