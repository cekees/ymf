"""A diffusive flux with a nonlinear potential: div(-a grad(phi(u))).

The Kirchhoff form of a nonlinear diffusion equation puts the nonlinearity in
the potential, phi(u) = u^2 + 10u, rather than in the tensor a. Space.d splits
grad(u^2 + 10u) into grad(u^2) + 10 grad(u); the classifier puts the terms of
one unknown's gradients back together as one tensor and one potential, and
the ADR form carries the potential beside the tensor (``potential``: phi, its
derivative, depends_on), as Proteus's nonlinear ``potential`` flags take it.
A potential that is the unknown itself reads exactly as before.
"""

import pytest

pytest.importorskip("sympy")
pytest.importorskip("strictyaml")

import sympy  # noqa: E402

from ymf.symbolic import SymbolicError, Space, adr_form, adr_problem, parse_equation  # noqa: E402
from ymf.symbolic.adr import classify  # noqa: E402


def _form(equation, unknowns=None, coefficients=None, dim=2, more=()):
    space = Space(dim, unknowns or {"u": 0}, coefficients or {})
    residuals = parse_equation(equation, space)
    for text in more:
        residuals.extend(parse_equation(text, space))
    return adr_form(residuals, space)


def _same(code, text):
    """A code string and an expression are the same function (10.0 and 10 alike)."""
    return sympy.simplify(sympy.nsimplify(sympy.sympify(code)) - sympy.sympify(text)) == 0


def _spec(equation, coefficients):
    return {
        "model": {
            "name": "Kirchhoff",
            "strong_form": {
                "provenance": "human_specified",
                "unknowns": ["u"],
                "equation_formulation": "steady nonlinear diffusion",
                "equations": [equation],
                "domain": "Ω = [0, 1] × [0, 1]",
                "boundary_regions": [{"name": "boundary", "where": "∂Ω"}],
                "boundary_conditions": [{"region": "boundary", "variable": "u",
                                         "type": "dirichlet", "value": 0.0}],
                "coefficients": coefficients,
            },
        },
        "solution_paths": {"discretizations": []},
    }


def test_a_polynomial_potential_is_one_tensor_and_one_potential():
    diffusion = _form("∇·(−∇(u² + 10u)) + sin(πx) = 0")["equations"][0]["diffusion"]
    assert list(diffusion) == ["u"]
    entry = diffusion["u"]
    assert entry["a"] == ["1", "1"] and entry["rowptr"] == [0, 1, 2] and entry["colind"] == [0, 1]
    assert entry["depends_on"] == {}
    potential = entry["potential"]
    assert _same(potential["phi"], "u**2 + 10*u")
    assert _same(potential["dphi"]["u"], "2*u + 10")
    assert potential["depends_on"] == {"u": "nonlinear"}


def test_the_leading_term_of_the_potential_has_coefficient_one():
    """div(-3 grad(u^2 + 10u)): a = 3 I and phi = u^2 + 10u, not a = I and phi = 3u^2 + 30u."""
    entry = _form("∇·(−3∇(u² + 10u)) = 0")["equations"][0]["diffusion"]["u"]
    assert entry["a"] == ["3", "3"]
    assert _same(entry["potential"]["phi"], "u**2 + 10*u")


def test_a_tensor_that_depends_on_the_unknown_keeps_its_own_flags():
    entry = _form("∇·(−(1 + u²)∇(u³)) = 0")["equations"][0]["diffusion"]["u"]
    assert entry["depends_on"] == {"u": "nonlinear"}
    assert _same(entry["potential"]["phi"], "u**3")
    assert _same(entry["potential"]["dphi"]["u"], "3*u**2")


def test_a_potential_that_is_the_unknown_reads_as_before():
    entry = _form("∇·(−(1 + u²)∇u) = 0")["equations"][0]["diffusion"]["u"]
    assert "potential" not in entry
    assert entry["depends_on"] == {"u": "nonlinear"}


def test_a_linear_potential_with_a_number_is_the_unknown():
    """div(-2 grad(u)): Space.d takes the number out, so the potential is u and a = 2 I."""
    entry = _form("∇·(−∇(2u)) = 0")["equations"][0]["diffusion"]["u"]
    assert "potential" not in entry
    assert entry["a"] == ["2", "2"]


def test_tensors_that_are_not_multiples_by_a_number_are_refused():
    with pytest.raises(SymbolicError, match="not multiples of each other by a number"):
        _form("∇·(−∇(u²) − x∇u) = 0")


def test_a_potential_of_two_unknowns_is_refused():
    with pytest.raises(SymbolicError, match="not a function of one unknown"):
        _form("∇·(−∇(u v)) = 0", unknowns={"u": 0, "v": 0}, more=["v = 1"])


def test_a_potential_that_depends_on_the_coordinates_is_refused():
    with pytest.raises(SymbolicError, match=r"also depends on x"):
        _form("∇·(−∇(x u²)) = 0")


def test_classify_keeps_the_potential_by_unknown():
    space = Space(2, {"u": 0})
    eq = classify(parse_equation("∇·(−∇(u² + 10u)) = 0", space)[0], space)
    assert set(eq.diffusion) == {"u"} and set(eq.potential) == {"u"}
    assert sympy.simplify(eq.potential["u"] - sympy.sympify("u**2 + 10*u",
                                                            locals={"u": space.component_symbols[0]})) == 0


@pytest.mark.parametrize("coefficients,equation", [
    ({"a": 1.0, "φ": {"formula": "u² + 10u"}, "r": {"formula": "−16π sin(4πx) sin(4πy)"}},
     "∇·(−a∇φ) + r = 0"),
    ({"r": {"formula": "−16π sin(4πx) sin(4πy)"}}, "∇·(−∇(u² + 10u)) + r = 0"),
])
def test_the_problem_is_the_same_through_a_coefficient_of_the_state_or_written_out(coefficients, equation):
    """phi as a coefficient that names the unknown, or written under the gradient: the same ADR form. Through a
    coefficient, the named form (the slots over the coefficients of the state) keeps the potential too."""
    problem = adr_problem(_spec(equation, coefficients))
    entry = problem["adr"]["equations"][0]["diffusion"]["u"]
    assert _same(entry["potential"]["phi"], "u**2 + 10*u")
    assert [float(a) for a in entry["a"]] == [1.0, 1.0]
    assert "reaction" in problem["adr"]["equations"][0]
    if "named" in problem:
        named = problem["named"]["equations"][0]["diffusion"]["u"]
        assert _same(named["phi"], "u**2 + 10*u")
        assert "reaction" in problem["named"]["equations"][0]
