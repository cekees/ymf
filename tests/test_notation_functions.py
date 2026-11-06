"""heaviside, sign, min, max, hypot, norm and sym in the notation.

The two-phase level-set model writes χ = heaviside(φ); its initial signed
distance is a min/max/sqrt formula; its momentum equation uses the strain
rate sym(∇v); shallow-water friction uses norm(q). heaviside and sign are
sharp: their derivatives are Dirac deltas, which a solver cannot evaluate,
so a sharp function of an unknown is refused where a derivative is needed
until a discretization regularizes it.
"""

import math

import pytest

sympy = pytest.importorskip("sympy")

from ymf.symbolic import Space, SymbolicError, adr_form, parse_equation  # noqa: E402
from ymf.symbolic import parse_expression, to_code, to_notation  # noqa: E402


def _space(**unknowns):
    return Space(2, unknowns or {"u": 0})


def _evaluate(code, **names):
    numpy = pytest.importorskip("numpy")       # not in ymf[symbolic] alone
    return eval(code, {"__builtins__": {}, "numpy": numpy}, names)  # noqa: S307


@pytest.mark.parametrize("text", [
    "heaviside(u)", "sign(u - 1)", "min(x, y)", "max(x, y, 0)",
    "min(max(x - 1.2, y - 0.6), 0)", "hypot(x, y - 1)",
])
def test_the_new_functions_round_trip(text):
    space = _space()
    value = parse_expression(text, space)
    assert parse_expression(to_notation(value), space) == value


def test_min_and_max_are_numpy_alone():
    code = to_code(parse_expression("max(x, y, 0) + min(x, 2)", _space()))
    assert "functools" not in code
    assert _evaluate(code, x=-1.0, y=-2.0) == 0.0 + -1.0


def test_the_box_signed_distance_of_the_dam_break_column():
    """The note's φ0 equals damBreakLevelSet.py's phi_init in all three branches."""
    text = "sqrt(max(x − 1.2, 0)² + max(y − 0.6, 0)²) + min(max(x − 1.2, y − 0.6), 0)"
    code = to_code(parse_expression(text, _space()))

    def phi_init(x, y):                       # test/TwoPhaseFlow/damBreakLevelSet.py
        dx, dy = x - 1.2, y - 0.6
        if dx < 0 and dy < 0:
            return max(dx, dy)
        if dx > 0 and dy > 0:
            return math.hypot(dx, dy)
        return max(dx, dy)

    for x, y in [(0.3, 0.2), (1.1, 0.55), (2.0, 1.4), (1.5, 0.1), (0.4, 1.0), (3.2, 1.8)]:
        assert math.isclose(_evaluate(code, x=x, y=y), phi_init(x, y), abs_tol=1e-15)


def test_norm_of_a_vector_and_of_a_scalar():
    space = _space(q=1)
    q0, q1 = space.component_symbols
    assert parse_expression("norm(q)", space) == sympy.sqrt(q0 ** 2 + q1 ** 2)
    assert parse_expression("norm(x - 1)", space) == sympy.Abs(space.x[0] - 1)


def test_sym_is_the_symmetric_part():
    space = _space(v=1)
    strain = parse_expression("sym(∇v)", space)
    assert strain == strain.T
    with pytest.raises(SymbolicError, match="square tensor"):
        parse_expression("sym(v)", space)


def test_the_symmetric_stress_is_a_diffusion_with_cross_terms():
    """∇·(2μ sym(∇v)) = ∇·(μ(∇v + ∇vᵀ)): v_0's equation diffuses v_1 too."""
    space = Space(2, {"v": 1, "p": 0}, {"μ": sympy.Float(0.5)})
    residuals = []
    for text in ["−∇·(2μ sym(∇v)) + ∇p = 0", "∇·v = 0"]:
        residuals.extend(parse_equation(text, space))
    momentum = adr_form(residuals, space)["equations"][0]
    assert set(momentum["diffusion"]) == {"v_0", "v_1"}


def test_a_sharp_function_of_an_unknown_has_no_derivative_to_give():
    space = _space()
    residuals = parse_equation("∂(heaviside(u))/∂t = 0", space)
    with pytest.raises(SymbolicError, match="regularizes it"):
        adr_form(residuals, space)


def test_a_sharp_function_of_the_coordinates_is_data():
    """A zone indicator in x, y: a reaction coefficient, never differentiated."""
    space = _space()
    residuals = parse_equation("∂u/∂t = ∇·∇u − heaviside(x − 0.5) u", space)
    reaction = adr_form(residuals, space)["equations"][0]["reaction"]
    assert reaction["depends_on"] == {"u": "linear"}
    assert "numpy.select" in reaction["r"]


def test_hypot_is_the_dam_break_column_bit_for_bit():
    """With hypot, φ0 is damBreakLevelSet.py's phi_init to the last bit, which the sqrt form is not
    (it rounds twice: at 2 of the 109 nodes of the deck's stored mesh it differs by one ulp), so a run
    from the spec can reproduce the deck bitwise. hypot(a, 0) is |a| exactly, so one formula covers the
    three branches."""
    numpy = pytest.importorskip("numpy")
    text = "hypot(max(x − 1.2, 0), max(y − 0.6, 0)) + min(max(x − 1.2, y − 0.6), 0)"
    value = parse_expression(text, _space())
    code = to_code(value)
    assert "numpy.hypot(" in code
    assert parse_expression(to_notation(value), _space()) == value

    def phi_init(x, y):                       # test/TwoPhaseFlow/damBreakLevelSet.py
        dx, dy = x - 1.2, y - 0.6
        if dx < 0 and dy < 0:
            return max(dx, dy)
        if dx > 0 and dy > 0:
            return numpy.hypot(dx, dy)
        return max(dx, dy)

    points = numpy.random.default_rng(0).uniform([0.0, 0.0], [3.22, 1.8], (20000, 2))
    found = _evaluate(code, x=points[:, 0], y=points[:, 1])
    assert numpy.array_equal(found, [phi_init(x, y) for x, y in points])


def test_hypot_differentiates():
    space = _space()
    value = parse_expression("hypot(x, y)", space)
    assert sympy.simplify(sympy.diff(value, space.x[0]) - space.x[0] / value) == 0
