"""The spec notation: Unicode mathematics, parsed by ymf.symbolic.notation."""

import pytest

sympy = pytest.importorskip("sympy")

from ymf.symbolic import Space, SymbolicError, parse_equation, parse_expression  # noqa: E402
from ymf.symbolic.notation import NotationError, Parser  # noqa: E402

X, Y = sympy.symbols("x y", real=True)
T = sympy.Symbol("t", real=True)


def expr(text, **coefficients):
    return parse_expression(text, Space(2, {}, {k: sympy.sympify(v) for k, v in coefficients.items()}))


@pytest.mark.parametrize("text, expected", [
    ("π² sin(πx) sin(πy)", sympy.pi**2 * sympy.sin(sympy.pi * X) * sympy.sin(sympy.pi * Y)),
    ("2π²", 2 * sympy.pi**2),
    ("x³ + y^2 + x**4", X**3 + Y**2 + X**4),
    ("1/(2π²)", 1 / (2 * sympy.pi**2)),
    ("-x²", -X**2),                                   # unary minus binds looser than ^
    ("sin(πx)/2", sympy.sin(sympy.pi * X) / 2),
    ("x − y", X - Y),                                 # U+2212 minus sign
    ("1.0e-3 x", sympy.Float("1.0e-3") * X),
])
def test_scalar_expressions(text, expected):
    assert sympy.simplify(expr(text) - expected) == 0


def test_implicit_multiplication_and_greek_coefficients():
    value = expr("300 + 10 exp(-2π²κt) sin(πx) sin(πy)", κ=sympy.Float("0.001"))
    assert sympy.simplify(value - (300 + 10 * sympy.exp(-0.002 * sympy.pi**2 * T)
                                   * sympy.sin(sympy.pi * X) * sympy.sin(sympy.pi * Y))) == 0


def test_dot_is_a_product_between_scalars():
    assert expr("G / (2μ) · y (H - y)", G=8, μ=1, H=1) == 4 * Y * (1 - Y)


def test_tuples_are_vectors():
    v = expr("(x, 2y)")
    assert isinstance(v, sympy.MatrixBase) and list(v) == [X, 2 * Y]


# --- the differential operators -------------------------------------------------


def residuals(text, unknowns, **coefficients):
    space = Space(2, unknowns, {k: sympy.sympify(v) for k, v in coefficients.items()})
    return space, parse_equation(text, space)


def test_unicode_and_spelled_out_operators_agree():
    _, a = residuals("∂T/∂t = ∇·(κ∇T)  in Ω", {"T": 0}, κ=0.5)
    _, b = residuals("dt(T) - div(κ*grad(T)) = 0", {"T": 0}, κ=0.5)
    assert sympy.simplify(a[0] - b[0]) == 0


def test_laplacian_is_div_grad():
    _, a = residuals("-Δu = 1", {"u": 0})
    _, b = residuals("-div(grad(u)) = 1", {"u": 0})
    assert a == b


def test_navier_stokes_momentum_in_unicode():
    space, r = residuals("ρ ∂v/∂t + ∇·(ρ v⊗v) − μΔv + ∇p = f", {"v": 1, "p": 0},
                         ρ=1, μ=0.1, f=sympy.Matrix([0, 0]))
    _, s = residuals("rho*dt(v) + div(rho*outer(v, v)) - div(mu*grad(v)) + grad(p) = f",
                     {"v": 1, "p": 0}, rho=1, mu=0.1, f=sympy.Matrix([0, 0]))
    assert len(r) == 2 and all(sympy.simplify(a - b) == 0 for a, b in zip(r, s))


def test_partial_derivative_in_space():
    _, r = residuals("∂u/∂x = 0", {"u": 0})
    assert str(r[0]) == "D(u, 0)"


def test_the_in_clause_is_ignored():
    assert residuals("-Δu = 1 in Ω", {"u": 0})[1] == residuals("-Δu = 1", {"u": 0})[1]


# --- domains, predicates, assignments -------------------------------------------


def test_domains_with_and_without_time():
    p = Parser(Space(2, {}, {"H": sympy.Integer(2)}))
    assert p.domain("Ω = [0, 4] × [0, H]") == ([0.0, 0.0], [4.0, 2.0], None)
    assert p.domain("Ω = [0, 1] × [0, 1],  t ∈ (0, 50]") == ([0.0, 0.0], [1.0, 1.0], (0.0, 50.0))


def test_an_empty_interval_is_refused():
    with pytest.raises(NotationError, match="empty interval"):
        Parser(Space(2, {})).domain("Ω = [1, 0] × [0, 1]")


def test_predicates():
    p = Parser(Space(2, {}, {"H": sympy.Integer(1)}))
    assert p.predicate("∂Ω") is None
    assert p.predicate("y = H") == "(abs(y - 1) <= tol)"
    assert p.predicate("x = 0 and y = 0 or x = 4 and y = 0") == \
        "(abs(x) <= tol and abs(y) <= tol) or (abs(x - 4) <= tol and abs(y) <= tol)"


def test_assignments_with_local_definitions():
    p = Parser(Space(2, {"v": 1, "p": 0}, {"μ": sympy.Rational(1, 40)}))
    out = p.assignments("λ = 1/(2μ) - sqrt(1/(4μ²) + 4π²)\n"
                        "v = (1 - exp(λx) cos(2πy), λ/(2π) exp(λx) sin(2πy))\n"
                        "p = (1 - exp(2λx))/2", {"v", "p"})
    assert set(out) == {"v", "p"}
    lam = 20 - sympy.sqrt(400 + 4 * sympy.pi**2)
    assert sympy.simplify(out["p"] - (1 - sympy.exp(2 * lam * X)) / 2) == 0


def test_function_style_left_hand_side():
    out = Parser(Space(2, {"u": 0})).assignments("u(x, y) = sin(πx) sin(πy) / 2", {"u"})
    assert out["u"] == sympy.sin(sympy.pi * X) * sympy.sin(sympy.pi * Y) / 2


# --- errors say what and where ------------------------------------------------------


@pytest.mark.parametrize("text, message", [
    ("2 kapa", r"unknown name 'kapa'"),
    ("∇·(κ∇Q)", r"unknown name 'Q'"),
    ("sin πx", r"sin needs \(arguments\)"),
    ("(x, y) * (x, y)", r"cannot multiply two vectors"),
    ("x +", r"expected a number, a name or '\(' at position 3"),
    ("x @ y", r"unexpected character '@'"),
])
def test_errors_name_the_problem(text, message):
    with pytest.raises(SymbolicError, match=message):
        expr(text, κ=1)
