"""ymf.symbolic.printer: sympy back to the notation, as an exact round trip."""

import random
from pathlib import Path

import pytest

sympy = pytest.importorskip("sympy")
pytest.importorskip("strictyaml")

from ymf.compose import load_composed
from ymf.symbolic import Space, equation_notation, parse_equation, to_notation
from ymf.symbolic.notation import NotationError, Parser

EXAMPLES = Path(__file__).parent.parent / "examples"
SPACE = Space(2, {"u": 0, "v": 1}, {"κ": sympy.Symbol("κ", real=True)})


def round_trip(expr, space=SPACE):
    text = to_notation(expr)
    back = Parser(space).expression(text)
    assert back == expr, "%s\n  printed %r\n  read back %s" % (expr, text, back)
    return text


@pytest.mark.parametrize("written, printed", [
    ("-Δu", "-∂(∂u/∂x)/∂x - ∂(∂u/∂y)/∂y"),
    ("∇·(κ∇u)", "∂(κ ∂u/∂x)/∂x + ∂(κ ∂u/∂y)/∂y"),
    ("∂u/∂t + ∇·(u v)", "∂(u v_0)/∂x + ∂(u v_1)/∂y + ∂u/∂t"),
    ("2π² sin(πx) sin(πy)", "2 π² sin(π x) sin(π y)"),
    ("sin(πx)/(2π²)", "sin(π x)/(2 π²)"),
    ("u/(x y) - 3/4", "u/(x y) - 3/4"),
    ("x^(1/3) + 1.0e-10 u", "1e-10 u + x^(1/3)"),
    ("(u, x²)", "(u, x²)"),
    ("exp(1) + 0.025 u³ - u^5", "-u^5 + 0.025 u³ + exp(1)"),
    ("(∂u/∂x)² + 1/sqrt(x)", "(∂u/∂x)² + 1/sqrt(x)"),
])
def test_the_notation_reads_back_as_written(written, printed):
    expr = Parser(SPACE).expression(written)
    assert round_trip(expr) == printed


def test_an_equation_is_written_as_lhs_equals_rhs():
    u = SPACE._fields["u"]
    text = equation_notation(SPACE.dt(u) - SPACE.div(SPACE.grad(u)), 1)
    assert text == "-∂(∂u/∂x)/∂x - ∂(∂u/∂y)/∂y + ∂u/∂t = 1"
    (residual,) = Parser(SPACE).equation(text)
    assert residual == SPACE.dt(u) - SPACE.div(SPACE.grad(u)) - 1


def test_what_the_notation_cannot_say_is_refused():
    with pytest.raises(NotationError, match="tensor"):
        to_notation(SPACE.grad(SPACE._fields["v"]))
    with pytest.raises(NotationError, match="no notation"):
        to_notation(sympy.atan(sympy.Symbol("x")))


def _space_of(doc):
    """The Space adr_problem parses a spec in (ymf.symbolic.problem)."""
    from ymf.symbolic.problem import _coefficients, _domain_dimension, box_bounds
    strong = doc["model"]["strong_form"]
    box = box_bounds(strong["domain"])
    dim = len(box[0]) if box else _domain_dimension(strong["domain"])
    ranks = {(u if isinstance(u, str) else u["name"]):
             (0 if isinstance(u, str) else int(u.get("rank", 0)))
             for u in strong["unknowns"]}
    return Space(dim, ranks, _coefficients(strong.get("coefficients") or {}, dim))


SPECS = sorted(p for p in EXAMPLES.rglob("*.ymf")
               if "ymf_state: abstract" not in p.read_text(encoding="utf-8"))


@pytest.mark.parametrize("path", SPECS, ids=lambda p: p.name)
def test_every_example_spec_round_trips(path):
    from ymf.compose import YmfCompositionError
    try:
        doc = load_composed(path)
    except YmfCompositionError as exc:
        if "does not exist" in str(exc):
            pytest.skip("extends an archive a run writes: %s" % (exc,))
        raise
    strong = doc["model"]["strong_form"]
    if not strong.get("equations"):
        pytest.skip("no equations in notation")
    space = _space_of(doc)
    count = 0
    for text in strong["equations"]:
        for residual in parse_equation(text, space):
            round_trip(residual, space)
            count += 1
    unknowns = list(space.unknowns)
    for solution in (doc.get("solution_paths") or {}).get("analytical") or []:
        formula = (solution.get("solution") or {}).get("formula")
        if formula:
            for value in Parser(space).assignments(formula, unknowns).values():
                round_trip(value, space)
                count += 1
    assert count


def _random_expression(rng, space, depth):
    x, y = space.x
    leaves = [x, y, space.t, sympy.pi, space._fields["u"], *space._fields["v"],
              sympy.Integer(rng.randint(1, 9)), sympy.Rational(rng.randint(1, 5), rng.randint(2, 7)),
              sympy.Float(rng.choice(["0.025", "1.5", "1e-10", "3.75"]))]
    if depth == 0:
        return rng.choice(leaves)
    a = _random_expression(rng, space, depth - 1)
    b = _random_expression(rng, space, depth - 1)
    choice = rng.randrange(10)
    if choice == 0:
        return a + b
    if choice == 1:
        return a - b
    if choice == 2:
        return a * b
    if choice == 3:
        return a / b if b != 0 else a
    if choice == 4:
        return a ** rng.choice([2, 3, 4, 5, -1, -2, sympy.Rational(1, 2), sympy.Rational(1, 3)])
    if choice == 5:
        return rng.choice([sympy.sin, sympy.cos, sympy.exp, sympy.tanh, sympy.Abs])(a)
    if choice == 6:
        return space.d(a, rng.randrange(2))
    if choice == 7:
        return space.dt(a)
    if choice == 8:
        return -a
    return a * space.d(b, rng.randrange(2))


def test_random_expressions_round_trip():
    rng = random.Random(20261018)
    for _ in range(400):
        expr = _random_expression(rng, SPACE, rng.randint(1, 4))
        if expr.has(sympy.zoo, sympy.nan, sympy.oo):
            continue
        round_trip(expr)
