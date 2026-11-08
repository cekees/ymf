"""Decimals under a power are kept as written (C5 finding Y2): sympy before
1.14 factors the largest coefficient out of a power of a sum with a decimal
in it ((x − 47.5)² -> 2256.25 (0.0210526 x − 1)²), which is the same number
by other arithmetic. dam3Bumps' bed is the case that showed it."""

import pytest

pytest.importorskip("sympy")

from ymf.symbolic.adr import to_code  # noqa: E402
from ymf.symbolic.language import Space, parse_expression  # noqa: E402


def _code(text):
    return to_code(parse_expression(text, Space(2, {}, {})))


@pytest.mark.parametrize("text, code", [
    ("(x − 47.5)²", "(x - 47.5)**2"),
    ("sqrt((x − 47.5)² + (y − 15)²)", "numpy.sqrt((x - 47.5)**2 + (y - 15)**2)"),
    ("3 − 3 sqrt((x − 47.5)² + (y − 15)²)/10", "3 - 3/10*numpy.sqrt((x - 47.5)**2 + (y - 15)**2)"),
    ("(2.5 x + 1)^(1/3)", "(2.5*x + 1)**(1/3)"),
    ("norm((x − 0.5, y))", "numpy.sqrt(y**2 + (x - 0.5)**2)"),
])
def test_a_power_of_a_sum_keeps_its_decimals(text, code):
    assert _code(text) == code


def test_integers_and_rationals_evaluate_as_before():
    assert _code("(x + 1)²") == "(x + 1)**2"
    assert _code("(x + 1)^0") == "1"
    assert _code("sqrt(4)") == "2"
