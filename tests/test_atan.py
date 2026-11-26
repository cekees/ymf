"""``atan`` in the notation (C42 round 7): santos_step's step, 0.1 (0.5 + arctan(x/0.075)/π), can be stated.

It parses to sympy's ``atan``, prints back as ``atan`` (the printer's round trip), is ``numpy.arctan`` in the code a
solver evaluates, and differentiates (the bed slope a shallow-water kernel reads is ∂Z/∂x = 0.1/(0.075 π (1 +
(x/0.075)²))).
"""

import math

import pytest

sympy = pytest.importorskip("sympy")
pytest.importorskip("strictyaml")

from _formulation import SHALLOW_WATER, roundtrip  # noqa: E402
from ymf import closure  # noqa: E402
from ymf.names import FUNCTION_NAMES  # noqa: E402
from ymf.references import resolve  # noqa: E402
from ymf.symbolic import split_problem  # noqa: E402
from ymf.symbolic.adr import to_code  # noqa: E402
from ymf.symbolic.language import Space  # noqa: E402
from ymf.symbolic.notation import Parser  # noqa: E402
from ymf.symbolic.printer import to_notation  # noqa: E402

SANTOS = SHALLOW_WATER / "santos_step.ymf"
STEP = "0.1 (0.5 + atan(x/0.075)/π)"


def _code_at(code, x):
    import numpy
    return eval(code, {"__builtins__": {}, "numpy": numpy}, {"x": x})


def test_atan_is_a_function_of_the_notation():
    assert "atan" in FUNCTION_NAMES
    space = Space(2, {})
    x = space.x[0]
    assert Parser(space).expression("atan(2x)") == sympy.atan(2 * x)


def test_atan_prints_back_and_reads_the_same():
    space = Space(2, {})
    parser = Parser(space)
    for text in ("atan(x)", STEP, "atan(x/y)²", "sin(atan(x))"):
        e = parser.expression(text)
        assert parser.expression(to_notation(e)) == e, text
    assert to_notation(sympy.atan(space.x[0])) == "atan(x)"


@pytest.mark.parametrize("x", [-15.0, -0.075, 0.0, 0.03, 15.0])
def test_the_step_is_the_decks_in_code(x):
    code = to_code(Parser(Space(2, {})).expression(STEP))
    assert "numpy.arctan" in code
    deck = 0.1 * (0.5 + 1.0 / math.pi * math.atan(x / 0.075))
    assert _code_at(code, x) == pytest.approx(deck, rel=1e-15, abs=1e-17)


def test_atan_differentiates():
    space = Space(2, {})
    x = space.x[0]
    step = Parser(space).expression(STEP)
    slope = sympy.diff(step, x)
    assert sympy.simplify(slope - 0.1 / (0.075 * sympy.pi * (1 + (x / 0.075) ** 2))) == 0


def test_santos_step_can_be_stated():
    spec, _, _ = closure.load(SANTOS)
    resolve(spec)
    split = split_problem(spec, "relaxed_convex_limited")
    (gn,) = split["subproblems"]
    # the step reaches the depth's initial value and the bed slope in the momentum equation
    assert "numpy.arctan" in split["initial"]["h"]
    reaction = next(e for e in gn["adr"]["equations"] if e["component"] == "q_0")["reaction"]["r"]
    assert "numpy.pi*(177.7777777777778*x**2 + 1)" in reaction        # ∂Z/∂x = (0.1/0.075/π)/(1 + (x/0.075)²)
    # the wave over still water: the depth's initial value at the crest (x = xs = -3, Z there ~ 0.0008)
    z = 0.1 * (0.5 + math.atan(-3 / 0.075) / math.pi)
    assert _code_at(split["initial"]["h"], -3.0) == pytest.approx(0.2 + 0.0365 - z, rel=1e-14)
    hw = next(a for a in gn["dispersion"]["auxiliaries"] if a["role"] == "relaxed_divergence")
    assert not hw["initial_is_default"] and "tanh" in hw["initial"]


def test_santos_step_round_trips_through_an_archive(tmp_path):
    spec, again, digest = roundtrip(tmp_path, SANTOS, "relaxed_convex_limited", {"cells": "80x8"})
    assert again["model"]["strong_form"]["coefficients"]["Z"]["formula"] == STEP
    assert split_problem(again, "relaxed_convex_limited")["initial"] == split_problem(
        spec, "relaxed_convex_limited")["initial"]
