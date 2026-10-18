"""Write sympy back in the notation of a spec: the inverse of :mod:`.notation`.

A tool that derives something from a spec -- ``pro weak`` writing the weak
form it will solve, say -- writes it in the notation a person reads, so the
record says what was solved in the same terms as the spec says what was
posed. The contract is a round trip: for an expression ``e`` built by
:class:`~ymf.symbolic.notation.Parser` (or by the
:class:`~ymf.symbolic.language.Space` operators it applies),
``Parser(space).expression(to_notation(e)) == e``, structurally, not just
mathematically.

What can be written is what the notation can read:

* numbers: integers, rationals (``1/2``), floats in their shortest exact
  form (``0.025``, ``1e-10``); π;
* names: the coordinates, t, the unknowns and their components
  (``u``, ``v_0``), and any other symbol by its name;
* ``+ - * /`` and powers, with ``² ³ ⁴`` for small integer exponents and
  ``^`` otherwise; ``sqrt`` for a square root;
* the functions sin cos tan exp log sqrt sinh cosh tanh abs;
* derivatives: ``D(e, axis)`` as ``∂e/∂x``, ``Dt(e)`` as ``∂e/∂t``;
* a column vector as a tuple ``(a, b)``.

Anything else (a tensor, an integral, a function the notation lacks, or
what sympy forms only from complex values -- ``I``, ``sign`` -- out of a
root of a negative quantity) is refused with :class:`NotationError` rather
than written in a form the parser would read differently.

Choices that keep the round trip exact: factors are separated by spaces,
so names never run together into one (``x y``, not ``xy``); a compound
denominator is parenthesized, since ``/`` binds left to right
(``x/(2 y)``); and a derivative raised to a power is parenthesized,
``(∂u/∂x)²``, which also reads unambiguously.
"""

from __future__ import annotations

from typing import List

import sympy

from ymf.symbolic.language import D, Dt
from ymf.symbolic.notation import NotationError

__all__ = ["to_notation", "equation_notation"]

_SUPERSCRIPTS = {2: "²", 3: "³", 4: "⁴"}
_AXES = "xyz"
_FUNCTIONS = {
    sympy.sin: "sin", sympy.cos: "cos", sympy.tan: "tan", sympy.exp: "exp",
    sympy.log: "log", sympy.sinh: "sinh", sympy.cosh: "cosh", sympy.tanh: "tanh",
    sympy.Abs: "abs",
}

# precedence of what a printed form is, for deciding on parentheses
_ATOM, _POWER, _PRODUCT, _SUM = 4, 3, 2, 1


def to_notation(expr) -> str:
    """``expr`` in the notation of a spec (see the module docstring)."""
    if isinstance(expr, sympy.MatrixBase):
        if expr.shape[1] != 1:
            raise NotationError("a %dx%d tensor has no notation; only vectors (a, b) do"
                                % expr.shape)
        return "(" + ", ".join(_text(e)[0] for e in expr) + ")"
    return _text(sympy.sympify(expr))[0]


def equation_notation(lhs, rhs=0) -> str:
    """``lhs = rhs`` in the notation."""
    return "%s = %s" % (to_notation(lhs), to_notation(rhs))


def _wrap(text: str, level: int, need: int) -> str:
    return "(%s)" % text if level < need else text


def _number(n) -> tuple:
    """A non-negative number's text and level (the sign is the caller's)."""
    if n.is_Integer:
        return str(int(n)), _ATOM
    if n.is_Rational:
        return "%d/%d" % (n.p, n.q), _PRODUCT
    if n.is_Float:
        return repr(float(n)), _ATOM
    raise NotationError("the number %s has no notation" % (n,))


def _text(e) -> tuple:
    """(text, level) for a scalar expression."""
    if e is sympy.pi:
        return "π", _ATOM
    if e is sympy.E:
        return "exp(1)", _ATOM
    if e.is_Number:
        if e < 0:
            inner, level = _number(-e)
            return "-" + _wrap(inner, level, _PRODUCT), _SUM
        return _number(e)
    if e.is_Symbol:
        return e.name, _ATOM
    if isinstance(e, D):
        operand, axis = e.args
        if int(axis) not in range(3):
            raise NotationError("D along axis %s has no notation" % (axis,))
        return "∂%s/∂%s" % (_operand(operand), _AXES[int(axis)]), _PRODUCT
    if isinstance(e, Dt):
        return "∂%s/∂t" % (_operand(e.args[0]),), _PRODUCT
    if e.is_Add:
        return _sum(e), _SUM
    if e.is_Mul:
        return _product(e)
    if e.is_Pow:
        return _power(e)
    for function, name in _FUNCTIONS.items():
        if isinstance(e, function):
            return "%s(%s)" % (name, _text(e.args[0])[0]), _ATOM
    raise NotationError("%s (%s) has no notation" % (e, type(e).__name__))


def _operand(e) -> str:
    """What follows ∂: a name bare, anything else parenthesized."""
    text, level = _text(e)
    return text if e.is_Symbol else "(%s)" % text


def _sum(e) -> str:
    out: List[str] = []
    for term in e.as_ordered_terms():
        text, _ = _text(term)
        if out:
            out.append("- " + text[1:] if text.startswith("-") else "+ " + text)
        else:
            out.append(text)
    return " ".join(out)


def _product(e) -> tuple:
    coeff, rest = e.as_coeff_Mul()
    negative = coeff < 0
    coeff = -coeff if negative else coeff
    factors = sympy.Mul.make_args(rest)
    if (negative or coeff != 1) and len(factors) > 1 and any(f.is_Add for f in factors):
        # sympy distributes a number over a sum only in a two-factor
        # product, so 2*(a + b)*c keeps its sum; read left to right,
        # "2 (a + b) c" would multiply 2 into the sum first. The number
        # multiplies the whole of the rest instead: -2 ((a + b) c)/3.
        body, _ = _product(rest)
        text = "(%s)" % body
        if coeff.is_Rational:
            if coeff.p != 1:
                text = "%d %s" % (coeff.p, text)
            if coeff.q != 1:
                text += "/%d" % coeff.q
        elif coeff != 1:
            text = "%s %s" % (_number(coeff)[0], text)
        return ("-" + text, _SUM) if negative else (text, _PRODUCT)
    numerator, denominator = [], []
    if coeff.is_Rational and not coeff.is_Integer:
        if coeff.p != 1:
            numerator.append(str(coeff.p))
        denominator.append(sympy.Integer(coeff.q))
    elif coeff != 1:
        numerator.append(_number(coeff)[0])
    for factor in sympy.Mul.make_args(rest):
        base, exponent = factor.as_base_exp()
        if exponent.is_Number and exponent < 0:
            denominator.append(base ** -exponent)
        else:
            numerator.append(_wrap(*_text(factor), _PRODUCT + 1)
                             if not isinstance(factor, (D, Dt)) else _text(factor)[0])
    text = " ".join(numerator) or "1"
    if denominator and any(f.is_Add or (f.is_Pow and f.base.is_Add) for f in denominator):
        # multiplying the denominators together would distribute a number
        # over a sum (2*(a + b) -> 2a + 2b); dividing by each in turn keeps
        # them apart: 9/2/(a + b), not 9/(2 (a + b))
        for factor in denominator:
            text += "/" + _wrap(*_text(factor), _ATOM)
    elif denominator:
        below = sympy.Mul(*denominator)
        below_text, below_level = _text(below)
        text += "/" + _wrap(below_text, below_level, _ATOM)
    if negative:
        return "-" + text, _SUM
    return text, _PRODUCT


def _power(e) -> tuple:
    base, exponent = e.args
    if exponent == sympy.Rational(1, 2):
        return "sqrt(%s)" % (_text(base)[0],), _ATOM
    if exponent.is_Number and exponent < 0:
        return _product(e)                    # 1/(...)
    base_text, base_level = _text(base)
    # a derivative or a negative number as a base reads ambiguously bare
    if isinstance(base, (D, Dt)) or base_level < _ATOM:
        base_text = "(%s)" % base_text
    if exponent.is_Integer and int(exponent) in _SUPERSCRIPTS:
        return base_text + _SUPERSCRIPTS[int(exponent)], _POWER
    exp_text, exp_level = _text(exponent)
    return "%s^%s" % (base_text, _wrap(exp_text, exp_level, _ATOM)), _POWER
