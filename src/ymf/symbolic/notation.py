"""The written notation of a YMF spec: Unicode mathematics, parsed.

A spec states each fact once, in the notation a person writes and reads:

    equations:
      - "∂T/∂t = ∇·(κ∇T)  in Ω"
      - "ρ ∂v/∂t + ∇·(ρ v⊗v) - μΔv + ∇p = f"
    domain: "Ω = [0, 1] × [0, 1],  t ∈ (0, 50]"
    boundary_regions:
      - {name: walls, geometry: "y = 0 or y = H"}
    coefficients:
      f: "π² sin(πx) sin(πy)"

This module turns those strings into sympy, through the operators of
:class:`ymf.symbolic.language.Space`. It is a small recursive-descent
parser, not a translation to Python, so the grammar is exactly what is
written here and anything outside it is an error that says where.

Expressions
    numbers ``2``, ``1.5``, ``1.0e-3``; names (coefficients, unknowns, the
    coordinates x y z, time t, π); ``+ - * / ^ **`` and the superscripts
    ``² ³ ⁴``; implicit multiplication (``2π²κt``, ``κ ΔT``, ``sin(πx)``);
    ``·`` (a dot product of vectors, a product otherwise) and ``⊗``
    (outer product, binding tighter than ``·`` and ``*``); tuples
    ``(a, b)`` are vectors; functions sin cos tan exp log sqrt sinh cosh
    tanh abs.
Differential operators, each applied to the operand that follows it
    ``∇u`` gradient, ``∇·F`` divergence, ``Δu`` Laplacian,
    ``∂u/∂t`` and ``∂u/∂x`` partial derivatives; also the spelled-out
    ``grad div lap dt dx dy dz dot outer``.
A run of letters that is not a known name is read as a product of
single-letter names (``πx`` is π·x) when every letter is known, and is an
error otherwise -- so a misspelt coefficient is reported, not multiplied.

Equations are ``lhs = rhs``, optionally followed by ``in <region>``, which
is ignored. Predicates (boundary regions) are ``a = b`` tests joined by
``and`` and ``or``. Domains are ``Ω = I × I [× I]`` with intervals
``[a, b]``, optionally followed by ``, t ∈ (t0, t1]``. Assignments
(analytical solutions) are ``name = expr`` lines; ``T(x, y, t) = ...`` is
accepted, and a name that is not an unknown becomes a local definition
usable on later lines (``λ = Re/2 - sqrt(Re²/4 + 4π²)``).
"""

from __future__ import annotations

import re
import unicodedata
from typing import Any, Dict, List, Optional, Tuple

import sympy

__all__ = ["NotationError", "Parser", "tokens"]


class NotationError(ValueError):
    """A string in a spec does not follow the notation."""


_SUPERSCRIPTS = {"²": 2, "³": 3, "⁴": 4}
_ALIASES = {"−": "-", "⋅": "·", "∙": "·", "–": "-", "∆": "Δ"}
_OPERATORS = ["**", "^", "+", "-", "*", "/", "(", ")", "[", "]", ",", "=",
              "·", "⊗", "∇", "Δ", "∂", "×", "∈", ";"] + list(_SUPERSCRIPTS)
_NUMBER = re.compile(r"\d+(\.\d*)?([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?")


def _is_letter(ch: str) -> bool:
    # Δ is a Greek capital letter to Unicode, but here it is the Laplacian
    return ch != "Δ" and (ch == "_" or unicodedata.category(ch).startswith("L"))


def tokens(text: str) -> List[Tuple[str, str, int]]:
    """(kind, value, position) for each token; kinds NUM, NAME, OP, NL."""
    out = []
    i = 0
    while i < len(text):
        ch = _ALIASES.get(text[i], text[i])
        if ch == "\n":
            out.append(("NL", "\n", i)); i += 1; continue
        if ch.isspace():
            i += 1; continue
        m = _NUMBER.match(text, i)
        if m and (ch in "0123456789" or ch == "."):
            out.append(("NUM", m.group(0), i)); i = m.end(); continue
        if _is_letter(ch):
            j = i + 1
            # ASCII digits only: '²'.isdigit() is True, and ² is an operator
            while j < len(text) and (_is_letter(text[j]) or text[j] in "0123456789"):
                j += 1
            out.append(("NAME", text[i:j], i)); i = j; continue
        for op in _OPERATORS:
            if text.startswith(op, i) or ch == op:
                out.append(("OP", op if text.startswith(op, i) else ch, i))
                i += len(op) if text.startswith(op, i) else 1
                break
        else:
            raise NotationError("unexpected character %r at %d in %r" % (text[i], i, text))
    return out


#: words that join predicates; never names, never implicit factors
_KEYWORDS = {"and", "or", "in"}

_FUNCTIONS = {
    "sin": sympy.sin, "cos": sympy.cos, "tan": sympy.tan, "exp": sympy.exp,
    "log": sympy.log, "sqrt": sympy.sqrt, "sinh": sympy.sinh, "cosh": sympy.cosh,
    "tanh": sympy.tanh, "abs": sympy.Abs,
}


class Parser:
    """Parse notation against a :class:`~ymf.symbolic.language.Space`."""

    def __init__(self, space, locals_: Optional[Dict[str, Any]] = None):
        self.space = space
        self.names: Dict[str, Any] = {"π": sympy.pi, "pi": sympy.pi, "t": space.t}
        for axis, n in enumerate("xyz"[:space.dim]):
            self.names[n] = space.x[axis]
        self.names.update(space.coefficients)
        self.names.update(space._fields)
        self.names.update(locals_ or {})
        self.functions = dict(_FUNCTIONS)
        self.functions.update({
            "div": space.div, "grad": space.grad, "lap": lambda e: space.div(space.grad(e)),
            "dt": space.dt, "dot": self._dot, "outer": space.outer,
            "vec": lambda *xs: sympy.Matrix([sympy.sympify(v) for v in xs]),
        })
        for axis, n in enumerate("xyz"[:space.dim]):
            self.functions["d" + n] = (lambda e, a=axis: space.d(e, a))

    # -- driving -------------------------------------------------------------

    def _start(self, text: str):
        self.text = text
        self.toks = [t for t in tokens(text) if t[0] != "NL"]
        self.pos = 0

    def _peek(self, k: int = 0):
        j = self.pos + k
        return self.toks[j] if j < len(self.toks) else ("END", "", len(self.text))

    def _take(self):
        tok = self._peek(); self.pos += 1; return tok

    def _expect(self, value: str):
        tok = self._take()
        if tok[1] != value:
            self._fail("expected %r" % value, tok)
        return tok

    def _fail(self, message: str, tok=None):
        tok = tok or self._peek()
        where = tok[2]
        raise NotationError("%s at position %d in %r (%s)"
                            % (message, where, self.text,
                               "end of text" if tok[0] == "END" else "near %r" % tok[1]))

    def _done(self):
        if self._peek()[0] != "END":
            self._fail("unexpected text")

    # -- grammar -------------------------------------------------------------

    def expression(self, text: str):
        self._start(text)
        value = self._sum()
        self._done()
        return value

    def _sum(self):
        value = self._product()
        while self._peek()[1] in ("+", "-"):
            op = self._take()[1]
            rhs = self._product()
            value = value + rhs if op == "+" else value - rhs
        return value

    def _starts_operand(self, tok) -> bool:
        if tok[0] == "NAME" and tok[1] in _KEYWORDS:
            return False
        return tok[0] in ("NUM", "NAME") or tok[1] in ("(", "∇", "Δ", "∂")

    def _product(self):
        value = self._unary()
        while True:
            tok = self._peek()
            if tok[1] in ("*", "/", "·"):
                self._take()
                rhs = self._unary()
                if tok[1] == "*":
                    value = self._mul(value, rhs)
                elif tok[1] == "/":
                    value = value / rhs if not isinstance(rhs, sympy.MatrixBase) else self._fail(
                        "division by a vector", tok)
                else:
                    value = self._dot(value, rhs)
            elif self._starts_operand(tok):
                value = self._mul(value, self._unary())      # implicit
            else:
                return value

    def _unary(self):
        if self._peek()[1] == "-":
            self._take(); return -self._unary()
        if self._peek()[1] == "+":
            self._take(); return self._unary()
        return self._outer()

    def _outer(self):
        value = self._power()
        while self._peek()[1] == "⊗":
            self._take()
            value = self.space.outer(value, self._power())
        return value

    def _power(self):
        base = self._primary()
        while True:
            tok = self._peek()
            if tok[1] in _SUPERSCRIPTS:
                self._take(); base = base ** _SUPERSCRIPTS[tok[1]]
            elif tok[1] in ("^", "**"):
                self._take(); base = base ** self._unary()
            else:
                return base

    def _primary(self):
        tok = self._take()
        kind, value = tok[0], tok[1]
        if kind == "NUM":
            return sympy.Integer(value) if re.fullmatch(r"\d+", value) else sympy.Float(value)
        if kind == "NAME":
            return self._name(tok)
        if value == "(":
            items = [self._sum()]
            while self._peek()[1] == ",":
                self._take(); items.append(self._sum())
            self._expect(")")
            return items[0] if len(items) == 1 else sympy.Matrix([sympy.sympify(v) for v in items])
        if value == "∇":
            if self._peek()[1] == "·":
                self._take(); return self.space.div(self._power())
            return self.space.grad(self._power())
        if value == "Δ":
            return self.space.div(self.space.grad(self._power()))
        if value == "∂":
            operand = self._power()
            self._expect("/")
            self._expect("∂")
            var = self._take()
            if var[1] == "t":
                return self.space.dt(operand)
            if var[1] in "xyz"[:self.space.dim] and len(var[1]) == 1:
                axis = "xyz".index(var[1])
                if isinstance(operand, sympy.MatrixBase):
                    return operand.applyfunc(lambda e: self.space.d(e, axis))
                return self.space.d(operand, axis)
            self._fail("∂/∂%s: only t, x, y, z" % var[1], var)
        self._fail("expected a number, a name or '('", tok)

    def _name(self, tok):
        name = tok[1]
        if self._peek()[1] == "(" and name in self.functions:
            self._take()
            args = []
            if self._peek()[1] != ")":
                args.append(self._sum())
                while self._peek()[1] == ",":
                    self._take(); args.append(self._sum())
            self._expect(")")
            try:
                return self.functions[name](*args)
            except (TypeError, ValueError) as exc:
                self._fail("%s(...): %s" % (name, exc), tok)
        if name in self.names:
            return self.names[name]
        if name in self.functions:
            self._fail("%s needs (arguments)" % name, tok)
        # a run of letters: a product of single-letter names, if all known
        parts = list(name)
        if all(p in self.names for p in parts):
            value = sympy.Integer(1)
            for p in parts:
                value = self._mul(value, self.names[p])
            return value
        unknown = name if len(name) == 1 else name
        self._fail("unknown name %r (not a coefficient, an unknown, a coordinate, t, "
                   "π or a function%s)" % (unknown, "" if len(name) == 1 else
                                            ", nor a product of single-letter names"), tok)

    # -- helpers -------------------------------------------------------------

    @staticmethod
    def _mul(a, b):
        if isinstance(a, sympy.MatrixBase) and isinstance(b, sympy.MatrixBase):
            raise NotationError("cannot multiply two vectors; use · (dot) or ⊗ (outer)")
        return a * b

    def _dot(self, a, b):
        if isinstance(a, sympy.MatrixBase) and isinstance(b, sympy.MatrixBase):
            return self.space.dot(a, b)
        return a * b

    # -- the larger forms ----------------------------------------------------

    def equation(self, text: str) -> List[sympy.Expr]:
        """``lhs = rhs [in region]`` -> scalar residuals lhs - rhs."""
        body = re.split(r"\s+in\s+", text, maxsplit=1)[0]
        self._start(body)
        left = self._sum()
        if self._peek()[1] == "=":
            self._take()
            right = self._sum()
        else:
            right = sympy.Integer(0)
        self._done()
        if isinstance(left, sympy.MatrixBase) and not isinstance(right, sympy.MatrixBase):
            if sympy.sympify(right) != 0:
                raise NotationError("%r sets a vector equal to a scalar" % text)
            right = sympy.zeros(*left.shape)
        elif isinstance(right, sympy.MatrixBase) and not isinstance(left, sympy.MatrixBase):
            raise NotationError("%r sets a scalar equal to a vector" % text)
        residual = left - right
        if isinstance(residual, sympy.MatrixBase):
            if residual.shape[1] != 1:
                raise NotationError("%r is a tensor equation" % text)
            return [sympy.sympify(r) for r in residual]
        return [sympy.sympify(residual)]

    def predicate(self, text: str) -> Optional[str]:
        """'x = 0 or y = H' -> numpy-free code using ``tol``; None for the whole boundary."""
        stripped = text.strip()
        if stripped in ("∂Ω", "∂ Ω", "boundary"):
            return None
        from ymf.symbolic.adr import to_code
        self._start(stripped)
        alternatives = []
        while True:
            tests = []
            while True:
                lhs = self._sum()
                self._expect("=")
                rhs = self._sum()
                tests.append("abs(%s) <= tol" % to_code(lhs - rhs))
                if self._peek()[1] == "and":
                    self._take(); continue
                break
            alternatives.append("(" + " and ".join(tests) + ")")
            if self._peek()[1] == "or":
                self._take(); continue
            break
        self._done()
        return " or ".join(alternatives)

    def _interval(self):
        open_ = self._take()
        if open_[1] not in ("[", "("):
            self._fail("expected an interval [a, b]", open_)
        a = self._sum()
        self._expect(",")
        b = self._sum()
        close = self._take()
        if close[1] not in ("]", ")"):
            self._fail("expected ] or ) to close the interval", close)
        return float(a), float(b)

    def domain(self, text: str):
        """'Ω = [0,1] × [0,1], t ∈ (0, 50]' -> (lower, upper, (t0, t1) or None)."""
        self._start(text)
        if self._peek()[0] == "NAME" and self._peek(1)[1] == "=":
            self._take(); self._take()
        boxes = [self._interval()]
        while self._peek()[1] == "×":
            self._take(); boxes.append(self._interval())
        time = None
        while self._peek()[1] == ",":
            self._take()
            name = self._take()
            self._expect("∈")
            if name[1] != "t":
                self._fail("only t ∈ (t0, t1] may follow the box", name)
            time = self._interval()
        self._done()
        lower = [a for a, _ in boxes]
        upper = [b for _, b in boxes]
        for a, b in boxes:
            if not b > a:
                raise NotationError("empty interval [%g, %g] in %r" % (a, b, text))
        return lower, upper, time

    def assignments(self, text: str, unknowns) -> Dict[str, Any]:
        """``name = expr`` lines -> {unknown: value}; other names become locals."""
        out: Dict[str, Any] = {}
        for line in re.split(r"[\n;]", text):
            line = line.strip()
            if not line:
                continue
            self._start(line)
            name = self._take()
            if name[0] != "NAME":
                self._fail("expected 'name = ...'", name)
            if self._peek()[1] == "(":            # T(x, y, t) = ...
                depth = 0
                while True:
                    tok = self._take()
                    depth += tok[1] == "("
                    depth -= tok[1] == ")"
                    if depth == 0:
                        break
            self._expect("=")
            value = self._sum()
            self._done()
            if name[1] in unknowns:
                out[name[1]] = value
            else:
                self.names[name[1]] = value
        return out
