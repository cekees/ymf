"""The names a formula uses, read without sympy.

The notation (``ymf.symbolic.notation``) parses a formula into sympy; this
module holds the part of it that needs no sympy -- the tokenizer, and the
rule by which a run of letters is a name -- so that the core (closure,
references) can say which names a formula uses in an environment without
the ``symbolic`` extra, and get the same answer the parser would.

The rule, as ``Parser._name`` applies it: a run of letters followed by
``(`` that is a function of the notation is that function; a declared name
(a coefficient, an unknown, a vector unknown's component ``v_0``, a
coordinate, ``t``, ``π``) is that name; any other run of letters is the
product of its letters, each of which must be a name.
"""

from __future__ import annotations

import re
import unicodedata
from typing import Dict, Iterable, List, Set, Tuple

__all__ = ["NotationError", "tokens", "FUNCTION_NAMES", "names_in", "state_dependence"]


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
            # The `else` of the `for` above: no operator matched.
            if _SUBSCRIPTS_START <= ch <= _SUBSCRIPTS_END:
                # The superscripts ² ³ ⁴ are powers, so a subscript looks as
                # if it ought to work too. It does not, and the fix is to
                # spell the index into the name and declare it.
                run = ch
                while i + len(run) < len(text) and \
                        _SUBSCRIPTS_START <= text[i + len(run)] <= _SUBSCRIPTS_END:
                    run += text[i + len(run)]
                stem = out[-1][1] if out and out[-1][0] == "NAME" else ""
                plain = "".join(_SUBSCRIPT_PLAIN.get(c, c) for c in run)
                raise NotationError(
                    "%r at %d in %r is a subscript, and the notation has none; "
                    "write %r and declare it (a coefficient or an unknown)"
                    % (run, i, text, stem + plain))
            raise NotationError("unexpected character %r at %d in %r" % (text[i], i, text))
    return out


#: Unicode subscripts, which the notation does not have. Reported on their
#: own rather than as "unexpected character", because the superscripts
#: ² ³ ⁴ *are* notation (powers), so a subscript is a reasonable guess.
_SUBSCRIPTS_START, _SUBSCRIPTS_END = "\u2080", "\u209c"
_SUBSCRIPT_PLAIN = {chr(0x2080 + d): str(d) for d in range(10)}
_SUBSCRIPT_PLAIN.update({
    "\u2090": "a", "\u2091": "e", "\u2092": "o", "\u2093": "x",
    "\u2095": "h", "\u2096": "k", "\u2097": "l", "\u2098": "m",
    "\u2099": "n", "\u209a": "p", "\u209b": "s", "\u209c": "t",
})


#: the notation's function names (``ymf.symbolic.notation``'s ``_FUNCTIONS``
#: and the parser's own; a test checks the two agree)
FUNCTION_NAMES = frozenset([
    "sin", "cos", "tan", "exp", "log", "sqrt", "sinh", "cosh", "tanh", "abs", "hypot",
    "heaviside", "sign", "min", "max",
    "div", "grad", "lap", "dt", "dot", "outer", "vec", "norm", "sym", "dx", "dy", "dz",
])

#: names every formula may use
BUILTIN_NAMES = frozenset(["x", "y", "z", "t", "π", "pi"])


def names_in(text: str, declared: Iterable[str]) -> Set[str]:
    """The declared names ``text`` uses, as the parser reads it.

    A function call is not a name; a run of letters that is not declared
    counts for each of its letters (``ρU`` is ``ρ·U``). Names that are not
    declared are left out: a formula that uses one is an error the parser
    reports, not this.
    """
    declared = set(declared)
    out: Set[str] = set()
    toks = [t for t in tokens(str(text)) if t[0] != "NL"]
    for i, (kind, value, _) in enumerate(toks):
        if kind != "NAME":
            continue
        following = toks[i + 1][1] if i + 1 < len(toks) else ""
        if following == "(" and value in FUNCTION_NAMES:
            continue
        if value in declared:
            out.add(value)
        elif value not in FUNCTION_NAMES:
            out.update(ch for ch in value if ch in declared)
    return out


def _component_of(name: str, unknowns: Iterable[str]) -> str:
    stem, _, index = name.rpartition("_")
    return stem if stem in unknowns and index.isdigit() else name


def state_dependence(coefficients: Dict[str, object], unknowns: Iterable[str]
                     ) -> Dict[str, List[str]]:
    """``{coefficient: [unknown, ...]}`` for each coefficient whose formula
    names an unknown (or a component of one), directly or through another
    coefficient; sorted, in the order the unknowns are declared.

    Read from the formulas' text alone (no sympy). A coefficient with a
    ``value`` is a number, whatever else it has.
    """
    unknowns = list(unknowns)
    formulas = {}
    for name, entry in (coefficients or {}).items():
        if isinstance(entry, dict):
            if entry.get("value") is None and entry.get("formula") is not None:
                formulas[name] = str(entry["formula"])
        elif isinstance(entry, str):
            formulas[name] = entry
    components = {"%s_%d" % (u, i) for u in unknowns for i in range(3)}
    declared = set(coefficients or {}) | set(unknowns) | components | BUILTIN_NAMES
    direct = {name: {_component_of(n, unknowns) for n in names_in(text, declared)}
              for name, text in formulas.items()}
    out: Dict[str, List[str]] = {}
    for name in formulas:
        seen, todo, state = set(), [name], set()
        while todo:
            current = todo.pop()
            if current in seen:
                continue
            seen.add(current)
            for used in direct.get(current, ()):
                if used in unknowns:
                    state.add(used)
                elif used in direct:
                    todo.append(used)
        if state:
            out[name] = [u for u in unknowns if u in state]
    return out
