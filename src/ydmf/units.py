"""Unit handling, dimensional consistency checking, and non-dimensionalization.

Implements the "immediate"/"short-term" recommendations from
``docs/ydmf-units-report.md``:

- pint-backed unit string parsing (report §1.1, §6.2), tolerant of the
  bare-exponent (``m2/s``), unicode-superscript (``m²``), caret (``m^2``),
  and dot-multiply (``kg/(m·s)``) idioms used throughout the YDMF docs and
  examples, none of which pint accepts natively.
- Cross-field dimensional consistency checks (report §4.4): boundary/initial
  condition units must be compatible with the unknown they apply to.
- Characteristic-scale-based non-dimensionalization (report §5).

Deliberately **non-intrusive** (report §2.1): every function accepts YDMF
documents with no units at all and simply reports "nothing to check" rather
than raising. Units are an enhancement, never a requirement.

Symbolic parsing of ``strong_form_expression`` term-by-term (report §4.2) is
explicitly **out of scope** here — it needs a real expression parser and is
deferred to the ibvp/sympy language layer per the v0.2 delta doc §7. What
*is* implemented is everything checkable from unit strings already present
on the parsed document.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import pint

from ydmf.normalize import normalize_unknowns

__all__ = [
    "UnitParseError",
    "UnitIssue",
    "UnitCheckResult",
    "normalize_unit_string",
    "to_pint_unit",
    "dimensionality",
    "units_compatible",
    "check",
    "characteristic_scale_values",
    "compute_dimensionless_numbers",
    "non_dimensionalize",
]

_UREG = pint.UnitRegistry()

# Unicode superscript digits / minus sign seen in report examples (m², m³, ⁻¹).
_SUPERSCRIPT_MAP = str.maketrans(
    {"²": "2", "³": "3", "⁴": "4", "¹": "1", "⁰": "0", "⁻": "-"}
)

# Matches a unit letter (incl. µ/ohm-ish greek-ish chars pint knows) directly
# followed by digits (optionally negative, for superscript-minus cases like
# "s⁻¹" -> "s-1" after translation), e.g. "m2", "s3", "s-1" — but not
# something already turned into "m**2" by an earlier substitution.
_BARE_EXPONENT_RE = re.compile(r"(?<!\*)(?<=[A-Za-zµ])(-?\d+)(?!\d)")


class UnitParseError(ValueError):
    """Raised when a YDMF unit string cannot be parsed by pint."""


def normalize_unit_string(unit_str: Optional[str]) -> Optional[str]:
    """Normalize a YDMF-style unit string into pint-parseable syntax.

    Handles the notational idioms used throughout the YDMF docs/examples
    that pint's default parser rejects outright:

    >>> normalize_unit_string("m2/s")
    'm**2/s'
    >>> normalize_unit_string("kg/m3")
    'kg/m**3'
    >>> normalize_unit_string("kg/(m·s)")
    'kg/(m*s)'
    >>> normalize_unit_string("m^2")
    'm**2'
    """
    if unit_str is None:
        return None
    s = unit_str.strip()
    s = s.translate(_SUPERSCRIPT_MAP)
    s = s.replace("·", "*").replace("×", "*")
    s = s.replace("^", "**")
    s = _BARE_EXPONENT_RE.sub(r"**\1", s)
    return s


def to_pint_unit(unit_str: Optional[str]) -> "pint.Unit":
    """Parse a YDMF unit string into a :class:`pint.Unit`.

    ``None``, ``""``, and ``"dimensionless"`` all map to pint's dimensionless
    unit. Raises :class:`UnitParseError` on anything pint can't parse even
    after :func:`normalize_unit_string`.
    """
    if not unit_str or unit_str.strip().lower() == "dimensionless":
        return _UREG.dimensionless
    normalized = normalize_unit_string(unit_str)
    try:
        return _UREG.parse_expression(normalized).units
    except Exception as exc:  # pint raises several distinct exception types
        raise UnitParseError(
            f"could not parse unit string {unit_str!r} (normalized: {normalized!r}): {exc}"
        ) from exc


def dimensionality(unit_str: Optional[str]) -> str:
    """Return the dimensionality string for a unit string, e.g. ``'[length] / [time]'``."""
    return str(to_pint_unit(unit_str).dimensionality)


def units_compatible(a: Optional[str], b: Optional[str]) -> bool:
    """True if two unit strings share the same dimensionality (i.e. are convertible)."""
    return to_pint_unit(a).dimensionality == to_pint_unit(b).dimensionality


@dataclass
class UnitIssue:
    severity: str  # "error" | "warning"
    location: str
    message: str

    def __str__(self) -> str:
        return f"[{self.severity}] {self.location}: {self.message}"


@dataclass
class UnitCheckResult:
    issues: List[UnitIssue] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not any(i.severity == "error" for i in self.issues)

    def __bool__(self) -> bool:
        return self.ok

    def __iter__(self):
        return iter(self.issues)

    def __len__(self) -> int:
        return len(self.issues)


def _collect_unknown_units(strong_form: Dict[str, Any]) -> Dict[str, str]:
    units: Dict[str, str] = {}
    for u in normalize_unknowns(strong_form.get("unknowns") or []):
        if u.get("units"):
            units[u["name"]] = u["units"]
    return units


def check(doc: Dict[str, Any]) -> UnitCheckResult:
    """Check dimensional consistency of a parsed YDMF document.

    Non-intrusive per report §2.1: any unknown, coefficient, boundary
    condition, or initial condition with no ``units`` field is silently
    skipped rather than flagged as an error.

    Checks performed (report §4.4):

    - Every declared unit string parses (bad strings -> ``error``).
    - Each ``boundary_conditions[].units`` is dimensionally compatible with
      the ``unknowns[].units`` of the variable it references.
    - Each ``initial_conditions[].units`` is dimensionally compatible with
      the ``unknowns[].units`` of the field it initializes.

    Returns a :class:`UnitCheckResult`; ``bool(result)`` is ``True`` iff there
    are no ``error``-severity issues (warnings don't affect truthiness).
    """
    result = UnitCheckResult()
    problem = doc.get("Problem") or {}
    strong_form = problem.get("strong_form") or {}
    unknown_units = _collect_unknown_units(strong_form)

    for name, unit_str in unknown_units.items():
        try:
            to_pint_unit(unit_str)
        except UnitParseError as exc:
            result.issues.append(UnitIssue("error", f"unknowns[{name}]", str(exc)))

    def _check_against_unknown(entries, var_key, location_fmt):
        for entry in entries:
            var = entry.get(var_key)
            entry_units = entry.get("units")
            if not entry_units or var not in unknown_units:
                continue
            location = location_fmt(entry)
            try:
                if not units_compatible(entry_units, unknown_units[var]):
                    result.issues.append(
                        UnitIssue(
                            "error",
                            location,
                            f"units {entry_units!r} ({dimensionality(entry_units)}) incompatible "
                            f"with {var!r}'s declared units {unknown_units[var]!r} "
                            f"({dimensionality(unknown_units[var])})",
                        )
                    )
            except UnitParseError as exc:
                result.issues.append(UnitIssue("error", location, str(exc)))

    _check_against_unknown(
        strong_form.get("boundary_conditions") or [],
        "variable",
        lambda bc: f"boundary_conditions[region={bc.get('region')!r}, variable={bc.get('variable')!r}]",
    )
    _check_against_unknown(
        strong_form.get("initial_conditions") or [],
        "field",
        lambda ic: f"initial_conditions[field={ic.get('field')!r}]",
    )

    coefficients = strong_form.get("coefficients") or {}
    for name, spec in coefficients.items():
        if isinstance(spec, dict) and spec.get("units"):
            try:
                to_pint_unit(spec["units"])
            except UnitParseError as exc:
                result.issues.append(UnitIssue("error", f"coefficients[{name}]", str(exc)))

    return result


def _eval_formula(formula: str, scope: Dict[str, float]) -> Optional[float]:
    """Evaluate a characteristic-scale/dimensionless-number formula string.

    Formulas come from the YDMF document author (trusted, same trust level
    as any other YAML content in the file) but evaluation is still
    restricted to arithmetic on the given scope: no builtins, no attribute
    or subscript access. Returns ``None`` (rather than raising) if the
    formula references a name not present in ``scope`` or otherwise fails to
    evaluate, since not all referenced scales may be available yet.
    """
    try:
        code = compile(formula, "<ydmf-formula>", "eval")
    except SyntaxError:
        return None
    for name in code.co_names:
        if name not in scope:
            return None
    try:
        return eval(code, {"__builtins__": {}}, dict(scope))  # noqa: S307
    except Exception:
        return None


def characteristic_scale_values(problem: Dict[str, Any]) -> Dict[str, float]:
    """Resolve ``Problem.characteristic_scales`` to a flat ``{name: value}`` dict.

    Direct (``value:``) scales are taken as-is. ``derived: true`` scales
    (e.g. ``time: {derived: true, formula: "L_char / U_char"}``, report
    §5.2) are resolved by substituting the standard ``*_char`` aliases
    (``L_char``, ``U_char``, ``rho_char``, ``mu_char``, ``T_char``) for
    ``length``, ``velocity``, ``density``, ``viscosity``, ``temperature``
    respectively. Unresolvable derived scales are omitted, not errored —
    non-dimensionalization is best-effort per report §5.3.
    """
    scales = problem.get("characteristic_scales") or {}
    values: Dict[str, float] = {}
    derived: Dict[str, str] = {}

    for name, spec in scales.items():
        if not isinstance(spec, dict):
            continue
        if spec.get("value") is not None:
            values[name] = spec["value"]
        elif spec.get("derived") and spec.get("formula"):
            derived[name] = spec["formula"]

    alias_of = {
        "length": "L_char",
        "velocity": "U_char",
        "density": "rho_char",
        "viscosity": "mu_char",
        "temperature": "T_char",
    }

    # Fixed-point resolution: a couple of passes handle chained derived
    # scales (e.g. pressure derived from a derived velocity), though the
    # common case (report examples) resolves in one pass.
    for _ in range(len(derived) + 1):
        scope = {alias: values[key] for key, alias in alias_of.items() if key in values}
        progressed = False
        for name, formula in list(derived.items()):
            v = _eval_formula(formula, scope)
            if v is not None:
                values[name] = v
                del derived[name]
                progressed = True
        if not progressed:
            break

    return values


def compute_dimensionless_numbers(problem: Dict[str, Any]) -> Dict[str, float]:
    """Evaluate ``Problem.dimensionless_numbers[].formula`` (report §5.2).

    Formulas are evaluated against a scope built from
    :func:`characteristic_scale_values` (using bare names ``L``, ``U``,
    ``rho``, ``mu``, ``T`` as well as the ``*_char`` aliases, matching the
    variable names used in the report's example formulas, e.g.
    ``"rho * U * L / mu"``) plus a fixed gravitational constant ``g`` and any
    scalar-valued ``strong_form.coefficients`` entries. Numbers whose formula
    can't be fully resolved (missing scale/coefficient) are omitted.
    """
    scales = characteristic_scale_values(problem)
    bare_names = {
        "length": "L",
        "velocity": "U",
        "density": "rho",
        "viscosity": "mu",
        "temperature": "T",
    }
    scope: Dict[str, float] = {"g": 9.81}
    for key, bare in bare_names.items():
        if key in scales:
            scope[bare] = scales[key]
    scope.update(scales)  # also expose raw scale names, e.g. "length"

    strong_form = problem.get("strong_form") or {}
    for name, spec in (strong_form.get("coefficients") or {}).items():
        if isinstance(spec, dict) and isinstance(spec.get("value"), (int, float)):
            scope.setdefault(name, spec["value"])
        elif isinstance(spec, (int, float)):
            scope.setdefault(name, spec)

    numbers: Dict[str, float] = {}
    for name, spec in (problem.get("dimensionless_numbers") or {}).items():
        formula = spec.get("formula") if isinstance(spec, dict) else None
        if not formula:
            continue
        v = _eval_formula(formula, scope)
        if v is not None:
            numbers[name] = v
    return numbers


# Dimensionalities (as pint UnitsContainer, via a representative unit string)
# used to classify unknowns for automatic substitution-relation generation.
_VELOCITY_DIM = _UREG.get_dimensionality("m/s")
_LENGTH_DIM = _UREG.get_dimensionality("m")
_TIME_DIM = _UREG.get_dimensionality("s")
_PRESSURE_DIM = _UREG.get_dimensionality("Pa")


def non_dimensionalize(doc: Dict[str, Any]) -> Dict[str, Any]:
    """Build the non-dimensionalization block described in report §5.2/§5.3.

    Given a document with ``Problem.characteristic_scales`` and
    units-annotated ``unknowns``, returns:

    - ``characteristic_scales``: resolved numeric values (incl. derived ones)
    - ``dimensionless_numbers``: evaluated ``{name: value}``
    - ``substitutions``: ``{"<name>*": "<formula>"}`` for each unknown whose
      dimensionality matches a known pattern (velocity, length, time,
      pressure), mirroring the report's worked examples (§5.4)

    This produces the substitution *relations*, not a symbolic substitution
    into ``strong_form_expression`` — doing that requires a full sympy-based
    PDE parser and is explicitly deferred (v0.2 delta doc §7: "automatic
    non-dimensionalization pipeline").
    """
    problem = doc.get("Problem") or {}
    scales = characteristic_scale_values(problem)
    numbers = compute_dimensionless_numbers(problem)

    strong_form = problem.get("strong_form") or {}
    substitutions: Dict[str, str] = {}
    for u in normalize_unknowns(strong_form.get("unknowns") or []):
        name = u["name"]
        unit_str = u.get("units")
        if not unit_str:
            continue
        try:
            dim = to_pint_unit(unit_str).dimensionality
        except UnitParseError:
            continue
        if dim == _VELOCITY_DIM:
            substitutions[f"{name}*"] = f"{name} / U_char"
        elif dim == _PRESSURE_DIM:
            substitutions[f"{name}*"] = f"{name} / (rho_char * U_char**2)"
        elif dim == _LENGTH_DIM:
            substitutions[f"{name}*"] = f"{name} / L_char"
        elif dim == _TIME_DIM:
            substitutions[f"{name}*"] = f"{name} * U_char / L_char"

    return {
        "characteristic_scales": scales,
        "dimensionless_numbers": numbers,
        "substitutions": substitutions,
    }
