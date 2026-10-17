"""Unit handling, dimensional consistency checking, and non-dimensionalization.

Implements the "immediate"/"short-term" recommendations from
``docs/ymf-units-report.md``:

- pint-backed unit string parsing (report §1.1, §6.2), tolerant of the
  bare-exponent (``m2/s``), unicode-superscript (``m²``), caret (``m^2``),
  and dot-multiply (``kg/(m·s)``) idioms used throughout the YMF docs and
  examples, none of which pint accepts natively.
- Cross-field dimensional consistency checks (report §4.4): boundary/initial
  condition units must be compatible with the unknown they apply to.
- Characteristic-scale-based non-dimensionalization (report §5).

Deliberately **non-intrusive** (report §2.1): every function accepts YMF
documents with no units at all and simply reports "nothing to check" rather
than raising. Units are an enhancement, never a requirement.

Symbolic parsing of ``strong_form_expression`` term-by-term (report §4.2) is
explicitly **out of scope** here — it needs a real expression parser and is
deferred to the ibvp/sympy language layer per the v0.2 delta doc §7. What
*is* implemented is everything checkable from unit strings already present
on the parsed document.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import pint

from ymf.normalize import normalize_unknowns

__all__ = [
    "UnitParseError",
    "UnitIssue",
    "UnitCheckResult",
    "normalize_unit_string",
    "to_pint_unit",
    "dimensionality",
    "units_compatible",
    "check",
    "MISSING_UNITS",
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
    """Raised when a YMF unit string cannot be parsed by pint."""


def normalize_unit_string(unit_str: Optional[str]) -> Optional[str]:
    """Normalize a YMF-style unit string into pint-parseable syntax.

    Handles the notational idioms used throughout the YMF docs/examples
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
    """Parse a YMF unit string into a :class:`pint.Unit`.

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


def check(doc: Dict[str, Any], missing_units: str = "skip") -> UnitCheckResult:
    """Check dimensional consistency of a parsed YMF document.

    Non-intrusive per report §2.1: any unknown, coefficient, boundary
    condition, or initial condition with no ``units`` field is silently
    skipped rather than flagged as an error.

    ``missing_units`` says what a boundary condition that cannot be checked
    for want of declared units becomes: ``"skip"`` (the default) says
    nothing, ``"warn"`` and ``"error"`` add an issue of that severity naming
    what is undeclared.

    Checks performed (report §4.4):

    - Every declared unit string parses (bad strings -> ``error``).
    - Each boundary condition's ``units``, the units of its right-hand side
      g, agree with its type (see :func:`_check_boundary_conditions`):
      ``[g] = [u]`` for dirichlet, ``[g] = [q.n]`` for neumann, and
      ``[alpha][u] = [beta][q.n] = [g]`` for robin, where ``q = -a grad u``
      is the diffusive flux of the variable's equation.
    - Each ``initial_conditions[].units`` is dimensionally compatible with
      the ``unknowns[].units`` of the field it initializes.
    - A ``density`` or ``viscosity`` scale given as a number while a ``rho``
      (``ρ``) or ``mu`` (``μ``) coefficient is also given: ``error`` if they differ,
      ``warning`` if they agree. Derive the scale from the coefficient
      instead, so the value is stated once.

    Returns a :class:`UnitCheckResult`; ``bool(result)`` is ``True`` iff there
    are no ``error``-severity issues (warnings don't affect truthiness).
    """
    if missing_units not in MISSING_UNITS:
        raise ValueError("missing_units must be one of %s, not %r"
                         % (", ".join(MISSING_UNITS), missing_units))
    result = UnitCheckResult()
    problem = doc.get("model") or {}
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

    _check_boundary_conditions(strong_form, result, missing_units)
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

    # A scale stated as a number beside a coefficient for the same quantity
    # is one value written twice, and the copies drift apart when someone
    # changes one (a parameter variant that sets mu but not the viscosity
    # scale). Matched by the conventional names. Copies that disagree are an
    # error; copies that agree are a warning, since that agreement is what
    # has to be maintained by hand.
    coefficient_values = _coefficient_values(problem)
    coefficient_units = _coefficient_units(problem)
    scales = problem.get("characteristic_scales") or {}
    paired = set()

    for scale, names in _SCALE_COEFFICIENT_NAMES.items():
        spec = scales.get(scale)
        if not isinstance(spec, dict) or spec.get("value") is None:
            continue
        coefficient = next((n for n in names if n in coefficient_values), None)
        if coefficient is None:
            continue
        stated, actual = spec["value"], coefficient_values[coefficient]
        agree = math.isclose(stated, actual, rel_tol=1e-12)
        agreement = "agrees with" if agree else "differs from"
        paired.add((scale, coefficient))
        result.issues.append(UnitIssue(
            "warning" if agree else "error", f"characteristic_scales[{scale}]",
            f"value {stated:g} {agreement} coefficients[{coefficient}] = {actual:g}; "
            f"state it once with {{derived: true, formula: \"{coefficient}\"}}"))

    # The pairings above are the conventional names. A coefficient under any
    # other name -- nu, eta, a diffusivity -- was invisible to them, so the
    # same value written twice went unreported for want of being called the
    # expected thing.
    #
    # Matching on the value instead needs something to tell a duplicate from
    # a coincidence: in a non-dimensionalized problem half the quantities
    # are 1.0. Dimension is that something. Equal numbers of the same kind
    # are one value written twice; equal numbers of different kinds (a
    # density of 1.0 kg/m3 and a viscosity of 1.0 Pa*s) are a coincidence,
    # and saying so would be noise. Where either side states no units there
    # is no way to tell, so nothing is said.
    for scale, spec in scales.items():
        if not isinstance(spec, dict) or spec.get("value") is None:
            continue
        scale_units = spec.get("units")
        if not scale_units:
            continue
        for name, value in coefficient_values.items():
            if (scale, name) in paired:
                continue
            if not math.isclose(spec["value"], value, rel_tol=1e-12):
                continue
            name_units = coefficient_units.get(name)
            if not name_units or not units_compatible(scale_units, name_units):
                continue
            paired.add((scale, name))
            result.issues.append(UnitIssue(
                "warning", f"characteristic_scales[{scale}]",
                f"value {spec['value']:g} equals coefficients[{name}], and both "
                f"are {dimensionality(scale_units)}; state it once with "
                f"{{derived: true, formula: \"{name}\"}}"))

    return result


#: What :func:`check` does with a boundary condition it cannot check for
#: want of declared units: say nothing, warn, or call it an error.
MISSING_UNITS = ("skip", "warn", "error")


class _Undetermined(Exception):
    """A dimension that depends on something undeclared, or can't be read."""

    def __init__(self, reasons):
        super().__init__(", ".join(reasons))
        self.reasons = list(reasons)


class _Inconsistent(Exception):
    """An expression whose own terms differ in dimension."""


#: pint units, each of a base dimension of its own, defined once per process
#: and handed out to each check as stand-ins
_STANDIN_UNITS: List[str] = []


class _Dimensions:
    """The dimension of each quantity a spec names.

    Something with declared units has their dimension. Something without
    them gets a dimension of its own, a stand-in that no real dimension
    equals. A comparison then decides itself: if the two sides' ratio is
    dimensionless they agree, whatever the stand-ins were (alpha u against
    beta q.n needs no units for u, because u cancels); if what is left is
    made of real dimensions they disagree; and if a stand-in is left, the
    answer depends on what was not declared.

    Only dimensions are compared, never scale: km and m agree.
    """

    def __init__(self):
        self._standins: Dict[str, str] = {}   # pint dimension -> what it stands for

    def standin(self, what: str):
        for dim, known in self._standins.items():
            if known == what:
                return _UREG.Quantity(1, dim[1:-1])
        index = len(self._standins)
        if index == len(_STANDIN_UNITS):
            name = "ymf_undeclared_%d" % index
            _UREG.define("%s = [%s]" % (name, name))
            _STANDIN_UNITS.append(name)
        name = _STANDIN_UNITS[index]
        self._standins["[%s]" % name] = what
        return _UREG.Quantity(1, name)

    def of(self, unit_str: Optional[str], what: str):
        """``unit_str``'s dimension, or a stand-in for ``what``."""
        if not unit_str:
            return self.standin(what)
        # in base units: an offset unit (degC) is otherwise refused in a product
        return _UREG.Quantity(1, to_pint_unit(unit_str)).to_base_units()

    def undeclared(self, quantity) -> List[str]:
        return ["nothing declares " + self._standins[d]
                for d in quantity.dimensionality if d in self._standins]

    def compare(self, a, b):
        """None if a and b agree; else ("error", text) or ("missing", reasons)."""
        ratio = (a / b).to_reduced_units()
        if ratio.dimensionless:
            return None
        missing = self.undeclared(ratio)
        if missing:
            return ("missing", missing)
        return ("error", "%s against %s" % (_show(a), _show(b)))


def _show(quantity) -> str:
    """A derived quantity's units for a message, in SI base units ("kg/s³")."""
    if quantity.dimensionless:
        return "dimensionless"
    return format(quantity.to_base_units().units, "~P")


def _expression_dimension(expr, names: Dict[str, Any], dims: _Dimensions):
    """The dimension of a sympy expression, given one for each free symbol."""
    import sympy
    from sympy.core.function import AppliedUndef
    from ymf.symbolic.language import D

    if expr.is_Number:
        return _UREG.Quantity(1, "dimensionless")
    if isinstance(expr, sympy.Symbol):
        if expr.name not in names:
            raise _Undetermined(["%s is not a coordinate, unknown or "
                                 "coefficient" % expr.name])
        return names[expr.name]
    if isinstance(expr, D):
        return (_expression_dimension(expr.args[0], names, dims)
                / names[("x", "y", "z")[int(expr.args[1])]])
    if expr.is_Mul:
        out = _UREG.Quantity(1, "dimensionless")
        for arg in expr.args:
            out = out * _expression_dimension(arg, names, dims)
        return out
    if expr.is_Pow:
        base, exponent = expr.args
        if exponent.is_Number:
            return _expression_dimension(base, names, dims) ** float(exponent)
        # x**y: both dimensionless, or nothing is known
        for arg in expr.args:
            if not _expression_dimension(arg, names, dims).dimensionless:
                raise _Undetermined(["the dimension of the power %s is not known" % expr])
        return _UREG.Quantity(1, "dimensionless")
    if expr.is_Add or isinstance(expr, (sympy.Max, sympy.Min, sympy.Abs)):
        terms = [_expression_dimension(arg, names, dims) for arg in expr.args]
        for term in terms[1:]:
            verdict = dims.compare(term, terms[0])
            if verdict is None:
                continue
            if verdict[0] == "missing":
                raise _Undetermined(verdict[1])
            raise _Inconsistent("the terms of %s differ in dimension: %s" % (expr, verdict[1]))
        return terms[0]
    if isinstance(expr, AppliedUndef):
        # a coefficient, k(x, y): see _symbolic_names
        name = type(expr).__name__
        if name not in names:
            raise _Undetermined(["%s is not a coefficient" % name])
        return names[name]
    if isinstance(expr, sympy.Function):
        # sin, exp, log, ...: a dimensionless argument and result
        return _UREG.Quantity(1, "dimensionless")
    raise _Undetermined(["ymf.units cannot take the dimension of %s" % expr])


def _check_boundary_conditions(strong_form: Dict[str, Any], result: "UnitCheckResult",
                               missing_units: str) -> None:
    """The dimension of each boundary condition's g, against its type.

    ``units`` on a condition is the units of g, its right-hand side.

    - dirichlet, u = g: ``[g] = [u]``.
    - neumann, q.n = g, where q = -a grad u is the diffusive flux of the
      variable's own equation: ``[g] = [a][u]/[L]``.
    - robin, alpha u + beta q.n = g: ``[alpha][u] = [beta][q.n] = [g]``.
      The first equality holds or fails without units for u or g, since u
      appears on both sides; it needs only alpha, beta, a and the domain's
      length.

    Anything undeclared -- the unknown's units, g's, a coefficient's, the
    domain's -- leaves the comparison it enters undecided, and
    ``missing_units`` says what that becomes. Neumann and robin need
    ymf[symbolic] to read alpha, beta and the equations; without it they
    are undecided too.
    """
    bcs = [bc for bc in (strong_form.get("boundary_conditions") or [])
           if bc.get("type") in ("dirichlet", "neumann", "robin")]
    if not bcs:
        return
    dims = _Dimensions()
    unknowns = normalize_unknowns(strong_form.get("unknowns") or [])
    ranks = {u["name"]: int(u.get("rank", 0)) for u in unknowns}
    unknown_dims = {}
    for u in unknowns:
        try:
            unknown_dims[u["name"]] = dims.of(u.get("units"), "the units of %s" % u["name"])
        except UnitParseError:
            unknown_dims[u["name"]] = dims.standin("the units of %s" % u["name"])
    domain = strong_form.get("domain")
    domain_units = domain.get("units") if isinstance(domain, dict) else None

    def report(location, verdict):
        if verdict is None:
            return
        kind, detail = verdict
        if kind == "error":
            result.issues.append(UnitIssue("error", location, detail))
        elif missing_units != "skip":
            reasons = detail if isinstance(detail, list) else [detail]
            result.issues.append(UnitIssue(
                "warning" if missing_units == "warn" else "error", location,
                "not checked: %s" % "; ".join(dict.fromkeys(reasons))))

    flux_cache: Dict[str, Any] = {}
    for bc in bcs:
        var, kind = bc.get("variable"), bc["type"]
        location = f"boundary_conditions[region={bc.get('region')!r}, variable={var!r}]"
        if var not in unknown_dims:
            continue
        try:
            g = dims.of(bc.get("units"), "the units of g on %s" % bc.get("region"))
        except UnitParseError as exc:
            result.issues.append(UnitIssue("error", location, str(exc)))
            continue
        u = unknown_dims[var]
        if kind == "dirichlet":
            verdict = dims.compare(g, u)
            if verdict and verdict[0] == "error":
                declared = _unknown_units(unknowns, var)
                verdict = ("error",
                           f"units {bc['units']!r} ({dimensionality(bc['units'])}) incompatible "
                           f"with {var!r}'s declared units {declared!r} "
                           f"({dimensionality(declared)})")
            report(location, verdict)
            continue
        try:
            names, components, space = _symbolic_names(
                strong_form, dims, unknown_dims, ranks, domain_units, flux_cache)
            fluxes = [_flux_dimension(strong_form, c, names, dims, space, flux_cache)
                      for c in components[var]]
            q = fluxes[0]
            for other in fluxes[1:]:
                verdict = dims.compare(other, q)
                if verdict:
                    raise (_Undetermined(verdict[1]) if verdict[0] == "missing"
                           else _Inconsistent("the components of %s's flux differ: %s"
                                              % (var, verdict[1])))
            if kind == "neumann":
                verdict = dims.compare(g, q)
                if verdict and verdict[0] == "error":
                    verdict = ("error", "g is %s but the flux q.n of %s is %s"
                               % (bc["units"], var, _show(q)))
                report(location, verdict)
                continue
            alpha_u = _bc_expression(bc, "alpha", names, dims, space) * u
            beta_q = _bc_expression(bc, "beta", names, dims, space) * q
        except _Undetermined as exc:
            report(location, ("missing", exc.reasons))
            continue
        except _Inconsistent as exc:
            result.issues.append(UnitIssue("error", location, str(exc)))
            continue
        # Three quantities, three pairwise relations; two that hold make the
        # third hold, so it is undecided only while fewer than two are known.
        relations = [(alpha_u, "alpha u", beta_q, "beta q.n"),
                     (g, "g", alpha_u, "alpha u"),
                     (g, "g", beta_q, "beta q.n")]
        holding, missing, failed = 0, [], None
        for a, a_text, b, b_text in relations:
            verdict = dims.compare(a, b)
            if verdict is None:
                holding += 1
            elif verdict[0] == "missing":
                missing += verdict[1]
            elif failed is None:
                show = lambda q, text: bc["units"] if text == "g" else _show(q)
                failed = "%s is %s but %s is %s" % (a_text, show(a, a_text),
                                                    b_text, show(b, b_text))
        if failed:
            report(location, ("error", failed))
        elif holding < 2:
            report(location, ("missing", missing))


def _unknown_units(unknowns, name):
    return next(u.get("units") for u in unknowns if u["name"] == name)


def _symbolic_names(strong_form, dims, unknown_dims, ranks, domain_units, cache):
    """({symbol name: dimension}, {unknown: [component, ...]}) for the
    expressions of a spec, read with ymf.symbolic."""
    if "names" in cache:
        return cache["names"]
    try:
        import sympy  # noqa: F401
        from ymf.symbolic.language import Space
    except ImportError:
        raise _Undetermined(["reading the flux q needs ymf[symbolic]"])
    length = (dims.of(domain_units, "the domain's length units") if domain_units
              else dims.standin("the domain's length units"))
    names: Dict[str, Any] = {"x": length, "y": length, "z": length,
                             "t": dims.standin("the units of t")}
    coefficients = strong_form.get("coefficients") or {}
    for name, entry in coefficients.items():
        declared = entry.get("units") if isinstance(entry, dict) else None
        try:
            names[name] = dims.of(declared, "the units of %s" % name)
        except UnitParseError:
            names[name] = dims.standin("the units of %s" % name)
    dim = _space_dimension(strong_form.get("domain"))
    components = {}
    for name, rank in ranks.items():
        if rank == 0:
            components[name] = [name]
            names[name] = unknown_dims[name]
        else:
            components[name] = ["%s_%d" % (name, i) for i in range(dim)]
            for c in components[name]:
                names[c] = unknown_dims[name]
    # Coefficients stay named, so their units survive the parse: k becomes
    # k(x, y), whose free symbols are coordinates, which classify accepts
    # (a plain symbol k would be an undefined name to it).
    space = Space(dim, ranks)
    space.coefficients = {n: sympy.Function(n)(*space.x) for n in coefficients}
    cache["names"] = (names, components, space)
    return cache["names"]


def _space_dimension(domain) -> int:
    from ymf.domain import box_bounds
    box = box_bounds(domain)
    if box is not None:
        return len(box[0])
    if isinstance(domain, dict) and domain.get("dimension"):
        return int(domain["dimension"])
    if isinstance(domain, str):
        from ymf.symbolic.problem import _domain_dimension
        try:
            return _domain_dimension(domain)
        except Exception:
            pass
    raise _Undetermined(["the domain's dimension is not known"])


def _bc_expression(bc, key, names, dims, space):
    from ymf.symbolic.language import SymbolicError, parse_expression
    text = bc.get(key)
    if text is None:
        raise _Undetermined(["the condition gives no %s" % key])
    try:
        expr = parse_expression(str(text), space)
    except SymbolicError as exc:
        raise _Undetermined(["%s does not read: %s" % (key, exc)])
    return _expression_dimension(expr, names, dims)


def _flux_dimension(strong_form, component, names, dims, space, cache):
    """[q] for one component's equation: [a][phi]/[L] over its diffusion terms."""
    from ymf.symbolic.adr import classify
    from ymf.symbolic.language import SymbolicError, parse_equation
    if "equations" not in cache:
        residuals = []
        try:
            for text in strong_form.get("equations") or []:
                residuals.extend(parse_equation(text, space))
        except SymbolicError as exc:
            raise _Undetermined(["the equations do not read: %s" % exc])
        cache["equations"] = residuals
    residuals = cache["equations"]
    comps = space.components
    if len(residuals) != len(comps) or component not in comps:
        raise _Undetermined(["no flux q: %d equations for %d unknowns"
                             % (len(residuals), len(comps))])
    try:
        eq = classify(residuals[comps.index(component)], space)
    except SymbolicError as exc:
        raise _Undetermined(["no flux q: %s" % exc])
    terms = []
    for phi, tensor in eq.diffusion.items():
        for entry in tensor:
            if entry != 0:
                terms.append(_expression_dimension(entry, names, dims) * names[phi] / names["x"])
    if not terms:
        raise _Undetermined(["the equation for %s has no diffusive flux" % component])
    for term in terms[1:]:
        verdict = dims.compare(term, terms[0])
        if verdict:
            if verdict[0] == "missing":
                raise _Undetermined(verdict[1])
            raise _Inconsistent("the diffusion terms of the equation for %s differ "
                                "in dimension: %s" % (component, verdict[1]))
    return terms[0]


def _coefficient_units(problem: Dict[str, Any]) -> Dict[str, str]:
    """The units each coefficient states, for the ones that state any."""
    strong_form = problem.get("strong_form") or {}
    units: Dict[str, str] = {}
    for name, entry in (strong_form.get("coefficients") or {}).items():
        if isinstance(entry, dict) and entry.get("units"):
            units[name] = entry["units"]
    return units


#: Characteristic scales whose physical quantity is conventionally also a
#: coefficient of the equations, by that coefficient's usual name.
_SCALE_COEFFICIENT_NAMES = {"density": ("rho", "ρ"), "viscosity": ("mu", "μ")}


def _eval_formula(formula: str, scope: Dict[str, float]) -> Optional[float]:
    """Evaluate a characteristic-scale/dimensionless-number formula string.

    Formulas come from the YMF document author (trusted, same trust level
    as any other YAML content in the file) but evaluation is still
    restricted to arithmetic on the given scope: no builtins, no attribute
    or subscript access. Returns ``None`` (rather than raising) if the
    formula references a name not present in ``scope`` or otherwise fails to
    evaluate, since not all referenced scales may be available yet.
    """
    try:
        code = compile(formula, "<ymf-formula>", "eval")
    except SyntaxError:
        # Not Python: the spec notation ("1/(2π²κ)", "ρ U L / μ"), which
        # needs ymf[symbolic] to read.
        return _eval_notation(formula, scope)
    for name in code.co_names:
        if name not in scope:
            return None
    try:
        return eval(code, {"__builtins__": {}}, dict(scope))  # noqa: S307
    except Exception:
        return None


def _eval_notation(formula: str, scope: Dict[str, float]) -> Optional[float]:
    """Evaluate a formula in the spec notation, or None if it can't be."""
    try:
        import sympy
        from ymf.symbolic.language import Space
        from ymf.symbolic.notation import NotationError, Parser
    except ImportError:
        return None
    space = Space(3, {}, {k: sympy.Float(v) for k, v in scope.items()})
    try:
        value = Parser(space).expression(formula)
    except (NotationError, TypeError, ValueError):
        return None
    if getattr(value, "free_symbols", set()):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def characteristic_scale_values(problem: Dict[str, Any]) -> Dict[str, float]:
    """Resolve ``Problem.characteristic_scales`` to a flat ``{name: value}`` dict.

    Direct (``value:``) scales are taken as-is. ``derived: true`` scales
    (e.g. ``time: {derived: true, formula: "L_char / U_char"}``, report
    §5.2) are resolved by substituting the standard ``*_char`` aliases
    (``L_char``, ``U_char``, ``rho_char``, ``mu_char``, ``T_char``) for
    ``length``, ``velocity``, ``density``, ``viscosity``, ``temperature``
    respectively, and the numeric ``strong_form.coefficients`` by name.

    Deriving a scale from a coefficient (``viscosity: {derived: true,
    formula: "mu"}``) keeps a physical value in one place: a problem that
    changes ``mu`` changes the scale, and every dimensionless number built
    on it, with it. Unresolvable derived scales are omitted, not errored —
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

    coefficients = _coefficient_values(problem)

    # Fixed-point resolution: a couple of passes handle chained derived
    # scales (e.g. pressure derived from a derived velocity), though the
    # common case (report examples) resolves in one pass.
    for _ in range(len(derived) + 1):
        scope = dict(coefficients)
        scope.update({alias: values[key] for key, alias in alias_of.items() if key in values})
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


def _as_number(value: Any) -> Optional[float]:
    """Return ``value`` as a finite number, or ``None`` if it isn't one.

    ``strong_form.coefficients`` is untyped in the schema (``Any()``), so
    strictyaml hands back scalars as strings: ``kappa: 1.0e-3`` loads as
    ``'1.0e-3'``. Numeric strings are coerced here; formula strings such as
    ``"π²"`` (and ``"nan"``/``"inf"``) are not numbers and yield ``None``.
    """
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, str):
        try:
            number = float(value)
        except ValueError:
            return None
        return number if math.isfinite(number) else None
    return None


def _coefficient_values(problem: Dict[str, Any]) -> Dict[str, float]:
    """The numeric ``strong_form.coefficients``, as ``{name: value}``.

    A coefficient written as an expression is a number too, once evaluated:
    ``Q: "π²"`` is 9.8696, and ``λ: "ρ/(2μ) - sqrt(ρ²/(4μ²) + 4π²)"`` is a
    number as soon as ρ and μ are known. Those go through
    :func:`ymf.symbolic.problem.coefficient_numbers`, the same resolver
    ``adr_problem`` uses, so a dimensionless number built on such a
    coefficient gets the value the solver would get rather than quietly
    going missing.

    ``sympy`` comes with the ``frontend`` extra. Without it -- a bare
    ``ymf[units]`` install -- only the coefficients that are already numbers
    are returned, and expression-valued ones are left out as before.
    """
    strong_form = problem.get("strong_form") or {}
    spec = strong_form.get("coefficients") or {}
    try:
        from ymf.symbolic.problem import coefficient_numbers
    except ImportError:            # no sympy: numbers only, as before
        pass
    else:
        return coefficient_numbers(spec)
    values: Dict[str, float] = {}
    for name, entry in spec.items():
        value = _as_number(entry.get("value") if isinstance(entry, dict) else entry)
        if value is not None:
            values[name] = value
    return values


def compute_dimensionless_numbers(problem: Dict[str, Any]) -> Dict[str, float]:
    """Evaluate ``Problem.dimensionless_numbers[].formula`` (report §5.2).

    Formulas are evaluated against a scope built from
    :func:`characteristic_scale_values` (using bare names ``L``, ``U``,
    ``rho``, ``mu``, ``T`` as well as the ``*_char`` aliases, matching the
    variable names used in the report's example formulas, e.g.
    ``"rho * U * L / mu"``) plus a fixed gravitational constant ``g`` and any
    numeric ``strong_form.coefficients`` entries (a bare number or a
    ``{value: ..., units: ...}`` map; numeric strings, as strictyaml returns
    them, are coerced). Numbers whose formula can't be fully resolved
    (missing scale/coefficient) are omitted.
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

    for name, value in _coefficient_values(problem).items():
        scope.setdefault(name, value)

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
    problem = doc.get("model") or {}
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
