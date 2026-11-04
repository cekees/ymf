"""A YMF problem specification -> an ADR problem a solver can run.

:func:`adr_problem` reads a spec -- ``strong_form.domain`` and
``equations``, the coefficients, the boundary regions' ``where``,
boundary and initial conditions, and the analytical solution's
``formula``, each written once in the notation of
:mod:`ymf.symbolic.notation` -- and returns a plain-data description:

    {"adr": {...},            # ymf.symbolic.adr.adr_form
     "unknowns": {name: [component, ...]},
     "domain": {"lower": [...], "upper": [...]},
     "dirichlet": [{"component", "region", "where", "value"}, ...],
     "periodic": [{"component", "region", "axes"}, ...],
     "initial": {component: code},
     "exact": {component: code},
     "time": [t0, t1] or None}

Every expression in it is a numpy code string in x, y, z, t (and, in the
ADR coefficients, the unknowns), so the solver side needs numpy and nothing
else. ``where`` is a boolean code string that also uses ``tol``, which the
consumer supplies (a length tolerance for testing points on a boundary).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import sympy

from ymf.domain import box_bounds
from ymf.symbolic.adr import adr_form, to_code
from ymf.symbolic.language import (D, Dt, Space, SymbolicError, parse_equation,
                                   parse_expression)

__all__ = ["adr_problem", "strong_form_equations", "state_coefficients"]


def _number(value):
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return sympy.Float(value)
    if isinstance(value, str):
        try:
            return sympy.Float(float(value))
        except ValueError:
            return None
    return None


def _coefficients(spec: Dict[str, Any], dim: int, strict: bool = True,
                  unknowns: Optional[Dict[str, int]] = None,
                  fields: Optional[Dict[str, int]] = None,
                  preset: Optional[Dict[str, Any]] = None,
                  functions: Optional[Dict[str, Any]] = None,
                  parameters: Sequence[sympy.Symbol] = ()) -> Dict[str, sympy.Expr]:
    """Resolve coefficient values and formulas; formulas may use each other.

    ``strict`` decides what an unresolvable formula means. Building a problem
    to solve, it is an error: a coefficient the equations use but nothing
    defines cannot be run. Reading the same spec to report on it -- which is
    what :func:`coefficient_numbers` does for ``ymf.units`` -- it is not, and
    the entry is simply left out.

    ``unknowns`` (name -> rank) lets a formula name the unknowns too: a
    density ``ρw (1 - χ) + ρa χ`` with ``χ`` an unknown, a conductivity
    ``1 + u²``. Such a coefficient is a function of the state, substituted
    into the equations like any other, so its dependence on the unknowns
    shows up in the ADR form's ``depends_on`` and derivatives. Without
    ``unknowns`` it is left unresolved, which is what a reader that wants
    numbers needs. A formula may name the unknowns, not their derivatives.

    ``fields`` (name -> rank) are known fields a formula may name too (a
    sub-problem's lagged reads, see ymf.symbolic.split), and ``preset`` gives
    coefficients their value outright, in place of their formula (a
    coefficient read through the sub-problem that carries it: χ -> χh).
    """
    unknowns = dict(unknowns or {})
    clash = sorted(set(unknowns) & set(spec))
    if clash:
        raise SymbolicError("coefficient %s has the name of an unknown; a name is one "
                            "or the other" % ", ".join(clash))
    state = Space(dim, unknowns)
    allowed = ({"x", "y", "z", "t"} | set(state.components)
               | {p.name for p in parameters})
    values: Dict[str, sympy.Expr] = dict(preset or {})
    formulas: Dict[str, str] = {}
    for name, entry in spec.items():
        if name in values:
            continue                  # preset: read from elsewhere, not computed
        if isinstance(entry, dict):
            if entry.get("value") is not None:
                v = _number(entry["value"])
                if v is None:
                    raise SymbolicError("coefficient %s: value %r is not a number"
                                        % (name, entry["value"]))
                values[name] = v
            elif entry.get("formula") is not None:
                formulas[name] = entry["formula"]
            # units only: declared, no value -- an error only if used
        else:
            v = _number(entry)
            if v is not None:
                values[name] = v
            else:
                formulas[name] = str(entry)
    for _ in range(len(formulas) + 1):
        progressed = False
        for name, text in list(formulas.items()):
            try:
                value = parse_expression(text, Space(dim, unknowns, values, fields,
                                                     functions=functions))
            except SymbolicError:
                continue
            entries = list(value) if isinstance(value, sympy.MatrixBase) else [value]
            stray = set().union(*[sympy.sympify(e).free_symbols for e in entries])
            if {s.name for s in stray} - allowed:
                continue          # depends on a coefficient not resolved yet
            if any(sympy.sympify(e).has(D, Dt) for e in entries):
                raise SymbolicError(
                    "coefficient %s: %r takes a derivative of an unknown; a coefficient "
                    "may depend on the unknowns, not (yet) on their derivatives"
                    % (name, text))
            values[name] = value
            del formulas[name]
            progressed = True
        if not progressed:
            break
    if formulas and strict:
        raise SymbolicError("coefficients %s could not be resolved (unknown names, "
                            "or a cycle)" % ", ".join(sorted(formulas)))
    return values


def coefficient_numbers(spec: Dict[str, Any], dim: int = 3) -> Dict[str, float]:
    """The coefficients that are numbers, with expressions evaluated.

    ``ymf.units`` needs the same answer about a coefficient that
    :func:`adr_problem` gets, or the same spec means two things: a
    dimensionless number built on a coefficient written as an expression
    would quietly go missing while the solver resolved it without
    difficulty. Both now come through this one resolver.

    Only coefficients that reduce to a real, finite scalar are returned. A
    vector (``f: {formula: "(G, 0)"}``), one that still depends on x or t,
    and one that nothing defines are all left out -- they are coefficients,
    but they are not numbers.
    """
    values = _coefficients(spec, dim, strict=False)
    numbers: Dict[str, float] = {}
    for name, value in values.items():
        if isinstance(value, sympy.MatrixBase) or value.free_symbols:
            continue
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue            # complex, or not evaluable to a float
        if math.isfinite(number):
            numbers[name] = number
    return numbers


def _given(value, space: Space, where: str):
    """``value``, refused if it depends on the unknowns.

    Boundary data, initial conditions and exact solutions are evaluated by
    the solver at points, in x, y, z and t alone. Since a coefficient may now
    be a function of the state, naming one there could slip an unknown in.
    """
    entries = list(value) if isinstance(value, sympy.MatrixBase) else [value]
    named = set().union(*[sympy.sympify(e).free_symbols for e in entries])
    state = sorted(s.name for s in named & set(space.component_symbols))
    if state:
        raise SymbolicError("%s depends on the unknowns (%s), directly or through a "
                            "coefficient; give it in x, y, z and t"
                            % (where, ", ".join(state)))
    return value


def state_coefficients(strong_form: Dict[str, Any], dim: int = 3) -> Dict[str, List[str]]:
    """The coefficients that are functions of the state, and of what.

    ``{coefficient: [component, ...]}`` for every coefficient whose formula
    names an unknown, directly or through another coefficient. Such a
    coefficient is not a number, so a scale or a dimensionless number built
    on it has no value; ``ymf.units`` says so rather than leave it out
    silently. Coefficients that do not resolve are left out, as in
    :func:`coefficient_numbers`.
    """
    ranks = {}
    for u in strong_form.get("unknowns") or []:
        if isinstance(u, str):
            ranks[u] = 0
        else:
            ranks[u["name"]] = int(u.get("rank", 0))
    values = _coefficients(strong_form.get("coefficients") or {}, dim, strict=False,
                           unknowns=ranks)
    components = {s.name: s for s in Space(dim, ranks).component_symbols}
    out = {}
    for name, value in values.items():
        entries = list(value) if isinstance(value, sympy.MatrixBase) else [value]
        named = set().union(*[sympy.sympify(e).free_symbols for e in entries])
        state = sorted(s.name for s in named if s.name in components)
        if state:
            out[name] = state
    return out


def _components(value, name: str, rank: int, space: Space, where: str) -> Dict[str, str]:
    """Per-component code strings for a scalar or vector value."""
    _given(value, space, where)
    if rank == 0:
        if isinstance(value, sympy.MatrixBase):
            raise SymbolicError("%s: %s is a scalar but got a vector" % (where, name))
        return {name: to_code(value)}
    if not isinstance(value, sympy.MatrixBase) or value.shape != (space.dim, 1):
        raise SymbolicError("%s: %s is a %d-vector; give it as (a, b%s)"
                            % (where, name, space.dim, ", c" if space.dim == 3 else ""))
    return {"%s_%d" % (name, i): to_code(value[i]) for i in range(space.dim)}


def _robin_part(bc, key, space):
    """alpha or beta of alpha u + beta (q.n) = g, as code."""
    text = bc.get(key)
    if text is None:
        raise SymbolicError(
            "robin condition on %s: give alpha and beta, the coefficients of "
            "u and of the outward flux in alpha u + beta (q.n) = g"
            % (bc["region"],))
    return to_code(_given(parse_expression(str(text), space), space,
                          "robin condition on %s: %s" % (bc["region"], key)))


def _robin_rhs(bc, space):
    """g, from formula or value like any other right-hand side."""
    if bc.get("formula") is not None:
        return to_code(_given(parse_expression(bc["formula"], space), space,
                              "robin condition on %s: g" % (bc["region"],)))
    if bc.get("value") is not None:
        return to_code(_number(bc["value"]))
    raise SymbolicError(
        "robin condition on %s: give g as value or formula" % (bc["region"],))


def _series(bc, name, rank):
    """``[[t, value], ...]``, checked.

    A value supplied as data rather than written as a formula: a measured or
    design record is the one thing an expression cannot be. It is a value
    form, not a kind of condition -- under ``neumann`` it is a discharge
    hydrograph, under ``dirichlet`` a stage hydrograph or a measured inlet
    value -- so the entry carries ``series`` where another carries
    ``value``, and everything downstream treats the two the same way.
    """
    where = "%s condition on %s" % (bc["type"], bc["region"])
    if rank != 0:
        raise SymbolicError(
            "%s: %s is a vector, and a series is one number per time. Give "
            "one condition per component." % (where, name))
    series = bc.get("series")
    if any(len(point) != 2 for point in series):
        raise SymbolicError("%s: each sample is [t, value]" % (where,))
    times = [float(t) for t, _ in series]
    if len(times) < 2:
        raise SymbolicError(
            "%s: one sample is a constant; give value instead" % (where,))
    if any(b <= a for a, b in zip(times, times[1:])):
        raise SymbolicError(
            "%s: times must increase, got %s"
            % (where, ", ".join(repr(t) for t in times)))
    return [[float(t), float(v)] for t, v in series]


def region_tags(region: Dict[str, Any], domain: Any) -> List[int]:
    """A tag region's boundary flags: its ints, and its names resolved
    through the domain's ``boundary_tags``."""
    names = (domain.get("boundary_tags") if isinstance(domain, dict) else None) or {}
    flags = []
    for tag in region["tags"]:
        if isinstance(tag, int) or (isinstance(tag, str) and tag.lstrip("-").isdigit()):
            flags.append(int(tag))
        elif tag in names:
            flags.append(int(names[tag]))
        else:
            raise SymbolicError("boundary region %s: tag %r is not one of the domain's "
                                "boundary_tags (have: %s)"
                                % (region["name"], tag, ", ".join(sorted(names)) or "none"))
    return flags


def _selector(region: Dict[str, Any], domain: Any, space: Space) -> Dict[str, Any]:
    """How a solver finds a region's boundary points.

    ``{"where": code}`` for a predicate (None: the whole boundary), or
    ``{"flags": [ints]}`` for a tag region, which the solver tests against the
    boundary flag of the point. A tag region has no ``where`` key at all, so a
    solver that does not know ``flags`` fails on it rather than reading a
    missing predicate as "the whole boundary".
    """
    if region.get("tags") is not None:
        return {"flags": region_tags(region, domain)}
    return {"where": _where(region["where"], space)}


def _where(text: str, space: Space) -> Optional[str]:
    """A region's ``where`` predicate -> boolean code using ``tol`` (None: the whole boundary)."""
    from ymf.symbolic.notation import NotationError, Parser
    try:
        return Parser(space).predicate(text)
    except NotationError as exc:
        raise SymbolicError("boundary region %r: %s" % (text, exc)) from exc


def _domain_dimension(text: str) -> int:
    """Count the space intervals of a domain before reading their bounds.

    Coefficient formulas need the dimension (they may use x, y, z) and a
    domain's bounds may name coefficients, so the dimension comes first,
    from the number of '×' factors before any ', t ∈ ...'.
    """
    from ymf.symbolic.notation import tokens
    count, depth = 1, 0
    for kind, value, _ in tokens(text):
        if value in "([":
            depth += 1
        elif value in ")]":
            depth -= 1
        elif value == "×" and depth == 0:
            count += 1
        elif value == "," and depth == 0:
            break
    return count


def _periodic_axes(where, lower, upper, space: Space, region: str) -> List[int]:
    """The axes whose two opposite faces both lie in a periodic region."""
    code = _where(where, space)
    if code is None:
        raise SymbolicError("periodic region %s must name its two faces, e.g. "
                            "'x = 0 or x = 4', not the whole boundary" % (region,))
    import math
    test = compile(code, "<where>", "eval")
    names = "xyz"[:space.dim]
    tol = 1e-8 * max(u - l for u, l in zip(upper, lower))
    axes = []
    for a in range(space.dim):
        hits = []
        for face in (lower[a], upper[a]):
            point = [(l + u) / 2 for l, u in zip(lower, upper)]
            point[a] = face
            ns = dict(zip(names, point), tol=tol)
            hits.append(bool(eval(test, {"__builtins__": {"abs": abs}, "numpy": math}, ns)))
        if all(hits):
            axes.append(a)
        elif any(hits):
            raise SymbolicError("periodic region %s contains only one face of axis %s"
                                % (region, names[a]))
    if not axes:
        raise SymbolicError("periodic region %s contains no pair of opposite faces" % (region,))
    return axes


def _frame(doc: Dict[str, Any]):
    """The strong form's setting, before any equation is read.

    Returns ``(strong, box, region, ranks, dim)``: the strong form, the
    domain's box (the bounding box of a region), the region itself (None
    for a box), each unknown's rank, and the space dimension. Shared by
    :func:`adr_problem` and the split reader, ymf.symbolic.split.
    """
    problem = doc["model"]
    strong = problem["strong_form"]
    domain = strong["domain"]
    box = box_bounds(domain)
    region = None
    if box is None and not isinstance(domain, str):
        if domain.get("kind") != "piecewise_linear":
            raise SymbolicError(
                "strong_form.domain: adr_problem can run a box or a "
                "piecewise_linear region; this one is kind %r. Other kinds "
                "validate and are carried, but nothing solves them yet."
                % (domain.get("kind"),))
        # A boundary representation: the solver meshes it (its
        # discretization gives mesh.h). The bounding box stands in for the
        # box where one is needed: the boundary tolerance's scale.
        from ymf.closure import plain
        region = plain(domain)
        vertices = [[float(c) for c in v] for v in region["vertices"]]
        box = [[min(c) for c in zip(*vertices)], [max(c) for c in zip(*vertices)]]
    dim = len(box[0]) if box else _domain_dimension(domain)
    if not strong.get("equations"):
        raise SymbolicError("strong_form.equations is needed to run a spec")
    ranks = {}
    for u in strong["unknowns"]:
        if isinstance(u, str):
            ranks[u] = 0
        else:
            ranks[u["name"]] = int(u.get("rank", 0))
    return strong, box, region, ranks, dim


def _posed(doc: Dict[str, Any]):
    """What adr_problem and strong_form_equations share.

    Returns ``(strong, box, region, ranks, coefficients, space, residuals)``.
    """
    strong, box, region, ranks, dim = _frame(doc)
    coefficients = _coefficients(strong.get("coefficients") or {}, dim, unknowns=ranks)
    space = Space(dim, ranks, coefficients)
    residuals: List[sympy.Expr] = []
    from ymf.normalize import equations_in_order
    try:
        texts = equations_in_order(strong)
    except ValueError as exc:
        raise SymbolicError(str(exc)) from exc
    for text in texts:
        residuals.extend(parse_equation(text, space))
    return strong, box, region, ranks, coefficients, space, residuals


def strong_form_equations(doc: Dict[str, Any]) -> Tuple[Space, List[sympy.Expr]]:
    """The Space a spec's strong form is written in, and its scalar residuals.

    Still in sympy, one residual per scalar component, as adr_problem
    classifies them (ymf.symbolic.adr.classify). For a tool that renders the
    problem rather than solving it -- writing its weak form in the notation
    (ymf.symbolic.to_notation), say.
    """
    _, _, _, _, _, space, residuals = _posed(doc)
    return space, residuals


def adr_problem(doc: Dict[str, Any], hamiltonian_gradients=()) -> Dict[str, Any]:
    """Build the plain-data ADR problem for a validated (composed) spec.

    ``hamiltonian_gradients``: see :func:`ymf.symbolic.adr.classify`.
    """
    from ymf.symbolic.notation import NotationError, Parser

    strong, box, region, ranks, coefficients, space, residuals = _posed(doc)
    domain = strong["domain"]
    dim = space.dim
    if box is not None:
        lower, upper = box
        interval = domain.get("time")
        time = tuple(interval) if interval else None
    else:
        try:
            lower, upper, time = Parser(Space(dim, {}, coefficients)).domain(domain)
        except NotationError as exc:
            raise SymbolicError("strong_form.domain: %s" % exc) from exc

    adr = adr_form(residuals, space, hamiltonian_gradients=tuple(hamiltonian_gradients))

    regions = {r["name"]: r for r in strong.get("boundary_regions", [])}
    dirichlet = []
    neumann = []
    robin = []
    periodic = []
    from ymf.normalize import boundary_conditions
    try:
        conditions = boundary_conditions(strong)
    except ValueError as exc:
        raise SymbolicError(str(exc)) from exc
    for bc in conditions:
        if bc["type"] not in ("dirichlet", "neumann", "periodic", "robin"):
            raise SymbolicError("boundary condition on %s: %s is not a kind of "
                                "boundary condition that runs in one ADR system (%s "
                                "conditions are for a split discretization's solver)"
                                % (bc["region"], bc["type"], bc["type"]))
        name = bc["variable"]
        if name not in ranks:
            raise SymbolicError("boundary condition on %s: %s is not an unknown%s"
                                % (bc["region"], name,
                                   " (a condition on a coefficient is for a split "
                                   "discretization's solver)"
                                   if name in (strong.get("coefficients") or {}) else ""))
        if bc["region"] not in regions:
            raise SymbolicError("boundary condition: no boundary region named %r" % bc["region"])
        selector = _selector(regions[bc["region"]], domain, space)
        where = regions[bc["region"]].get("where")
        if bc["type"] == "periodic":
            if where is None:
                raise SymbolicError(
                    "periodic condition on %s: periodic faces are found by position, "
                    "so the region needs where:, not tags:" % (bc["region"],))
            if region is not None:
                raise SymbolicError(
                    "periodic condition on %s: periodic faces are a box's; "
                    "this domain is piecewise_linear" % (bc["region"],))
            axes = _periodic_axes(where, lower, upper, space, bc["region"])
            comps = [name] if ranks[name] == 0 else ["%s_%d" % (name, i) for i in range(dim)]
            periodic.extend({"component": c, "region": bc["region"], "axes": axes} for c in comps)
            continue
        if bc["type"] == "robin":
            if ranks[name] != 0:
                raise SymbolicError(
                    "robin condition on %s: %s is a vector. alpha u + beta q.n = g "
                    "is one relation, so give one per component."
                    % (bc["region"], name))
            robin.append({"component": name, "region": bc["region"], **selector,
                          "alpha": _robin_part(bc, "alpha", space),
                          "beta": _robin_part(bc, "beta", space),
                          "g": _robin_rhs(bc, space)})
            continue
        entries = dirichlet if bc["type"] == "dirichlet" else neumann
        if bc.get("series") is not None:
            # Tabulated rather than written: a discharge hydrograph under
            # neumann, a stage hydrograph or a measured inlet value under
            # dirichlet. One number per time, so a vector wants one
            # condition per component.
            entries.append({"component": name, "region": bc["region"], **selector,
                            "series": _series(bc, name, ranks[name])})
            continue
        if bc.get("formula") is not None:
            value = parse_expression(bc["formula"], space)
        elif bc.get("value") is not None:
            value = _number(bc["value"])
            if ranks[name] == 1:
                value = sympy.Matrix([value] * dim)
        else:
            raise SymbolicError("boundary condition on %s: give value or formula" % bc["region"])
        for component, code in _components(value, name, ranks[name], space,
                                           "boundary condition on %s" % bc["region"]).items():
            entries.append({"component": component, "region": bc["region"], **selector,
                            "value": code})

    initial = {}
    for ic in strong.get("initial_conditions") or []:
        name = ic["field"]
        if ic.get("formula") is not None:
            value = parse_expression(ic["formula"], space)
        else:
            value = _number(ic.get("value"))
        initial.update(_components(value, name, ranks.get(name, 0), space,
                                   "initial condition for %s" % name))

    exact = {}
    for entry in (doc.get("solution_paths") or {}).get("analytical") or []:
        try:
            values = Parser(space).assignments(entry["solution"]["formula"], ranks)
        except NotationError as exc:
            raise SymbolicError("analytical solution %s: %s" % (entry["name"], exc)) from exc
        for name, value in values.items():
            exact.update(_components(value, name, ranks[name], space,
                                     "analytical solution %s" % entry["name"]))
        break

    has_mass = any("mass" in e for e in adr["equations"])
    up_to_constant = _up_to_a_constant(residuals, space, ranks, dirichlet)
    unknowns = {name: ([name] if rank == 0 else ["%s_%d" % (name, i) for i in range(dim)])
                for name, rank in ranks.items()}
    out = {"adr": adr, "unknowns": unknowns,
           "domain": dict({"lower": lower, "upper": upper},
                          **({"region": region} if region else {})),
           # A neumann entry's value is the OUTWARD diffusive flux q.n with
           # q = -a grad u, which is the sign Proteus's
           # diffusiveFluxBoundaryConditions wants (see
           # test/poisson_2d/poisson_het_2d_p.py, which checks it against an
           # exact solution).
           "dirichlet": dirichlet, "neumann": neumann, "robin": robin,
           "periodic": periodic, "initial": initial, "exact": exact,
           # The time interval the problem is posed on; None for a steady
           # problem. A time derivative in a model used steadily (Couette
           # from the Navier-Stokes model) is dropped by a steady solve.
           "time": list(time) if time else None, "has_mass": has_mass}
    if up_to_constant:
        # Scalar unknowns the problem determines only up to a constant:
        # they appear in the equations only through their gradients, and
        # no Dirichlet condition fixes them anywhere (the pressure of an
        # incompressible flow with the velocity given on the whole
        # boundary). A solver puts the constant in the null space of its
        # operator; errors in them are measured modulo constants. Absent
        # when there are none.
        out["up_to_constant"] = up_to_constant
    named = _state_coefficient_data(coefficients, space)
    if named:
        # The coefficients that follow the state, by name. The ADR form above
        # has them substituted, which is all a solver of one system needs; a
        # code generator or a dispatcher wants the name too, to emit rho as
        # one inline function of the state rather than its expansion in every
        # slot. Absent when there are none, so a problem without them reads
        # exactly as before.
        out["coefficients"] = named
    return out


def _state_coefficient_data(coefficients, space: Space) -> Dict[str, Any]:
    """``{name: {code, depends_on, d}}`` for each coefficient of the state.

    ``code`` is the numpy code of the coefficient (a list for a vector),
    ``depends_on`` maps each component it depends on to ``linear`` or
    ``nonlinear``, and ``d`` holds its derivative in each. A sharp function
    of an unknown (``heaviside(φ)``) has no derivative with a value at a
    point; it is marked ``sharp: true`` and its ``d`` is left out, for the
    discretization to regularize.
    """
    from sympy.core.function import AppliedUndef
    from ymf.symbolic.adr import _flag
    unknowns = set(space.component_symbols)
    out: Dict[str, Any] = {}
    for name, value in coefficients.items():
        entries = list(value) if isinstance(value, sympy.MatrixBase) else [value]
        entries = [sympy.sympify(e) for e in entries]
        fields = sorted({f.func.__name__ for e in entries for f in e.atoms(AppliedUndef)})
        if not (set().union(*[e.free_symbols for e in entries]) & unknowns or fields):
            continue
        flags: Dict[str, str] = {}
        for e in entries:
            for component, flag in _flag(e, space).items():
                if flags.get(component) != "nonlinear":
                    flags[component] = flag
        vector = isinstance(value, sympy.MatrixBase)
        record: Dict[str, Any] = {
            "code": [to_code(e) for e in entries] if vector else to_code(entries[0]),
            "depends_on": flags}
        if fields:
            # known fields it reads (a sub-problem's lagged data, ymf.symbolic.split)
            record["reads"] = fields
        try:
            record["d"] = {
                c: ([to_code(sympy.diff(e, s)) for e in entries] if vector
                    else to_code(sympy.diff(entries[0], s)))
                for c, s in zip(space.components, space.component_symbols) if c in flags}
        except SymbolicError:
            record["sharp"] = True
        out[name] = record
    return out


def _up_to_a_constant(residuals, space, ranks, dirichlet) -> List[str]:
    from ymf.symbolic.language import D
    fixed = {e["component"] for e in dirichlet}
    out = []
    for name, rank in ranks.items():
        if rank != 0 or name in fixed:
            continue
        symbol = sympy.Symbol(name, real=True)
        stand_in = sympy.Dummy()
        appears = False
        for r in residuals:
            gradients = {d: stand_in for d in sympy.sympify(r).atoms(D) if d.args[0] == symbol}
            rest = sympy.sympify(r).xreplace(gradients)
            if rest.has(symbol):
                break                  # appears itself, not only through gradients
            appears = appears or bool(gradients)
        else:
            if appears:
                out.append(name)
    return out
