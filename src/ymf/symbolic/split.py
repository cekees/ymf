"""A split discretization -> one plain-data problem per sub-problem.

A discretization with ``subproblems`` solves the model as a sequence of
smaller problems in each step (see the guide, "Splitting a discretization
into sub-problems"). The whole model need not be one ADR system -- two-phase
level-set flow is not: ``v·∇φ`` is a Hamiltonian in another unknown -- but
each sub-problem is, once what it reads from the others is data. So each one
is classified on its own:

* its **unknowns** are the unknowns of the model equations it ``solves``
  (their ``for``), and the unknown of its own ``equation`` (its ``for``, by
  default the unknown it ``introduces``);
* every other field it reads must be declared in ``lags``; it is a **known
  field**, a function of the coordinates and time, which the plain data
  names (``v_0``) and whose gradient it names ``grad_<field>_<axis>``; a
  coefficient lagged ``via`` the sub-problem that carries it (``χ`` via the
  one that introduces ``χh``, approximating ``χ``) is that carrier;
* the symbols the discretization ``binds`` (``h``, ``ε = 1.5 h``) are
  substituted, and the sharp functions it ``regularize``\\ s become smooth
  ones: Proteus's smoothed Heaviside of width ε,
  ``H_ε(φ) = (1 + φ/ε + sin(πφ/ε)/π)/2`` on ``|φ| ≤ ε`` (``gf.H`` in
  ``equivalent_polynomials.h``), whose derivative is its ``gf.D``, and the
  smoothed sign ``2 H_ε − 1`` of its redistancing.

The result is plain data, like :func:`ymf.symbolic.adr_problem`'s: what a
dispatcher maps onto a solver's modules, and what a kernel generator reads.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import sympy

from ymf.symbolic.adr import adr_form, to_code
from ymf.symbolic.language import Space, SymbolicError, parse_equation, parse_expression

__all__ = ["split_problem", "smoothed_heaviside", "smoothed_sign", "REGULARIZATIONS"]


def smoothed_heaviside(phi, eps):
    """Proteus's ``gf.H(eps, phi)``: 0 below -ε, 1 above ε, smooth between."""
    return sympy.Piecewise(
        (sympy.Integer(1), phi > eps),
        (sympy.Integer(0), phi < -eps),
        ((1 + phi / eps + sympy.sin(sympy.pi * phi / eps) / sympy.pi) / 2, True))


def smoothed_sign(phi, eps):
    """``2 H_ε(φ) − 1``, the sign Proteus's redistancing uses."""
    return 2 * smoothed_heaviside(phi, eps) - 1


#: sharp function -> {regularization name: f(argument, width)}
REGULARIZATIONS = {
    "heaviside": {"smoothed_heaviside": smoothed_heaviside},
    "sign": {"smoothed_sign": smoothed_sign},
}


def split_problem(doc: Dict[str, Any], discretization: str) -> Dict[str, Any]:
    """The plain-data problem of each sub-problem of ``discretization``.

    ``doc`` is a validated, composed spec (its names were checked on load,
    ymf.references). Returns::

        {"discretization": name, "dim": d,
         "domain": {"lower": [...], "upper": [...]}, "time": [t0, t1] | None,
         "bound": {symbol: meaning or code},
         "initial": {component: code},          # the model unknowns'
         "subproblems": [{"name", "unknowns", "adr", "reads", "coefficients",
                          "regularized", "initial", "boundary_conditions",
                          ...}, ...]}

    A sub-problem that solves nothing it can classify (one that only spins
    the run up, or enforces a constraint it gives no equation for) has
    ``"classified": False`` and the reason instead of ``adr``.
    """
    from ymf.conditions import expand
    try:
        doc = expand(doc)          # aliases and condition sets: only their expansion is solved
    except ValueError as exc:
        raise SymbolicError(str(exc)) from exc
    from ymf.discretizations import subproblem_solver
    from ymf.normalize import normalize_equations
    from ymf.symbolic.problem import (_coefficients, _frame, _named_form,
                                      _state_coefficient_data)
    entries = (doc.get("solution_paths") or {}).get("discretizations") or []
    entry = next((d for d in entries if d.get("name") == discretization), None)
    if entry is None:
        raise SymbolicError("no discretization named %r (have: %s)"
                            % (discretization, ", ".join(d.get("name", "?") for d in entries)))
    subproblems = entry.get("subproblems") or []
    if not subproblems:
        raise SymbolicError("discretization %s is not split into sub-problems; "
                            "adr_problem reads it as one system" % (discretization,))
    strong, box, region, ranks, dim = _frame(doc)
    equations = {e["name"]: e for e in normalize_equations(strong.get("equations") or [])}
    spec = strong.get("coefficients") or {}
    introduced = {s["introduces"]["name"]: s for s in subproblems if s.get("introduces")}
    approximates = {s["introduces"]["approximates"]: s["introduces"]["name"]
                    for s in subproblems if (s.get("introduces") or {}).get("approximates")}
    entry_bound = _bind(entry.get("binds") or {}, {}, dim, "the discretization's binds")

    numbers = _coefficients(spec, dim, strict=False)
    domain = strong["domain"]
    time = None
    if box is not None:
        lower, upper = box
        if isinstance(domain, dict) and domain.get("time"):
            time = list(domain["time"])
    else:
        from ymf.symbolic.notation import Parser
        lower, upper, interval = Parser(Space(dim, {}, numbers)).domain(domain)
        time = list(interval) if interval else None
    out_subs = []
    for sub in subproblems:
        where = "sub-problem %s" % (sub["name"],)
        bound = _bind(sub.get("binds") or {}, entry_bound, dim, where + " binds")
        record: Dict[str, Any] = {"name": sub["name"]}
        for key in ("solves", "enforces", "writes_back", "pseudo_time", "spin_up",
                    "finite_element", "time_integration", "stabilization", "proteus"):
            if key in sub:
                record[key] = sub[key]
        record["solver"] = subproblem_solver(entry, sub)

        # its unknowns, and the text of every equation it classifies
        own: Dict[str, int] = {}
        texts: List[Tuple[str, str]] = []
        for name in sub.get("solves") or []:
            equation = equations[name]
            if equation["for"] is None:
                raise SymbolicError("%s solves %s, which does not say which unknown it is "
                                    "for; give the model's equations 'for'" % (where, name))
            own[equation["for"]] = ranks[equation["for"]]
            texts.append((name, equation["equation"]))
        if sub.get("equation"):
            target = sub.get("for") or (sub.get("introduces") or {}).get("name")
            if target is None:
                raise SymbolicError("%s gives an equation but not the unknown it is for: "
                                    "give for:" % (where,))
            own[target] = ranks.get(target, 0)
            texts.append(("equation", sub["equation"]))

        # what it reads: known fields, and coefficients read through a carrier
        fields: Dict[str, int] = {}
        preset: Dict[str, Any] = {}
        reads: Dict[str, Any] = {}
        for name, lag in (sub.get("lags") or {}).items():
            lag = lag if isinstance(lag, dict) else {"at": lag}
            target = lag.get("of", name)
            read = {"at": lag["at"]}
            if "via" in lag:
                read["via"] = lag["via"]
            if "of" in lag:
                read["of"] = lag["of"]
            if target in ranks or target in introduced:
                fields[name] = ranks.get(target, 0)
            elif target in spec and target in approximates:
                carrier = approximates[target]
                fields[carrier] = 0
                read["carried_by"] = carrier
            reads[name] = read
        known = Space(dim, own, {}, fields).known
        for name, read in reads.items():
            field = read.get("carried_by", name)
            if field in known:
                if "carried_by" in read:
                    preset[name] = known[field]
                read["components"] = ([field] if fields[field] == 0 else
                                      ["%s_%d" % (field, i) for i in range(dim)])
        if not own:
            record["classified"] = False
            record["why"] = ("it spins the run up" if sub.get("spin_up") else
                             "it enforces %r and gives no equation for it" % sub.get("enforces"))
            record["reads"] = reads
            out_subs.append(record)
            continue

        smooth = _regularizations(entry, sub, bound, dim, where)
        functions = {fn: (lambda a, *rest, s=smoother, w=width: s(a, w))
                     for fn, (_, width, smoother) in smooth.items()}
        supplied = [s for s in bound.values() if isinstance(s, sympy.Symbol)]
        coefficients = _coefficients(spec, dim, strict=False, unknowns=own,
                                     fields=fields, preset=preset, functions=functions,
                                     parameters=supplied)
        declared = (set(spec) | set(ranks) | set(introduced) | set(bound)
                    | set(sub.get("lags") or {}))
        used = _used([t for _, t in texts], spec, declared, carried=set(preset))
        values = dict(coefficients)
        values.update(bound)
        space = Space(dim, own, values, fields, supplied, functions)
        residuals: List[sympy.Expr] = []
        for label, text in texts:
            try:
                residuals.extend(parse_equation(text, space))
            except SymbolicError as exc:
                raise SymbolicError(
                    "%s, %s: %s (a field a sub-problem reads must be one of its unknowns "
                    "or declared in its lags)" % (where, label, exc)) from exc
        written = [t for _, t in texts] + [
            (spec[c].get("formula") if isinstance(spec[c], dict) else spec[c]) or ""
            for c in used if c not in preset] + [sub.get("initial_condition") or ""]
        sharp = _written(written, smooth)
        try:
            record["adr"] = adr_form(residuals, space)
        except SymbolicError as exc:
            raise SymbolicError("%s: %s" % (where, exc)) from exc
        record["unknowns"] = {name: ([name] if rank == 0 else
                                     ["%s_%d" % (name, i) for i in range(dim)])
                              for name, rank in own.items()}
        record["reads"] = reads
        named = _state_coefficient_data(
            {k: v for k, v in values.items() if k in spec and k in used}, space)
        if named:
            # by name too: each one inline over the state, and the
            # sub-problem's ADR form over those names (problem._named_form)
            state = [c for c in spec if c in named]
            form, inline = _named_form([t for _, t in texts], spec, values, state, dim,
                                       own, fields, supplied, functions, preset)
            for name, extra in inline.items():
                named[name].update(extra)
            record["coefficients"] = named
            record["named"] = form
        if sub.get("initial_condition"):
            record["initial"] = _initial(sub, own, ranks, values, functions, dim, where)
        if sharp:
            record["regularized"] = {fn: {"by": by, "width": to_code(width)}
                                     for fn, (by, width, _) in smooth.items() if fn in sharp}
        if sub.get("enforces"):
            try:
                parse_equation(sub["enforces"], space)
            except SymbolicError as exc:
                raise SymbolicError("%s, enforces: %s" % (where, exc)) from exc
        carried = {(s.get("introduces") or {}).get("approximates"): sub["introduces"]["name"]
                   for s in [sub] if (s.get("introduces") or {}).get("approximates")}
        record["boundary_conditions"] = _conditions(strong, own, carried, space, dim,
                                                    lower, upper, where)
        out_subs.append(record)

    initial = {}
    for ic in strong.get("initial_conditions") or []:
        name = ic["field"]
        if ic.get("formula") is not None:
            value = parse_expression(ic["formula"], Space(dim, {}, numbers))
        else:
            value = sympy.Float(ic.get("value"))
        if ranks.get(name, 0) == 1 and not isinstance(value, sympy.MatrixBase):
            value = sympy.Matrix([value] * dim)
        if isinstance(value, sympy.MatrixBase):
            initial.update({"%s_%d" % (name, i): to_code(v) for i, v in enumerate(value)})
        else:
            initial[name] = to_code(value)
    return {"discretization": discretization, "dim": dim,
            "domain": dict({"lower": lower, "upper": upper},
                           **({"region": region} if region else {})), "time": time,
            "bound": {k: (v if isinstance(v, str) else to_code(v))
                      for k, v in _describe_bound(entry.get("binds") or {}, entry_bound).items()},
            "initial": initial, "subproblems": out_subs}


def _conditions(strong, own, carried, space, dim, lower, upper, where) -> List[Dict[str, Any]]:
    """The boundary conditions one sub-problem imposes, as plain data.

    Those on its own unknowns, and those on a coefficient its introduced
    unknown carries (atmosphere's χ inflow is VOF's χh inflow). Conditions on
    what it only reads are another sub-problem's. Each entry has ``type``,
    ``variable``, ``region``, ``where`` (boolean code using ``tol``, None for
    the whole boundary), ``set`` if it came from one, the ``components`` it
    acts on, and its data: ``value`` (code; per component for a vector,
    except slip's normal component), ``series``, robin's ``alpha``/``beta``/
    ``g``, or periodic's ``axes``.
    """
    from ymf.normalize import boundary_conditions
    from ymf.symbolic.problem import (_components, _given, _number, _periodic_axes,
                                      _robin_part, _robin_rhs, _series)
    from ymf.symbolic.problem import _selector
    regions = {r["name"]: r for r in strong.get("boundary_regions") or []}
    out = []
    try:
        conditions = boundary_conditions(strong)
    except ValueError as exc:
        raise SymbolicError("%s: %s" % (where, exc)) from exc
    for bc in conditions:
        variable = bc["variable"]
        if variable in own:
            target, rank = variable, own[variable]
        elif variable in carried:
            target, rank = carried[variable], 0
        else:
            continue
        here = "%s, boundary condition %s on %s" % (where, bc["type"], bc["region"])
        entry: Dict[str, Any] = {"type": bc["type"], "variable": variable,
                                 "region": bc["region"],
                                 **_selector(regions[bc["region"]], strong["domain"], space)}
        if target != variable:
            entry["carried_by"] = target
        for key in ("component", "only", "imposed", "expanded_from"):
            if bc.get(key) is not None:
                entry[key] = bc[key]
        entry["components"] = ([target] if rank == 0 else
                               ["%s_%d" % (target, i) for i in range(dim)])
        kind = bc["type"]
        if kind == "periodic":
            if lower is None:
                raise SymbolicError("%s: periodic faces are a box's" % (here,))
            if regions[bc["region"]].get("where") is None:
                raise SymbolicError("%s: periodic faces are found by position, so the "
                                    "region needs where:, not tags:" % (here,))
            entry["axes"] = _periodic_axes(regions[bc["region"]]["where"], lower, upper,
                                           space, bc["region"])
            domain = strong["domain"]
            if isinstance(domain, dict) and domain.get("kind") == "piecewise_linear":
                from ymf.symbolic.problem import periodic_faces_match
                periodic_faces_match(domain, entry["axes"], lower, upper, bc["region"])
        elif kind == "robin":
            entry.update(alpha=_robin_part(bc, "alpha", space),
                         beta=_robin_part(bc, "beta", space), g=_robin_rhs(bc, space))
        elif bc.get("series") is not None:
            entry["series"] = _series(bc, target, 0 if bc.get("component") else rank)
        elif bc.get("formula") is not None or bc.get("value") is not None:
            value = (parse_expression(bc["formula"], space) if bc.get("formula") is not None
                     else _number(bc["value"]))
            if bc.get("component") or rank == 0:      # a projection is one number
                entry["value"] = to_code(_given(value, space, here))
            else:
                if not isinstance(value, sympy.MatrixBase):
                    value = sympy.Matrix([value] * dim)
                entry["value"] = _components(value, target, rank, space, here)
        out.append(entry)
    return out


def _bind(binds: Dict[str, Any], outer: Dict[str, Any], dim: int, where: str) -> Dict[str, Any]:
    """Resolve ``binds``: a symbol with only a ``meaning`` stays a symbol the
    solver supplies (h, the element size); one with a formula or value is
    that, in the others."""
    values = dict(outer)
    formulas = {}
    for name, b in binds.items():
        if b.get("value") is not None:
            values[name] = sympy.Float(b["value"])
        elif b.get("formula") is not None:
            formulas[name] = b["formula"]
        else:
            values[name] = sympy.Symbol(name, positive=True)
    for _ in range(len(formulas) + 1):
        for name, text in list(formulas.items()):
            try:
                values[name] = parse_expression(text, Space(dim, {}, values))
            except SymbolicError:
                continue
            del formulas[name]
    if formulas:
        raise SymbolicError("%s: %s could not be resolved" % (where, ", ".join(sorted(formulas))))
    return values


def _describe_bound(binds, values):
    return {name: (b.get("meaning") if b.get("formula") is None and b.get("value") is None
                   else values[name]) for name, b in binds.items()}


def _regularizations(entry, sub, bound, dim, where):
    """{sharp function: (by, width, smoother)} for this sub-problem."""
    chosen = dict(entry.get("regularize") or {})
    chosen.update(sub.get("regularize") or {})
    out = {}
    for function, rule in chosen.items():
        smoothers = REGULARIZATIONS.get(function, {})
        if rule["by"] not in smoothers:
            raise SymbolicError("%s: %s is regularized by %r, which ymf does not know (have: %s)"
                                % (where, function, rule["by"], ", ".join(smoothers) or "none"))
        if rule.get("width") is None:
            raise SymbolicError("%s: regularizing %s needs a width" % (where, function))
        width = parse_expression(rule["width"], Space(dim, {}, bound))
        out[function] = (rule["by"], width, smoothers[rule["by"]])
    return out


def _written(texts, smooth) -> set:
    """The regularized sharp functions these texts use."""
    from ymf.symbolic.notation import tokens
    names = set()
    for text in texts:
        names |= {v for k, v, _ in tokens(str(text)) if k == "NAME"}
    return {fn for fn in smooth if fn in names}


def _initial(sub, own, ranks, values, functions, dim, where):
    """An introduced unknown's initial condition: an expression in the model
    unknowns at the start, which the solver reads as known fields."""
    target = sub.get("for") or sub["introduces"]["name"]
    model_fields = {n: r for n, r in ranks.items() if n not in own}
    space = Space(dim, {}, {k: v for k, v in values.items() if not _names_unknown(v, own)},
                  model_fields, functions=functions)
    try:
        value = parse_expression(sub["initial_condition"], space)
    except SymbolicError as exc:
        raise SymbolicError("%s, initial_condition: %s" % (where, exc)) from exc
    return {target: to_code(value),
            "reads": sorted({f.func.__name__ for f in sympy.sympify(value).atoms(
                sympy.core.function.AppliedUndef)})}


def _used(texts: List[str], spec: Dict[str, Any], declared: set, carried=()) -> set:
    """The coefficients the texts name, directly or through other coefficients'
    formulas; not through a ``carried`` one, which is read rather than
    computed. A run of letters that is not a declared name counts for each
    letter, as the notation reads an undeclared ``ρU`` as ``ρ·U``."""
    from ymf.symbolic.notation import tokens

    def names(text):
        out = set()
        for kind, value, _ in tokens(str(text)):
            if kind == "NAME":
                out.update([value] if value in declared or len(value) == 1 else value)
        return out & set(spec)

    used, todo = set(), set().union(*[names(t) for t in texts]) if texts else set()
    while todo:
        name = todo.pop()
        used.add(name)
        if name in carried:
            continue
        entry = spec[name]
        formula = entry.get("formula") if isinstance(entry, dict) else entry
        if formula is not None:
            todo |= names(formula) - used
    return used


def _names_unknown(value, own) -> bool:
    entries = list(value) if isinstance(value, sympy.MatrixBase) else [value]
    names = set().union(*[sympy.sympify(e).free_symbols for e in entries])
    return any(s.name in own or s.name.rsplit("_", 1)[0] in own for s in names)
