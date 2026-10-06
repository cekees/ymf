"""A YMF problem specification -> an ADR problem a solver can run.

:func:`adr_problem` reads the machine-readable parts of a spec --
``geometry``, ``strong_form.equations``, the coefficients, boundary
regions with ``where``, boundary and initial conditions, and analytical
``expressions`` -- and returns a plain-data description:

    {"adr": {...},            # ymf.symbolic.adr.adr_form
     "unknowns": {name: [component, ...]},
     "geometry": {"lower": [...], "upper": [...]},
     "dirichlet": [{"component", "region", "where", "value"}, ...],
     "periodic": [{"component", "region", "axes"}, ...],
     "initial": {component: code},
     "exact": {component: code},
     "transient": bool}

Every expression in it is a numpy code string in x, y, z, t (and, in the
ADR coefficients, the unknowns), so the solver side needs numpy and nothing
else. ``where`` is a boolean code string that also uses ``tol``, which the
consumer supplies (a length tolerance for testing points on a boundary).
"""

from __future__ import annotations

from typing import Any, Dict, List

import sympy

from ymf.symbolic.adr import adr_form, to_code
from ymf.symbolic.language import Space, SymbolicError, parse_equation, parse_expression

__all__ = ["adr_problem"]


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


def _coefficients(spec: Dict[str, Any], dim: int) -> Dict[str, sympy.Expr]:
    """Resolve coefficient values and formulas; formulas may use each other."""
    values: Dict[str, sympy.Expr] = {}
    formulas: Dict[str, str] = {}
    for name, entry in spec.items():
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
                value = parse_expression(text, Space(dim, {}, values))
            except SymbolicError:
                continue
            stray = (value.free_symbols if not isinstance(value, sympy.MatrixBase)
                     else set().union(*[e.free_symbols for e in value]))
            if {s.name for s in stray} - {"x", "y", "z", "t"}:
                continue          # depends on a coefficient not resolved yet
            values[name] = value
            del formulas[name]
            progressed = True
        if not progressed:
            break
    if formulas:
        raise SymbolicError("coefficients %s could not be resolved (unknown names, "
                            "or a cycle)" % ", ".join(sorted(formulas)))
    return values


def _components(value, name: str, rank: int, space: Space, where: str) -> Dict[str, str]:
    """Per-component code strings for a scalar or vector value."""
    if rank == 0:
        if isinstance(value, sympy.MatrixBase):
            raise SymbolicError("%s: %s is a scalar but got a vector" % (where, name))
        return {name: to_code(value)}
    if not isinstance(value, sympy.MatrixBase) or value.shape != (space.dim, 1):
        raise SymbolicError("%s: %s is a %d-vector; give it as (a, b%s)"
                            % (where, name, space.dim, ", c" if space.dim == 3 else ""))
    return {"%s_%d" % (name, i): to_code(value[i]) for i in range(space.dim)}


def _where(text: str, space: Space) -> str:
    """'x = 0 and y = 1 or x = 4' -> numpy-free boolean code using tol."""
    alternatives = []
    for alternative in text.split(" or "):
        tests = []
        for atom in alternative.split(" and "):
            if "=" not in atom:
                raise SymbolicError("boundary predicate %r: each test is 'a = b'" % (atom,))
            lhs, rhs = atom.split("=", 1)
            diff = parse_expression(lhs, space) - parse_expression(rhs, space)
            tests.append("abs(%s) <= tol" % to_code(diff))
        alternatives.append("(" + " and ".join(tests) + ")")
    return " or ".join(alternatives)


def _periodic_axes(where, lower, upper, space: Space, region: str) -> List[int]:
    """The axes whose two opposite faces both lie in a periodic region."""
    if not where:
        raise SymbolicError("periodic region %s needs a where predicate naming its "
                            "two faces, e.g. 'x = 0 or x = 4'" % (region,))
    import math
    test = compile(_where(where, space), "<where>", "eval")
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


def adr_problem(doc: Dict[str, Any]) -> Dict[str, Any]:
    """Build the plain-data ADR problem for a validated (composed) spec."""
    problem = doc["Problem"]
    if "geometry" not in problem:
        raise SymbolicError("Problem.geometry is needed to run a spec (a box)")
    box = problem["geometry"]["box"]
    lower, upper = [float(v) for v in box["lower"]], [float(v) for v in box["upper"]]
    if len(lower) != len(upper):
        raise SymbolicError("geometry.box: lower and upper differ in length")
    dim = len(lower)

    strong = problem["strong_form"]
    if not strong.get("equations"):
        raise SymbolicError("strong_form.equations is needed to run a spec")
    ranks = {}
    for u in strong["unknowns"]:
        if isinstance(u, str):
            ranks[u] = 0
        else:
            ranks[u["name"]] = int(u.get("rank", 0))

    coefficients = _coefficients(strong.get("coefficients") or {}, dim)
    space = Space(dim, ranks, coefficients)

    residuals: List[sympy.Expr] = []
    for text in strong["equations"]:
        residuals.extend(parse_equation(text, space))
    adr = adr_form(residuals, space)

    regions = {r["name"]: r.get("where") for r in strong.get("boundary_regions", [])}
    dirichlet = []
    periodic = []
    for bc in strong.get("boundary_conditions") or []:
        if bc["type"] not in ("dirichlet", "periodic"):
            raise SymbolicError("boundary condition on %s: only dirichlet and periodic "
                                "conditions are runnable so far, not %s" % (bc["region"], bc["type"]))
        name = bc["variable"]
        if name not in ranks:
            raise SymbolicError("boundary condition on %s: %s is not an unknown"
                                % (bc["region"], name))
        if bc["region"] not in regions:
            raise SymbolicError("boundary condition: no boundary region named %r" % bc["region"])
        where = regions[bc["region"]]
        if bc["type"] == "periodic":
            axes = _periodic_axes(where, lower, upper, space, bc["region"])
            comps = [name] if ranks[name] == 0 else ["%s_%d" % (name, i) for i in range(dim)]
            periodic.extend({"component": c, "region": bc["region"], "axes": axes} for c in comps)
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
            dirichlet.append({"component": component, "region": bc["region"],
                              "where": _where(where, space) if where else None,
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
        expressions = entry["solution"].get("expressions")
        if expressions:
            for name, text in expressions.items():
                exact.update(_components(parse_expression(text, space), name,
                                         ranks.get(name, 0), space,
                                         "analytical solution %s" % entry["name"]))
            break

    transient = any("mass" in e for e in adr["equations"])
    unknowns = {name: ([name] if rank == 0 else ["%s_%d" % (name, i) for i in range(dim)])
                for name, rank in ranks.items()}
    return {"adr": adr, "unknowns": unknowns,
            "geometry": {"lower": lower, "upper": upper},
            "dirichlet": dirichlet, "periodic": periodic, "initial": initial, "exact": exact,
            "transient": transient}
