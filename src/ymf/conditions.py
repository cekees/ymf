"""Boundary-condition aliases, and their expansion into the abstract forms.

Dirichlet, Neumann and Robin conditions (and periodic pairing) are the
abstract forms; each is complete with its expression. The other kinds a spec
may write -- ``slip``, ``open``, ``inflow``, ``characteristic`` -- and the
named condition sets of a model (``free_slip``, ``atmosphere``) are
shorthand for particular combinations of them. A spec may use the shorthand;
what is solved, and what a closed spec and its archive record, is the
expansion (cekees, 2026-11-04): otherwise two solvers could close the same
spec differently and the archive could not show which, and the run would not
be reproducible from its record.

An expanded condition is one of the abstract forms, on

* the whole variable, or a ``component`` of a vector: its ``normal`` or
  ``tangential`` part on the boundary, or a Cartesian component ``x``, ``y``,
  ``z``;
* the whole region, or a ``part`` of it given by a predicate (``x = 0 and
  0.4 ≤ y ≤ 0.6``: a state imposed on part of a boundary), and ``only`` the
  part of that where the flow enters (``inflow``) or where characteristics
  enter (``incoming_characteristics``);

with its value (``value``, ``formula`` or ``series``), optionally how it is
``imposed`` (``strong`` or ``weak``: C12 makes that a property of each
condition), and ``expanded_from``, the alias and the condition set it came
from.

:data:`ALIASES` is the definition, as data. :func:`expand` applies it to a
document; it is idempotent, and a document with no alias and no condition set
comes back unchanged, so the outputs it keys keep their keys.

**A part, and the rest of its region.** A condition with a ``part`` holds
on the points of its region where the predicate holds, and there it replaces
the conditions on the same region and variable that have no part and are on
a component that overlaps its own (:func:`overlap`: the whole variable
overlaps every component; ``x`` and ``y``, or ``normal`` and ``tangential``,
do not): those hold on the rest of the region. So ``boundary_sets: [{region: all, set: wall}]``
and ``{region: all, part: "x = 0", variable: q, type: dirichlet, formula:
"(0.75, 0)"}`` is a wall everywhere but on x = 0, where the discharge is
given; ``{region: all, variable: q, type: dirichlet, component: y, value:
0}`` with ``{region: all, part: "x = 0", variable: q, type: dirichlet,
component: x, value: 0.18}`` holds q's y component everywhere, x = 0
included. A part is not a region of its own, so a condition with a part does
not replace a condition set's member either: the member is the rest.
:func:`excluded` gives the parts a condition without one does not hold on.

This module is stdlib only, like the rest of the spec core.
"""

from __future__ import annotations

import copy
from typing import Any, Dict, List, Tuple

__all__ = ["ALIASES", "ABSTRACT", "AXES", "expand", "expand_condition", "has_shorthand", "describe", "parts",
           "excluded", "overlap"]

#: the Cartesian components a condition may be on, and the index of each
AXES = {"x": 0, "y": 1, "z": 2}

#: the abstract forms: what an expansion is made of
ABSTRACT = ("dirichlet", "neumann", "robin", "periodic")

#: the keys that carry a condition's value
_VALUE_KEYS = ("value", "formula", "series", "units")

#: the order an expanded condition's keys are written in
_ORDER = ("region", "part", "variable", "type", "component", "only", "value", "formula", "series",
          "alpha", "beta", "units", "imposed", "penalty", "form", "source", "source_url", "temporal",
          "expanded_from")


def _ordered(entry: Dict[str, Any]) -> Dict[str, Any]:
    head = {k: entry[k] for k in _ORDER if k in entry}
    head.update({k: v for k, v in entry.items() if k not in head})
    return head

#: alias -> what it means. Each component is an abstract condition; ``takes``
#: says it receives the alias's value (value, formula or series), and
#: ``default`` what it is when the alias gives none. ``needs`` restricts the
#: alias to a vector variable (a normal and a tangential part).
ALIASES: Dict[str, Dict[str, Any]] = {
    # a free-slip wall: no flow through it, no tangential traction on it
    "slip": {
        "needs": "vector",
        "means": "the normal component is the value (0 if none is given) and the "
                 "tangential traction (the diffusive flux of the tangential part) is zero",
        "components": [
            {"type": "dirichlet", "component": "normal", "takes": True, "default": 0.0},
            {"type": "neumann", "component": "tangential", "value": 0.0},
        ],
    },
    # a declared open boundary (C12): the natural condition, stated
    "open": {
        "means": "zero diffusive flux, the natural condition, declared: the advective "
                 "flux leaves with the solution's own trace",
        "components": [
            {"type": "neumann", "value": 0.0},
        ],
    },
    # a value imposed where the flow comes in, nothing where it leaves
    "inflow": {
        "means": "the value, imposed only where the flow enters",
        "components": [
            {"type": "dirichlet", "only": "inflow", "takes": True},
        ],
    },
    # a state imposed through the characteristics that enter (a tailwater)
    "characteristic": {
        "means": "the value, imposed only through the incoming characteristics (so a "
                 "supercritical outflow ignores it)",
        "components": [
            {"type": "dirichlet", "only": "incoming_characteristics", "takes": True},
        ],
    },
}


def expand_condition(bc: Dict[str, Any], rank: int = 0) -> List[Dict[str, Any]]:
    """One condition as the abstract conditions it means.

    An abstract condition comes back as itself; an alias as its components,
    each with the condition's region and variable, its value where the
    component takes one, and ``expanded_from: {type: <alias>}`` (and the set,
    when the condition came from one). ``rank`` is the variable's (1 for a
    vector), which an alias such as ``slip`` needs.
    """
    kind = bc.get("type")
    origin = {}
    if bc.get("set") is not None:
        origin["set"] = bc["set"]
    if kind in ABSTRACT or kind not in ALIASES:
        out = {k: v for k, v in bc.items() if k != "set"}
        if origin:
            out["expanded_from"] = dict(out.get("expanded_from") or {}, **origin)
        return [out]
    alias = ALIASES[kind]
    if alias.get("needs") == "vector" and rank != 1:
        raise ValueError("a %s condition is on a vector, and %s is not one"
                         % (kind, bc.get("variable")))
    given = {k: bc[k] for k in _VALUE_KEYS if bc.get(k) is not None}
    shared = {k: v for k, v in bc.items()
              if k not in _VALUE_KEYS and k not in ("type", "set", "expanded_from")}
    # how a condition is imposed is the Dirichlet part's (ymf.boundaries): the natural parts stay natural
    kind_keys = {k: shared.pop(k) for k in ("imposed", "penalty", "form") if k in shared}
    out = []
    for part in alias["components"]:
        entry = dict(shared)
        entry.update({k: v for k, v in part.items() if k not in ("takes", "default")})
        if part["type"] == "dirichlet":
            entry.update(kind_keys)
        if part.get("takes"):
            if given:
                entry.update(given)
            elif "default" in part:
                entry["value"] = part["default"]
        entry["expanded_from"] = dict({"type": kind}, **origin)
        out.append(entry)
    return out


def has_shorthand(strong_form: Dict[str, Any]) -> bool:
    """Whether a strong form uses an alias or a condition set."""
    if strong_form.get("boundary_sets") or strong_form.get("condition_sets"):
        return True
    return any(isinstance(bc, dict) and bc.get("type") in ALIASES
               for bc in strong_form.get("boundary_conditions") or [])


def expand(doc: Dict[str, Any]) -> Dict[str, Any]:
    """``doc`` with every alias and condition set expanded (see the module).

    The strong form's ``boundary_conditions`` become the abstract conditions,
    in order (written ones first, then the applied sets', each alias replaced
    by its components), and ``condition_sets`` and ``boundary_sets`` go, their
    meaning now written out and their names kept in ``expanded_from``. A
    document without either comes back as the same object, unchanged.
    Raises ``ValueError`` for what does not expand (an undefined set, two sets
    deciding one region's variable, a vector alias on a scalar).
    """
    model = doc.get("model") if isinstance(doc, dict) else None
    strong = (model or {}).get("strong_form") if isinstance(model, dict) else None
    if not isinstance(strong, dict) or not has_shorthand(strong):
        return doc
    from ymf.normalize import boundary_conditions, normalize_unknowns
    ranks: Dict[str, int] = {}
    for u in strong.get("unknowns") or []:
        if isinstance(u, dict):
            ranks[u["name"]] = int(u.get("rank") or 0)
    for u in normalize_unknowns(strong.get("unknowns") or []):
        ranks.setdefault(u["name"], 0)
    expanded: List[Dict[str, Any]] = []
    for bc in boundary_conditions(strong):
        expanded.extend(_ordered(e) for e in
                        expand_condition(bc, ranks.get(bc.get("variable"), 0)))
    out = copy.copy(doc)
    out["model"] = dict(model)
    new_strong = {k: v for k, v in strong.items()
                  if k not in ("condition_sets", "boundary_sets", "boundary_conditions")}
    # keep the strong form's key order: the conditions where they were
    ordered: Dict[str, Any] = {}
    for k in strong:
        if k == "boundary_conditions" or (k == "boundary_sets"
                                          and "boundary_conditions" not in strong):
            ordered["boundary_conditions"] = expanded
        elif k in new_strong:
            ordered[k] = new_strong[k]
    ordered.setdefault("boundary_conditions", expanded)
    out["model"]["strong_form"] = ordered
    return out


def describe() -> List[Tuple[str, str, List[Dict[str, Any]]]]:
    """``(alias, what it means, its components)`` for the documentation."""
    return [(name, a["means"], copy.deepcopy(a["components"])) for name, a in ALIASES.items()]



def overlap(a: Any, b: Any) -> bool:
    """Whether conditions on components ``a`` and ``b`` of one variable (None: the whole variable) can constrain the
    same thing: not when they are two different Cartesian components, nor the normal and the tangential part."""
    if a is None or b is None or a == b:
        return True
    if a in AXES and b in AXES:
        return False
    return {a, b} != {"normal", "tangential"}


def parts(conditions: List[Dict[str, Any]]) -> Dict[Tuple[Any, Any], List[Tuple[str, Any]]]:
    """``{(region, variable): [(part, component), ...]}``: the parts of each region that conditions on a variable claim,
    with the component each is on, in the order written, each once."""
    out: Dict[Tuple[Any, Any], List[Tuple[str, Any]]] = {}
    for bc in conditions:
        if isinstance(bc, dict) and bc.get("part") is not None:
            claimed = out.setdefault((bc.get("region"), bc.get("variable")), [])
            if (bc["part"], bc.get("component")) not in claimed:
                claimed.append((bc["part"], bc.get("component")))
    return out


def excluded(bc: Dict[str, Any], claimed: Dict[Tuple[Any, Any], List[Tuple[str, Any]]]) -> List[str]:
    """The parts a condition without a part does not hold on: those of its region that conditions on the same variable,
    on a component that overlaps its own, claim (``claimed`` is :func:`parts`). Empty for a condition with a part."""
    if bc.get("part") is not None:
        return []
    out: List[str] = []
    for part, component in claimed.get((bc.get("region"), bc.get("variable")), []):
        if overlap(component, bc.get("component")) and part not in out:
            out.append(part)
    return out
