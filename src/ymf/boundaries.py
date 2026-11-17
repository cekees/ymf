"""Boundary kind: how a boundary condition reaches the discrete problem (strong, weak, none), and a weak one's parameters.

A condition says *what* holds on a boundary (a value, a flux). Its ``imposed`` says *how the discrete problem enforces
it*, which is what picks the boundary term a kernel evaluates there (the engine's face kinds ``Strong``,
``WeakDirichlet`` and ``Flux``; C19's gap note G8)::

    - {region: inlet,  variable: φ, type: dirichlet, only: inflow, imposed: strong}   # rows of the system replaced
    - {region: inlet,  variable: φ, type: dirichlet, only: inflow, imposed: weak,     # penalty and consistency terms
       penalty: 100.0, form: symmetric}
    - {region: outlet, variable: φ, type: neumann,   value: 0.0}                      # none: the natural condition

* ``strong``: the unknown's degrees of freedom on the boundary are set to the value (a Dirichlet condition).
* ``weak``: the value is enforced through terms of the weak form (Nitsche): a penalty ``penalty / h`` and, by ``form``, the
  adjoint consistency term: ``symmetric`` (SIPG), ``nonsymmetric`` (NIPG) or ``incomplete`` (IIPG, none). Dirichlet only.
* ``none``: nothing constrains the solution on the boundary; the boundary integral carries the given flux (Neumann, Robin)
  or leaves with the solution's own trace (``open``).

The kind is what a condition's *type* allows: a Dirichlet condition (and the Dirichlet part of ``slip``, ``inflow``,
``characteristic``) is ``strong`` or ``weak`` (default ``strong``); Neumann, Robin and ``open`` are ``none``;
``periodic`` is none of them (the faces are identified). :func:`kind_of` gives the kind of a condition, stated or the
default; :func:`check_condition` refuses what does not make sense (a strong Neumann condition, a penalty on a strong
one); :func:`kinds` lists every condition of a document with its kind.

When an alias is expanded (``ymf.conditions``) the kind goes to the Dirichlet components only: ``slip, imposed: weak`` is
a weak normal condition and a natural tangential one. ``split_problem`` gives each condition's ``kind``, and for a weak
one its ``penalty`` and ``form``, so a generator picks the face kernel from the declaration and not from the model's name.

Stdlib only, like the rest of the spec core.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

__all__ = ["BoundaryError", "BOUNDARY_KINDS", "WEAK_FORMS", "kind_of", "check_condition", "kinds"]


class BoundaryError(ValueError):
    """A boundary kind a spec states that does not make sense."""


BOUNDARY_KINDS: Dict[str, str] = {
    "strong": "the degrees of freedom on the boundary are set to the value",
    "weak": "the value is enforced through penalty and consistency terms of the weak form (Nitsche)",
    "none": "no constraint on the solution: the boundary integral carries the given flux, or the solution's own trace",
}

WEAK_FORMS: Dict[str, str] = {
    "symmetric": "SIPG: the adjoint consistency term with the sign that keeps the form symmetric",
    "nonsymmetric": "NIPG: the adjoint consistency term with the opposite sign",
    "incomplete": "IIPG: no adjoint consistency term",
}

#: the types a condition's kind is chosen for (the Dirichlet part), and the types that are always ``none``
_DIRICHLET = ("dirichlet", "slip", "inflow", "characteristic")
_NATURAL = ("neumann", "robin", "open")


def kind_of(condition: Dict[str, Any]) -> Optional[str]:
    """The boundary kind of a condition (an alias, or an abstract one): what it states, or its type's default. None for
    ``periodic``."""
    kind = condition.get("type")
    if kind == "periodic":
        return None
    stated = condition.get("imposed")
    if stated is not None:
        return str(stated)
    return "strong" if kind in _DIRICHLET else "none"


def check_condition(condition: Dict[str, Any], where: str) -> None:
    """``imposed``, ``penalty`` and ``form`` of a condition agree with its type. Raises :class:`BoundaryError`."""
    kind = condition.get("type")
    stated = condition.get("imposed")
    has_weak_keys = [k for k in ("penalty", "form") if condition.get(k) is not None]
    if stated is not None and str(stated) not in BOUNDARY_KINDS:
        raise BoundaryError("%s: imposed %s is not one of %s" % (where, stated, ", ".join(BOUNDARY_KINDS)))
    if kind == "periodic":
        if stated is not None or has_weak_keys:
            raise BoundaryError("%s: a periodic condition identifies faces; it has no kind" % (where,))
        return
    if kind in _NATURAL and stated is not None and str(stated) != "none":
        raise BoundaryError(
            "%s: a %s condition is imposed none (the boundary integral carries it); %s is for a Dirichlet condition"
            % (where, kind, stated))
    if kind in _DIRICHLET and stated is not None and str(stated) == "none":
        raise BoundaryError(
            "%s: a %s condition imposes a value, strong or weak; to impose nothing give no condition, or open"
            % (where, kind))
    weak = kind_of(condition) == "weak"
    if has_weak_keys and not weak:
        raise BoundaryError("%s: %s belong to a weak condition (imposed: weak); this one is %s" % (
            where, " and ".join(has_weak_keys), kind_of(condition)))
    if condition.get("penalty") is not None:
        try:
            ok = float(condition["penalty"]) > 0
        except (TypeError, ValueError):
            ok = False
        if not ok:
            raise BoundaryError("%s: penalty %s is not a positive number" % (where, condition["penalty"]))
    if condition.get("form") is not None and str(condition["form"]) not in WEAK_FORMS:
        raise BoundaryError("%s: form %s is not one of %s" % (where, condition["form"], ", ".join(WEAK_FORMS)))


def kinds(doc: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Every boundary condition of the document, expanded, as ``{region, variable, type, [component], [only], kind,
    [penalty], [form]}`` in order: what a solver reads to choose its boundary terms."""
    from ymf.conditions import expand
    from ymf.normalize import boundary_conditions
    strong = ((expand(doc).get("model") or {}).get("strong_form") or {})
    out = []
    for bc in boundary_conditions(strong):
        entry = {k: bc[k] for k in ("region", "variable", "type", "component", "only") if bc.get(k) is not None}
        entry["kind"] = kind_of(bc)
        for k in ("penalty", "form"):
            if bc.get(k) is not None:
                entry[k] = bc[k]
        out.append(entry)
    return out
