"""The dispersive shallow-water family (Serre-Green-Naghdi) solved by hyperbolic relaxation: the closure, its auxiliary
unknowns and their roles.

C5's dispatcher maps a spec onto proteus's ``GN_SW2DCV``, which solves six unknowns (h, hu, hv, hη, hw, hβ): the model's
depth and discharge, and three auxiliaries the closure introduces. The spec could not say them: the model's family
(``dispersive``, round 3) was a word, with no way to state the closure, the auxiliaries, what each is, or what the
discretization supplies of the model. A sub-problem now states::

    dispersion:
      closure: hyperbolic_relaxation
      closes: [Pd, Qd]                 # the model's coefficients the closure supplies
      auxiliaries:
        - {name: hη, role: relaxed_depth,      units: m2}
        - {name: hw, role: relaxed_divergence, units: m2/s, initial: "0"}        # optional
        - {name: hβ, role: relaxed_bed_flux,   units: m2/s}

* ``closure``: how the dispersive terms are discretized. ``hyperbolic_relaxation`` (Guermond, Kees, Popov, Tovar): the
  elliptic (non-hydrostatic) problem of Serre-Green-Naghdi is replaced by three transport equations that relax toward it
  over a length of the element size, so the whole system is hyperbolic and one explicit scheme solves it.
* ``closes``: the model's coefficients whose values the closure computes from the auxiliaries, by role: the dispersive
  pressure (``dispersive_pressure``, the non-hydrostatic part of the pressure) and the dispersive bed force
  (``dispersive_bed_force``, the topography-dispersion force along ∇Z). Each must be a coefficient of the model with that
  role, so the model's equations are the same whichever closure solves them.
* ``auxiliaries``: one per role of :data:`AUXILIARY_ROLES`, each once. ``name`` is the unknown's name, ``units`` its units,
  ``initial`` (optional) its initial value as a formula in the model's unknowns and coordinates (a solitary wave's hw is
  not 0: proteus's santos_step deck gives one); unstated, the role's
  default (:data:`AUXILIARY_ROLES`) applies (proteus's decks: ``hη = h²``, ``hw = hβ = 0``). Boundary data is the role's
  ``boundary_like`` (hη and hw as h, hβ as q): the auxiliary carries no condition of its own in the spec.

Only a model of family ``dispersive_shallow_water`` may state it (checked on load by ``ymf.vocabulary``).

Stdlib only, like the rest of the spec core.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

__all__ = ["DispersionError", "AUXILIARY_ROLES", "resolve", "check"]


class DispersionError(ValueError):
    """A dispersion block a spec states that does not make sense."""


#: role -> what the auxiliary is, the value it has when the relaxation is exact (``consistent``, as proteus's decks write
#: it: "hη -> h², hβ -> q·∇Z, hw -> h² ∇·u + 3/2 hβ"), its default initial value (the decks', which take hw and hβ as 0:
#: "often okay to take hβ = 0"), and whose boundary data it takes (``h`` or ``q``).
AUXILIARY_ROLES: Dict[str, Dict[str, str]] = {
    "relaxed_depth": {
        "means": "h·η, the depth's relaxed copy",
        "consistent": "h^2", "initial": "h^2", "boundary_like": "h"},
    "relaxed_divergence": {
        "means": "h·w, the relaxed vertical-velocity moment",
        "consistent": "h^2 ∇·u + (3/2) hβ", "initial": "0", "boundary_like": "h"},
    "relaxed_bed_flux": {
        "means": "h·β, the relaxed bed-slope flux",
        "consistent": "q·∇Z", "initial": "0", "boundary_like": "q"},
}

#: the coefficient roles a closure supplies, by closure
CLOSURES: Dict[str, Dict[str, Any]] = {
    "hyperbolic_relaxation": {
        "means": "the non-hydrostatic problem replaced by three transport equations that relax toward it over a length "
                 "of the element size; one explicit hyperbolic system",
        "closes": ("dispersive_pressure", "dispersive_bed_force"),
        "roles": tuple(AUXILIARY_ROLES),
    },
}

_KEYS = {"closure", "closes", "auxiliaries"}
_AUX_KEYS = {"name", "role", "units", "initial"}


def resolve(block: Any, where: str = "dispersion") -> Dict[str, Any]:
    """The block with ymf's words, every auxiliary's default ``initial`` filled in and ``boundary_like`` added:
    ``{closure, closes, auxiliaries: [{name, role, units?, initial, initial_is_default, boundary_like}]}``. Raises
    :class:`DispersionError`."""
    if not isinstance(block, dict) or "closure" not in block:
        raise DispersionError("%s: give closure: (%s)" % (where, ", ".join(CLOSURES)))
    unknown = sorted(set(block) - _KEYS)
    if unknown:
        raise DispersionError("%s: no key %s; its keys are %s" % (where, ", ".join(unknown), ", ".join(sorted(_KEYS))))
    closure = str(block["closure"])
    if closure not in CLOSURES:
        raise DispersionError("%s.closure: %s is not one of %s" % (where, closure, ", ".join(CLOSURES)))
    spec = CLOSURES[closure]
    closes = [str(c) for c in (block.get("closes") or [])]
    if not closes:
        raise DispersionError("%s: closes: names the model's coefficients the closure supplies (%s)" % (
            where, " and ".join(spec["closes"])))
    if len(set(closes)) != len(closes):
        raise DispersionError("%s.closes: a coefficient twice" % (where,))
    auxiliaries = block.get("auxiliaries") or []
    if not isinstance(auxiliaries, list) or not auxiliaries:
        raise DispersionError("%s.auxiliaries: one per role (%s)" % (where, ", ".join(spec["roles"])))
    out: List[Dict[str, Any]] = []
    seen_roles: Dict[str, str] = {}
    names: List[str] = []
    for i, aux in enumerate(auxiliaries):
        here = "%s.auxiliaries[%d]" % (where, i)
        if not isinstance(aux, dict) or "name" not in aux or "role" not in aux:
            raise DispersionError("%s: give name: and role:" % (here,))
        unread = sorted(set(aux) - _AUX_KEYS)
        if unread:
            raise DispersionError("%s: no key %s; its keys are %s" % (here, ", ".join(unread), ", ".join(sorted(_AUX_KEYS))))
        role = str(aux["role"])
        if role not in AUXILIARY_ROLES:
            raise DispersionError("%s.role: %s is not one of %s" % (here, role, ", ".join(AUXILIARY_ROLES)))
        if role not in spec["roles"]:
            raise DispersionError("%s.role: %s is not an auxiliary of %s" % (here, role, closure))
        if role in seen_roles:
            raise DispersionError("%s.role: %s is already the role of %s" % (here, role, seen_roles[role]))
        name = str(aux["name"])
        if name in names:
            raise DispersionError("%s.name: %s is already an auxiliary" % (here, name))
        seen_roles[role] = name
        names.append(name)
        entry: Dict[str, Any] = {"name": name, "role": role}
        if aux.get("units") is not None:
            entry["units"] = str(aux["units"])
        stated = aux.get("initial")
        entry["initial"] = str(stated) if stated is not None else AUXILIARY_ROLES[role]["initial"]
        entry["initial_is_default"] = stated is None
        entry["boundary_like"] = AUXILIARY_ROLES[role]["boundary_like"]
        out.append(entry)
    absent = [r for r in spec["roles"] if r not in seen_roles]
    if absent:
        raise DispersionError("%s.auxiliaries: %s needs one for each of %s; missing %s" % (
            where, closure, ", ".join(spec["roles"]), ", ".join(absent)))
    return {"closure": closure, "closes": closes, "auxiliaries": out}


def check(block: Any, doc: Dict[str, Any], family: Optional[str], where: str = "dispersion") -> Dict[str, Any]:
    """:func:`resolve`, and the block agrees with the model: its family is dispersive, the auxiliaries' names are not
    the model's, and each ``closes`` is a coefficient of the model with the role the closure supplies."""
    out = resolve(block, where)
    if family != "dispersive_shallow_water":
        raise DispersionError("%s: only a model of family dispersive_shallow_water has a dispersion; this one is %s" % (
            where, family or "of no family"))
    strong = ((doc.get("model") or {}).get("strong_form") or {})
    taken = set((strong.get("coefficients") or {}))
    for u in strong.get("unknowns") or []:
        taken.add(u if isinstance(u, str) else u.get("name"))
    for aux in out["auxiliaries"]:
        if aux["name"] in taken:
            raise DispersionError("%s: the auxiliary %s is already a name in the model" % (where, aux["name"]))
    from ymf.vocabulary import coefficient_roles
    by_role = {v: k for k, v in coefficient_roles(doc).items()}
    wanted = list(CLOSURES[out["closure"]]["closes"])
    for name in out["closes"]:
        role = by_role.get(name)
        if role not in wanted:
            raise DispersionError("%s.closes: %s has role %s; %s supplies %s" % (
                where, name, role or "none", out["closure"], " and ".join(wanted)))
    got = sorted(by_role[n] for n in out["closes"])
    if got != sorted(wanted):
        raise DispersionError("%s.closes: %s supplies %s, and these are %s" % (
            where, out["closure"], " and ".join(wanted), " and ".join(got)))
    return out
