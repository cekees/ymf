"""The words a spec states and a solver maps: ymf's vocabulary for model
families, sub-problem roles, discretization options and coefficient roles.

A dispatcher (proteus's ``DispatchTable``) maps a spec onto modules and
flags. It used to coin the words it matched on -- a sub-problem's role, a
stabilization method, a time integration -- which meant the spec could not
say them and a second solver could not read them. They live here instead
(C42 round 3, cekees 2026-11-06, decision 16): the spec states them, ymf
checks them on load, the closed spec writes the ones it derives, and a
solver's table only maps them.

Each kind is a mapping ``{word: meaning}``. Where the dispatch table's word
was a solver's name for the thing (``vbdf``, ``SSP_RK3``),
ymf's word replaces it and the old one is accepted as an alias
(:data:`ALIASES`): a spec using it still loads, and :func:`canonical` gives
ymf's word. Aliases are not rewritten on load or in the closed spec yet, so
that a reader keyed on the old word keeps working until it maps the new one.

What is here:

* :data:`FAMILIES` -- model families, by the ``physical_model`` processes and
  assumptions that identify them; :func:`family_of`.
* :data:`ROLES` -- what a sub-problem of a split discretization does, by
  family; :func:`role_of` reads it from the sub-problem's structure.
* :data:`TIME_INTEGRATION`, :data:`STABILIZATION`, :data:`SCHEMES`,
  :data:`STEP_CONTROL`, :data:`REGULARIZATIONS`, :data:`BINDING_MEANINGS`,
  :data:`WRITES_BACK`, :data:`INERTIA_FORMS` -- the discretization's words.
* :data:`LINEAR_SOLVERS` -- ``direct`` and ``iterative``, and the packages a
  spec may still name; :func:`linear_solver_kind`.

:func:`check` is called on load (ymf.references). Stdlib only, like the rest
of the spec core.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Iterable, List, Optional

__all__ = ["FAMILIES", "ROLES", "TIME_INTEGRATION", "STABILIZATION", "SCHEMES",
           "STEP_CONTROL", "REGULARIZATIONS", "BINDING_MEANINGS", "WRITES_BACK",
           "INERTIA_FORMS", "LINEAR_SOLVERS", "ALIASES",
           "canonical", "family_of", "role_of", "roles_of", "linear_solver_kind",
           "check", "VocabularyError"]


class VocabularyError(ValueError):
    """A word the vocabulary does not have, or one used where it cannot be."""


# ---------------------------------------------------------------------------
# Model families
# ---------------------------------------------------------------------------

#: a model belongs to a family when its physical_model lists every process
#: and assumption the family names
FAMILIES: Dict[str, Dict[str, Any]] = {
    "two_phase_level_set": {
        "processes": ["incompressible_flow", "two_phase_flow", "interface_transport"],
        "means": "incompressible two-phase flow with a sharp interface carried by a level set",
    },
    "shallow_water": {
        "assumptions": ["hydrostatic", "depth_averaged", "nondispersive"],
        "means": "the hydrostatic shallow-water equations (Saint-Venant)",
    },
    # The dispersive model (Serre-Green-Naghdi) is depth averaged and not
    # hydrostatic. The dispatch table keyed it on "hyperbolic_relaxation",
    # which is how GN_SW2DCV discretizes it, not an assumption of the model;
    # that spelling is still recognized (FAMILY_ALIASES).
    "dispersive_shallow_water": {
        "assumptions": ["depth_averaged", "dispersive"],
        "means": "the dispersive (Serre-Green-Naghdi) shallow-water equations",
    },
}

#: older ways of identifying a family, still recognized
FAMILY_ALIASES: Dict[str, List[Dict[str, List[str]]]] = {
    "dispersive_shallow_water": [{"assumptions": ["depth_averaged", "hyperbolic_relaxation"]}],
}


def _physical(doc: Dict[str, Any]) -> Dict[str, Any]:
    return ((doc.get("model") or {}).get("physical_model") or {})


def _matches(needs: Dict[str, Any], have: Dict[str, set]) -> bool:
    return all(set(needs.get(key) or []) <= have[key] for key in ("processes", "assumptions"))


def families_of(doc: Dict[str, Any]) -> List[str]:
    """Every family the model's physical_model identifies, in table order."""
    physical = _physical(doc)
    have = {"processes": set(physical.get("processes") or []),
            "assumptions": set(physical.get("assumptions") or [])}
    return [name for name, needs in FAMILIES.items()
            if _matches(needs, have)
            or any(_matches(alt, have) for alt in FAMILY_ALIASES.get(name, []))]


def family_of(doc: Dict[str, Any]) -> Optional[str]:
    """The model's family: the one it states (``physical_model.family``), or
    the one its processes and assumptions identify; None if none does.
    Raises :class:`VocabularyError` if several do and the model states none."""
    stated = _physical(doc).get("family")
    if stated is not None:
        return str(stated)
    found = families_of(doc)
    if len(found) > 1:
        raise VocabularyError(
            "model.physical_model identifies the families %s; say which (family: ...)"
            % ", ".join(found))
    return found[0] if found else None


# ---------------------------------------------------------------------------
# Sub-problem roles
# ---------------------------------------------------------------------------

#: what a sub-problem does, by family. role_of reads it from what the
#: sub-problem solves, enforces, introduces and spins up.
ROLES: Dict[str, Dict[str, str]] = {
    "two_phase_level_set": {
        "velocity_pressure": "solves the equations for the velocity and the pressure together",
        "velocity": "solves the momentum equation alone (a projection scheme's velocity step)",
        "pressure_increment": "enforces a divergence-free velocity (the projection)",
        "pressure": "solves the continuity equation for the pressure alone (the pressure update)",
        "pressure_spin_up": "a spin-up solve for the initial pressure",
        "volume_fraction": "introduces a conservative carrier of heaviside(level set) and transports it",
        "level_set": "solves the interface equation for the level set",
        "level_set_redistanced": "solves the interface equation and enforces |∇(level set)| = 1 "
                                 "in one solve",
        "redistancing": "enforces |∇(level set)| = 1 to steady state, writing the level set back",
        "mass_correction": "introduces a correction of the level set so that heaviside(level set) "
                           "matches the carrier",
    },
    "shallow_water": {
        "shallow_water": "solves the mass and momentum equations together",
    },
    "dispersive_shallow_water": {
        "shallow_water": "solves the mass and momentum equations together",
    },
}

_HEAVISIDE = re.compile(r"^\s*heaviside\(\s*([^\s()]+)\s*\)\s*$")


def _strong(doc: Dict[str, Any]) -> Dict[str, Any]:
    return ((doc.get("model") or {}).get("strong_form") or {})


def _ranks(strong: Dict[str, Any]) -> Dict[str, int]:
    out = {}
    for u in strong.get("unknowns") or []:
        if isinstance(u, str):
            out[u] = 0
        elif isinstance(u, dict):
            out[u["name"]] = int(u.get("rank") or 0)
    return out


def heaviside_of(doc: Dict[str, Any]) -> Dict[str, str]:
    """``{coefficient: unknown}`` for each coefficient that is heaviside of an
    unknown (``χ := heaviside(φ)`` gives ``{χ: φ}``): a phase indicator."""
    strong = _strong(doc)
    ranks = _ranks(strong)
    out = {}
    for name, entry in (strong.get("coefficients") or {}).items():
        formula = entry.get("formula") if isinstance(entry, dict) else entry
        m = _HEAVISIDE.match(str(formula or ""))
        if m and m.group(1) in ranks:
            out[name] = m.group(1)
    return out


def _normal(text: Any) -> str:
    return re.sub(r"\s+", "", str(text or ""))


def _is_eikonal(enforces: Any, level_sets: Iterable[str]) -> bool:
    text = _normal(enforces)
    return any(text in ("norm(∇%s)=1" % phi, "|∇%s|=1" % phi) for phi in level_sets)


def _is_divergence_free(enforces: Any, ranks: Dict[str, int]) -> bool:
    text = _normal(enforces)
    return any(rank == 1 and text == "∇·%s=0" % name for name, rank in ranks.items())


def role_of(sub: Dict[str, Any], doc: Dict[str, Any],
            family: Optional[str] = None) -> Optional[str]:
    """What one sub-problem does (a word of :data:`ROLES`), read from its
    structure; None if no role of the family fits. ``family`` defaults to
    :func:`family_of`. A role the sub-problem states is not consulted: this
    is what its structure says (:func:`check` compares the two)."""
    family = family if family is not None else family_of(doc)
    strong = _strong(doc)
    ranks = _ranks(strong)
    equations = {e.get("name"): e for e in strong.get("equations") or [] if isinstance(e, dict)}
    solved = [equations[e].get("for") for e in sub.get("solves") or [] if e in equations]
    vectors = [u for u in solved if ranks.get(u) == 1]
    scalars = [u for u in solved if ranks.get(u) == 0]
    phase = heaviside_of(doc)
    level_sets = set(phase.values())
    new = sub.get("introduces") or {}
    enforces = sub.get("enforces")
    if family == "two_phase_level_set":
        if solved:
            if len(vectors) == 1 and len(scalars) == 1 and not set(scalars) & level_sets:
                return "velocity_pressure"
            if len(vectors) == 1 and not scalars:
                return "velocity"
            if not vectors and len(scalars) == 1 and scalars[0] in level_sets:
                return "level_set_redistanced" if _is_eikonal(enforces, level_sets) else "level_set"
            if not vectors and len(scalars) == 1:
                return "pressure"
            return None
        if new.get("approximates") in phase:
            return "volume_fraction"
        if new.get("correction_of") in level_sets:
            return "mass_correction"
        if _is_eikonal(enforces, level_sets):
            return "redistancing"
        if _is_divergence_free(enforces, ranks):
            return "pressure_increment"
        if sub.get("spin_up"):
            return "pressure_spin_up"
        return None
    if family in ("shallow_water", "dispersive_shallow_water"):
        if solved and set(solved) == set(ranks):
            return "shallow_water"
    return None


def roles_of(doc: Dict[str, Any], entry: Dict[str, Any]) -> Dict[str, Optional[str]]:
    """``{sub-problem: role}`` for one split discretization entry."""
    family = family_of(doc)
    return {s["name"]: role_of(s, doc, family) for s in entry.get("subproblems") or []}


# ---------------------------------------------------------------------------
# The discretization's words
# ---------------------------------------------------------------------------

TIME_INTEGRATION: Dict[str, str] = {
    "backward_euler": "backward (implicit) Euler, first order",
    "bdf2": "the second-order backward differentiation formula, variable step",
    "ssp_rk3": "the three-stage, third-order strong-stability-preserving Runge-Kutta method",
    "none": "the sub-problem has no time derivative (a constraint, a projection)",
}

STABILIZATION: Dict[str, str] = {
    "ASGS": "algebraic subgrid scale stabilization (variational multiscale)",
    "SUPG": "streamline-upwind Petrov-Galerkin",
    "residual_shock_capturing": "residual-based shock-capturing viscosity alone",
    "entropy_viscosity": "entropy viscosity (Guermond, Pasquetti, Popov)",
    "graph_viscosity_smoothness": "graph viscosity scaled by a smoothness indicator",
    "graph_viscosity_entropy_residual": "graph viscosity scaled by the entropy residual",
}

#: an algebraic scheme stated in a discretization's ``discretization:`` block:
#: its name, and the options it reads (option -> meaning). ``time_integration``,
#: ``outputs`` and ``dt`` are read for every scheme.
SCHEMES: Dict[str, Dict[str, Any]] = {
    "graph_viscosity_convex_limiting": {
        "means": "low-order graph viscosity with convex limiting toward the high-order "
                 "(entropy-viscosity) solution",
        "options": {
            "entropy_viscosity_coefficient": "the high-order solution's entropy-viscosity "
                                             "coefficient (cE)",
            "lumped_mass": "a lumped (true) or consistent (false) mass matrix",
            "bounds_relaxation": "relax the convex-limiting bounds (true) or not",
        },
    },
}

#: read in every scheme block
SCHEME_COMMON = {"scheme": "the scheme's name", "time_integration": "a TIME_INTEGRATION word",
                 "outputs": "the number of output intervals", "dt": "the time step"}

STEP_CONTROL: Dict[str, str] = {
    "min_cfl": "each sub-problem proposes a step from its CFL number; the system takes the "
               "smallest",
    "fixed": "a fixed step, dt",
}

REGULARIZATIONS: Dict[str, str] = {
    "smoothed_heaviside": "H_ε(φ) = (1 + φ/ε + sin(πφ/ε)/π)/2 on |φ| ≤ ε",
    "smoothed_sign": "2 H_ε − 1",
    "wet_dry_cutoff": "1/h -> 2h/(h² + max(h, hε)²)",
}

BINDING_MEANINGS: Dict[str, str] = {
    "element_size": "the size of each element of the mesh",
    "dry_depth": "the depth below which a node is dry (set by the solver)",
}

WRITES_BACK: Dict[str, str] = {
    "replace": "the sub-problem's solution replaces the field",
    "add": "the sub-problem's solution is added to the field (a correction)",
    "mass_quadrature": "the field's values at the quadrature points are reset from it",
}

#: weak-form labels with a meaning dispatchers share: the inertia form of a
#: momentum equation. Both are formulations of the same strong form when the
#: velocity is divergence free.
INERTIA_FORMS: Dict[str, str] = {
    "conservative_inertia": "∂(ρv)/∂t + ∇·(ρ v⊗v)",
    "nonconservative_inertia": "ρ (∂v/∂t + v·∇v)",
}

# ---------------------------------------------------------------------------
# Linear solvers
# ---------------------------------------------------------------------------

#: ``direct`` and ``iterative`` say what is asked of the linear solve; each
#: solver maps them to its own options and refuses what it cannot honour. A
#: package name says which implementation, and is a request of its kind.
LINEAR_SOLVERS: Dict[str, Dict[str, Any]] = {
    "direct": {"kind": "direct", "means": "a direct factorization; the solver picks the package "
                                          "(serial or distributed)"},
    "iterative": {"kind": "iterative", "means": "a Krylov method with a preconditioner; "
                                                "petsc_options refine it"},
    "superlu": {"kind": "direct", "means": "SuperLU (SuperLU_DIST in parallel)"},
    "mumps": {"kind": "direct", "means": "MUMPS"},
    "umfpack": {"kind": "direct", "means": "UMFPACK (serial)"},
    "petsc": {"kind": "iterative", "means": "a PETSc KSP; direct if petsc_options ask for "
                                            "pc_type lu or cholesky with ksp_type preonly"},
}

_DIRECT_PC = ("lu", "cholesky")


def linear_solver_kind(solver: Optional[Dict[str, Any]]) -> Optional[str]:
    """``direct``, ``iterative``, or None when the solver block names no
    linear solver."""
    if not solver or solver.get("linear_solver") is None:
        return None
    word = str(solver["linear_solver"])
    if word not in LINEAR_SOLVERS:
        raise VocabularyError("linear_solver %s: ymf knows %s" % (word, ", ".join(LINEAR_SOLVERS)))
    if word == "petsc":
        options = solver.get("petsc_options") or {}
        if (str(options.get("pc_type")) in _DIRECT_PC
                and str(options.get("ksp_type", "preonly")) == "preonly"):
            return "direct"
    return LINEAR_SOLVERS[word]["kind"]


#: the words of each kind, for check() and canonical()
_KINDS: Dict[str, Dict[str, Any]] = {
    "family": FAMILIES, "role": {r: m for fam in ROLES.values() for r, m in fam.items()},
    "time_integration": TIME_INTEGRATION, "stabilization": STABILIZATION, "scheme": SCHEMES,
    "step_control": STEP_CONTROL, "regularization": REGULARIZATIONS,
    "binding_meaning": BINDING_MEANINGS, "writes_back": WRITES_BACK,
    "linear_solver": LINEAR_SOLVERS,
    "scheme_option": {o: m for s in SCHEMES.values() for o, m in s["options"].items()},
}

#: kind -> {older word: ymf's word}. The older words are the dispatch table's
#: (proteus's class or flag names); a spec may still use them.
ALIASES: Dict[str, Dict[str, str]] = {
    "time_integration": {"vbdf": "bdf2", "SSP_RK3": "ssp_rk3"},
    "stabilization": {"shock_capturing": "residual_shock_capturing"},
    "role": {"pressure_initial": "pressure_spin_up"},
    "scheme_option": {"entropy_viscosity_cE": "entropy_viscosity_coefficient"},
}


def canonical(kind: str, word: Any) -> str:
    """ymf's word for ``word`` of ``kind`` (itself, or what its alias means);
    raises :class:`VocabularyError` if ymf does not know it."""
    words = _KINDS[kind]
    word = str(word)
    if word in words:
        return word
    alias = ALIASES.get(kind, {}).get(word)
    if alias is not None:
        return alias
    raise VocabularyError("%s %s: ymf knows %s" % (
        kind.replace("_", " "), word,
        ", ".join(list(words) + ["%s (= %s)" % kv for kv in ALIASES.get(kind, {}).items()])))


# ---------------------------------------------------------------------------
# The check on load
# ---------------------------------------------------------------------------

def _word(kind: str, word: Any, where: str) -> str:
    try:
        return canonical(kind, word)
    except VocabularyError as exc:
        raise VocabularyError("%s: %s" % (where, exc)) from None


def _stabilization_method(s: Any) -> Any:
    if isinstance(s, str):
        return s
    if isinstance(s, dict):
        return s.get("method")
    return None


def _check_solver(solver: Any, where: str) -> None:
    if isinstance(solver, dict) and solver.get("linear_solver") is not None:
        _word("linear_solver", solver["linear_solver"], where + ".solver.linear_solver")


def _check_entry(doc: Dict[str, Any], entry: Dict[str, Any], family: Optional[str]) -> None:
    where = "solution_paths.discretizations[%s]" % (entry.get("name", "?"),)
    _check_solver(entry.get("solver"), where)
    step = entry.get("step_control")
    if isinstance(step, dict) and step.get("type") is not None:
        _word("step_control", step["type"], where + ".step_control.type")
    block = entry.get("discretization")
    if isinstance(block, dict):
        if block.get("time_integration") is not None:
            _word("time_integration", block["time_integration"],
                  where + ".discretization.time_integration")
        if block.get("scheme") is not None:
            scheme = _word("scheme", block["scheme"], where + ".discretization.scheme")
            known = set(SCHEME_COMMON) | set(SCHEMES[scheme]["options"]) | {
                a for a, w in ALIASES["scheme_option"].items() if w in SCHEMES[scheme]["options"]}
            unread = sorted(set(block) - known)
            if unread:
                raise VocabularyError(
                    "%s.discretization: the scheme %s has no option %s; its options are %s"
                    % (where, scheme, ", ".join(unread),
                       ", ".join(list(SCHEMES[scheme]["options"]) + list(SCHEME_COMMON))))
    for label, binding in list((entry.get("binds") or {}).items()):
        if isinstance(binding, dict) and binding.get("meaning") is not None:
            _word("binding_meaning", binding["meaning"], "%s.binds.%s.meaning" % (where, label))
    for function, rule in (entry.get("regularize") or {}).items():
        if isinstance(rule, dict) and rule.get("by") is not None:
            _word("regularization", rule["by"], "%s.regularize.%s.by" % (where, function))
    roles = ROLES.get(family or "", {})
    for sub in entry.get("subproblems") or []:
        here = "%s.subproblems[%s]" % (where, sub.get("name", "?"))
        if sub.get("time_integration") is not None:
            _word("time_integration", sub["time_integration"], here + ".time_integration")
        method = _stabilization_method(sub.get("stabilization"))
        if method is not None:
            _word("stabilization", method, here + ".stabilization.method")
        _check_solver(sub.get("solver"), here)
        for label, binding in (sub.get("binds") or {}).items():
            if isinstance(binding, dict) and binding.get("meaning") is not None:
                _word("binding_meaning", binding["meaning"], "%s.binds.%s.meaning" % (here, label))
        for function, rule in (sub.get("regularize") or {}).items():
            if isinstance(rule, dict) and rule.get("by") is not None:
                _word("regularization", rule["by"], "%s.regularize.%s.by" % (here, function))
        for field, how in (sub.get("writes_back") or {}).items():
            _word("writes_back", how, "%s.writes_back.%s" % (here, field))
        if sub.get("role") is not None:
            stated = _word("role", sub["role"], here + ".role")
            if family is None:
                raise VocabularyError("%s.role: %s is a role of a model family, and the model "
                                      "belongs to none" % (here, stated))
            if stated not in roles:
                raise VocabularyError("%s.role: %s is not a role in family %s; its roles are %s"
                                      % (here, stated, family, ", ".join(roles)))
            derived = role_of(sub, doc, family)
            if derived != stated:
                raise VocabularyError(
                    "%s.role says %s, but what it %s makes it %s" % (
                        here, stated, "solves, enforces and introduces",
                        derived if derived else "no role of family %s" % family))


def check(doc: Dict[str, Any]) -> None:
    """Every vocabulary word the document states is ymf's (or an alias of
    one); a stated family or sub-problem role agrees with what the model and
    the sub-problem's structure say.
    Raises :class:`VocabularyError`."""
    physical = _physical(doc)
    family = None
    if physical.get("family") is not None:
        family = _word("family", physical["family"], "model.physical_model.family")
        found = families_of(doc)
        if family not in found:
            raise VocabularyError(
                "model.physical_model.family is %s, but its processes and assumptions %s"
                % (family, "identify %s" % ", ".join(found) if found else
                   "identify no family (%s needs %s)" % (family, _needs(family))))
    else:
        found = families_of(doc)
        family = found[0] if len(found) == 1 else None
    for entry in ((doc.get("solution_paths") or {}).get("discretizations") or []):
        if isinstance(entry, dict):
            _check_entry(doc, entry, family)


def _needs(family: str) -> str:
    needs = FAMILIES[family]
    return "; ".join("%s %s" % (k, ", ".join(needs[k])) for k in ("processes", "assumptions")
                     if needs.get(k))
