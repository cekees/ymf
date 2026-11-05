"""What a closed spec records: the input, with its shorthand written out.

A spec may leave things to rules: a boundary condition may be an alias, a
sub-problem's solver refines the entry's, its regularization replaces the
entry's function by function, a lag may be a bare ``previous_step``, a
coefficient read through the sub-problem that carries it is that carrier.
Each rule has one meaning, but if it lives only in a reader's code, two
readers can close the same spec differently and the archive cannot show
which. So the closed spec writes each one out (cekees, 2026-11-04, for the
boundary conditions; C42 round 2 for the rest):

* boundary-condition aliases and condition sets, expanded
  (:func:`ymf.conditions.expand`);
* every **coefficient of the state** -- one whose formula names an unknown,
  directly or through another coefficient (``ρ = ρw (1 − χ) + ρa χ`` with
  ``χ = heaviside(φ)``) -- in the map form, with ``depends_on``: the
  unknowns it follows. A solver of one system substitutes such a
  coefficient; a kernel generator emits it by name, as one inline function
  of the state, and this is where the closed spec says which names those
  are;
* every **sub-problem** of a split discretization with what it solves for
  (``for``, by default the unknown it introduces), its solver (the entry's,
  refined by its own), its regularization (the entry's, function by
  function replaced by its own), and each lag as ``{at, via, of,
  carried_by}``, where ``carried_by`` names the introduced unknown that
  carries a lagged coefficient (``χ`` via VOF is read as ``χh``).

Applying it twice changes nothing, and a spec with none of the shorthand
(no aliases, no coefficients of the state, no sub-problems) is returned
unchanged. Readers may read either form; what they solve is the same.
"""

from __future__ import annotations

import copy
from typing import Any, Dict, List

__all__ = ["explicit", "without_derived", "state_coefficients", "split_out", "carriers",
           "lag_of", "DERIVED_KEYS"]

#: keys explicit() derives from what the input already says, which therefore
#: do not count towards an output's key (ymf.closure.input_slice)
DERIVED_KEYS = ("depends_on",)


def _strong(doc: Dict[str, Any]) -> Dict[str, Any]:
    return ((doc.get("model") or {}).get("strong_form") or {})


def _unknowns(strong: Dict[str, Any]) -> List[str]:
    return [u if isinstance(u, str) else u.get("name")
            for u in strong.get("unknowns") or []]


def state_coefficients(doc: Dict[str, Any]) -> Dict[str, List[str]]:
    """``{coefficient: [unknown, ...]}``: the coefficients of the state, read
    from the formulas' text (ymf.names; no sympy)."""
    from ymf.names import state_dependence
    strong = _strong(doc)
    return state_dependence(strong.get("coefficients") or {}, _unknowns(strong))


def _name_state(doc: Dict[str, Any]) -> Dict[str, Any]:
    state = state_coefficients(doc)
    if not state:
        return doc
    doc = dict(doc)
    model = dict(doc["model"])
    strong = dict(model["strong_form"])
    coefficients = {}
    for name, entry in strong["coefficients"].items():
        if name in state:
            entry = dict(entry) if isinstance(entry, dict) else {"formula": str(entry)}
            entry["depends_on"] = list(state[name])
        coefficients[name] = entry
    strong["coefficients"] = coefficients
    model["strong_form"] = strong
    doc["model"] = model
    return doc


def carriers(entry: Dict[str, Any]) -> Dict[str, str]:
    """``{model quantity: introduced unknown}`` for one split entry: what each
    introduced carrier approximates (``χ`` -> ``χh``)."""
    out = {}
    for sub in entry.get("subproblems") or []:
        new = sub.get("introduces") or {}
        if new.get("approximates"):
            out[new["approximates"]] = new["name"]
    return out


def lag_of(name: str, lag: Any, carried: Dict[str, str]) -> Dict[str, Any]:
    """One lag written out: ``{of?, via?, at, carried_by?}``."""
    lag = dict(lag) if isinstance(lag, dict) else {"at": lag}
    out = {key: lag[key] for key in ("of", "via") if lag.get(key) is not None}
    out["at"] = lag["at"]
    target = lag.get("of", name)
    if target in carried:
        out["carried_by"] = carried[target]
    elif lag.get("carried_by") is not None:
        out["carried_by"] = lag["carried_by"]     # refused on load if it is wrong
    return out


def split_out(entry: Dict[str, Any]) -> Dict[str, Any]:
    """One discretization entry with every sub-problem written out (see the
    module docstring); an entry without sub-problems is returned as it is."""
    if not entry.get("subproblems"):
        return entry
    from ymf.discretizations import subproblem_solver
    entry = dict(entry)
    carried = carriers(entry)
    out = []
    for sub in entry["subproblems"]:
        sub = dict(sub)
        if sub.get("equation") and not sub.get("for") and sub.get("introduces"):
            sub["for"] = sub["introduces"]["name"]
        if sub.get("lags"):
            sub["lags"] = {name: lag_of(name, lag, carried)
                           for name, lag in sub["lags"].items()}
        regularize = dict(entry.get("regularize") or {})
        regularize.update(sub.get("regularize") or {})
        if regularize:
            sub["regularize"] = regularize
        solver = subproblem_solver(entry, sub)
        if solver:
            sub["solver"] = solver
        out.append(sub)
    entry["subproblems"] = out
    return entry


def explicit(doc: Dict[str, Any]) -> Dict[str, Any]:
    """``doc`` as a closed spec records it (see the module docstring)."""
    from ymf.conditions import expand
    doc = _name_state(expand(doc))
    entries = (doc.get("solution_paths") or {}).get("discretizations") or []
    if any(isinstance(e, dict) and e.get("subproblems") for e in entries):
        doc = dict(doc)
        paths = dict(doc["solution_paths"])
        paths["discretizations"] = [split_out(e) if isinstance(e, dict) else e
                                    for e in entries]
        doc["solution_paths"] = paths
    return doc


def without_derived(doc: Dict[str, Any]) -> Dict[str, Any]:
    """``doc`` with the keys :func:`explicit` derives from the rest removed."""
    strong = _strong(doc)
    coefficients = strong.get("coefficients") or {}
    if not any(isinstance(e, dict) and any(k in e for k in DERIVED_KEYS)
               for e in coefficients.values()):
        return doc
    doc = copy.deepcopy(doc)
    for entry in doc["model"]["strong_form"]["coefficients"].values():
        if isinstance(entry, dict):
            for key in DERIVED_KEYS:
                entry.pop(key, None)
    return doc
