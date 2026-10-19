"""How far a YMF document has got: its ``ymf_state``.

A document is the record of a progression towards a closure, and
``ymf_state`` says where along it the document stands::

    abstract      part of a model: it must be extended before it can be
                  solved (equations and unknowns, say, with no domain,
                  conditions or values); validated against the partial
                  schema only
    well-posed    the model block is complete: it passes the full schema
    discretized   well-posed, and every discretization gives a mesh, so
                  every output it plans has a key
    closed        discretized, and every planned output is present and
                  still keyed to the input that produced it

Each state includes the ones before it. A document may declare its state;
the loader checks that the content supports it (declaring less than the
content reaches is allowed). Tools that advance a document -- a runner
writing an archive -- write the state it has reached.

``ymf_state`` describes the file, like ``extends``: it is not inherited
through composition.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

__all__ = ["STATES", "YmfStateError", "Branch", "branches", "reached", "check"]

STATES = ("abstract", "well-posed", "discretized", "closed")


class YmfStateError(ValueError):
    """A document declares a state its content does not reach."""


def _discretizations(doc: Dict[str, Any]) -> List[Dict[str, Any]]:
    return list(((doc.get("solution_paths") or {}).get("discretizations")) or [])


@dataclass(frozen=True)
class Branch:
    """One discretization: the state it reaches, and why not the next one.

    ``reasons`` say what keeps the branch from the next state, each with a
    hint where one helps: "P2 is missing a mesh: give mesh: {h: ...}".
    ``missing`` lists the planned outputs not yet present (discretized).
    """
    name: str
    state: str
    reasons: Tuple[str, ...] = ()
    planned: Tuple[str, ...] = ()
    missing: Tuple[str, ...] = ()


def branches(doc: Dict[str, Any]) -> List[Branch]:
    """Each discretization of a well-posed document, as a :class:`Branch`."""
    from ymf import closure
    from ymf.discretizations import missing as lacks
    kept, _ = closure.partition(doc, doc.get("outputs") or {})
    present = {closure._base_key(k) for k in kept}
    out = []
    for entry in _discretizations(doc):
        name = entry.get("name")
        try:
            lacking = lacks(entry)
        except ValueError as exc:          # an unknown method
            out.append(Branch(name, "well-posed", (str(exc),)))
            continue
        if lacking:
            out.append(Branch(name, "well-posed", tuple(
                "%s is missing %s: %s" % (name, m.what, m.hint) for m in lacking)))
            continue
        planned = tuple(closure.output_key(doc, name, r)
                        for r in closure.realizations(entry))
        owed = tuple(k for k in planned if k not in present)
        if owed:
            out.append(Branch(name, "discretized", (), planned, owed))
        else:
            out.append(Branch(name, "closed", (), planned))
    return out


def reached(doc: Dict[str, Any], well_posed: bool = True) -> Tuple[str, str]:
    """``(state, why not the next one)`` for a loaded document.

    ``well_posed`` says whether the document passed the full schema (the
    loader knows; a document that did not is ``abstract``). Otherwise the
    document is as far as its least advanced branch (:func:`branches`), and
    the reason gathers what holds those branches back. The reason is empty
    for ``closed``.
    """
    if not well_posed:
        return "abstract", "it does not pass the full schema"
    every = branches(doc)
    if not every:
        return "well-posed", "it has no discretizations"
    state = min((b.state for b in every), key=STATES.index)
    if state == "well-posed":
        return state, "; ".join(r for b in every if b.state == state for r in b.reasons)
    if state == "discretized":
        planned = [k for b in every for k in b.planned]
        missing = [k for b in every for k in b.missing]
        return state, "%d of %d planned outputs %s missing: %s" % (
            len(missing), len(planned), "is" if len(missing) == 1 else "are",
            ", ".join(missing))
    return "closed", ""


def check(declared: str, doc: Dict[str, Any], label: str,
          well_posed: bool = True) -> None:
    """Refuse a declared state the content of ``doc`` does not reach."""
    if declared not in STATES:
        raise YmfStateError("%s: ymf_state %r is not one of %s"
                            % (label, declared, ", ".join(STATES)))
    state, why = reached(doc, well_posed)
    if STATES.index(declared) > STATES.index(state):
        raise YmfStateError(
            "%s: declares ymf_state: %s, but it is only %s: %s"
            % (label, declared, state, why))
