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

from typing import Any, Dict, List, Tuple

__all__ = ["STATES", "YmfStateError", "reached", "check"]

STATES = ("abstract", "well-posed", "discretized", "closed")


class YmfStateError(ValueError):
    """A document declares a state its content does not reach."""


def _discretizations(doc: Dict[str, Any]) -> List[Dict[str, Any]]:
    return list(((doc.get("solution_paths") or {}).get("discretizations")) or [])


def reached(doc: Dict[str, Any], well_posed: bool = True) -> Tuple[str, str]:
    """``(state, why not the next one)`` for a loaded document.

    ``well_posed`` says whether the document passed the full schema (the
    loader knows; a document that did not is ``abstract``). The reason is
    empty for ``closed``.
    """
    from ymf import closure
    if not well_posed:
        return "abstract", "it does not pass the full schema"
    discretizations = _discretizations(doc)
    if not discretizations:
        return "well-posed", "it has no discretizations"
    meshless = [d.get("name") for d in discretizations if not closure.realizations(d)]
    if meshless:
        return "well-posed", "discretization %s gives no mesh (cells or h)" % (
            ", ".join(repr(n) for n in meshless),)
    planned = [key for key, _, _ in closure.planned_outputs(doc)]
    kept, _ = closure.partition(doc, doc.get("outputs") or {})
    present = {closure._base_key(k) for k in kept}
    missing = [k for k in planned if k not in present]
    if missing:
        return "discretized", "%d of %d planned outputs %s missing: %s" % (
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
