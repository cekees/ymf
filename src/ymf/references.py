"""The reference graph: the edges between a document's branches.

A YMF document holds competing approaches to one problem and records how
they relate. Those relations are cross-references, and there are four:

    solution_paths.discretizations[].from_weak_form  -> weak_forms[].label
    vvuq.verification.reference.analytical           -> analytical[].name
    vvuq.verification.branches[]                     -> discretizations[].name
    vvuq.verification.expected_rates (keys)          -> discretizations[].name

Until now none of them was checked: a reference that named nothing still
validated, so a typo in a verification plan was simply ignored. They are
also where a document says what it *is* -- which weak form a
discretization implements is the difference between two branches of the
same problem and two different problems -- so leaving them unchecked made
the structure advisory.

Two of the four are single-valued, and those may be omitted when the graph
does not branch: a discretization in a document with one weak form can only
implement that one, and a verification plan in a document with one
analytical solution can only compare against that one. The omission is
filled in during composition, at the point the file is merged, and written
into the composed document, so the result has no implicit bindings and an
archive carries none either.

Filling in at that point rather than at the end is deliberate. A model with
one weak form is unambiguous when written; a problem extending it may add a
second. Resolving against the finished document would make the model's
omission retroactively ambiguous -- punishing a file for what a later file
did.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

__all__ = [
    "YmfReferenceError",
    "apply_defaults",
    "resolve",
    "weak_form_labels",
    "analytical_names",
    "discretization_names",
]


class YmfReferenceError(ValueError):
    """A reference names a branch the document does not have."""


def _problem(doc: Dict[str, Any]) -> Dict[str, Any]:
    return doc.get("Problem") or {}


def _paths(doc: Dict[str, Any]) -> Dict[str, Any]:
    return doc.get("solution_paths") or {}


def _verification(doc: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    block = (doc.get("vvuq") or {}).get("verification")
    return block if isinstance(block, dict) else None


def weak_form_labels(doc: Dict[str, Any]) -> List[str]:
    return [w["label"] for w in _problem(doc).get("weak_forms") or []
            if isinstance(w, dict) and "label" in w]


def analytical_names(doc: Dict[str, Any]) -> List[str]:
    return [a["name"] for a in _paths(doc).get("analytical") or []
            if isinstance(a, dict) and "name" in a]


def discretization_names(doc: Dict[str, Any]) -> List[str]:
    return [d["name"] for d in _paths(doc).get("discretizations") or []
            if isinstance(d, dict) and "name" in d]


def _only(candidates: List[str]) -> Optional[str]:
    """The single candidate, or None if there are none or several."""
    return candidates[0] if len(candidates) == 1 else None


# ---------------------------------------------------------------------------
# defaulting
# ---------------------------------------------------------------------------


def apply_defaults(doc: Dict[str, Any], source: str) -> List[Dict[str, Any]]:
    """Fill in omitted references where the graph does not branch.

    Returns one record per reference filled in, for ``composition.defaults``.
    Nothing is filled in where there is no candidate (a document may have
    discretizations and no weak forms) or more than one -- there, omitting
    the reference is a question only the author can answer, and
    :func:`resolve` asks it.
    """
    filled: List[Dict[str, Any]] = []

    label = _only(weak_form_labels(doc))
    if label is not None:
        for entry in _paths(doc).get("discretizations") or []:
            if isinstance(entry, dict) and "from_weak_form" not in entry:
                entry["from_weak_form"] = label
                filled.append({
                    "path": "solution_paths.discretizations[%s].from_weak_form"
                            % entry.get("name", "?"),
                    "value": label,
                    "because": "the document has one weak form",
                    "defaulted_by": source,
                })

    verification = _verification(doc)
    if verification is not None:
        reference = verification.get("reference")
        if isinstance(reference, dict) and "analytical" not in reference:
            name = _only(analytical_names(doc))
            if name is not None:
                reference["analytical"] = name
                filled.append({
                    "path": "vvuq.verification.reference.analytical",
                    "value": name,
                    "because": "the document has one analytical solution",
                    "defaulted_by": source,
                })
    return filled


# ---------------------------------------------------------------------------
# resolution
# ---------------------------------------------------------------------------


def _check(named: Any, candidates: List[str], path: str, what: str) -> None:
    if not isinstance(named, str):
        raise YmfReferenceError(
            "%s must name a %s; found %r" % (path, what, named))
    if named not in candidates:
        raise YmfReferenceError(
            "%s names %r, which is not a %s in this document (have: %s)"
            % (path, named, what, ", ".join(candidates) or "none"))


def resolve(doc: Dict[str, Any]) -> None:
    """Check every reference, raising :class:`YmfReferenceError` on the first bad one."""
    labels = weak_form_labels(doc)
    discretizations = discretization_names(doc)

    for entry in _paths(doc).get("discretizations") or []:
        if not isinstance(entry, dict):
            continue
        named = entry.get("from_weak_form")
        where = "solution_paths.discretizations[%s].from_weak_form" % (
            entry.get("name", "?"),)
        if named is None:
            if labels:
                raise YmfReferenceError(
                    "%s is omitted, but this document has %d weak forms (%s), "
                    "so which one it implements has to be said"
                    % (where, len(labels), ", ".join(labels)))
            continue              # no weak forms at all: nothing to name
        _check(named, labels, where, "weak form label")

    verification = _verification(doc)
    if verification is None:
        return

    if "branch_ref" in (verification.get("reference") or {}):
        raise YmfReferenceError(
            "vvuq.verification.reference.branch_ref was renamed to "
            "vvuq.verification.reference.analytical")

    reference = verification.get("reference")
    if isinstance(reference, dict) and "analytical" in reference:
        _check(reference["analytical"], analytical_names(doc),
               "vvuq.verification.reference.analytical", "analytical solution")

    for i, named in enumerate(verification.get("branches") or []):
        _check(named, discretizations,
               "vvuq.verification.branches[%d]" % i, "discretization")

    for named in (verification.get("expected_rates") or {}):
        _check(named, discretizations,
               "vvuq.verification.expected_rates.%s" % named, "discretization")
