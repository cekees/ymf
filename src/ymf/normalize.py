"""Normalization helpers for YMF's backward-compatible schema evolution.

See ``docs/history/schema-v0.2-delta.md`` §2 for the rationale: ``unknowns``
widened from ``Seq(Str())`` (v0.1) to accept either bare strings or
structured maps with ``units``/``std_name``/``provenance`` (v0.2). Downstream code should
call :func:`normalize_unknowns` immediately after parsing so it only ever
sees one uniform shape.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

UnknownEntry = Union[str, Dict[str, Any]]


def normalize_unknowns(
    unknowns: List[UnknownEntry], provenance: Optional[str] = None
) -> List[Dict[str, Any]]:
    """Normalize v0.1 (bare str) or v0.2 (map) unknowns to a uniform list.

    Each returned entry has keys ``name``, ``units``, ``std_name`` and
    ``provenance`` (each ``None`` if absent), regardless of which input form
    was used. Pass the strong form's ``provenance`` as ``provenance`` to
    have it apply to every unknown that doesn't state its own.

    >>> normalize_unknowns(["u", "p"])
    [{'name': 'u', 'units': None, 'std_name': None, 'provenance': None}, {'name': 'p', 'units': None, 'std_name': None, 'provenance': None}]

    >>> normalize_unknowns([{"name": "u", "units": "m/s"},
    ...                     {"name": "p", "units": "Pa", "provenance": "human_edited"}],
    ...                    provenance="llm_derived")
    [{'name': 'u', 'units': 'm/s', 'std_name': None, 'provenance': 'llm_derived'}, {'name': 'p', 'units': 'Pa', 'std_name': None, 'provenance': 'human_edited'}]
    """
    result: List[Dict[str, Any]] = []
    for entry in unknowns:
        if isinstance(entry, str):
            result.append({"name": entry, "units": None, "std_name": None,
                           "provenance": provenance})
        elif isinstance(entry, dict):
            result.append(
                {
                    "name": entry["name"],
                    "units": entry.get("units"),
                    "std_name": entry.get("std_name"),
                    "provenance": entry.get("provenance", provenance),
                }
            )
        else:
            raise TypeError(
                f"Unsupported unknowns entry type: {type(entry)!r} (value={entry!r})"
            )
    return result


EquationEntry = Union[str, Dict[str, Any]]


def normalize_equations(equations: List[EquationEntry]) -> List[Dict[str, Any]]:
    """Normalize ``strong_form.equations`` to ``{name, for, equation}`` entries.

    An equation is written either as a bare string, ordered like the
    unknowns, or as a mapping that names it and, with ``for``, the unknown
    it determines::

        - name: momentum
          for: v
          equation: "ρ ∂v/∂t + ∇·(ρ v⊗v) − ∇·(μ∇v) + ∇p = ρ g  in Ω"

    A name is what a discretization's sub-problem refers to (``solves:
    [momentum]``), and what an extending file merges onto. ``name`` and
    ``for`` are ``None`` for a bare string.

    >>> normalize_equations(["-Δu = f  in Ω"])
    [{'name': None, 'for': None, 'equation': '-Δu = f  in Ω'}]
    """
    result: List[Dict[str, Any]] = []
    for entry in equations:
        if isinstance(entry, str):
            result.append({"name": None, "for": None, "equation": entry})
        elif isinstance(entry, dict):
            result.append({"name": entry.get("name"), "for": entry.get("for"),
                           "equation": entry.get("equation")})
        else:
            raise TypeError(
                f"Unsupported equations entry type: {type(entry)!r} (value={entry!r})"
            )
    return result


def equations_in_order(strong_form: Dict[str, Any], complete: bool = True) -> List[str]:
    """The equations' text, one per unknown, in the order of the unknowns.

    Bare strings are taken in the order written, as they always were. When
    every equation says which unknown it is ``for``, the order written does
    not matter: they are put in the unknowns' order. Raises ``ValueError``
    when some say and some do not, when one names something that is not an
    unknown, or when two are for the same unknown -- each is a document
    that does not say which equation determines which unknown.
    ``complete=False`` allows unknowns no equation is for yet (an abstract
    file the rest of the chain completes); the result then holds only the
    equations given.
    """
    entries = normalize_equations(strong_form.get("equations") or [])
    stated = [e for e in entries if e["for"] is not None]
    if not stated:
        return [e["equation"] for e in entries]
    if len(stated) != len(entries):
        raise ValueError(
            "equations: give 'for' on every equation or on none; %s %s but %s %s not"
            % (_listing([e["name"] for e in stated]), "says" if len(stated) == 1 else "say",
               _listing([e["name"] or e["equation"] for e in entries if e["for"] is None]),
               "does" if len(entries) - len(stated) == 1 else "do"))
    names = [u["name"] for u in normalize_unknowns(strong_form.get("unknowns") or [])]
    by_unknown: Dict[str, Dict[str, Any]] = {}
    for e in entries:
        if e["for"] not in names:
            raise ValueError("equations[%s] is for %r, which is not an unknown (have: %s)"
                             % (e["name"], e["for"], ", ".join(names)))
        if e["for"] in by_unknown:
            raise ValueError("equations %s and %s are both for %s"
                             % (by_unknown[e["for"]]["name"], e["name"], e["for"]))
        by_unknown[e["for"]] = e
    missing = [n for n in names if n not in by_unknown]
    if missing and complete:
        raise ValueError("equations: none is for %s" % ", ".join(missing))
    return [by_unknown[n]["equation"] for n in names if n in by_unknown]


def _listing(items: List[Any]) -> str:
    return ", ".join(str(i) for i in items)
