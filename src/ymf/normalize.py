"""Normalization helpers for YMF's backward-compatible schema evolution.

See ``docs/ymf-schema-v0.2-delta.md`` §2 for the rationale: ``unknowns``
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
