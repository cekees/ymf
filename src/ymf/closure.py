"""The archive as the closure of a spec.

A spec is the input; an archive is that spec, as realized, followed by its
``outputs``. Each output is one discretization run on one mesh (one
*realization*), and is keyed by what produced it::

    ymf: 0.3.0
    Problem: ...
    solution_paths: ...
    outputs:
      P1/cells=8/3f2a9c01d4e7:
        input: {sha256: 3f2a9c01d4e7...}
        realization: {cells: 8, levels: 1}
        transformations: ...
        approximation: {...}          # the XDMF Domain

The key hashes the output's slice of the input, whole: the Problem, the
analytical solutions it is checked against, the one discretization entry
(with its studied parameters replaced by the realization's values) and the
realization. Prose is included, so a key is a fingerprint of the input text
as well as of its numbers. Anything that changes the slice changes the key;
editing another discretization does not.

A discretization's ``mesh.cells`` and ``mesh.levels`` (and its
``discretization.dt``) may be lists: a study. Each combination is one
realization, and one output.
"""

from __future__ import annotations

import copy
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

from ymf.archive import write_document

__all__ = [
    "KEY_DIGITS",
    "realizations",
    "realize",
    "input_slice",
    "input_digest",
    "output_key",
    "planned_outputs",
    "load",
    "write",
]

#: hex digits of the digest shown in a key (the full digest is in the output)
KEY_DIGITS = 12

# The studied parameters: (realization name, path in the discretization entry)
_STUDIED = (
    ("cells", ("mesh", "cells")),
    ("levels", ("mesh", "levels")),
    ("dt", ("discretization", "dt")),
)
_DEFAULTS = {"levels": 1}


def _get(entry: Dict[str, Any], path: Tuple[str, ...]) -> Any:
    for part in path:
        if not isinstance(entry, dict) or part not in entry:
            return None
        entry = entry[part]
    return entry


def _set(entry: Dict[str, Any], path: Tuple[str, ...], value: Any) -> None:
    for part in path[:-1]:
        entry = entry[part]
    entry[path[-1]] = value


def realizations(discretization: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The realizations a discretization entry asks for, in study order.

    Each is a mapping such as ``{"cells": 8, "levels": 1}`` (and ``"dt"``
    when the entry steps in time). An entry without ``mesh`` asks for none:
    the cells must then come from the command line, as an override.
    """
    axes = []
    for name, path in _STUDIED:
        value = _get(discretization, path)
        if value is None:
            if name in _DEFAULTS and _get(discretization, ("mesh",)) is not None:
                axes.append((name, [_DEFAULTS[name]]))
            continue
        axes.append((name, list(value) if isinstance(value, (list, tuple)) else [value]))
    if not any(name == "cells" for name, _ in axes):
        return []
    names = [name for name, _ in axes]
    return [dict(zip(names, combo)) for combo in itertools.product(*(v for _, v in axes))]


def realize(discretization: Dict[str, Any], realization: Dict[str, Any]) -> Dict[str, Any]:
    """The discretization entry with each studied list replaced by one value."""
    entry = copy.deepcopy(_plain(discretization))
    for name, path in _STUDIED:
        if name in realization and _get(entry, path[:-1]) is not None:
            _set(entry, path, realization[name])
    return entry


def _plain(value: Any) -> Any:
    """strictyaml's mappings and sequences, as plain dicts and lists."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value


def _discretization(doc: Dict[str, Any], name: str) -> Dict[str, Any]:
    for entry in doc["solution_paths"]["discretizations"]:
        if entry["name"] == name:
            return entry
    raise KeyError("no discretization named %r" % (name,))


def input_slice(doc: Dict[str, Any], discretization: str,
                realization: Dict[str, Any]) -> Dict[str, Any]:
    """The part of the input one output depends on (see the module docstring)."""
    paths = doc.get("solution_paths") or {}
    return _plain({
        "Problem": doc["Problem"],
        "analytical": paths.get("analytical") or [],
        "discretization": realize(_discretization(doc, discretization), realization),
        "realization": realization,
    })


def input_digest(doc: Dict[str, Any], discretization: str,
                 realization: Dict[str, Any]) -> str:
    """sha256 of the input slice in a canonical JSON form."""
    text = json.dumps(input_slice(doc, discretization, realization), sort_keys=True,
                      ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def output_key(doc: Dict[str, Any], discretization: str,
               realization: Dict[str, Any]) -> str:
    """``<discretization>/cells=N[/levels=L][/dt=T]/<digest prefix>``."""
    parts = [discretization, "cells=%s" % (realization["cells"],)]
    if realization.get("levels", 1) != 1:
        parts.append("levels=%s" % (realization["levels"],))
    if "dt" in realization:
        parts.append("dt=%s" % (realization["dt"],))
    parts.append(input_digest(doc, discretization, realization)[:KEY_DIGITS])
    return "/".join(parts)


def planned_outputs(doc: Dict[str, Any]) -> Iterator[Tuple[str, str, Dict[str, Any]]]:
    """Every (key, discretization name, realization) the spec asks for."""
    for entry in doc["solution_paths"]["discretizations"]:
        for realization in realizations(entry):
            yield output_key(doc, entry["name"], realization), entry["name"], realization


def _base_key(key: str) -> str:
    """A key without its ``~N`` suffix (a kept, differing rerun)."""
    return key.split("~", 1)[0]


def load(path: str | Path) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, str]]:
    """Load a spec or archive: ``(spec, outputs, dropped)``.

    ``spec`` is the composed, validated input. ``outputs`` are those -- the
    file's own and any inherited through ``extends`` -- whose key, recomputed
    against the composed input, is still theirs. ``dropped`` maps every other
    key to why: an inherited output whose input the overrides changed belongs
    to the old archive, not to this one.
    """
    from ymf.compose import load_composed
    doc = load_composed(path)
    candidates = doc.pop("outputs", None) or {}
    spec = doc
    outputs, dropped = {}, {}
    for key, output in candidates.items():
        base = _base_key(key)
        name = base.split("/", 1)[0]
        realization = (output or {}).get("realization")
        try:
            current = output_key(spec, name, realization) if realization else None
        except KeyError:
            dropped[key] = "the discretization %r is not in this input" % (name,)
            continue
        if current == base:
            outputs[key] = output
        else:
            dropped[key] = "its input changed (it would now be %s)" % (current,)
    return spec, outputs, dropped


def write(path: str | Path, spec: Dict[str, Any], outputs: Dict[str, Any]) -> None:
    """Write an archive: ``ymf:``, the realized spec, then ``outputs``."""
    document = dict(_plain(spec))
    document.pop("outputs", None)
    if outputs:
        document["outputs"] = outputs
    write_document(path, document)
