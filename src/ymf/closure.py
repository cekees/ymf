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

The key hashes the output's slice of the input, whole: the model without
its weak forms, the one weak form the discretization uses, the analytical
solutions it is checked against, the one discretization entry (with its
studied parameters replaced by the realization's values) and the
realization. Prose is included, so a key is a fingerprint of the input text
as well as of its numbers. Anything that changes the slice changes the key;
editing another discretization, or another weak form, does not.

Until 2026-10-18 (between releases: no release wrote archives) the slice
held the whole model block, every weak form included, so adding a weak form
re-keyed every output. An output keyed that way is recognised by
:func:`partition`, which re-keys it rather than drop it when its input is
otherwise unchanged (see :func:`legacy_input_slice`).

A discretization's ``mesh.cells`` or ``mesh.h``, its ``mesh.levels`` (and its
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
    "plain",
    "used_weak_form",
    "input_slice",
    "legacy_input_slice",
    "input_digest",
    "output_key",
    "planned_outputs",
    "partition",
    "file_stem",
    "archive_stem",
    "load",
    "write",
]

#: hex digits of the digest shown in a key (the full digest is in the output)
KEY_DIGITS = 12

# The studied parameters: (realization name, path in the discretization entry)
_STUDIED = (
    ("cells", ("mesh", "cells")),
    ("h", ("mesh", "h")),
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

    Each is a mapping such as ``{"cells": 8, "levels": 1}`` or
    ``{"h": 0.1, "levels": 1}`` (and ``"dt"`` when the entry steps in time).
    An entry without ``mesh`` asks for none: the mesh must then come from
    the command line, as an override. ``cells`` and ``h`` are alternatives;
    an entry giving both is refused.
    """
    mesh = _get(discretization, ("mesh",))
    if isinstance(mesh, dict) and "cells" in mesh and "h" in mesh:
        raise ValueError("discretization %r: mesh gives both cells and h; "
                         "cells is for a box, h for any other region"
                         % (discretization.get("name"),))
    axes = []
    for name, path in _STUDIED:
        value = _get(discretization, path)
        if value is None:
            if name in _DEFAULTS and _get(discretization, ("mesh",)) is not None:
                axes.append((name, [_DEFAULTS[name]]))
            continue
        axes.append((name, list(value) if isinstance(value, (list, tuple)) else [value]))
    if not any(name in ("cells", "h") for name, _ in axes):
        return []
    names = [name for name, _ in axes]
    return [dict(zip(names, combo)) for combo in itertools.product(*(v for _, v in axes))]


def realize(discretization: Dict[str, Any], realization: Dict[str, Any]) -> Dict[str, Any]:
    """The discretization entry with each studied list replaced by one value."""
    entry = copy.deepcopy(plain(discretization))
    for name, path in _STUDIED:
        if name in realization and _get(entry, path[:-1]) is not None:
            _set(entry, path, realization[name])
    return entry


def plain(value: Any) -> Any:
    """strictyaml's mappings and sequences, as plain dicts and lists."""
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    return value


def _discretization(doc: Dict[str, Any], name: str) -> Dict[str, Any]:
    for entry in doc["solution_paths"]["discretizations"]:
        if entry["name"] == name:
            return entry
    raise KeyError("no discretization named %r" % (name,))


def used_weak_form(doc: Dict[str, Any], discretization: Dict[str, Any]
                   ) -> Optional[Dict[str, Any]]:
    """The weak form a discretization entry solves, or None if it names none.

    Its ``from_weak_form``; the composer fills that in where the document has
    one weak form, and this does the same for a document built without it.
    """
    weak_forms = list((doc.get("model") or {}).get("weak_forms") or [])
    label = discretization.get("from_weak_form")
    if label is None:
        return weak_forms[0] if len(weak_forms) == 1 else None
    for weak_form in weak_forms:
        if str(weak_form.get("label")) == str(label):
            return weak_form
    raise KeyError("discretization %r names weak form %r, which the model does not have"
                   % (discretization.get("name"), label))


def input_slice(doc: Dict[str, Any], discretization: str,
                realization: Dict[str, Any]) -> Dict[str, Any]:
    """The part of the input one output depends on (see the module docstring)."""
    paths = doc.get("solution_paths") or {}
    entry = _discretization(doc, discretization)
    return plain({
        "model": {k: v for k, v in doc["model"].items() if k != "weak_forms"},
        "weak_form": used_weak_form(doc, entry),
        "analytical": paths.get("analytical") or [],
        "discretization": realize(entry, realization),
        "realization": realization,
    })


def legacy_input_slice(doc: Dict[str, Any], discretization: str,
                       realization: Dict[str, Any]) -> Dict[str, Any]:
    """The slice until 2026-10-18: the whole model block, every weak form in it.

    Kept only to recognise outputs keyed that way, so they can be re-keyed
    instead of recomputed.
    """
    paths = doc.get("solution_paths") or {}
    return plain({
        "model": doc["model"],
        "analytical": paths.get("analytical") or [],
        "discretization": realize(_discretization(doc, discretization), realization),
        "realization": realization,
    })


def _digest(slice_: Dict[str, Any]) -> str:
    text = json.dumps(slice_, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def input_digest(doc: Dict[str, Any], discretization: str,
                 realization: Dict[str, Any]) -> str:
    """sha256 of the input slice in a canonical JSON form."""
    return _digest(input_slice(doc, discretization, realization))


def _key(discretization: str, realization: Dict[str, Any], digest: str) -> str:
    size = "cells" if "cells" in realization else "h"
    parts = [discretization, "%s=%s" % (size, realization[size])]
    if realization.get("levels", 1) != 1:
        parts.append("levels=%s" % (realization["levels"],))
    if "dt" in realization:
        parts.append("dt=%s" % (realization["dt"],))
    parts.append(digest[:KEY_DIGITS])
    return "/".join(parts)


def output_key(doc: Dict[str, Any], discretization: str,
               realization: Dict[str, Any]) -> str:
    """``<discretization>/cells=N[/levels=L][/dt=T]/<digest prefix>``,
    with ``h=H`` in place of ``cells=N`` for an unstructured mesh."""
    return _key(discretization, realization, input_digest(doc, discretization, realization))


def planned_outputs(doc: Dict[str, Any]) -> Iterator[Tuple[str, str, Dict[str, Any]]]:
    """Every (key, discretization name, realization) the spec asks for."""
    for entry in doc["solution_paths"]["discretizations"]:
        for realization in realizations(entry):
            yield output_key(doc, entry["name"], realization), entry["name"], realization


def _base_key(key: str) -> str:
    """A key without its ``~N`` suffix (a kept, differing rerun)."""
    return key.split("~", 1)[0]


def load(path: str | Path, check_state: bool = True
         ) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, str]]:
    """Load a spec or archive: ``(spec, outputs, dropped)``.

    ``spec`` is the composed, validated input. ``outputs`` are those -- the
    file's own and any inherited through ``extends`` -- whose key, recomputed
    against the composed input, is still theirs. ``dropped`` maps every other
    key to why: an inherited output whose input the overrides changed belongs
    to the old archive, not to this one. ``check_state`` is
    :func:`ymf.compose.load_composed`'s: False for a runner that brings an
    over-claiming archive up to date.
    """
    from ymf.compose import load_composed
    spec = load_composed(path, check_state=check_state)
    outputs, dropped = partition(spec, spec.pop("outputs", None) or {})
    return spec, outputs, dropped


def partition(spec: Dict[str, Any], outputs: Dict[str, Any]
              ) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """Split ``outputs`` into those ``spec`` still produces, and the rest.

    Returns ``(kept, dropped)``; ``dropped`` maps each other key to why.

    An output keyed by the slice used until 2026-10-18 (:func:`legacy_input_slice`)
    whose input is otherwise unchanged is kept under its current key, with
    ``input.sha256`` brought up to date and ``input.rekeyed_from`` naming
    the key it had. Its files are not renamed: an output records its own
    (``storage``), and only new outputs are named by their key.
    """
    kept, dropped = {}, {}
    for key, output in outputs.items():
        base = _base_key(key)
        name = base.split("/", 1)[0]
        realization = (output or {}).get("realization")
        if not realization:
            dropped[key] = "it records no realization"
            continue
        try:
            current = output_key(spec, name, realization)
        except KeyError:
            dropped[key] = "the discretization %r is not in this input" % (name,)
            continue
        if current == base:
            kept[key] = output
        elif _key(name, realization, _digest(legacy_input_slice(spec, name, realization))) == base:
            output = copy.deepcopy(output)
            output["input"] = dict(output.get("input") or {},
                                   sha256=input_digest(spec, name, realization),
                                   rekeyed_from=key)
            kept[current + key[len(base):]] = output
        else:
            dropped[key] = "its input changed (it would now be %s)" % (current,)
    return kept, dropped


def file_stem(stem: str, key: str) -> str:
    """The name of an output's files: ``poisson`` + key -> ``poisson_3f2a9c01d4e7``.

    A kept, differing rerun (``key~2``) gets ``poisson_3f2a9c01d4e7_2``.
    """
    base, _, n = key.partition("~")
    return "%s_%s%s" % (stem, base.rsplit("/", 1)[-1], "_" + n if n else "")


def archive_stem(path: str | Path) -> str:
    """``poisson.ymf`` and ``poisson.archive.ymf`` -> ``poisson``."""
    name = Path(path).name
    for suffix in (".ymf", ".archive"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def write(path: str | Path, spec: Dict[str, Any], outputs: Dict[str, Any]) -> None:
    """Write an archive: ``ymf:``, its ``ymf_state``, the realized spec, then
    ``outputs``. The state is the one the archive reaches (see ymf.state),
    whatever its input declared."""
    from ymf.state import reached
    document = dict(plain(spec))
    document.pop("outputs", None)
    document.pop("ymf_state", None)
    document.pop("generated", None)     # it describes the input file, not this one
    if outputs:
        document["outputs"] = outputs
    state, _ = reached(document)
    write_document(path, dict({"ymf_state": state}, **document))
