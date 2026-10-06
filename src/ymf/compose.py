"""Compose a YMF document from several files with ``extends``.

A document may name another document it builds on::

    extends: navier_stokes_model.yaml     # relative to this file

The usual shape is a **model** file that states the physics and the
equations but is not yet a well-posed problem (no domain, boundary
conditions or coefficient values), marked ``kind: model``, and several
**problem** files that each extend it into a complete problem. Chains are
allowed: a problem can extend another problem to vary a parameter.

Validation
----------
Every file in a chain is validated on its own as it is read, so an error is
reported against the file and line that holds it:

- a ``kind: model`` file against :data:`ymf.schema.PARTIAL_YMF_SCHEMA`,
  in which any section or field may be missing but whatever is present
  must be well formed;
- every other file against the partial schema too (a problem file holds
  only what it adds), and the composed result against the full
  :data:`ymf.schema.YMF_SCHEMA`.

Merge rules
-----------
The child is merged onto the parent:

- mappings merge key by key;
- lists of mappings that carry an identity (``label``, else ``name``)
  merge entry by entry on it, so a child can add a weak form or fill in a
  field of an existing unknown; new entries are appended;
- lists of strings (``processes``, ``assumptions``, ...) are unioned,
  keeping order;
- anything else the child sets replaces the parent's value.

When a child changes a value its parent had set, that is an **override**,
and it is recorded rather than refused. The composed document gets a
``composition`` block naming every source file and every override::

    composition:
      sources: [navier_stokes_model.yaml, plane_poiseuille.yaml, ...]
      overrides:
        - path: Problem.characteristic_scales.viscosity.value
          was: 1.0
          now: 0.01
          set_by: plane_poiseuille_re100.yaml

``extends`` and ``kind`` describe a file, not the problem, and are not
merged: the composed document is self-contained. A composed *model* is
returned with ``kind: model``, so that whoever receives it can tell it is
not a problem yet.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import strictyaml
import yaml

from ymf.schema import _ALLOW_FLOW_STYLE, PARTIAL_YMF_SCHEMA, YMF_SCHEMA

__all__ = ["YmfCompositionError", "load_composed"]

#: Keys that describe a file rather than the problem; never merged.
_FILE_KEYS = ("extends", "kind", "composition")
#: Keys that identify an entry in a list of mappings, in order of preference.
_IDENTITY_KEYS = ("label", "name")


class YmfCompositionError(ValueError):
    """An ``extends`` chain cannot be followed (missing file, cycle)."""


def _identity(entry: Any) -> Optional[Tuple[str, Any]]:
    if isinstance(entry, dict):
        for key in _IDENTITY_KEYS:
            if key in entry:
                return key, entry[key]
    return None


def _merge(parent: Any, child: Any, path: str, source: str,
           overrides: List[Dict[str, Any]]) -> Any:
    """Merge ``child`` onto ``parent``, appending any overrides found."""
    if isinstance(parent, dict) and isinstance(child, dict):
        merged = dict(parent)
        for key, value in child.items():
            where = "%s.%s" % (path, key) if path else key
            merged[key] = (_merge(parent[key], value, where, source, overrides)
                           if key in parent else value)
        return merged

    if isinstance(parent, list) and isinstance(child, list):
        if all(isinstance(x, str) for x in parent + child):
            return parent + [x for x in child if x not in parent]
        if parent and child and all(_identity(x) for x in parent + child):
            merged = list(parent)
            index = {_identity(x)[1]: i for i, x in enumerate(parent)}
            for entry in child:
                key = _identity(entry)[1]
                if key in index:
                    i = index[key]
                    merged[i] = _merge(parent[i], entry, "%s[%s]" % (path, key),
                                       source, overrides)
                else:
                    merged.append(entry)
            return merged

    if parent != child:
        overrides.append({"path": path, "was": parent, "now": child, "set_by": source})
    return child


def _read(path: Path) -> Dict[str, Any]:
    """Validate one file against the partial schema, reporting its own name."""
    return strictyaml.dirty_load(
        path.read_text(encoding="utf-8"), PARTIAL_YMF_SCHEMA,
        label=str(path), allow_flow_style=_ALLOW_FLOW_STYLE,
    ).data


def _compose(path: Path, chain: Tuple[Path, ...]) -> Tuple[Dict[str, Any], List[Path], List[Dict[str, Any]]]:
    """Return ``(merged data, source files root-first, overrides)``."""
    if path in chain:
        cycle = " -> ".join(str(p) for p in chain + (path,))
        raise YmfCompositionError("extends cycle: %s" % (cycle,))
    if not path.is_file():
        raise YmfCompositionError(
            "%s extends %s, which does not exist" % (chain[-1], path) if chain
            else "no such file: %s" % (path,))
    data = _read(path)
    own = {k: v for k, v in data.items() if k not in _FILE_KEYS}
    if "extends" not in data:
        return own, [path], []
    parent, sources, overrides = _compose(
        (path.parent / data["extends"]).resolve(), chain + (path,))
    merged = _merge(parent, own, "", path.name, overrides)
    return merged, sources + [path], overrides


def load_composed(path: str | Path) -> Dict[str, Any]:
    """Load a YMF document, following ``extends``, and validate the result.

    A file marked ``kind: model`` is returned validated against the partial
    schema only, since a model is not a complete problem. Anything else must
    compose into a document that passes the full schema.
    """
    path = Path(path).resolve()
    if not path.is_file():
        raise YmfCompositionError("no such file: %s" % (path,))
    data = _read(path)
    kind = data.get("kind", "problem")
    if "extends" not in data and kind == "problem":
        # A self-contained problem: validate the text itself, so errors
        # cite its own lines exactly as they appear in the file.
        return strictyaml.dirty_load(
            path.read_text(encoding="utf-8"), YMF_SCHEMA,
            label=str(path), allow_flow_style=_ALLOW_FLOW_STYLE,
        ).data
    merged, sources, overrides = _compose(path, ())
    if len(sources) > 1:
        base = path.parent
        merged["composition"] = {
            "sources": [os.path.relpath(p, base) for p in sources]}
        # Omitted when empty: strictyaml cannot write an empty list under
        # an untyped (Any) key.
        if overrides:
            merged["composition"]["overrides"] = overrides
    if kind == "model":
        # Kept, so that whoever receives this knows it is not a problem yet.
        merged["kind"] = "model"
        return merged
    # The composed document has no file of its own, so it is written out as
    # YAML and validated as text, under a label naming the chain. (Not
    # strictyaml.as_document, which cannot write an empty list or mapping
    # even where the schema allows one.) Errors in any single file were
    # already reported against that file by _read.
    label = "%s (composed from %s)" % (path.name, ", ".join(p.name for p in sources))
    text = yaml.safe_dump(merged, allow_unicode=True, sort_keys=False,
                          default_flow_style=False, width=1000)
    data = strictyaml.dirty_load(text, YMF_SCHEMA, label=label,
                                 allow_flow_style=_ALLOW_FLOW_STYLE).data
    # The schema types "composition" as Any, which reads every scalar back
    # as a string. The loader wrote it, so return its own typed copy.
    data["composition"] = merged["composition"]
    return data
