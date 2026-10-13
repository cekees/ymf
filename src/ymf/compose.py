"""Compose a YMF document from several files with ``extends``.

A document may name another document it builds on::

    extends: navier_stokes_model.ymf     # relative to this file

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
      sources: [navier_stokes_model.ymf, plane_poiseuille.ymf, ...]
      overrides:
        - path: Problem.characteristic_scales.viscosity.value
          was: 1.0
          now: 0.01
          set_by: plane_poiseuille_re100.ymf

``extends`` and ``kind`` describe a file, not the problem, and are not
merged: the composed document is self-contained. A composed *model* is
returned with ``kind: model``, so that whoever receives it can tell it is
not a problem yet.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import strictyaml
import yaml

from ymf.domain import validate_domain_block
from ymf.references import apply_defaults, resolve
from ymf.schema import _ALLOW_FLOW_STYLE, PARTIAL_YMF_SCHEMA, YMF_SCHEMA

__all__ = ["YmfCompositionError", "load_composed"]

#: Keys that describe a file rather than the problem; never merged.
_FILE_KEYS = ("extends", "kind", "composition", "remove")
#: Keys that identify an entry in a list of mappings, in order of preference.
_IDENTITY_KEYS = ("label", "name")


class YmfCompositionError(ValueError):
    """An ``extends`` chain cannot be followed (missing file, cycle)."""


class YmfRemovalError(YmfCompositionError):
    """A ``remove:`` path names nothing in the document it composes."""


def _parse_path(path: str) -> List[Tuple[str, str]]:
    """``a.b[c].d`` -> ``[("key","a"), ("key","b"), ("item","c"), ("key","d")]``.

    The same shape ``_merge`` writes into the override log, so a path that
    appears in ``composition.overrides`` can be pasted into ``remove:``.
    A key containing a literal ``.`` or ``[`` cannot be addressed; no key in
    the schema does.
    """
    steps: List[Tuple[str, str]] = []
    for part in path.split("."):
        if not part:
            raise YmfRemovalError("empty segment in remove path %r" % (path,))
        while "[" in part:
            head, _, rest = part.partition("[")
            item, bracket, part = rest.partition("]")
            if not bracket:
                raise YmfRemovalError("unclosed '[' in remove path %r" % (path,))
            if head:
                steps.append(("key", head))
            steps.append(("item", item))
        if part:
            steps.append(("key", part))
    if not steps:
        raise YmfRemovalError("empty remove path")
    return steps


def _remove_at(data: Any, path: str) -> Any:
    """Delete what ``path`` names, returning it. Raise if it names nothing.

    A removal that matches nothing is an error rather than a no-op: it means
    the document being composed has moved under the file doing the removing,
    and silently composing the unremoved thing is how a variant ends up
    stating the opposite of what its author wrote.
    """
    steps = _parse_path(path)
    node = data
    for i, (kind, token) in enumerate(steps[:-1]):
        where = ".".join(t for _, t in steps[: i + 1])
        if kind == "key":
            if not isinstance(node, dict) or token not in node:
                raise YmfRemovalError(
                    "remove %s: nothing at %s" % (path, where))
            node = node[token]
        else:
            node = _find_item(node, token, path, where)
    kind, token = steps[-1]
    if kind == "key":
        if not isinstance(node, dict) or token not in node:
            raise YmfRemovalError("remove %s: nothing to remove" % (path,))
        return node.pop(token)
    index = _index_of(node, token, path)
    return node.pop(index)


def _find_item(node: Any, token: str, path: str, where: str) -> Any:
    index = _index_of(node, token, path, where)
    return node[index]


def _index_of(node: Any, token: str, path: str, where: str = "") -> int:
    """The position of ``token`` in a list, by value or by its identity key."""
    if not isinstance(node, list):
        raise YmfRemovalError(
            "remove %s: %s is not a list" % (path, where or path))
    for i, entry in enumerate(node):
        if entry == token:
            return i                      # a plain-string list
        identity = _identity(entry)
        if identity and identity[1] == token:
            return i                      # a list of labelled mappings
    raise YmfRemovalError("remove %s: no entry %r" % (path, token))


def _apply_removals(data: Dict[str, Any], paths: List[str],
                    source: str) -> List[Dict[str, Any]]:
    """Apply a file's ``remove:`` paths, returning the record of what went.

    Removals run after the file has been merged onto what it extends, so a
    file can drop an inherited item and add its own replacement in one step:
    removing the ``incompressible`` assumption while adding ``compressible``
    leaves only the second, because the string lists union first and the
    removal then takes the one named.
    """
    record = []
    for path in paths:
        was = _remove_at(data, path)
        record.append({"path": path, "was": was, "removed_by": source})
    return record


def _identity(entry: Any) -> Optional[Tuple[str, Any]]:
    if isinstance(entry, dict):
        for key in _IDENTITY_KEYS:
            if key in entry:
                return key, entry[key]
    return None


def _note_writer(provenance: Dict[str, List[str]], path: str,
                 source: str) -> None:
    """Record ``source`` as a writer of ``path``, keeping chain order."""
    writers = provenance.setdefault(path, [])
    if source not in writers:
        writers.append(source)


def _record_paths(data: Any, path: str, source: str,
                  provenance: Dict[str, List[str]]) -> None:
    """Note ``source`` as a writer of every mapping key under ``data``."""
    if isinstance(data, dict):
        for key, value in data.items():
            where = "%s.%s" % (path, key) if path else key
            _note_writer(provenance, where, source)
            _record_paths(value, where, source, provenance)
    elif isinstance(data, list):
        for entry in data:
            identity = _identity(entry)
            if identity:
                where = "%s[%s]" % (path, identity[1])
                _note_writer(provenance, where, source)
                _record_paths(entry, where, source, provenance)


def _merge(parent: Any, child: Any, path: str, source: str,
           overrides: List[Dict[str, Any]],
           provenance: Optional[Dict[str, List[str]]] = None) -> Any:
    """Merge ``child`` onto ``parent``, appending any overrides found.

    ``provenance`` collects the file that last wrote each path, which is
    what lets a failure in the composed document name a file rather than a
    line of a document that exists nowhere on disk.
    """
    if provenance is None:
        provenance = {}
    if isinstance(parent, dict) and isinstance(child, dict):
        merged = dict(parent)
        for key, value in child.items():
            where = "%s.%s" % (path, key) if path else key
            _note_writer(provenance, where, source)
            if key in parent:
                merged[key] = _merge(parent[key], value, where, source,
                                     overrides, provenance)
            else:
                merged[key] = value
                _record_paths(value, where, source, provenance)
        return merged

    if isinstance(parent, list) and isinstance(child, list):
        if all(isinstance(x, str) for x in parent + child):
            return parent + [x for x in child if x not in parent]
        if parent and child and all(_identity(x) for x in parent + child):
            merged = list(parent)
            index = {_identity(x)[1]: i for i, x in enumerate(parent)}
            for entry in child:
                key = _identity(entry)[1]
                where = "%s[%s]" % (path, key)
                _note_writer(provenance, where, source)
                if key in index:
                    i = index[key]
                    merged[i] = _merge(parent[i], entry, where, source,
                                       overrides, provenance)
                else:
                    merged.append(entry)
                    _record_paths(entry, where, source, provenance)
            return merged

    if parent != child:
        overrides.append({"path": path, "was": parent, "now": child, "set_by": source})
    # Replaced wholesale -- a scalar, or a list the merge could not pair up
    # entry by entry (an empty parent, say). Whatever is inside it is this
    # file's now, so attribute it.
    _record_paths(child, path, source, provenance)
    return child


_KEY_LINE = re.compile(r"^(\s*)(- )?([^\s#][^:]*?):(?:\s+(.*))?$")


def _path_at_line(text: str, line: int) -> Optional[str]:
    """The dotted path of the mapping key at ``line`` of a dumped document.

    The document a composed problem is validated against is one this module
    dumped, so a validation error's line number can be turned back into a
    path, and from there into the file that wrote it.

    A list entry is named the way ``_merge`` and ``remove:`` name it --
    ``discretizations[new_path]`` -- rather than by the key that happens to
    carry its identity, so the path points at the entry and matches what
    the override and removal records print.
    """
    stack: List[Tuple[int, str, bool]] = []     # indent, token, is_item
    lines = text.split("\n")
    if not 0 <= line < len(lines):
        return None
    for current in lines[: line + 1]:
        match = _KEY_LINE.match(current)
        if not match:
            continue
        pad, dash, key, value = match.groups()
        indent = len(pad) + (2 if dash else 0)
        while stack and stack[-1][0] >= indent:
            stack.pop()
        key = key.strip()
        if dash and key in _IDENTITY_KEYS and value:
            # "- name: new_path" identifies the entry, not a field in it
            stack.append((indent, value.strip().strip("'\""), True))
        else:
            stack.append((indent, key, False))
    if not stack:
        return None
    out = ""
    for _, token, is_item in stack:
        out += ("[%s]" % token) if is_item else (("." + token) if out else token)
    return out or None


def _attribute(exc, text: str, provenance: Dict[str, List[str]],
               sources: List[Path]) -> None:
    """Add to ``exc`` the file responsible for the place it failed.

    A composed problem has no file of its own, so strictyaml can only point
    at a line of a document that exists nowhere on disk. Here that line is
    turned back into a path, and the path into the file that last wrote it.

    Which of the error's two marks to read depends on what went wrong. A
    missing required key is about the *mapping* that lacks it, which is the
    context mark; the problem mark points at the end of that mapping, which
    is some unrelated line. Anything else -- a value of the wrong type, a
    key the schema does not have -- is about the problem mark itself.

    The case the gaps doc named is a required key that no file in the chain
    sets. There the useful information is precisely that the chain is
    silent, so that is said outright, along with the file that writes the
    mapping it should have gone in.
    """
    problem = getattr(exc, "problem", "") or ""
    missing = re.findall(r"required key\(s\) (.+?) not found", problem)
    try:
        mark = exc.context_mark if missing else exc.problem_mark
        line = mark.line
    except Exception:                       # pragma: no cover - defensive
        return
    path = _path_at_line(text, line)
    if path is None:
        return

    if missing:
        keys = [k.strip().strip("'\"") for k in missing[0].split(",")]
        holders = provenance.get(path) or []
        where = ("%s is written by %s" % (path, ", ".join(holders)) if holders
                 else "%s comes from the chain %s"
                 % (path, " <- ".join(p.name for p in sources)))
        note = ("no file in the chain sets %s under %s; %s"
                % (", ".join(sorted(keys)), path, where))
    else:
        writers = provenance.get(path) or []
        note = ("at %s, last written by %s" % (path, writers[-1]) if writers
                else "at %s, which no file in the chain sets" % (path,))

    exc.note = (exc.note + "\n" + note) if getattr(exc, "note", None) else note


_ARCHIVE_KEYS = ("ymf", "outputs")


class _NoAliasDumper(yaml.SafeDumper):
    # strictyaml refuses anchors, and a merged document can share objects
    def ignore_aliases(self, data):
        return True


def _dump(data: Any, **kwargs: Any) -> str:
    return yaml.dump(data, Dumper=_NoAliasDumper, **kwargs)


def _is_archive_text(text: str) -> bool:
    import re
    return re.search(r"^(ymf|outputs):", text, re.M) is not None


#: Spellings accepted on input and folded to the canonical one as the file is
#: read. Keys are lower_case_with_underscores (see "Naming" in
#: docs/problem-specification.md); these two were capitalized before that rule
#: was written down.
#:
#: Folding happens before anything is merged, deliberately. Composition pairs
#: keys by exact match, so a model writing "Problem:" and a variant writing
#: "problem:" would not merge -- they would become siblings, and the result
#: would be a document with two problems. An alias eases a rename; it is not
#: a standing second spelling.
_KEY_ALIASES = {"Problem": "model", "problem": "model", "Xdmf": "xdmf"}


def _read_text(path: Path) -> str:
    """The file's text, with deprecated top-level spellings folded."""
    text = path.read_text(encoding="utf-8")
    for old, new in _KEY_ALIASES.items():
        # top level only: a nested key is indented, and a path inside a
        # string value ("Problem.strong_form...") never starts a line
        text = re.sub(r"^%s:" % old, "%s:" % new, text, flags=re.M)
    return text


def _read(path: Path) -> Tuple[Dict[str, Any], Optional[Dict[str, Any]]]:
    """Validate one file against the partial schema; return (data, outputs).

    A plain spec is validated as written, so errors cite its own lines. An
    archive (a spec with ``ymf:`` and ``outputs:``) is read with pyyaml --
    its outputs hold arrays, far too much for strictyaml -- and only its
    specification part is validated. Its outputs come back separately.
    """
    text = _read_text(path)
    if not _is_archive_text(text):
        return strictyaml.dirty_load(text, PARTIAL_YMF_SCHEMA, label=str(path),
                                     allow_flow_style=_ALLOW_FLOW_STYLE).data, None
    from ymf.archive import read_document
    document = read_document(path)
    outputs = document.pop("outputs", None)
    document.pop("ymf", None)
    spec_text = _dump(document, allow_unicode=True, sort_keys=False,
                               default_flow_style=False, width=1000)
    data = strictyaml.dirty_load(spec_text, PARTIAL_YMF_SCHEMA, label=str(path),
                                 allow_flow_style=_ALLOW_FLOW_STYLE).data
    if "composition" in document:
        data["composition"] = document["composition"]   # typed, as written
    return data, outputs


def _compose(path: Path, chain: Tuple[Path, ...]):
    """Return ``(merged, sources, overrides, removals, provenance, outputs, record)``."""
    if path in chain:
        cycle = " -> ".join(str(p) for p in chain + (path,))
        raise YmfCompositionError("extends cycle: %s" % (cycle,))
    if not path.is_file():
        raise YmfCompositionError(
            "%s extends %s, which does not exist" % (chain[-1], path) if chain
            else "no such file: %s" % (path,))
    data, outputs = _read(path)
    if outputs and chain:
        # an inherited output's files are relative to its own archive
        from ymf.archive import rebase_domain
        outputs = {k: dict(v, approximation=rebase_domain(
                       v["approximation"], path.parent, chain[0].parent))
                   if v.get("approximation") else v for k, v in outputs.items()}
    own = {k: v for k, v in data.items() if k not in _FILE_KEYS}
    drop = list(data.get("remove") or ())
    if "extends" not in data:
        provenance: Dict[str, List[str]] = {}
        _record_paths(own, "", path.name, provenance)
        removals = _apply_removals(own, drop, path.name) if drop else []
        defaults = apply_defaults(own, path.name)
        return (own, [path], [], removals, defaults, provenance,
                dict(outputs or {}), data.get("composition"))
    (parent, sources, overrides, removals, defaults, provenance, inherited,
     record) = _compose((path.parent / data["extends"]).resolve(),
                        chain + (path,))
    merged = _merge(parent, own, "", path.name, overrides, provenance)
    # after the merge, so a file can drop an inherited item and add its
    # replacement in the same step
    removals = removals + _apply_removals(merged, drop, path.name)
    # and after the removals: a file that drops a weak form changes what is
    # left to default to. "At that point" is this point, for this file.
    defaults = defaults + apply_defaults(merged, path.name)
    inherited.update(outputs or {})          # a file's own outputs win
    return (merged, sources + [path], overrides, removals, defaults,
            provenance, inherited, record)


def _check_document(doc: Dict[str, Any]) -> None:
    """The checks a finished document gets, whichever path produced it."""
    resolve(doc)
    strong = (doc.get("model") or {}).get("strong_form") or {}
    if "domain" in strong:
        validate_domain_block(strong["domain"], "model.strong_form.domain")


def load_composed(path: str | Path) -> Dict[str, Any]:
    """Load a YMF document, following ``extends``, and validate the result.

    A file marked ``kind: model`` is returned validated against the partial
    schema only, since a model is not a complete problem. Anything else must
    compose into a document that passes the full schema. An archive's
    ``outputs`` (its own, and any it inherits through ``extends``) are
    attached to the result unvalidated; ymf.closure decides which of them
    still belong to the composed input.
    """
    path = Path(path).resolve()
    if not path.is_file():
        raise YmfCompositionError("no such file: %s" % (path,))
    text = _read_text(path)
    archive = _is_archive_text(text)
    data, _ = _read(path)
    kind = data.get("kind", "problem")
    if "extends" not in data and kind == "problem" and not archive:
        # A self-contained problem: validate the text itself, so errors
        # cite its own lines exactly as they appear in the file.
        document = strictyaml.dirty_load(
            text, YMF_SCHEMA, label=str(path), allow_flow_style=_ALLOW_FLOW_STYLE,
        ).data
        filled = apply_defaults(document, path.name)
        if filled:
            # recorded even for one file, so the document still says how it
            # came to hold a reference its author did not write
            document["composition"] = {"defaults": filled}
        _check_document(document)
        return document
    (merged, sources, overrides, removals, defaults, provenance, outputs,
     record) = _compose(path, ())
    composition = None
    if len(sources) > 1:
        base = path.parent
        composition = {"sources": [os.path.relpath(p, base) for p in sources]}
        # Omitted when empty: strictyaml cannot write an empty list under
        # an untyped (Any) key.
        if overrides:
            composition["overrides"] = overrides
        if record:
            composition["inherited"] = record
    elif record:
        composition = record                  # an archive's own record, kept
    # A single file can remove from itself and can have a reference filled
    # in for it, so these are recorded whether or not anything was composed.
    if removals or defaults:
        composition = dict(composition or {})
        if removals:
            composition["removals"] = removals
        if defaults:
            composition["defaults"] = defaults
    if composition:
        merged["composition"] = composition
    if kind == "model":
        # Kept, so that whoever receives this knows it is not a problem yet.
        merged["kind"] = "model"
        _check_document(merged)
        return merged
    # The composed document has no file of its own, so it is written out as
    # YAML and validated as text, under a label naming the chain. (Not
    # strictyaml.as_document, which cannot write an empty list or mapping
    # even where the schema allows one.) Errors in any single file were
    # already reported against that file by _read.
    label = "%s (composed from %s)" % (path.name, ", ".join(p.name for p in sources))
    text = _dump(merged, allow_unicode=True, sort_keys=False,
                          default_flow_style=False, width=1000)
    try:
        data = strictyaml.dirty_load(text, YMF_SCHEMA, label=label,
                                     allow_flow_style=_ALLOW_FLOW_STYLE).data
    except strictyaml.YAMLValidationError as exc:
        # The composed document has no file of its own, so the line this
        # error points at exists nowhere on disk. Name the file instead.
        _attribute(exc, text, provenance, sources)
        raise
    # The schema types "composition" as Any, which reads every scalar back
    # as a string. The loader wrote it, so return its own typed copy.
    if composition:
        data["composition"] = composition
    if outputs:
        data["outputs"] = outputs
    _check_document(data)
    return data
