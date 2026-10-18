"""Files a tool wrote: the ``generated`` record, and whether it still holds.

A step such as ``pro weak`` or ``pro disc`` writes its output as a new file
that ``extends`` its input, so the person who reviews it reads one file and
the steps stay separate. The file says what wrote it::

    extends: poisson.ymf
    generated:
      by: pro weak (proteus 2.0.0.dev0, ymf 0.2.1.dev)
      reads: model_without_weak_forms
      from: 9c41e0a2b7d3       # sha256[:12] of what the step read
      content: 51f0c2aa19e8    # sha256[:12] of what it wrote
    model:
      weak_forms: [...]

Two questions then have cheap answers:

* **stale**: the input changed since the step ran. ``from`` hashes the part
  of the composed parent the step read (``reads`` names which part), so an
  edit upstream that the step did not read does not count.
* **edited**: someone changed what the step wrote. ``content`` hashes the
  file's own content -- what it contributes, its ``remove`` list included. An edited file
  is legitimate (it is the person's now); a tool must not overwrite it.

Like ``extends`` and ``ymf_state``, ``generated`` describes a file, not what
it specifies: it is not inherited, so a file that extends a generated one
is not itself generated. An archive does not carry it either.
"""

from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

__all__ = ["READS", "DIGITS", "digest", "slice_of", "content_of", "record",
           "write", "Status", "check"]

#: hex digits of a digest kept in a record (as in output keys)
DIGITS = 12


def _model_without_weak_forms(doc: Dict[str, Any]) -> Dict[str, Any]:
    return {"model": {k: v for k, v in (doc.get("model") or {}).items()
                      if k != "weak_forms"}}


def _model_and_analytical(doc: Dict[str, Any]) -> Dict[str, Any]:
    return {"model": doc.get("model") or {},
            "analytical": (doc.get("solution_paths") or {}).get("analytical") or []}


#: What a step may read, by name: the part of its parent document whose
#: change makes the step's output stale. pro weak reads the posed model
#: (its weak forms are what it writes); pro disc reads the model with its
#: weak forms, and the analytical solutions its outputs are checked against.
READS: Dict[str, Callable[[Dict[str, Any]], Dict[str, Any]]] = {
    "model_without_weak_forms": _model_without_weak_forms,
    "model_and_analytical": _model_and_analytical,
}

#: the file-level keys that are not content. ``remove`` is a file-level key
#: too, but what it drops is part of what the file contributes, so it counts.
_NOT_CONTENT = ("extends", "ymf_state", "composition", "generated")


def _plain(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value


def digest(value: Any) -> str:
    """sha256[:DIGITS] of plain data, in a canonical JSON form."""
    text = json.dumps(_plain(value), sort_keys=True, ensure_ascii=False,
                      separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:DIGITS]


def slice_of(doc: Dict[str, Any], reads: str) -> Dict[str, Any]:
    """The part of ``doc`` a step that ``reads`` this depends on."""
    if reads not in READS:
        raise ValueError("generated.reads %r is not one of %s"
                         % (reads, ", ".join(sorted(READS))))
    return READS[reads](doc)


def content_of(data: Dict[str, Any]) -> Dict[str, Any]:
    """A file's own content: what it contributes, ``remove`` included."""
    return {k: v for k, v in data.items() if k not in _NOT_CONTENT}


def record(by: str, reads: str, parent: Dict[str, Any],
           content: Dict[str, Any]) -> Dict[str, str]:
    """The ``generated`` block for a step that read ``parent`` and wrote ``content``."""
    return {"by": by, "reads": reads, "from": digest(slice_of(parent, reads)),
            "content": digest(content_of(content))}


def write(path: str | Path, extends: str, by: str, reads: str,
          parent: Dict[str, Any], content: Dict[str, Any],
          state: Optional[str] = None) -> Dict[str, str]:
    """Write a generated file: ``extends``, its ``generated`` record, its content.

    ``parent`` is the composed document the step read (``extends`` names
    its file, relative to ``path``). ``content`` is what the step adds, as
    plain data; it is checked against the partial schema by reading the
    file back, and ``content`` in the record is the digest of what was read
    back, so a later :func:`check` compares like with like. ``state``, when
    given, is written as the file's ``ymf_state``: the state the step says
    the composed document reaches (``pro disc`` writes ``discretized``),
    which the loader then checks. Like ``extends``, it is not content.
    Returns the record.
    """
    from ymf.compose import _dump, _read
    path = Path(path)
    own = content_of(copy.deepcopy(_plain(content)))

    if state is not None:
        from ymf.state import STATES
        if state not in STATES:
            raise ValueError("ymf_state %r is not one of %s" % (state, ", ".join(STATES)))

    def text(generated: Optional[Dict[str, str]]) -> str:
        head = {"extends": extends}
        if state is not None:
            head["ymf_state"] = state
        if generated is not None:
            head["generated"] = generated
        return _dump(dict(head, **own), allow_unicode=True, sort_keys=False,
                     default_flow_style=False, width=1000)

    path.write_text(text(None), encoding="utf-8")
    data, _ = _read(path)
    generated = {"by": by, "reads": reads,
                 "from": digest(slice_of(parent, reads)),
                 "content": digest(content_of(data))}
    path.write_text(text(generated), encoding="utf-8")
    return generated


@dataclass
class Status:
    """What :func:`check` found for one file."""

    path: Path
    record: Optional[Dict[str, str]]          # None: not a generated file
    stale: bool = False                       # its input changed since
    edited: bool = False                      # someone changed what it holds
    reasons: List[str] = field(default_factory=list)

    @property
    def current(self) -> bool:
        return self.record is not None and not (self.stale or self.edited)


def check(path: str | Path) -> Status:
    """Whether a generated file still holds: neither stale nor edited.

    Reads the file, composes what it ``extends`` (without checking that
    document's declared state), and compares both digests with the record.
    A file with no ``generated`` record gives a Status whose ``record`` is
    None.
    """
    from ymf.compose import _read, load_composed
    path = Path(path).resolve()
    data, _ = _read(path)
    generated = data.get("generated")
    status = Status(path=path, record=None if generated is None else dict(generated))
    if generated is None:
        return status
    missing = [k for k in ("by", "reads", "from", "content") if not generated.get(k)]
    if missing:
        raise ValueError("%s: generated lacks %s" % (path, ", ".join(missing)))
    if "extends" not in data:
        raise ValueError("%s: a generated file extends the document it was "
                         "generated from, and this one extends nothing" % (path,))
    parent = load_composed(path.parent / data["extends"], check_state=False)
    parent.pop("outputs", None)
    now = digest(slice_of(parent, generated["reads"]))
    if now != generated["from"]:
        status.stale = True
        status.reasons.append(
            "%s changed in what %s read (%s): %s, it was %s"
            % (data["extends"], generated["by"], generated["reads"], now, generated["from"]))
    held = digest(content_of(data))
    if held != generated["content"]:
        status.edited = True
        status.reasons.append(
            "its content was edited since %s wrote it: %s, it was %s"
            % (generated["by"], held, generated["content"]))
    return status
