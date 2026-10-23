"""Command line entry points.

Only ``ymf2xmf`` so far: it turns a ``.ymf`` archive into the ``.xmf`` that
ParaView and VisIt read.

This exists because of the campaign's decision that a solver writes YMF and
nothing else. XDMF is a *consumer* of a YMF archive, not a second format
the solver has to keep in step -- so the conversion happens here, on
demand, rather than in the writer. A solver that emits both has two
representations to keep consistent and will eventually fail to.

Deliberately in the ``pyyaml``-only core: this reads YAML and writes XML,
both of which :mod:`ymf.archive` and :mod:`ymf.xdmf` do without the
validation stack. Running the converter must not require ``ymf[spec]``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Sequence

from ymf.archive import (YmfArchiveError, canonicalize_domain, read_ymf, select_domain,
                         validate_domain)
from ymf.xdmf import write_xdmf

__all__ = ["archive2xmf", "describe", "ymf2xmf", "ymf2xmf_main"]


def _default_output(source: Path) -> Path:
    """``run.ymf`` -> ``run.xmf``, keeping any other suffixes intact."""
    return source.with_suffix(".xmf")


def ymf2xmf(
    source: str | Path,
    output: Optional[str | Path] = None,
    *,
    validate: bool = True,
    collections: Optional[Sequence[str]] = None,
    fields: Optional[Sequence[str]] = None,
    steps: Optional[str] = None,
) -> Path:
    """Write the ``.xmf`` for a ``.ymf`` archive. Returns the output path.

    Parameters
    ----------
    source:
        The ``.ymf`` archive.
    output:
        Where to write. Defaults to ``source`` with a ``.xmf`` suffix.
    validate:
        Check the archive's structure before converting. On by default: a
        malformed archive is better reported here, with the offending
        field named, than as a file a viewer silently mis-draws.
    collections, fields, steps:
        Convert only these: see :func:`ymf.archive.select_domain`. The
        ``.xmf`` then describes a part of the archive, which is what a
        viewer opening one space, one field or one instant of a long run
        wants; the ``.ymf`` is untouched.

    The archive's ``extra`` payload -- ``Problem``, ``solution_paths``,
    ``vvuq`` and anything else with no XDMF equivalent -- is carried into
    the ``.xmf`` as the base64 ``<Information Name="YMF">`` element that
    :func:`ymf.xdmf.read_xdmf` can recover. Generic XDMF consumers parse
    and ignore it, so nothing is lost by converting and nothing is
    imposed on the viewer.
    """
    source = Path(source)
    domain, extra = read_ymf(source)
    if not domain:
        raise YmfArchiveError(
            "%s holds no domain, so there is nothing to convert" % (source,))
    domain = validate_domain(domain) if validate else canonicalize_domain(domain)
    if collections is not None or fields is not None or steps is not None:
        domain = select_domain(domain, collections=collections, fields=fields, steps=steps)
    destination = Path(output) if output is not None else _default_output(source)
    write_xdmf(domain, destination, ymf_extra=extra)
    return destination


def archive2xmf(
    source: str | Path,
    outdir: Optional[str | Path] = None,
    keys: Optional[Sequence[str]] = None,
    *,
    validate: bool = True,
) -> List[Path]:
    """Write one ``.xmf`` per output of an archive. Returns the paths.

    Each is named ``<stem>_<digest>.xmf`` after the output's key (beside
    its ``<stem>_<digest>.h5``, if it has one), and
    carries the realized spec and that output's record (everything but
    its arrays) as the ``<Information Name="YMF">`` element.
    """
    from ymf.archive import read_document
    from ymf.closure import archive_stem, file_stem
    source = Path(source)
    document = read_document(source)
    outputs = document.get("outputs") or {}
    if not outputs:
        raise YmfArchiveError("%s holds no outputs, so there is nothing to convert" % (source,))
    unknown = [k for k in keys or () if k not in outputs]
    if unknown:
        raise YmfArchiveError("%s has no output %s" % (source, ", ".join(unknown)))
    spec = {k: v for k, v in document.items() if k != "outputs"}
    outdir = Path(outdir) if outdir is not None else source.parent
    written = []
    for key in keys or list(outputs):
        output = outputs[key]
        domain = output.get("approximation")
        if not domain:
            continue
        domain = validate_domain(domain) if validate else canonicalize_domain(domain)
        record = {k: v for k, v in output.items() if k != "approximation"}
        destination = outdir / (file_stem(archive_stem(source), key) + ".xmf")
        write_xdmf(domain, destination, ymf_extra=dict(spec, output=key, record=record))
        written.append(destination)
    return written


def _is_archive(source: Path) -> bool:
    """An archive in the current form: a spec with ``outputs``."""
    import re
    with open(source, encoding="utf-8") as f:
        for line in f:
            if re.match(r"(ymf|outputs):", line):
                return True
            if line.startswith("ymf_archive_version:"):
                return False
    return False


def _fields_of(collection) -> List[str]:
    """Every field name any grid of the collection carries, in first-seen order."""
    names: List[str] = []
    for step in collection.get("data", []):
        grids = step.get("spatial_collection", [step])
        for g in grids:
            for a in g.get("attributes", []):
                if a["name"] not in names:
                    names.append(a["name"])
    return names


def describe(domain) -> List[str]:
    """What an archive holds, one line per collection, for ``--show`` and ``--verbose``.

    Each collection is numbered so the interactive prompt can refer to it.
    """
    lines = []
    for i, collection in enumerate(domain.get("time_collections", []), start=1):
        steps = collection["data"]
        times = [step["time"] for step in steps]
        span = ("no steps" if not times
                else "1 step at t=%g" % times[0] if len(times) == 1
                else "%d steps, t=%g..%g" % (len(times), times[0], times[-1]))
        detail = ""
        if steps:
            first = steps[0]
            if "spatial_collection" in first:
                detail = ", %d subdomains" % len(first["spatial_collection"])
                first = first["spatial_collection"][0]
            detail += ", %s topology" % first["topology"]["type"]
            fields = _fields_of(collection)
            if fields:
                detail += ", fields: %s" % ", ".join(fields)
        lines.append("  %d. %s: %s%s" % (i, collection["name"], span, detail))
    return lines


_describe = describe


def _pick(prompt: str, options: Sequence[str], ask=None) -> Optional[List[str]]:
    """Ask for some of ``options`` by number; blank keeps them all (None)."""
    ask = ask or input
    for i, option in enumerate(options, start=1):
        print("  %d. %s" % (i, option))
    answer = ask(prompt).strip()
    if not answer:
        return None
    chosen = []
    for token in answer.replace(",", " ").split():
        try:
            index = int(token)
        except ValueError:
            raise YmfArchiveError("%r is not a number" % (token,)) from None
        if not 1 <= index <= len(options):
            raise YmfArchiveError("%d is not one of 1..%d" % (index, len(options)))
        if options[index - 1] not in chosen:
            chosen.append(options[index - 1])
    return chosen


def interactive_selection(domain, ask=None):
    """Ask which collections, fields and steps to convert.

    Returns ``(collections, fields, steps)`` as :func:`ymf2xmf` takes
    them, each ``None`` where everything was kept. ``ask`` defaults to
    ``input`` (looked up at call time, so a test can replace it); tests
    pass their own.
    """
    ask = ask or input
    collections = domain.get("time_collections", [])
    names = [c["name"] for c in collections]
    print("Collections (blank for all):")
    chosen_names = _pick("collections> ", names, ask)
    kept = [c for c in collections if chosen_names is None or c["name"] in chosen_names]
    fields: List[str] = []
    for c in kept:
        for f in _fields_of(c):
            if f not in fields:
                fields.append(f)
    chosen_fields = None
    if fields:
        print("Fields (blank for all):")
        chosen_fields = _pick("fields> ", fields, ask)
    n_steps = max((len(c.get("data", [])) for c in kept), default=0)
    steps = None
    if n_steps > 1:
        answer = ask("Steps of %d, as indices or slices, e.g. -1, 0:10, ::5 (blank for all)> "
                     % (n_steps,)).strip()
        steps = answer or None
    return chosen_names, chosen_fields, steps


def ymf2xmf_main(argv: Optional[Sequence[str]] = None) -> int:
    """``ymf2xmf`` console entry point."""
    parser = argparse.ArgumentParser(
        prog="ymf2xmf",
        description="Convert a YMF archive into an XDMF (.xmf) file that "
                    "ParaView and VisIt can open. The heavy data is not "
                    "copied -- the .xmf references the same HDF5 datasets "
                    "the .ymf does, so it must sit beside them.")
    parser.add_argument("archive", help="the .ymf archive to convert")
    parser.add_argument("-o", "--output",
                        help="output path (default: the archive with a .xmf suffix); "
                             "for an archive with outputs, the directory to write to")
    parser.add_argument("--key", action="append",
                        help="convert only this output (repeatable; default: all)")
    parser.add_argument("--no-validate", action="store_true",
                        help="skip the structural check before converting")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="describe what the archive holds")
    parser.add_argument("--show", action="store_true",
                        help="describe what the archive holds and stop; convert nothing")
    parser.add_argument("-c", "--collection", action="append", metavar="NAME",
                        help="convert only this time collection, by name or part of one "
                             "(repeatable; default: all)")
    parser.add_argument("-f", "--field", action="append", metavar="NAME",
                        help="convert only this field (repeatable; default: all)")
    parser.add_argument("-s", "--steps", metavar="SPEC",
                        help="convert only these steps: indices and slices, e.g. -1, "
                             "0:10, ::5 (default: all)")
    parser.add_argument("-i", "--interactive", action="store_true",
                        help="ask which collections, fields and steps to convert")
    args = parser.parse_args(argv)

    source = Path(args.archive)
    if not source.exists():
        parser.error("no such file: %s" % (source,))

    try:
        if _is_archive(source):
            for destination in archive2xmf(source, args.output, args.key,
                                           validate=not args.no_validate):
                print("wrote %s" % (destination,))
            return 0
        collections, fields, steps = args.collection, args.field, args.steps
        if args.verbose or args.show or args.interactive:
            domain, _ = read_ymf(source)
            domain = canonicalize_domain(domain)
            print("%s" % (source,))
            if not args.interactive:
                for line in describe(domain):
                    print(line)
            if args.show:
                return 0
            if args.interactive:
                collections, fields, steps = interactive_selection(domain)
        destination = ymf2xmf(source, args.output,
                               validate=not args.no_validate,
                               collections=collections, fields=fields, steps=steps)
    except YmfArchiveError as exc:
        # A structural problem in the archive, named. Not a traceback:
        # the user asked to convert a file, and the file is the problem.
        print("ymf2xmf: %s" % (exc,), file=sys.stderr)
        return 1
    except OSError as exc:
        print("ymf2xmf: %s" % (exc,), file=sys.stderr)
        return 1

    print("wrote %s" % (destination,))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(ymf2xmf_main())
