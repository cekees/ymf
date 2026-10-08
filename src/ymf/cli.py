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

from ymf.archive import YmfArchiveError, canonicalize_domain, read_ymf, validate_domain
from ymf.xdmf import write_xdmf

__all__ = ["archive2xmf", "ymf2xmf", "ymf2xmf_main"]


def _default_output(source: Path) -> Path:
    """``run.ymf`` -> ``run.xmf``, keeping any other suffixes intact."""
    return source.with_suffix(".xmf")


def ymf2xmf(
    source: str | Path,
    output: Optional[str | Path] = None,
    *,
    validate: bool = True,
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


def _describe(domain) -> List[str]:
    """A short summary of what an archive holds, for --verbose."""
    lines = []
    for collection in domain.get("TimeCollections", []):
        steps = collection["Data"]
        times = [step["Time"] for step in steps]
        span = ("no steps" if not times
                else "1 step at t=%g" % times[0] if len(times) == 1
                else "%d steps, t=%g..%g" % (len(times), times[0], times[-1]))
        detail = ""
        if steps:
            first = steps[0]
            if "SpatialCollection" in first:
                detail = ", %d subdomains" % len(first["SpatialCollection"])
                first = first["SpatialCollection"][0]
            fields = [a["Name"] for a in first.get("Attributes", [])]
            detail += ", %s topology" % first["Topology"]["Type"]
            if fields:
                detail += ", fields: %s" % ", ".join(fields)
        lines.append("  %s: %s%s" % (collection["Name"], span, detail))
    return lines


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
        if args.verbose:
            domain, _ = read_ymf(source)
            print("%s" % (source,))
            for line in _describe(canonicalize_domain(domain)):
                print(line)
        destination = ymf2xmf(source, args.output,
                               validate=not args.no_validate)
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
