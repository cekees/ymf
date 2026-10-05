"""Read a YMF archive back: metadata from the .ymf, arrays from the HDF5.

    python examples/write_archive.py --outdir /tmp/ymf-demo
    python examples/read_archive.py /tmp/ymf-demo/heat.ymf

This is the post-processing side. It does three things a collaborator's
script would do with an archive:

1. Recover the problem specification the run carried in ``extra``.
2. Follow each DataItem's ``file.h5:/dataset`` reference to the array,
   for a uniform archive or a per-subdomain one alike.
3. Check the data against the spec: the peak temperature excess should
   decay as exp(-2 pi^2 kappa t), so fitting that rate to the archived
   fields recovers kappa, which is compared with the value in the spec.
   write_archive.py stores the exact solution, so the two agree to
   round-off; run against a real solver's archive, the gap measures
   discretization error.

It finishes by reading the derived ``.xmf`` back and confirming it holds
exactly the same domain and problem as the ``.ymf``.

Needs ``numpy`` and ``h5py`` beside ymf's core.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import h5py
import numpy as np

from ymf.archive import canonicalize_domain, read_ymf
from ymf.xdmf import read_xdmf


def load_array(data_item, base_dir):
    """Fetch the array a DataItem points at.

    ``Data`` is ``"<file>:/<dataset>"`` with the file relative to the
    archive's own directory -- that is what lets an archive be moved as a
    directory. ``rsplit`` keeps a colon inside the path (a Windows drive
    letter) from being mistaken for the separator.
    """
    if data_item["Format"] != "HDF":
        raise NotImplementedError("this example only follows HDF references")
    filename, dataset = data_item["Data"].rsplit(":", 1)
    with h5py.File(base_dir / filename, "r") as h5:
        return h5[dataset][()]


def grids_of(step):
    """A uniform step is one grid; a spatial step is several."""
    return step["SpatialCollection"] if "SpatialCollection" in step else [step]


def field(grid, name):
    for attr in grid["Attributes"]:
        if attr["Name"] == name:
            return attr
    raise KeyError("no attribute %r in this grid" % (name,))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("archive", type=Path, help="a .ymf written by write_archive.py")
    args = parser.parse_args(argv)
    base_dir = args.archive.parent

    domain, extra = read_ymf(args.archive)
    problem = extra["Problem"]
    kappa_spec = problem["strong_form"]["coefficients"]["kappa"]["value"]
    t_wall = problem["strong_form"]["boundary_conditions"][0]["value"]
    print("problem: %s" % problem["name"])
    print("kappa in the spec: %g %s"
          % (kappa_spec, problem["strong_form"]["coefficients"]["kappa"]["units"]))

    collection = domain["TimeCollections"][0]
    times, peaks = [], []
    print("\n  time   subdomains   max T [K]")
    for step in collection["Data"]:
        grids = grids_of(step)
        # The peak over the whole domain is the peak over its subdomains.
        peak = max(load_array(field(g, "T")["DataItem"], base_dir).max()
                   for g in grids)
        times.append(step["Time"])
        peaks.append(peak)
        print("  %5.1f   %10d   %9.4f" % (step["Time"], len(grids), peak))

    # ln(T_max - T_wall) = ln(A) - 2 pi^2 kappa t, so the slope gives kappa.
    slope = np.polyfit(times, np.log(np.array(peaks) - t_wall), 1)[0]
    kappa_fit = -slope / (2.0 * math.pi ** 2)
    print("\nkappa fitted from the archived fields: %g" % kappa_fit)
    print("relative difference from the spec:     %.1e"
          % (abs(kappa_fit - kappa_spec) / kappa_spec))

    xmf_path = args.archive.with_suffix(".xmf")
    if xmf_path.exists():
        xmf_domain, xmf_extra = read_xdmf(xmf_path)
        same = (canonicalize_domain(xmf_domain) == canonicalize_domain(domain)
                and xmf_extra == extra)
        print("\n%s holds the same domain and problem: %s" % (xmf_path.name, same))
        if not same:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
