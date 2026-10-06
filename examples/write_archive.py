"""Write a YMF archive the way a solver does, then convert it for ParaView.

    python examples/write_archive.py --outdir /tmp/ymf-demo
    python examples/write_archive.py --outdir /tmp/ymf-inline --inline
    paraview /tmp/ymf-demo/heat.xmf

There is no solver here. The temperature at each output time is the exact
solution from ``examples/heat_equation.ymf``, evaluated at the mesh nodes,
so the example is about the archive and nothing else. A real solver writes
its own arrays in exactly the same way.

Two archives are written, holding the same data:

``heat.ymf``
    One grid per time step -- what a serial run, or a parallel run that
    gathers to global arrays, produces.
``heat_split.ymf``
    Each time step is a collection of two subdomain grids, as if two MPI
    ranks had each written their own piece. ParaView shows the same field.

By default each archive keeps its arrays in an HDF5 file beside it, and
the ``.ymf`` is the metadata that says what each dataset is. With
``--inline`` the arrays go into the archive itself: no HDF5 file, and the
``.ymf`` and ``.xmf`` are each self-contained. That suits small problems,
and solvers that would rather not depend on HDF5. Each archive also carries the
heat-equation problem specification in its ``extra`` section, so the
output records the problem that produced it.

Needs ``numpy`` in addition to ymf's core (``pyyaml``), and ``h5py``
unless ``--inline`` is given. It does *not* need the validation extras:
writing an archive never does.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import yaml

from ymf.archive import (
    add_spatial_step,
    add_uniform_step,
    attribute,
    data_item_for,
    geometry,
    grid,
    new_domain,
    topology,
    write_ymf,
)
from ymf.cli import ymf2xmf

EXAMPLES = Path(__file__).resolve().parent

# From heat_equation.ymf. Kept as constants rather than parsed from the
# formula strings, because the spec's formulas are for people and LLMs to
# read; turning them into code is the symbolic layer's job, which ymf
# doesn't have yet.
KAPPA = 1.0e-3     # m^2/s
T_WALL = 300.0     # K
AMPLITUDE = 10.0   # K


def unit_square_mesh(n):
    """An n x n grid of squares on [0,1]^2, each cut into two triangles.

    Returns ``(nodes, elements)``: ``nodes`` is ``((n+1)^2, 3)`` float64
    with z = 0, because XDMF's ``XYZ`` geometry wants three coordinates
    even in 2D; ``elements`` is ``(2 n^2, 3)`` int32 node indices.
    """
    x = np.linspace(0.0, 1.0, n + 1)
    X, Y = np.meshgrid(x, x, indexing="xy")
    nodes = np.column_stack([X.ravel(), Y.ravel(), np.zeros(X.size)])
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="xy")
    a = (j * (n + 1) + i).ravel()          # lower-left corner of each square
    b, c, d = a + 1, a + n + 2, a + n + 1  # lower-right, upper-right, upper-left
    elements = np.vstack([np.column_stack([a, b, c]),
                          np.column_stack([a, c, d])]).astype(np.int32)
    return nodes, elements


def exact_temperature(nodes, t):
    """T(x,y,t) = 300 + 10 exp(-2 pi^2 kappa t) sin(pi x) sin(pi y)."""
    x, y = nodes[:, 0], nodes[:, 1]
    decay = np.exp(-2.0 * np.pi ** 2 * KAPPA * t)
    return T_WALL + AMPLITUDE * decay * np.sin(np.pi * x) * np.sin(np.pi * y)


def exact_heat_flux(nodes, t):
    """q = -kappa grad T, with a zero z-component so it is an XDMF Vector."""
    x, y = nodes[:, 0], nodes[:, 1]
    scale = -KAPPA * AMPLITUDE * np.pi * np.exp(-2.0 * np.pi ** 2 * KAPPA * t)
    q = np.zeros_like(nodes)
    q[:, 0] = scale * np.cos(np.pi * x) * np.sin(np.pi * y)
    q[:, 1] = scale * np.sin(np.pi * x) * np.cos(np.pi * y)
    return q


class Hdf5Store:
    """Puts each array in an HDF5 file; DataItems reference it."""

    def __init__(self, path):
        import h5py  # only this store needs it

        self.name = path.name
        self.file = h5py.File(path, "w")

    def item(self, dataset, array):
        self.file[dataset] = array
        # data_item_for reads shape and dtype off the array, so Dimensions,
        # DataType and Precision cannot disagree with what is in the file.
        return data_item_for(array, "%s:/%s" % (self.name, dataset))

    def close(self):
        self.file.close()


class InlineStore:
    """Puts each array in the DataItem itself: no HDF5, no other files."""

    def item(self, dataset, array):
        return data_item_for(array, inline=True)

    def close(self):
        pass


def open_store(outdir, stem, inline):
    return InlineStore() if inline else Hdf5Store(outdir / (stem + ".h5"))


def write_serial(outdir, nodes, elements, times, problem, inline=False):
    """One grid per time step. Returns the path of the .ymf written."""
    store = open_store(outdir, "heat", inline)
    domain = new_domain("Mesh Spatial_Domain")
    # The mesh does not move, so it is stored once and every step refers
    # to the same arrays.
    topo = topology("Triangle", len(elements), store.item("elements", elements))
    geom = geometry(store.item("nodes", nodes))

    # A cell-centred field: which half of the square each triangle is in.
    centroid_x = nodes[elements, 0].mean(axis=1)
    side = (centroid_x > 0.5).astype(np.int32)
    side_attr = attribute("side", store.item("side", side), center="Cell")

    for k, t in enumerate(times):
        T = exact_temperature(nodes, t)
        q = exact_heat_flux(nodes, t)
        add_uniform_step(domain, t, topo, geom, [
            attribute("T", store.item("T_t%d" % k, T)),
            attribute("q", store.item("q_t%d" % k, q), attribute_type="Vector"),
            side_attr,
        ])
    store.close()

    ymf_path = outdir / "heat.ymf"
    write_ymf(domain, ymf_path, extra=problem)
    return ymf_path


def split_mesh(nodes, elements):
    """Cut the mesh into left and right halves, renumbering each half's nodes.

    That is what a subdomain looks like to the rank that owns it: its own
    elements, and only the nodes they touch, numbered from zero. Nodes on
    the cut appear in both halves.
    """
    centroid_x = nodes[elements, 0].mean(axis=1)
    parts = []
    for owned in (centroid_x <= 0.5, centroid_x > 0.5):
        local_elements = elements[owned]
        used = np.unique(local_elements)
        renumber = np.full(len(nodes), -1, dtype=np.int32)
        renumber[used] = np.arange(len(used), dtype=np.int32)
        parts.append((nodes[used], renumber[local_elements]))
    return parts


def write_split(outdir, nodes, elements, times, problem, inline=False):
    """Each step is a SpatialCollection of per-subdomain grids."""
    store = open_store(outdir, "heat_split", inline)
    parts = split_mesh(nodes, elements)
    domain = new_domain("Mesh Spatial_Domain")
    meshes = [
        (topology("Triangle", len(local_elements),
                  store.item("elements_p%d" % rank, local_elements)),
         geometry(store.item("nodes_p%d" % rank, local_nodes)))
        for rank, (local_nodes, local_elements) in enumerate(parts)
    ]
    for k, t in enumerate(times):
        grids = []
        for rank, (local_nodes, _) in enumerate(parts):
            T = exact_temperature(local_nodes, t)
            topo, geom = meshes[rank]
            grids.append(grid(topo, geom, [
                attribute("T", store.item("T_p%d_t%d" % (rank, k), T)),
            ], name="subdomain_%d" % rank))
        add_spatial_step(domain, t, grids)
    store.close()

    ymf_path = outdir / "heat_split.ymf"
    write_ymf(domain, ymf_path, extra=problem)
    return ymf_path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--outdir", type=Path, default=Path("ymf-example-output"),
                        help="directory to write into (created if missing)")
    parser.add_argument("-n", type=int, default=16,
                        help="squares per side of the mesh (default 16)")
    parser.add_argument("--inline", action="store_true",
                        help="store the arrays inline in the archive: no HDF5 "
                             "file, and h5py is not needed")
    args = parser.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    # The problem spec rides along in the archive. Plain yaml.safe_load,
    # not ymf.load_ymf: the writer trusts its input and so does not need
    # the strictyaml front-end. Validate specs before a run, not during.
    problem = yaml.safe_load((EXAMPLES / "heat_equation.ymf").read_text())

    nodes, elements = unit_square_mesh(args.n)
    times = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]

    for ymf_path in (write_serial(args.outdir, nodes, elements, times, problem, args.inline),
                     write_split(args.outdir, nodes, elements, times, problem, args.inline)):
        xmf_path = ymf2xmf(ymf_path)
        print("wrote %s and %s" % (ymf_path, xmf_path.name))
    return 0


if __name__ == "__main__":
    sys.exit(main())
