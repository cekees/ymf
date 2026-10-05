# YMF examples

Everything here runs, and `tests/test_examples.py` runs it, so the examples
stay in step with the code.

```bash
pip install -e '.[frontend,examples]'
python examples/validate_spec.py                          # check every spec below
python examples/write_archive.py --outdir /tmp/ymf-demo   # write two archives
python examples/read_archive.py /tmp/ymf-demo/heat.ymf    # read one back
```

## Problem specifications

| File | Shows |
|---|---|
| [poisson.yaml](poisson.yaml) | The smallest complete spec: one weak form, one discretization, an exact solution. v0.1-style bare unknowns. |
| [poisson_v02_enriched.yaml](poisson_v02_enriched.yaml) | The same problem with v0.2 additions: units, a standard name, characteristic scales, a typed boundary condition. |
| [heat_equation.yaml](heat_equation.yaml) | A time-dependent problem: initial and boundary conditions with units, a time integrator, an analytical solution. The archive examples embed this spec. |
| [kovasznay_flow.yaml](kovasznay_flow.yaml) | Branching: Navier–Stokes with two weak forms (Taylor–Hood, stabilized equal-order), a discretization for each, a Reynolds number computed from the scales, and a verification plan against the exact solution. |

## Scripts

**[validate_spec.py](validate_spec.py)**: runs each spec through the
front-end. It checks the schema, prints the branch tree from weak forms to
discretizations, runs the unit check and resolves dimensionless numbers.
At the end it breaks one spec on purpose to show what a unit error looks
like. Needs `ymf[frontend]`.

**[write_archive.py](write_archive.py)**: writes an archive the way a solver
would, on a 16×16 triangulated unit square over six output times. There is
no solver: the fields are the heat equation's exact solution, so the
example is about the archive format alone. It writes:

```text
heat.h5  heat.ymf  heat.xmf                    one grid per step (serial, or gathered)
heat_split.h5  heat_split.ymf  heat_split.xmf  two subdomain grids per step (per-rank)
```

Open either `.xmf` in ParaView or VisIt to see the temperature decay. The
script needs `numpy` and `h5py` (`ymf[examples]`) but not the front-end:
writing an archive never does.

**[read_archive.py](read_archive.py)**: the post-processing side. It
recovers the problem spec from the archive, follows each field's HDF5
reference, fits the decay rate of the peak temperature to recover κ, and
compares κ with the value in the spec. It then confirms the `.xmf` holds
exactly what the `.ymf` does.

## Notebooks

[notebooks/](notebooks/) holds the original design notebook and its port to
the `ymf` package, with [COMPARISON.md](notebooks/COMPARISON.md) recording
that both produce byte-identical XDMF. They predate the archive core
(`ymf.archive`) and show the earlier `write_xdmf` path. For new work, start
from the scripts above.
