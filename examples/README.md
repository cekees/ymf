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

### Composed from several files: [navier_stokes/](navier_stokes/)

One model, extended into three problems with `extends`. See
[Composing documents](../docs/problem-specification.md#composing-documents-models-and-problems).

| File | Shows |
|---|---|
| [navier_stokes_model.yaml](navier_stokes/navier_stokes_model.yaml) | A `kind: model` file: the physics, the equations, the unknowns, and the units of the coefficients, with no domain, conditions or values. Not a problem on its own. |
| [planar_couette.yaml](navier_stokes/planar_couette.yaml) | Extends the model into shear-driven flow between plates: geometry, periodic BCs, coefficient values, the exact (linear) solution, a Galerkin weak form and Taylor–Hood, which should reproduce the exact solution on any grid. |
| [plane_poiseuille.yaml](navier_stokes/plane_poiseuille.yaml) | Extends the model into pressure-driven flow: the parabolic exact solution, a stabilized weak form and equal-order P1, with a grid-convergence plan. |
| [plane_poiseuille_re100.yaml](navier_stokes/plane_poiseuille_re100.yaml) | Extends Poiseuille to change the viscosity and the driving gradient only (Re 1 → 100). Loading it records each overridden value in `composition.overrides`. |

## Scripts

**[validate_spec.py](validate_spec.py)**: runs each spec through the
front-end. It checks the schema, prints the branch tree from weak forms to
discretizations, runs the unit check and resolves dimensionless numbers.
For a composed spec it lists the source files and any overrides.
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

Open either `.xmf` in ParaView or VisIt to see the temperature decay. With
`--inline` the arrays go into the archives themselves: only the `.ymf` and
`.xmf` files are written, each self-contained, and h5py isn't needed. The
script needs `numpy`, and `h5py` without `--inline` (`ymf[examples]`), but
not the front-end: writing an archive never does.

**[read_archive.py](read_archive.py)**: the post-processing side. It
recovers the problem spec from the archive, gets each field's values
(following its HDF5 reference, or reading them inline), fits the decay rate of the peak temperature to recover κ, and
compares κ with the value in the spec. It then confirms the `.xmf` holds
exactly what the `.ymf` does.

## Notebooks

[notebooks/](notebooks/) holds the original design notebook and its port to
the `ymf` package, with [COMPARISON.md](notebooks/COMPARISON.md) recording
that both produce byte-identical XDMF. They predate the archive core
(`ymf.archive`) and show the earlier `write_xdmf` path. For new work, start
from the scripts above.
