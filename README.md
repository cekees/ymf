# YMF — YAML Modeling Format

YMF is a plain-YAML record of a computational physics problem: what it
models and why, the equations, how they were discretized and solved, and
what came out. It is meant to be written and reviewed by people and LLMs
together, and read and written by solvers.

It has two parts that share one format:

| | Problem specification | Archive |
|---|---|---|
| **What** | The problem from physics down to solver settings: processes and assumptions, strong form, weak forms, discretizations, verification plan | A solver's output: which mesh, time and field each HDF5 dataset holds |
| **Written by** | a person, an LLM, or both | the solver |
| **Checked by** | `ymf.schema` (strictyaml) and `ymf.units` (pint): schema, units, dimensional consistency | `ymf.archive`: structure and array shapes only, cheaply |
| **Depends on** | `ymf[frontend]` | `pyyaml` alone |
| **Guide** | [docs/problem-specification.md](docs/problem-specification.md) | [docs/archive-format.md](docs/archive-format.md) |

An archive can carry the specification that produced it, so a run's output
records the problem it solved.

## Why

Between "simulate this system" and "run this solver" sit a dozen decisions:
which physics to keep, which unknowns, which weak form, which elements,
which stabilization. They usually live in someone's head, a notebook, or a
solver input file that records the outcome but not the reasoning. YMF makes
each stage a structured, reviewable block, and tags each one with its
**provenance**: whether an LLM proposed it, a person wrote it, or a person
edited an LLM's proposal. One document can hold several **branches**
(competing weak forms and discretizations of the same problem) alongside
the plan for verifying them against each other. Branches can also be
separate files: a **model** (the physics and equations, often most of the
work for multiphysics) is written once and extended, file by file, into
many well-posed problems.

On the output side, solver archives have usually been XDMF XML written by
hand-built element trees, with no data model behind them and nothing
checking that the metadata matches the arrays. The YMF archive is a small,
explicit data model that converts one-to-one to XDMF, so ParaView and VisIt
still read it.

## Status

Pre-release (0.2.x). The repository is private while the format settles.

| Piece | State |
|---|---|
| Archive data model, YAML I/O, structural checks (`ymf.archive`) | **Implemented.** Tested against real Proteus output; used by Proteus's `ymf_io` branch, aimed at Proteus 2.0 |
| `.ymf` → `.xmf` conversion (`ymf2xmf`, `ymf.xdmf`) | **Implemented.** Lossless round trip; output opens in ParaView. Reads XDMF from other tools too, refusing by name what the model can't hold ([docs/xdmf-model.yaml](docs/xdmf-model.yaml)) |
| Self-contained archives: arrays inline, no HDF5 | **Implemented.** `data_item_for(array, inline=True)`; `write_archive.py --inline` |
| Problem-spec schema v0.2 (`ymf.schema`) | **Implemented.** See the [known gaps](docs/problem-specification.md#known-gaps) |
| Composing a spec from several files: a model extended into problems (`ymf.compose`) | **Implemented.** Overrides are recorded in the composed document |
| Units, scales, dimensionless numbers (`ymf.units`) | **Implemented**, apart from rewriting the PDE in dimensionless form |
| Manning's n lookup and data-source catalog (`ymf.data_sources`) | Implemented; data fetching is a stub |
| Symbolic layer: parsing the strong and weak forms (sympy/ibvp) | **Planned.** The mathematics in a spec is free text today |
| Dispatch to several solver backends | Planned |
| Multi-component coupling through BMI (`System:` block) | Deferred to v0.3 |

## Quick start

```bash
git clone git@github.com:cekees/ymf.git && cd ymf
pip install -e '.[dev]'        # front-end, examples and tests
pytest

python examples/validate_spec.py
python examples/write_archive.py --outdir /tmp/ymf-demo
python examples/read_archive.py /tmp/ymf-demo/heat.ymf
ymf2xmf -v /tmp/ymf-demo/heat.ymf      # then open /tmp/ymf-demo/heat.xmf in ParaView
```

A solver that only writes archives installs ymf with no extras, which
brings in `pyyaml` and nothing else. See [examples/README.md](examples/README.md) for
what each example demonstrates.

## A first look

A specification, abridged from
[examples/kovasznay_flow.yaml](examples/kovasznay_flow.yaml):

```yaml
Problem:
  name: "Kovasznay flow, Re = 40"
  physical_model:
    provenance: llm_derived
    processes: [incompressible_flow]
    assumptions: [steady, newtonian, constant_density, two_dimensional]
  characteristic_scales:
    length: {value: 1.0, units: m}
    velocity: {value: 1.0, units: m/s}
    density: {value: 1.0, units: kg/m3}
    viscosity: {value: 0.025, units: Pa*s}
  dimensionless_numbers:
    Re: {formula: "rho * U * L / mu"}           # resolved to 40.0
  strong_form:
    unknowns: [{name: v, units: m/s}, {name: p, units: Pa}]
    strong_form_expression: |-
      ρ (v·∇)v - μ Δv + ∇p = 0,   ∇·v = 0   in Ω
    ...
  weak_forms:
    - label: "mixed"                            # needs an inf-sup stable pair
    - label: "stabilized_equal_order"           # SUPG/PSPG
solution_paths:
  discretizations:
    - name: "taylor_hood"
      from_weak_form: "mixed"
      finite_element: {velocity: {family: CG, order: 2}, pressure: {family: CG, order: 1}}
    - name: "equal_order_p1"
      from_weak_form: "stabilized_equal_order"
      finite_element: {velocity: {family: CG, order: 1}, pressure: {family: CG, order: 1}}
```

Writing an archive, from
[examples/write_archive.py](examples/write_archive.py):

```python
from ymf.archive import (new_domain, add_uniform_step, topology, geometry,
                         attribute, data_item_for, write_ymf)

domain = new_domain("Mesh Spatial_Domain")
topo = topology("Triangle", len(elements), data_item_for(elements, "heat.h5:/elements"))
geom = geometry(data_item_for(nodes, "heat.h5:/nodes"))
for k, t in enumerate(times):
    add_uniform_step(domain, t, topo, geom,
                     [attribute("T", data_item_for(T[k], "heat.h5:/T_t%d" % k))])
write_ymf(domain, "heat.ymf", extra={"Problem": problem})
```

`data_item_for` reads shape, type and precision off the array, so the
metadata can't disagree with the data.

## Design principles

- **Composable.** A new branch adds weak forms or discretizations and never
  edits existing ones.
- **Provenance everywhere.** Each decision records whether it was
  `llm_derived`, `human_specified`, `human_edited` or `proteus_derived`.
- **A trust boundary.** Validation of human/LLM input is thorough and runs
  once, in the front-end. The solver trusts validated input and checks
  only structure, cheaply. The archive core therefore depends on `pyyaml`
  alone. strictyaml measured ~27× slower than libyaml-backed pyyaml
  on archive metadata, too slow for a solver's write path.
- **XDMF is a consumer, not a dependency.** The archive maps one-to-one
  onto XDMF, and `.xmf` files are derived from it on demand, carrying the
  YMF-only content in a standard `<Information>` element that viewers
  ignore.
- **Additive evolution.** Schema changes add optional fields where
  possible. The one breaking change so far (structured `unknowns`) comes
  with `normalize_unknowns()`.

## Repository layout

```text
src/ymf/
  archive.py        archive data model, YAML I/O, structural checks   (pyyaml)
  xdmf.py           domain <-> XDMF conversion                       (stdlib)
  cli.py            the ymf2xmf command
  normalize.py      normalize_unknowns                               (stdlib)
  schema.py         problem-spec schema                              (strictyaml, ymf[spec])
  compose.py        `extends`: models, problems built on them        (strictyaml, ymf[spec])
  units.py          unit checks, scales, dimensionless numbers       (pint, ymf[units])
  data_sources.py   Manning's n lookup, data-source catalog          (stdlib)
docs/               guides, design drafts, background reports -- start at docs/README.md
examples/           specs and runnable scripts -- start at examples/README.md
tests/              pytest suite, including real Proteus .xmf fixtures and the examples
```

`import ymf` loads only the pyyaml/stdlib modules. `ymf.load_ymf`,
`ymf.check_units` and the other front-end names load strictyaml or pint on
first use, so an archive-only install never needs them.

## License

MIT (see `LICENSE`).
