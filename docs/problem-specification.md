# Writing a YMF problem specification

A problem specification is a YAML document that describes a PDE problem
from its physics down to a solver configuration, so that a person, an LLM
and a solver can all work from the same record. This guide covers the
current schema (v0.2) as `ymf.schema` implements it, using the documents in
[`examples/`](../examples/) to illustrate. For the design history and
rationale, see [ymf-schema.md](ymf-schema.md) (v0.1) and
[ymf-schema-v0.2-delta.md](ymf-schema-v0.2-delta.md).

- [The idea](#the-idea)
- [Document structure](#document-structure)
- [Branching: one problem, several solution paths](#branching-one-problem-several-solution-paths)
- [Provenance](#provenance)
- [Units, scales and dimensionless numbers](#units-scales-and-dimensionless-numbers)
- [Validating a document](#validating-a-document)
- [What is checked and what is free text](#what-is-checked-and-what-is-free-text)
- [Known gaps](#known-gaps)

## The idea

Getting from "simulate this physical system" to "run this solver" passes
through decisions that usually live only in someone's head or a notebook:

```text
Stage 0  physical processes      which physics to include, which assumptions
Stage 1  strong form             the PDEs, domain, boundary and initial conditions
Stage 2  solution variables      which unknowns (u, p  vs  ρu, p ...)
Stage 3  weak form               the variational statement and its function spaces
Stage 4  discretization          element family and order, stabilization
Stage 5  solver                  linear/nonlinear, time stepping, tolerances
```

A YMF spec records each stage as a structured block, and records **who
made each decision** (`provenance`). The intended workflow has an LLM
propose the early stages (model selection, formulation, choice of
unknowns) and draft the later ones, with a person reviewing at each step.
The document is the reviewable artifact.

## Document structure

```yaml
Problem:                     # required
  name: ...                  # required
  physical_model: ...        # stage 0
  strong_form: ...           # stages 1-2
  weak_forms: [...]          # stage 3, one or more
  characteristic_scales: ... # optional: for non-dimensionalization
  dimensionless_numbers: ... # optional: formulas over the scales
  units: ...                 # optional: {convention: SI|CGS|US|dimensionless}
  unit_systems: ...          # optional
  mesh_generation: ...       # optional: where a mesh comes from (DEM, filters)

solution_paths:              # required
  analytical: [...]          # exact solutions, possibly empty
  discretizations: [...]     # stages 4-5, possibly empty

vvuq: ...                    # optional, free-form: verification/validation plans and results
archive: ...                 # optional, free-form: reserved for output configuration
Xdmf: ...                    # optional, free-form: legacy, from the v0.1 draft
```

The smallest complete example is
[`examples/poisson.yaml`](../examples/poisson.yaml). The key blocks:

**`physical_model`**: `provenance` (required), plus lists of `processes`
and `assumptions`, and optionally a `source_document`.

**`strong_form`**: `unknowns`, `equation_formulation`,
`strong_form_expression`, `domain` and `boundary_regions` are required.
Unknowns are either bare names (`[u, p]`) or structured entries with
`units` and an optional `std_name`. `ymf.normalize_unknowns()` turns either
form into the structured one. Optional: `initial_conditions`,
`boundary_conditions`, `coefficients`, `known_analytical_solution`,
`dimensional_check`.

```yaml
initial_conditions:
  - field: T
    type: function          # constant | function | data_source
    formula: "300 + 10 sin(πx) sin(πy)"
    units: K
boundary_conditions:
  - region: walls           # a boundary_regions[].name
    variable: T
    type: dirichlet         # dirichlet | neumann | robin | periodic | hydrograph
    value: 300.0
    units: K
```

**`weak_forms`**: each has a `label`, `provenance`, a `derivation`, the
trial and test `solution_spaces`, and the `bilinear` and `linear` forms.
Optional notes record inf-sup stability, existence, uniqueness and
regularity, plus the stabilization method.

**`discretizations`**: each has a `name`, `provenance`, `solver`, and
usually `from_weak_form` and `finite_element`. `finite_element` takes
`fields`, or `velocity` and `pressure` for mixed problems, each as
`{family: CG|DG|RT|BDM|RTF|N1curl, order: n}`. The `solver` block takes
`type` (`linear | nonlinear | time_marching | saddle_point`),
`linear_solver` (`petsc | mumps | umfpack | superlu`), `nonlinear_solver`
(`newton | broyden | line_search`), `tolerance` and `max_iterations`.
`discretization` is a free-form block for anything else, such as the time
integrator in `heat_equation.yaml`.

## Branching: one problem, several solution paths

The schema exists to let one document hold *competing* approaches to the
same problem, and to keep track of how they relate. Each discretization
names the weak form it implements with `from_weak_form` (a label, or an
integer index into `weak_forms`), and the `vvuq` block names the branches
it compares. [`examples/kovasznay_flow.yaml`](../examples/kovasznay_flow.yaml)
has this tree:

```text
strong form: incompressible Navier-Stokes, primitive variables
  weak form 'mixed' ─────────────────── discretization 'taylor_hood'     velocity CG2, pressure CG1
  weak form 'stabilized_equal_order' ── discretization 'equal_order_p1'  velocity CG1, pressure CG1
  analytical 'kovasznay_exact' ◄─────── vvuq.verification.reference
```

Adding a branch adds entries; it never edits existing ones. That is what
makes a document safe for several contributors, or an LLM and a person, to
extend.

## Provenance

Every decision block carries `provenance`, one of:

| value | meaning |
|---|---|
| `llm_derived` | proposed by an LLM, not yet confirmed by a person |
| `human_specified` | written by a person |
| `human_edited` | LLM-proposed, then changed by a person |
| `proteus_derived` | derived by a solver (e.g. from a Proteus model) |

It is required on `physical_model`, `strong_form` (as
`unknown_provenance`), every weak form, and every solution path. Its
purpose is to calibrate trust: a reviewer can see at a glance which parts
of a model a person has looked at.

## Units, scales and dimensionless numbers

Units are optional everywhere. A document with none still passes. Where
they are given, `ymf.check_units(doc)` (pint-backed) checks that every unit
string parses, and that each boundary and initial condition's units are
dimensionally compatible with the unknown it constrains:

```text
[error] boundary_conditions[region='walls', variable='T']: units 'm' ([length]) incompatible with 'T''s declared units 'K' ([temperature])
```

Unit strings accept the common notations (`m2/s`, `m²/s`, `m^2/s`,
`kg/(m·s)`, `Pa*s`).

`characteristic_scales` takes `length`, `velocity`, `time`, `density`,
`viscosity` and `temperature`, each either a value with units or `derived:
true` with a formula over `L_char`, `U_char`, `rho_char`, `mu_char` and
`T_char`. `dimensionless_numbers` are formulas over `L`, `U`, `rho`, `mu`,
`T` and `g`:

```yaml
characteristic_scales:
  length: {value: 1.0, units: m}
  velocity: {value: 1.0, units: m/s}
  density: {value: 1.0, units: kg/m3}
  viscosity: {value: 0.025, units: Pa*s}
  time: {derived: true, formula: "L_char / U_char"}
dimensionless_numbers:
  Re: {formula: "rho * U * L / mu"}
```

`ymf.non_dimensionalize(doc)` resolves these (`Re = 40.0` here) and derives
the substitution relations for unknowns whose units it recognizes
(`v* = v / U_char`, `p* = p / (rho_char * U_char**2)`). It does not rewrite
the PDE itself. That needs the symbolic layer, which doesn't exist yet.

## Validating a document

```python
from ymf import load_ymf, check_units, non_dimensionalize

doc = load_ymf("examples/kovasznay_flow.yaml")   # raises on a schema violation
result = check_units(doc)                        # truthy if no errors
for issue in result:
    print(issue)
print(non_dimensionalize(doc))
```

Or run [`examples/validate_spec.py`](../examples/validate_spec.py), which
does all of this for every example and prints the branch tree. These need
the front-end extras: `pip install 'ymf[frontend]'`.

`load_ymf` uses strictyaml, which keeps comments and checks types strictly.
YMF re-enables YAML flow style (`[u, p]`, `{family: CG, order: 2}`), which
strictyaml rejects by default.

## What is checked and what is free text

The schema checks structure, enumerations and types. The mathematics is
free text: `strong_form_expression`, `bilinear`, `linear`, `formula` and
`derivation` are strings that people and LLMs read, and nothing parses them
yet. A sympy-based layer that would turn them into something checkable
(and into solver input) is the main piece of planned work. Until then, a
discretization's correctness rests on review, and on the verification
results recorded in `vvuq`.

`vvuq`, `archive` and `Xdmf` are accepted with any content. The intended
shape of `vvuq` is in [ymf-schema.md §3](ymf-schema.md#3-vvuq-results).
`archive` is reserved for output configuration. The archive format itself
is implemented separately; see [archive-format.md](archive-format.md).

## Known gaps

These are real limitations of the current schema and code, recorded here so
nobody has to rediscover them:

- **Branch references aren't resolved.** A `from_weak_form` or `branch_ref`
  that names nothing still validates.
- **Coefficients come back as strings.** `coefficients` is untyped in the
  schema, so after `load_ymf` a value like `kappa: 1.0e-3` is the string
  `'1.0e-3'`. `compute_dimensionless_numbers` only uses numeric
  coefficients, so a formula that refers to a coefficient (a Péclet or
  Fourier number) is silently dropped. Use characteristic scales for now.
- **Boundary conditions are constants.** `boundary_conditions[].value` is a
  number. Spatially varying data, such as Kovasznay's exact velocity on the
  boundary, can only be stated in a weak form's free-text
  `boundary_conditions`.
- **Write `1.0e-10`, not `1e-10`.** strictyaml accepts both, but PyYAML
  (YAML 1.1), which reads a spec carried inside an archive, treats `1e-10`
  as a string.
- **`characteristic_scales` is a fixed set** of six. There is no slot for a
  diffusivity or a pressure scale.
