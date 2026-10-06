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
- [Composing documents: models and problems](#composing-documents-models-and-problems)
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

extends: model.yaml          # optional: the document this one builds on
kind: model                  # optional: model | problem (default)
composition: ...             # written by the loader into a composed document
```

The last three describe the file rather than the problem; see
[Composing documents](#composing-documents-models-and-problems).

The smallest complete example is
[`examples/poisson.yaml`](../examples/poisson.yaml). The key blocks:

**`physical_model`**: `provenance` (required), plus lists of `processes`
and `assumptions`, and optionally a `source_document`.

**`strong_form`**: `provenance`, `unknowns`, `equation_formulation`,
`strong_form_expression`, `domain` and `boundary_regions` are required.
Unknowns are either bare names (`[u, p]`) or structured entries with
`units` and optional `std_name` and `provenance`. `ymf.normalize_unknowns()`
turns either form into the structured one. Optional: `initial_conditions`,
`boundary_conditions`, `coefficients`, `known_analytical_solution`,
`dimensional_check`.

```yaml
initial_conditions:
  - field: T
    type: function          # constant | function | data_source
    formula: |-
      300 + 10 sin(πx) sin(πy)
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

## Composing documents: models and problems

Branches can also live in separate files. A document names the one it
builds on with `extends`, and `load_ymf` merges the two. The common shape
is a **model** that many **problems** extend:

```text
examples/navier_stokes/
  navier_stokes_model.yaml          kind: model -- physics, equations, unknowns,
  │                                 coefficient units; no domain, conditions or values
  ├── planar_couette.yaml           geometry, BCs, values, exact solution,
  │                                 a Galerkin weak form, Taylor-Hood
  └── plane_poiseuille.yaml         the same, driven by a pressure gradient,
      │                             a stabilized weak form, equal-order P1
      └── plane_poiseuille_re100.yaml   changes the viscosity and gradient only
```

For multiphysics, getting the model right (the processes, the coupled
equations, the unknowns) is often most of the work. Writing it once as a
model means it is reviewed once, and every well-posed problem built on it
shares it rather than copying it.

A model is marked `kind: model`. It is validated on its own against a
**partial** schema, in which any section or field may be missing but
whatever is present must be well formed. A list entry that is present (an
unknown, a weak form) must still be complete. A problem file holds only
what it adds:

```yaml
extends: navier_stokes_model.yaml   # relative to this file
Problem:
  name: "Planar Couette flow, Re = 10"
  physical_model:
    provenance: human_specified
    assumptions: [steady, fully_developed]   # added to the model's
  strong_form:
    domain: |-
      Ω = [0, 4] × [0, 1]  (m)
    coefficients:
      mu: {value: 0.1}                       # the model gave the units
  weak_forms: [...]
solution_paths: {...}
```

**Merge rules.** Mappings merge key by key. Lists of entries with a
`label` or `name` (unknowns, weak forms, discretizations, boundary regions)
merge entry by entry on it, and new entries are appended. Lists of strings
(`processes`, `assumptions`) are unioned. Anything else the child sets
replaces the parent's value. Every file in a chain is checked against the
partial schema as it is read, so an error is reported against the file and
line that holds it. The composed result must then pass the full schema,
unless the file being loaded is itself a model.

**Overrides are recorded, not refused.** When a child changes a value its
parent set, the composed document's `composition` block records it, and
also lists every source file:

```yaml
composition:
  sources: [navier_stokes_model.yaml, plane_poiseuille.yaml, plane_poiseuille_re100.yaml]
  overrides:
    - path: Problem.strong_form.coefficients.mu.value
      was: '1.0'
      now: '0.01'
      set_by: plane_poiseuille_re100.yaml
    - path: Problem.strong_form.coefficients.G.value
      ...
```

(Coefficient values are recorded as strings because the schema leaves
coefficients untyped; see [Known gaps](#known-gaps).)

This makes a parameter variant a three-line file whose changes a reviewer
can read off directly. `extends` and `kind` describe a file, not the
problem, and are dropped from a composed problem. The result is a single
self-contained document, which is what an archive should carry in `extra`.

## Provenance

Every decision block carries `provenance`, one of:

| value | meaning |
|---|---|
| `llm_derived` | proposed by an LLM, not yet confirmed by a person |
| `human_specified` | written by a person |
| `human_edited` | LLM-proposed, then changed by a person |
| `proteus_derived` | derived by a solver (e.g. from a Proteus model) |

It is required on `physical_model`, `strong_form`, every weak form, and
every solution path. Its purpose is to calibrate trust: a reviewer can see
at a glance which parts of a model a person has looked at.

The rule is the same at every level: a block's `provenance` covers
everything inside it, and a nested entry may state its own to override it.
The strong form's `provenance` covers the equations and the choice of
unknowns. A structured unknown that someone else chose says so:

```yaml
strong_form:
  provenance: llm_derived          # the equations, and v
  unknowns:
    - name: v
      units: m/s
    - name: p
      units: Pa
      provenance: human_edited     # a person changed this one
```

`normalize_unknowns(unknowns, provenance=strong_form["provenance"])` fills
in each unknown's effective provenance.

Earlier drafts, including [ymf-schema.md](ymf-schema.md), called the
strong form's key `unknown_provenance`. It was renamed because it read as
"the provenance is unknown", and the old key is now rejected.

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
`T_char`, and over the numeric `coefficients` by name.
`dimensionless_numbers` are formulas over `L`, `U`, `rho`, `mu`, `T`, `g`
and the numeric coefficients:

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

**State each value once.** Viscosity and density are both coefficients of
the equations and characteristic scales. Don't write the number twice:
derive the scale from the coefficient, and a problem that changes the
coefficient changes the scale and every dimensionless number with it.
From [`plane_poiseuille.yaml`](../examples/navier_stokes/plane_poiseuille.yaml):

```yaml
characteristic_scales:
  length: {value: 1.0, units: m}
  velocity: {derived: true, formula: "G * L_char**2 / (8 * mu)", units: m/s}
  density: {derived: true, formula: "rho", units: kg/m3}
  viscosity: {derived: true, formula: "mu", units: Pa*s}
strong_form:
  coefficients:
    rho: {value: 1.0}
    mu: {value: 1.0}
    G: {value: 8.0, units: Pa/m}
```

Its Re = 100 variant then changes `mu` and `G` and nothing else. The unit
check enforces this for the conventional names: a `density` or `viscosity`
scale given as a number beside a `rho` or `mu` coefficient is an **error**
if the two differ and a **warning** if they agree:

```text
[error] characteristic_scales[viscosity]: value 1 differs from coefficients[mu] = 0.01; state it once with {derived: true, formula: "mu"}
```

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
  `'1.0e-3'`. `compute_dimensionless_numbers` coerces numeric strings, both
  bare and as `{value: ..., units: ...}`, so a Péclet or Fourier formula
  that refers to `kappa` now evaluates. Every other consumer of
  `coefficients` still sees strings, and a non-numeric coefficient (a
  formula such as `"π²"`) is still left out of the formula scope, which
  silently drops any number that refers to it.
- **Boundary conditions are constants.** `boundary_conditions[].value` is a
  number. Spatially varying data, such as Kovasznay's exact velocity on the
  boundary, can only be stated in a weak form's free-text
  `boundary_conditions`.
- **Write `1.0e-10`, not `1e-10`.** strictyaml accepts both, but PyYAML
  (YAML 1.1), which reads a spec carried inside an archive, treats `1e-10`
  as a string.
- **Composition can add and change, but not remove.** A problem cannot
  drop an assumption or a field its model states. That takes a new model.
- **An error in a composed problem that no single file holds** (say, a
  required key that none of the files sets) is reported against a line of
  the merged document, labelled with the chain of files, not against a
  source file.
- **Duplicated values are only caught by name.** The check above covers a
  `density` or `viscosity` scale against a `rho` or `mu` coefficient. A
  coefficient under another name (`nu`, `eta`), or a velocity scale that
  restates a boundary speed, isn't checked. Derive scales from coefficients
  wherever possible.
- **`characteristic_scales` is a fixed set** of six. There is no slot for a
  diffusivity or a pressure scale.
