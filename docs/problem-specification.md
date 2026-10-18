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
- [Writing it once: the notation](#writing-it-once-the-notation)
- [Running a spec](#running-a-spec)
- [Branching: one problem, several solution paths](#branching-one-problem-several-solution-paths)
- [Document state](#document-state)
- [Composing documents](#composing-documents)
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
model:                     # required
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

extends: model.ymf          # optional: the document this one builds on
ymf_state: abstract          # optional: abstract | well-posed | discretized | closed
composition: ...             # written by the loader into a composed document
```

The last three describe the file rather than what it specifies; see
[Document state](#document-state) and
[Composing documents](#composing-documents).

The smallest complete example is
[`examples/poisson.ymf`](../examples/poisson.ymf). The key blocks:

**`physical_model`**: `provenance` (required), plus lists of `processes`
and `assumptions`, and optionally a `source_document`.

**`strong_form`**: `provenance`, `unknowns`, `equation_formulation`,
`domain` and `boundary_regions` are required, and `equations` holds the
equations themselves, in the [notation](#writing-it-once-the-notation).
Unknowns are either bare names (`[u, p]`) or structured entries with
`units`, an optional `std_name`, `provenance`, and `rank` (1 for a vector).
`ymf.normalize_unknowns()` turns either form into the structured one.
Optional: `initial_conditions`, `boundary_conditions`, `coefficients`,
`dimensional_check`. (`strong_form_expression` and
`known_analytical_solution` are accepted from older documents; they
duplicate `equations` and the analytical solution path, so new documents
leave them out.)

```yaml
initial_conditions:
  - field: T
    type: function          # constant | function | data_source
    formula: "300 + 10 sin(πx) sin(πy)"
    units: K
boundary_conditions:
  - region: walls           # a boundary_regions[].name
    variable: T
    type: dirichlet         # dirichlet | neumann | robin | periodic
    value: 300.0
    units: K
```

**`weak_forms`**: each has a `label`, `provenance`, a `derivation`, the
trial and test `solution_spaces`, and the `bilinear` and `linear` forms.
Optional notes record inf-sup stability, existence, uniqueness and
regularity, plus the stabilization method.

**`discretizations`**: each has a `name`, `provenance`, `solver`, and
usually `finite_element` and `from_weak_form` (a weak form's label, which
may be omitted when the document has only one — see
[Branching](#branching-one-problem-several-solution-paths)). `finite_element` takes
`fields`, or `velocity` and `pressure` for mixed problems, each as
`{family: CG|DG|RT|BDM|RTF|N1curl, order: n}`. The `solver` block takes
`type` (`linear | nonlinear | time_marching | saddle_point`),
`linear_solver` (`petsc | mumps | umfpack | superlu`), `petsc_options`,
`nonlinear_solver` (`newton | broyden | line_search`), `tolerance` and
`max_iterations`. `petsc_options` are PETSc options without the leading
dash, as strings (`{ksp_type: cg, pc_type: gamg, ksp_rtol: "1e-10"}`; an
option that takes no value is given `""`), for a `linear_solver` of
`petsc` or `mumps`. Proteus honours `linear_solver`: `superlu` is its own
serial LU, `mumps` a direct solve through PETSc, and `petsc` a PETSc
Krylov solve set by `petsc_options` (PETSc's defaults without them);
`umfpack` is refused, as it has no Proteus solver. Omitted, Proteus picks
its serial LU, or MUMPS when the operator has a null space.
`mesh: {cells, levels}` gives the meshes: cells along each side of the
domain's box, and levels of uniform refinement (the solution is on the
finest). A region that is not a box (`kind: piecewise_linear`) gives
`mesh: {h, levels}` instead: `h` is the element size of an unstructured
mesh, the h of convergence rates, and the solver turns it into its
mesh generator's own constraint. A mesh gives `cells` or `h`, never both.
Any of them may be a list, a study with one output per combination.
`discretization` is a free-form block for anything else, such as the time
integrator in `heat_equation.ymf`; its `dt` may be a list too.

## Naming

**Keys are `lower_case_with_underscores`. Values that name a type are
capitalized.**

```yaml
model:
  strong_form:
    unknowns: ["u"]
    coefficients:
      κ: {value: 1.0e-3, units: m2/s}
solution_paths:
  discretizations:
    - name: "p1"
      finite_element:
        fields: {family: CG, order: 1}     # CG is a type, not a key
```

So `strong_form`, `boundary_regions`, `characteristic_scales`, `data_item`,
`number_of_elements` — and `CG`, `DG`, `RT`, `SI`, `Triangle`,
`Tetrahedron`, `Scalar`, `Node`, `Float`, `XYZ`, which name kinds of thing
rather than being keys. A value that is not naming a type is the author's
to spell as they like: a branch called `taylor_hood` and one called
`TaylorHood` are both fine, and neither means anything to the schema.

This holds for the archive as well as the specification. The archive's
vocabulary came from XDMF, which capitalizes everything — elements
(`Topology`, `DataItem`) and their attributes (`NumberOfElements`,
`Precision`) alike — and YMF followed by analogy, including for keys XDMF
does not have at all (`TimeCollections`, `SpatialCollection`, `Attributes`
as a list). Translation to XDMF is a mapping either way, so the archive
uses the same rule as everything else and
[`ymf2xmf`](../src/ymf/cli.py) supplies the capitalization XDMF wants. What
comes back from a round trip is the lower-case form.

### Why the rule, and not case-insensitivity

Composition pairs keys by exact match. A model file writing `Model:` and a
variant writing `model:` would not merge — they would become siblings, and
the result would be a document with two problems, failing the schema
somewhere unhelpful. The same exact-string assumption runs through the
identity keys that pair list entries, `remove:` paths, the override log and
the reference graph.

So where a spelling has to be accepted for compatibility it is folded to the
canonical one as the file is read, before anything is merged. Aliases ease a
rename; they are not a standing second spelling.

## Writing it once: the notation

Every fact in a spec is written once, in the notation a person writes and
reads, and the same string is what `ymf.symbolic` parses. There is no
separate machine-readable copy to keep in step. From the examples:

```yaml
strong_form:
  equations:
    - "∂T/∂t = ∇·(κ∇T)  in Ω"
    - "ρ ∂v/∂t + ∇·(ρ v⊗v) − μΔv + ∇p = f  in Ω"
  domain: "Ω = [0, 1] × [0, 1],  t ∈ (0, 50]"
  boundary_regions:
    - {name: walls, where: "y = 0 or y = H"}
    - {name: pressure_datum, where: "x = 0 and y = 0"}
    - {name: all, where: "∂Ω"}
  boundary_conditions:
    - {region: walls, variable: v, type: dirichlet, formula: "(U, 0)"}
  coefficients:
    κ: {value: 1.0e-3, units: m2/s}
    f: "π² sin(πx) sin(πy)"
    λ: "ρ/(2μ) - sqrt(ρ²/(4μ²) + 4π²)"
solution_paths:
  analytical:
    - solution:
        formula: |-
          v = (1 - exp(λx) cos(2πy), λ/(2π) exp(λx) sin(2πy))
          p = (1 - exp(2λx))/2
```

| written | means |
|---|---|
| `∇u`, `∇·F`, `Δu` | gradient, divergence, Laplacian of the operand that follows |
| `∂u/∂t`, `∂u/∂x` | partial derivatives |
| `a·b`, `a⊗b` | dot product (a product of scalars), outer product |
| `(a, b)` | a vector |
| `x² π³ x^n x**n` | powers |
| `2π²κt`, `κ ΔT`, `sin(πx)` | implicit multiplication |
| `sin cos tan exp log sqrt sinh cosh tanh abs` | functions, with parentheses |
| `grad div lap dt dx dy dz dot outer` | the operators spelled out, for ASCII-only text |
| `lhs = rhs  in Ω` | an equation; the `in ...` clause is ignored |
| `y = 0 or x = 0 and y = H` | a boundary region; `∂Ω` is the whole boundary |
| `Ω = [a, b] × [c, d],  t ∈ (t0, t1]` | the domain, a box, and the time interval of a transient problem |

Names are the coordinates `x y z`, time `t`, `π`, the unknowns and the
coefficients. A coefficient may be any name, Greek or not (`κ`, `ρ`, `μ`,
`G`), and a coefficient's value may itself be a formula in the others and
in `x`, `y`, `z`, `t`. A run of letters that is not a known name is read as
a product of single-letter names (`πx` is π·x) if every letter is known,
and is an error otherwise, so a misspelling is reported, not multiplied:

```text
unknown name 'kapa' (not a coefficient, an unknown, a coordinate, t, π or a
function, nor a product of single-letter names) at position 2 in '2 kapa T'
```

Tensors follow the continuum-mechanics convention `grad(v)[i, j] =
∂v_j/∂x_i`, so `v·∇v` is (v·∇)v and `∇·(μ∇v)` is the vector Laplacian.

The equations are sorted into the advection-diffusion-reaction form that
transport codes solve (mass, advective flux, diffusion, reaction,
Hamiltonian), so they must be in a form that has one: a second derivative
appears as the divergence of a flux, and a nonlinear first-order term
depends on the gradients of the unknowns, not the unknowns themselves.
That is why the Navier-Stokes model writes its advection as `∇·(ρ v⊗v)`,
which equals ρ(v·∇)v when ∇·v = 0. The convective form is refused, with
the reason.

## Running a spec

`ymf.symbolic.adr_problem(doc)` turns a validated, possibly composed, spec
into a plain-data problem: each equation's coefficients and their
derivatives as numpy code strings, the domain, the Dirichlet and periodic
boundary data, initial conditions and the exact solution. A solver consumes
that without sympy. For Proteus, `scripts/ymf_run` does the whole run:

```bash
ymf_run examples/poisson.ymf --outdir out              # every discretization, every mesh
ymf_run out/poisson.archive.ymf --check                # rerun it all, and compare
```

The meshes are part of the spec. Each discretization gives them as a study,

```yaml
  discretizations:
    - name: "P1_linear"
      from_weak_form: "global"
      finite_element: {fields: {family: CG, order: 1}}
      mesh: {cells: [4, 8, 16]}      # cells per side; levels: refines each uniformly
      solver: {type: linear, tolerance: 1.0e-12}
```

and `ymf_run` runs each discretization on each mesh, printing the L2
error of every unknown against the exact solution, with observed rates.
The solver's `tolerance` and `max_iterations` reach Newton. A
discretization whose weak form asks for a stabilization the runner cannot
provide is skipped and says so, rather than run without it. `--cells`,
`--levels` and `--discretization` narrow or change the study from the
command line, and the archive records them as overrides.

Everything goes into one archive per input, `poisson.archive.ymf`: the
spec as realized, then one output per discretization and mesh, keyed by a
hash of the part of the spec that produced it. Each output records the
transformations that produced it (strong form → ADR form → Proteus's
discrete problem → the solve, with software versions), the errors and the
approximation. Like make, `ymf_run` computes only the outputs the archive
lacks; `--check` reruns the rest and confirms that each reproduces its
record bitwise. See [The archive of a spec](archive-format.md#the-archive-of-a-spec).

An unknown that appears in the equations only through its gradient, and
that no Dirichlet condition fixes anywhere, is determined only up to a
constant: the pressure of an incompressible flow with the velocity given
on the whole boundary, say. `adr_problem` lists such unknowns under
`up_to_constant`. Nothing in the spec states it separately, so it cannot
disagree with the equations and boundary conditions. Proteus then keeps the
constant in the null space of its linear operator (PETSc GMRES around a
MUMPS LU that detects the null pivot), rather than pinning a value at a
point, and the errors in those unknowns are measured modulo constants.
[`kovasznay_flow.ymf`](../examples/kovasznay_flow.ymf) is solved this way.
The runs pin BLAS and OpenMP to one thread, since threaded BLAS inside
MUMPS sums in a different order on each run and breaks bitwise
reproducibility; the archive records the thread settings.

How a gradient term is sorted is a choice too: by default `∇p` becomes the
flux `∇·(pI)`, which the weak form integrates by parts;
`adr_problem(doc, hamiltonian_gradients=("p",))` makes it a Hamiltonian,
`H = ∂p/∂x_j`, which it does not. The two formulations differ where the
boundary flux of momentum matters (outflow), not with Dirichlet velocity.

`--emit-pn` also writes each output's problem as the classic Proteus files
`<stem>_<digest>_p.py`, `_n.py` and `_so.py`: what `proteus.ADRProblem`
builds in memory, as readable Python, with the coefficients' code strings
copied verbatim. `parun <stem>_<digest>_so.py` runs them, and its solution
is bitwise the archived one. They are a way to read exactly what was
solved, or to start a hand-written Proteus model from it.

## The domain

`model.strong_form.domain` is the region the model is posed on. It is
either the notation short form, which means a box:

```yaml
domain: "Ω = [0, 1] × [0, 1],  t ∈ (0, 50]"
```

or a mapping that says what **kind** of domain it is:

| kind | holds | for |
|---|---|---|
| `box` | `lower`, `upper`, optional `time` | the same thing the string says |
| `piecewise_linear` | `vertices`, `segments` (2-D) or `facets`/`volumes` (3-D), `holes`, `regions` | a boundary representation — Proteus's PSLG and PLC |
| `discrete` | `topology`, `geometry` | already a mesh; the one kind XDMF holds directly |
| `file` | `format`, `path` | someone else's language — a `.geo`, an STL, a CAD model |

```yaml
domain:
  kind: piecewise_linear
  dimension: 2
  vertices: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
  segments: [[0, 1], [1, 2], [2, 3], [3, 0]]
  boundary_tags: {bottom: 1, right: 2, top: 3, left: 4}
```

`kind` is written, not inferred: a reader should know what it is holding
before it reads it, and whoever writes a domain should have to say what
they are making. A kind is checked against what it requires, and against
what belongs to *another* kind — a `kind: box` carrying `vertices` is a
confusion worth catching, and the message says which kind that key belongs
to.

Many domains are technically expressible in YAML and should still be
`kind: file`. Inlining a CAD model helps nobody, and half of Proteus's own
domain classes are already a path and a format.

**What a domain is not.** It is not how to mesh the region — that is
`discretizations[].mesh` — and it is not the state of whatever object
holds it. Proteus's `Domain` classes carry both, and neither is part of
what a domain *is*.

**`domain` in an archive** is this same domain, discretized: XDMF's
`<Domain>` holding grids of `topology` and `geometry`. The same region at
two levels, as with strong and weak forms.

`kind: box` and `kind: piecewise_linear` can be solved today: a box is
meshed by `cells` per side, any other region by `h`, the element size of an
unstructured mesh (see the discretization's `mesh`). Boundary regions are
still named by `where` predicates, which hold on any mesh of the region;
periodic conditions need a box. `discrete` and `file` domains validate and
are carried, and `adr_problem` says so rather than failing obscurely.

## Branching: one problem, several solution paths

The schema exists to let one document hold *competing* approaches to the
same problem, and to keep track of how they relate.
[`examples/kovasznay_flow.ymf`](../examples/kovasznay_flow.ymf) has this
tree:

```text
strong form: incompressible Navier-Stokes, primitive variables
  weak form 'mixed' ─────────────────── discretization 'taylor_hood'     velocity CG2, pressure CG1
  weak form 'stabilized_equal_order' ── discretization 'equal_order_p1'  velocity CG1, pressure CG1
  analytical 'kovasznay_exact' ◄─────── vvuq.verification.reference
```

### The edges

Those relations are cross-references, and there are four of them:

| from | names |
|---|---|
| `solution_paths.discretizations[].from_weak_form` | a `weak_forms[].label` |
| `vvuq.verification.reference.analytical` | an `analytical[].name` |
| `vvuq.verification.branches[]` | a `discretizations[].name` |
| `vvuq.verification.expected_rates` keys | a `discretizations[].name` |

Each is checked after composition. A reference that names nothing is an
error, with the available names listed — before, it simply validated and
was ignored, so a typo in a verification plan did nothing and said nothing.

A weak form is named by its **label**. Integer indices into `weak_forms`
used to be accepted too, and are not any more: `remove:` can delete a weak
form, which would silently renumber every index after it and repoint
discretizations at a form they do not implement.

An analytical entry has no `from_weak_form`. It solves the strong form; the
weak formulation belongs to theory and numerics. (Non-classical solutions
constructed directly from the integral form of a conservation law are a
real thing, but they are not what that key meant.)

### What may be left unsaid

Two of the edges are single-valued, and where the graph does not branch
they may be omitted. A discretization in a document with one weak form can
only implement that one; a verification plan in a document with one
analytical solution can only compare against that one. So
[`examples/poisson.ymf`](../examples/poisson.ymf) states neither, and
[`examples/kovasznay_flow.ymf`](../examples/kovasznay_flow.ymf), which
branches, states both.

An omitted reference is filled in and **written into the composed
document**, and recorded:

```yaml
composition:
  defaults:
    - path: solution_paths.discretizations[P1_backward_euler].from_weak_form
      value: galerkin
      because: the document has one weak form
      defaulted_by: heat_equation.ymf
```

so the document that comes out holds no implicit bindings, and neither does
an archive carrying it.

Filling in happens as each file is composed, not once at the end. A model
with one weak form is unambiguous when it is written; a problem extending
it may add a second. Resolving against the finished document would make the
model's omission retroactively ambiguous — punishing a file for what a
later file did.

Adding a branch adds entries; it never edits existing ones. That is what
makes a document safe for several contributors, or an LLM and a person, to
extend.

## Document state

A YMF document is the record of a progression: from a model, through a
well-posed specification and its discretizations, to the closure that holds
what was computed and exactly what it was computed from. `ymf_state` says
where along that progression a document stands:

| `ymf_state` | the document | checked by |
|---|---|---|
| `abstract` | part of a model; it must be extended before it can be solved (equations and unknowns, say, with no domain, conditions or values) | the partial schema only |
| `well-posed` | a complete `model:` block | the full schema |
| `discretized` | well-posed, and every discretization gives a mesh (`cells` or `h`), so every output it plans has a key | its planned outputs |
| `closed` | discretized, and every planned output is present and still keyed to the input that produced it | its outputs' keys |

Each state includes the ones before it. A document may declare its state,
and the loader checks that its content reaches it; declaring less is
allowed. Tools that advance a document write the state they reach: a
runner writing an archive records `closed` when every output it plans is
there, and `discretized` while some are still owed. An archive whose input
was edited no longer reaches the `closed` it declares, and loading it says
so, naming the outputs that are now stale. (A runner bringing it up to
date loads it anyway, and recomputes.) A document that declares nothing
must be well-posed, as before.

## Generated documents

A tool that advances a document -- `pro weak` adding weak forms, `pro disc`
adding discretizations -- writes a new file that `extends` its input, so
each step can be reviewed as one file, and records that it did:

```yaml
extends: poisson.ymf
generated:
  by: pro weak (proteus 2.0.0.dev0, ymf 0.2.1.dev)
  reads: model_without_weak_forms   # or model_and_analytical
  from: 9c41e0a2b7d3                # what it read, hashed
  content: 51f0c2aa19e8             # what it wrote, hashed
model:
  weak_forms: [...]
```

`ymf.generated.check` answers two questions from that record. **Stale**:
the part of the input the step read (`reads`) has changed since, so the
step should run again; an edit to a part it did not read does not count.
**Edited**: someone changed what the step wrote. An edited file is
legitimate -- it is the editor's now -- but a tool must not overwrite it.
`ymf.generated.write` writes such a file. Like `extends` and `ymf_state`,
`generated` describes a file: it is not inherited, and an archive of a
generated file does not carry it.

## Composing documents

Branches can also live in separate files. A document names the one it
builds on with `extends`, and `load_ymf` merges the two. The common shape
is an **abstract** document holding the model, which many well-posed
specs extend:

```text
examples/navier_stokes/
  navier_stokes_model.ymf          ymf_state: abstract -- physics, equations, unknowns,
  │                                 coefficient units; no domain, conditions or values
  ├── planar_couette.ymf           geometry, BCs, values, exact solution,
  │                                 a Galerkin weak form, Taylor-Hood
  └── plane_poiseuille.ymf         the same, driven by a pressure gradient,
      │                             a stabilized weak form, equal-order P1
      └── plane_poiseuille_re100.ymf   changes the viscosity and gradient only
```

For multiphysics, getting the model right (the processes, the coupled
equations, the unknowns) is often most of the work. Writing it once as a
model means it is reviewed once, and every well-posed spec built on it
shares it rather than copying it.

An abstract document is marked `ymf_state: abstract`. It is validated on
its own against a **partial** schema, in which any section or field may be
missing but whatever is present must be well formed. A list entry that is
present (an unknown, a weak form) must still be complete. A file that
extends another holds only what it adds:

```yaml
extends: navier_stokes_model.ymf   # relative to this file
model:
  name: "Planar Couette flow, Re = 10"
  physical_model:
    provenance: human_specified
    assumptions: [steady, fully_developed]   # added to the model's
  strong_form:
    domain: "Ω = [0, 4] × [0, H]"
    coefficients:
      μ: {value: 0.1}                        # the model gave the units
      H: {value: 1.0, units: m}              # and this spec adds one
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
unless the file being loaded declares `ymf_state: abstract`.

To fill in a coefficient the abstract document declares, use the mapping
form (`f: {formula: "(G, 0)"}`), which merges with its `{units: ...}`. A
bare string (`f: "(G, 0)"`) replaces the whole entry, units and all, and is
recorded as an override.

**Overrides are recorded, not refused.** When a child changes a value its
parent set, the composed document's `composition` block records it, and
also lists every source file:

```yaml
composition:
  sources: [navier_stokes_model.ymf, plane_poiseuille.ymf, plane_poiseuille_re100.ymf]
  overrides:
    - path: model.strong_form.coefficients.μ.value
      was: 1.0
      now: 0.01
      set_by: plane_poiseuille_re100.ymf
    - path: model.strong_form.coefficients.G.value
      ...
```

**A file may also remove.** `remove:` lists paths this file drops from what
it extends, applied after the merge, so a file can drop an inherited item
and add its replacement in one step:

```yaml
extends: navier_stokes_model.ymf
remove:
  - "model.physical_model.assumptions[incompressible]"
model:
  physical_model:
    assumptions: [compressible]
```

A path is the same shape the override log prints: `.` for a mapping key,
`[...]` for an entry in a list — by its own text in a list of strings
(`assumptions[incompressible]`), by its `label` or `name` in a list of
mappings (`weak_forms[stabilized]`). So a path read out of `overrides` can
be pasted into `remove:`.

This exists because the moves are real ones: incompressible flow solved as
compressible, inviscid Euler as viscous, or in the other direction an
unknown field dropped or replaced by an analytical expression. Before, each
took a new abstract file.

Removals are recorded beside the overrides, with what they took:

```yaml
composition:
  removals:
    - path: model.strong_form.coefficients.μ
      was: {value: 1.0, units: Pa*s}
      removed_by: euler_variant.ymf
```

**A path that names nothing is an error**, not a quiet no-op. It means the
document being composed has moved under the file doing the removing, and
composing the unremoved thing silently is how a variant ends up stating the
opposite of what its author wrote.

This makes a parameter variant a three-line file whose changes a reviewer
can read off directly. `extends` describes a file, not what it specifies,
and is dropped from a composed document; `ymf_state` is kept as the top
file declared it, and is not inherited from the files it extends. The
result is a single self-contained document, which is what an archive
should carry in `extra`.

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
string parses, that each initial condition's units are dimensionally
compatible with the unknown it initializes, and that each boundary
condition's `units` — the units of its right-hand side `g` — fit its type:

- `dirichlet`, `u = g`: `[g] = [u]`.
- `neumann`, `q·n = g` with `q = -a∇u` the diffusive flux of the variable's
  own equation: `[g] = [a][u]/[L]`.
- `robin`, `αu + β q·n = g`: `[α][u] = [β][q·n] = [g]`. The first equality
  needs no units for `u` or `g`, since `u` is a factor of both sides; `α`,
  `β`, `a` and the domain's length decide it.

```text
[error] boundary_conditions[region='walls', variable='T']: units 'm' ([length]) incompatible with 'T''s declared units 'K' ([temperature])
[error] boundary_conditions[region='wall', variable='T']: g is K but the flux q.n of T is kg/s³
```

Only dimensions are compared, never scale, so `km` and `m` agree. Reading
`α`, `β` and the flux needs ymf[symbolic]. A condition that cannot be
checked because something it involves declares no units — the unknown, `g`,
a coefficient, the domain — is skipped silently by default.
`check_units(doc, missing_units="warn")` or `"error"` reports it instead,
naming what is undeclared:

```text
[warning] boundary_conditions[region='wall', variable='T']: not checked: nothing declares the units of k
```

Unit strings accept the common notations (`m2/s`, `m²/s`, `m^2/s`,
`kg/(m·s)`, `Pa*s`).

`characteristic_scales` takes `length`, `velocity`, `time`, `density`,
`viscosity` and `temperature`, each either a value with units or `derived:
true` with a formula over `L_char`, `U_char`, `rho_char`, `mu_char` and
`T_char`, and over the numeric `coefficients` by name.
`dimensionless_numbers` are formulas over `L`, `U`, `rho`, `mu`, `T`, `g`
and the numeric coefficients. Both may be written in the
[notation](#writing-it-once-the-notation) (`"1/(2π²κ)"`, `"ρ U L / μ"`),
which needs `ymf[symbolic]` to evaluate:

```yaml
characteristic_scales:
  length: {value: 1.0, units: m}
  velocity: {value: 1.0, units: m/s}
  density: {derived: true, formula: "ρ", units: kg/m3}
  viscosity: {derived: true, formula: "μ", units: Pa*s}
  time: {derived: true, formula: "L_char / U_char"}
dimensionless_numbers:
  Re: {formula: "ρ U L / μ"}
```

`ymf.non_dimensionalize(doc)` resolves these (`Re = 40.0` here) and derives
the substitution relations for unknowns whose units it recognizes
(`v* = v / U_char`, `p* = p / (rho_char * U_char**2)`). It does not yet
rewrite the PDE itself in dimensionless form.

**State each value once.** Viscosity and density are both coefficients of
the equations and characteristic scales. Don't write the number twice:
derive the scale from the coefficient, and a problem that changes the
coefficient changes the scale and every dimensionless number with it.
From [`plane_poiseuille.ymf`](../examples/navier_stokes/plane_poiseuille.ymf):

```yaml
characteristic_scales:
  length: {derived: true, formula: "H", units: m}
  velocity: {derived: true, formula: "G H²/(8μ)", units: m/s}
  density: {derived: true, formula: "ρ", units: kg/m3}
  viscosity: {derived: true, formula: "μ", units: Pa*s}
strong_form:
  coefficients:
    ρ: {value: 1.0}
    μ: {value: 1.0}
    G: {value: 8.0, units: Pa/m}
    H: {value: 1.0, units: m}
```

Its Re = 100 variant then changes `μ` and `G` and nothing else. The unit
check enforces this for the conventional names: a `density` or `viscosity`
scale given as a number beside a `ρ`/`rho` or `μ`/`mu` coefficient is an **error**
if the two differ and a **warning** if they agree:

```text
[error] characteristic_scales[viscosity]: value 1 differs from coefficients[mu] = 0.01; state it once with {derived: true, formula: "mu"}
```

## Validating a document

```python
from ymf import load_ymf, check_units, non_dimensionalize

doc = load_ymf("examples/kovasznay_flow.ymf")   # raises on a schema violation
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

The schema checks structure, enumerations and types. With `ymf[symbolic]`,
the mathematics is checked too, when a spec is run: the equations, the
domain, the boundary regions, every coefficient, condition and solution
formula are parsed, and a string that does not parse, or names something
undefined, is an error that says where. Equations must also have an
advection-diffusion-reaction form (see
[the notation](#writing-it-once-the-notation)).

Still free text, for people: `derivation`, the weak forms' `bilinear` and
`linear`, `solution_spaces`, and notes. Weak forms are not parsed yet: the
solver derives the weak form from the strong one (Proteus integrates the
flux terms by parts), and the weak form's `stabilization_method` selects
what the solver adds.

`vvuq`, `archive` and `Xdmf` are accepted with any content. The intended
shape of `vvuq` is in [ymf-schema.md §3](ymf-schema.md#3-vvuq-results).
`archive` is reserved for output configuration. The archive format itself
is implemented separately; see [archive-format.md](archive-format.md).

## Known gaps

These are real limitations of the current schema and code, recorded here so
nobody has to rediscover them:

- *(Closed.)* Branch references were not resolved: a `from_weak_form` or a
  `branch_ref` that named nothing still validated. All four edges of the
  reference graph are now checked, `branch_ref` has become
  `reference.analytical`, integer indices into `weak_forms` are gone, and a
  single-valued reference may be omitted where the graph does not branch.
  See [Branching](#branching-one-problem-several-solution-paths).
- *(Closed.)* `coefficients` was untyped, so every scalar arrived as a
  string, and a coefficient written as an expression
  (`λ: "ρ/(2μ) - sqrt(ρ²/(4μ²) + 4π²)"`, `Q: "π²"`) was dropped from
  formula scope, silently taking any number that referred to it. Numbers
  are now typed in the schema, and expressions are evaluated by the one
  resolver `adr_problem` uses, so `ymf.units` and the solver give the same
  answer for the same coefficient. `sympy` moved into the `frontend` extra
  for this: a front-end that cannot evaluate an expression would report
  different numbers from the solver reading the same file. A bare
  `ymf[units]` install still works, and still sees only the coefficients
  that are already numbers.
- *(Closed.)* A boundary condition's units were compared with its
  unknown's whatever its type, which is right only for `dirichlet`: a
  `neumann` flux in W/m² failed and one in K passed. Each type now has its
  own rule (see [Units](#units-scales-and-dimensionless-numbers)). What is
  left: `neumann` and `robin` see the diffusive flux only, and `t` has no
  declared units anywhere, so an `α` or `β` that depends on time is never
  decided.
- **Dirichlet, Neumann, Robin and periodic conditions run**, with values
  given as a number, a formula or a table, and only on a box.
  A `neumann` condition gives the **outward** diffusive flux `q·n` with
  `q = -a∇u`, which is the sign Proteus's
  `diffusiveFluxBoundaryConditions` wants; where no condition claims a
  face, the natural condition applies, which for a diffusion term is zero
  flux. A Neumann-only problem is detected as determined up to a constant,
  like any other with no Dirichlet condition.
- **A value may be given as a table**, `series: [[t, value], ...]` with `t`
  ascending — a measured or design record, which is the one thing an
  expression cannot be. Linear between samples; outside the range the first
  and last values are held, because extrapolating a record past its end is
  not something to do quietly.

  This is a **value form**, not a kind of condition: it is orthogonal to
  `type`, exactly as `value` and `formula` are. Under `neumann` it is a
  discharge hydrograph; under `dirichlet` a stage hydrograph — a gauge
  reading against time — or a measured inlet temperature. `hydrograph` used
  to be a `type`, which tangled the two and made the Dirichlet case
  impossible to write. A series is one number per time, so it is refused on
  a vector unknown: give one condition per component.
- **`robin` runs too**, `αu + β q·n = g`, with `alpha` and `beta` beside
  the condition and `g` given as `value` or `formula` like any other
  right-hand side. All three are expressions, so they may vary along the
  boundary. It is one relation, so a vector unknown wants one condition per
  component. `beta: 0` is refused — that is a Dirichlet condition, and
  belongs where it can be imposed strongly.

  This needed a numerical flux written for Proteus
  (`NumericalFlux.RobinFactory`), not wiring: the flux depends on the
  solution's trace, and a Proteus flux boundary condition is a function of
  `(x, t)` alone. That the *dependency set* of a boundary expression is the
  thing to get right, rather than the list of condition types, is the point
  taken up in `YMF_TAG_PROPAGATION.md`.
- **Only `box` and `piecewise_linear` domains are solved.** A domain of any
  kind validates and is carried — see [The domain](#the-domain) — but
  `adr_problem` refuses `discrete` and `file`, saying which kind it was
  given rather than failing obscurely.
- **The notation's limits**, each of which now reports the fix as well as
  the fault. Functions need parentheses — `sin πx` says to write `sin(x)`,
  not `sin x`. There are no subscripts — `p₀` says to write `p0` and declare
  it, since the superscripts `² ³ ⁴` *are* notation (powers) and a subscript
  is a reasonable guess. And a multi-letter name must be declared to be a
  name rather than a product — `Re` says to declare it, and that it is not
  read as `R·e` because `R` and `e` are not names.
- **SUPG/PSPG with equal-order P1 on Kovasznay flow is pre-asymptotic at
  reachable mesh sizes**, not wrong. Measured rates sit below the a priori
  values on coarse meshes and climb towards them as the mesh refines (1.88
  for velocity at 512 cells). This was investigated and closed: it is a
  resolution limit, not a defect in `NavierStokesASGS_velocity_pressure`.
  Taylor-Hood on the same problem reaches its rates sooner, so a
  side-by-side comparison at one mesh size misleads. Reaching the
  asymptotic range needs an iterative solver (512 cells is already 24 GB
  with MUMPS) or a lower-Reynolds companion case.
- *(Closed.)* `1e-10` used to be a string to PyYAML and a number to
  strictyaml, so a spec meant different things read as a file and read
  inside an archive. The archive reader now resolves YAML 1.2 scientific
  notation, so `1e-10`, `1E-10`, `1e+10`, `-1e-10` and `1.5e3` are all
  numbers to both. The resolver is registered on a subclass, so pyyaml's own
  loaders are unchanged for anything else in the process.
- *(Closed.)* Composition could add and change but not remove, so dropping
  an assumption or a field a model stated took a whole new model file.
  `remove:` now does it, recorded in `composition.removals`. See
  [Composing documents](#composing-documents-models-and-problems).
- *(Closed.)* An error that no single file holds -- a required key none of
  them sets -- was reported against a line of the merged document, which
  exists nowhere on disk. The message now also names the path and the files
  that write it: *no file in the chain sets `equation_formulation` under
  `model.strong_form`; `model.strong_form` is written by
  `model.ymf`, `model.ymf`*. Where one file added an incomplete entry,
  only that file is named. The path is the same shape `remove:` and the
  override log use, so it can be pasted between them.
- *(Mostly closed.)* Duplicated values were caught by name only, so a
  coefficient called `nu` or `eta`, or any scale outside the conventional
  six, was invisible. A scale whose value *and dimension* match a
  coefficient's is now reported whatever either is called. Dimension is what
  separates a duplicate from a coincidence: in a non-dimensionalized problem
  half the quantities are 1.0, and 1.0 kg/m3 beside 1.0 m is not one value
  written twice. Where either side states no units, nothing is said, because
  there is no way to tell. Still unchecked: a velocity scale that restates a
  boundary speed, which would need the boundary condition's expression
  evaluated. Derive scales from coefficients wherever possible.
- *(Closed.)* `characteristic_scales` was a fixed set of six, with no slot
  for a diffusivity or a pressure scale. It now takes any name. The
  original six keep their aliases (`L`, `U`, `rho`, `mu`, `T`) in formula
  scope; any other name resolves under itself.
