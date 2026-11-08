# Writing a YMF problem specification

A problem specification is a YAML document that describes a PDE problem
from its physics down to a solver configuration, so that a person, an LLM
and a solver can all work from the same record. This guide covers the
current schema (v0.2) as `ymf.schema` implements it, using the documents in
[`examples/`](../examples/) to illustrate. For the design history and
rationale, see [history/schema-v0.1.md](history/schema-v0.1.md) and
[history/schema-v0.2-delta.md](history/schema-v0.2-delta.md); for the gaps
that have been closed and why, [changes.md](changes.md).

- [The idea](#the-idea)
- [Document structure](#document-structure)
- [Writing it once: the notation](#writing-it-once-the-notation)
- [Running a spec](#running-a-spec)
- [Branching: one problem, several solution paths](#branching-one-problem-several-solution-paths)
- [Gauges](#gauges)
- [Splitting a discretization into sub-problems](#splitting-a-discretization-into-sub-problems)
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
  mesh_generation: ...       # optional: where a mesh comes from (DEM, filters)

solution_paths:              # required
  analytical: [...]          # exact solutions, possibly empty
  discretizations: [...]     # stages 4-5, possibly empty

vvuq: ...                    # optional, free-form: verification/validation plans and results

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
Optional: `initial_conditions`, `boundary_conditions`, `coefficients`.
(`strong_form_expression`, `known_analytical_solution` and
`dimensional_check` used to be accepted and are refused now; see
[changes.md](changes.md).)

An equation is a string, and the strings are taken in the order of the
unknowns. Or it is named, and says which unknown it is `for`:

```yaml
equations:
  - name: continuity
    for: p
    equation: "∇·v = 0  in Ω"
  - name: momentum
    for: v
    equation: "ρ ∂v/∂t + ∇·(ρ v⊗v) − ∇·(μ∇v) + ∇p = f  in Ω"
```

With `for` on every equation the order written no longer matters (it is
`for` on all of them or on none). A name is unique, is what a file that
extends this one merges onto (replacing `momentum`'s text replaces that
equation and is recorded as an override), and is how a discretization will
say which equations each of its sub-problems solves. An abstract model may
leave an unknown without an equation for the files that extend it to
supply; a well-posed document may not.

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

**Boundary conditions** of all four types run. A `neumann` condition gives
the **outward** diffusive flux `q·n` with `q = -a∇u`, which is the sign
Proteus's `diffusiveFluxBoundaryConditions` wants; where no condition claims
a face, the natural condition applies, which for a diffusion term is zero
flux. A Neumann-only problem is detected as determined up to a constant,
like any other with no Dirichlet condition. A `robin` condition is
`αu + β q·n = g`, with `alpha` and `beta` beside the condition and `g` given
as `value` or `formula` like any other right-hand side. All three are
expressions, so they may vary along the boundary. It is one relation, so a
vector unknown wants one condition per component. `beta: 0` is refused:
that is a Dirichlet condition, and belongs where it can be imposed strongly.
`periodic` conditions need a box, or a region whose two faces are
translations of each other (see [The domain](#the-domain)).

**Aliases, and what is solved.** Dirichlet, Neumann and Robin (and
periodic pairing) are the abstract forms; each is complete with its
expression. Four more kinds are **aliases**, shorthand for particular
combinations of them, defined as data in `ymf.conditions.ALIASES`:

| alias | expands to | value |
|---|---|---|
| `slip` (a vector) | `dirichlet` on the `normal` component, the value (0 if none is given); `neumann` on the `tangential` component, 0 (no tangential traction) | optional |
| `open` | `neumann` 0: the natural condition, declared (C12's open boundary) | none |
| `inflow` | `dirichlet`, the value, `only: inflow` (where the flow enters) | required |
| `characteristic` | `dirichlet`, the value, `only: incoming_characteristics` (a tailwater stage that a supercritical outflow ignores) | required |

A spec may write the aliases. **What is solved, and what the closure of
the spec records, is the expansion**: `ymf.conditions.expand` replaces
every alias, and every applied condition set, with its abstract
conditions, each marked with what it came from:

```yaml
# boundary_sets: [{region: walls, set: free_slip}] closes as
boundary_conditions:
  - {region: walls, variable: v, type: dirichlet, component: normal, value: 0.0,
     expanded_from: {type: slip, set: free_slip}}
  - {region: walls, variable: v, type: neumann, component: tangential, value: 0.0,
     expanded_from: {type: slip, set: free_slip}}
```

The archive of a spec holds the expansion, its output keys are computed
from it (so a spec and its archive key alike, and a spec without aliases
keeps its keys), and `adr_problem`, `split_problem` and `ymf.units` read
only it. Otherwise two solvers could close the same spec differently and
the archive could not show which. An expanded condition may say how it is
`imposed`, `strong` or `weak`: C12 makes that a property of each
condition. One ADR system (`adr_problem`) runs Dirichlet, Neumann and
Robin on whole variables; a condition on a component, restricted with
`only`, or imposed weakly is for a split discretization's solver, and
`adr_problem` says so.

A condition may also name a coefficient of the state rather than an
unknown: `{variable: χ, type: inflow, value: 1.0}` says the fluid entering
is air, whatever carries `χ` in a given discretization. Checked on load:
the region exists, the variable is an unknown or a coefficient, `slip` and
a `component` are on a vector, and a kind gets a value exactly when it
takes one.

**Condition sets.** What a free-slip wall or an open top means is part of
the model: written once there, over its own unknowns and coefficients, and
applied to regions by each problem.

```yaml
# the model
condition_sets:
  free_slip:
    - {variable: v, type: slip}
  atmosphere:
    - {variable: p, type: dirichlet, value: 0.0}
    - {variable: v, type: open}
    - {variable: χ, type: inflow, value: 1.0}     # what enters is air
# a problem extending it
boundary_sets:
  - {region: walls, set: free_slip}
  - {region: top, set: atmosphere}
```

`ymf.boundary_conditions(strong_form)` is every condition with the sets
expanded, each member carrying its region and the set it came from. A
condition written out under `boundary_conditions` for the same region and
variable replaces the set's member, so a case can apply `atmosphere` and
give the pressure its own value; two sets that both decide one region's
variable are an error, as is a set nobody defines. A set's members are
checked where they are defined, as any condition is. Like an alias, a set is
input shorthand: the closure writes its members out (above).

A value is given as a number (`value`), an expression (`formula`), or a
table, `series: [[t, value], ...]` with `t` ascending: a measured or design
record, which is the one thing an expression cannot be. A series is linear
between samples; outside its range the first and last values are held,
because extrapolating a record past its end is not something to do quietly.
It is a **value form**, not a kind of condition, orthogonal to `type`
exactly as `value` and `formula` are. Under `neumann` it is a discharge
hydrograph; under `dirichlet` a stage hydrograph, a gauge reading against
time, or a measured inlet temperature. A series is one number per time, so
it is refused on a vector unknown: give one condition per component.

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
`linear_solver` (`direct | iterative`, or a package: `petsc | mumps | umfpack
| superlu`), `petsc_options`,
`nonlinear_solver` (`newton | broyden | line_search`), `tolerance` and
`max_iterations`. `petsc_options` are PETSc options without the leading
dash, as strings (`{ksp_type: cg, pc_type: gamg, ksp_rtol: "1e-10"}`; an
option that takes no value is given `""`), for a `linear_solver` of
`petsc` or `mumps`. Proteus honours `linear_solver`: `superlu` is its own
serial LU, `mumps` a direct solve through PETSc, and `petsc` a PETSc
Krylov solve set by `petsc_options` (PETSc's defaults without them);
`umfpack` is refused, as it has no Proteus solver. Omitted, Proteus picks
its serial LU, or MUMPS when the operator has a null space.
`direct` and `iterative` say what is asked of the linear solve rather than
which package does it: a direct factorization (serial or distributed, the
solver's choice), or a Krylov method with a preconditioner that
`petsc_options` refine. Every solver path of a solver reads the same two
words and maps them to its own options, and refuses what it cannot honour
instead of overriding the spec. A package name is a request of its kind
(`superlu`, `mumps` and `umfpack` direct; `petsc` iterative, unless
`petsc_options` give `pc_type: lu` or `cholesky` with `ksp_type: preonly`):
`ymf.vocabulary.linear_solver_kind` says which.
`mesh: {cells, levels}` gives the meshes: cells along each side of the
domain's box, and levels of uniform refinement (the solution is on the
finest). A region that is not a box (`kind: piecewise_linear`) gives
`mesh: {h, levels}` instead: `h` is the element size of an unstructured
mesh, the h of convergence rates, and the solver turns it into its
mesh generator's own constraint. A mesh gives `cells` or `h`, never both.
Any of them may be a list, a study with one output per combination.
`discretization` is a free-form block for anything else, such as the time
integrator in `heat_equation.ymf`; its `dt` may be a list too.

A mesh is one way of reducing the model to a computable problem, and
`method` says which way an entry takes:

| `method` | its block | what it needs | studied |
|---|---|---|---|
| `mesh` (default) | `mesh:` | `cells` (a box) or `h` | `cells`, `h`, `levels` |
| `neural` | `neural:` | `width`, `depth`, `points` | all three |
| `particle` | `particle:` | `spacing`, `kernel` | `spacing` |

Every method also studies `dt`. `neural` (PINNs, weak adversarial
networks) and `particle` (RKPM, SPH) are schemas so far: their fields
validate, and a missing one is named with a hint, but no outputs are
planned for them yet, so a spec using one reaches `well-posed` and says
why it goes no further. Whether a solver can run a method is the solver's
question; the document's state does not depend on it.

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
| `hypot(a, b)` | `sqrt(a² + b²)` rounded once (C99 `hypot`, `numpy.hypot` in code): a distance that matches code computing it with `hypot`, bit for bit |
| `heaviside(φ)`, `sign(e)`, `min(a, b, ...)`, `max(a, b, ...)` | sharp functions; a solver needs `heaviside` and `sign` of an unknown regularized, since their derivatives are Dirac deltas |
| `norm(v)`, `sym(∇v)` | `sqrt(v·v)`, and the symmetric part `(A + Aᵀ)/2` of a square tensor |
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

**A coefficient may be a function of the state.** Its formula may name the
unknowns as well: a conductivity that grows with temperature, a density
that follows a phase indicator, a friction law in the depth and discharge.

```yaml
strong_form:
  unknowns: [T]
  equations:
    - "∂T/∂t = ∇·(κ∇T)  in Ω"
  coefficients:
    κ: {formula: "κ0 (1 + T²)", units: m2/s}
    κ0: {value: 1.0e-3, units: m2/s}
```

It is substituted into the equations like any other coefficient, so the ADR
form carries the dependence: `κ` here makes the diffusion of `T` nonlinear
in `T`, with its derivative, exactly as if `κ0 (1 + T²)` were written in the
equation. Write a variable coefficient inside the divergence (`∇·(κ∇T)`,
`∇·(μ∇v)`); `κΔT` is the divergence form only while `κ` is constant, and is
refused otherwise. A formula may name the unknowns, not their derivatives.

Such a coefficient is not a number, and three things that need one say so
rather than go quiet. Boundary data, initial conditions and exact solutions
are evaluated at points in `x`, `y`, `z`, `t`, so one that names it is an
error (the heat example's exact solution names `κ`, and stops being exact
the moment `κ` depends on `T`). A characteristic scale or a dimensionless
number built on it has no value, so `ymf.units` warns and asks for a
constant coefficient instead (`ρw`, not `ρ = ρw (1 − χ) + ρa χ`). And
`coefficient_numbers` leaves it out. A coefficient may not share its name
with an unknown.

`adr_problem` substitutes such a coefficient, and also lists it by name
under `coefficients` -- its code, what it depends on, and its derivative
in each unknown (`sharp: true` and no derivative for `heaviside` or `sign`
of an unknown) -- so a code generator can emit `ρ` once as a function of
the state rather than its expansion in every term. Each also has its
`inline` code -- its formula over the state, constants substituted, other
coefficients of the state kept by name (`uses`) -- and the problem gains
`named`: the `order` to evaluate those coefficients in (each names only
earlier ones) and the ADR form's slot values over the names (`m: v_0*ρ`,
`a: [2*μ, μ]`). The equations are read a second time for it, with each
coefficient of the state a known field, so a gradient of one is kept (as
`grad_ρ_0`) rather than lost as a constant's. The derivatives and
`depends_on` are only in the substituted form, which is exact; in `named` a
name is held fixed. A problem with no such coefficient has neither key and
reads exactly as before.

A closed spec (an archive) names them too: closure writes each coefficient
of the state in the map form with `depends_on`, the unknowns it follows,
directly or through another coefficient (`ρ: {formula: "ρw (1 − χ) + ρa χ",
depends_on: [φ]}`). It is read from the formulas' text (`ymf.names`, no
sympy, so every environment closes a spec alike), checked on load (a
`depends_on` that is not what the formula says is refused), and does not
count towards an output's key: the formula it comes from already does.

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

**Written back.** `ymf.symbolic.to_notation` is the inverse: it writes a
sympy expression in this notation so that parsing it gives back exactly the
same expression. That is how a tool records what it derived -- the weak
form `pro weak` will solve, say -- in the same terms the spec uses. Two
details make the round trip exact. A numeral is an IEEE double, as the
solver computes with, however many digits it has. And the components of a
vector unknown are names too (`v_0`, `v_1`), because what works on them,
the scalarized form a solver reads, has to be writable.

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
unstructured mesh (see the discretization's `mesh`). A periodic condition
on a `piecewise_linear` region holds when the region's boundary pieces on
the two faces of each periodic axis are translations of each other: the
segments (2-D) or facets (3-D) on the lower face, moved by the length of
the bounding box along the axis, must be exactly those on the upper face.
That is checked at posing, with a count of the pieces that have no partner;
the mesher must then add no points to them, and the solver checks after
meshing that every node on one face has its partner on the other.

A `file` domain is carried whole for the solver to read. A mesh file
(`2dm`, `3dm`, `hex`, `tetgen`, `ymf`) fixes the mesh, so its discretization
gives `mesh: {levels: ...}` only, the levels of uniform refinement of the
file's mesh, and a mesh size there is reported as what it lacks; a geometry
(`gmsh_geo`) is meshed by `h` like a region. `bounds: {lower, upper}` is
optional: it sets the scale of the boundary tolerance, and without it the
solver takes the bounds from the mesh, so for this kind alone posing
follows meshing; the dimension must then be given. A file's boundary
regions are naturally its boundary tags. A file domain cannot be periodic
(Gmsh's `Periodic` constraint is not carried).

A `discrete` domain is a mesh as ymf stores one (`topology`, `geometry`, and
optionally tagged `facets` and `provenance`; see `ymf.mesh`): a mesh is a
type of domain. It is posed whole, like a file. Its bounds are its nodes'
when the geometry holds them inline; when they are in an HDF5 file the
solver takes them, or `bounds` gives them, and the dimension comes from the
geometry's type (`XY`, `XYZ`). Its boundary regions are its facets' tags, a
discretization gives `mesh: {levels}` only, and it cannot be periodic. The
design is [design/domains-beyond-box.md](design/domains-beyond-box.md).

**Boundary regions by position or by tag.** A region is named by `where`, a
predicate on the point that holds on any mesh of the domain, or by `tags`,
the mesh's boundary tags: ints, or names the domain's `boundary_tags`
resolve. Exactly one of the two:

```yaml
domain:
  kind: piecewise_linear
  vertices: [[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]]
  segments: [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]
  segment_flags: [1, 2, 3, 3, 4, 5]
  boundary_tags: {bottom: 1, right: 2, notch: 3, top: 4, left: 5}
boundary_regions:
  - {name: dry, tags: [bottom, left]}
  - {name: outlet, where: "x = 2"}
```

A tag region reaches a solver as `flags: [1, 5]`, which it tests against the
boundary flag of the point, and carries no `where`, so a solver that does
not read `flags` fails on it rather than taking a missing predicate for the
whole boundary. A box carries no tags (name its faces by position), and a
periodic region is found by position, so it needs `where`.

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

## Gauges

What a run measures as it goes is part of the problem, beside its
conditions: `model.gauges`, one entry per gauge file.

```yaml
gauges:
  - name: heights                   # the level set along four verticals
    kind: line
    fields: [φ]
    lines: [[[2.724, 0.0], [2.724, 1.8]], [[2.228, 0.0], [2.228, 1.8]]]
    file: height1.csv
  - name: pressure
    kind: point
    fields: [p]
    points: [[3.22, 0.16], [3.22, 0.584]]
  - kind: boundary_integral         # forces and moments per boundary region
    tags: [cylinder]
    quantities: [force, moment]
    moment_point: [0.2, 0.2, 0.0]
    sample_rate: 0.01
```

`point`, `line` and `line_integral` gauges sample `fields` (unknowns or
coefficients, `ρ` as well as `p`) at `points` or along `lines`, each line
its two ends. A `boundary_integral` gauge integrates `quantities`
(`wetted_area`, `force_pressure`, `force_viscous`, `force`, `moment`, the
moment about `moment_point`) over each boundary region in `tags`. Checked
on load: what each kind needs and may not have, that fields and tags
exist, and that a gauge's points share one dimension. A solver reads them;
Proteus's reader is C22's next step.

## Splitting a discretization into sub-problems

A model with several unknowns is often solved as a sequence of smaller
problems in each time step. Two-phase level-set flow is the standard case:
the model is momentum, continuity and the interface transport of a level
set `φ`, with `ρ` and `μ` coefficients of `χ = heaviside(φ)`; Proteus solves
it as flow, then a volume-of-fluid carrier of `χ`, then the level set, then
redistancing, then a mass correction. That splitting is a property of the
discretization, not of the model, so another discretization of the same
model can split it differently (CLSVOF and a projection scheme) or not at
all. [`tests/fixtures/two_phase`](../tests/fixtures/two_phase) holds the
model and the dam break split this way; [`tests/fixtures/shallow_water`](../tests/fixtures/shallow_water)
holds the shallow-water family the same way (one sub-problem, as Proteus's
SW2DCV solves the whole system at once), with a lake at rest and the
Torino flume as cases.

```yaml
discretizations:
  - name: level_set_five
    provenance: proteus_derived
    mesh: {h: 0.0322}
    binds:                          # symbols the model leaves to the discretization
      h: {meaning: element_size}
      ε: {formula: "1.5 h"}
    regularize:                     # sharp model functions made smooth
      heaviside: {by: smoothed_heaviside, width: ε}
    step_control: {type: min_cfl, cfl: 0.25}
    solver: {type: nonlinear, nonlinear_solver: newton, linear_solver: petsc}
    subproblems:                    # solved in this order every step
      - name: flow
        solves: [momentum, continuity]
        lags: {φ: previous_step, χ: {via: vof, at: previous_step}}
        finite_element: {velocity: {family: CG, order: 1}, pressure: {family: CG, order: 1}}
        solver: {max_iterations: 100}           # refines the entry's solver
      - name: vof
        introduces: {name: χh, approximates: χ}
        equation: "∂χh/∂t + ∇·(v χh) = 0"
        lags: {v: {via: flow, at: current_step}}
      - name: ncls
        solves: [interface]
        lags: {v: {via: flow, at: current_step}}
      - name: rdls
        enforces: "norm(∇φ) = 1"
        pseudo_time: steady
        writes_back: {φ: replace}
      - name: mcorr
        enforces: "heaviside(φ) = χh"
        introduces: {name: φc, correction_of: φ}
        writes_back: {φ: add}
```

A sub-problem **solves** model equations by name (so the model's equations
must be [named](#document-structure)), or **enforces** a constraint the
others lose, or **introduces** an unknown of its own (an auxiliary carrier
of a model quantity, or a correction to a model unknown), or spins the run
up. It says what it **lags** -- reads from elsewhere -- and when: as the
previous step left it, or as an earlier sub-problem of this step left it
(`via`). What it **writes back** into another field it says too. Its
`solver` refines the entry's, which stays required as the default.
`binds` names the symbols the discretization supplies and the model may not
use (the element size, a smoothing width), and `regularize` the sharp
functions it replaces with smooth ones (`heaviside`, `sign`), and the
reciprocal of an unknown that may vanish: `inverse: {by: wet_dry_cutoff,
of: h, width: hε}` reads every negative power of the depth `h` (`q⊗q/h`,
Manning's `h^(-7/3)`) through Proteus SW2DCV's `one_over_hReg`,
`2h/(h² + max(h, hε)²)`, which is `1/h` wherever `h ≥ hε` and 0 at `h = 0`.
It is applied before anything is differentiated, so the derivatives are the
regularized expression's; `of` must name a scalar unknown the sub-problem
solves for.

Checked on load: sub-problem names are unique; every model equation is
solved by exactly one sub-problem; every name a sub-problem uses (an
equation, a lagged field, a `via`, what it introduces, approximates,
corrects or writes back) exists, and an introduced or bound name does not
shadow one that does. The split is part of what an output's key covers, so
changing a sub-problem's solver plans a new output.

A sub-problem may give an `equation` of its own, `for` the unknown it
determines (by default the one it introduces); a lag with `of` reads a
copy of a field under a name of its own (redistancing reads `φs`, the
level set before it changes `φ`); and a sub-problem's own `regularize`
replaces the entry's for the functions it names (redistancing smooths its
sign over `0.75 h`, the rest their Heaviside over `1.5 h`).

**A closed spec writes the split out.** What the rules above leave implicit
is written into each sub-problem when a spec is closed (`ymf.explicit`,
for the reason the boundary-condition aliases are expanded: two readers
must not close one spec differently): its `for` (the unknown it
introduces, when it gives an equation and no `for`); its `solver`, the
entry's refined by its own; its `regularize`, the entry's with its own
replacing function by function; and every lag as `{of, via, at,
carried_by}`, where `carried_by` names the introduced unknown a lagged
coefficient is read as (`χ: {via: vof, at: previous_step, carried_by:
χh}`). A `carried_by` that is not the carrier is refused on load. Written
out or not, the split classifies the same (`split_problem` gives equal
results for both) and keys its outputs alike.

**Each sub-problem is classified on its own.** `adr_problem` reads a model
as one system, which this one is not (`v·∇φ` is a Hamiltonian in another
unknown). `ymf.symbolic.split.split_problem(doc, "level_set_five")` gives
each sub-problem's ADR form instead: its unknowns are those of the
equations it solves and of its own equation; every other field it reads is
a known field, declared in its `lags` (a read that is not is an error that
says so), which the plain data names (`v_0`) and whose gradient it names
`grad_v_0_0`; a coefficient lagged through the sub-problem that carries it
(`χ` via `vof`, which introduces `χh`) is that carrier. Bound symbols are
substituted, and the sharp functions are smoothed as Proteus smooths them:
`heaviside` becomes `gf.H` of `equivalent_polynomials.h`, whose derivative
is its `gf.D`, and `sign` becomes `2 H_ε − 1`. Each sub-problem also
carries what it reads and when, its coefficients of the state by name
(with `inline` code and its ADR form over them in `named`, as for
`adr_problem`), the
regularizations it used, an introduced unknown's initial condition, its
effective solver, and its **boundary conditions**: those on the unknowns it
solves for, and those on a coefficient its introduced unknown carries
(atmosphere's `χ` inflow is VOF's `χh` inflow), with condition sets
expanded. A condition on a field a sub-problem only reads is another's.
A condition on a **phase indicator no sub-problem carries** (`χ :=
heaviside(φ)` with no VOF, as when CLSVOF solves the level set alone) is a
condition on the sign of the level set, given to the sub-problem that
solves it: `heaviside(φ) = 1` says `φ > 0` and no more. Its entry has
`variable: φ`, `sign` (`1` for χ = 1, `-1` for 0, `0` for 1/2), the
indicator's own `{variable, value}` and no `value`; a solver that imposes a
Dirichlet value picks one with that sign. So the atmosphere's `χ` inflow
of 1 (what enters is air) is VOF's `χh` inflow in `level_set_five` and the
positive `φ` inflow of CLSVOF in `level_set_clsvof`
([`dambreak_clsvof.ymf`](../tests/fixtures/two_phase/dambreak_clsvof.ymf)).
Anything else on such an indicator (a formula, another value) is refused.
A solver maps a sub-problem onto its own modules (a dispatcher); `proteus:`
blocks are interim hints for it, carried as text.

### The words a dispatcher maps

A solver's dispatcher maps the spec onto its modules by the words the spec
states, and those words are ymf's ([`ymf.vocabulary`](../src/ymf/vocabulary.py)),
checked on load, so that a second solver reads the same ones:

| kind | words |
|---|---|
| model family (`physical_model.family`) | `two_phase_level_set`, `shallow_water`, `dispersive_shallow_water` |
| sub-problem role (`role`) | two-phase: `velocity_pressure`, `velocity`, `pressure_increment`, `pressure`, `pressure_spin_up`, `volume_fraction`, `level_set`, `level_set_redistanced`, `redistancing`, `mass_correction`; shallow water: `shallow_water` |
| `time_integration` | `backward_euler`, `bdf2`, `ssp_rk3`, `none` |
| `stabilization.method` | `ASGS`, `SUPG`, `residual_shock_capturing`, `entropy_viscosity`, `graph_viscosity_smoothness`, `graph_viscosity_entropy_residual` |
| `discretization.scheme` and its options | `graph_viscosity_convex_limiting`: `entropy_viscosity_coefficient`, `lumped_mass`, `bounds_relaxation` (and `time_integration`, `outputs`, `dt` for every scheme) |
| `step_control.type` | `min_cfl`, `fixed` |
| `regularize.*.by` | `smoothed_heaviside`, `smoothed_sign`, `wet_dry_cutoff` |
| `binds.*.meaning` | `element_size`, `dry_depth` |
| `writes_back` | `replace`, `add`, `mass_quadrature` |
| `linear_solver` | `direct`, `iterative`; `superlu`, `mumps`, `umfpack`, `petsc` |
| coefficient `role` | `density`, `dynamic_viscosity`, `kinematic_viscosity`, `gravity`, `phase_indicator`, `surface_tension`, `bed_elevation`, `manning_coefficient`, `porosity`, `linear_drag`, `quadratic_drag`, `bed_friction`, `canopy_drag`, `diffusivity` |

A family is identified by the processes and assumptions of
`physical_model` (two-phase: `incompressible_flow`, `two_phase_flow`,
`interface_transport`; shallow water: `hydrostatic`, `depth_averaged`,
`nondispersive`; dispersive: `depth_averaged`, `dispersive`). A role is
read from the sub-problem's structure: what it solves (the unknowns its
equations are for), enforces, introduces and spins up. Both may be stated
(`family:`, `role:`), and are then checked against what the model and the
structure say; a closed spec writes both out (`ymf.explicit`), and neither
counts towards an output's key, since each follows from what does.

A **coefficient role** says what a coefficient means whatever the model
calls it: `ρ: {formula: "ρw (1 − χ) + ρa χ", units: kg/m3, role: density}`.
A solver that needs the density, the gravity or the bed finds it by its
role (`ymf.vocabulary.coefficient_by_role(doc, "density")`), not by the
name the family file happens to use. A role is used by one coefficient of
a model, and `ymf.units.check` refuses units of another dimension than the
role's.

Some words were first the dispatch table's, named after a solver's class
or flag; ymf's word replaces them and the old one is still read as an
alias (`ymf.vocabulary.ALIASES`, `canonical`): `vbdf` is `bdf2`, `SSP_RK3`
is `ssp_rk3`, `shock_capturing` (the method) is `residual_shock_capturing`,
`pressure_initial` is `pressure_spin_up`, `entropy_viscosity_cE` is
`entropy_viscosity_coefficient`, and the coefficient roles `mannings`,
`bathymetry`, `drag_alpha`, `drag_beta` are `manning_coefficient`,
`bed_elevation`, `linear_drag`, `quadratic_drag`. A spec using an alias
loads, and is not rewritten: a reader keyed on the old word keeps working
until it maps the new one. `dispersive_shallow_water` was identified by
the assumption `hyperbolic_relaxation`, which is how GN_SW2DCV discretizes
the model, not an assumption of it; that is still recognized.

## Document state

A YMF document is the record of a progression: from a model, through a
well-posed specification and its discretizations, to the closure that holds
what was computed and exactly what it was computed from. `ymf_state` says
where along that progression a document stands:

| `ymf_state` | the document | checked by |
|---|---|---|
| `abstract` | part of a model; it must be extended before it can be solved (equations and unknowns, say, with no domain, conditions or values) | the partial schema only |
| `well-posed` | a complete `model:` block | the full schema |
| `discretized` | well-posed, and every discretization has what its `method` needs (a mesh: `cells` or `h`), so every output it plans has a key | its planned outputs |
| `closed` | discretized, and every planned output is present and still keyed to the input that produced it | its outputs' keys |

Each state includes the ones before it. A document is as far as its least
advanced discretization, and whatever reports the state names each branch
that holds it back, with a hint: "P2_quadratic is missing a mesh: give
mesh: {h: ...} (or {cells: ...} on a box), or --h (--cells) on the command
line" (`ymf.state.branches` gives the same per branch).

A document may declare its state,
and the loader checks that its content reaches it; declaring less is
allowed. Tools that advance a document write the state they reach: a
runner writing an archive records `closed` when every output it plans is
there, and `discretized` while some are still owed. An archive whose input
was edited no longer reaches the `closed` it declares, and loading it says
so, naming the outputs that are now stale. (A runner bringing it up to
date loads it anyway, and recomputes.) A document that declares nothing
must be well-posed, as before.

## Generated documents

A tool that advances a document -- `pro weak` adding weak forms, `pro dis`
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

Earlier drafts, including [schema v0.1](history/schema-v0.1.md), called the
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

`vvuq` takes `verification` and `validation`, each with any content; the
reference graph checks `vvuq.verification.reference`. The intended shape of
the rest is in [schema v0.1 § 3](history/schema-v0.1.md#3-vvuq-results). The
top-level `archive` and `Xdmf` keys of the v0.1 draft are refused now: the
archive is its own format, [archive-format.md](archive-format.md).

## Known gaps

These are real limitations of the current schema and code, recorded here so
nobody has to rediscover them. Gaps that have been closed, and how, are in
[changes.md](changes.md).

- **A solver may not read every domain kind yet.** ymf poses all four;
  Proteus reads `discrete` domains of triangles or tetrahedra (in memory,
  cekees/proteus#200), and takes the bounds of a bounds-free file only from
  tetgen, 2dm and 3dm files.
- **Periodic conditions only along the axes.** A periodic region names
  two faces of the domain's bounding box, on a box or on a region whose two
  faces match; a rotation, or a shift that is not along an axis, is not
  expressible, and periodic conditions on a `file` domain are not
  supported (Gmsh's `Periodic` constraint is not carried).
- **The units of a flux condition are partly checked.** `neumann` and
  `robin` see the diffusive flux only, and `t` has no declared units
  anywhere, so an `α` or `β` that depends on time is never decided.
- **A duplicated value is found only between a scale and a coefficient.** A
  velocity scale that restates a boundary speed is not caught, because that
  would need the boundary condition's expression evaluated. Derive scales
  from coefficients wherever possible.
- **The notation's limits**, each of which reports the fix as well as the
  fault. Functions need parentheses: `sin πx` says to write `sin(x)`, not
  `sin x`. There are no subscripts: `p₀` says to write `p0` and declare it,
  since the superscripts `² ³ ⁴` *are* notation (powers) and a subscript is
  a reasonable guess. And a multi-letter name must be declared to be a name
  rather than a product: `Re` says to declare it, and that it is not read as
  `R·e` because `R` and `e` are not names.
