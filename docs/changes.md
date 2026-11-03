# Closed gaps and the decisions behind them

The guides describe what YMF does. This page records the gaps that have been
closed, newest first, with the commit that closed each, so that nobody has to
rediscover why a thing is the way it is. Open gaps are in
[problem-specification.md § Known gaps](problem-specification.md#known-gaps).

## 2026-11

- **2026-11-03: `slip`, `open`, `inflow`, `characteristic`; conditions on
  a coefficient.** The boundary kinds two-phase flow and shallow water need:
  a free-slip wall, an open boundary (C12's declared open boundary), a value
  imposed where the flow enters, a stage imposed through the incoming
  characteristics. A condition may name a coefficient of the state (`χ`).
  Conditions are now checked on load: the region exists, the variable is an
  unknown or a coefficient, `slip` is on a vector, and a value is given
  exactly when the kind takes one. `ymf.units` checks `slip`, `inflow` and
  `characteristic` like `dirichlet`. `adr_problem` refuses the four kinds
  and conditions on coefficients, naming the solver they are for.
- **2026-11-03: each sub-problem classified on its own
  (`ymf.symbolic.split.split_problem`).** The plain data of a split
  discretization, one ADR form per sub-problem. Its own unknowns are
  unknowns; what it reads from the others is a known field (a function of
  the coordinates and time, written by name, its gradient
  `grad_<field>_<axis>`), and must be declared in `lags`; a coefficient
  lagged via the sub-problem that carries it is replaced by the carrier.
  `binds` are substituted (a symbol with only a meaning, `h`, stays a
  parameter the solver supplies), and `regularize` smooths `heaviside` and
  `sign` exactly as Proteus does (`gf.H`, `gf.D`, `2 H_ε − 1`). Schema:
  a sub-problem's `for` and `regularize`, a lag's `of`. The two-phase dam
  break classifies in all five sub-problems: flow as variable-density
  Navier-Stokes with ρ and μ from VOF's carrier, VOF and the level set as
  transport in the lagged velocity (`v·∇φ` is a Hamiltonian in `φ` alone
  there), redistancing as an eikonal Hamiltonian, mass correction as a
  regularized reaction-diffusion.
- **2026-11-03: coefficients of the state keep their names.**
  `adr_problem` lists each coefficient that follows the state under
  `coefficients`: its numpy code, `depends_on`, and its derivative in each
  unknown it depends on (`sharp: true` instead, for `heaviside` or `sign`
  of an unknown). The ADR form keeps them substituted, so a solver of one
  system is unaffected; the names are for a generator or a dispatcher. A
  problem without such coefficients has no `coefficients` key, so its
  record is unchanged.
- **2026-11-03: a discretization may be split into sub-problems.**
  `subproblems:` (solved in order each step: what each `solves`,
  `enforces`, `introduces`, `lags`, `writes_back`, with its own space,
  time integration, stabilization and a `solver` that refines the
  entry's), `binds:` (symbols such as `h` and `ε = 1.5 h` that the
  discretization supplies), `regularize:` (sharp functions made smooth) and
  `step_control:`. The structure and every name are checked on load
  (`ymf.references`); `ymf.discretizations.subproblem_solver` gives a
  sub-problem's effective solver. Classification per sub-problem is not
  done yet. Splitting is not a new `method`: a split discretization is
  still a mesh one, and its studies and output keys work as before.
- **2026-11-03: `heaviside`, `sign`, `min`, `max`, `norm`, `sym`.** The
  two-phase model writes `χ = heaviside(φ)`, its initial signed distance
  with `min`, `max` and `sqrt`, its viscous stress with `sym(∇v)`;
  shallow-water friction uses `norm(q)`. They parse, print back
  (`to_notation`), and become numpy code (`min`/`max` as nested
  `numpy.minimum`/`numpy.maximum`, not `functools.reduce`). The derivative
  of `heaviside` or `sign` is a Dirac delta, which has no value at a point:
  where the ADR form needs one (a sharp function of an unknown), it is
  refused with a message saying the discretization must regularize it.
- **2026-11-03: named equations.** `strong_form.equations` takes
  `{name, for, equation}` beside bare strings. Bare strings are still
  ordered like the unknowns; with `for` on every equation the order written
  stops mattering. Names merge across `extends` like every other named
  list, and are what a discretization's sub-problems will cite. Checked on
  load: a name used twice, `for` naming something that is not an unknown,
  two equations for one unknown, `for` on some but not all, and (outside an
  abstract file) an unknown no equation is for. `ymf.equations_in_order`
  is the one reader, shared by `ymf.symbolic` and `ymf.units`.
- **2026-11-03: a coefficient's formula may name the unknowns.** A
  coefficient could be a formula in the coordinates, time and the other
  coefficients; one that named an unknown was "could not be resolved". Now
  it may: a state-dependent density, conductivity or friction law is a
  coefficient whose formula names the unknowns, substituted into the
  equations so the ADR form carries the dependence and its derivatives.
  Decided by cekees for the model-family work (two-phase level-set flow:
  `χ = heaviside(φ)`, `ρ(χ)`, `μ(χ)`; shallow water: Manning friction,
  canopy drag), in place of a separate `relations:` block: one kind of
  thing, a coefficient, whether it is constant, varies in space or follows
  the state. What came with it: boundary data, initial conditions and exact
  solutions that lean on such a coefficient are refused (they are evaluated
  at points in x, y, z, t); `ymf.units` warns when a characteristic scale or
  a dimensionless number is built on one; a coefficient may not take an
  unknown's name; a formula may not take a derivative of an unknown (yet).
  `ymf.symbolic.state_coefficients` lists which coefficients follow the
  state, and on what.
- **2026-11-01: one name per thing: `model` and `approximation`.** The spec's
  top-level key has been `model` since a66ea1c, and the archive of a spec
  stores the realized spec under the key the spec uses, so it has written
  `model` since then too. The README, the archive guide and several
  docstrings still showed `Problem:`; they now say `model:`. `Problem:` and
  `problem:` still read as `model:`.
  The functions that build and check the tree a solver writes were named
  for a *domain*, which is also the region a problem is posed on
  (`model.strong_form.domain`, `ymf.domain`). They are now named for an
  *approximation*, which is what the archive of a spec already called that
  tree (`outputs[key].approximation`): `new_approximation`,
  `validate_approximation`, `canonicalize_approximation`,
  `select_approximation`, `rebase_approximation`, `inline_approximation`,
  `approximation_arrays`. The old names are the same functions, kept for one
  release; a release should warn before they go. The solver file's
  `domain:` key and the functions' parameter names are unchanged: the key
  is part of `ymf_archive_version: 2`.
- **2026-11-01: eight keys that nothing read are refused.** The schema
  accepted them from the v0.1 and v0.2 drafts, and no code in ymf or
  Proteus read or wrote them, so a document carrying one said something
  nobody heard. A document that still carries one is now told which key:
  - `strong_form.strong_form_expression`: superseded by `equations`, which
    ymf.symbolic reads.
  - `strong_form.known_analytical_solution`: superseded by
    `solution_paths.analytical` and `vvuq.verification.reference`.
  - `strong_form.dimensional_check`: the unit check is `ymf.units.check`,
    run by whoever wants it; a switch in the document turned nothing on.
  - `model.unit_systems`: nothing defined what a named system selects; units
    are given on each unknown, coefficient and condition.
  - `discretizations[].bmi_interface`: no BMI adapter exists.
  - top-level `archive` and `xdmf` (and the v0.1 spelling `Xdmf`, which is
    no longer folded): the archive is its own format, written by a solver
    and by `ymf.closure`, not a block of the spec.

  `mesh_generation` and `roughness_source` stay for now: `ymf.data_sources`
  reads them, and whether that module stays is a separate decision.

## 2026-10

- **2026-10-27: a hot start repeating a time replaces.** The two open
  questions of [design/time-outside-space.md](design/time-outside-space.md)
  § 8 were answered on the Proteus side (cekees/proteus#121): a hot start
  replaces the steps after the one it resumes from, and a space first
  archived mid-run is covered by a writer-level test.
- **2026-10-26: time outside, space inside (archive format 2).** The domain
  used to hold one temporal collection per function space, because XDMF
  does; it now holds one list of instants with the spaces inside each.
  Format 1 is still read and upgraded in memory. See
  [design/time-outside-space.md](design/time-outside-space.md)
  (proteustoolkit/ymf#22, #23; `a3eac96`).
- **2026-10-17: a boundary condition's units follow its type** (`54765f8`,
  proteustoolkit/ymf#9). They used to be compared with the unknown's units
  whatever the type, which is right only for `dirichlet`: a `neumann` flux in
  W/m² failed and one in K passed. Each type now has its own rule (see
  [Units](problem-specification.md#units-scales-and-dimensionless-numbers)).
  What is left is in Known gaps.
- **2026-10-16: `piecewise_linear` domains are solved** (`54a9f8f`). Any
  region that is not a box is meshed by `h`, the element size of an
  unstructured mesh. Periodic conditions still need a box.
- **2026-10-15: a table is a value form, not a kind of condition**
  (`6de07fe`). `hydrograph` used to be a `type`, which tangled the value
  with the condition and made a Dirichlet stage hydrograph impossible to
  write. `series:` is now orthogonal to `type`, as `value` and `formula` are.
- **2026-10-14: Neumann and Robin conditions run** (`600e38e`, `ad10e70`).
  Robin needed a numerical flux written for Proteus
  (`NumericalFlux.RobinFactory`), not wiring: the flux depends on the
  solution's trace, and a Proteus flux boundary condition is a function of
  `(x, t)` alone. That the *dependency set* of a boundary expression is the
  thing to get right, rather than the list of condition types, is the point
  taken up in the tag-propagation design notes.
- **2026-10-08: the reference graph is checked** (`e467ff1`). A
  `from_weak_form` or `branch_ref` that named nothing used to validate. All
  four edges are now checked, `branch_ref` became
  `reference.analytical`, integer indices into `weak_forms` are gone, and a
  single-valued reference may be omitted where the graph does not branch.
- **2026-10-08: numbers are numbers, and expressions are evaluated once**
  (`487182b`, `c86adf1`). `coefficients` was untyped, so every scalar arrived
  as a string, and a coefficient written as an expression was dropped from
  formula scope. Numbers are now typed in the schema, and expressions are
  evaluated by the one resolver `adr_problem` uses, so `ymf.units` and the
  solver give the same answer for the same coefficient. `sympy` moved into
  the `frontend` extra for this. The decision behind it (cekees): **numbers
  and expressions are treated uniformly**, because a front-end that
  evaluates expressions differently from the solver would report one thing
  and run another, and a `markdown | tex -> ymf` path is an intended
  direction.
- **2026-10-08: `1e-10` is a number to every reader** (`487182b`). It used
  to be a string to PyYAML and a number to strictyaml, so a spec meant
  different things read as a file and read inside an archive. The archive
  reader resolves YAML 1.2 scientific notation, on a loader subclass so
  pyyaml's own loaders are unchanged.
- **2026-10-08: `characteristic_scales` takes any name** (`487182b`). It
  was a fixed set of six. The original six keep their aliases (`L`, `U`,
  `rho`, `mu`, `T`) in formula scope.
- **2026-10-08: composition may remove** (`8452cd6`). Decided by cekees:
  removal is a real modelling move (incompressible solved as compressible,
  inviscid as viscous, dropping an unknown or replacing it with an
  analytical expression), and it must be logged. `remove:` names paths and
  is recorded in `composition.removals`.
- **2026-10-08: an error names the files behind it** (`d187db1`). An error
  that no single file holds, such as a required key none of them sets, was
  reported against a line of the merged document, which exists nowhere on
  disk. The message now names the path and the files that write it.
- **2026-10-08: duplicates are found by value and dimension** (`08773b8`).
  They were caught by name only. A scale whose value *and dimension* match a
  coefficient's is now reported whatever either is called. Dimension is what
  separates a duplicate from a coincidence: in a non-dimensionalized problem
  half the quantities are 1.0. What is left is in Known gaps.
- **2026-10-08: SUPG/PSPG with equal-order P1 on Kovasznay flow is
  pre-asymptotic, not wrong.** Investigated and closed: measured rates sit
  below the a priori values on coarse meshes and climb towards them as the
  mesh refines (1.88 for velocity at 512 cells). It is a resolution limit,
  not a defect in `NavierStokesASGS_velocity_pressure`. Taylor-Hood on the
  same problem reaches its rates sooner, so a side-by-side comparison at one
  mesh size misleads. Reaching the asymptotic range needs an iterative
  solver (512 cells is already 24 GB with MUMPS) or a lower-Reynolds
  companion case.
