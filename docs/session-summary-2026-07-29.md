# YMF Session Summary — 2026-07-29

**Copy-paste this into a new session** to restore context and continue the work.

---

## Who we are / What we're doing

User (cekees) is the developer behind Proteus (a PDE solver framework) and the original ibvp project (Andreas Kloeckner's symbolic language for PDEs). We're designing **YMF** — a YAML-based problem description format for computational physics that supports LLM-assisted model selection, multi-backend dispatch, and built-in VVUQ.

## Core architecture

### The pipeline (6 stages, from physics to solution)

```
Stage 0: Physical processes → choose which to model
Stage 1: Equation formulation → strong form PDEs (conservative/non-conservative, dimensionless)
Stage 2: Solution variables → unknowns (u,p vs ρu,p etc.) — separate decision from equation
Stage 3: Weak form → theoretical reformulation (C²→H1, integration by parts)
Stage 4: Discretization → FEM/FD/FV, element type, order, mesh
Stage 5: Solver → linear/nonlinear, time stepping, tolerances
```

**LLM's role**: assists at stages 0-2 (model selection, formulation, unknown choice), then generates boilerplate for stages 3-5. Human reviews and approves at each step.

**YMF's role**: captures every stage as a structured, verifiable artifact. Each stage is a branch-point — new discretizations or weak forms can be added without modifying existing ones.

## YMF Schema (draft at `ymf-schema.md` in workspace)

### Structure (composable, append-only branches)

```yaml
Problem:
  physical_model: ...  # processes, assumptions, source docs
  strong_form:         # equation, unknowns, domain, BCs, coefficients
  weak_forms:          # theoretical reformulations (shared parent of FEM branches)
    - label: "global"  # human/LLM name; from_weak_form uses this or int index
      bilinear: |...|
      linear: |...|

solution_paths:
  analytical:          # strong form → analytical solution (no weak form needed)
  discretizations:     # numerical methods
    - name: "taylor_hood"
      from_weak_form: "global"  # int (0) or str ("global") — both resolve same entry
      finite_element: {velocity: {family: CG, order: 2}, ...}
      solver: {type: nonlinear, nonlinear_solver: newton, ...}
    - name: "fd_central"
      from_weak_form: null     # bypasses weak form — FD goes directly from strong form
      strategy: finite_difference
      discretization: {spatial: central_difference, order: 4, ...}

vvuq:  # computed after solving all branches
  verification: {grid_convergence, exact_comparison, MMS}
  validation: {comparison against real-world data}

archive:  # output config (mappable to XDMF but XDMF not embedded in YMF)
  time_collection: append | batch
  bundle: true
  include_ymf: true
```

### Key design decisions

1. **sympy over pymbolic** — sympy is the language layer (CAS, dsolve, geometry, integration). pymbolic is used only for UFL emission (FEniCS/Firedrake, since UFL is built on pymbolic). Conversion: sympy → pymbolic → UFL.

2. **Weak form is Phase 2**, not a discretization strategy. It's a theoretical reformulation (C²→H1, existence/uniqueness/regularity). Multiple FEM branches can share the same weak form with different discretizations. FD/FV bypass weak form entirely.

3. **Two-tier referencing**: `from_weak_form` accepts `int` (programmatic) or `str` (human/LLM label). Both resolve to the same array position.

4. **Provenance at every level**: `llm_derived`, `human_specified`, `proteus_derived`, `human_edited`. Used for trust calibration and debugging.

5. **XDMF is a consumer, not a dependency**. Archive data maps to XDMF but doesn't embed XDMF references. ParaView reads the archive; XDMF format isn't in the YMF schema.

6. **YAML + strictyaml** — YAML for human/LLM authoring (comments, flow style, block scalars), strictyaml for deterministic validation.

## Proteus 2.0 integration

- **Archiver module**: current state is a mess — multiple HDF5 approaches, XDMF extension pain. New state: archive module that simply writes XDMF output + embeds YMF provenance.
- **XML → YAML output**: YAML natively supports streaming time collection (append timestep-by-timestep). No more "save, overwrite, read-back, reconstruct" kludge.
- **UFL emission**: new target path — YMF → sympy → pymbolic → UFL → FEniCS/Firedrake kernel compilation.
- **ibvp resurrected**: original language layer (Andreas Kloeckner's sympy-based DSL) is revived. `target/proteus` transpiler existed; new `target/ufl` transpiler is added.

## ibvp original architecture (what was revived)

The original ibvp repo had:
- `ibvp/language/symbolic/` — sympy/pymbolic expression tree for PDEs
- `ibvp/language/geometry.py` — domain/geometry definitions
- `ibvp/target/proteus/` — transpiler to Proteus input
- Used pymbolic primitives (`Variable`, `d(x,y)`, etc.)

New architecture:
- Language layer: **sympy** (not pymbolic) — CAS capabilities needed for analytical solutions
- Targets: `proteus` (existing), `ufl` (new), `analysis` (new), `mms` (new)
- YAML input validated by strictyaml, converted to sympy expression tree

## VVUQ integration

YMF embeds VVUQ results natively — not an afterthought. The archive captures:
- **Verification**: grid convergence, MMS, exact comparison — error rates per branch per grid level
- **Validation**: comparison against real-world datasets with quantitative metrics (RMSE, R², etc.)
- **Uncertainty**: parameter distributions, propagation methods, sensitivity analysis

Multiple discretization branches are compared against each other and against data, all in one document.

## Files in workspace

- `ymf-schema.md` — full schema draft (v0.1)
- `memory/2026-07-29.md` — today's session notes
- User's Jupyter notebook (from earlier session) — Poisson XDMF data

## Open questions / next steps

1. **YMF name**: still not finalized as of this session. The candidate then was
   "YAML Data Model and Format" (YDMF), which felt forced. *Settled later:* the
   project is **YMF — YAML Modeling Format**, renamed throughout on 2026-09-02.
2. **sympy expression system**: design the language layer (how to parse strong form expressions from YAML, handle unknowns vs coefficients, boundary regions as sympy domains).
3. **UFL transpiler**: sympy → pymbolic → UFL code generation.
4. **Proteus 2.0 refactoring**: Archiver rewrite, UFL emission, ibvp integration.
5. **YAML → sympy parser**: how to convert `strong_form_expression` text (unicode, LaTeX) to a sympy expression tree reliably.

---

**End of session summary.**
