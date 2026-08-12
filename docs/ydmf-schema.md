# YDMF — Yet Another Data Model for Physics
# Draft Schema (v0.1)

A structured, human- and LLM-editable format for describing physical
models, numerical discretizations, and VVUQ results, with support for
branching across formulations, solution paths, and solver backends.

## Design Principles

- **Composability**: each branch adds discretizations or weak forms
  without modifying existing entries.
- **Provenance**: every block records how it was generated
  (`llm_derived`, `human_specified`, `proteus_derived`, `human_edited`).
- **Two-tier referencing**: integer index for programmatic use,
  string label for human/LLM readability.  Both resolve to the same
  array position.
- **XDMF preservation**: every valid XDMF 2.0/3.x document is a valid
  YDMF document (embedded under the `Xdmf` key).  ParaView and other
  XDMF consumers skip unknown top-level keys.

---

## Top-Level Structure

```yaml
# ============================================================
# YDMF — Yet Another Data Model for Physics (working name)
# ============================================================

# === 1. Problem definition (shared by all branches) ===
Problem:
  name: "Flow over flat plate"

# === 2. Discretizations (numerical approximation paths) ===
solution_paths:
  analytical: [...]
  discretizations: [...]

# === 3. VVUQ results (computed after solving) ===
vvuq:
  verification: { ... }
  validation: { ... }

# === 4. Output / archive configuration ===
archive:
  format: xdmf | vtk | hdf5
  time_collection: append | batch
  include_auxiliary: [residual_norm, solver_iterations, ...]

# === 5. Archive output (easily mappable to XDMF) ===
# A valid archive is a data tree plus metadata.  The data tree
# maps cleanly to XDMF without XDMF being a dependency — XDMF is
# only a consumer of the archive, not a building block of YDMF.
Archive:
  Domain: { ... }          # XDMF 2.0/3.x data tree
  Extension:               # non-XDMF metadata, skipped by standard consumers
    YDMF_provenance: "..."
```

---

## 1. Problem Definition

### 1.1 Physical Model (Phase 0)

```yaml
Problem:
  name: "Flow over flat plate"
  physical_model:
    provenance: llm_derived  | human_specified | proteus_derived
    processes:
      - incompressible_NS
    assumptions:
      - steady
      - constant_properties
      - no_body_forces
    source_document: "paper_2024_smith.pdf"
    notes: "RANS turbulence model not yet included"
```

### 1.2 Strong Form (Phase 1)

The strong form is the equation as it appears in literature — a
statement about PDEs, unknowns, domain, and boundary regions.

```yaml
  strong_form:
    unknowns: [u, p]                    # primary unknowns (Evan's u convention)
    unknown_provenance: llm_derived
    equation_formulation: "incompressible_NS_primitive"
    strong_form_expression: |
      ∂u/∂t + u·∇u = -∇p + νΔu + f
      ∇·u = 0
    domain: |
      Ω = [0,10] × [0,1]
      T ∈ [0, 1.0]
    boundary_regions:
      - name: "inlet"
        geometry: "x=0, 0≤y≤1"
      - name: "outlet"
        geometry: "x=10, 0≤y≤1"
      - name: "walls"
        geometry: "y=0 ∪ y=1, 0≤x≤10"
    # Optional: analytical solution (if known)
    known_analytical_solution:
      formula: "..."
      reference: "..."
    # Coefficients referenced in the strong form expression
    coefficients:
      ν: 1.0e-5                          # kinematic viscosity
      f: "sin(πx)cos(πt)"               # body force
```

Key: **`strong_form_expression`** is free-form text (unicode, LaTeX,
or sympy-compatible notation).  The parser/LLM converts it to a
sympy expression tree.

### 1.3 Weak Forms (Phase 2)

Weak forms are **theoretical reformulations** of the strong form.
They live under `weak_forms` (a list).  Each entry has a `label`
(string name) for human/LLM reference and an implicit array index.
Discretizations reference by **either** the label string or the
integer index.

```yaml
  weak_forms:
    - label: "supg_standard"            # ← human/LLM name
      provenance: llm_derived           # ← how this was generated
      derivation: |
        Multiply by test functions (v,q) ∈ H1×H1,
        integrate over Ω, integrate by parts on diffusion term.
        Apply essential BCs on trial space, natural BCs on test space.
      solution_spaces:
        trial: "H1 × H1"                # [u,p] ∈ H1(Ω) × H1(Ω)
        test:  "H1_0 × H1"              # v ∈ H1_0(Ω), q ∈ H1(Ω)
      inf_sup_stability: true           # LBB condition verified
      existence_note: "Lax-Milgram applies (with inf-sup condition)"
      uniqueness_note: "guaranteed under inf-sup stability"
      regularity_note: "u ∈ H2, p ∈ H1 (on convex domain)"
      bilinear: |
        B([(u,p)], [(v,q)]) = 
          ∫Ω ∇u:∇v dx - ∫Ω p(∇·v) dx + ∫Ω q(∇·u) dx
      linear: |
        F([(v,q)]) = ∫Ω f·v dx
      boundary_conditions: |
        u|Γ_inlet = u_inlet   (essential)
        p|Γ_outlet = 0              (natural, embedded in weak form)
        u|Γ_wall = 0                (essential)
      stabilization_method: SUPG
      notes: |
        Standard Galerkin with SUPG stabilization.
        Taylor-Hood elements (CG2/CG1) are LBB-stable for this form.

    - label: "cps"                       # ← human/LLM name
      provenance: human_specified        # ← expert correction
      derivation: |
        Same as weak form 0, but with Continuous Petrov-Galerkin
        (CPS) stabilization instead of SUPG.
      solution_spaces:
        trial: "H1 × H1"
        test:  "H1_0 × H1"
      inf_sup_stability: true
      bilinear: |
        B([(u,p)], [(v,q)]) = 
          ∫Ω ∇u:∇v dx - ∫Ω p(∇·v) dx + ∫Ω q(∇·u) dx
          - τ_CPS ∫Ω (u·∇u + ∇p):∇(v·∇v + ∇q) dx
      linear: |
        F([(v,q)]) = ∫Ω f·v dx
      stabilization_method: CPS
```

**Referencing**:

```yaml
# Integer index (programmatic):
from_weak_form: 0    # → weak_forms[0] = "supg_standard"

# String label (human/LLM readable):
from_weak_form: "supg_standard"   # → same weak_forms[0]

# Mixed: both are valid and resolve to the same entry.
```

### 1.4 Discretizations (Phase 3)

Discretizations are grouped under `solution_paths` into three
top-level categories:

1. **`analytical`** — strong form → analytical solution (classical
   C²), no weak form required.
2. **`discretizations`** — numerical approximation methods.  Each
   entry references a weak form by index or label (or `null` for
   FD/FV which bypass weak forms entirely).
3. **(future)** — additional path types (e.g., `model_reduction`)

#### 2.1 Analytical Path

```yaml
solution_paths:
  analytical:
    - name: "exact_poisson_reference"
      provenance: llm_derived
      from_weak_form: null        # bypasses weak form entirely
      method: "separation_of_variables"
      solution:
        formula: "u(x,y) = sin(πx)sin(πy)/(2π²)"
        notes: "used for verification only; not a Navier-Stokes solution"
        derivation_source: "sympy.dsolve"    # sympy derived this
        domain_restriction: "unit square only"
```

#### 2.2 Numerical Discretizations

```yaml
  discretizations:

    # FEM branch 0: from weak_form[0] = "supg_standard"
    - name: "taylor_hood_cg2"
      provenance: llm_derived
      from_weak_form: 0           # integer index
      finite_element:
        velocity: {family: CG, order: 2}
        pressure: {family: CG, order: 1}
      discretization_details: "Taylor-Hood P2/P1 on triangle mesh"
      solver:
        type: nonlinear
        nonlinear_solver: newton
        linear_solver: petsc
        tolerance: 1e-10
        max_iterations: 50

    # FEM branch 1: from weak_form[0] = "supg_standard"
    # Same weak form, different FEM spaces
    - name: "dg_upwind_2"
      provenance: llm_derived
      from_weak_form: "supg_standard"   # string label
      finite_element:
        velocity: {family: DG, order: 2}
        pressure: {family: DG, order: 1}
      discretization_details: "DG2/DG1 with upwind flux"
      flux_method: "upwind"
      numerical_flux: "Rusanov"
      solver:
        type: nonlinear
        nonlinear_solver: newton
        linear_solver: petsc

    # FEM branch 2: from weak_form[1] = "cps"
    - name: "taylor_hood_cps"
      provenance: human_specified
      from_weak_form: "cps"
      finite_element:
        velocity: {family: CG, order: 2}
        pressure: {family: CG, order: 1}
      discretization_details: "Taylor-Hood P2/P1 with CPS stabilization"
      solver:
        type: nonlinear
        nonlinear_solver: newton
        linear_solver: petsc

    # FD branch: from strong form directly (no weak form)
    - name: "fd_central_4th_order"
      provenance: llm_derived
      from_weak_form: null
      strategy: finite_difference
      discretization:
        spatial: "central_difference"
        order: 4
        time: "runge_kutta_4"
        grid:
          type: structured
          nx: 256
          ny: 64
          extent: [[0, 10], [0, 1]]
        staggered: false
      boundary_discretization:
        inlet: "extrapolate_u"
        outlet: "zero_gradient"
        wall: "no_slip_second_order"
      solver:
        type: time_marching
        dt: 1e-5
        max_time: 1.0
        CFL: 0.5
        max_iterations: 1000000

    # FV branch: from strong form directly (no weak form)
    - name: "fv_weno5_hllc"
      provenance: human_specified
      from_weak_form: null
      strategy: finite_volume
      discretization:
        reconstruction: "weno5"
        flux: "HLLC"
        quadrature: "Gauss_2pt"
        grid_type: unstructured
        mesh_file: "mesh.h5"
      solver:
        type: time_marching
        dt: 1e-6
        max_time: 0.5
        CFL: 0.5

    # Proteus-derived branch (deterministic from _p.py / _n.py)
    - name: "proteus_default_NS"
      provenance: proteus_derived
      from_weak_form: "supg_standard"
      proteus_modules:
        physics: "_p.py"
        numerics: "_n.py"
      finite_element:
        velocity: {family: DG, order: 2}
        pressure: {family: DG, order: 1}
      solver:
        type: nonlinear
        nonlinear_solver: newton
        linear_solver: petsc
```

#### 2.3 Weak-Form Discretization Details

The `finite_element` block describes the discretization of the weak
form itself:

```yaml
finite_element:
  # For FEM: element family + order per field
  velocity: {family: CG, order: 2}
  pressure: {family: CG, order: 1}
  # Alternative: mixed methods
  # velocity: {family: RT, order: 1}
  # pressure: {family: DG, order: 0}
  # Alternative: DG (same family for all fields)
  # fields: {family: DG, order: 2}

  # For non-FEM: strategy is finite_difference, finite_volume, spectral, ...
  # and discretization block replaces finite_element.

# Stabilization that modifies the weak form at the discrete level
stabilization:
  method: "SUPG" | "CPS" | "DG_upwind_flux" | "GLS" | "none"
  parameters:
    tau: "1/(2|u|/h + ν/h²)"  # stabilization parameter formula
```

---

## 3. VVUQ Results

Computed **after** solving all branches.

```yaml
vvuq:
  verification:
    method: "grid_convergence"   # or "mms" | "exact_comparison"
    grid_levels: [64, 128, 256]
    reference:
      type: "analytical"
      branch_ref: "exact_poisson_reference"
      # type: "manufactured"       # MMS
      # formula: "..."
    results:
      - branch: "taylor_hood_cg2"
        grids:
          - n: 64    error_L2: 1.2e-2  error_H1: 3.4e-2  rate_L2: null  rate_H1: null
          - n: 128   error_L2: 3.1e-3  error_H1: 8.7e-3  rate_L2: 2.0  rate_H1: 2.0
          - n: 256   error_L2: 7.8e-4  error_H1: 2.2e-3  rate_L2: 2.0  rate_H1: 2.0
        convergence_order: "second_order_L2"
      - branch: "taylor_hood_cps"
        grids:
          - n: 64    error_L2: 1.1e-2  error_H1: 3.2e-2  rate_L2: null  rate_H1: null
          - n: 128   error_L2: 2.8e-3  error_H1: 7.8e-3  rate_L2: 2.1  rate_H1: 2.0
          - n: 256   error_L2: 7.0e-4  error_H1: 2.0e-3  rate_L2: 2.0  rate_H1: 2.0

  validation:
    dataset:
      source: "NASA Turbine Blade Cooling Data"
      year: 1976
      doi: "10.xxxx/xxxxx"
    quantities_of_interest:
      - name: "NusseltNumber"
        location: "turbine blade surface"
      - name: "HeatTransferCoefficient"
        location: "turbine blade surface"
    comparison_method: "profile"  # or "pointwise" | "integral"
    results:
      - branch: "taylor_hood_cg2"
        metrics:
          RMSE: 0.034
          max_deviation: 0.12
          R2: 0.96
        notes: "Trend agreement good; peak values underpredicted by ~8%"
```

---

## 4. Archive Configuration

```yaml
archive:
  format: xdmf | vtk | hdf5
  time_collection:
    mode: append | batch        # append: stream timestep-by-timestep
                                # batch: write all timesteps at end
  output_interval: 0.1          # write every 0.1 time units (when mode=append)
  auxiliary_variables:
    include:
      - residual_norm
      - solver_iterations
      - condition_number
      - time_step_accepted
  bundle: true                  # single HDF5 file for all data
  compression: gzip             # or lzf, szip, none
  include_ydmf: true            # embed full YDMF as provenance
  version: "ydmf-0.1"
  software_versions:
    proteus: "2.0.0"
    fenics: "2023.1.0"
    firedrake: "2023.10.12"
  provenance:
    created_by: "..."
    git_commit: "..."
    compute_node: "..."
    date: "2026-07-29"
```

---

## 5. XDMF Preservation

```yaml
Xdmf:
  # Valid XDMF 2.0 or 3.x XML structure, embedded here.
  # Every valid XDMF document is a valid YDMF document.
  # ParaView and other XDMF consumers read the Xdmf key
  # and ignore all others (unknown top-level keys).
  Domain:
    Grid:
      - Name: "solution_field_u"
        GridType: "Uniform"
        Time: {Value: "0.1"}
        Topology: {Type: "Tetrahedron", NumberOfElements: "687737"}
        Geometry: {Type: "XYZ"}
        Attribute: {Name: "u", AttributeType: "Scalar", ...}

  # Extension elements: arbitrary XML that XDMF consumers skip.
  Extension:
    YDMF_provenance: "supg_standard → taylor_hood_cg2"
    weak_form_label: "supg_standard"
    discretization_index: 0
```

---

## Full Example

A complete working example (combining all sections):

```yaml
# ============================================================
# YDMF Example: Poisson on unit square
# ============================================================

Problem:
  name: "Poisson unit square"
  physical_model:
    provenance: llm_derived
    processes: [heat_conduction]
    assumptions: [steady, constant_properties]

  strong_form:
    unknowns: [u]
    unknown_provenance: llm_derived
    equation_formulation: "Poisson"
    strong_form_expression: "-Δu = f  in Ω"
    domain: |
      Ω = [0,1] × [0,1]
    boundary_regions:
      - name: "all_boundaries"
        geometry: "∂Ω"
    coefficients:
      f: "π²sin(πx)sin(πy)"

weak_forms:
  - label: "global"
    provenance: llm_derived
    derivation: |
      Multiply by test function v ∈ H1_0(Ω),
      integrate over Ω, integrate by parts.
    solution_spaces:
      trial: "H1"
      test:  "H1_0"
    inf_sup_stability: true
    existence_note: "Lax-Milgram applies (coercive bilinear form)"
    uniqueness_note: "guaranteed by coercivity"
    regularity_note: "u ∈ H2(Ω) on convex domain"
    bilinear: |
      a(u,v) = ∫Ω ∇u·∇v dx
    linear: |
      l(v) = ∫Ω f·v dx
    boundary_conditions: |
      u|∂Ω = 0  (essential)
    stabilization_method: "none"

solution_paths:
  analytical:
    - name: "exact_Poisson"
      provenance: llm_derived
      from_weak_form: null
      method: "separation_of_variables"
      solution:
        formula: "u(x,y) = sin(πx)sin(πy) / (2π²)"
        derivation_source: "sympy.dsolve"

  discretizations:
    - name: "P1_linear"
      provenance: llm_derived
      from_weak_form: "global"   # string label
      finite_element:
        fields: {family: CG, order: 1}
      solver:
        type: linear
        linear_solver: petsc
        tolerance: 1e-12

    - name: "P2_quadratic"
      provenance: llm_derived
      from_weak_form: "global"   # same weak form
      finite_element:
        fields: {family: CG, order: 2}
      solver:
        type: linear
        linear_solver: petsc
        tolerance: 1e-12

    - name: "FD_central_2nd"
      provenance: llm_derived
      from_weak_form: null            # bypasses weak form
      strategy: finite_difference
      discretization:
        spatial: "central_difference"
        order: 2
        grid:
          type: structured
          nx: 32
          ny: 32
      solver:
        type: linear
        linear_solver: petsc
        tolerance: 1e-12

vvuq:
  verification:
    method: "exact_comparison"
    reference:
      type: "analytical"
      formula: "u(x,y) = sin(πx)sin(πy) / (2π²)"
    results:
      - branch: "P1_linear"
        error_L2: 0.0412
        error_H1: 0.287
        convergence_rate: "first_order_L2, first_order_H1_semi"
      - branch: "P2_quadratic"
        error_L2: 0.0031
        error_H1: 0.062
        convergence_rate: "second_order_L2, second_order_H1_semi"
      - branch: "FD_central_2nd"
        error_L2: 0.0387
        error_H1: 0.301
        convergence_rate: "second_order_L2 (expected)"

archive:
  format: xdmf
  time_collection: batch
  bundle: true
  compression: gzip
  include_ydmf: true
  version: "ydmf-0.1"

Archive:
  Format:
    Grid:
      - Name: "solution_P1"
        GridType: "Uniform"
        Topology: {Type: "Triangle", NumberOfElements: "578"}
        Geometry: {Type: "XYZ"}
        Attribute:
          Name: "u"
          AttributeType: "Scalar"
          Center: "Node"
          DataItem:
            Format: "HDF"
            DataType: "Float"
            Precision: "8"
            Dimensions: "324 1"
            Hyperslab: "0 324 0 1 0 1 0 324 0 1"
```

---

## Schema Validation (strictyaml)

The YDMF file is validated with `strictyaml` before processing:

```python
from strictyaml import (
    Int, Str, Seq, Map, Enum, Optional, Float, Any,
    validator_with_error,
)

# Core validators
ProvenanceEnum = Enum(["llm_derived", "human_specified",
                        "proteus_derived", "human_edited"])

FEFamilyEnum = Enum(["CG", "DG", "RT", "BDM", "RTF", "N1curl"])

WeakFormDef = Map({
    "label": Str(),
    "provenance": ProvenanceEnum,
    "derivation": Str(),
    "solution_spaces": Map({
        "trial": Str(),
        "test": Str(),
    }),
    "inf_sup_stability": Optional(Str(), default=""),
    "existence_note": Optional(Str(), default=""),
    "uniqueness_note": Optional(Str(), default=""),
    "regularity_note": Optional(Str(), default=""),
    "bilinear": Str(),
    "linear": Str(),
    "boundary_conditions": Optional(Str(), default=""),
    "stabilization_method": Optional(Str(), default="none"),
    "notes": Optional(Str(), default=""),
})

FiniteElementDef = Map({
    Optional("velocity"): Map({
        "family": FEFamilyEnum,
        "order": Int(minimum=0),
    }),
    Optional("pressure"): Map({
        "family": FEFamilyEnum,
        "order": Int(minimum=0),
    }),
    Optional("fields"): Map({
        "family": FEFamilyEnum,
        "order": Int(minimum=0),
    }),
})

SolverConfig = Map({
    Optional("type"): Enum(["linear", "nonlinear", "time_marching", "saddle_point"]),
    Optional("nonlinear_solver"): Enum(["newton", "broyden", "line_search"]),
    Optional("linear_solver"): Enum(["petsc", "mumps", "umfpack", "superlu"]),
    Optional("tolerance"): Float(),
    Optional("max_iterations"): Int(minimum=1),
})

AnalyticalEntry = Map({
    "name": Str(),
    "provenance": ProvenanceEnum,
    "from_weak_form": Optional(Any(), default=None),
    "method": Str(),
    "solution": Map({
        "formula": Str(),
        Optional("notes"): Str(),
        Optional("derivation_source"): Str(),
        Optional("domain_restriction"): Str(),
    }),
})

DiscretizationEntry = Map({
    "name": Str(),
    "provenance": ProvenanceEnum,
    "from_weak_form": Optional(Any(), default=None),  # Int or Str or None
    "finite_element": Optional(FiniteElementDef),
    "discretization": Optional(Any()),  # strategy-specific
    "solver": SolverConfig,
    Optional("flux_method"): Str(),
    Optional("numerical_flux"): Str(),
    Optional("discretization_details"): Str(),
    Optional("proteus_modules"): Map({
        "physics": Str(),
        "numerics": Str(),
    }),
})

DiscretizationCategory = Map({
    "analytical": Seq(AnalyticalEntry),
    "discretizations": Seq(DiscretizationEntry),
})

StrongFormDef = Map({
    "unknowns": Seq(Str()),
    "unknown_provenance": ProvenanceEnum,
    "equation_formulation": Str(),
    "strong_form_expression": Str(),
    "domain": Str(),
    "boundary_regions": Seq(Map({
        "name": Str(),
        "geometry": Str(),
    })),
    Optional("known_analytical_solution"): Map({
        "formula": Str(),
        "reference": Str(),
    }),
    Optional("coefficients"): Map({
        Any(): Any(),
    }),
})

PhysicalModelDef = Map({
    "provenance": ProvenanceEnum,
    Optional("processes"): Seq(Str()),
    Optional("assumptions"): Seq(Str()),
    Optional("source_document"): Str(),
    Optional("notes"): Str(),
})

# Top-level schema
YDMF = Map({
    "Problem": Map({
        "name": Str(),
        Optional("physical_model"): PhysicalModelDef,
        Optional("strong_form"): StrongFormDef,
        Optional("weak_forms"): Seq(WeakFormDef),
    }),
    "solution_paths": DiscretizationCategory,
    Optional("vvuq"): Map({
        Optional("verification"): Any(),
        Optional("validation"): Any(),
    }),
    Optional("archive"): Any(),
    Optional("Xdmf"): Any(),  # preserve any valid XDMF
})
```

---

## Notes

- **`from_weak_form: null`** means the discretization does NOT
  reference any weak form.  This is valid for: analytical solutions,
  finite-difference, finite-volume, spectral, and any method that
  operates directly on the strong form.
- **`from_weak_form` can be `int` or `str`.**  Both resolve to the
  same array position.  Integers are programmatic; strings are
  human/LLM-readable labels.
- **Weak forms are a parent of FEM discretizations, not a
  sibling.**  Adding a new FEM branch is appending to the
  `discretizations` list — the shared weak form definition stays
  untouched.
- **Archive mapping**: the Archive section is a data tree that is
  easily mappable to XDMF (or other formats) without requiring XDMF
  as a dependency.  The mapping preserves semantics; XDMF is
  *consumed* via ParaView or similar tools, not embedded.
