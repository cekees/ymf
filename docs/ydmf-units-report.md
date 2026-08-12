# Physical Units Support for YDMF

A report on adding physical unit support to YDMF in a non-intrusive way that enables dimensional consistency checks and non-dimensionalization.

---

## 1. Current Landscape

### 1.1 Python Unit Libraries

| Library | Strengths | Weaknesses | Best For |
|---|---|---|---|
| **pint** | Widely used, flexible unit definitions, auto-conversion, registry of ~10,000 units, supports custom units, arithmetic with units | Slower than bare numbers, sympy integration is separate | General Python unit handling |
| **unyt** (yt project) | Built for scientific computing, Numpy integration, C/Python interoperability, built-in astropy integration | Less flexible than pint, narrower adoption | Astrophysics, fluid dynamics |
| **astropy.units** | Standard in astronomy, well-tested, supports physical constants, C/Python interoperability | Domain-specific (astronomy), not general-purpose | Astropy ecosystem |
| **sympy.physics.units** | Native sympy integration, symbolic manipulation with units, dimensional analysis built-in | Slow, limited unit definitions, incomplete coverage | Symbolic computation with units |
| **quantities** (pyquantities) | Numpy-style API, unit conversion, simple syntax | Minimal documentation, small community | Simple scientific scripts |

**Recommendation**: Use **pint** as the primary unit library for YDMF, with **sympy.physics.units** for symbolic dimensional analysis. Pint's extensive unit registry and conversion capabilities complement sympy's symbolic manipulation.

### 1.2 Scientific Formats

| Format | Units Support | Notes |
|---|---|---|
| **CF Conventions** | Rich — `units`, `standard_name`, `axis`, `positive`, `units` attributes on every variable | The gold standard for climate/weather data. Defines canonical units and standard names. |
| **NetCDF** | Built-in variable attributes for units (`units`, `long_name`, `standard_name`) | Widely supported, CF conventions extend it |
| **XDMF** | Optional metadata in `DataItem` elements | No built-in unit system; units stored as XML attributes |
| **OASIS** | Unit conversion between coupled models | Focused on model coupling, not problem specification |

---

## 2. Design Principles

### 2.1 Non-Intrusive by Default

The primary design constraint is **non-intrusiveness**. YDMF should work perfectly without any unit annotation:

```yaml
# This is valid YDMF today — no units needed
strong_form:
  unknowns: [u, p]
  coefficients:
    ν: 1.0e-5

# Adding units is optional
strong_form:
  unknowns:
    - name: u
      units: m/s   # optional annotation
    - name: p
      units: Pa
  coefficients:
    ν:
      value: 1.0e-5
      units: m2/s  # optional annotation
```

**Guiding principle**: Unit annotation is an *enhancement*, not a requirement. A YDMF file without units is still a valid, usable problem specification.

### 2.2 Dimensional Compatibility

When units ARE provided:

- All terms in an equation term must have the same dimensions
- Boundary/initial conditions must be dimensionally compatible with the equation
- Discretization parameters (grid spacing, time step) must have appropriate dimensions
- The system validates these constraints but does not *require* them

### 2.3 Multi-Unit System Support

YDMF should support multiple unit systems:

| System | Units | Typical Use |
|---|---|---|
| **SI (MKS)** | m, kg, s, K, mol, A, cd | Scientific computing, general physics |
| **CGS** | cm, g, s, erg | Astrophysics, some chemistry |
| **Imperial/US** | ft, lb, s, °F | US engineering, environmental modeling |
| **Dimensionless** | No units | Non-dimensional analysis, benchmark cases |

A single YDMF file should be able to specify which unit system it uses, and the system can convert between them automatically.

---

## 3. Unit Annotation Strategies

Three approaches, ranked from least to most intrusive:

### Strategy A: Global Units Table (Minimal change)

```yaml
# YDMF stays mostly unchanged
Problem:
  name: "Flow over flat plate"
  strong_form:
    unknowns: [u, p]
    coefficients:
      ν: 1.0e-5

# Separate global units table
units:
  convention: SI
  variables:
    u: m/s
    p: Pa
    ν: m2/s
    t: s
    x: m
    y: m
```

**Pros**: Minimal change to existing YDMF structure; backwards compatible; clear separation of values and units.

**Cons**: Loose coupling between variable and unit definitions; no per-variable unit checking at the YAML level.

### Strategy B: Annotated Variable Names (Moderate change)

```yaml
Problem:
  strong_form:
    unknowns: ["u [m/s]", "p [Pa]"]
    coefficients:
      ν_m2_s: 1.0e-5
```

**Pros**: Units travel with variable names; no separate table needed.

**Cons**: Variable names become harder to use programmatically; not LLM-friendly (hard to parse variable + unit from string).

### Strategy C: Inline Unit Specification (Most explicit)

```yaml
Problem:
  strong_form:
    unknowns:
      - name: u
        units: m/s
      - name: p
        units: Pa
    coefficients:
      ν:
        value: 1.0e-5
        units: m2/s
```

**Pros**: Most explicit; enables automatic validation; LLM-friendly (structured fields); integrates naturally with CSDMS standard names (which also have units).

**Cons**: More verbose; slightly more complex schema.

### 3.1 Recommendation: Hybrid (A + C)

Use **Strategy C** (inline unit specification) for explicit annotation, with **Strategy A** (global units table) for bulk specification:

```yaml
Problem:
  strong_form:
    unknowns:
      - name: u
        units: m/s   # or reference to global
      - name: p
        units: Pa

  # Optional global units fallback
  units:
    convention: SI
    default: true    # all bare numbers use these units
    variables:
      ν: m2/s

  # Or per-domain unit sets
  unit_systems:
    SI:
      length: m
      mass: kg
      time: s
      temperature: K
    US:
      length: ft
      mass: lb
      time: s
      temperature: F
```

This approach:
- Maintains backwards compatibility (bare numbers are still valid)
- Enables full dimensional analysis when units are specified
- Is LLM-friendly (structured fields with clear semantics)
- Integrates with pint/unyt for automatic conversion
- Can reference CSDMS standard names (which carry units)

---

## 4. Dimensional Consistency Checking

### 4.1 Dimensional Representation

YDMF should track dimensions using the standard SI dimensional notation:

```
[M]^a [L]^b [T]^c [Θ]^d [mol]^e [A]^f [cd]^g
```

For most computational physics, the relevant dimensions are:

| Dimension | Symbol | Typical Quantity |
|---|---|---|
| Mass | M | density, force |
| Length | L | velocity, grid spacing |
| Time | T | time step, time, kinematic viscosity |
| Temperature | Θ | temperature, thermal expansion |
| Current | A | electromagnetic problems |
| Substance | mol | chemistry problems |

### 4.2 Equation-Level Validation

The strong form expression should be checked for dimensional consistency:

```yaml
strong_form:
  strong_form_expression: "-Δu = f"
  dimensional_check:
    method: "automatic"   # or "manual", "off"
    tolerance: 1e-10       # for floating-point dimensionless checks
```

The validation process:
1. Parse the `strong_form_expression` to extract terms (e.g., `-Δu` and `f`)
2. Look up the units of each variable from the YDMF specification
3. Check that all terms have the same dimensions
4. Report any mismatches

### 4.3 Discretization Parameter Validation

```yaml
discretizations:
  - name: "fd_central"
    from_weak_form: null
    strategy: finite_difference
    discretization:
      spatial: "central_difference"
      order: 4
      grid:
        nx: 256
        ny: 64
        extent: [[0, 10], [0, 1]]
        units: [m, m]        # ← dimension check: extent matches spatial units
      time: "runge_kutta_4"
      dt: 1e-5               # ← dimension check: dt has time units
```

The validation checks:
- Grid extent has length dimensions [L]
- Grid spacing (derived from extent/nx, ny) has length dimensions [L]
- Time step dt has time dimensions [T]
- dt/grid_spacing has appropriate dimensions for the equation (e.g., for advection: [T]/[L] = [T⁻¹L⁻¹] must be compatible with velocity [L/T])
- CFL number is dimensionless (implicit check: dt * velocity / grid_spacing is dimensionless)

### 4.4 Boundary/Initial Condition Validation

```yaml
strong_form:
  unknowns:
    - name: u
      units: m/s

boundary_conditions:
  - region: "inlet"
    type: "dirichlet"
    value: 5.0
    units: m/s    # ← must match unknown's units
```

The validation checks:
- BC value has the same dimensions as the unknown it applies to
- IC value has the same dimensions as the unknown at t=0
- BC type is dimensionally appropriate (e.g., Neumann BC value has dimensions of derivative: [u]/[L])

### 4.5 Sympy Integration

Use `sympy.physics.units` for dimensional analysis of the strong form expression:

```python
from sympy.physics.units import Dimension, Quantity, dimensions
from sympy.physics.units import velocity, length, time

# Check if an expression is dimensionally consistent
u = Quantity('u')
f = Quantity('f')
lap_u = Dimension(velocity / length)  # dimensions of ∇²u

# If both terms have the same dimensions, the equation is consistent
assert dimensions(lap_u) == dimensions(f)
```

This enables:
- **Automatic dimensional analysis** of the strong form expression
- **Unit propagation** through the equation (e.g., if u has units m/s and the equation involves u·∇u, the term has units m²/s²)
- **Dimensionless number calculation** (Reynolds, Froude, Grashof numbers, etc.)

### 4.6 Pint Integration

Use pint for practical unit conversion and validation:

```python
import pint

# Define the unit registry
ureg = pint.UnitRegistry()

# Validate that all variables have compatible units
u_units = ureg.meter / ureg.second        # m/s
p_units = ureg.pascal                     # Pa = kg/(m·s²)
nu_units = ureg.meter**2 / ureg.second    # m²/s

# Check that the equation terms are dimensionally consistent
# For Laplacian: Δu has units u/L → (m/s)/m = 1/s
# For f: should also have units 1/s
# Check: ν·Δu → (m²/s)·(1/s) = m²/s² ≠ 1/s → dimension mismatch!
# (This is why the weak form, not strong form, is typically used for NS)
```

---

## 5. Non-Dimensionalization

### 5.1 Why Non-Dimensionalize?

- **Numerical stability**: dimensionless numbers of order 1 are easier for solvers
- **Parameter space reduction**: fewer parameters to explore
- **Generalizability**: results apply to entire classes of problems (e.g., "high Reynolds number flow" applies to many geometries and fluids)
- **Scaling**: results from one simulation can predict behavior at different scales

### 5.2 Characteristic Scales in YDMF

```yaml
Problem:
  name: "Flow over flat plate"

  # Characteristic scales for non-dimensionalization
  characteristic_scales:
    length:
      value: 1.0
      units: m
    velocity:
      value: 5.0
      units: m/s
    time:                        # optional — derived from length/velocity
      derived: true
      formula: "L_char / U_char"
    pressure:                    # optional — derived from ρU²
      derived: true
      formula: "rho_char * U_char^2"
    density:
      value: 1.225
      units: kg/m3
    viscosity:
      value: 1.81e-5
      units: kg/(m·s)

  # Derived dimensionless numbers
  dimensionless_numbers:
    Reynolds:                     # Re = ρUL/μ
      value: (1.225 * 5.0 * 1.0) / 1.81e-5
      formula: "rho * U * L / mu"
    Froude:                       # optional
      value: 5.0 / (9.81 * 1.0)**0.5
      formula: "U / (g * L)**0.5"
    Mach:                         # optional
      value: 5.0 / 343.0
      formula: "U / c"
```

### 5.3 Automatic Non-Dimensionalization

YDMF can store both the dimensional and dimensionless forms:

```yaml
# Dimensional form (for input to solver)
strong_form:
  strong_form_expression: "∂u/∂t + u·∇u = -∇p + νΔu + f"
  coefficients:
    ν:
      value: 1.81e-5
      units: m2/s

# Dimensionless form (for analysis/verification)
strong_form_dimensionless:
  strong_form_expression: "∂u*/∂t* + u*·∇*u* = -∇*p* + (1/Re)Δ*u* + f*"
  coefficients:
    Reynolds:
      value: 276243
      units: dimensionless    # or just omit units

  non_dimensionalization:
    variables:
      u*: "u / U_char"
      t*: "t * U_char / L_char"
      x*: "x / L_char"
      p*: "p / (rho * U_char^2)"
      f*: "f * L_char / (U_char^2)"
    parameters:
      Re: "rho * U_char * L_char / mu"
```

The system can:
1. **Convert dimensional → dimensionless**: automatically substitute the non-dimensionalization relations
2. **Convert dimensionless → dimensional**: multiply by characteristic scales
3. **Run simulations in either form**: solver can operate on the dimensional or dimensionless form
4. **Compare results**: verification/error metrics can be computed in either form

### 5.4 Non-Dimensionalization Examples

**Poisson Equation** (already dimensionless if u is a scalar):
```
Dimensional:     -Δu = f    [u] = m²/s (e.g., heat equation with u = temperature)
Non-dimensional: -Δu* = f*
  u* = (u - u_ref) / Δu_ref
  x* = x / L_char
  f* = f / (f_ref)
```

**Navier-Stokes**:
```
Dimensional:     ∂u/∂t + u·∇u = -∇p/ρ + νΔu + f
Non-dimensional: ∂u*/∂t* + u*·∇*u* = -∇*p* + (1/Re)Δ*u* + f*
  where: u* = u/U_char, t* = t·U_char/L_char, p* = p/(ρU²)
  and:   Re = ρUL/μ
```

**Diffusion Equation**:
```
Dimensional:     ∂u/∂t = D·Δu
Non-dimensional: ∂u*/∂t* = Δ*u*
  where: u* = u/u_ref, t* = t·D/L_char²
  and:   the dimensionless number is the Fourier number (Fo = D·t/L²)
```

---

## 6. Implementation Recommendations

### 6.1 Minimal Viable Implementation (Version 1)

```yaml
# Schema additions to YDMF
units:
  convention: SI | CGS | US | dimensionless   # optional, default SI
  variables:                                    # optional mapping
    variable_name: "unit_string"                # e.g., "u": "m/s"
  characteristic_scales:                        # optional
    length: {value: 1.0, units: m}
    velocity: {value: 5.0, units: m/s}
    time: {derived: true, formula: "L/U"}
```

```python
# Python module for dimensional analysis
import pint

def check_dimensional_consistency(ydmf):
    """Check that all terms in strong_form have compatible dimensions."""
    ureg = pint.UnitRegistry()

    # Parse unknowns and coefficients
    unknown_units = {k: ureg(v) for k, v in ydmf['units']['variables'].items()}

    # For each equation term, check dimensionality
    for equation in ydmf['Problem']['strong_form']['equations']:
        terms = extract_terms(equation['expression'])
        term_dims = [look_up_dimension(term, unknown_units) for term in terms]

        if not all_dims_consistent(term_dims):
            raise DimensionalInconsistencyError(f"Equation {equation['name']} has dimensionally inconsistent terms")
```

### 6.2 What Requires Infrastructure

1. **Pint integration**: Install `pip install pint` (already a common dependency)
2. **Sympy unit integration**: Use `sympy.physics.units` for symbolic dimensional analysis
3. **Dimension checker module**: A `ydmf.check_units()` function that validates dimensional consistency
4. **Non-dimensionalization module**: A `ydmf.non_dimensionalize()` function that generates the dimensionless form
5. **Unit conversion in output**: XDMF output should include unit metadata (via CF Conventions-style attributes)

### 6.3 Making It LLM-Friendly

- **Unit auto-completion**: If an LLM generates `ν: 1.0e-5`, the system suggests `units: m2/s` (kinematic viscosity)
- **Unit validation feedback**: If an LLM writes `dt: 1e-5` with `units: m` (wrong dimension), the system warns "expected time dimension for time step"
- **Unit suggestion from context**: If the equation involves velocity, suggest `[L/T]` for all velocity-related variables
- **Dimensionless number suggestions**: When the LLM generates characteristic scales, suggest relevant dimensionless numbers (Re, Fr, Gr, etc.)

### 6.4 XDMF Output Integration

XDMF doesn't have a built-in unit system, but it can carry unit metadata:

```xml
<!-- XDMF output with unit metadata (CF Conventions style) -->
<Attribute Name="velocity" AttributeType="Vector" Center="Node">
  <DataItem Format="HDF" DataType="Float" Precision="8" Dimensions="118084 3">
    solution.h5:/velocity
  </DataItem>
  <Extension>
    <!-- CF Conventions attributes -->
    <units>m/s</units>
    <long_name>velocity</long_name>
    <standard_name>surface_water_velocity</standard_name>
  </Extension>
</Attribute>
```

The XDMF output module can automatically add unit attributes from the YDMF specification.

### 6.5 Integration with CSDMS Standard Names

CSDMS standard names carry expected units. The dimensional checker can validate YDMF units against CSDMS expectations:

```python
# If YDMF has:
#   unknowns:
#     - name: u
#       units: m/s
# And CSDMS standard_name is "surface_water_velocity":
#   Then check: does "m/s" match the expected units for "surface_water_velocity"?
#   CSDMS says: "surface_water_velocity" has units "m/s" → VALID
#   CSDMS says: "surface_water_velocity" has units "ft/s" → WARN (convertible but not SI)
```

This creates a **three-layer unit system**:
1. **YDMF explicit units** (from the YAML file)
2. **CSDMS standard name units** (from the registry)
3. **Physical reality** (what the units should actually be)

---

## 7. Concrete Examples

### 7.1 Poisson Problem (Simple)

```yaml
# Minimal — no units
Problem:
  name: "Poisson unit square"
  strong_form:
    unknowns: [u]
    equation: "-Δu = f in Ω"
    coefficients:
      f: "π²sin(πx)sin(πy)"

# With units
Problem:
  name: "Poisson heat conduction"
  strong_form:
    unknowns:
      - name: T                  # temperature
        units: K
    equation: "-kΔT = Q in Ω"
    coefficients:
      k:                          # thermal conductivity
        value: 1.0
        units: W/(m·K)
      Q:                          # heat source
        value: "π²"               # in W/m³
        units: W/m3

  units:
    characteristic_scales:
      length: {value: 1.0, units: m}
      temperature: {value: 300, units: K}   # reference temperature
```

### 7.2 Navier-Stokes (Complex)

```yaml
Problem:
  name: "Navier-Stokes, incompressible"
  strong_form:
    unknowns:
      - name: u
        units: m/s
      - name: p
        units: Pa
    equation: |
      ∂u/∂t + u·∇u = -∇p/ρ + νΔu + f
      ∇·u = 0

    coefficients:
      ρ:
        value: 1.225
        units: kg/m3
      ν:
        value: 1.5e-5
        units: m2/s
      f:
        value: 0.0
        units: m/s2

    characteristic_scales:
      length: {value: 1.0, units: m}
      velocity: {value: 10.0, units: m/s}

    # Derived dimensionless number
    dimensionless_numbers:
      Reynolds:
        formula: "ρ * U * L / μ"   # Note: ν = μ/ρ, so Re = U·L/ν
        value: 666667              # 10 * 1 / 1.5e-5

    # Dimensional consistency check
    dimensional_check:
      enabled: true
      tolerance: 1e-10

discretizations:
  - name: "FD_central_4th_order"
    discretization:
      grid:
        nx: 256
        ny: 64
        extent: [[0, 10], [0, 1]]
        units: [m, m]              # length dimensions
      dt: 1e-5
      dt_units: s                  # time dimension

# Output units
output:
  fields:
    velocity:
      units: m/s
    pressure:
      units: Pa
    time:
      units: s
```

**Dimensional consistency checks that pass**:
- `∂u/∂t` has units: [L/T]/[T] = [L/T²]
- `u·∇u` has units: [L/T]·[L/T]/[L] = [L/T²] ✓
- `∇p/ρ` has units: [M/(L·T²)]/[M/L³] = [L²/T²]/[L] = [L/T²] ✓
- `νΔu` has units: [L²/T]·[L/T]/[L²] = [L/T²] ✓
- `f` has units: [L/T²] ✓

All terms have the same dimensions [L/T²] — **consistent**.

**Dimensional inconsistency example** (user error):
```yaml
# User accidentally wrote:
velocity: 10.0
units: m     # should be m/s!
```

The dimensional checker would flag: `velocity` has units [L] but should have [L/T] (from the context of the Navier-Stokes equation).

---

## 8. Summary of Recommendations

### For YDMF Schema (immediate)

1. **Add optional `units` field** to unknowns, coefficients, and boundary conditions
2. **Add optional `characteristic_scales`** section for non-dimensionalization
3. **Add optional `dimensional_check`** flag for validation

### For Python Infrastructure (short-term)

1. **Install pint** — `pip install pint`
2. **Create `ydmf.units.check()`** — validates dimensional consistency of a YDMF file
3. **Create `ydmf.units.non_dimensionalize()`** — converts dimensional YDMF to dimensionless form
4. **Create `ydmf.units.auto_suggest()`** — suggests units based on context (LLM-friendly)

### For XDMF Output (short-term)

1. **Add unit metadata** to XDMF output (CF Conventions-style attributes)
2. **Include dimensionless numbers** in archive metadata

### For Long-Term

1. **Integrate with CSDMS Standard Names** — validate YDMF units against expected units for each standard name
2. **Automatic non-dimensionalization** — generate dimensionless forms from dimensional YDMF automatically
3. **Unit-aware solvers** — some solvers benefit from knowing units; pass units through the pipeline

---

## Appendix: SI Units Reference

| Dimension | Symbol | Base Units | Common Derived Units |
|---|---|---|---|
| Length | L | m | km, mm, μm, ft, in |
| Mass | M | kg | g, tonne, lb |
| Time | T | s | min, h, day, yr |
| Temperature | Θ | K | °C, °F |
| Current | A | A | — |
| Amount | mol | mol | — |
| Luminosity | cd | cd | — |

| Quantity | SI Unit | Symbol | Dimensions |
|---|---|---|---|
| Velocity | m/s | | [L/T] |
| Acceleration | m/s² | | [L/T²] |
| Force | N = kg·m/s² | | [M·L/T²] |
| Pressure | Pa = N/m² | | [M/(L·T²)] |
| Energy | J = N·m | | [M·L²/T²] |
| Power | W = J/s | | [M·L²/T³] |
| Viscosity (dynamic) | Pa·s = kg/(m·s) | | [M/(L·T)] |
| Viscosity (kinematic) | m²/s | | [L²/T] |
| Thermal conductivity | W/(m·K) | | [M·L/(T³·Θ)] |
| Heat capacity | J/(kg·K) | | [L²/(T²·Θ)] |
