# CSDMS Standard Names and BMI Integration for YDMF

A report on integrating CSDMS Standard Names and the Basic Model Interface (BMI) into the YAML Data Model and Format (YDMF) for computational physics.

---

## 1. Overview

CSDMS (the Center for Subsurface Modeling and Simulation at the University of Colorado Boulder) has developed two complementary standards that are highly relevant to YDMF's goals:

- **CSDMS Standard Names** — a registry of ~4,000+ standardized variable names for earth system and subsurface modeling variables, enabling interoperability across codes regardless of how each developer chose to name their variables.

- **BMI (Basic Model Interface)** — a component-based API specification that allows any model to be treated as a plug-and-play component within a larger coupled modeling framework.

Both standards align well with YDMF's dual goals: (1) a structured, human/LLM-friendly problem description format, and (2) multi-backend interoperability across solvers (Proteus, FEniCS, Firedrake, ParFlow, Delft3D, WRF, etc.).

---

## 2. CSDMS Standard Names

### 2.1 Registry Structure

The CSDMS Standard Names registry is a hierarchical vocabulary organized by domain:

| Domain | Example Variables |
|---|---|
| **atmosphere** | `air_temperature`, `air_density`, `surface_downwelling_shortwave_flux`, `surface_wind_speed`, `specific_humidity` |
| **land_surface** | `land_surface_temperature`, `soil_moisture_content`, `leaf_area_index`, `vegetation_fraction` |
| **surface_water** | `surface_water_depth`, `water_surface_elevation`, `surface_water_velocity_x`, `surface_water_velocity_y`, `water_surface_speed`, `surface_water_dimensionless_velocity_x` |
| **subsurface** | `subsurface_water_pressure`, `groundwater_head`, `soil_water_content`, `subsurface_water_velocity_x`, `porosity`, `hydraulic_conductivity` |
| **ocean** | `ocean_salinity`, `sea_surface_temperature`, `sea_surface_height_above_geoid`, `ocean_mixed_layer_thickness` |
| **ice** | `sea_ice_area_fraction`, `sea_ice_thickness`, `glacier_surface_elevation` |
| **chemistry** | `dissolved_oxygen_concentration`, `nitrate_concentration`, `phosphate_concentration` |

### 2.2 Naming Convention

Standard names follow a consistent pattern:

```
[domain]_[subdomain]_[property]_[modifier]
```

Examples:
- `surface_water_velocity_x` — velocity in the x-direction for surface water
- `land_surface_temperature` — temperature at the land surface
- `subsurface_water_pressure` — pressure in the subsurface water domain

### 2.3 Registry API

The registry is accessible programmatically at:
- **Web API**: https://csdms.colorado.edu/wiki/Standard_Names_API
- **Python client**: `csdms_stdname` package (`pip install csdms_stdname`)
- **JSON download**: The full registry is available as a structured JSON file

Each entry includes:
- `standard_name` (string): the canonical name
- `units` (string): expected SI units
- `description` (string): description of the quantity
- `long_name` (string): human-readable long form
- `is_valid` (boolean): whether the entry is active

### 2.4 Integration with YDMF

YDMF can integrate CSDMS Standard Names at multiple levels:

#### Level 1: Direct replacement (most explicit)

```yaml
# Current (bare names)
unknowns: [u, p]
coefficients:
  ν: 1.0e-5

# With CSDMS Standard Names
unknowns:
  - std_name: surface_water_velocity_x
    alias: u
  - std_name: water_surface_speed
    alias: v
  - std_name: surface_water_pressure
    alias: p
coefficients:
  kinematic_viscosity:
    std_name: null  # not in CSDMS registry
    value: 1.0e-5
    units: m2/s
```

#### Level 2: Metadata field (retains bare names)

```yaml
strong_form:
  unknowns: [u, p]
  std_names:
    u: surface_water_velocity_x
    p: surface_water_pressure
  coefficients:
    ν:
      alias: kinematic_viscosity
      std_name: null
      units: m2/s
```

#### Level 3: Validation-only (minimal change)

```yaml
# YDMF stays exactly as-is today
unknowns: [u, p]
strong_form_expression: "∂u/∂t + u·∇u = -∇p + νΔu + f"

# Separate validation layer (optional)
validation:
  std_name_check:
    unknowns:
      u: surface_water_velocity_x  # validates u maps to this
      p: surface_water_pressure
```

### 2.5 Benefits

1. **Interoperability**: Any solver reading YDMF can look up the standard name and map it to their internal variable representation, regardless of naming conventions.

2. **Unit awareness**: CSDMS entries carry expected SI units, which feeds into the units system (Report 2).

3. **Automatic documentation**: The `description` and `long_name` fields provide human-readable documentation for each variable.

4. **Cross-domain linkage**: Standard names explicitly encode the physical domain (surface_water, subsurface, etc.), which can drive automatic coupling between models in different domains.

### 2.6 Limitations

1. **Incomplete coverage**: Many PDE variable names (e.g., "u", "v", "p" for fluid dynamics, "φ" for electromagnetics) don't have direct CSDMS entries because the registry is earth-science focused.

2. **Dimensionality ambiguity**: `surface_water_velocity_x` is inherently 3D but only specifies one component. Multi-component vector fields need multiple entries.

3. **Not a substitute for physics**: A standard name says *what* a variable is, not *how* it behaves in the PDE. YDMF still needs the equation-level specification.

4. **Not all variables fit**: Some YDMF variables (e.g., numerical artifacts, solver diagnostics) have no natural standard name.

---

## 3. Basic Model Interface (BMI)

### 3.1 Overview

BMI is a component-based API that allows any model to be treated as a "plug-and-play" component in a coupled modeling framework. It defines a minimal set of methods that every model must implement, enabling model interoperability without requiring shared codebases or file formats.

### 3.2 Core API

Every BMI-compliant model implements these methods:

| Category | Method | Purpose |
|---|---|---|
| **Lifecycle** | `initialize(config)` | Set up the model from configuration |
| | `finalize()` | Clean up resources |
| | `update()` | Advance the model one time step (legacy) |
| | `update_every(n_steps)` | Advance n steps (newer, more efficient) |
| | `take_step()` | Advance one time step (legacy alias) |
| | `update_connection()` | Update data exchange between components |
| **Time** | `get_start_time()` | Return model start time |
| | `get_end_time()` | Return model end time |
| | `get_current_time()` | Return current simulation time |
| | `get_time_step()` | Return time step size (Δt) |
| | `get_time_units()` | Return time units (e.g., "seconds") |
| **Variables** | `get_var_type(name)` | Return numpy dtype for a variable |
| | `get_var_type_string(name)` | Return type as string |
| | `get_var_units(name)` | Return SI units for a variable |
| | `get_var_location(name)` | Return location ("surface", "near_surface", "subsurface") |
| | `get_var_grid(name)` | Return grid ID for a variable |
| | `get_var_ndims(name)` | Return number of dimensions |
| | `get_var_size(name)` | Return size in elements |
| | `get_var_itemsize(name)` | Return size in bytes per element |
| | `get_var_offset(name)` | Return byte offset in memory |
| | `get_input_var_ndims(name)` | Return ndims for an input variable |
| | `get_output_var_ndims(name)` | Return ndims for an output variable |
| | `get_input_var_grid(name)` | Return grid for an input variable |
| | `get_output_var_grid(name)` | Return grid for an output variable |
| | `get_input_var_type(name)` | Return type for an input variable |
| | `get_output_var_type(name)` | Return type for an output variable |
| | `get_input_var_units(name)` | Return units for an input variable |
| | `get_output_var_units(name)` | Return units for an output variable |
| | `get_input_var_location(name)` | Return location for an input variable |
| | `get_output_var_location(name)` | Return location for an output variable |
| **Data** | `get_value(name)` | Get a copy of the variable's data |
| | `get_value_at_indices(name, indices)` | Get values at specific indices |
| | `get_value_slice(name, start, stop)` | Get a slice of the variable |
| | `get_value_ptr(name)` | Get pointer to variable data (in-place modification) |
| | `set_value(name, value)` | Set the variable's data |
| | `set_value_at_indices(name, indices, values)` | Set values at indices |
| **Grid** | `get_grid_type(name)` | Return grid type ("rectilinear", "structured", "unstructured", etc.) |
| | `get_grid_size(name)` | Return number of nodes |
| | `get_grid_ndims(name)` | Return number of dimensions |
| | `get_grid_shape(name, shape)` | Fill shape array |
| | `get_grid_spacing(name, spacing)` | Fill spacing array |
| | `get_grid_origin(name, origin)` | Fill origin array |
| | `get_grid_x(name, x)` | Fill x-coordinate array |
| | `get_grid_y(name, y)` | Fill y-coordinate array |
| | `get_grid_z(name, z)` | Fill z-coordinate array |
| | `get_grid_edges(name, edges)` | Fill edge arrays |
| | `get_grid_nodes(name, nodes)` | Fill node arrays |
| | `get_grid_faces(name, faces)` | Fill face arrays |
| | `get_grid_center_coordinates(name, coords)` | Fill center coordinate arrays |
| **Registration** | `get_component_name()` | Return model name |
| | `add_var(name, var_type, is_input, is_grid_var, location)` | Register a variable |
| | `add_grid(name, grid_type, size, ndims, spacing, origin, x, y, z, is_regular)` | Register a grid |
| | `remove_grid(name)` | Remove a grid |

### 3.3 BMI Implementation Languages

BMI is implemented in:
- **Python**: `pymt.bmi` (https://github.com/csdms/bmi)
- **C**: `libbmi` (https://github.com/csdms/bmi-c)
- **Fortran**: `bmi-fortran` (https://github.com/csdms/bmi-fortran)
- **Rust**: (community implementations)

### 3.4 Integration with YDMF

YDMF can serve as the **configuration layer** for BMI-compliant models. The YDMF file describes the problem specification, and a BMI adapter reads the YDMF file and exposes the problem to the solver as a BMI component.

#### Example: YDMF → BMI Adapter

```python
from pymt.bmi import BMI

class YDMFBMIAdapter(BMI):
    """Adapter that exposes a YDMF problem specification as a BMI component."""

    def initialize(self, config_path):
        # Parse the YDMF file
        self.ydmf = load_ydmf(config_path)

        # Register variables from YDMF weak_form definition
        self._register_variables()

        # Register grids from YDMF discretization
        self._register_grids()

        # Set up the solver using YDMF specification
        self.solver = build_solver(self.ydmf)
        self.solver.initialize()

    def get_var_units(self, name):
        # Look up unit from YDMF coefficient/unknown definition
        return self._get_unit_from_ydmf(name)

    def get_var_location(self, name):
        # Look up location from YDMF boundary condition definition
        return self._get_location_from_ydmf(name)

    def get_grid_spacing(self, name, spacing):
        # Look up grid spacing from YDMF discretization grid definition
        return self._get_spacing_from_ydmf(name, spacing)

    def take_step(self):
        self.solver.step()

    def get_value(self, name):
        return self.solver.get_field(name)

    def set_value(self, name, value):
        self.solver.set_field(name, value)
```

#### YDMF → BMI Variable Mapping

| YDMF Section | BMI Method | Mapping |
|---|---|---|
| `strong_form.unknowns` | `add_var()` | Register each unknown as an output variable |
| `strong_form.coefficients` | `add_var()` | Register each coefficient as an input variable |
| `weak_forms[].solution_spaces` | `add_grid()` | Register grid for each solution space (trial/test) |
| `discretizations[].grid` | `add_grid()` | Register discretization grid |
| `discretizations[].boundary_conditions` | `get_var_location()` | Register boundary condition locations |
| `strong_form.domain` | `add_grid()` | Register computational domain grid |
| `solution_paths[].solver.*` | `initialize()` | Pass solver configuration |

#### Example YDMF → BMI Data Flow

```yaml
# YDMF input
Problem:
  strong_form:
    unknowns: [u, p]
    coefficients:
      ν: 1.0e-5
    domain: |
      Ω = [0,1] × [0,1]

discretizations:
  - name: "P1_linear"
    from_weak_form: "global"
    finite_element: {fields: {family: CG, order: 1}}
    solver: {type: linear}
```

The BMI adapter would register:
- **Output variables**: `u`, `p` (the unknowns)
- **Input variables**: `ν` (kinematic viscosity)
- **Grids**: one grid for the computational domain Ω
- **Time methods**: derived from the solver configuration

### 3.5 Multi-Component Coupling

BMI's real power is in coupling multiple models. YDMF can describe a **system of coupled components**:

```yaml
# Multi-component YDMF (future extension) — illustrative example only
System:
  components:
    - name: "fluid_solver"
      bmi_interface: true
      ydmf_file: "fluid.ydmf"
      coupling:
        inputs: [fluid_velocity, fluid_pressure]
        outputs: [fluid_force]
    - name: "solid_solver"
      bmi_interface: true
      ydmf_file: "solid.ydmf"
      coupling:
        inputs: [fluid_force, solid_stress]
        outputs: [solid_displacement]
    - name: "boundary_manager"
      bmi_interface: true
      ydmf_file: "boundary.ydmf"
      coupling:
        inputs: [solid_displacement]
        outputs: [boundary_velocity]
  coupling_strategy: "strong"  # or "weak" (partitioned)
  coupling_interval: 0.01
```

This approach would enable YDMF to drive multi-physics simulations (fluid-structure interaction, heat transfer in flowing fluids, etc.) by connecting BMI-compatible components.

### 3.6 Existing BMI-Compatible Models

Many earth system models are already BMI-compliant or have adapters:
- **ParFlow** (subsurface flow)
- **Delft3D** (surface water hydrodynamics)
- **WRF** (weather/atmospheric modeling)
- **CLM** (land surface model)
- **CISM** (ice sheet model)
- **ParMOSEK** (optimization)
- **PyMT** (the framework itself, provides BMI wrappers)

YDMF could serve as a **problem specification format** that works across all of these.

---

## 4. Integration Strategies

### 4.1 Minimal Viable Integration

**What can be done first (low effort, high value):**

1. **Add `std_name` field** to YDMF unknowns and coefficients, alongside existing bare names:
   ```yaml
   unknowns:
     - name: u
       std_name: surface_water_velocity_x
   ```

2. **Add `units` field** to YDMF coefficients and unknowns:
   ```yaml
   coefficients:
     ν:
       value: 1.0e-5
       units: m2/s
   ```

3. **Create a `std_names.py` module** that maps YDMF variable names to CSDMS standard names using the `csdms_stdname` package.

### 4.2 Intermediate Integration

**What requires more infrastructure:**

1. **BMI adapter layer** that reads YDMF files and exposes them to BMI-compatible solvers.

2. **Automated variable lookup** that suggests CSDMS standard names from the YDMF's `strong_form_expression` (parse the expression, identify variables, propose standard names).

3. **Unit validation** that checks YDMF coefficients against CSDMS expected units.

### 4.3 Long-Term Vision

1. **YDMF as a universal problem specification format** that works with any BMI-compliant model framework, not just Proteus or FEniCS.

2. **Automatic coupling** where YDMF describes multi-component systems and the BMI adapter handles data exchange between components.

3. **Registry integration** where YDMF variable lookups against the CSDMS registry happen at load time, providing warnings if a variable name doesn't match any standard.

---

## 5. Challenges and Limitations

### 5.1 Scope Mismatch

The CSDMS Standard Names registry is earth-science focused. YDMF aims to support general computational physics (including non-earth-science domains like structural mechanics, electromagnetics, chemistry, etc.). **Not all PDE variables have CSDMS standard names.**

**Mitigation**: Use CSDMS standard names *where available* but allow arbitrary names for variables outside the registry.

### 5.2 BMI Grid Assumption

BMI assumes a grid-based data structure. YDMF supports grid-based methods (FEM on meshes, FD on structured grids) but also non-grid methods (spectral methods, meshfree methods, global bases).

**Mitigation**: BMI's `get_grid_type()` method returns `"unstructured"` or other types, which can accommodate various discretizations. For truly non-grid methods, the BMI adapter could create a "virtual grid" (e.g., spectral coefficients arranged on a conceptual grid).

### 5.3 Unknown vs. Input/Output

BMI has a clear separation between input variables and output variables. YDMF's unknown/coefficient distinction doesn't map cleanly:
- Some variables are both input and output (e.g., velocity is input to the momentum equation but output from the solver)
- Coefficients can change over time (e.g., viscosity as a function of temperature)
- Boundary conditions can be time-dependent

**Mitigation**: Use the YDMF `solution_paths[].solver.type` to determine which variables are inputs vs outputs. For steady-state problems, all variables are effectively outputs. For time-dependent problems, coefficients are inputs and unknowns are outputs.

### 5.4 Performance

BMI's Python API has overhead compared to direct calls to solver libraries. For very large-scale simulations, the BMI adapter layer could add non-trivial overhead.

**Mitigation**: BMI is designed for model coupling, not for maximizing throughput within a single solver. Use BMI for coupling YDMF-specified problems; for single-solver performance, bypass BMI and call the solver directly.

---

## 6. Recommendations

### For YDMF Schema (immediate)

1. **Add `std_name` and `units` fields** to YDMF's unknowns, coefficients, and boundary conditions. These are optional, backward-compatible additions.

2. **Add `bmi_interface` flag** to discretizations, indicating whether the solver should be exposed via BMI.

3. **Document the CSDMS standard name lookup convention** — e.g., "if a variable name matches a CSDMS standard name, use it; otherwise use the bare name."

### For BMI Integration (short-term)

1. **Create a `ydmf_bmi.py` module** that loads a YDMF file and returns a BMI-compliant adapter object.

2. **Support the most common BMI methods** first: `initialize`, `finalize`, `take_step`, `get_var_units`, `get_var_location`, `get_grid_spacing`.

3. **Test with at least one existing BMI-compatible model** (ParFlow or Delft3D) to validate the adapter.

### For Long-Term

1. **Consider a "components" section** in YDMF for multi-component coupled problems, using BMI for inter-component communication.

2. **Integrate with PyMT** (Python Modeling Toolkit) as the coupling framework, since it already provides BMI wrappers for many models.

3. **Explore BMI's grid API** for automated mesh generation from YDMF's domain/geometry specifications.

---

## 7. Conclusion

CSDMS Standard Names and BMI are highly complementary to YDMF's goals. The standard names provide a semantic layer for variable identification and unit awareness, while BMI provides a mechanical layer for model interoperability and coupling.

The integration doesn't have to be all-or-nothing. YDMF can adopt these standards incrementally: first adding `std_name` and `units` fields to the schema, then creating a BMI adapter layer, and finally supporting multi-component coupling.

The biggest challenge is scope mismatch — CSDMS standard names are earth-science focused, while YDMF aims for general computational physics. However, this is a manageable limitation: use standard names where available, fall back to arbitrary names otherwise.

**Bottom line**: CSDMS and BMI are not dependencies of YDMF, but they are natural extensions that would significantly improve YDMF's interoperability without compromising its core design goals.
