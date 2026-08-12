# YDMF Schema v0.2 — Consolidation Delta

**Date**: 2026-08-12
**Status**: Draft — proposed additive/breaking changes to `ydmf-schema.md` (v0.1), synthesized from four background reports:

- `csdms-bmi-report.md` — CSDMS Standard Names + BMI component interoperability
- `ydmf-units-report.md` — physical units, dimensional consistency, non-dimensionalization
- `data-sources-report.md` — topography/bathymetry/vegetation/roughness/remote-sensing data sourcing
- `session-summary-2026-07-29.md` — original v0.1 design rationale

This document is a **delta**, not a replacement. It lists exactly what changes in the v0.1 schema, flags the one breaking change, and gives the merged strictyaml validators. Once reviewed, fold this into `ydmf-schema.md` proper and bump the version string.

---

## 1. Summary of Changes

| # | Section | Change Type | Source Report(s) |
|---|---|---|---|
| 1 | `unknowns` | **Breaking** (type widened) | units, csdms-bmi |
| 2 | `coefficients` | Additive (already loose enough) | units, csdms-bmi |
| 3 | `initial_conditions` (new) | Additive (new key) | data-sources, units |
| 4 | `boundary_conditions` (new, typed) | Additive (new key, supersedes bare `boundary_regions` for value-carrying BCs) | data-sources, units |
| 5 | `characteristic_scales` / `dimensionless_numbers` (new) | Additive | units |
| 6 | `dimensional_check` (new) | Additive | units |
| 7 | `mesh_generation` (new, under `Problem`) | Additive | data-sources |
| 8 | `roughness_source` (new, under discretization) | Additive | data-sources |
| 9 | `bmi_interface` flag (new, under discretization) | Additive | csdms-bmi |
| 10 | `units` / `unit_systems` block (new, under `Problem`) | Additive | units |
| 11 | `std_name` field (new, on unknowns/coefficients) | Additive | csdms-bmi |

**Only #1 is breaking.** Everything else is additive (`Optional(...)` in strictyaml terms) and preserves every valid v0.1 document.

---

## 2. The One Breaking Change: `unknowns`

### v0.1 (current)

```yaml
strong_form:
  unknowns: [u, p]              # Seq(Str())
```

### v0.2 (proposed)

`unknowns` must accept **either** a bare string (v0.1 compatibility) **or** a map with `name` + optional `units`/`std_name` (v0.2 enrichment):

```yaml
strong_form:
  # Still valid (v0.1 style):
  unknowns: [u, p]

  # Also valid (v0.2 enriched style):
  unknowns:
    - name: u
      units: m/s
      std_name: surface_water_velocity_x
    - name: p
      units: Pa
      std_name: surface_water_pressure
```

**Why this is breaking**: a strict `Seq(Str())` validator rejects the map form outright. Any parser doing `for name in unknowns: ...` assuming bare strings will break on v0.2 documents unless it normalizes first.

**Mitigation**: provide a `normalize_unknowns()` helper that always returns `List[{"name": str, "units": Optional[str], "std_name": Optional[str]}]` regardless of which form was used, so downstream code only ever sees one shape.

```python
def normalize_unknowns(unknowns):
    """Normalize v0.1 (bare str) or v0.2 (map) unknowns to a uniform list of dicts."""
    result = []
    for u in unknowns:
        if isinstance(u, str):
            result.append({"name": u, "units": None, "std_name": None})
        else:
            result.append({"name": u["name"], "units": u.get("units"), "std_name": u.get("std_name")})
    return result
```

---

## 3. New Sections (Additive)

### 3.1 `initial_conditions` / `boundary_conditions`

v0.1 only had `boundary_regions` (name + geometry, no value). v0.2 adds typed, value-carrying conditions with provenance for external data sources:

```yaml
strong_form:
  boundary_regions:             # unchanged from v0.1 — pure geometry
    - name: "inlet"
      geometry: "x=0, 0≤y≤1"

  initial_conditions:           # NEW
    - field: "u"
      type: "constant" | "function" | "data_source"
      value: 0.0                          # for constant
      formula: "sin(pi*x)*sin(pi*y)"      # for function
      source: "CHIRPS"                    # for data_source
      source_url: "https://..."
      temporal: "daily"                   # temporal resolution, if time-varying
      units: "m/s"

  boundary_conditions:          # NEW — references boundary_regions by name
    - region: "inlet"
      variable: "u"
      type: "dirichlet" | "neumann" | "robin" | "periodic" | "hydrograph"
      value: 5.0
      units: "m/s"
      source: null                        # or external source, see below
      source_url: null
    - region: "upstream"
      variable: "discharge"
      type: "hydrograph"
      source: "USGS NWIS"
      source_url: "https://waterdata.usgs.gov/nwis/"
      temporal: "hourly"
```

### 3.2 `characteristic_scales` / `dimensionless_numbers`

```yaml
Problem:
  characteristic_scales:
    length:   {value: 1.0, units: m}
    velocity: {value: 5.0, units: m/s}
    time:     {derived: true, formula: "L_char / U_char"}
    density:  {value: 1.225, units: kg/m3}
    viscosity: {value: 1.81e-5, units: "kg/(m*s)"}

  dimensionless_numbers:
    Reynolds: {formula: "rho * U * L / mu", value: 276243}
    Froude:   {formula: "U / (g * L)**0.5"}
```

### 3.3 `dimensional_check`

```yaml
strong_form:
  dimensional_check:
    enabled: true
    method: "automatic" | "manual" | "off"
    tolerance: 1e-10
```

### 3.4 `units` / `unit_systems` (Problem-level)

```yaml
Problem:
  units:
    convention: SI | CGS | US | dimensionless
    default: true
  unit_systems:
    SI: {length: m, mass: kg, time: s, temperature: K}
    US: {length: ft, mass: lb, time: s, temperature: F}
```

### 3.5 `mesh_generation` (Problem-level)

```yaml
Problem:
  mesh_generation:
    source: "SRTM 1arcsec"
    source_url: "https://www.usgs.gov/centers/eros/data-archive/srtm-90m-digital-elevation-database-v4"
    format: "GeoTIFF"
    vertical_datums: {horizontal: "WGS84", vertical: "EGM96"}
    resolution_target: 30
    filter: {method: "smoothing", kernel: "gaussian", sigma: 2}
```

### 3.6 `roughness_source` (per-discretization)

```yaml
discretizations:
  - name: "fd_with_manning"
    roughness_source:
      type: "land_cover_mapping"
      land_cover_data:
        source: "USGS NLCD"
        resolution: "30 m"
        url: "https://www.mrlc.gov/"
      reference: "Arcement and Schneider (1989)"
      lookup_table: "manning_n_lookup.csv"
```

### 3.7 `bmi_interface` (per-discretization)

```yaml
discretizations:
  - name: "taylor_hood_cg2"
    bmi_interface: true
    # ... existing fields unchanged
```

---

## 4. Merged strictyaml Validators (v0.2)

Only the changed/new validators are shown; everything else from v0.1 (`ProvenanceEnum`, `FEFamilyEnum`, `WeakFormDef`, `FiniteElementDef`, `SolverConfig`, `AnalyticalEntry`, `DiscretizationCategory`, `PhysicalModelDef`) is unchanged and still applies.

```python
from strictyaml import (
    Int, Str, Seq, Map, Enum, Optional, Float, Any, Bool,
)

# --- 1. Unknowns: widened to accept str OR map (BREAKING) ---
UnknownDef = Str() | Map({
    "name": Str(),
    Optional("units"): Str(),
    Optional("std_name"): Str(),
})

# --- 2. Initial conditions (NEW) ---
InitialConditionDef = Map({
    "field": Str(),
    "type": Enum(["constant", "function", "data_source"]),
    Optional("value"): Float(),
    Optional("formula"): Str(),
    Optional("source"): Str(),
    Optional("source_url"): Str(),
    Optional("temporal"): Str(),
    Optional("units"): Str(),
})

# --- 3. Boundary conditions (NEW) ---
BoundaryConditionDef = Map({
    "region": Str(),
    "variable": Str(),
    "type": Enum(["dirichlet", "neumann", "robin", "periodic", "hydrograph"]),
    Optional("value"): Float(),
    Optional("units"): Str(),
    Optional("source"): Str(),
    Optional("source_url"): Str(),
    Optional("temporal"): Str(),
})

# --- 4. Characteristic scales / dimensionless numbers (NEW) ---
ScaleDef = Map({
    Optional("value"): Float(),
    Optional("units"): Str(),
    Optional("derived"): Bool(),
    Optional("formula"): Str(),
})

CharacteristicScalesDef = Map({
    Optional("length"): ScaleDef,
    Optional("velocity"): ScaleDef,
    Optional("time"): ScaleDef,
    Optional("density"): ScaleDef,
    Optional("viscosity"): ScaleDef,
    Optional("temperature"): ScaleDef,
})

DimensionlessNumberDef = Map({
    "formula": Str(),
    Optional("value"): Float(),
})

# --- 5. Dimensional check (NEW) ---
DimensionalCheckDef = Map({
    Optional("enabled"): Bool(),
    Optional("method"): Enum(["automatic", "manual", "off"]),
    Optional("tolerance"): Float(),
})

# --- 6. Units / unit systems (NEW, Problem-level) ---
UnitSystemDef = Map({
    Optional("length"): Str(),
    Optional("mass"): Str(),
    Optional("time"): Str(),
    Optional("temperature"): Str(),
})

UnitsDef = Map({
    Optional("convention"): Enum(["SI", "CGS", "US", "dimensionless"]),
    Optional("default"): Bool(),
})

# --- 7. Mesh generation (NEW, Problem-level) ---
MeshGenerationDef = Map({
    "source": Str(),
    Optional("source_url"): Str(),
    Optional("format"): Str(),
    Optional("vertical_datums"): Map({
        Optional("horizontal"): Str(),
        Optional("vertical"): Str(),
    }),
    Optional("resolution_target"): Float(),
    Optional("filter"): Map({
        Optional("method"): Str(),
        Optional("kernel"): Str(),
        Optional("sigma"): Float(),
    }),
})

# --- 8. Roughness source (NEW, per-discretization) ---
RoughnessSourceDef = Map({
    "type": Str(),
    Optional("land_cover_data"): Map({
        "source": Str(),
        Optional("resolution"): Str(),
        Optional("url"): Str(),
    }),
    Optional("reference"): Str(),
    Optional("lookup_table"): Str(),
})

# --- 9. Updated StrongFormDef (widened unknowns, + new sections) ---
StrongFormDef = Map({
    "unknowns": Seq(UnknownDef),                      # CHANGED: was Seq(Str())
    "unknown_provenance": ProvenanceEnum,
    "equation_formulation": Str(),
    "strong_form_expression": Str(),
    "domain": Str(),
    "boundary_regions": Seq(Map({
        "name": Str(),
        "geometry": Str(),
    })),
    Optional("initial_conditions"): Seq(InitialConditionDef),     # NEW
    Optional("boundary_conditions"): Seq(BoundaryConditionDef),   # NEW
    Optional("known_analytical_solution"): Map({
        "formula": Str(),
        "reference": Str(),
    }),
    Optional("coefficients"): Map({
        Any(): Any(),   # unchanged — already accepts bare float or {value, units, std_name} map
    }),
    Optional("dimensional_check"): DimensionalCheckDef,           # NEW
})

# --- 10. Updated DiscretizationEntry (+ bmi_interface, roughness_source) ---
DiscretizationEntry = Map({
    "name": Str(),
    "provenance": ProvenanceEnum,
    "from_weak_form": Optional(Any(), default=None),
    "finite_element": Optional(FiniteElementDef),
    "discretization": Optional(Any()),
    "solver": SolverConfig,
    Optional("flux_method"): Str(),
    Optional("numerical_flux"): Str(),
    Optional("discretization_details"): Str(),
    Optional("proteus_modules"): Map({
        "physics": Str(),
        "numerics": Str(),
    }),
    Optional("bmi_interface"): Bool(),                            # NEW
    Optional("roughness_source"): RoughnessSourceDef,             # NEW
})

# --- 11. Updated top-level Problem (+ characteristic_scales, units, mesh_generation) ---
Problem = Map({
    "name": Str(),
    Optional("physical_model"): PhysicalModelDef,
    Optional("strong_form"): StrongFormDef,
    Optional("weak_forms"): Seq(WeakFormDef),
    Optional("characteristic_scales"): CharacteristicScalesDef,   # NEW
    Optional("dimensionless_numbers"): Map({Any(): DimensionlessNumberDef}),  # NEW
    Optional("units"): UnitsDef,                                  # NEW
    Optional("unit_systems"): Map({Any(): UnitSystemDef}),        # NEW
    Optional("mesh_generation"): MeshGenerationDef,                # NEW
})

# --- Top-level YDMF (unchanged structurally, Problem definition updated above) ---
YDMF = Map({
    "Problem": Problem,
    "solution_paths": DiscretizationCategory,
    Optional("vvuq"): Map({
        Optional("verification"): Any(),
        Optional("validation"): Any(),
    }),
    Optional("archive"): Any(),
    Optional("Xdmf"): Any(),
})
```

---

## 5. Cross-Report Reconciliation Notes

1. **Units string format** is consistent across `csdms-bmi-report.md` and `ydmf-units-report.md` — both use bare unit strings like `"m/s"`, `"m2/s"`, `"Pa"`, parseable by `pint`. No conflict.

2. **CSDMS `std_name` and units report's `units` field coexist** on the same `UnknownDef`/coefficient map — they are complementary, not overlapping. The units report's §4.6/§6.5 explicitly proposes validating `units` against the CSDMS-expected units for a given `std_name`. This is a **post-parse validation step**, not a schema constraint (strictyaml can't do cross-field semantic validation on its own) — implement as a `validate_std_name_units()` Python function.

3. **BMI adapter now has a concrete units source**: `csdms-bmi-report.md` §3.2's `get_var_units()` BMI method reads directly from the new `units` field on `UnknownDef`/coefficients (previously undefined — the BMI report's own example was inconsistent about this, see reconciliation note below).

4. **data-sources report's BC/IC provenance fields** (`source`, `source_url`, `temporal`) integrate directly into the new `InitialConditionDef`/`BoundaryConditionDef` — no conflict, this report essentially specified what those new sections needed to contain.

5. **Known remaining inconsistency (flagged, not yet resolved)**: `csdms-bmi-report.md` §3.4's example YDMF→BMI data flow still shows bare `unknowns: [u, p]` and `coefficients: {ν: 1.0e-5}` without units/std_name — that example predates this consolidation and should be updated to the enriched v0.2 form when the BMI report is folded into permanent docs.

6. **Multi-component `System:` coupling block** (`csdms-bmi-report.md` §3.5) is explicitly marked "illustrative example only" / future extension in that report and is **not included** in this v0.2 delta — it doesn't fit under the current `Problem`/`solution_paths` structure and needs its own design pass (likely a new top-level key, e.g. `System`, sibling to `Problem`, only present for multi-component YDMF documents). Deferred to v0.3.

---

## 6. Migration Path

1. v0.1 documents remain valid under v0.2 (only `unknowns` type is widened, not narrowed — old bare-string lists still parse).
2. New optional sections default to absent/`None` — no behavior change for documents that don't use them.
3. `normalize_unknowns()` (§2) should be called immediately after YAML parsing, before any code touches `unknowns`, so internal code never has to branch on str-vs-map.
4. Recommend bumping `archive.version` to `"ydmf-0.2"` once a document uses any v0.2-only field.

---

## 7. Deferred to v0.3 (Not in This Delta)

- Multi-component `System:` coupling block for BMI-based multi-physics (csdms-bmi §3.5) — **deferred, on hold per user request (2026-08-12)**
- ~~Automatic non-dimensionalization pipeline (`ydmf.units.non_dimensionalize()`, units report §5.3/§6)~~ — **implemented 2026-08-12**, see `src/ydmf/units.py`. Produces resolved characteristic scales, dimensionless numbers, and substitution relations for units-tagged unknowns; does *not* symbolically substitute into `strong_form_expression` (still needs the sympy/ibvp language layer) — that piece remains deferred.
- ~~`ydmf.data_sources` module for automated DEM/land-cover fetching (data-sources report §7)~~ — **partially implemented 2026-08-12**, see `src/ydmf/data_sources.py`. Manning's n land-cover lookup (incl. NLCD class codes) and a known-source URL validator are done and tested; actual raster/NetCDF fetching (SRTM/CHIRPS/ERA5/NLCD tile downloads) is **not** implemented — needs rasterio/xarray + network access, neither available in this sandbox. `fetch()` is a minimal stdlib-only URL downloader as a placeholder extension point.
- XDMF output unit-metadata embedding (units report §6.4)
