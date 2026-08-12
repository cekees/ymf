# YDMF — YAML Data Model and Format for Physics

A structured, human- and LLM-editable YAML format for describing physical
models, numerical discretizations, and VVUQ (verification, validation,
uncertainty quantification) results, with support for branching across
formulations, solution paths, and solver backends.

## Status

Early design/implementation phase (schema v0.2 draft). Not yet released.

## Design Principles

- **Composability** — each branch adds discretizations or weak forms without
  modifying existing entries.
- **Provenance** — every block records how it was generated (`llm_derived`,
  `human_specified`, `proteus_derived`, `human_edited`).
- **Two-tier referencing** — integer index for programmatic use, string label
  for human/LLM readability; both resolve to the same array position.
- **XDMF as a consumer, not a dependency** — archive output maps cleanly to
  XDMF for tools like ParaView, without embedding XDMF as a schema
  requirement.
- **Backward compatible evolution** — schema changes are additive wherever
  possible; the one breaking change so far (`unknowns` widened from bare
  strings to optional structured entries) ships with a normalization helper.

## Repository Layout

```text
ydmf/
  docs/                      Design docs and background reports
    ydmf-schema.md           v0.1 schema draft (original)
    ydmf-schema-v0.2-delta.md  v0.2 consolidation delta (current)
    csdms-bmi-report.md      CSDMS Standard Names + BMI integration report
    ydmf-units-report.md     Physical units / dimensional analysis report
    data-sources-report.md   Geospatial/environmental data source catalog
    session-summary-2026-07-29.md  Original design session notes
  src/ydmf/                  Python package (schema validation, helpers)
  tests/                     pytest test suite
```

## Pipeline (conceptual)

```
Stage 0: Physical processes   → choose which to model
Stage 1: Equation formulation → strong form PDEs
Stage 2: Solution variables   → unknowns (u,p vs ρu,p etc.)
Stage 3: Weak form            → theoretical reformulation (C² → H1)
Stage 4: Discretization       → FEM/FD/FV, element type, order, mesh
Stage 5: Solver                → linear/nonlinear, time stepping, tolerances
```

An LLM assists at stages 0-2 (model selection, formulation, unknown choice)
and generates boilerplate for stages 3-5; a human reviews and approves at
each step. YDMF captures every stage as a structured, verifiable artifact.

## Installation (development)

```bash
pip install -e ".[dev]"
pytest
```

## Quick Example

```python
from ydmf import load_ydmf, normalize_unknowns

doc = load_ydmf("examples/poisson.yaml")
unknowns = normalize_unknowns(doc["Problem"]["strong_form"]["unknowns"])
```

## Units and Dimensional Analysis

`ydmf.units` (pint-backed) provides:

- `check(doc)` — validates dimensional consistency between `unknowns`,
  `boundary_conditions`, and `initial_conditions` unit strings. Documents
  with no units at all pass trivially (units are optional everywhere).
- `characteristic_scale_values(problem)` / `compute_dimensionless_numbers(problem)`
  — resolve `characteristic_scales` (incl. `derived: true` formulas) and
  evaluate `dimensionless_numbers` formulas (e.g. Reynolds number).
- `non_dimensionalize(doc)` — resolves characteristic scales/dimensionless
  numbers and derives substitution relations (`u* = u / U_char`, etc.) for
  units-tagged unknowns. Does not symbolically rewrite
  `strong_form_expression` — that needs the sympy/ibvp language layer
  (still deferred, see roadmap).

Unit strings accept the notational variants used throughout the docs
(`m2/s`, `m²/s`, `m^2/s`, `kg/(m·s)`) via `normalize_unit_string()`.

## Data Sources

`ydmf.data_sources` provides:

- `resolve_roughness_source(roughness_source, land_cover)` /
  `manning_n_for_land_cover(...)` — Manning's n lookup from land cover
  class (internal keys or raw USGS NLCD class codes), per
  Chow (1959) / Arcement & Schneider (1989) reference ranges. Supports a
  CSV `lookup_table` override.
- `validate_source_record(record)` — checks a `source`/`source_url` pair
  against a small built-in catalog (SRTM, CHIRPS, ERA5, NLCD, etc.) and
  returns warnings (never errors — sources outside the catalog are fine).
- `fetch(url, dest)` — minimal stdlib-only downloader; a placeholder
  extension point, not a real geospatial fetch pipeline (no
  rasterio/xarray dependency assumed).

## Roadmap

See `docs/ydmf-schema-v0.2-delta.md` §7 for what's deferred to v0.3:
multi-component BMI `System:` coupling block (on hold per user request),
symbolic non-dimensionalization of `strong_form_expression`, automated
raster/NetCDF data fetching, XDMF unit-metadata embedding.

## License

MIT (see `LICENSE`).
