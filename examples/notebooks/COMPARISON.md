# YDMF-original.ipynb vs YDMF-updated.ipynb — executed comparison

Both notebooks were executed for real through a live Jupyter kernel
(`jupyter nbconvert --to notebook --execute --inplace ...`), not just
inspected or run as extracted scripts. Results below are from those actual
runs, on 2026-08-12.

## Bottom line

**Both notebooks run cleanly (no errors) and produce byte-identical final
XDMF output** (`poisson.xmf`). The updated notebook is a faithful,
correctness-preserving port — same worked example, same result, built on
package code instead of inline/ad hoc code.

```diff
--- original poisson.xmf output
+++ updated poisson.xmf output
(no differences)
```

## Structural differences

| | Original | Updated |
|---|---|---|
| Code cells | 8 | 3 |
| Markdown cells | 7 | 5 |
| Total stream output lines | 175 | 28 |
| Schema used | Bespoke inline `Model`/`ModelDomain`/`Domain`/`Grid`/`Topology`/`Geometry`/`DataItem` strictyaml `Map`s, defined in cell 3 of the notebook itself | `ydmf.YDMF_SCHEMA` from `src/ydmf/schema.py` — versioned, tested, shared across the whole package |
| XML-writing code | ~50 lines of `ElementTree`/`SubElement` calls + a hand-rolled `indentXML()` pretty-printer, both pasted directly into the notebook (cells 11-12) | One call: `ydmf.write_xdmf(domain, "poisson.xmf")` — implementation lives in `src/ydmf/xdmf.py`, unit-tested independently of any notebook |
| Flow-style YAML handling | Called `strictyaml.dirty_load(..., allow_flow_style=True)` directly in the notebook, with a comment explaining why plain `load()` wouldn't work | Same underlying need, but handled once inside `ydmf.load_ydmf`/`ydmf.validate_ydmf` — every caller gets it automatically, no per-notebook boilerplate or tribal knowledge required |
| What the schema models | One schema conflates the PDE definition (`Model.Equation`, `Model.BoundaryConditions`) with the mesh archive (`Domain.TimeCollection`) and even includes an unused `ModelDomain.Grid` (single hexahedron, never fed into the XDMF conversion — dead code in the original) | Two clear, separate concerns: a `Problem`/`solution_paths` document (physics + weak form + discretization + provenance) validated against the real package schema, and a plain `Archive.Domain` dict for the mesh time-series, matching `docs/ydmf-schema.md` §5's stated design ("XDMF is a consumer of the archive, not a schema dependency of YDMF") |
| Redundant re-print of the YAML | Cell 6 prints raw `poisson_data`; cell 9 prints `poisson_yaml.as_yaml()` — same content, printed twice, once before and once after strictyaml round-trips it (mostly to demonstrate strictyaml's re-serialization, e.g. dropping trailing whitespace/reformatting flow lists) | Not needed — `validate_ydmf()` is called once and the result is used directly; no round-trip demonstration since the schema/loader behavior is already covered by the package's own test suite (`tests/test_schema.py`) |

## What's genuinely different (not just refactored)

1. **Original models a toy `ModelDomain.Grid`** (a single unit-cube hexahedron with 8 hardcoded nodes) **that is never actually used** by the XML-conversion code in cells 11-12 — only `Domain.TimeCollection` gets converted. This is either leftover scaffolding or a demonstration of the schema accepting it without ever exercising the conversion path. The updated notebook doesn't carry this dead weight forward.

2. **Original's `Model.Equation` is a free-text string** (`"Δu=f\nΔv=g"` — note: two equations for a problem that's actually about a *single* scalar unknown `u`, likely a copy/paste leftover from a coupled-system template) with no connection to weak forms, discretization, or solver choice. The updated notebook's `Problem.strong_form` + `weak_forms` + `solution_paths.discretizations` sections make that connection explicit and structured (equation → weak form → FE space → solver), matching how the rest of the YDMF schema/package actually works.

3. **Original has no discretization/solver metadata at all** — just the raw mesh. The updated notebook's `Problem` document declares `P1_tetrahedral`/`CG`/order 1/`petsc`/`proteus_derived` provenance, which is exactly the kind of information `docs/ydmf-schema.md` was designed to carry and that downstream tooling (BMI adapters, VVUQ, non-dimensionalization) in this package now consumes.

## Verification method

```bash
jupyter nbconvert --to notebook --execute --inplace YDMF-original.ipynb
jupyter nbconvert --to notebook --execute --inplace YDMF-updated.ipynb
```

Both completed with exit code 0, no `error` output cells. Final XDMF text
output (last code cell of each) compared with a plain string diff — no
differences found.
