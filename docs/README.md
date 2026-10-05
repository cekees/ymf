# YMF documentation

Start with the two guides. They describe what exists today. The rest are
design drafts and background research, kept for their reasoning.

## Guides (current)

| Document | What it covers |
|---|---|
| [problem-specification.md](problem-specification.md) | Writing and validating a problem spec: structure, branching, provenance, units, known gaps |
| [archive-format.md](archive-format.md) | The archive a solver writes: files on disk, the domain model, parallel output, XDMF conversion, use in Proteus |

## Design drafts

| Document | Status |
|---|---|
| [ymf-schema.md](ymf-schema.md) | The original v0.1 schema draft. Still the best statement of the design principles and of the intended `vvuq` shape. Its §4 "Archive Configuration" and §5 "XDMF Preservation" were superseded by the implemented [archive format](archive-format.md), and the title expands the old working name. |
| [ymf-schema-v0.2-delta.md](ymf-schema-v0.2-delta.md) | The v0.2 changes on top of v0.1, which `ymf.schema` implements. §7 lists what was deferred to v0.3. |

## Background research

These were drafted with an LLM research agent, and have not all been
checked against their sources. Treat specific numbers and API names in
them as leads to verify, not facts.

| Document | Covers | Caution |
|---|---|---|
| [ymf-units-report.md](ymf-units-report.md) | Unit libraries, dimensional analysis, non-dimensionalization | The basis for `ymf.units`. |
| [data-sources-report.md](data-sources-report.md) | Topography, bathymetry, land cover, roughness and forcing datasets | The basis for `ymf.data_sources`; dataset details not re-verified. |
| [csdms-bmi-report.md](csdms-bmi-report.md) | CSDMS Standard Names and the Basic Model Interface | **Contains known errors.** CSDMS is the *Community Surface Dynamics Modeling System*, not what §1 says. The registry fetched at the time listed 2,653 names, not "~4,000+". The real naming grammar is `object__quantity` (double underscore), not the pattern the report gives. Its "Registry API" section and example variable names were never confirmed against a source. |

## History

[session-summary-2026-07-29.md](session-summary-2026-07-29.md) is the notes
from the session where YMF was designed, under its earlier name YDMF. It
records the original intent, including ideas that have not been built (the
`ibvp`/sympy symbolic layer, multi-backend dispatch).
