# YMF documentation

Start with the two guides. They describe what exists today. Design notes say
what was decided and whether it is built. History and research are kept for
their reasoning, not as a description of the code.

## Guides (current)

| Document | What it covers |
|---|---|
| [problem-specification.md](problem-specification.md) | Writing and validating a problem spec: structure, branching, composition, provenance, units, known gaps |
| [archive-format.md](archive-format.md) | Archives: the archive of a spec (outputs keyed by their input, studies, reruns, `extends`), and the archive a solver writes: files on disk, the domain model, parallel output, XDMF conversion, use in Proteus |
| [xdmf-model.yaml](xdmf-model.yaml) | The archive model laid out in the order of XDMF's own [Model and Format](https://xdmf.org/index.php/XDMF_Model_and_Format.html) page: each XDMF element and attribute, the `.ymf` key that holds it, and whether it is supported. For readers who know XDMF. Checked by `tests/test_xdmf_model.py`. |
| [changes.md](changes.md) | Gaps that have been closed, dated, with the commit that closed each. The guides describe the result; this says why it is that way. |

## Design notes

A design note is written before the code, records the decisions taken on it,
and lands in `design/` with the code that implements it.

| Note | Status |
|---|---|
| [design/time-outside-space.md](design/time-outside-space.md) | **Implemented** 2026-10-26: archive format 2, one list of instants with the spaces inside each. |

## History

Superseded drafts, kept for the reasoning behind decisions still in force.
Each starts with a note saying what replaced it.

| Document | What it was |
|---|---|
| [history/schema-v0.1.md](history/schema-v0.1.md) | The original v0.1 schema draft. Its § 3 is still the only description of the intended `vvuq` shape. |
| [history/schema-v0.2-delta.md](history/schema-v0.2-delta.md) | The v0.2 changes on top of v0.1, which `ymf.schema` implements. |
| [history/session-2026-07-29.md](history/session-2026-07-29.md) | Notes from the session where YMF was designed, under its earlier name YDMF. |

## Background research

Drafted with an LLM research agent, and not all checked against their
sources. Treat specific numbers and API names in them as leads to verify,
not facts. Each says so in its own header.

| Document | Covers | Caution |
|---|---|---|
| [research/units.md](research/units.md) | Unit libraries, dimensional analysis, non-dimensionalization | The basis for `ymf.units`. |
| [research/data-sources.md](research/data-sources.md) | Topography, bathymetry, land cover, roughness and forcing datasets | The basis for `ymf.data_sources`; dataset details not re-verified. |
| [research/csdms-bmi.md](research/csdms-bmi.md) | CSDMS Standard Names and the Basic Model Interface | **Contains known errors**, listed in its header. |
