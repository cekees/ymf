# Design note: time outside, space inside

**Status: proposed, for cekees's approval. No code changes accompany this note.**
Closes the open question at the end of [archive-format.md](archive-format.md)
(recorded 2026-10-23). Drafted 2026-10-25.

Nothing here touches the Proteus writer until the decisions in
[§ Decisions needed](#decisions-needed) are answered.

## 1. The proposal

Today a domain is a list of **time collections, one per finite-element space**,
each holding its own list of instants. The proposal transposes it: **one list
of instants, with the spaces inside each instant.**

```yaml
# format 1 (today)                       # format 2 (proposed)
time_collections:                         times:
- name: Mesh_c0p1_Lagrange                - time: 0.0
  data: [{time: 0.0, ...},                  spaces:
         {time: 0.5, ...}]                    Mesh_c0p1_Lagrange: {topology: ..., geometry: ..., attributes: ...}
- name: Mesh_c0p2_Lagrange                    Mesh_c0p2_Lagrange: {...}
  data: [{time: 0.0, ...},                - time: 0.5
         {time: 0.5, ...}]                  spaces: {...}
```

A space's value is exactly what a step is today, minus `time`: a uniform grid
`{topology, geometry, attributes}` or a spatial one `{spatial_collection: [...]}`.
Everything below the space level (grids, topology, geometry, attributes,
data items, inline values, canonicalization, validation) is **unchanged**.

## 2. Why it is worth a format version

Measured against the code on `main` (ymf `ef44cbe`, proteus `9a3f1160`):

1. **The writer already requires it.** `allGatherIncremental` raises
   `RuntimeError` if any collection lacks a grid for the instant being
   archived ("every collection is written at every archived instant").
   The only thing that varies is *which* spaces exist: a model that first
   archives at step 40 adds its space then. The proposal encodes exactly that
   (§ 4.3). Format 1 spends a whole second level of indexing on a
   redundancy the writer enforces anyway.
2. **An instant becomes the unit of atomicity.** Proteus's live append
   (`_appendLiveStep`) writes one record per collection per step. A run killed
   between two of those leaves collection A at step N and collection B at
   step N-1, and the reader keeps both. With one record per instant a kill
   loses a whole instant or none of it. This is a correctness gain, not
   tidiness.
3. **"Which instants does this archive hold?"** is `[t["time"] for t in
   domain["times"]]`. Today `Archiver.archived_times` walks the first
   collection and has a comment that "every collection covers the same
   instants, so one is enough", a fact the format does not state.
4. **`append_step` stops needing collection names** (§ 4.2), which is the
   thing the lazily created collections forced.

The objection (operator splitting: no instant at which every space is
current) is answered in archive-format.md: the archive samples at
synchronization points, and the archiving step is a choice separate from the
numerical step. It also never applied to what Proteus does today, per item 1.
A space that is genuinely absent at an instant is simply absent from `spaces`.

## 3. Two version numbers, and which of them moves

"v1/v2" is ambiguous in this codebase. There are two independent numbers:

| number | where | values today | meaning |
|---|---|---|---|
| `ymf_archive_version` | top of a `.ymf` document (`ymf.archive.ARCHIVE_FORMAT_VERSION`) | 1 | shape of the `domain` tree |
| `METADATA_FORMAT_VERSION` | HDF5 root attribute `ymf_archive_metadata_version` (`proteus.Archiver.AR_base`) | 1 = XDMF `<Grid>` fragments (Proteus <= 1.9.x), 2 = YAML grid dicts | how per-step, per-rank grid metadata is stashed in the `.h5` |

**Proposal: `ymf_archive_version` goes 1 -> 2; `METADATA_FORMAT_VERSION`
stays at 2.** The HDF5 stash is datasets named `<collection>_<step>` and the
collection names live in `ymf_archive_collections`. That is a storage layout
for the *writer's* working state, and the transposition happens when the
document is assembled from it (`assemble_domain`) and when live records are
appended. Leaving it alone means:

- no change to the collective MPI path (`allGatherIncremental`, the
  `allreduce(MAX)` width agreement, the per-rank gather);
- hot start (`_collection_names`, `_seed_domain`) is untouched;
- every existing Proteus 2.0-line archive's `.h5` assembles into a format 2
  document with no migration, and `scripts/gatherTimes` can regenerate the
  `.ymf` from any `.h5` (it already does).

The cost is that the `.h5` keeps `n_spaces x n_steps` small metadata
datasets. That is the existing cost and can be revisited separately.

The closure document (`ymf: 0.2.x`, an output's `approximation`) carries the
domain too. It has no archive-version key; it is guarded by the package
version in `ymf:` (`check_version` refuses a newer writer). A format 2
`approximation` is written by a newer ymf, so an old reader refuses it
instead of misreading it. See § 7 for the version number this implies.

## 4. Changes, by component

### 4.1 Data model (`ymf.archive`)

New canonical in-memory shape: `domain = {"times": [instant, ...]}` with
`instant = {"time": float, "spaces": {name: grid_or_spatial}}`.

- `new_domain()` -> `{"times": []}`. The `name` argument goes (a domain no
  longer starts with a named collection). Callers that passed a name keep
  working for one release with a `DeprecationWarning` (§ 7).
- `add_instant(domain, time, spaces)`, with `spaces` a mapping built from
  `uniform_grid(...)` / `spatial_grid(...)` (thin wrappers over the existing
  `grid()` and `{"spatial_collection": ...}`).
- `add_uniform_step` / `add_spatial_step` / `add_collection` keep their
  signatures and map onto the new shape: they find the instant whose time
  equals `time` (or append one) and set `spaces[collection]`. They exist so
  the Proteus writer, `examples/write_archive.py` and the tests can migrate in
  steps. `collection=` as an integer index is only defined while a domain was
  built through these shims; a name is the supported form.
- **Space order is significant** and is the order of the first instant that
  contains the space, then insertion order. `pyyaml`'s dumper already runs
  `sort_keys=False`. Hot start needs the creation order (`EBMesh` would sort
  ahead of `Mesh`), so the format *states* it: the key order of `spaces` is
  document order.
- Two views for consumers that think per space (the `MeshTools` mesh reader,
  `scripts/clearh5.py`, `ymf2xmf --show`):
  `space_names(domain)` and `series(domain, name)` -> `[(time, grid), ...]`.
  `findArchiveGrid(domain, MeshTag, id_in_collection)` becomes a three-line
  use of `series`.
- `canonicalize_domain` / `validate_domain` gain: `times` is a list,
  `time` is a float, instants are not required to be strictly increasing (a
  hot-started run may repeat a time; see § 8), `spaces` is a non-empty
  mapping, and the existing grid checks run per space with the path
  `times[i].spaces[name]...` in errors.
- `select_domain(domain, collections=, fields=, steps=)`: `collections`
  filters `spaces`; `fields` unchanged; **`steps` now indexes instants, not
  per-collection steps**. This is a behavior change only when collections have
  different step counts, which the writer forbids except for lazily added
  spaces. Keep the keyword name `collections=` (it is the CLI's `-c`).

### 4.2 `append_step` without collection names

```python
append_step(path, instant, *, fsync=False)     # instant = {"time": t, "spaces": {...}}
```

One YAML document per instant, same `---` ... `...` framing, same
truncation rule (an instant is complete or absent). The record is
`{"step": instant}`; the `"collection"` key disappears. `read_ymf` folds a
record by appending it to `times`, creating nothing by name because the
spaces are *in* the record. A space that first appears at step 40 is just a
key that first appears in the 41st record.

The old three-argument form `append_step(path, collection, step)` is kept for
format 1 writers during the deprecation window and writes a **format 1**
record; a file is never mixed (§ 4.4 `read_ymf` refuses it).

### 4.3 Writer (Proteus `Archiver.py`) — after approval

Territory note: this is the only proteus change, and it is contained.

- `allGatherIncremental` keeps its loop over collections and its HDF5 stash
  unchanged. It builds `live_steps` as today, then `_appendLiveStep` groups
  them into **one** `{"time": t, "spaces": {name: grid_or_spatial}}` record
  and makes **one** `append_step` call per instant instead of one per
  collection. Its header writes `new_domain()` (empty `times`).
- `assemble_domain` iterates steps outermost and collections innermost,
  building each instant from the datasets that exist at that step. Its
  `n_steps is None` end-detection ("count up until a step is missing") becomes
  "until no collection has the step". The v1-metadata path (`_v1_grids`)
  feeds the same builder.
- `_inlineStep` / `_inlineDomain` build instants directly.
- `ArchiveCollection` / `ArchiveDomain` (the *running*, one-step-deep objects
  that the ~25 `write*Xdmf_*` writers use) are **not changed**. They are a
  writer-side convenience and never see the document's shape.
- `archived_times` becomes `[i["time"] for i in domain["times"]]`.
- Readers in proteus that index `time_collections`: `MeshTools.findArchiveGrid`
  and `readArchiveDomain` callers (`MeshTools.py:7003, 7138`),
  `scripts/clearh5.py:44`, `scripts/gatherTimes:43`,
  `proteus/cli/ymfrun.py` (via `domain_arrays`, which walks the tree and
  needs only the new walk), and the tests listed in § 6.

### 4.4 Reader (`read_ymf`)

- `ymf_archive_version` 1 **or** 2 is accepted; anything else is refused as
  today. A version 1 document is **upgraded in memory** by
  `transpose_collections(domain)` (§ 5) and returned in the format 2 shape.
  Callers never see `time_collections` again.
- Appended records: a stream whose header is format 1 may contain
  `{"collection", "step"}` records (an in-flight run from before the
  upgrade); they are folded into the header's collections *then* the whole
  domain is transposed. A format 2 header with a format 1 record, or the
  reverse, is an error naming the line.
- The truncation rule is unchanged.

### 4.5 `ymf2xmf` and `ymf.xdmf`

XDMF's model is the old one (a Temporal `GridCollection` of grids of one
type), and every viewer expects it. The conversion therefore **transposes
back**:

- `build_xdmf_tree`: for each space name (in document order), one
  `Grid CollectionType=Temporal` holding that space's instants, exactly the
  tree emitted today. The XML is byte-identical for an archive that round
  trips, which makes the change testable: convert a format 1 archive and its
  format 2 upgrade, diff the `.xmf`.
- `parse_xdmf_domain`: reads per-space collections, then
  `transpose_collections` joins them on time (§ 5). `read_xdmf` and
  `round_trip_equal` keep their contracts.
- `ymf.xdmf.collection_names` (used by Proteus hot start for <= 1.9
  archives) is about the `.xmf` and does not change.
- CLI: `--show` prints one line per space as now (`series()` provides it)
  plus an instants line (`6 instants, t=0..50, 2 spaces`). `-s/--steps`
  selects instants. `-c/--collection` selects spaces. The `-i` prompt is
  unchanged in appearance. No flag is renamed.

I considered emitting XDMF as a Temporal collection of Spatial collections
(one per instant, spaces as the members). Viewers would then treat the
linear and quadratic meshes at an instant as one multi-block dataset, which is
what ParaView does with per-rank subdomain grids, but the two meshes overlap
cell for cell and are not pieces of one mesh. Keeping them as separate
datasets is the current behavior and the right one.

## 5. Transposing format 1: the migration

One function, used by `read_ymf`, `parse_xdmf_domain`, and the Proteus
v1-metadata path:

```
transpose_collections(collections) -> times
```

1. Collections are taken in list order; that order becomes the order of
   `spaces`.
2. Each collection's steps are keyed by `(time, occurrence)`, where
   `occurrence` is how many earlier steps of the **same collection** had that
   exact `time`. This pairs repeats correctly if a time legitimately recurs.
3. Instants are the union of those keys, ordered by first appearance when the
   collections are walked step-major (step 0 of every collection, then
   step 1, ...). When all collections have equal-length, equal-time series
   (the writer's normal case) this is just "zip".
4. A space missing at a key is omitted from that instant's `spaces`.
5. Times are compared exactly. Proteus stores them with `repr(float)` and
   YAML round trips a float exactly, so equality is the right test; a
   tolerance would merge instants that are really distinct.

A collection with a single step (a Uniform or Spatial top-level grid read
from foreign XDMF) is one key at its time (0.0 if it has none), as
`parse_xdmf_domain` does today.

**Archives on disk need no rewrite.** Format 1 files stay readable
indefinitely. `ymf2xmf --upgrade run.ymf` (new, optional) rewrites a format 1
`.ymf` as format 2, and `scripts/gatherTimes` regenerates the document from
the `.h5` in whatever format the installed ymf writes. The `.h5` is
format-agnostic (§ 3).

Cases and where each is covered:

| archive | handled by |
|---|---|
| Proteus <= 1.9.x: metadata v1 (XDMF fragments), no `.ymf` | unchanged read path; `assemble_domain` builds instants |
| Proteus 2.0-line, metadata v2, `.ymf` format 1 | `read_ymf` upgrades in memory |
| In-flight / killed run, format 1 header + per-collection records | `read_ymf` folds then transposes (§ 4.4) |
| Closure `approximation: {time_collections: [...]}` | the same `transpose_collections`, applied on load; written as `times` |
| Foreign `.xmf` | `parse_xdmf_domain` |

## 6. Test plan

New tests, then the existing suites that mention the shape.

- `transpose_collections`: zip case, lazily added space, repeated time,
  single-step collection, empty collection, order preservation.
- Round trip: format 1 fixture -> read -> write (format 2) -> read equals;
  format 1 -> `ymf2xmf` and format 2 -> `ymf2xmf` give identical `.xmf`.
- `append_step`: a killed-mid-instant file drops the whole instant; a space
  appearing at record 41; format mismatch errors.
- `select_domain` with `steps=-1` and unequal step counts.
- The shipped examples and `tests/fixtures` in both formats.
- Proteus (after approval): `test_archive_dtypes`, `test_archive_no_xml`,
  `test_adr_problem` archive checks, `test/ci/test_systemStep.py`
  (`archived_times`), a hot-start test from a format 1 archive, and a
  two-space run (c0p1 + c0p2) at 1 and 2 ranks, global and per-subdomain,
  comparing `.xmf` from `ymf2xmf` to the pre-change output.
- Because proteus CI installs ymf from `@main`, the ymf PR merges first;
  until the Proteus PR lands the new ymf must still serve the old Proteus
  writer, which is what the shims in § 4.1 are for.

## 7. Compatibility and release order

- ymf: the new reader accepts both formats; the new writer writes format 2.
  `ymf.__version__` goes to the next minor (0.3.0 proposed) since the closure
  `approximation` shape changes and `check_version` is what protects old
  readers.
- Shims (`new_domain(name)`, `add_collection`, `add_uniform_step(...,
  collection=)`, three-argument `append_step`) keep the current Proteus
  writer working against the new ymf unchanged, but with the format 1 record
  for `append_step`. One release later they warn; one after that they go.
- Order: (1) ymf PR implementing § 4.1, 4.2, 4.4, 4.5 and § 5, merged;
  (2) Proteus PR changing `Archiver.py` (§ 4.3) and its readers, merged;
  (3) shims deprecated. Between 1 and 2 nothing in Proteus changes behavior.
- External: any script that does `domain["time_collections"]` directly breaks
  on format 2. In this tree that is `scripts/clearh5.py`, `scripts/gatherTimes`,
  `MeshTools.py`. I have not searched user code outside these repositories;
  `air-water-vv` is the one worth checking before release.

## 8. Things I could not verify

- **Hot start repeating a time.** Whether a hot-started run can write an
  instant at a time already in the archive. The occurrence index in § 5
  handles it; whether the writer should *replace* or *append* is a Proteus
  question this note does not answer.
- **Lazily added spaces in a real run.** The writer's own check forbids a
  collection with no grid at an instant, but I found no case in the test suite
  that adds a space mid-run. The docstring of `append_step` says it happens.
  The design handles it; no fixture exercises it yet.
- **`approximation` consumers outside these repos.**
- **Cost.** I have not measured parse time of a transposed document against
  the current one. The record count drops from `n_spaces x n_steps` to
  `n_steps`, so the append path should get cheaper; the one-pass transpose on
  load is O(total steps).

## Decisions needed

1. **Approve the transposition** and `ymf_archive_version` 1 -> 2.
2. **Key name**: `times` / `spaces` as drawn, or something else
   (`instants`)? Rename is free until code exists.
3. **`spaces` as a mapping** (lookup by name; order is document order, which
   the format must state) **or a list** of `{name, ...}` like `attributes`
   (order is structural; one extra key per entry). Recommended: mapping.
4. **Leave `METADATA_FORMAT_VERSION` at 2** and transpose only at assembly
   (§ 3). Recommended: yes.
5. **Release number** for ymf (0.3.0?) and whether the deprecation window is
   one release or two.
6. **Optional `ymf2xmf --upgrade`**: wanted, or is "format 1 reads forever"
   enough?
