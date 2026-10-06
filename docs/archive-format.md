# The YMF archive format

A YMF archive is how a simulation records its output: the heavy arrays go
in HDF5, and a small YAML document (`.ymf`) says what each array is — which
mesh, which time, which field, what shape and type. This half of YMF is
implemented, tested against real solver output, and is what Proteus writes
on its `ymf_io` branch.

The code is `ymf.archive` (data model and YAML I/O), `ymf.xdmf` (conversion
to XDMF), and the `ymf2xmf` command. All three depend on `pyyaml` alone.

If you know XDMF, [xdmf-model.yaml](xdmf-model.yaml) goes through XDMF's
model element by element, in the order of its specification, and says
which parts a `.ymf` archive represents.

- [Files on disk](#files-on-disk)
- [The `.ymf` document](#the-ymf-document)
- [The domain model](#the-domain-model)
- [Writing an archive](#writing-an-archive)
- [Self-contained archives: inline values](#self-contained-archives-inline-values)
- [Parallel output](#parallel-output)
- [What the writer checks, and what it doesn't](#what-the-writer-checks-and-what-it-doesnt)
- [Viewing: converting to XDMF](#viewing-converting-to-xdmf)
- [Why pyyaml and not strictyaml](#why-pyyaml-and-not-strictyaml)
- [Use in Proteus](#use-in-proteus)
- [Limits](#limits)

## Files on disk

```text
run.h5     heavy data: node coordinates, connectivity, field values
run.ymf    the archive of record: YAML metadata pointing into run.h5
run.xmf    optional, derived from run.ymf by `ymf2xmf`, for ParaView/VisIt
```

The `.ymf` refers to datasets as `run.h5:/T_t3`, with the file path
relative to the `.ymf` itself, so the three files move together as a
directory. The `.xmf` refers to the same datasets; it holds no data of its
own.

An archive can also hold its arrays **inline**, with no HDF5 file at all;
see [Self-contained archives](#self-contained-archives-inline-values).

The `.xmf` is a *view* of the archive. A solver writes `.ymf` and nothing
else, and the `.xmf` is produced on demand. A solver that writes both
formats has to keep two representations consistent, and eventually won't.

## The `.ymf` document

Three top-level keys:

```yaml
ymf_archive_version: 1        # refused on read if it isn't the version this ymf knows
domain:                       # the mesh/time/field tree -- see below
  TimeCollections: [...]
extra:                        # optional: anything with no XDMF equivalent
  Problem: {...}              #   typically the problem spec that produced the run
  solution_paths: {...}
```

`extra` is free-form, and is where an archive records its provenance. The
example archives carry the whole `examples/heat_equation.ymf` problem
specification, so the output says what problem it solved. The archive
writer doesn't interpret `extra`; it stores it.

Here is one time step of `heat.ymf` from
[`examples/write_archive.py`](../examples/write_archive.py):

```yaml
ymf_archive_version: 1
domain:
  TimeCollections:
  - Name: Mesh Spatial_Domain
    Data:
    - Time: 0.0
      Topology:
        Type: Triangle
        NumberOfElements: 512
        DataItem: {Format: HDF, DataType: Int, Precision: 4,
                   Dimensions: [512, 3], Data: 'heat.h5:/elements'}
      Geometry:
        Type: XYZ
        DataItem: {Format: HDF, DataType: Float, Precision: 8,
                   Dimensions: [289, 3], Data: 'heat.h5:/nodes'}
      Attributes:
      - Name: T
        AttributeType: Scalar
        Center: Node
        DataItem: {Format: HDF, DataType: Float, Precision: 8,
                   Dimensions: [289], Data: 'heat.h5:/T_t0'}
      - Name: q
        AttributeType: Vector
        Center: Node
        DataItem: {Format: HDF, DataType: Float, Precision: 8,
                   Dimensions: [289, 3], Data: 'heat.h5:/q_t0'}
    - Time: 10.0
      ...
extra:
  Problem:
    name: Transient heat conduction, unit square
    ...
```

(The file on disk is block-style YAML throughout; flow style is used here
only to keep the listing short.)

## The domain model

A domain is plain dicts and lists — no classes — so it costs nothing to
build, gathers across MPI natively, and serializes directly. The names
follow XDMF on purpose, so the conversion is one-to-one.

```text
domain
└── TimeCollections: [collection, ...]       one per finite-element space
    └── collection: {Name, Data: [step, ...]}
        └── step: one instant, in one of two forms
            ├── uniform:  {Time, Topology, Geometry, Attributes}
            └── spatial:  {Time, SpatialCollection: [grid, ...]}   one grid per subdomain
                          grid = {Name?, Topology, Geometry, Attributes}

Topology   {Type, NumberOfElements, NodesPerElement?, DataItem}
Geometry   {Type, DataItem}
Attribute  {Name, AttributeType, Center, DataItem}
DataItem   {Format: HDF, DataType, Precision, Dimensions, Data: "file.h5:/dataset"}
         | {Format: XML, DataType, Precision, Dimensions, Include: "sidecar.txt"}
         | {Format: XML, DataType, Precision, Dimensions, Values: [1.0, 2.5, ...]}
```

**Why a list of time collections.** A solver writes one collection per
finite-element space: a Proteus run with quadratic velocity holds both
`Mesh Spatial_Domain` (the linear base mesh) and `Mesh_c0p2_Lagrange` (same
elements, more nodes). 195 of the 814 Proteus archives surveyed when this
model was designed had two. `{"TimeCollection": {...}}` (singular) is
accepted as input shorthand and normalized to the list.

**Closed sets.** These fields take a fixed set of values, and anything else
is rejected when the archive is written:

| Field | Legal values |
|---|---|
| `Topology.Type` | `Polyvertex Polyline Polygon Triangle Quadrilateral Tetrahedron Pyramid Wedge Hexahedron`, the quadratic `Edge_3 Tri_6 Triangle_6 Quad_8 Quadrilateral_8 Tet_10 Tetrahedron_10 Pyramid_13 Wedge_15 Hex_20 Hexahedron_20`, `Mixed`, and the structured `2DSMesh 2DRectMesh 2DCoRectMesh 3DSMesh 3DRectMesh 3DCoRectMesh` |
| `Geometry.Type` | `XYZ XY X_Y_Z X_Y VXVYVZ ORIGIN_DXDYDZ ORIGIN_DXDY` |
| `AttributeType` | `Scalar Vector Tensor Tensor6 Matrix GlobalID` |
| `Center` | `Node Cell Face Edge Grid` |
| `DataType` | `Float Int UInt Char UChar` |
| `Format` | `HDF XML Binary` |

**Defaults.** `canonicalize_domain()` fills every defaultable field
explicitly, and the file on disk is always canonical. The one default worth
knowing is that XDMF's default `Precision` is **4** bytes, which is wrong
for `float64`. Use `data_item_for()`, which reads the precision off the
array.

## Writing an archive

The constructor helpers build the dicts and apply defaults consistently.
`data_item_for()` is the one to reach for: it takes the array and derives
`Dimensions`, `DataType` and `Precision` from it, so the metadata can't
disagree with what's in the file.

```python
import h5py, numpy as np
from ymf.archive import (new_domain, add_uniform_step, topology, geometry,
                         attribute, data_item_for, write_ymf)

domain = new_domain("Mesh Spatial_Domain")
with h5py.File("heat.h5", "w") as h5:
    h5["nodes"], h5["elements"] = nodes, elements          # mesh written once
    topo = topology("Triangle", len(elements), data_item_for(elements, "heat.h5:/elements"))
    geom = geometry(data_item_for(nodes, "heat.h5:/nodes"))
    for k, t in enumerate(times):
        h5["T_t%d" % k] = T = solve_to(t)
        add_uniform_step(domain, t, topo, geom, [
            attribute("T", data_item_for(T, "heat.h5:/T_t%d" % k)),
        ])

write_ymf(domain, "heat.ymf", extra={"Problem": problem})
```

[`examples/write_archive.py`](../examples/write_archive.py) is the
complete, runnable version, including a vector field, a cell-centred field,
and a per-subdomain archive. [`examples/read_archive.py`](../examples/read_archive.py)
reads one back.

Reading is `domain, extra = read_ymf("heat.ymf")`. Each DataItem's `Data`
splits on its last `:` into a file (relative to the `.ymf`) and a dataset
path.

## Self-contained archives: inline values

A DataItem can hold its values itself instead of pointing at them, as
XDMF's XML format does. The `.ymf` and the `.xmf` derived from it are then
each complete on their own: no HDF5 file, no sidecar, and no HDF5
dependency for the solver that writes them. This suits small problems,
test cases, and solvers that would rather not depend on HDF5. For large
arrays, HDF5 remains the right choice; inline values are text.

```python
topo = topology("Triangle", len(cells), data_item_for(cells, inline=True))
geom = geometry(data_item_for(nodes, inline=True))
# or, without an array library:
geom = geometry(data_item([4, 3], values=[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, ...], precision=8))
```

```yaml
DataItem:
  Format: XML
  DataType: Int
  Precision: 4
  Dimensions: [2, 3]
  Values: [0, 1, 3, 0, 3, 2]
```

`Values` is flat, in row-major order, as in XDMF's element text; nested
lists and arrays are flattened on the way in. The values are coerced to
the declared `DataType` (a non-integer under `Int` is refused) and their
count is checked against `Dimensions`. `ymf2xmf` writes them as the
DataItem's text, one row per line, with floats written exactly, so a
`.ymf` → `.xmf` → `.ymf` round trip loses nothing. Inline and referenced
DataItems can be mixed freely in one archive.
`python examples/write_archive.py --inline` writes the examples this way,
and runs without h5py installed.

## Parallel output

There are two ways to write from many MPI ranks, and the archive represents
both:

- **Global arrays** (`uniform` steps): ranks write their slices of
  one global array per field, collectively. The archive has one grid per
  step, describing the assembled array. A rank holding only its own slice
  calls `data_item_for(local, ref, dimensions=global_shape, check=False)`:
  the array supplies the dtype, and the dimension check is skipped because
  the global shape is larger than the slice on purpose.
- **Per-subdomain** (`spatial` steps): each rank writes its own piece with
  its own local numbering, and the step is a `SpatialCollection` holding
  one grid per rank. Nodes on subdomain boundaries appear in more than one
  grid. `heat_split.ymf` in the examples is laid out this way.

For the per-subdomain case, `dump_grid()` / `load_grid()` turn one grid into
a YAML string and back. That string is what travels: Proteus gathers the
per-rank grid metadata as these strings and stores them in fixed-width HDF5
datasets. (Parallel HDF5 rejects variable-length strings, so the width is
agreed across ranks with an `allreduce` of the maximum length.)

## What the writer checks, and what it doesn't

YMF splits checking across a trust boundary:

- A **front-end** (`ymf.schema` + `ymf.units`, see
  [problem-specification.md](problem-specification.md)) validates input a
  human or an LLM wrote: schema, units, dimensional consistency. It runs
  once, before the simulation, and can afford to be slow.
- A **solver** trusts the validated input and writes the archive. On that
  path the checks are only structural and cost O(1) per field.

`write_ymf()` runs `validate_domain()` by default. It checks required keys,
the closed sets above, and that no grid has two fields with the same name
and centering (viewers key fields by name within point data and within
cell data, so a duplicate silently hides the first; the same name cell-
and node-centred is fine).
`data_item_for()` and `check_dimensions()` compare the declared
`Dimensions` against the actual array while it is still in hand. A
flattened declaration (`[N*k]` for an `(N, k)` array) is accepted; a
genuine mismatch is not. Errors name the offending field:

```text
data_item_for().Dimensions: declared [81] (81 values) but the array has shape [25] (25 values)
TimeCollections[0].Data[0].Topology.Type: 'Triangel' is not a known XDMF topology type. ...
```

The first of these is a real bug this check found in Proteus: a field
sized per mesh vertex (25 values) archived through a quadratic space that
declared 81. The archive had been written that way without any error.

The writer does **not** check units, physical meaning, or standard names.
That is the front-end's job.

## Viewing: converting to XDMF

```bash
ymf2xmf run.ymf                 # writes run.xmf beside it
ymf2xmf -v run.ymf -o view.xmf  # -v summarizes what the archive holds
```

```text
run.ymf
  Mesh Spatial_Domain: 6 steps, t=0..50, Triangle topology, fields: T, q, side
wrote view.xmf
```

ParaView and VisIt open the `.xmf` directly. `extra` is carried into it as
a base64-encoded JSON `<Information Name="YMF" Value="..."/>` element, a
child of `<Domain>`, which the XDMF DTD permits
(`<!ELEMENT Domain (Information*, Grid+)>`). Viewers parse it and ignore
it. `ymf.xdmf.read_xdmf()` recovers both the domain and `extra` from an
`.xmf`, and the round trip is exact. `examples/read_archive.py` confirms
this on every run, and the example `.xmf` files have been opened in
ParaView 6.0.1, both the uniform and the per-subdomain one.

## Why pyyaml and not strictyaml

`strictyaml` is the right tool for documents people write (the problem
spec), and the wrong one for metadata a machine writes. Serializing a
realistic archive document (100 time steps × 16 ranks × 10 fields):

| serializer | time |
|---|---|
| `json.dumps` | 0.017 s |
| `xml.etree` (what Proteus did before) | 0.078 s |
| `pyyaml` `CSafeDumper` (libyaml) | 0.704 s |
| `pyyaml` `safe_dump` (pure Python) | 2.258 s |
| `strictyaml.as_document` | 18.778 s |

So the archive core uses pyyaml and prefers the libyaml C extension.
`require_libyaml()` lets a caller on a hot path fail loudly rather than
silently run 3× slower. YAML is also about 1.5× the bytes of the equivalent
XML. Proteus keeps the per-step cost O(1) by serializing each step's grids
as they are written and the full document once, at close.

## Use in Proteus

Proteus's `ymf_io` branch (aimed at the 2.0 line) replaces its hand-built
XDMF XML with this model:

- Every field goes through one `write_field` call, which builds its
  `Attribute` + `DataItem`. Fields beyond the solution itself are declared
  by the physics class, which yields `ArchiveField`s from an
  `archiveFields` hook, rather than written by hand-built XML.
- Per-step, per-rank grid metadata is stored in the HDF5 file as YAML
  (`dump_grid`). At close, `assemble_domain()` builds one domain from it and
  writes `run.ymf`.
- The HDF5 root carries a metadata format version. Archives from Proteus
  ≤ 1.9, which stored XDMF XML fragments instead, are still read, so hot
  start and post-processing keep working on old output.
- Proteus currently still writes an `.xmf` beside the `.ymf`. Dropping it
  in favour of `ymf2xmf` is the remaining step.

## Limits

- Fields carry no units. The archive says what shape and type a field is,
  not what physical quantity it is. Today that information lives in the
  problem spec in `extra`.
- One level of spatial nesting: a step is a single grid or a flat
  collection of subdomain grids.
- `extra` is stored as-is and not validated.
- Format version 1 is the only one. A reader refuses any other version
  rather than guess.
