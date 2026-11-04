# Domains beyond a box: `file`, `discrete`, and periodic on a non-box

Status: proposal 2026-10-26 (branch `domains-beyond-box`, bed1608); decided by
cekees 2026-11-04 (§ 5); built in cekees/ymf#8 and, for Proteus,
cekees/proteus#188, except Proteus's reading of a `discrete` domain, which
waits on cekees/proteus#186 (§ 3). Sections 2-4 are the proposal as written,
each followed by what was decided and built. Scope is what was left of C27
Phase 5: Neumann, Robin, tabulated hydrographs and `piecewise_linear` regions
already ran.

## 1. Where things stand

`adr_problem` (`symbolic/problem.py:_posed`) solves a `box`, or a
`piecewise_linear` region, and raises `SymbolicError` for any other kind,
naming it. It refuses periodic conditions when the domain is a region
(`problem.py`, "periodic faces are a box's"). Three facts shape the options.

- **Proteus can already build these domains.** `DomainYMF.from_ymf` turns a
  `kind: file` block into the matching `Domain` class (`gmsh_geo` and the
  mesh-file formats in `_FILE_CLASSES`). A `discrete` block becomes a
  `TriangulatedSurfaceDomain` only, and only for `Triangle` topology: a
  surface, not a volume to solve on. The refusals are in ymf's layer.
- **Everything downstream of the domain assumes a box's bounds.**
  `ADRProblem.physics` sets `L` and `x0` from `domain.lower/upper` (for a
  region, the bounding box of its vertices), and that box sets the boundary
  tolerance. A `file` domain has no vertices to take a box from until it is
  meshed, and `_posed` runs before meshing.
- **Boundary regions are coordinate predicates only.**
  `boundary_regions: [{name, where}]`, evaluated on the point. `getDBC(x, flag)`
  is already handed the mesh's boundary flag, but nothing in the spec can
  select by it. For a region that is fine (the predicate holds on any mesh of
  it). For a CAD or Gmsh model, whose physical groups are the natural names,
  writing coordinates is the wrong tool.

## 2. `file` domains

**Proposal: solve them, with boundary regions selectable by tag.**

1. `boundary_regions` entries take `tags: [ints or names]` as an alternative to
   `where` (exactly one of the two). Names resolve through the domain's
   `boundary_tags` (Gmsh physical groups come through
   `getGmshPhysicalGroups`, which `Domain.gmsh2proteus` already uses). Tag
   selection is by the `flag` argument Proteus passes, so it needs no new
   plumbing in `ADRProblem._dirichlet/_neumann/_robin`: a test of
   `flag in tags`, beside the predicate.
   This is a change to the schema and to composition (a region now has two
   forms), so it belongs in the same PR as the first solve.
2. Bounds. A `file` block gains optional `bounds: {lower, upper}`. With it,
   the box machinery is unchanged. Without it, `adr_problem` refuses with the
   fix ("give bounds, or run `pro msh` first") rather than guessing. Reading
   the bounds back from the generated mesh is nicer but moves a step that is
   pure today (`_posed` builds the problem from the document alone) to after
   meshing; do not do that until someone needs it.
3. The mesh. `discretizations[].mesh` is `cells` (box) or `h` (region). A
   file's mesh is set by the file, so `mesh: {levels}` only, and `h`/`cells`
   on a `file` domain is an error that says so. Multilevel from a file is
   whatever the Gmsh route supports; check before promising `levels > 1`.

Out of scope: formats Proteus has no `Domain` class for (`from_ymf` already
reports these by name), and STL surfaces, which are boundaries, not solvable
volumes.

**Decided and built.** D1: `tags` beside `where`, exactly one per region,
names through `boundary_tags`; a tag region reaches the solver as `flags` and
no `where` (an older solver fails on it, `KeyError: 'where'`, rather than
taking it for the whole boundary). D2: `bounds` is optional; without it the
solver takes the bounds from the mesh, so posing follows meshing for this
kind alone (Proteus reads them from a tetgen `.node` or 2dm/3dm `ND` lines;
other formats give `bounds`). The mesh: a mesh file (`2dm`, `3dm`, `hex`,
`tetgen`, `ymf`) takes `mesh: {levels}` only, keyed `levels=L`; a geometry
(`gmsh_geo`) takes `h`, which Proteus's gmsh route already uses (`-clmax
0.5 he`), so item 3 above was half wrong. `levels > 1` works from a file:
every mesh-file class refines its coarse mesh uniformly
(`generateFromExistingCoarseMesh`), and a tetgen cube read as a file
converges at rate 1.94, 1.96 over levels 1-3.

## 3. `discrete` domains

**Superseded by D4 (cekees, 2026-11-04): a `discrete` domain is solved.** It
is a mesh as ymf stores one (C7.1, mesh decision 3: a mesh is a type of
domain): `topology`, `geometry` and, from cekees/ymf#6, tagged `facets` and
`provenance`. ymf poses it whole, with its bounds from its nodes when the
document holds them inline (from the solver, or `bounds`, when they are in an
HDF5 file), its boundary regions by its facet tags, and `mesh: {levels}`
only. Proteus reads it through `MeshYmfDomain` (cekees/proteus#186, the C7.1
session's): that PR adds `ymf` as a file format; what it does not yet do is
turn an inline `discrete` block into a `MeshYmfDomain`, and the dispatch is
in `DomainYMF.py`, which #186 owns. Until it lands, Proteus refuses a
`discrete` block it cannot read, as before.

The proposal, kept for the record: **keep them carried and refused, and make the refusal useful.**
The only consumer is `TriangulatedSurfaceDomain`, a surface. A volume mesh in
the document would have to be written out as `.node/.ele` (or `.msh`) and read
through `MeshTetgenDomain`/`Mesh2DMDomain`; that is a second path to a mesh
that `file` already gives for the price of one export step. Revisit if an
archive's own grids are wanted as a domain (the "replay on the mesh I got"
case), which is a closure question, not a Phase 5 one. The change now is the
message: say "write the mesh to a file and use `kind: file`".

## 4. Periodic on a non-box domain

Today's mechanism (`ADRProblem._periodic`) is a key function: a node on either
face of a periodic axis maps to its coordinates with that axis snapped to the
lower face, and Proteus gives nodes with equal keys one dof. That is purely
geometric, so it does not care whether the domain is a box. What it needs is a
mesh where **every node on one face has a partner on the other**.

- `box` meshes (structured) have that for free.
- A `piecewise_linear` region is meshed by triangle/tetgen. If the boundary
  segments on the two faces are translations of each other and the generator
  adds no boundary points (triangle's `-Y`; tetgen's equivalent), the nodes
  match. If it does add them, the faces get different nodes and periodicity
  fails *silently* (unmatched nodes are simply not merged).

**Proposal:** allow `periodic` on a `piecewise_linear` region when all hold,
and check them at validation, not by trusting the mesh:

1. the periodic faces are bounding-box faces (the `_periodic_axes` test
   already finds them, on the region's bounding box);
2. the boundary segments (2-D) or facets (3-D) on the lower face translate
   onto those on the upper face, within tolerance;
3. the mesher is told not to refine the boundary.

After meshing, assert that every keyed node has a partner, and fail with the
count if not; this check is cheap and is the one that catches (3) being wrong.
Non-axis-aligned or rotational periodicity is out of scope: it needs a shift
vector per pair rather than a snapped coordinate, and nobody has asked.
Periodic on a `file` domain is out of scope for the same reason plus Gmsh's
`Periodic` constraint, which Proteus's reader does not carry.

**Decided and built.** ymf checks condition 2 at posing
(`periodic_faces_match`: the pieces on the lower face, moved by the box's
length, as point sets to 1e-8 of the box, with a count of those without a
partner). Proteus does 3 and the check after meshing: `Y` for triangle and
tetgen, and, since `Y` alone leaves a face with only its own vertices (2
nodes on the notched channel's x face), every 2-D boundary segment cut into
equal pieces no longer than `h` before meshing, so translated faces get
translated nodes; `run` stops with the count of unpaired face nodes. On the
notched channel the old options left 0, 3 and 6 unpaired nodes at
h = 0.2, 0.1, 0.05; now 0, with P1 and P2 converging at rates 2 and 3 on a
solution whose x-derivative is nonzero on the periodic faces. 3-D faces are
not cut (tetgen's `Y` keeps the facets' own vertices); the check reports if
that is not enough. Periodic on `file` and `discrete` domains is refused.

## 5. Order, and what to decide

1. Non-box periodic (section 4): smallest, no schema change, one test with a
   notched or L-shaped region that keeps two opposite faces.
2. `file` domains with tag selection (section 2): schema change plus
   composition; the largest.
3. `discrete`: the message only.

Decisions for cekees (answered 2026-11-04; the order was kept, periodic first,
in one round with `file` domains, D3):

- **Tag selection in `boundary_regions`** (`tags:` beside `where`): yes/no. It is
  the one schema change, and the whole value of `file` domains rests on it.
  **D1: yes, exactly one of the two.**
- **`bounds` required on `file` domains** (refuse without), or read from the
  mesh later. **D2: optional; without it, read from the mesh; posing stays
  pure whenever bounds are given.**
- **Order**: periodic first, or `file` first. **D3: one round, periodic
  first.**
- **D4 (new):** `discrete` domains solved through `MeshYmfDomain` and the
  facet tags, bounds from the nodes; § 3 rewritten.

## 6. Not verified (updated 2026-11-04)

- ~~Whether Proteus passes `-Y`~~: it did not (`VApq30Dena<area>`,
  `VApq1.35q12feena<vol>`); it does now for a periodic region.
- ~~That `levels > 1` works from a `file` domain~~: it does (tetgen cube,
  levels 1-3, rates 1.94 and 1.96).
- ~~Which classes `_FILE_CLASSES` covers~~: 2dm, 3dm, hex and tetgen, plus
  `ymf` with cekees/proteus#186; Proteus reads bounds from tetgen, 2dm and 3dm.
- The cost of the partner check on a large 3-D mesh (it is two set
  constructions over the face nodes; not timed).
- The `gmsh_geo` route with `h` and without `bounds` (Gmsh must make the mesh
  before the bounds can be read; not exercised).
- Proteus solving a `discrete` domain: waits on cekees/proteus#186.
