# Domains beyond a box: `file`, `discrete`, and periodic on a non-box

Status: proposal, 2026-10-26. Nothing here is built. Scope is what is left of
C27 Phase 5 after the check recorded in the roster: Neumann, Robin, tabulated
hydrographs and `piecewise_linear` regions already run.

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

## 3. `discrete` domains

**Proposal: keep them carried and refused, and make the refusal useful.**
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

## 5. Order, and what to decide

1. Non-box periodic (section 4): smallest, no schema change, one test with a
   notched or L-shaped region that keeps two opposite faces.
2. `file` domains with tag selection (section 2): schema change plus
   composition; the largest.
3. `discrete`: the message only.

Decisions for cekees:

- **Tag selection in `boundary_regions`** (`tags:` beside `where`): yes/no. It is
  the one schema change, and the whole value of `file` domains rests on it.
- **`bounds` required on `file` domains** (refuse without), or read from the
  mesh later.
- **Order**: periodic first, or `file` first.

## 6. Not verified

- Whether the triangle/tetgen invocation in `proteus` passes `-Y`, or whether
  `MeshOptions` exposes it.
- That `levels > 1` works from a `file` domain.
- Which `Domain` classes `_FILE_CLASSES` covers beyond `gmsh_geo` (I read the
  dispatch, not each class).
- The cost of the post-mesh partner check on a large 3-D mesh.
