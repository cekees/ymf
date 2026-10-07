"""The YMF archive core: the domain data model, plus YAML read/write.

This is the module simulation codes depend on. Its **only** third-party
dependency is ``pyyaml`` -- deliberately, so that a solver like Proteus can
write YMF archives without inheriting ``strictyaml``, ``pint`` or ``sympy``
(see ``ymf.schema`` / ``ymf.units``, which are behind the ``ymf[spec]`` and
``ymf[units]`` extras).

The trust boundary
------------------
Validation and writing are different jobs with different dependency
budgets:

- A **front-end** (``ymf.schema`` + ``ymf.units``, plus whatever tooling
  wraps them) validates YMF *input* authored by a human or an LLM. It does
  schema validation, unit coherence, dimensional analysis, standard-name
  checking. It runs once, off the compute path, and may be as slow and as
  thorough as it likes.
- A **solver** reads validated input and writes archive *output*. It
  trusts the input's semantics and does not re-derive them.

So nothing here checks units, dimensional consistency or standard names.
What it *does* check is **structure**: that required keys are present, that
closed-set fields hold legal values, and that declared ``Dimensions`` match
the arrays they claim to describe. Those checks are O(1) per field and
fail loudly, naming the field. Trust, but verify cheaply.

Why not strictyaml here
-----------------------
``strictyaml`` is the right tool for human-authored documents -- schema
validation, comment preservation, exact round-tripping -- and the wrong
tool for machine-written archive metadata. Measured on a realistic archive
document (100 timesteps x 16 ranks x 10 attributes):

======================================  =========
serializer                              time
======================================  =========
``json.dumps``                            0.017 s
``xml.etree`` (what Proteus did before)    0.078 s
``pyyaml`` ``CSafeDumper`` (libyaml)       0.704 s
``pyyaml`` ``safe_dump`` (pure Python)     2.258 s
``strictyaml.as_document``                18.778 s
======================================  =========

That is ~240x the ElementTree code it replaces. This module therefore uses
``pyyaml``, and prefers the libyaml-backed C dumper/loader -- the 3.2x gap
between ``CSafeDumper`` and the pure-Python fallback is large enough to be
worth knowing about, so :func:`require_libyaml` is provided for callers on
a hot path to assert on it rather than silently degrade.

The domain model
----------------
A *domain* is a plain dict -- no classes, no schema objects, so it costs
nothing to build, gathers across MPI natively, and serializes directly::

    {"TimeCollections": [{"Name": str, "Data": [step, ...]}, ...]}

There is a *list* of time collections, not one, because a solver writes one
per finite-element space: a real Proteus archive commonly holds both
``Mesh Spatial_Domain`` (the linear base mesh) and ``Mesh_c0p2_Lagrange``
(the quadratic space, on the same elements with more nodes). 195 of the 814
archives surveyed when this model was designed had two. A single-collection
document may be written as ``{"TimeCollection": {...}}`` for convenience;
:func:`canonicalize_domain` normalizes it to the list form, which is the
only form written to disk.

Each ``step`` is one instant in time and takes exactly one of two forms:

*uniform* -- a single grid, which is what a serial run or a globally
synchronized parallel write produces::

    {"Time": float, "Topology": T, "Geometry": G, "Attributes": [A, ...]}

*spatial* -- a collection of grids, one per MPI rank, which is what a
subdomain-per-rank write produces::

    {"Time": float, "SpatialCollection": [grid, ...]}

and the pieces are::

    grid = {"Name": str?, "Topology": T, "Geometry": G, "Attributes": [A, ...]}
    T    = {"Type": str, "NumberOfElements": int,
            "NodesPerElement": int?, "DataItem": D}
    G    = {"Type": str, "DataItem": D}
    A    = {"Name": str, "AttributeType": str, "Center": str, "DataItem": D}
    D    = {"Format": str, "DataType": str, "Precision": int,
            "Dimensions": [int, ...], "Data": str}      # Format="HDF"
         | {"Format": "XML", ..., "Include": str}       # text sidecar file
         | {"Format": "XML", ..., "Values": [num, ...]} # inline, row-major

Inline ``Values`` make an archive self-contained -- no HDF5 file, no
sidecar -- which suits small problems and solvers that avoid an HDF5
dependency. ``ymf2xmf`` writes them as the DataItem's text, as XDMF does.

Build these with the constructor helpers below (:func:`data_item`,
:func:`attribute`, :func:`topology`, :func:`geometry`, :func:`grid`,
:func:`new_domain`, :func:`add_uniform_step`, :func:`add_spatial_step`)
rather than by hand -- they apply defaults consistently, which is what
makes :func:`canonicalize_domain` able to state a complete list of
defaultable fields instead of guessing.

``ymf.xdmf`` converts a domain built here into an ``.xmf`` file that
ParaView and VisIt read natively. Nothing in this module knows about XML.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

try:  # pragma: no cover - depends on how pyyaml was built
    from yaml import CSafeDumper as _YamlDumper, CSafeLoader as _YamlLoader

    HAVE_LIBYAML = True
except ImportError:  # pragma: no cover
    from yaml import SafeDumper as _YamlDumper, SafeLoader as _YamlLoader

    HAVE_LIBYAML = False


__all__ = [
    "YmfArchiveError",
    "HAVE_LIBYAML",
    "require_libyaml",
    "KNOWN_TOPOLOGY_TYPES",
    "ATTRIBUTE_TYPES",
    "CENTERINGS",
    "DATA_TYPES",
    "DATA_ITEM_FORMATS",
    "GEOMETRY_TYPES",
    "data_item",
    "data_item_for",
    "attribute",
    "topology",
    "geometry",
    "grid",
    "new_domain",
    "add_collection",
    "add_uniform_step",
    "add_spatial_step",
    "canonicalize_domain",
    "validate_domain",
    "check_dimensions",
    "write_ymf",
    "load_array",
    "inline_domain",
    "domain_arrays",
    "read_ymf",
    "dump_grid",
    "load_grid",
    "ARCHIVE_FORMAT_VERSION",
]


class YmfArchiveError(ValueError):
    """A YMF archive document is structurally invalid.

    Raised by :func:`validate_domain` and :func:`check_dimensions`. The
    message always names the offending path within the document (e.g.
    ``TimeCollections[0].Data[3].Attributes[1].Center``) so a failure points at
    the field that caused it rather than at the serializer.
    """


#: Bumped when the on-disk archive layout changes incompatibly. Written
#: into every document by :func:`write_ymf` and checked by :func:`read_ymf`,
#: so a version mismatch produces a clear error instead of a confusing
#: parse failure somewhere downstream.
ARCHIVE_FORMAT_VERSION = 1

#: XDMF 2.0 topology types. Closed on purpose: getting one of these wrong
#: yields an ``.xmf`` that a viewer rejects (or worse, silently misdraws),
#: and the set a solver can legitimately emit is knowable. Structured-mesh
#: types are included because XDMF defines them, even though the current
#: Proteus writers only emit unstructured ones.
KNOWN_TOPOLOGY_TYPES = frozenset(
    {
        # unstructured, linear
        "Polyvertex",
        "Polyline",
        "Polygon",
        "Triangle",
        "Quadrilateral",
        "Tetrahedron",
        "Pyramid",
        "Wedge",
        "Hexahedron",
        # unstructured, quadratic
        "Edge_3",
        "Tri_6",
        "Triangle_6",
        "Quad_8",
        "Quadrilateral_8",
        "Tet_10",
        "Tetrahedron_10",
        "Pyramid_13",
        "Wedge_15",
        "Hex_20",
        "Hexahedron_20",
        # heterogeneous
        "Mixed",
        # structured
        "2DSMesh",
        "2DRectMesh",
        "2DCoRectMesh",
        "3DSMesh",
        "3DRectMesh",
        "3DCoRectMesh",
    }
)

#: Legal ``AttributeType`` values -- the rank of a field.
ATTRIBUTE_TYPES = frozenset({"Scalar", "Vector", "Tensor", "Tensor6", "Matrix", "GlobalID"})

#: Legal ``Center`` values -- what mesh entity a field is attached to.
CENTERINGS = frozenset({"Node", "Cell", "Face", "Edge", "Grid"})

#: Legal ``DataType`` values for a DataItem.
DATA_TYPES = frozenset({"Float", "Int", "UInt", "Char", "UChar"})

#: Legal ``Format`` values. ``HDF`` references an HDF5 dataset; ``XML``
#: means the values live in a text sidecar file pulled in by reference;
#: ``Binary`` references a raw binary file.
DATA_ITEM_FORMATS = frozenset({"HDF", "XML", "Binary"})

#: Legal ``Geometry`` types.
GEOMETRY_TYPES = frozenset({"XYZ", "XY", "X_Y_Z", "X_Y", "VXVYVZ", "ORIGIN_DXDYDZ", "ORIGIN_DXDY"})

# XDMF's own documented default for DataItem/@Precision. Faithful to the
# format rather than to what a solver usually writes -- use data_item_for()
# to infer the right value from an actual array instead of relying on this.
_DEFAULT_PRECISION = 4
_DEFAULT_TIME_COLLECTION_NAME = "TimeCollection"


def require_libyaml() -> None:
    """Raise unless ``pyyaml`` was built with its libyaml C extension.

    The pure-Python fallback is ~3.2x slower to dump and load. That is
    tolerable for a one-off read but not for a solver writing an archive
    every timestep, so a caller on that path should assert on it at
    startup rather than discover it as unexplained slowness later.
    """
    if not HAVE_LIBYAML:
        raise YmfArchiveError(
            "pyyaml is installed without its libyaml C extension, which makes "
            "YMF archive read/write about 3.2x slower. Install a pyyaml built "
            "with libyaml (the conda-forge and manylinux wheels both are), or "
            "call write_ymf/read_ymf without require_libyaml() to accept the "
            "slower pure-Python path."
        )


# ---------------------------------------------------------------------------
# constructors
#
# These exist so that defaults are applied in exactly one place. That is
# what lets canonicalize_domain() enumerate every defaultable field rather
# than covering only the ones that happened to surface in a test.
# ---------------------------------------------------------------------------


def data_item(
    dimensions: Sequence[int],
    data: Optional[str] = None,
    *,
    include: Optional[str] = None,
    values: Optional[Sequence[Any]] = None,
    data_type: str = "Float",
    precision: int = _DEFAULT_PRECISION,
    fmt: Optional[str] = None,
) -> Dict[str, Any]:
    """Describe where an array lives.

    Parameters
    ----------
    dimensions:
        The array's shape, as a sequence of ints.
    data:
        An HDF5 reference of the form ``"file.h5:/dataset_name"``. Give
        this for ``Format="HDF"`` (the default).
    include:
        Path to a text sidecar file holding the values. Give this instead
        of ``data`` for ``Format="XML"``, the no-HDF5 fallback.
    values:
        The values themselves, flattened in row-major order, for an inline
        ``Format="XML"`` DataItem: the archive then needs no other file.
        Exactly one of ``data`` / ``include`` / ``values`` must be given.
    data_type:
        One of :data:`DATA_TYPES`.
    precision:
        Bytes per value. Defaults to XDMF's own default of 4, which is
        very often *not* what you want -- ``float64`` needs 8. Prefer
        :func:`data_item_for`, which reads it off the array.
    fmt:
        One of :data:`DATA_ITEM_FORMATS`. Inferred from which of
        ``data``/``include``/``values`` was given, so it rarely needs passing.
    """
    if sum(x is not None for x in (data, include, values)) != 1:
        raise YmfArchiveError(
            "data_item() needs exactly one of data= (an HDF5 reference like "
            "'out.h5:/u_t3'), include= (a path to a text sidecar file) or "
            "values= (the values inline); "
            f"got data={data!r}, include={include!r}, "
            f"values={'<%d values>' % len(values) if values is not None else None}"
        )
    if fmt is None:
        fmt = "HDF" if data is not None else "XML"
    item: Dict[str, Any] = {
        "Format": fmt,
        "DataType": data_type,
        "Precision": int(precision),
        "Dimensions": _FlowList(int(d) for d in dimensions),
    }
    if data is not None:
        item["Data"] = str(data)
    elif include is not None:
        item["Include"] = str(include)
    else:
        item["Values"] = _inline_values(values, data_type, "data_item()")
        check_dimensions(item, [len(item["Values"])], where="data_item()")
    return item


class _FlowList(list):
    """A list written as one ``[a, b, c]`` line rather than one per item.

    ``Dimensions`` and inline ``Values`` use it, so a shape reads
    ``[8, 3]`` and an array costs one (wrapped) line rather than a line
    per number. It is a plain list in every other respect.
    """


#: DataTypes whose values are integers.
_INTEGER_DATA_TYPES = frozenset({"Int", "UInt", "Char", "UChar"})


def _flatten(values: Any) -> List[Any]:
    """Flatten nested sequences (or anything with ``.tolist()``) row-major."""
    if hasattr(values, "tolist"):  # numpy arrays, without importing numpy
        values = values.tolist()
    if not isinstance(values, (list, tuple)):
        return [values]
    flat: List[Any] = []
    for v in values:
        flat.extend(_flatten(v) if isinstance(v, (list, tuple)) else [v])
    return flat


def _inline_values(values: Any, data_type: str, where: str) -> "_FlowList":
    """Inline values as a flat list of numbers of the declared type."""
    integral = data_type in _INTEGER_DATA_TYPES
    result = _FlowList()
    for i, v in enumerate(_flatten(values)):
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            raise YmfArchiveError(
                f"{where}.Values[{i}]: {v!r} is not a number")
        if integral:
            if v != int(v):
                raise YmfArchiveError(
                    f"{where}.Values[{i}]: {v!r} is not an integer, but the "
                    f"DataType is {data_type}")
            v = int(v)
        else:
            v = float(v)
        result.append(v)
    return result


#: numpy dtype kind -> XDMF DataType. Covers the kinds that can appear in a
#: mesh or field array; anything else is rejected rather than guessed at.
_DTYPE_KIND_TO_DATA_TYPE = {"f": "Float", "i": "Int", "u": "UInt", "b": "UChar"}


def data_item_for(
    array: Any,
    data: Optional[str] = None,
    *,
    include: Optional[str] = None,
    inline: bool = False,
    dimensions: Optional[Sequence[int]] = None,
    check: bool = True,
) -> Dict[str, Any]:
    """Build a DataItem by reading shape and dtype off ``array``.

    With ``inline=True`` the values are copied into the DataItem itself
    (``Format="XML"``, ``Values``) instead of being referenced, so the
    archive needs no HDF5 file or sidecar.

    Prefer this over :func:`data_item` whenever the array is in hand. It
    removes the three fields a caller can most easily get wrong --
    ``Dimensions``, ``DataType`` and ``Precision`` -- by deriving all three
    from the array itself, so a ``float64`` field cannot be described as
    4-byte, and a connectivity array cannot claim the wrong element count.

    ``array`` needs a numpy-style ``.shape`` and ``.dtype``; numpy is not
    imported here, so anything with those attributes works.

    Pass ``dimensions`` only to declare a deliberately different shape from
    the array's own -- flattening an ``(N, k)`` connectivity array to
    ``[N*k]`` is the normal case. It is still checked for consistency, so a
    genuine mismatch is rejected.

    Pass ``check=False`` when ``array`` is deliberately *not* the whole of
    what the DataItem describes, and is present only to supply the dtype.
    That is the case for a collective parallel write: the DataItem covers
    the assembled global array while this caller holds one rank's slice, so
    the declared dimensions are larger than the array by design and
    checking them would reject a correct write.
    """
    try:
        shape, dtype = array.shape, array.dtype
    except AttributeError as exc:  # pragma: no cover - defensive
        raise YmfArchiveError(
            "data_item_for() needs an array with .shape and .dtype; got "
            f"{type(array).__name__}. Use data_item() to state Dimensions, "
            "DataType and Precision explicitly."
        ) from exc

    kind = dtype.kind
    if kind not in _DTYPE_KIND_TO_DATA_TYPE:
        raise YmfArchiveError(
            f"data_item_for(): dtype {dtype!r} (kind {kind!r}) has no XDMF "
            f"DataType equivalent. Supported kinds: "
            f"{sorted(_DTYPE_KIND_TO_DATA_TYPE)}"
        )

    if inline and (data is not None or include is not None):
        raise YmfArchiveError(
            "data_item_for(): inline=True copies the values in, so it takes "
            "no data= or include= reference")
    item = data_item(
        list(shape) if dimensions is None else dimensions,
        data,
        include=include,
        values=array if inline else None,
        data_type=_DTYPE_KIND_TO_DATA_TYPE[kind],
        precision=dtype.itemsize,
    )
    if dimensions is not None and check:
        check_dimensions(item, shape, where="data_item_for()")
    return item


def attribute(
    name: str,
    item: Dict[str, Any],
    *,
    center: str = "Node",
    attribute_type: str = "Scalar",
) -> Dict[str, Any]:
    """Describe one field defined over a grid.

    ``center`` and ``attribute_type`` are the XDMF ``Center`` and
    ``AttributeType``; both are validated against closed sets, so a
    malformed field is not representable.
    """
    return {
        "Name": str(name),
        "AttributeType": attribute_type,
        "Center": center,
        "DataItem": item,
    }


def topology(
    topology_type: str,
    n_elements: int,
    item: Dict[str, Any],
    *,
    nodes_per_element: Optional[int] = None,
) -> Dict[str, Any]:
    """Describe a grid's element connectivity.

    ``nodes_per_element`` is only meaningful for the types that don't imply
    it -- ``Polyvertex`` and ``Polyline`` -- and is omitted otherwise.
    """
    topo: Dict[str, Any] = {
        "Type": topology_type,
        "NumberOfElements": int(n_elements),
        "DataItem": item,
    }
    if nodes_per_element is not None:
        topo["NodesPerElement"] = int(nodes_per_element)
    return topo


def geometry(item: Dict[str, Any], *, geometry_type: str = "XYZ") -> Dict[str, Any]:
    """Describe a grid's node coordinates."""
    return {"Type": geometry_type, "DataItem": item}


def grid(
    topo: Dict[str, Any],
    geom: Dict[str, Any],
    attributes: Optional[Iterable[Dict[str, Any]]] = None,
    *,
    name: Optional[str] = None,
) -> Dict[str, Any]:
    """Assemble one uniform grid: connectivity, coordinates, and fields."""
    g: Dict[str, Any] = {"Topology": topo, "Geometry": geom}
    if name is not None:
        g["Name"] = str(name)
    g["Attributes"] = list(attributes) if attributes is not None else []
    return g


def new_domain(name: str = _DEFAULT_TIME_COLLECTION_NAME) -> Dict[str, Any]:
    """Create a domain with one empty, named time collection.

    Add further collections with :func:`add_collection` -- one per
    finite-element space that needs its own mesh.
    """
    return {"TimeCollections": [{"Name": str(name), "Data": []}]}


def add_collection(domain: Dict[str, Any], name: str) -> Dict[str, Any]:
    """Append a new named time collection. Returns it.

    Use one collection per finite-element space: the nodes of a quadratic
    space are not the nodes of the linear mesh, so they cannot share a
    grid, and a viewer needs to see them as separate meshes over the same
    elements.
    """
    collections = domain.setdefault("TimeCollections", [])
    if any(c.get("Name") == name for c in collections):
        raise YmfArchiveError(
            f"domain already has a time collection named {name!r}; collection "
            "names are how a reader tells the meshes apart, so they must be unique"
        )
    collection: Dict[str, Any] = {"Name": str(name), "Data": []}
    collections.append(collection)
    return collection


def _resolve_collection(domain: Dict[str, Any], collection: Any) -> Dict[str, Any]:
    """Find a collection by index or by name, with a legible failure."""
    collections = domain.get("TimeCollections")
    if not collections:
        raise YmfArchiveError(
            "domain has no time collections; build it with new_domain()"
        )
    if isinstance(collection, int):
        try:
            return collections[collection]
        except IndexError:
            raise YmfArchiveError(
                f"domain has {len(collections)} time collection(s); no index "
                f"{collection}"
            ) from None
    for c in collections:
        if c.get("Name") == collection:
            return c
    names = [c.get("Name") for c in collections]
    raise YmfArchiveError(
        f"domain has no time collection named {collection!r}; it has {names}"
    )


def add_uniform_step(
    domain: Dict[str, Any],
    time: float,
    topo: Dict[str, Any],
    geom: Dict[str, Any],
    attributes: Optional[Iterable[Dict[str, Any]]] = None,
    *,
    collection: Any = 0,
) -> Dict[str, Any]:
    """Append a timestep holding a single grid. Returns the new step.

    ``collection`` selects which time collection to append to, by index or
    by name; it defaults to the first, which is the only one in a
    single-mesh archive.
    """
    step: Dict[str, Any] = {"Time": float(time), "Topology": topo, "Geometry": geom}
    step["Attributes"] = list(attributes) if attributes is not None else []
    _resolve_collection(domain, collection)["Data"].append(step)
    return step


def add_spatial_step(
    domain: Dict[str, Any],
    time: float,
    grids: Iterable[Dict[str, Any]],
    *,
    collection: Any = 0,
) -> Dict[str, Any]:
    """Append a timestep holding one grid per subdomain. Returns the new step."""
    step = {"Time": float(time), "SpatialCollection": list(grids)}
    _resolve_collection(domain, collection)["Data"].append(step)
    return step


# ---------------------------------------------------------------------------
# canonicalization
# ---------------------------------------------------------------------------


def _canonicalize_data_item(item: Dict[str, Any], where: str) -> Dict[str, Any]:
    """Fill in every defaultable DataItem field.

    Serialization loses the distinction between "field omitted, assume the
    default" and "field explicitly set to the default" -- both come back
    explicit. So comparing a sparse hand-written dict against one that has
    been through a write/read cycle only works if the sparse one is
    canonicalized first. This is the single place that decides what those
    defaults are.
    """
    if "Dimensions" not in item:
        raise YmfArchiveError(f"{where}: DataItem has no 'Dimensions'")
    present = [k for k in ("Data", "Include", "Values") if k in item]
    if len(present) != 1:
        raise YmfArchiveError(
            f"{where}: DataItem needs exactly one of 'Data' (HDF5 reference), "
            f"'Include' (text sidecar path) or 'Values' (inline values); got "
            f"{' and '.join(present) if present else 'none'}"
        )
    has_data = present == ["Data"]
    fmt = item.get("Format", "HDF" if has_data else "XML")
    canonical: Dict[str, Any] = {
        "Format": fmt,
        "DataType": item.get("DataType", "Float"),
        "Precision": int(item.get("Precision", _DEFAULT_PRECISION)),
        "Dimensions": _FlowList(int(d) for d in item["Dimensions"]),
    }
    if has_data:
        canonical["Data"] = str(item["Data"])
    elif "Include" in item:
        canonical["Include"] = str(item["Include"])
    else:
        canonical["Values"] = _inline_values(
            item["Values"], canonical["DataType"], where)
    return canonical


def _canonicalize_attribute(attr: Dict[str, Any], where: str) -> Dict[str, Any]:
    for key in ("Name", "DataItem"):
        if key not in attr:
            raise YmfArchiveError(f"{where}: Attribute has no {key!r}")
    return {
        "Name": str(attr["Name"]),
        "AttributeType": attr.get("AttributeType", "Scalar"),
        "Center": attr.get("Center", "Node"),
        "DataItem": _canonicalize_data_item(attr["DataItem"], f"{where}.DataItem"),
    }


def _canonicalize_topology(topo: Dict[str, Any], where: str) -> Dict[str, Any]:
    for key in ("Type", "NumberOfElements", "DataItem"):
        if key not in topo:
            raise YmfArchiveError(f"{where}: Topology has no {key!r}")
    canonical: Dict[str, Any] = {
        "Type": topo["Type"],
        "NumberOfElements": int(topo["NumberOfElements"]),
        "DataItem": _canonicalize_data_item(topo["DataItem"], f"{where}.DataItem"),
    }
    # NodesPerElement is genuinely optional rather than defaultable: there is
    # no correct value to invent for a Tetrahedron, so an absent key stays
    # absent and round-trips as absent.
    if "NodesPerElement" in topo:
        canonical["NodesPerElement"] = int(topo["NodesPerElement"])
    return canonical


def _canonicalize_geometry(geom: Dict[str, Any], where: str) -> Dict[str, Any]:
    if "DataItem" not in geom:
        raise YmfArchiveError(f"{where}: Geometry has no 'DataItem'")
    return {
        "Type": geom.get("Type", "XYZ"),
        "DataItem": _canonicalize_data_item(geom["DataItem"], f"{where}.DataItem"),
    }


def _canonicalize_grid(g: Dict[str, Any], where: str) -> Dict[str, Any]:
    canonical: Dict[str, Any] = {}
    if "Name" in g:
        canonical["Name"] = str(g["Name"])
    canonical["Topology"] = _canonicalize_topology(g["Topology"], f"{where}.Topology")
    canonical["Geometry"] = _canonicalize_geometry(g["Geometry"], f"{where}.Geometry")
    canonical["Attributes"] = [
        _canonicalize_attribute(a, f"{where}.Attributes[{i}]")
        for i, a in enumerate(g.get("Attributes", []))
    ]
    return canonical


def _canonicalize_step(step: Dict[str, Any], where: str) -> Dict[str, Any]:
    if "Time" not in step:
        raise YmfArchiveError(f"{where}: step has no 'Time'")
    is_spatial = "SpatialCollection" in step
    is_uniform = "Topology" in step or "Geometry" in step
    if is_spatial and is_uniform:
        raise YmfArchiveError(
            f"{where}: step has both 'SpatialCollection' and inline "
            "'Topology'/'Geometry'; a step is either one uniform grid or a "
            "collection of subdomain grids, not both"
        )
    if not is_spatial and not is_uniform:
        raise YmfArchiveError(
            f"{where}: step has neither 'SpatialCollection' nor "
            "'Topology'/'Geometry'"
        )
    canonical: Dict[str, Any] = {"Time": float(step["Time"])}
    if is_spatial:
        canonical["SpatialCollection"] = [
            _canonicalize_grid(g, f"{where}.SpatialCollection[{i}]")
            for i, g in enumerate(step["SpatialCollection"])
        ]
    else:
        for key in ("Topology", "Geometry"):
            if key not in step:
                raise YmfArchiveError(f"{where}: uniform step has no {key!r}")
        canonical["Topology"] = _canonicalize_topology(step["Topology"], f"{where}.Topology")
        canonical["Geometry"] = _canonicalize_geometry(step["Geometry"], f"{where}.Geometry")
        canonical["Attributes"] = [
            _canonicalize_attribute(a, f"{where}.Attributes[{i}]")
            for i, a in enumerate(step.get("Attributes", []))
        ]
    return canonical


def canonicalize_domain(domain: Dict[str, Any]) -> Dict[str, Any]:
    """Return ``domain`` with every defaultable field filled in explicitly.

    Use this before comparing a hand-written sparse domain against one read
    back from a file -- see :func:`_canonicalize_data_item` for why a direct
    ``==`` can otherwise fail without anything having been lost.

    The defaults applied are exactly: ``DataItem.Format`` (inferred from
    ``Data`` vs ``Include``), ``DataItem.DataType`` (``Float``),
    ``DataItem.Precision`` (4, XDMF's default),
    ``Attribute.AttributeType`` (``Scalar``),
    ``Attribute.Center`` (``Node``), ``Geometry.Type`` (``XYZ``),
    ``Attributes`` (empty list), and ``TimeCollections[i].Name``. A singular
    ``TimeCollection`` is normalized into a one-element ``TimeCollections``.
    Every other
    field is required, and ``Topology.NodesPerElement`` is optional-absent
    rather than defaulted.
    """
    collections = domain.get("TimeCollections")
    if collections is None:
        single = domain.get("TimeCollection")
        if single is None:
            return {}
        # Singular form is input sugar; the canonical output is always the list.
        collections = [single]
    result_collections = []
    seen: Dict[str, int] = {}
    for ci, collection in enumerate(collections):
        name = collection.get("Name", _DEFAULT_TIME_COLLECTION_NAME)
        if name in seen:
            raise YmfArchiveError(
                f"TimeCollections[{ci}].Name: duplicate collection name {name!r} "
                f"(also at TimeCollections[{seen[name]}]); names are how a reader "
                "tells the meshes apart, so they must be unique"
            )
        seen[name] = ci
        result_collections.append(
            {
                "Name": name,
                "Data": [
                    _canonicalize_step(step, f"TimeCollections[{ci}].Data[{i}]")
                    for i, step in enumerate(collection.get("Data", []))
                ],
            }
        )
    return {"TimeCollections": result_collections}


# ---------------------------------------------------------------------------
# structural validation -- the cheap half of the trust boundary
# ---------------------------------------------------------------------------


def _validate_data_item(item: Dict[str, Any], where: str) -> None:
    fmt = item["Format"]
    if fmt not in DATA_ITEM_FORMATS:
        raise YmfArchiveError(
            f"{where}.Format: {fmt!r} is not one of {sorted(DATA_ITEM_FORMATS)}"
        )
    if item["DataType"] not in DATA_TYPES:
        raise YmfArchiveError(
            f"{where}.DataType: {item['DataType']!r} is not one of {sorted(DATA_TYPES)}"
        )
    if item["Precision"] <= 0:
        raise YmfArchiveError(f"{where}.Precision: must be positive, got {item['Precision']}")
    dims = item["Dimensions"]
    if not dims:
        raise YmfArchiveError(f"{where}.Dimensions: must not be empty")
    if any(d < 0 for d in dims):
        raise YmfArchiveError(f"{where}.Dimensions: must be non-negative, got {dims}")
    if fmt == "XML" and "Include" not in item and "Values" not in item:
        raise YmfArchiveError(f"{where}: Format='XML' requires 'Include' or 'Values'")
    if "Values" in item:
        if fmt != "XML":
            raise YmfArchiveError(
                f"{where}: inline 'Values' need Format='XML', not {fmt!r}")
        check_dimensions(item, [len(item["Values"])], where=where)
    if fmt == "HDF" and "Data" not in item:
        raise YmfArchiveError(f"{where}: Format='HDF' requires 'Data'")


def _validate_grid_like(g: Dict[str, Any], where: str) -> None:
    topo = g["Topology"]
    if topo["Type"] not in KNOWN_TOPOLOGY_TYPES:
        raise YmfArchiveError(
            f"{where}.Topology.Type: {topo['Type']!r} is not a known XDMF "
            f"topology type. Known types: {sorted(KNOWN_TOPOLOGY_TYPES)}"
        )
    _validate_data_item(topo["DataItem"], f"{where}.Topology.DataItem")

    geom = g["Geometry"]
    if geom["Type"] not in GEOMETRY_TYPES:
        raise YmfArchiveError(
            f"{where}.Geometry.Type: {geom['Type']!r} is not one of {sorted(GEOMETRY_TYPES)}"
        )
    _validate_data_item(geom["DataItem"], f"{where}.Geometry.DataItem")

    seen: Dict[Tuple[str, str], int] = {}
    for i, attr in enumerate(g["Attributes"]):
        at = f"{where}.Attributes[{i}]"
        if attr["AttributeType"] not in ATTRIBUTE_TYPES:
            raise YmfArchiveError(
                f"{at}.AttributeType: {attr['AttributeType']!r} is not one of "
                f"{sorted(ATTRIBUTE_TYPES)}"
            )
        if attr["Center"] not in CENTERINGS:
            raise YmfArchiveError(
                f"{at}.Center: {attr['Center']!r} is not one of {sorted(CENTERINGS)}"
            )
        # Viewers keep fields on different entities apart (point data and
        # cell data are separate namespaces in VTK), so the same name may
        # appear once per centering -- exporters commonly write T both
        # cell- and node-centred. Twice on the same entity is a duplicate.
        key = (attr["Name"], attr["Center"])
        if key in seen:
            raise YmfArchiveError(
                f"{at}.Name: duplicate {attr['Center']}-centred attribute name "
                f"{attr['Name']!r} (also at {where}.Attributes[{seen[key]}]); "
                "viewers key fields by name, so duplicates silently shadow "
                "each other"
            )
        seen[key] = i
        _validate_data_item(attr["DataItem"], f"{at}.DataItem")


def validate_domain(domain: Dict[str, Any]) -> Dict[str, Any]:
    """Check ``domain`` is structurally sound; return it canonicalized.

    This is the solver-side half of the trust boundary: it verifies
    structure, not semantics. It checks that required keys are present,
    that closed-set fields (``AttributeType``, ``Center``, ``DataType``,
    ``Format``, topology and geometry types) hold legal values, and that no
    grid has two attributes of the same name. It does **not** check units,
    dimensional consistency, or standard names -- see the module docstring.

    Raises :class:`YmfArchiveError` naming the offending path.
    """
    canonical = canonicalize_domain(domain)
    for ci, collection in enumerate(canonical.get("TimeCollections", [])):
        for i, step in enumerate(collection["Data"]):
            where = f"TimeCollections[{ci}].Data[{i}]"
            if "SpatialCollection" in step:
                if not step["SpatialCollection"]:
                    raise YmfArchiveError(f"{where}.SpatialCollection: must not be empty")
                for j, g in enumerate(step["SpatialCollection"]):
                    _validate_grid_like(g, f"{where}.SpatialCollection[{j}]")
            else:
                _validate_grid_like(step, where)
    return canonical


def check_dimensions(item: Dict[str, Any], shape: Sequence[int], where: str = "DataItem") -> None:
    """Assert a DataItem's declared ``Dimensions`` match an array's shape.

    The cheapest high-value check at the boundary: a mismatch here is the
    failure mode where an archive is written without error and then loads
    as garbage, or not at all, in a viewer. Call it where the array is
    still in hand.

    ``shape`` is compared flattened-length-wise as well as element-wise, so
    declaring ``[N*k]`` for an ``(N, k)`` array is accepted -- that is a
    normal and intentional thing to do for XDMF connectivity arrays.
    """
    declared = [int(d) for d in item["Dimensions"]]
    actual = [int(s) for s in shape]
    if declared == actual:
        return
    n_declared, n_actual = 1, 1
    for d in declared:
        n_declared *= d
    for s in actual:
        n_actual *= s
    if n_declared != n_actual:
        raise YmfArchiveError(
            f"{where}.Dimensions: declared {declared} ({n_declared} values) but "
            f"the array has shape {actual} ({n_actual} values)"
        )


# ---------------------------------------------------------------------------
# YAML serialization
# ---------------------------------------------------------------------------


class _ArchiveDumper(_YamlDumper):
    """The libyaml dumper, writing inline ``Values`` on one line."""


_ArchiveDumper.add_representer(
    _FlowList,
    lambda dumper, data: dumper.represent_sequence(
        "tag:yaml.org,2002:seq", data, flow_style=True))


def _dump(obj: Any) -> str:
    if not HAVE_LIBYAML:
        warnings.warn(
            "pyyaml has no libyaml C extension; YMF archive writes are about "
            "3.2x slower than they need to be. Call "
            "ymf.archive.require_libyaml() at startup to make this an error.",
            RuntimeWarning,
            stacklevel=3,
        )
    return yaml.dump(
        obj, Dumper=_ArchiveDumper, default_flow_style=False, sort_keys=False,
        allow_unicode=True,
    )


def _load(text: str) -> Any:
    return yaml.load(text, Loader=_YamlLoader)


def write_ymf(
    domain: Dict[str, Any],
    path: str | Path,
    extra: Optional[Any] = None,
    *,
    validate: bool = True,
) -> None:
    """Write a domain as a ``.ymf`` YAML document.

    Parameters
    ----------
    domain:
        The domain dict. Canonicalized on the way out, so the file always
        has every defaultable field explicit.
    extra:
        Optional non-archive YMF content (``Problem``, ``solution_paths``,
        ``vvuq``, ...) to carry alongside. Unlike the XDMF path, which has
        to base64 it into an ``<Information>`` element, YAML holds it
        directly and readably.
    validate:
        Run :func:`validate_domain` first. Leave this on unless you are
        deliberately writing a document you know to be partial.
    """
    checked = validate_domain(domain) if validate else canonicalize_domain(domain)
    document: Dict[str, Any] = {"ymf_archive_version": ARCHIVE_FORMAT_VERSION, "domain": checked}
    if extra is not None:
        document["extra"] = extra
    Path(path).write_text(_dump(document), encoding="utf-8")


def read_ymf(path: str | Path) -> Tuple[Dict[str, Any], Optional[Any]]:
    """Read a ``.ymf`` document written by :func:`write_ymf`.

    Returns a ``(domain, extra)`` tuple. Raises :class:`YmfArchiveError` on
    a format-version mismatch, so an archive from an incompatible writer
    fails with a clear message rather than a confusing downstream
    ``KeyError``.
    """
    document = _load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(document, dict):
        raise YmfArchiveError(f"{path}: expected a YAML mapping at the top level")
    version = document.get("ymf_archive_version")
    if version != ARCHIVE_FORMAT_VERSION:
        raise YmfArchiveError(
            f"{path}: archive format version {version!r}, but this ymf writes "
            f"and reads version {ARCHIVE_FORMAT_VERSION}"
        )
    return document.get("domain", {}), document.get("extra")


def load_array(item: Dict[str, Any], base_dir: str | Path = "."):
    """The array a DataItem holds or points at, as a numpy array.

    Inline ``Values`` are reshaped to ``Dimensions``; an HDF5 reference
    ``"file.h5:/dataset"`` is read with the file relative to ``base_dir``
    (the archive's directory). Needs numpy, and h5py for references --
    imported here, so the archive core still needs neither.
    """
    import numpy
    if "Values" in item:
        dtype = "f%d" % item["Precision"] if item["DataType"] == "Float" else (
            "u%d" if item["DataType"] in ("UInt", "UChar") else "i%d") % item["Precision"]
        return numpy.array(item["Values"], dtype=dtype).reshape(item["Dimensions"])
    if item.get("Format") == "HDF" and "Data" in item:
        import h5py
        filename, dataset = str(item["Data"]).rsplit(":", 1)
        with h5py.File(Path(base_dir) / filename, "r") as f:
            return f[dataset][()]
    raise YmfArchiveError("load_array: no inline Values or HDF5 reference in %r" % (item,))


def _data_items(domain: Dict[str, Any]):
    """Every DataItem dict of a canonical domain, with a readable path."""
    for collection in domain.get("TimeCollections", []):
        for i, step in enumerate(collection["Data"]):
            for j, g in enumerate(step.get("SpatialCollection", [step])):
                where = "%s/%d%s" % (collection["Name"], i, "" if "Topology" in step else "/%d" % j)
                yield where + "/Topology", g["Topology"]
                yield where + "/Geometry", g["Geometry"]
                for a in g["Attributes"]:
                    yield "%s/%s/%s" % (where, a.get("Center", "Node"), a["Name"]), a


def inline_domain(domain: Dict[str, Any], base_dir: str | Path = ".") -> Dict[str, Any]:
    """A copy of ``domain`` with every referenced array copied in as ``Values``.

    The result needs no HDF5 file or sidecar: the ``.ymf`` (and the ``.xmf``
    derived from it) is self-contained.
    """
    import copy
    out = copy.deepcopy(canonicalize_domain(domain))
    for _, owner in _data_items(out):
        item = owner["DataItem"]
        if "Values" in item:
            continue
        array = load_array(item, base_dir)
        owner["DataItem"] = data_item_for(array, inline=True, dimensions=item["Dimensions"])
    return validate_domain(out)


def domain_arrays(domain: Dict[str, Any], base_dir: str | Path = "."):
    """{path: array} for every DataItem, however each is stored."""
    return {where: load_array(owner["DataItem"], base_dir)
            for where, owner in _data_items(canonicalize_domain(domain))}


def dump_grid(g: Dict[str, Any]) -> str:
    """Serialize one grid to a YAML string.

    For metadata that travels rather than lands in a file: a per-subdomain
    grid fragment stashed in an HDF5 dataset, or gathered across MPI. Pairs
    with :func:`load_grid`.
    """
    return _dump(g)


def load_grid(text: str) -> Dict[str, Any]:
    """Inverse of :func:`dump_grid`."""
    return _load(text)
