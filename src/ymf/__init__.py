"""YMF — YAML Modeling Format for computational physics.

Importing this package pulls in the **core** only: :mod:`ymf.archive` (the
archive data model and YAML I/O, needing just ``pyyaml``),
:mod:`ymf.xdmf` (the ``.xmf`` converter, stdlib-only), and
:mod:`ymf.normalize` (stdlib-only).

The validation stack is *not* imported eagerly. ``ymf.schema`` needs
``strictyaml`` (``pip install ymf[spec]``) and ``ymf.units`` needs ``pint``
(``ymf[units]``); both are resolved on first attribute access via PEP 562,
so a solver that only writes archives never loads them and never needs them
installed. That is deliberate and load-bearing -- see the trust-boundary
section of :mod:`ymf.archive` -- so please don't turn these back into
module-level imports for tidiness.
"""

from importlib.metadata import PackageNotFoundError, version as _version
from typing import Any

from ymf.archive import (
    ARCHIVE_FORMAT_VERSION,
    YmfArchiveError,
    add_spatial_step,
    add_uniform_step,
    attribute,
    canonicalize_domain,
    check_dimensions,
    data_item,
    data_item_for,
    dump_grid,
    geometry,
    grid,
    load_grid,
    new_domain,
    read_ymf,
    require_libyaml,
    topology,
    validate_domain,
    write_ymf,
)
from ymf.normalize import normalize_unknowns
from ymf.xdmf import (
    build_xdmf_tree,
    grid_element_time,
    parse_grid_element,
    read_xdmf,
    round_trip_equal,
    write_xdmf,
)

#: Attribute name -> (module, name within that module). Resolved lazily by
#: :func:`__getattr__` so that the heavy dependencies behind these modules
#: are only imported if something actually asks for them.
_LAZY_ATTRS = {
    # ymf.schema -- requires strictyaml (ymf[spec])
    "YMF_SCHEMA": ("ymf.schema", "YMF_SCHEMA"),
    "load_ymf": ("ymf.schema", "load_ymf"),
    "validate_ymf": ("ymf.schema", "validate_ymf"),
    # ymf.units -- requires pint (ymf[units])
    "check_units": ("ymf.units", "check"),
    "non_dimensionalize": ("ymf.units", "non_dimensionalize"),
    # ymf.data_sources -- stdlib only, but a front-end concern
    "resolve_roughness_source": ("ymf.data_sources", "resolve_roughness_source"),
    "validate_source_record": ("ymf.data_sources", "validate_source_record"),
}

#: Which extra to suggest when a lazy import fails for a missing dependency.
_EXTRA_FOR_MODULE = {"ymf.schema": "spec", "ymf.units": "units"}

__all__ = [
    # --- core: archive model and YAML I/O ---
    "ARCHIVE_FORMAT_VERSION",
    "YmfArchiveError",
    "new_domain",
    "add_uniform_step",
    "add_spatial_step",
    "grid",
    "topology",
    "geometry",
    "attribute",
    "data_item",
    "data_item_for",
    "canonicalize_domain",
    "validate_domain",
    "check_dimensions",
    "write_ymf",
    "read_ymf",
    "dump_grid",
    "load_grid",
    "require_libyaml",
    # --- core: XDMF conversion ---
    "build_xdmf_tree",
    "write_xdmf",
    "read_xdmf",
    "round_trip_equal",
    "parse_grid_element",
    "grid_element_time",
    # --- core: stdlib helpers ---
    "normalize_unknowns",
    # --- lazy: front-end validation (see _LAZY_ATTRS) ---
    *_LAZY_ATTRS,
]


def __getattr__(name: str) -> Any:
    """Resolve the front-end attributes in :data:`_LAZY_ATTRS` on demand."""
    target = _LAZY_ATTRS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = target
    import importlib

    try:
        module = importlib.import_module(module_name)
    except ImportError as exc:
        extra = _EXTRA_FOR_MODULE.get(module_name)
        if extra is None:
            raise
        raise ImportError(
            f"ymf.{name} lives in {module_name}, which needs an optional "
            f"dependency that isn't installed. Install it with: "
            f"pip install 'ymf[{extra}]'. (The YMF core deliberately depends "
            f"on pyyaml alone, so a solver writing archives doesn't have to "
            f"carry the validation stack.)"
        ) from exc
    value = getattr(module, attr_name)
    globals()[name] = value  # cache, so this path runs at most once per name
    return value


def __dir__() -> list:
    return sorted(__all__)


try:
    __version__ = _version("ymf")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"
