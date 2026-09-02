"""YMF — YAML Modeling Format for computational physics."""

from importlib.metadata import PackageNotFoundError, version as _version

from ymf.schema import YMF_SCHEMA, load_ymf, validate_ymf
from ymf.normalize import normalize_unknowns
from ymf.units import check as check_units
from ymf.units import non_dimensionalize
from ymf.data_sources import resolve_roughness_source, validate_source_record
from ymf.xdmf import (
    build_xdmf_tree,
    write_xdmf,
    read_xdmf,
    round_trip_equal,
    canonicalize_domain,
)

__all__ = [
    "YMF_SCHEMA",
    "load_ymf",
    "validate_ymf",
    "normalize_unknowns",
    "check_units",
    "non_dimensionalize",
    "resolve_roughness_source",
    "validate_source_record",
    "build_xdmf_tree",
    "write_xdmf",
    "read_xdmf",
    "round_trip_equal",
    "canonicalize_domain",
]

try:
    __version__ = _version("ymf")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"
