"""YDMF — YAML Data Model and Format for computational physics."""

from importlib.metadata import PackageNotFoundError, version as _version

from ydmf.schema import YDMF_SCHEMA, load_ydmf, validate_ydmf
from ydmf.normalize import normalize_unknowns
from ydmf.units import check as check_units
from ydmf.units import non_dimensionalize
from ydmf.data_sources import resolve_roughness_source, validate_source_record
from ydmf.xdmf import (
    build_xdmf_tree,
    write_xdmf,
    read_xdmf,
    round_trip_equal,
    canonicalize_domain,
)

__all__ = [
    "YDMF_SCHEMA",
    "load_ydmf",
    "validate_ydmf",
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
    __version__ = _version("ydmf")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"
