"""Geospatial/environmental data-sourcing helpers for YDMF.

Implements the "short-term integration" items from
``docs/data-sources-report.md`` §7:

- A land-cover -> Manning's n lookup table (Arcement & Schneider 1989 /
  Chow 1959 reference ranges, report §4.1), driving ``roughness_source``
  resolution for a YDMF discretization entry.
- A thin, dependency-light fetch helper for the named remote-sensing/DEM
  sources (SRTM, CHIRPS, ERA5, NLCD, ...) referenced by ``source``/
  ``source_url`` fields in ``mesh_generation``, ``initial_conditions``,
  ``boundary_conditions``, and ``roughness_source``.

**Scope note**: this sandbox has no ``numpy``/``rasterio``/``requests`` and
no outbound network access at the time of writing (see the units-report and
data-sources-report background sessions), so the actual raster I/O and
tile/service-specific download logic (a real "auto-fetch" per the report's
long-term recommendation) is **not implemented** here — that needs a real
geospatial stack (rasterio/xarray) and network access, and is left as a
documented extension point (:func:`fetch`, which is a thin
``urllib``-only downloader plus format dispatch, deliberately kept minimal).
What *is* fully implemented and tested: the Manning's n lookup (pure
stdlib, no deps) and YDMF-side metadata validation for source records
(``validate_source_record``).
"""

from __future__ import annotations

import csv
import io
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

__all__ = [
    "ManningLookupError",
    "DEFAULT_MANNING_N_TABLE",
    "load_manning_n_table",
    "manning_n_for_land_cover",
    "resolve_roughness_source",
    "KNOWN_DATA_SOURCES",
    "validate_source_record",
    "fetch",
]


class ManningLookupError(KeyError):
    """Raised when a land cover class has no Manning's n entry and no default."""


# Reference Manning's n ranges per report §4.1 (Chow 1959; Arcement &
# Schneider 1989 for the full methodology). Values are the midpoint of the
# cited typical range, in SI units of Manning's n (dimensionless-ish,
# s/m^(1/3), but per SI convention just given as bare n here).
DEFAULT_MANNING_N_TABLE: Dict[str, float] = {
    "open_water": 0.0325,
    "grassy_channel": 0.040,
    "wooded": 0.105,
    "urban_grass": 0.025,
    "urban_pavement": 0.0135,
    "residential": 0.0225,
    "cropland": 0.030,
    "forest_dense": 0.175,
    "desert_sand": 0.0325,
}

# NLCD (USGS National Land Cover Database) class code -> our internal
# land-cover key, for convenience when roughness_source.land_cover_data
# points at NLCD (report §4.2). Not exhaustive; covers the classes with an
# obvious Manning's n mapping. Falls back to KeyError -> ManningLookupError
# for anything unmapped, so callers know to extend rather than guess.
NLCD_CLASS_TO_LAND_COVER: Dict[int, str] = {
    11: "open_water",
    21: "urban_grass",  # developed, open space
    22: "residential",  # developed, low intensity
    23: "urban_pavement",  # developed, medium intensity
    24: "urban_pavement",  # developed, high intensity
    41: "forest_dense",  # deciduous forest
    42: "forest_dense",  # evergreen forest
    43: "forest_dense",  # mixed forest
    52: "wooded",  # shrub/scrub
    71: "grassy_channel",  # grassland/herbaceous
    81: "cropland",  # pasture/hay
    82: "cropland",  # cultivated crops
    90: "wooded",  # woody wetlands
    95: "grassy_channel",  # emergent herbaceous wetlands
}


def load_manning_n_table(path: Optional[str | Path] = None) -> Dict[str, float]:
    """Load a Manning's n lookup table.

    With no ``path``, returns :data:`DEFAULT_MANNING_N_TABLE`. With a path
    to a CSV (matching the ``lookup_table`` field on
    ``roughness_source``, report §3.6/§4.2), expects two columns
    ``land_cover,n`` (header optional, sniffed) and returns a merged table
    (CSV entries override defaults, so a partial CSV can extend defaults
    without repeating every row).
    """
    if path is None:
        return dict(DEFAULT_MANNING_N_TABLE)

    table = dict(DEFAULT_MANNING_N_TABLE)
    text = Path(path).read_text()
    reader = csv.reader(io.StringIO(text))
    for row in reader:
        if len(row) < 2:
            continue
        key, value = row[0].strip(), row[1].strip()
        if key.lower() in ("land_cover", "name", "class"):
            continue  # header row
        try:
            table[key] = float(value)
        except ValueError:
            continue
    return table


def manning_n_for_land_cover(
    land_cover: str,
    table: Optional[Dict[str, float]] = None,
    nlcd_class: Optional[int] = None,
) -> float:
    """Return Manning's n for a land cover key, or NLCD numeric class code.

    Exactly one of ``land_cover`` (internal key, e.g. ``"forest_dense"``) or
    ``nlcd_class`` (raw NLCD integer code, e.g. ``42``) should be meaningful;
    if both are given, ``nlcd_class`` is resolved first and ``land_cover``
    is only used as a fallback label for the error message.

    Raises :class:`ManningLookupError` if the class isn't in the table.
    """
    table = table if table is not None else DEFAULT_MANNING_N_TABLE
    key = land_cover
    if nlcd_class is not None:
        try:
            key = NLCD_CLASS_TO_LAND_COVER[nlcd_class]
        except KeyError as exc:
            raise ManningLookupError(
                f"no Manning's n mapping for NLCD class {nlcd_class!r}"
            ) from exc
    try:
        return table[key]
    except KeyError as exc:
        raise ManningLookupError(
            f"no Manning's n entry for land cover {key!r}"
        ) from exc


def resolve_roughness_source(
    roughness_source: Dict[str, Any],
    land_cover: str | int,
) -> float:
    """Resolve a YDMF ``roughness_source`` block (report §3.6) to a Manning's n value.

    ``land_cover`` is either an internal land-cover key (``str``) or a raw
    NLCD class code (``int``), matching the convention of
    ``roughness_source.land_cover_data.source`` (``"USGS NLCD"`` implies
    integer class codes; anything else implies internal string keys).
    """
    table = load_manning_n_table(roughness_source.get("lookup_table"))
    land_cover_data = roughness_source.get("land_cover_data") or {}
    source = (land_cover_data.get("source") or "").upper()

    if isinstance(land_cover, int) or "NLCD" in source:
        return manning_n_for_land_cover(str(land_cover), table, nlcd_class=int(land_cover))
    return manning_n_for_land_cover(str(land_cover), table)


@dataclass(frozen=True)
class DataSourceInfo:
    name: str
    category: str  # topography | bathymetry | vegetation | precipitation | wind | temperature
    resolution: str
    temporal: Optional[str]
    url: str
    license: str


# Catalog transcribed from data-sources-report.md so source names/URLs used
# in YDMF ``source``/``source_url`` fields can be validated against a known
# registry (report §7: "Add `url` field to reference data sources directly
# in the YAML"). Not exhaustive relative to the report — covers every
# source that has a URL cited in the report's appendix.
KNOWN_DATA_SOURCES: Dict[str, DataSourceInfo] = {
    "SRTM": DataSourceInfo(
        "SRTM", "topography", "30-90 m", None,
        "https://www.usgs.gov/centers/eros/data-archive/srtm-90m-digital-elevation-database-v4",
        "public_domain",
    ),
    "ASTER GDEM": DataSourceInfo(
        "ASTER GDEM", "topography", "30 m", None,
        "https://lpdaac.usgs.gov/", "free_registration",
    ),
    "ALOS AW3D": DataSourceInfo(
        "ALOS AW3D", "topography", "30 m (5 m select)", None,
        "https://www.eorc.jaxa.jp/ALOS/en/aw3d30/index.htm", "free",
    ),
    "Copernicus DEM": DataSourceInfo(
        "Copernicus DEM", "topography", "30-90 m", None,
        "https://spaceml.org/products/merder_dem", "free",
    ),
    "NASADEM": DataSourceInfo(
        "NASADEM", "topography", "30 m", None,
        "https://lpdaac.usgs.gov/products/nasademv001/", "free",
    ),
    "NOAA ETOPO": DataSourceInfo(
        "NOAA ETOPO", "bathymetry", "150 m - 2 km", None,
        "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO/", "free",
    ),
    "CHIRPS": DataSourceInfo(
        "CHIRPS", "precipitation", "5 km", "daily",
        "https://www.chc.ucsb.edu/data/chirps", "free",
    ),
    "GPM IMERG": DataSourceInfo(
        "GPM IMERG", "precipitation", "10 km", "30 min",
        "https://gpm1.gesdisc.eosdis.nasa.gov/", "free",
    ),
    "ERA5": DataSourceInfo(
        "ERA5", "wind", "25-30 km", "hourly",
        "https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels",
        "free",
    ),
    "MERRA-2": DataSourceInfo(
        "MERRA-2", "wind", "50 km", "hourly",
        "https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/", "free",
    ),
    "CCMP": DataSourceInfo(
        "CCMP", "wind", "25 km", "6-hourly",
        "https://podaac.jpl.nasa.gov/CCMP", "free",
    ),
    "WorldClim": DataSourceInfo(
        "WorldClim", "temperature", "1 km", "monthly",
        "https://www.worldclim.org/", "free",
    ),
    "CHELSA": DataSourceInfo(
        "CHELSA", "temperature", "1 km", "monthly",
        "https://chelsa-climate.org/", "free",
    ),
    "USGS NLCD": DataSourceInfo(
        "USGS NLCD", "vegetation", "30 m", None,
        "https://www.mrlc.gov/", "free",
    ),
    "ESA CCI": DataSourceInfo(
        "ESA CCI", "vegetation", "300 m", None,
        "https://www.esa-cci.org/", "free",
    ),
    "MODIS": DataSourceInfo(
        "MODIS", "vegetation", "250-500 m", "daily",
        "https://lpdaac.usgs.gov/", "free",
    ),
    "USGS NWIS": DataSourceInfo(
        "USGS NWIS", "hydrology", "station-based", "hourly",
        "https://waterdata.usgs.gov/nwis/", "free",
    ),
}


def validate_source_record(record: Dict[str, Any]) -> List[str]:
    """Validate a YDMF ``source``/``source_url`` pair against the known catalog.

    ``record`` is any dict with optional ``source`` and ``source_url`` keys
    (matches the shape of ``initial_conditions[]``, ``boundary_conditions[]``,
    and ``mesh_generation`` entries). Returns a list of warning strings
    (empty list = no issues). Never raises — an unrecognized source is a
    warning, not an error, since YDMF explicitly supports sources outside
    this catalog (report is a starting catalog, not an allowlist).
    """
    warnings: List[str] = []
    source = record.get("source")
    source_url = record.get("source_url")
    if not source:
        return warnings

    known = KNOWN_DATA_SOURCES.get(source)
    if known is None:
        warnings.append(
            f"source {source!r} is not in the known data source catalog "
            "(not an error — just unverified against docs/data-sources-report.md)"
        )
        return warnings

    if source_url and source_url != known.url:
        warnings.append(
            f"source_url {source_url!r} for {source!r} does not match the catalog URL "
            f"{known.url!r} (may be a mirror or a more specific endpoint — verify)"
        )
    return warnings


def fetch(url: str, dest: str | Path, timeout: float = 30.0) -> Path:
    """Download a URL to ``dest`` using only the standard library.

    Deliberately minimal: no retry/backoff policy, no format-specific
    parsing (GeoTIFF/HDF5/NetCDF decoding needs rasterio/xarray/h5py, none
    of which are assumed available — see module docstring). This exists so
    ``mesh_generation.source_url`` and similar fields have *something*
    runnable end-to-end without a heavy geospatial dependency; swap in a
    proper client (e.g. an authenticated Copernicus/CDS API client) when
    integrating a specific source.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:  # noqa: S310
            dest.write_bytes(response.read())
    except urllib.error.URLError as exc:
        raise ConnectionError(f"failed to fetch {url!r}: {exc}") from exc
    return dest
