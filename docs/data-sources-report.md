# Geospatial and Environmental Data Sources for Physical Modeling

A report on available data sources for topography, bathymetry, vegetation, roughness, and remote sensing — both free and commercial — relevant to physical modeling workflows.

---

## 1. Topography

### 1.1 Free Datasets

#### SRTM (Shuttle Radar Topography Mission)
- **Resolution**: 3 arcsec (~90 m), 1 arcsec (~30 m) in select areas
- **Coverage**: Global (±60° latitude)
- **Source**: USGS EarthExplorer, USGS Earthdata
- **Format**: GeoTIFF, SRTM1/2/3 DEM
- **Quality**: Good for large-scale modeling; gaps in tropical regions
- **License**: Public domain (for most uses)

#### ASTER GDEM (Advanced Spaceborne Thermal Emission and Reflection Radiometer Global Digital Elevation Model)
- **Resolution**: 30 m
- **Coverage**: Global (±60° latitude)
- **Source**: USGS EarthExplorer, NASA Earthdata
- **Format**: GeoTIFF, HDF5
- **Quality**: Better vertical accuracy than SRTM in some regions; noisy in others
- **License**: Free with registration

#### ALOS World 3D (AW3D30)
- **Resolution**: 30 m, 5 m (select areas)
- **Coverage**: Global
- **Source**: JAXA, Alaska Aerospace Dataset
- **Format**: GeoTIFF, ASTER L1a
- **Quality**: Better than SRTM; newer, more accurate
- **License**: Free

#### TanDEM-X DEM
- **Resolution**: 12 m
- **Coverage**: Global
- **Source**: DLR (German Aerospace Center)
- **Format**: GeoTIFF, HDF5
- **Quality**: Excellent; ~1–2 m vertical accuracy
- **License**: Free for non-commercial use

### 1.2 Commercial Datasets

#### DigitalGlobe/Maxar
- **Resolution**: < 0.5 m
- **Coverage**: Select areas worldwide
- **Format**: GeoTIFF, proprietary
- **Quality**: Very high; sub-meter accuracy
- **License**: Commercial, per-area

#### WorldDEM (NovarGIS/Trimble)
- **Resolution**: 12.5 m global, 5 m regional, 1 m local
- **Coverage**: Global
- **Quality**: High; ~2–3 m vertical accuracy
- **License**: Commercial

### 1.3 Specialized

#### Copernicus DEM (European Space Agency)
- **Resolution**: 30 m, 90 m
- **Coverage**: Global
- **Source**: ESA Copernicus Programme
- **License**: Free, no restrictions

#### NASADEM
- **Resolution**: 1 arcsec (~30 m)
- **Coverage**: Global
- **Source**: NASA/USGS
- **Quality**: Updated SRTM with improved processing
- **License**: Free

---

## 2. Bathymetry

### 2.1 Free Datasets

#### EMODnet (Europe)
- **Coverage**: European waters
- **Resolution**: Variable (100 m to 10 m)
- **Source**: European Commission
- **License**: Free

#### NOAA ETOPO
- **Resolution**: 1 arcmin (~2 km), 5 arcsec (~150 m)
- **Coverage**: Global
- **Source**: NOAA National Centers for Environmental Information
- **License**: Free

#### GMRT (General Bathymetric Chart of the Oceans)
- **Resolution**: Variable
- **Coverage**: Global oceans
- **Source**: IHO/IOC
- **License**: Free

#### SRTM30+
- **Resolution**: 30 arcsec (~900 m)
- **Coverage**: Global (includes bathymetry)
- **Source**: USGS/US Navy
- **License**: Free

### 2.2 Regional

#### NOAA NGS (National Geodetic Survey)
- **Coverage**: US waters and territories
- **Resolution**: High (sub-meter in some areas)
- **License**: Free

#### UK Hydrographic Office (UKHO)
- **Coverage**: UK waters, select international waters
- **Resolution**: High
- **License**: Free for many applications

---

## 3. Vegetation and Ecosystem Type

### 3.1 Free Datasets

#### MODIS Vegetation Index
- **Resolution**: 250 m (MOD13Q1), 500 m (MOD13A3)
- **Coverage**: Global, daily
- **Source**: NASA LP DAAC
- **License**: Free
- **Product**: NDVI, EVI, LAI, FVC

#### ESA CCI Land Cover
- **Resolution**: 300 m
- **Coverage**: Global
- **Source**: ESA Climate Change Initiative
- **License**: Free

#### GlobCover
- **Resolution**: 300 m
- **Coverage**: Global
- **Source**: ESA/ESA-CCI
- **License**: Free

#### USGS NLCD (National Land Cover Database)
- **Resolution**: 30 m
- **Coverage**: US only
- **Source**: USGS
- **License**: Free (for most uses)

#### LandScan Population
- **Resolution**: 1 km
- **Coverage**: Global
- **Source**: Oak Ridge National Laboratory
- **License**: Free for research

### 3.2 Commercial Datasets

#### ESRI Land Cover
- **Resolution**: Variable (30 m global)
- **Coverage**: Global
- **License**: Commercial

#### Maxar WorldView
- **Resolution**: 0.3–0.5 m
- **Coverage**: Global
- **License**: Commercial

---

## 4. Manning/Chezy Roughness Coefficients

### 4.1 Free Datasets

#### Land Cover → Manning's n Mapping

This is the most common approach: derive roughness from land cover data.

| Land Cover Type | Manning's n (typical range) | Reference |
|---|---|---|
| Open water | 0.025–0.040 | Chow, 1959 |
| Grassy channels | 0.030–0.050 | Chow, 1959 |
| Wooded areas | 0.060–0.150 | Chow, 1959 |
| Urban grass | 0.020–0.030 | Chow, 1959 |
| Urban pavement | 0.012–0.015 | Chow, 1959 |
| Residential areas | 0.015–0.030 | Chow, 1959 |
| Cropland | 0.020–0.040 | Chow, 1959 |
| Forest (dense) | 0.100–0.250 | Chow, 1959 |
| Desert/sand | 0.025–0.040 | Chow, 1959 |

**Sources for Manning's n reference values**:
- Chow (1959) — "Open Channel Hydraulics"
- Barnes (1967) — "Roughness Characteristics of Natural Channels"
- Arcement and Schneider (1989) — "Guide for Selecting Manning's Roughness Coefficients"
- USGS Technical Report 41

### 4.2 Derived Datasets

#### USGS National Map (topographic)
- Can derive roughness from terrain data
- Combines with land cover for more accurate Manning's n
- **Resolution**: 1/9 arcsec (~3 m) for Lidar areas

#### ESA CCI Land Cover + Manning's n lookup table
- Automatically maps land cover classes to Manning's n ranges
- **Resolution**: 300 m (from land cover source)

#### USGS NLCD + Manning's n lookup table
- Higher resolution (30 m) than CCI
- US-specific only
- **Resolution**: 30 m

### 4.3 Commercial Datasets

#### ESRI Land Cover + Manning's n mapping
- Commercial land cover data mapped to roughness coefficients
- **Resolution**: 30 m

#### Topobathy.com / Bathymetry.com
- Can derive roughness from bathymetry data
- **License**: Commercial

---

## 5. Remote Sensing — Precipitation, Wind, Temperature

### 5.1 Precipitation

#### FREE

**GPM (Global Precipitation Measurement) IMERG**
- **Resolution**: 0.1° × 0.1° (~10 km)
- **Temporal**: 30-minute intervals, 3-hourly final
- **Coverage**: Global
- **Source**: NASA
- **License**: Free

**CHIRPS (Climate Hazards Group InfraRed Precipitation with Station)**
- **Resolution**: 0.05° (~5 km)
- **Temporal**: Daily, 30+ years
- **Coverage**: Global (60°S–60°N)
- **Source**: UC Santa Barbara Climate Hazard Group
- **License**: Free (for most uses)

**TRMM (Tropical Rainfall Measuring Mission) 3B42**
- **Resolution**: 0.25° (~25 km)
- **Temporal**: 3-hourly
- **Coverage**: Global (50°S–50°N)
- **Source**: NASA/GSFC
- **License**: Free

**ERA5-Land (ECMWF Reanalysis)**
- **Resolution**: 0.1° (~9 km)
- **Temporal**: Hourly, 40+ years
- **Coverage**: Global
- **Source**: ECMWF
- **License**: Free

**NOAA CPC Morphing Technique (CMORPH)**
- **Resolution**: 8 km
- **Temporal**: 30-minute intervals
- **Coverage**: Global
- **Source**: NOAA
- **License**: Free

**CMFD (China Meteorological Forcing Dataset)**
- **Resolution**: 0.1° × 0.1°
- **Temporal**: 3-hourly, 40+ years
- **Coverage**: China and surrounding areas
- **Source**: Chinese Academy of Sciences
- **License**: Free

#### COMMERCIAL

**QuantArc**
- **Resolution**: High (sub-kilometer)
- **Temporal**: Near real-time
- **License**: Commercial

**RainViewer**
- **Resolution**: Radar-based
- **Coverage**: Europe, select regions
- **License**: Commercial/free tier

### 5.2 Wind Speed and Direction

#### FREE

**ERA5 (ECMWF Reanalysis v5)**
- **Resolution**: 0.25° (~30 km)
- **Temporal**: Hourly, 40+ years
- **Coverage**: Global
- **Source**: ECMWF
- **License**: Free (via Copernicus Climate Data Store)
- **Variables**: 10m wind speed, 10m wind direction, surface wind, geostrophic wind, etc.

**MERRA-2 (Modern-Era Retrospective Analysis)**
- **Resolution**: 0.5° × 0.625° (~50 km)
- **Temporal**: Hourly
- **Coverage**: Global
- **Source**: NASA/GSFC
- **License**: Free

**CCMP (Cross-Calibrated Multi-Platform)**
- **Resolution**: 0.25° (~25 km)
- **Temporal**: 6-hourly
- **Coverage**: Global ocean
- **Source**: NASA/JPL
- **License**: Free

**NCEP/NCAR Reanalysis**
- **Resolution**: 2.5° × 2.5° (~250 km)
- **Temporal**: 6-hourly, 40+ years
- **Coverage**: Global
- **Source**: NOAA/ESRL
- **License**: Free

#### COMMERCIAL

**Hindcast Design**
- **Resolution**: 0.25°–1 km (model-dependent)
- **Temporal**: Historical to near real-time
- **License**: Commercial

### 5.3 Temperature

#### FREE

**ERA5**
- **Resolution**: 0.25° × 0.25°
- **Temporal**: Hourly
- **Coverage**: Global
- **License**: Free

**NCEP/NCAR Reanalysis**
- **Resolution**: 2.5°
- **Temporal**: 6-hourly
- **Coverage**: Global
- **License**: Free

**GMFD (Global Meteorological Forcing Dataset)**
- **Resolution**: 0.5° × 0.5°
- **Temporal**: Daily, 30+ years
- **Coverage**: Global
- **License**: Free

**WorldClim**
- **Resolution**: 30 arcsec (~1 km)
- **Temporal**: Monthly averages
- **Coverage**: Global
- **Source**: WorldClim
- **License**: Free

**CHELSA (Climatologies at High Resolution for the Earth's Land Surface Areas)**
- **Resolution**: 30 arcsec (~1 km)
- **Temporal**: Monthly, 1979–present
- **Coverage**: Global (land)
- **License**: Free

---

## 6. Integration with YMF

### 6.1 Remote Sensing Fields as Input

YMF can include remote sensing fields as initial/boundary conditions:

```yaml
strong_form:
  initial_conditions:
    - field: precipitation
      source: "CHIRPS"
      temporal: "daily"
      source_url: "https://www.chc.ucsb.edu/data/chirps"
    - field: temperature
      source: "ERA5"
      temporal: "hourly"
      source_url: "https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels"
    - field: wind_speed
      source: "ERA5"
      temporal: "hourly"

  boundary_conditions:
    - region: "upstream"
      type: "dirichlet"
      variable: "discharge"
      source: "USGS NLDI"   # National Hydrologic Data for the US
      source_url: "https://www.usgs.gov/national-hydrography/"
    - region: "downstream"
      type: "neumann"
      variable: "stage"
      value: 0
```

### 6.2 Land Cover → Roughness Mapping

```yaml
discretizations:
  - name: "FD_with_manning"
    roughness_source:
      type: "land_cover_mapping"
      land_cover_data:
        source: "USGS NLCD"   # or "ESA CCI" for global
        resolution: "30 m"
        url: "https://www.mrlc.gov/"
      reference: "Arcement and Schneider (1989)"
      lookup_table: "manning_n_lookup.csv"   # embedded or external
    discretization:
      # ...
```

### 6.3 Topography/Bathymetry as Mesh Generation

```yaml
Problem:
  strong_form:
    domain: |
      Ω defined from SRTM DEM
    mesh_generation:
      source: "SRTM 1arcsec"
      source_url: "https://www.usgs.gov/centers/eros/data-archive/srtm-90m-digital-elevation-database-v4"
      format: "GeoTIFF"
      vertical_datums:
        horizontal: "WGS84"
        vertical: "EGM96"
      resolution_target: 30   # meters (from DEM resolution)
      filter:
        method: "smoothing"
        kernel: "gaussian"
        sigma: 2
```

### 6.4 Environmental Forcing as Time-Varying BCs

```yaml
discretizations:
  - name: "storm_flood_simulation"
    boundary_conditions:
      upstream:
        type: "hydrograph"
        source: "USGS NWIS"   # National Water Information System
        source_url: "https://waterdata.usgs.gov/nwis/"
        temporal_resolution: "hourly"
        variables: ["discharge", "stage", "turbidity"]
      atmospheric:
        precipitation:
          source: "CHIRPS"
          temporal: "daily"
        wind:
          source: "ERA5"
          temporal: "hourly"
        temperature:
          source: "CHELSA"
          temporal: "monthly"
```

### 6.5 Data Accessibility Summary

| Data Type | Primary Free Source | Primary Commercial Source | YMF Integration |
|---|---|---|---|
| Topography | SRTM, Copernicus DEM, NASADEM | Maxar, WorldDEM | Mesh generation, initial conditions |
| Bathymetry | NOAA ETOPO, EMODnet | Seabed 2030 consortium | Domain definition, initial conditions |
| Vegetation | MODIS, NLCD, ESA CCI | ESRI Land Cover | Roughness mapping, land cover BCs |
| Manning's n | Derived from land cover | Commercial land cover → n mapping | BC parameter |
| Precipitation | GPM IMERG, CHIRPS | QuantArc | Time-varying BCs, initial conditions |
| Wind | ERA5, CCMP, NCEP | Hindcast Design | Atmospheric BCs, surface stress |
| Temperature | ERA5, WorldClim, CHELSA | DWD | Atmospheric BCs, thermal BCs |

---

## 7. Recommendations

### For YMF Schema (immediate)

1. **Add `data_source` field** to initial conditions, boundary conditions, and domain specification
2. **Add `url` field** to reference data sources directly in the YAML
3. **Add `temporal` field** to specify the time resolution of remote sensing data
4. **Add `source_url` field** to point to the data source for reproducibility

### For Integration (short-term)

1. **Create `ymf.data_sources` module** that can load and process topography, bathymetry, and land cover data from YMF-specified sources
2. **Add land cover → Manning's n lookup table** (standardized reference: Arcement & Schneider, 1989)
3. **Create mesh generation from DEM data** using YMF's domain specification
4. **Add environmental forcing** (precipitation, wind, temperature) as time-varying BCs

### For Long-Term

1. **Automated data sourcing**: YMF automatically downloads and validates data from specified sources
2. **Multi-resolution support**: YMF handles coarse (ERA5) to fine (NLCD) data seamlessly
3. **Temporal consistency checks**: Ensure temporal resolution of different forcing data is compatible
4. **Spatial registration**: Ensure topography, bathymetry, land cover, and weather data are all in the same coordinate reference system

---

## Appendix: Key Data Source URLs

| Source | URL |
|---|---|
| USGS EarthExplorer | https://earthexplorer.usgs.gov/ |
| USGS National Map | https://apps.nationalmap.gov/downloader/ |
| Copernicus DEM | https://spaceml.org/products/merder_dem |
| NOAA ETOPO | https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO/ |
| CHIRPS | https://www.chc.ucsb.edu/data/chirps |
| GPM IMERG | https://gpm1.gesdisc.eosdis.nasa.gov/ |
| ERA5 | https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels |
| MERRA-2 | https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/ |
| CCMP | https://podaac.jpl.nasa.gov/CCMP |
| MODIS | https://lpdaac.usgs.gov/ |
| NLCD | https://www.mrlc.gov/ |
| ESA CCI | https://www.esa-cci.org/ |
| WorldClim | https://www.worldclim.org/ |
| CHELSA | https://chelsa-climate.org/ |
| SRTM | https://www.usgs.gov/centers/eros/data-archive/srtm-90m-digital-elevation-database-v4 |
| NASADEM | https://lpdaac.usgs.gov/products/nasademv001/ |
| ALOS AW3D | https://www.eorc.jaxa.jp/ALOS/en/aw3d30/index.htm |
