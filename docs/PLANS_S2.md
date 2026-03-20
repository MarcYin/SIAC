# SIAC v2 — Sentinel-2 Sensor Plan

This document covers the **Sentinel-2–specific** design for SIAC v2. It implements and
extends the sensor-agnostic core architecture defined in [PLANS.md](PLANS.md).

All Sentinel-2 work follows the core plan's modular contracts:
- M1: `Sentinel2Preprocessor` class → returns `ObservationBundle`
- M3: S2-specific surface prior options → returns `SurfacePrior`
- `S2DataBackend` protocol for data access (search + download) — sits before M1
- `SensorBand` / `SensorConfig` for band definitions
- AOI-scoped data access rules
- `sensor_to_reference()` / `reference_to_sensor()` functions for reference-band mapping

When adding support for a new sensor, use this document as a **reference
implementation** — the structure here mirrors what any sensor-specific plan
should cover (see [Core Requirements for Adding New Sensors](PLANS.md#8-core-requirements-for-adding-new-sensors) in the core plan).

---

## Table of Contents

1. [S2 Module Mapping (How S2 Implements Core Contracts)](#1-s2-module-mapping)
2. [Sentinel-2 Data Access Module](#2-sentinel-2-data-access-module)
3. [S2 Spectral Band Definitions](#3-s2-spectral-band-definitions)
4. [S2 Preprocessor (M1) — ObservationBundle Producer](#4-s2-preprocessor-m1--observationbundle-producer)
5. [S2 Multi-Resolution Handling (M4 / M6 Concerns)](#5-s2-multi-resolution-handling-m4--m6-concerns)
6. [S2-Specific Custom Overrides](#6-s2-specific-custom-overrides)
7. [Files Summary](#7-files-summary)
8. [Usage Examples](#8-usage-examples)
9. [Verification Plan](#9-verification-plan)

---

## 1. S2 Module Mapping

This section maps the core plan's six modules (M1–M6) to their S2-specific implementations, showing what is sensor-specific and what is reused from core.

| Core Module | S2 Implementation | Output Contract | S2-Specific? |
|-------------|-------------------|-----------------|-------------|
| **M1** Preprocessor | `Sentinel2Preprocessor` class | `ObservationBundle` | **Yes** — reads SAFE format, handles S2 angles, multi-res bands |
| **M2** Atmo Provider | `CAMSProvider` / `MERRA2Provider` (generic) | `AtmosphericState` | **No** — sensor-agnostic, reused as-is |
| **M3** Surface Prior | `BRDFDerivedPriorProvider` or `PrebuiltPriorStore` | `SurfacePrior` | Partially — needs S2 SensorConfig for band convolution |
| **M4** Grid Assembler | Generic `assemble_grids()` with S2 multi-res awareness | `SolverInputBundle` | Partially — S2 has 10/20/60 m bands that need special resampling |
| **M5** Solver | `solve_aerosol()` (generic) | `SolvedAtmosphere` | **No** — sensor-agnostic |
| **M6** Corrector | `correct_atmosphere()` (generic) with S2 multi-res upsample | `CorrectionResult` | Partially — correction at native 10/20/60 m per-band |

**Key insight**: only M1 is fully S2-specific. M4 and M6 have minor S2 awareness (multi-resolution), but the contract types are identical to the core.

---

## 2. Sentinel-2 Data Access Module

**Status**: Partially implemented (core data-access path now wired)

Implemented in code:
- `S2Query`, `S2Product`, `S2DataAccess` orchestration (`python/siac/io/s2_data_source.py`)
- `CopernicusDataspaceBackend` (`python/siac/io/copernicus_dataspace.py`)
- `GCSSentinel2Backend` with real listing/download + retries (`python/siac/io/gcs_sentinel2.py`)
- `S2DataAccessConfig` in `SIACConfig` (`python/siac/core/config.py`)
- S2 convenience entrypoints `resolve_s2_input()`, `search_sentinel2()`, `siac_process_s2()` (`python/siac/siac.py`)

### Context

The current `preprocess_s2()` reads **local SAFE directories only** — there's no discovery/download step. Users must manually obtain S2 data before running SIAC. This plan adds a **data access layer** that searches, selects, and fetches S2 L1C products, then hands local paths to M1.

This module sits **before M1** in the pipeline — it resolves a user query into a local path, then M1 produces the `ObservationBundle`.

**AOI interaction**: If an AOI is provided, the downloaded S2 TOA data will be cropped to the AOI extent by M1 during `preprocess()`. If no AOI is provided, the full S2 tile extent becomes the AOI. All other modules (M2, M3) only receive the resolved bounds.

### Data Sources

1. **Copernicus Data Space Ecosystem (CDSE)** — S3-compatible API: https://documentation.dataspace.copernicus.eu/APIs/S3.html
2. **Google Cloud Storage public dataset** — `gs://gcp-public-data-sentinel-2` (no auth, MGRS-based paths)

#### GCS Path Convention

```python
mgrs_tile = safe_id.split('_')[-2][1:]
mgrs_tile_path = f"{mgrs_tile[:2]}/{mgrs_tile[2:3]}/{mgrs_tile[3:]}"
url = f"gs://gcp-public-data-sentinel-2/tiles/{mgrs_tile_path}/{safe_id}.SAFE"
```

### Input Modes

| Mode | Input Example | Behavior |
|------|--------------|----------|
| **Product ID** | `S2B_MSIL1C_20210801T103629_N0500_R008_T31UDQ_20230731T120241.SAFE` | Direct fetch, no search |
| **MGRS + date shorthand** | `T31UDQ_20210801` | Search for tile+date, pick newest baseline |
| **MGRS tile + time range** | `tile="31UDQ", start="2021-08-01", end="2021-08-05"` | Search by tile+temporal window |
| **Location + time range** | `bbox=(1.5, 50.0, 2.5, 51.0)` or GeoJSON/AOI | Spatial+temporal search |
| **All above** with filters | `+ cloud_cover_max=30` | Filter results by cloud % |

### Processing Baseline Gotcha

Multiple copies of the same tile/date can exist with different processing baselines (N0301, N0400, N0500...). The module **always selects the newest baseline** (highest `N####` number) when duplicates exist.

### Architecture

```
siac/
  io/
    s2_data_source.py          ← NEW: search + download orchestrator
    copernicus_dataspace.py    ← NEW: CDSE S3 backend
    gcs_sentinel2.py           ← NEW: GCS public bucket backend
  satellite/
    sentinel2.py               ← EXISTING (no changes needed)
  core/
    config.py                  ← MODIFY: add S2DataAccessConfig
```

#### Key Abstraction

```python
class S2DataBackend(Protocol):
    """Backend for fetching S2 SAFE directories."""
    def search(self, query: S2Query) -> list[S2Product]: ...
    def download(self, product: S2Product, dest_dir: Path) -> Path: ...
```

Two concrete backends:
1. **`CopernicusDataspaceBackend`** — CDSE S3 (primary, full catalog search)
2. **`GCSSentinel2Backend`** — GCS public bucket (no auth, MGRS-based listing only)

### Phase 1: Data Types & Query Model

#### New file: `python/siac/io/s2_data_source.py`

```python
@dataclass
class S2Query:
    """Flexible query for S2 products."""
    product_id: str | None = None           # Direct product ID
    mgrs_tile: str | None = None            # e.g. "31UDQ"
    date: date | None = None                # Single date (for tile+date shorthand)
    start_date: date | None = None          # Time range start
    end_date: date | None = None            # Time range end
    bbox: tuple[float,float,float,float] | None = None  # (W,S,E,N) WGS84
    aoi: AOI | None = None                  # From existing AOI class
    max_cloud_cover: float = 100.0          # Max cloud cover %
    processing_level: str = "L1C"           # L1C or L2A

    @classmethod
    def from_product_id(cls, product_id: str) -> S2Query: ...

    @classmethod
    def from_tile_date(cls, tile_date: str) -> S2Query:
        """Parse 'T31UDQ_20210801' or '31UDQ_20210801'."""
        ...

    def validate(self) -> None:
        """Ensure at least one spatial constraint is set."""
        ...

@dataclass
class S2Product:
    """Metadata for a discovered S2 product."""
    product_id: str                          # Full SAFE name
    mgrs_tile: str                           # e.g. "31UDQ"
    sensing_date: datetime
    processing_baseline: str                 # e.g. "N0500"
    cloud_cover: float                       # 0-100
    satellite: str                           # "S2A" or "S2B"
    orbit_number: int
    source_url: str                          # Backend-specific URI
    size_mb: float | None = None

    @property
    def baseline_number(self) -> int:
        """Extract numeric baseline for comparison, e.g. 'N0500' → 500."""
        return int(self.processing_baseline.replace("N", ""))


class S2DataAccess:
    """
    Unified S2 data access — search, select, download.

    This is one of the few classes in the codebase — justified because it
    holds a configured backend + cache directory (genuinely coupled state).
    All public methods are thin wrappers that delegate to plain functions.
    """

    def __init__(self, backend: S2DataBackend, cache_dir: Path | None = None):
        self._backend = backend
        self._cache_dir = cache_dir or Path.home() / ".cache" / "siac" / "s2"

    def get(self, query: S2Query | str, dest_dir: Path | None = None) -> Path:
        """Main entry point. Returns path to local SAFE directory (input to M1)."""
        return fetch_s2(self._backend, query, dest_dir or self._cache_dir)

    def search(self, query: S2Query | str) -> list[S2Product]:
        """Search only (no download). Returns deduplicated list."""
        return search_s2(self._backend, query)


# ── Plain helper functions (the real logic) ────────────────────────────

def fetch_s2(backend: S2DataBackend, query: S2Query | str, dest_dir: Path) -> Path:
    """Search, select best product, download. Returns local SAFE path."""
    ...

def search_s2(backend: S2DataBackend, query: S2Query | str) -> list[S2Product]:
    """Search and deduplicate. No download."""
    ...

def deduplicate_products(products: list[S2Product]) -> list[S2Product]:
    """Group by (mgrs_tile, sensing_date), keep highest baseline_number."""
    ...

def select_best_product(products: list[S2Product]) -> S2Product:
    """Pick single best product (newest baseline, then newest sensing)."""
    ...
```

### Phase 2: Copernicus Data Space Ecosystem Backend

#### New file: `python/siac/io/copernicus_dataspace.py`

Uses CDSE's **S3-compatible API** via `boto3`:

```python
# Constants
CDSE_ENDPOINT = "https://eodata.dataspace.copernicus.eu"
CDSE_BUCKET = "eodata"
CDSE_PREFIX = "Sentinel-2/MSI/L1C"
CDSE_ODATA_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1"

def search_cdse(query: S2Query, access_key=None, secret_key=None) -> list[S2Product]:
    """Search CDSE OData catalog for S2 products."""
    ...

def download_cdse(product: S2Product, dest_dir: Path, access_key=None, secret_key=None) -> Path:
    """Download S2 SAFE dir from CDSE S3 bucket."""
    ...
```

**Why CDSE S3 over HTTP zip download**: S3 allows downloading individual bands/files (for future AOI-scoped partial reads), supports resume, and avoids the need to unzip. The OData catalog API is free and doesn't require tokens for searching.

### Phase 3: Google Cloud Storage Backend

#### New file: `python/siac/io/gcs_sentinel2.py`

```python
# Constants
GCS_BUCKET = "gcp-public-data-sentinel-2"

def search_gcs(query: S2Query) -> list[S2Product]:
    """List S2 products via GCS prefix listing (MGRS tile or product ID required)."""
    ...

def download_gcs(product: S2Product, dest_dir: Path) -> Path:
    """Download S2 SAFE dir from GCS public bucket."""
    ...
```

**Trade-offs**: GCS is public (no auth), fast, but has no catalog search — only prefix listing. Good when you know the MGRS tile. CDSE has full catalog search but requires credentials.

**S2DataBackend protocol** (kept as a protocol because backends bundle search + download):
```python
class S2DataBackend(Protocol):
    def search(self, query: S2Query) -> list[S2Product]: ...
    def download(self, product: S2Product, dest_dir: Path) -> Path: ...
```

The `search_cdse`/`download_cdse` and `search_gcs`/`download_gcs` plain functions
are wrapped into backend objects only when passed to `S2DataAccess`. Users working
directly with functions never need the protocol.

### Phase 4: Config Integration

#### Modify: `python/siac/core/config.py`

```python
class S2DataAccessConfig(BaseModel):
    """Configuration for Sentinel-2 data access."""

    backend: Literal["cdse", "gcs", "local"] = Field(
        default="local",
        description="Data source: 'cdse' (Copernicus DataSpace), 'gcs' (Google Cloud), 'local'",
    )
    cache_dir: Path | None = Field(default=None)
    max_cloud_cover: float = Field(default=80.0, ge=0, le=100)
    prefer_newest_baseline: bool = Field(default=True)
    cdse_access_key: str | None = None
    cdse_secret_key: str | None = None
```

Add to `SIACConfig`:
```python
s2_data: S2DataAccessConfig = Field(
    default_factory=S2DataAccessConfig,
    description="Sentinel-2 data access configuration",
)
```

### Phase 5: Pipeline Integration

#### Modify: `python/siac/siac.py`

```python
def siac_process_s2(
    config: SIACConfig,
    query: str | S2Query | Path,
    *,
    output_path: Path | None = None,
    **overrides,
) -> CorrectionResult:
    """Process S2 with flexible input (path, product ID, tile+date, or query).

    Resolves the query to a local SAFE path, then calls siac_process()
    with the standard pipeline — the same contracts apply.
    """
    input_path = resolve_s2_input(query, config)
    return siac_process(config, input_path, **overrides)


def resolve_s2_input(query: str | S2Query | Path, config: SIACConfig) -> Path:
    """Resolve flexible query to a local SAFE directory path."""
    ...
```

New convenience functions:
```python
def search_sentinel2(tile=None, date=None, start_date=None, end_date=None,
                     bbox=None, max_cloud_cover=80.0, backend="cdse") -> list[S2Product]:
    """Search for available S2 products without downloading."""
    ...
```

---

## 3. S2 Spectral Band Definitions

Sentinel-2 MSI has 13 bands across three native resolutions. The following
`SensorBand` definitions map S2 into the core plan's sensor-agnostic
spectral model (see [Core Plan §9.2](PLANS.md#92-spectral-band-descriptor)).

Authoritative RSRF source:
- [SentiWiki S2 Mission](https://sentiwiki.copernicus.eu/web/s2-mission)
- the implementation should follow the linked `Sentinel-2 Spectral Response Functions (S2-SRF)` document from that page / its linked documents page, rather than hard-coding a transient attachment URL

### S2A/S2B Band Table

| Band | Name | Center λ (nm) S2A | Center λ (nm) S2B | FWHM (nm) | Resolution (m) | Spectral Region | Role in SIAC |
|------|------|--------------------|--------------------|-----------|----------------|-----------------|-------------|
| B01 | Coastal | 443.9 | 442.3 | 27 | 60 | Coastal/Deep blue | Aerosol retrieval (high AOT sensitivity) |
| B02 | Blue | 496.6 | 492.1 | 98 | 10 | Blue | Aerosol retrieval (primary) |
| B03 | Green | 560.0 | 559.0 | 45 | 10 | Green | Surface prior target |
| B04 | Red | 664.5 | 665.0 | 38 | 10 | Red | Surface prior target |
| B05 | VRE-1 | 703.9 | 703.8 | 19 | 20 | Red edge | Vegetation |
| B06 | VRE-2 | 740.2 | 739.1 | 18 | 20 | Red edge | Vegetation |
| B07 | VRE-3 | 782.5 | 779.7 | 28 | 20 | NIR | Vegetation |
| B08 | NIR | 835.1 | 833.0 | 145 | 10 | NIR | Broad NIR |
| B8A | NIR-n | 864.8 | 864.0 | 33 | 20 | NIR | Climatology query (narrow) |
| B09 | WV | 945.0 | 943.2 | 26 | 60 | NIR/WV | Water vapour |
| B10 | Cirrus | 1373.5 | 1376.9 | 75 | 60 | Cirrus | Cirrus detection (excluded from surface) |
| B11 | SWIR-1 | 1613.7 | 1610.4 | 143 | 20 | SWIR-1 | Climatology query (aerosol-transparent) |
| B12 | SWIR-2 | 2202.4 | 2185.7 | 242 | 20 | SWIR-2 | Climatology query (aerosol-transparent) |

### SensorConfig Construction

The S2 preprocessor (M1) builds a `SensorConfig` from the SAFE metadata
(specifically from `MTD_MSIL1C.xml`). This maps directly into the core plan's
`SensorConfig` contract:

```python
def build_s2_sensor_config(satellite_id: str = "S2A") -> SensorConfig:
    """
    Construct SensorConfig for Sentinel-2 MSI.

    Uses actual S2A/S2B center wavelengths and FWHMs. Full RSRFs can
    optionally be loaded from the official SentiWiki `S2 Mission` page and its
    linked `Sentinel-2 Spectral Response Functions (S2-SRF)` document.
    """
    # S2A wavelengths shown; S2B has slightly shifted centers
    bands = (
        SensorBand("B01", 443.9, 27.0, 60.0),
        SensorBand("B02", 496.6, 98.0, 10.0),
        SensorBand("B03", 560.0, 45.0, 10.0),
        SensorBand("B04", 664.5, 38.0, 10.0),
        SensorBand("B05", 703.9, 19.0, 20.0),
        SensorBand("B06", 740.2, 18.0, 20.0),
        SensorBand("B07", 782.5, 28.0, 20.0),
        SensorBand("B08", 835.1, 145.0, 10.0),
        SensorBand("B8A", 864.8, 33.0, 20.0),
        SensorBand("B09", 945.0, 26.0, 60.0),
        SensorBand("B10", 1373.5, 75.0, 60.0),
        SensorBand("B11", 1613.7, 143.0, 20.0),
        SensorBand("B12", 2202.4, 242.0, 20.0),
    )
    return SensorConfig(
        sensor_id="MSI",
        satellite_id=satellite_id,
        bands=bands,
        default_ref_scale=1.0 / 10000.0,
        default_ref_offset=0.0,
    )
```

### Band Selection for the Core Pipeline

Using the core plan's wavelength-based selection, the S2 bands map to pipeline
roles as follows:

| Core Pipeline Role | Wavelength Criterion | S2 Band(s) Selected |
|-------------------|---------------------|---------------------|
| `sensor.vis_bands` (surface prior targets) | 400–700 nm | B01, B02, B03, B04 |
| `sensor.nir_bands` (climatology query) | 750–1000 nm | B07, B08, B8A, B09 |
| `sensor.swir_bands` (climatology query) | 1000–2500 nm excl. WV/cirrus | B11, B12 |
| Aerosol-sensitive (solver) | 400–520 nm | B01, B02 |
| Excluded (cirrus/WV) | 1350–1420 nm, 1800–1950 nm | B10 |

### S2 → MODIS Reference Band Mapping

When using `SpectralConvolver.sensor_to_reference()`, the nearest-band matching
for S2 (multispectral, ≤ 30 bands) produces:

| MODIS Reference Band | λ (nm) | Nearest S2 Band | S2 λ (nm) | Δ (nm) |
|---------------------|--------|-----------------|-----------|--------|
| Band 1 (Red) | 645 | B04 | 664.5 | 19.5 |
| Band 2 (NIR) | 858 | B8A | 864.8 | 6.8 |
| Band 3 (Blue) | 469 | B02 | 496.6 | 27.6 |
| Band 4 (Green) | 555 | B03 | 560.0 | 5.0 |
| Band 5 (SWIR-0) | 1240 | — | — | No S2 band near 1240 nm |
| Band 6 (SWIR-1) | 1640 | B11 | 1613.7 | 26.3 |
| Band 7 (SWIR-2) | 2130 | B12 | 2202.4 | 72.4 |

**Note**: MODIS Band 5 (1240 nm) has no S2 equivalent. The `SpectralConvolver`
handles this via interpolation between B8A (865 nm) and B11 (1614 nm), or by
excluding Band 5 from the query dimensions.

---

## 4. S2 Preprocessor (M1) — ObservationBundle Producer

The `Sentinel2Preprocessor` class is the S2-specific implementation of M1. It reads SAFE
directories and produces an `ObservationBundle` (the M1 output contract). Each method
can be overridden independently for fine-grained customisation.

### Method Mapping to Core Contracts

| Core Contract | S2 Method | Output / Notes |
|---------------------|-------------------|----------------|
| `preprocess(path, aoi)` | `Sentinel2Preprocessor.preprocess()` | Composes all methods below → **`ObservationBundle`** |
| detect sensor | `detect_sensor_s2()` (module-level function) | Check for `MTD_MSIL1C.xml` + `S2*_MSIL1C_*` naming |
| get metadata | `.get_metadata()` | Parse `MTD_MSIL1C.xml` → dict with `observation_time`, orbit, baseline, tile ID |
| extract geometry | `.extract_geometry()` | Parse angle grids from `MTD_TL.xml` → `GeometryAngles` (radians) |
| load TOA | `.load_toa()` | Read JP2 bands from `GRANULE/*/IMG_DATA/` → `xr.Dataset` |
| extract cloud mask | `.extract_cloud_mask()` | Read SCL or QA60 → `xr.DataArray` (bool) |
| sensor config | `.build_sensor_config()` | Construct from metadata + known band table → `SensorConfig` |

### How `Sentinel2Preprocessor` Produces `ObservationBundle`

```python
class Sentinel2Preprocessor:
    """S2 preprocessor. Each method can be overridden independently."""

    def get_metadata(self, input_path: Path) -> dict[str, Any]:
        """Parse MTD_MSIL1C.xml for observation time, orbit, baseline, tile ID."""
        ...

    def extract_geometry(self, input_path: Path) -> GeometryAngles:
        """Parse angle grids from MTD_TL.xml, mosaic detectors, convert to radians."""
        ...

    def load_toa(self, input_path: Path) -> xr.Dataset:
        """Read JP2 bands from GRANULE/*/IMG_DATA/ into xr.Dataset."""
        ...

    def extract_cloud_mask(self, input_path: Path) -> xr.DataArray:
        """Read QA60 or SCL and produce boolean cloud mask."""
        ...

    def build_sensor_config(self, metadata: dict) -> SensorConfig:
        """Construct SensorConfig from metadata + known S2 band table."""
        return build_s2_sensor_config(metadata.get("satellite_id", "S2A"))

    def preprocess(self, input_path: Path, aoi: AOI | None = None) -> ObservationBundle:
        """Default implementation: calls the methods above and assembles the bundle."""
        metadata = self.get_metadata(input_path)        # fast
        geometry = self.extract_geometry(input_path)     # fast (5000m grids, interpolated)
        toa = self.load_toa(input_path)                 # slow (reads JP2 bands)
        cloud_mask = self.extract_cloud_mask(input_path)
        sensor_config = self.build_sensor_config(metadata)

        # Derive CRS and bounds from TOA
        crs = str(toa.attrs.get("crs", toa.rio.crs))
        bounds = tuple(toa.rio.bounds())

        # Optionally clip to AOI
        if aoi is not None:
            toa = aoi.clip(toa)
            bounds = aoi.get_bounds(crs)

        return ObservationBundle(
            toa=toa,
            geometry=geometry,
            cloud_mask=cloud_mask,
            sensor_config=sensor_config,
            metadata=metadata,
            crs=crs,
            bounds=bounds,
        )
```

### S2-Specific Concerns Inside M1 (Not Exposed to Core)

These are implementation details internal to the S2 preprocessor class. The core pipeline
and downstream modules never see them — they only see the `ObservationBundle`.

1. **Multi-resolution bands**: S2 has 10 m, 20 m, and 60 m bands. The preprocessor
   stores them in the `xr.Dataset` with per-variable resolution metadata. M4 handles
   resampling to the aux grid.
2. **Per-detector angle grids**: S2 provides angles per detector (12 detectors per band).
   These overlap at swath edges and must be mosaicked. This is handled inside
   `extract_geometry()`, which returns a single `GeometryAngles` object.
3. **Processing baseline**: Baseline N0400+ changed the reflectance offset from 0 to
   -1000. The preprocessor accounts for this via `default_ref_offset` in `SensorConfig`.
4. **SAFE directory structure**: Nested structure with `GRANULE/`, `DATASTRIP/`, etc.
   All path traversal is encapsulated in the preprocessor.

### User Override Points for S2 M1

Users can customise parts of S2 preprocessing by **subclassing** and overriding
individual methods. The IDE shows all available methods, and stack traces remain
clear (e.g. `MyS2Preprocessor.extract_cloud_mask`).

```python
# Override just the cloud mask (e.g. use an external cloud model)
class MyS2Preprocessor(Sentinel2Preprocessor):
    def extract_cloud_mask(self, input_path: Path) -> xr.DataArray:
        return my_model.predict(input_path)

result = siac_process(config, input_path, preprocessor=MyS2Preprocessor())
```

Or replace M1 entirely with a standalone function:

```python
def my_s2_loader(input_path: Path, aoi=None) -> ObservationBundle:
    # Custom S2 loading (e.g. from COG, STAC, or xarray-sentinel)
    ...

result = siac_process(config, input_path, preprocessor=my_s2_loader)
```

---

## 5. S2 Multi-Resolution Handling (M4 / M6 Concerns)

S2 is the primary example of a **multi-resolution** sensor in SIAC. The core plan's
multi-resolution data assembly applies as follows for S2.

### Resolution Mapping

```
S2 Native Resolutions (inside ObservationBundle from M1):

  10 m:   B02, B03, B04, B08          (4 bands)
  20 m:   B05, B06, B07, B8A, B11, B12  (6 bands)
  60 m:   B01, B09, B10               (3 bands)

Pipeline Resolution Hierarchy (handled by M4 and M6):

  TOA resolution     → per-band native (10/20/60 m)     — M1 output, M6 input
  Aux resolution     → configurable, default 500 m      — M4 output, M5 input
  Aerosol resolution → configurable, default 1000 m     — M5 output
```

### Band Resolution Handling in M4 (Grid Assembler)

For S2, all bands are resampled to the **aux resolution** (e.g. 500 m) by M4
before entering the `SolverInputBundle`. The resampling method depends on band type:

| Band Type | Native (m) | Resampling to Aux (500 m) | Method |
|-----------|-----------|--------------------------|--------|
| Reflectance (B02–B12) | 10/20/60 | Average | Bilinear or area-weighted mean |
| Cloud mask (QA60) | 60 | Nearest | Conservative (if any pixel is cloud → cloud) |
| Angles (SZA, VZA, SAA, VAA) | 5000 | Bilinear interpolation | Already coarse; interpolate to aux grid |

### Final Correction Resolution (M6)

The atmospheric correction (TOA → BOA) in M6 is applied at **each band's native
resolution**. M6 receives the original `ObservationBundle` from M1 (at native
resolution) and the `SolvedAtmosphere` from M5. Solved AOT/TCWV (at aerosol
resolution, e.g. 1000 m) are upsampled to each band's native resolution using
bilinear interpolation:

```
Solved AOT (1000 m)  →  upsampled to 10 m  →  correct B02, B03, B04, B08
                     →  upsampled to 20 m  →  correct B05–B07, B8A, B11, B12
                     →  upsampled to 60 m  →  correct B01, B09
```

This preserves the native spatial detail of each band in the corrected output
(`CorrectionResult.boa`).

---

## 6. S2-Specific Custom Overrides

This section shows how users can inject custom functions at each module boundary
for S2 processing, while still producing the correct output contracts.

### 6.1 Custom Cloud Mask (Override Part of M1)

```python
class Sen2CorCloudPreprocessor(Sentinel2Preprocessor):
    """S2 preprocessor with Sen2Cor SCL-based cloud masking."""
    def extract_cloud_mask(self, input_path: Path) -> xr.DataArray:
        scl = load_sen2cor_scl(input_path)
        return scl.isin([3, 8, 9, 10])  # cloud shadow, cirrus, cloud high/med

result = siac_process(config, input_path, preprocessor=Sen2CorCloudPreprocessor())
```

### 6.2 Custom Atmospheric Prior (Override M2)

```python
# Use MERRA-2 reanalysis instead of CAMS
from siac.adapters.atmo.merra2 import MERRA2Provider
result = siac_process(config, input_path, atmo_provider=MERRA2Provider().get_prior)
```

### 6.3 Custom Surface Prior from Landsat Composites (Override M3)

```python
def landsat_surface_prior(bounds, crs, obs_time, sensor_config, geometry, resolution) -> SurfacePrior:
    """Build surface prior from Landsat 30-day composite."""
    composite = load_landsat_composite(bounds, crs, obs_time)
    # Convolve Landsat bands to S2 sensor bands via spectral functions
    boa = reference_to_sensor(composite, sensor_config.vis_bands, sensor_config)
    unc = estimate_uncertainty(composite)
    mask = np.isfinite(boa)
    return SurfacePrior(boa=boa, boa_unc=unc, mask=mask)

result = siac_process(config, input_path, surface_prior_provider=landsat_surface_prior)
```

### 6.4 Pre-Solved AOT (Override M5)

```python
def my_aot_retrieval(inputs: SolverInputBundle, config: SolverConfig) -> SolvedAtmosphere:
    """Use AOT from AERONET or external retrieval."""
    aot_field = interpolate_aeronet(station_data, inputs.geometry)
    return SolvedAtmosphere(
        atmo_state=inputs.atmo_prior.with_updated_aot_tcwv(aot_field, inputs.atmo_prior.tcwv),
        aot=aot_field, tcwv=inputs.atmo_prior.tcwv,
        aot_unc=inputs.atmo_prior.aot_unc, tcwv_unc=inputs.atmo_prior.tcwv_unc,
        cost_final=0.0, n_iterations=0, converged=True,
    )

result = siac_process(config, input_path, solver=my_aot_retrieval)
```

### 6.5 Custom Output Format (Override After M6)

```python
result = siac_process(config, Path("/path/to/S2_SAFE/"))
# CorrectionResult has standard fields — user can write to any format
write_to_stac(result.boa, result.metadata, output_dir)
```

---

## 7. Files Summary

### S2-Specific Files (New)

| File | Action | Description | Produces |
|------|--------|-------------|----------|
| `python/siac/io/s2_data_source.py` | **NEW** | `S2Query`, `S2Product`, `fetch_s2()`, `search_s2()` | Local `Path` (input to M1) |
| `python/siac/io/copernicus_dataspace.py` | **NEW** | `search_cdse()`, `download_cdse()` functions | `list[S2Product]` / local `Path` |
| `python/siac/io/gcs_sentinel2.py` | **NEW** | `search_gcs()`, `download_gcs()` functions | `list[S2Product]` / local `Path` |

### Modified Core Files

| File | Modification | Description |
|------|-------------|-------------|
| `python/siac/core/config.py` | Add `S2DataAccessConfig` | S2 data access configuration sub-model |
| `python/siac/siac.py` | Add `siac_process_s2()`, `resolve_s2_input()` | S2 convenience entry points (pre-M1 resolution) |
| `python/siac/satellite/sentinel2.py` | Update to export `Sentinel2Preprocessor` class | Return `ObservationBundle` via `.preprocess()` |

---

## 8. Usage Examples

```python
from siac import siac_process, siac_process_s2
from siac.adapters.data.s2_data_source import S2Query, search_sentinel2
from siac.config import SIACConfig
from datetime import date
from pathlib import Path

config = SIACConfig(sensor="s2")

# 1. Local SAFE (current behavior, unchanged)
#    M1 reads SAFE → ObservationBundle → M2/M3/M4/M5/M6 → CorrectionResult
result = siac_process(config, Path("/data/S2B_MSIL1C_...SAFE"))

# 2. Product ID → auto-download → M1 → ... → CorrectionResult
result = siac_process_s2(config, "S2B_MSIL1C_20210801T103629_N0500_R008_T31UDQ_20230731T120241")

# 3. Tile + date shorthand → search → newest baseline → download → M1 → ...
result = siac_process_s2(config, "T31UDQ_20210801")

# 4. Complex query
query = S2Query(
    bbox=(1.5, 50.0, 2.5, 51.0),
    start_date=date(2021, 8, 1),
    end_date=date(2021, 8, 15),
    max_cloud_cover=30.0,
)
result = siac_process_s2(config, query)

# 5. Custom modules for S2
result = siac_process(
    config,
    Path("/data/S2_SAFE/"),
    preprocessor=MyS2Preprocessor(),    # custom M1 subclass
    atmo_provider=my_custom_atmo,       # custom M2 (function or class.get_prior)
    surface_prior_provider=my_prior,    # custom M3 (function or class.get_surface_prior)
)

# 6. Search without processing
products = search_sentinel2(
    bbox=(1.5, 50.0, 2.5, 51.0),
    start_date="2021-08-01",
    end_date="2021-08-15",
)
for p in products:
    print(f"{p.product_id}  baseline={p.processing_baseline}  cloud={p.cloud_cover:.0f}%")
```

---

## 9. Verification Plan

### Unit Tests

1. **S2Query parsing**: product ID, tile+date, validation
2. **S2Product baseline comparison**: deduplication logic
3. **Band mapping**: `build_s2_sensor_config()` produces correct wavelengths
4. **SpectralConvolver**: `sensor_to_reference()` maps S2 → MODIS reference bands correctly
5. **ObservationBundle construction**: `Sentinel2Preprocessor().preprocess()` returns a valid
   `ObservationBundle` that passes contract validation

### Integration Tests

1. **Mock OData/S3 responses**: `CopernicusDataspaceBackend.search()` returns correctly parsed `S2Product` list
2. **Baseline dedup**: given products with baselines N0301, N0400, N0500 for same tile+date → only N0500 survives
3. **Multi-resolution**: verify TOA loading at native resolutions, M4 resampling to aux grid, M6 correction at native resolution
4. **Custom provider injection**: verify that a user-supplied M2/M3/M5 function is called
   and its output flows through the pipeline correctly

### Contract Validation Tests

1. **ObservationBundle validation**: missing `observation_time` → clear error
2. **AtmosphericState validation**: negative uncertainties → clear error
3. **SurfacePrior validation**: mismatched shapes → clear error
4. **SolverInputBundle validation**: misaligned grids → clear error

### Regression Tests

All 107 existing tests must still pass (no modification to core solver/corrector contracts).
