# SIAC v2 — Code Structure & Module Dependency Graph

## Package Tree

```
siac/
├── __init__.py
├── siac.py                          # Main entry point & _resolve_* helpers
├── pipeline.py                      # Pipeline orchestrator (run_pipeline)
│
├── core/                            # Shared types, protocols, config
│   ├── __init__.py
│   ├── types.py                     # All data types & sensor configs
│   ├── protocols.py                 # @runtime_checkable Protocol definitions
│   ├── config.py                    # SIACConfig (pydantic_settings)
│   ├── exceptions.py                # Error hierarchy
│   ├── spectral.py                  # SpectralBandDescriptor, band convolution
│   ├── validation.py                # Contract validators (_validate_*)
│   └── aoi.py                       # AOI (area of interest) handling
│
├── satellite/                       # M1 — Satellite Preprocessors
│   ├── __init__.py
│   ├── base.py                      # BaseSatellitePreprocessor, registry
│   └── sentinel2.py                 # Sentinel-2 SAFE reader
│
├── priors/                          # M2 & M3 — Prior Providers
│   ├── __init__.py
│   ├── atmospheric/                 # M2 — Atmospheric priors
│   │   ├── __init__.py
│   │   ├── cams.py                  # ECMWF CAMS provider
│   │   └── merra2.py               # NASA MERRA-2 provider
│   ├── brdf/                        # BRDF kernel data & computation
│   │   ├── __init__.py
│   │   ├── kernels.py               # Ross-Thick Li-Sparse kernels
│   │   ├── mcd43_earthaccess.py     # MODIS MCD43 via earthaccess
│   │   └── gee_stub.py             # Google Earth Engine (stub)
│   └── surface/                     # M3 — Surface prior derivation
│       ├── __init__.py
│       ├── kernel_model.py          # BRDF → SurfacePrior deriver
│       └── prior_store.py           # Pre-built Zarr prior store
│
├── rt/                              # RT Model Backends
│   ├── __init__.py
│   ├── emulator/
│   │   ├── __init__.py
│   │   └── two_nn.py               # Two-layer NN emulator
│   ├── lut/
│   │   ├── __init__.py
│   │   ├── constants.py            # Default LUT URLs and coordinate constants
│   │   ├── http_zip_store.py       # ReadOnlyZipFileSystem-style ZIP access (local/HTTP/S3)
│   │   ├── store.py                # Local/remote/ZIP store resolution
│   │   ├── backend.py              # ZarrLUTBackend interpolation engine
│   │   ├── create.py               # LUT generation from Py6S
│   │   └── zarr_lut.py             # Compatibility facade (re-exports)
│   └── direct/
│       └── __init__.py              # (placeholder)
│
├── grid/                            # M4 — Grid Assembly
│   ├── __init__.py
│   └── assembler.py                 # assemble_grids() resampler
│
├── solver/                          # M5 — Aerosol Solver
│   ├── __init__.py
│   ├── cost.py                      # Cost function & gradient
│   └── multigrid.py                 # Multi-grid optimiser
│
├── correction/                      # M6 — Atmospheric Correction
│   ├── __init__.py
│   └── atmospheric.py               # TOA → BOA corrector
│
└── io/                              # I/O & Geospatial Utilities
    ├── __init__.py
    ├── readers.py                   # Raster/HDF/NetCDF/Zarr readers
    ├── writers.py                   # COG/Zarr/NetCDF writers
    ├── reprojection.py              # CRS transforms, resampling, clipping
    ├── geometry.py                  # AOI loading, bounds math
    ├── earthaccess_source.py        # NASA Earthdata wrapper
    ├── s2_data_source.py            # S2Query, S2Product, deduplication
    ├── copernicus_dataspace.py      # CDSE backend (stub)
    └── gcs_sentinel2.py             # GCS public bucket backend (stub)
```

---

## Pipeline Data Flow (M1 → M6)

The six functional modules form a DAG. Arrows show which contract type
flows between them.

```
                        ┌─────────────────────────────┐
                        │       siac.siac              │
                        │  siac_process() / SIAC       │
                        └──────────┬──────────────────┘
                                   │ delegates to
                                   ▼
                        ┌─────────────────────────────┐
                        │      siac.pipeline           │
                        │      run_pipeline()          │
                        └──────────┬──────────────────┘
                                   │
              ┌────────────────────┼────────────────────┐
              │                    │                     │
              ▼                    ▼                     ▼
   ┌──────────────────┐ ┌──────────────────┐ ┌──────────────────┐
   │ M1: Preprocessor  │ │ M2: Atmo Prior   │ │ M3: Surface Prior│
   │                   │ │                  │ │                  │
   │ siac.satellite    │ │ siac.priors      │ │ siac.priors      │
   │   .sentinel2      │ │   .atmospheric   │ │   .surface       │
   │                   │ │   .cams          │ │   .kernel_model  │
   │  → Observation    │ │   .merra2        │ │   .prior_store   │
   │    Bundle         │ │                  │ │                  │
   │                   │ │  → Atmospheric   │ │  → SurfacePrior  │
   └────────┬──────────┘ │    State         │ │                  │
            │            └────────┬─────────┘ └────────┬─────────┘
            │                     │   (concurrent)      │
            │                     │                     │
            └─────────┬───────────┴─────────────────────┘
                      │
                      ▼
           ┌──────────────────────┐
           │ M4: Grid Assembler   │
           │                      │
           │ siac.grid.assembler  │
           │ assemble_grids()     │
           │                      │       ┌──────────────────┐
           │  → SolverInputBundle │◄──────│  RT Model Backend│
           └──────────┬───────────┘       │  siac.rt         │
                      │                   │    .emulator      │
                      ▼                   │    .lut           │
           ┌──────────────────────┐       └──────────────────┘
           │ M5: Aerosol Solver   │              │
           │                      │              │
           │ siac.solver          │              │
           │   .multigrid         │              │
           │   .cost              │              │
           │                      │              │
           │  → SolvedAtmosphere  │              │
           └──────────┬───────────┘              │
                      │                          │
                      ▼                          │
           ┌──────────────────────┐              │
           │ M6: Corrector        │◄─────────────┘
           │                      │
           │ siac.correction      │
           │   .atmospheric       │
           │                      │
           │  → CorrectionResult  │
           └──────────────────────┘
```

---

## Module Dependency Graph

Each arrow reads **"imports from"**. Leaf nodes have no internal siac
dependencies.

```
siac.siac
 ├── siac.core.config
 ├── siac.core.types
 ├── siac.core.aoi
 ├── siac.pipeline
 ├── siac.satellite
 ├── siac.priors.atmospheric.cams
 ├── siac.priors.surface.kernel_model
 ├── siac.rt.emulator
 ├── siac.rt.lut
 ├── siac.solver
 ├── siac.correction
 └── siac.grid.assembler

siac.pipeline
 ├── siac.core.aoi
 ├── siac.core.types
 └── siac.core.validation
      └── siac.core.types

siac.satellite.sentinel2
 ├── siac.core.types
 ├── siac.satellite.base
 │    ├── siac.core.types
 │    └── siac.io.reprojection          ← leaf
 └── siac.io                            ← leaf (readers/writers)

siac.priors.atmospheric.cams
 └── siac.core.types

siac.priors.atmospheric.merra2
 ├── siac.core.types
 └── siac.io.earthaccess_source         ← leaf

siac.priors.brdf.mcd43_earthaccess
 ├── siac.core.types
 └── siac.io.earthaccess_source         ← leaf

siac.priors.brdf.kernels                ← leaf (+ optional Rust)

siac.priors.surface.kernel_model
 ├── siac.core.types
 └── siac.priors.brdf.kernels           ← leaf

siac.priors.surface.prior_store
 └── siac.core.types

siac.grid.assembler
 └── siac.core.types

siac.solver.multigrid
 ├── siac.core.types
 ├── siac.core.protocols
 │    └── siac.core.types               (TYPE_CHECKING)
 └── siac.solver.cost
      ├── siac.core.types
      └── siac.core.protocols

siac.correction.atmospheric
 ├── siac.core.types
 └── siac.core.protocols

siac.rt.emulator.two_nn
 └── siac.core.types                    (+ optional Rust)

siac.rt.lut.backend
 ├── siac.core.types
 └── siac.rt.lut.store
      └── siac.rt.lut.http_zip_store

siac.rt.lut.create                   ← leaf (Py6S optional)
siac.rt.lut.zarr_lut                 ← compatibility re-export facade

siac.core.aoi
 ├── siac.io.geometry                   ← leaf
 └── siac.io.reprojection               ← leaf

siac.core.spectral                      ← leaf
siac.core.config                        ← leaf (pydantic_settings)
siac.core.exceptions                    ← leaf
siac.core.types                         ← leaf (foundation)
siac.io.readers                         ← leaf
siac.io.writers                         ← leaf
siac.io.reprojection                    ← leaf
siac.io.geometry                        ← leaf
siac.io.earthaccess_source              ← leaf
siac.io.s2_data_source                  ← leaf
siac.io.copernicus_dataspace
 └── siac.io.s2_data_source             ← leaf
siac.io.gcs_sentinel2
 └── siac.io.s2_data_source             ← leaf
```

---

## Contract Types & Where They Live

All contract types are defined in **`siac.core.types`** and validated by
**`siac.core.validation`**.

| Contract | Produced by (module) | Consumed by | Import path |
|---|---|---|---|
| `ObservationBundle` | M1 `siac.satellite.*` | M4, M6 | `siac.core.types.ObservationBundle` |
| `AtmosphericState` | M2 `siac.priors.atmospheric.*` | M4 | `siac.core.types.AtmosphericState` |
| `SurfacePrior` | M3 `siac.priors.surface.*` | M4 | `siac.core.types.SurfacePrior` |
| `SolverInputBundle` | M4 `siac.grid.assembler` | M5 | `siac.core.types.SolverInputBundle` |
| `SolvedAtmosphere` | M5 `siac.solver.multigrid` | M6 | `siac.core.types.SolvedAtmosphere` |
| `CorrectionResult` | M6 `siac.correction.atmospheric` | User | `siac.core.types.CorrectionResult` |

Supporting types used within contracts:

| Type | Import path | Used by |
|---|---|---|
| `GeometryAngles` | `siac.core.types.GeometryAngles` | ObservationBundle, SolverInputBundle |
| `SensorConfig` | `siac.core.types.SensorConfig` | ObservationBundle, SolverInputBundle |
| `SensorBand` | `siac.core.types.SensorBand` | SensorConfig, SolverInputBundle |
| `BRDFKernelWeights` | `siac.core.types.BRDFKernelWeights` | SurfacePrior |
| `RTCoefficients` | `siac.core.types.RTCoefficients` | RT backends |
| `SpectralBandDescriptor` | `siac.core.spectral.SpectralBandDescriptor` | Spectral convolution |

---

## Protocol Interfaces (`siac.core.protocols`)

Protocols define the plug-in points. Any class implementing the required
methods is accepted — no base class inheritance needed.

| Protocol | Required methods | Implemented by |
|---|---|---|
| `RTModelBackend` | `compute_coefficients()`, `supports_jacobian()`, `backend_name`, `is_available_for_sensor()` | `siac.rt.emulator.two_nn.TwoLayerNNEmulator`, `siac.rt.lut.backend.ZarrLUTBackend` |
| `SatellitePreprocessor` | `load_toa()`, `extract_geometry()`, `extract_cloud_mask()`, `get_metadata()` | `siac.satellite.sentinel2.Sentinel2Preprocessor` |
| `AtmosphericPriorProvider` | `get_prior()`, `source_name` | `siac.priors.atmospheric.cams.CAMSProvider`, `siac.priors.atmospheric.merra2.MERRA2Provider` |
| `BRDFProductProvider` | `get_brdf_parameters()`, `source_name` | `siac.priors.brdf.mcd43_earthaccess.MCD43EarthAccessProvider` |
| `SurfacePriorDeriver` | `compute_surface_prior()` | `siac.priors.surface.kernel_model.KernelModelDeriver` |
| `AerosolSolver` | `solve()` | `siac.solver.multigrid.MultiGridSolver` |
| `OutputWriter` | `write_boa()`, `write_auxiliary()` | (user-provided) |

---

## Callable Type Aliases (`siac.pipeline`)

The pipeline accepts plain functions or bound methods for each module slot.

| Alias | Signature | Default resolver |
|---|---|---|
| `PreprocessorFn` | `(Path, AOI\|None) -> ObservationBundle` | `siac.siac._resolve_preprocessor` |
| `AtmoPriorFn` | `(bounds, crs, datetime, float) -> AtmosphericState` | `siac.siac._resolve_atmo_provider` |
| `SurfacePriorFn` | `(bounds, crs, datetime, SensorConfig, GeometryAngles, float) -> SurfacePrior` | `siac.siac._resolve_surface_prior_provider` |
| `GridAssemblerFn` | `(ObservationBundle, AtmosphericState, SurfacePrior, RT, float, float) -> SolverInputBundle` | `siac.siac._resolve_grid_assembler` |
| `SolverFn` | `(SolverInputBundle, config) -> SolvedAtmosphere` | `siac.siac._resolve_solver` |
| `CorrectorFn` | `(ObservationBundle, SolvedAtmosphere, RT) -> CorrectionResult` | `siac.siac._resolve_corrector` |

---

## I/O & Data Access Layer

```
siac.io
 ├── readers.py          read_raster(), read_jp2(), read_zarr_array(), ...
 ├── writers.py          write_cog(), write_zarr(), write_dataset(), ...
 ├── reprojection.py     reproject_to_crs(), resample(), clip_to_bounds(), ...
 ├── geometry.py         load_aoi(), bounds_intersect(), create_grid_from_bounds(), ...
 ├── earthaccess_source  EarthAccessSource (NASA Earthdata wrapper)
 │
 ├── s2_data_source.py   S2Query, S2Product, S2DataBackend protocol
 │    ├── copernicus_dataspace.py   CopernicusDataspaceBackend (stub)
 │    └── gcs_sentinel2.py         GCSSentinel2Backend (stub)
 │
 └── (used by)
      ├── siac.satellite.sentinel2   (reads SAFE via siac.io.readers)
      ├── siac.priors.atmospheric    (CAMS NetCDF, MERRA-2 via earthaccess)
      ├── siac.priors.brdf           (MCD43 HDF via earthaccess)
      └── siac.core.aoi              (geometry + reprojection)
```

---

## Test Structure

```
tests/
├── conftest.py                      # Shared fixtures (MockRTModel, bundles, callables)
├── unit/
│   ├── test_types.py                # GeometryAngles, AtmosphericState, SensorConfig
│   ├── test_contracts.py            # ObservationBundle, SolverInputBundle, etc.
│   ├── test_validation.py           # _validate_* functions
│   ├── test_spectral.py             # SpectralBandDescriptor, convolution
│   ├── test_grid_assembler.py       # assemble_grids()
│   ├── test_prior_store.py          # PrebuiltPriorStore (Zarr fixtures)
│   ├── test_s2_data_source.py       # S2Query, S2Product, deduplication
│   ├── test_resolve.py              # _resolve_* helpers (requires pydantic_settings)
│   ├── test_correction.py           # AtmosphericCorrector
│   ├── test_solver.py               # Cost function, multi-grid solver
│   ├── test_satellite.py            # Preprocessor registry, sensor detection
│   ├── test_io.py                   # Readers, writers
│   └── test_config.py               # SIACConfig (requires pydantic_settings)
├── integration/
│   ├── test_pipeline.py             # run_pipeline() orchestration & concurrency
│   ├── test_injection.py            # Custom callable injection & contract violation
│   └── test_e2e_synthetic.py        # Full M1→M6 with synthetic data
└── benchmarks/
    └── test_perf.py                 # Grid assembler performance (256x256, 512x512)
```
