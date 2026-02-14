# SIAC v2 — Core Architecture Plan

This document defines the **sensor-agnostic, modular** architecture for SIAC v2.

- It specifies the core contracts (data model + interfaces) that any sensor integration must follow.
- Every module has a **defined output contract** — a typed dataclass that satisfies the exact input requirements of downstream consumers (solver, corrector).
- Users can replace any module with a **custom callable** as long as it returns the correct output type.
- Sensor-specific designs (e.g. Sentinel-2) live in dedicated plans (see [PLANS_S2.md](PLANS_S2.md)).

The intent is to keep this file stable as the "core", so adding a new sensor or swapping a data provider mostly means:
1) implement a callable/class that returns the expected output contract, and
2) register it with the pipeline,
without changing the solver or correction stages.

---

## Table of Contents

1. [Design Goals](#1-design-goals)
2. [Module Boundary Model](#2-module-boundary-model)
3. [Module Output Contracts (The Currency Types)](#3-module-output-contracts-the-currency-types)
4. [Module Definitions](#4-module-definitions)
5. [Custom Provider Pattern (User-Injected Functions)](#5-custom-provider-pattern-user-injected-functions)
6. [AOI Scoping Rules (Global)](#6-aoi-scoping-rules-global)
7. [Pipeline Orchestration (Sequential → Parallel)](#7-pipeline-orchestration-sequential--parallel)
8. [Core Requirements for Adding New Sensors](#8-core-requirements-for-adding-new-sensors)
9. [Sensor-Agnostic Spectral Model & Surface Prior](#9-sensor-agnostic-spectral-model--surface-prior)
10. [Pluggable Data Providers (Atmosphere / BRDF / Surface Prior)](#10-pluggable-data-providers-atmosphere--brdf--surface-prior)
11. [Centralised Authentication](#11-centralised-authentication)
12. [Earthaccess Rollout Plan (Pre-Implementation)](#12-earthaccess-rollout-plan-pre-implementation)

---

## 1. Design Goals

### Goals

- **Right tool for the job**: use plain functions for stateless transforms (grid assembly, solver, corrector, spectral math). Use small classes for modules that hold configured state (data directory paths, auth credentials, cached connections). Never build deep class hierarchies — keep classes small and flat.
- **Explicit module contracts**: every module declares a typed output that exactly matches the input requirements of downstream consumers. No implicit coupling.
- **User-replaceable modules**: users can inject their own callable (function or class instance) for any module as long as the return type satisfies the contract. Internal data sourcing is an implementation detail; the pipeline only cares about the output type.
- **Sensor-agnostic correction pipeline**: the solver and correction stages must not depend on sensor-specific band names or file formats.
- **Pluggable data providers**: atmospheric priors (CAMS/MERRA-2/etc.) and surface priors (BRDF / prior stores) are interchangeable.
- **AOI-scoped I/O by default**: all auxiliary data reads and remote fetches are bounded by the AOI extent.
- **Composable performance**: independent I/O tasks (satellite bands, priors, DEM, surface prior) should be parallelisable.
- **Extensibility**: adding a new sensor should require adding a new preprocessor + config, not rewriting the core.

### Non-goals

- Implementing every sensor's data access in core (those are sensor plans).
- Hard-coding any sensor-specific band names or SAFE/MTL layout in the solver.
- Prescribing how a module acquires its data internally — only its output shape matters.
- Building deep class hierarchies with abstract base classes that require understanding multiple layers to debug.

---

## 2. Module Boundary Model

### 2.1 Architecture Overview

SIAC v2 is composed of **six functional modules**, each with a strictly defined output contract. The pipeline orchestrator assembles these outputs and feeds them into the solver and corrector.

```
┌─────────────────────────────────────────────────────────────────────┐
│                        SIAC Pipeline                                │
│                                                                     │
│  ┌──────────────┐   ┌──────────────┐   ┌──────────────┐            │
│  │ M1: Satellite │   │ M2: Atmo     │   │ M3: Surface  │            │
│  │ Preprocessor  │   │ Prior        │   │ Prior        │            │
│  │              │   │ Provider     │   │ Provider     │            │
│  │ → ObsBundle  │   │ → AtmoState  │   │ → SurfPrior  │            │
│  └──────┬───────┘   └──────┬───────┘   └──────┬───────┘            │
│         │                  │                  │                     │
│  ┌──────┴──────────────────┴──────────────────┴───────┐            │
│  │              M4: Grid Assembler                     │            │
│  │              → SolverInputBundle                    │            │
│  └──────────────────────┬──────────────────────────────┘            │
│                         │                                           │
│  ┌──────────────────────▼──────────────────────────────┐            │
│  │              M5: Aerosol Solver                      │            │
│  │              → SolvedAtmosphere                      │            │
│  └──────────────────────┬──────────────────────────────┘            │
│                         │                                           │
│  ┌──────────────────────▼──────────────────────────────┐            │
│  │              M6: Atmospheric Corrector               │            │
│  │              → CorrectionResult                      │            │
│  └─────────────────────────────────────────────────────┘            │
└─────────────────────────────────────────────────────────────────────┘
```

### 2.2 Module Dependency Graph

Each module's output feeds exactly into the inputs of downstream modules. No cross-cutting data flows.

```
  M1 (ObservationBundle) ───────────────┐
                                        │
  M2 (AtmosphericState) ───────────────┤
                                        ├──▶ M4 (SolverInputBundle) ──▶ M5 (SolvedAtmosphere) ──▶ M6 (CorrectionResult)
  M3 (SurfacePrior) ───────────────────┤
                                        │
  RT Backend (local or remote store) ──┘
```

### 2.3 LUT Backend Internal Split

To keep LUT code maintainable, the Zarr LUT backend is separated by concern:

- `siac.rt.lut.backend`: RT interpolation logic (`ZarrLUTBackend`)
- `siac.rt.lut.store`: local/remote/S3/ZIP store resolution
- `siac.rt.lut.http_zip_store`: ReadOnlyZipFileSystem-style ZIP access for local/HTTP/S3 LUT archives
- `siac.rt.lut.create`: LUT generation utilities (`create_lut_from_py6s`)
- `siac.rt.lut.constants`: default public LUT URL and LUT coordinate constants
- `siac.rt.lut.zarr_lut`: compatibility facade that re-exports the public API

### 2.4 Key Design Rule

> **A module is any callable** that accepts its documented inputs and returns its documented output contract type. The pipeline does not care whether the callable is a function, a bound method, or a class instance with `__call__` — only that the return value satisfies the frozen dataclass contract.

This rule enables:
- Replacing any built-in module with a user-provided function *or* class
- Testing modules in isolation with mock inputs
- Composing modules in novel ways (e.g. skip the solver and supply pre-solved AOT)

**When to use a function vs. a class:**

| Use a **function** when… | Use a **small class** when… |
|---|---|
| The transform is stateless (inputs → output, no config) | The module needs configured state (directory paths, auth, cached connections) |
| There's nothing to initialise | There's meaningful setup (load LUTs, open sessions, validate paths) |
| Examples: `assemble_grids()`, `solve_aerosol()`, `correct_atmosphere()`, `sensor_to_reference()` | Examples: `CAMSProvider(cams_dir)`, `Sentinel2Preprocessor()`, `EmulatorBackend(weights_path)` |

Either way, keep classes **small and flat** — no deep inheritance, no abstract method chains. A class with one public method that takes no constructor args should be a function instead.

---

## 3. Module Output Contracts (The Currency Types)

These frozen dataclasses are the **only types** that cross module boundaries. They live in `python/siac/core/types.py`.

### 3.1 `ObservationBundle` — Output of M1 (Satellite Preprocessor)

Everything extracted from the satellite input, in a single package:

```python
@dataclass(frozen=True)
class ObservationBundle:
    """Complete observation data from satellite preprocessing.

    This is the single output contract of M1 (Satellite Preprocessor).
    It contains everything the pipeline needs from the satellite input.
    """
    toa: xr.Dataset                  # TOA reflectance, one var per band
    geometry: GeometryAngles         # SZA/SAA/VZA/VAA in radians
    cloud_mask: xr.DataArray         # bool, True = cloudy/invalid
    sensor_config: SensorConfig      # band definitions + scale/offset
    metadata: dict[str, Any]         # must include 'observation_time': datetime
    crs: str                         # e.g. "EPSG:32632"
    bounds: tuple[float, float, float, float]  # (xmin, ymin, xmax, ymax) in crs
```

**Validation rules** (checked at construction):
- `metadata["observation_time"]` must exist and be a `datetime`
- `toa` must have spatial dimensions and a CRS
- `geometry` arrays must be broadcastable to `toa` spatial dimensions
- `cloud_mask` must match `toa` spatial shape

### 3.2 `AtmosphericState` — Output of M2 (Atmospheric Prior Provider)

Already defined in the existing codebase. Repeated here for completeness:

```python
@dataclass(frozen=True)
class AtmosphericState:
    """Atmospheric parameters for radiative transfer.

    Output contract of M2. Also reused as the output of M5 (solver)
    via `with_updated_aot_tcwv()`.
    """
    aot: xr.DataArray           # AOT at 550 nm
    tcwv: xr.DataArray          # Total Column Water Vapour (g/cm²)
    tco3: xr.DataArray          # Total Column Ozone (atm-cm)
    aot_unc: xr.DataArray       # AOT uncertainty
    tcwv_unc: xr.DataArray      # TCWV uncertainty
    tco3_unc: xr.DataArray      # TCO3 uncertainty
    elevation: xr.DataArray     # Surface elevation (km)
```

**Validation rules**:
- All arrays must be spatially aligned (same grid or broadcastable)
- Uncertainty arrays must be non-negative
- `elevation` is in km above sea level

### 3.3 `SurfacePrior` — Output of M3 (Surface Prior Provider)

```python
@dataclass(frozen=True)
class SurfacePrior:
    """Surface reflectance prior for the solver.

    Output contract of M3. Provides expected BOA and uncertainty
    that the solver uses as a constraint.
    """
    boa: xr.DataArray           # Prior BOA reflectance (bands, y, x) or (y, x)
    boa_unc: xr.DataArray       # Uncertainty in BOA
    mask: xr.DataArray          # bool, True = valid pixel
```

**Validation rules**:
- `boa` and `boa_unc` must have matching shapes
- `mask` must broadcast to `boa` spatial dimensions
- Band dimension (if present) must match the expected solver band count

### 3.4 `SolverInputBundle` — Output of M4 (Grid Assembler)

The grid assembler resamples and aligns all upstream outputs to the solver's
working resolution. This is the **complete, validated** input to the solver —
the solver never fetches or aligns data itself.

```python
@dataclass(frozen=True)
class SolverInputBundle:
    """All inputs to the aerosol solver, resampled to solver grids.

    Output contract of M4. Everything in this bundle is spatially
    aligned and ready for the solver to consume directly.
    """
    # Observation (resampled to aux resolution)
    toa: xr.DataArray               # (bands, y, x) at aux resolution
    geometry: GeometryAngles        # at aux resolution
    cloud_mask: xr.DataArray        # (y, x) at aux resolution
    sensor_config: SensorConfig
    bands: list[SensorBand]         # solver bands (wavelength-selected)

    # Atmospheric prior (resampled to aux resolution)
    atmo_prior: AtmosphericState    # all fields at aux resolution

    # Surface prior (resampled to aux resolution)
    surface_prior: SurfacePrior     # at aux resolution

    # RT backend (not resampled — it's a model, not raster data)
    rt_model: RTModelBackend

    # Grid metadata
    aux_resolution_m: float         # e.g. 500.0
    aerosol_resolution_m: float     # e.g. 1000.0
```

**Validation rules**:
- All raster fields (`toa`, `geometry.*`, `cloud_mask`, `atmo_prior.*`, `surface_prior.*`) must share the same spatial grid at `aux_resolution_m`
- `bands` must be a subset of `sensor_config.bands`
- `rt_model` must support all bands in `bands`

### 3.5 `SolvedAtmosphere` — Output of M5 (Aerosol Solver)

```python
@dataclass(frozen=True)
class SolvedAtmosphere:
    """Solver output: retrieved atmospheric parameters + diagnostics.

    Output contract of M5. Contains the solved AOT/TCWV fields at the
    aerosol retrieval resolution, plus the full AtmosphericState
    (with solved values merged in) for use by the corrector.
    """
    atmo_state: AtmosphericState    # full state with solved AOT/TCWV replacing priors
    aot: xr.DataArray               # solved AOT at aerosol resolution
    tcwv: xr.DataArray              # solved TCWV at aerosol resolution
    aot_unc: xr.DataArray           # posterior AOT uncertainty
    tcwv_unc: xr.DataArray          # posterior TCWV uncertainty

    # Diagnostics
    cost_final: float               # final cost function value
    n_iterations: int               # total optimizer iterations
    converged: bool                 # did the solver converge?
```

### 3.6 `CorrectionResult` — Output of M6 (Atmospheric Corrector)

```python
@dataclass(frozen=True)
class CorrectionResult:
    """Final output of atmospheric correction.

    Output contract of M6.
    """
    boa: xr.Dataset                 # BOA reflectance, one var per band, at native resolution
    boa_unc: xr.Dataset | None      # per-band uncertainty (optional)
    aot: xr.DataArray               # solved AOT map
    tcwv: xr.DataArray              # solved TCWV map
    cloud_mask: xr.DataArray        # final cloud mask
    metadata: dict[str, Any]        # processing metadata (timings, versions, etc.)
```

### 3.7 Contract Hierarchy Summary

| Contract Type | Produced by | Consumed by | Key fields |
|--------------|-------------|-------------|------------|
| `ObservationBundle` | M1 Preprocessor | M4 Grid Assembler, M6 Corrector | toa, geometry, cloud_mask, sensor_config |
| `AtmosphericState` | M2 Atmo Provider | M4 Grid Assembler | aot, tcwv, tco3, elevation + uncertainties |
| `SurfacePrior` | M3 Surface Provider | M4 Grid Assembler | boa, boa_unc, mask |
| `SolverInputBundle` | M4 Grid Assembler | M5 Solver | all upstream data, resampled + aligned |
| `SolvedAtmosphere` | M5 Solver | M6 Corrector | solved aot/tcwv + merged AtmosphericState |
| `CorrectionResult` | M6 Corrector | User / Output Writer | boa, uncertainties, diagnostics |

### 3.8 Contract Tests

**File**: `tests/unit/test_contracts.py`

Every contract type must have paired tests covering construction, immutability,
properties, and validation. When a new contract field or validation rule is added,
a corresponding test must be added in the same commit.

| Contract | Test | Assert |
|----------|------|--------|
| `ObservationBundle` | `test_obs_bundle_construction` | All fields accessible, types correct |
| | `test_obs_bundle_frozen` | Attribute assignment raises `FrozenInstanceError` |
| | `test_obs_bundle_metadata_time` | `metadata["observation_time"]` is `datetime` |
| | `test_obs_bundle_bounds_tuple` | `bounds` is 4-tuple of floats |
| `AtmosphericState` | `test_atmo_state_construction` | All 7 fields populated |
| | `test_atmo_state_with_updated_aot_tcwv` | New object returned, original unchanged, tco3 preserved |
| | `test_atmo_state_shape_consistency` | All fields share spatial shape |
| `SurfacePrior` | `test_surface_prior_construction` | `boa`, `boa_unc`, `mask` fields match |
| | `test_surface_prior_mask_dtype` | `mask` is boolean |
| | `test_surface_prior_boa_unc_nonneg` | `boa_unc >= 0` |
| `SolverInputBundle` | `test_sib_construction` | All upstream fields present |
| | `test_sib_bands_subset` | `bands ⊂ sensor_config.bands` |
| | `test_sib_resolution_metadata` | `aux_resolution_m > 0`, `aerosol_resolution_m > 0` |
| `SolvedAtmosphere` | `test_solved_atmo_construction` | `converged`, `cost_final`, `n_iterations` have correct types |
| | `test_solved_atmo_state_updated` | `atmo_state.aot` matches `aot` field |
| `CorrectionResult` | `test_correction_result_construction` | All fields present |
| | `test_correction_result_optional_unc` | `boa_unc=None` is valid |
| | `test_correction_result_boa_bands` | BOA dataset has expected band variables |

**File**: `tests/unit/test_validation.py`

Each `_validate_*` function gets a happy-path test (valid fixture passes) and
negative tests (one broken field each). When a new validation rule is added,
a matching negative test must be added.

| Validator | Test | Corrupted field | Expected error |
|-----------|------|----------------|----------------|
| `_validate_observation_bundle` | `test_validate_obs_happy` | None | No error |
| | `test_validate_obs_missing_time` | Remove `metadata["observation_time"]` | `"must include 'observation_time'"` |
| | `test_validate_obs_wrong_time_type` | `observation_time = "string"` | Type assertion |
| | `test_validate_obs_empty_toa` | `toa = xr.Dataset()` | `"must have spatial dimensions"` |
| | `test_validate_obs_cloud_shape` | `cloud_mask` shape ≠ TOA | Shape assertion |
| `_validate_atmospheric_state` | `test_validate_atmo_happy` | None | No error |
| | `test_validate_atmo_neg_aot_unc` | `aot_unc` contains -0.1 | `"non-negative"` |
| | `test_validate_atmo_neg_tcwv_unc` | `tcwv_unc` contains -0.5 | `"non-negative"` |
| `_validate_surface_prior` | `test_validate_prior_happy` | None | No error |
| | `test_validate_prior_shape_mismatch` | `boa.shape ≠ boa_unc.shape` | Shape assertion |
| | `test_validate_prior_mask_broadcast` | Non-broadcastable `mask` | Broadcast error |
| `_validate_solver_input_bundle` | `test_validate_sib_happy` | None | No error |
| | `test_validate_sib_misaligned` | `toa` grid ≠ `atmo_prior.aot` grid | Grid alignment error |
| | `test_validate_sib_bands_invalid` | `bands` not in `sensor_config` | Subset violation |

---

## 4. Module Definitions

Each module is defined by:
1. **Its callable signature** — what inputs it takes
2. **Its output contract** — the exact type it must return
3. **A Protocol** (for class-based implementations) or a `Callable` type alias (for functional implementations)

### 4.1 M1: Satellite Preprocessor

**Purpose**: hide sensor file formats behind a shared interface. Produce an `ObservationBundle`.

**Callable signature** (what the pipeline calls):
```python
PreprocessorFn = Callable[[Path, AOI | None], ObservationBundle]
```

The pipeline accepts anything that satisfies this signature — a function, a bound method,
or a class instance. **Built-in preprocessors are small classes** because they hold
sensor-specific configuration and expose overridable methods for fine-grained customisation:

```python
class Sentinel2Preprocessor:
    """S2 preprocessor. Each method can be overridden independently."""

    def load_toa(self, input_path: Path) -> xr.Dataset: ...
    def extract_geometry(self, input_path: Path) -> GeometryAngles: ...
    def extract_cloud_mask(self, input_path: Path) -> xr.DataArray: ...
    def get_metadata(self, input_path: Path) -> dict[str, Any]: ...
    def build_sensor_config(self, metadata: dict) -> SensorConfig: ...

    def preprocess(self, input_path: Path, aoi: AOI | None = None) -> ObservationBundle:
        """Default implementation: calls the methods above and assembles the bundle."""
        metadata = self.get_metadata(input_path)
        toa = self.load_toa(input_path)
        geometry = self.extract_geometry(input_path)
        cloud_mask = self.extract_cloud_mask(input_path)
        sensor_config = self.build_sensor_config(metadata)
        crs = str(toa.rio.crs)
        bounds = tuple(toa.rio.bounds())
        if aoi is not None:
            toa = aoi.clip(toa)
            bounds = aoi.get_bounds(crs)
        return ObservationBundle(toa, geometry, cloud_mask, sensor_config, metadata, crs, bounds)
```

**Why a class here?** The preprocessor has multiple cooperating methods that users
want to override independently (e.g. swap cloud mask, keep everything else).
An IDE can show all overridable methods. Stack traces show
`Sentinel2Preprocessor.extract_cloud_mask`, not an opaque `functools.partial`.

**Design rules**:
- The pipeline calls `preprocessor.preprocess(path, aoi)` and receives a complete `ObservationBundle`.
- Users can subclass and override individual methods (e.g. custom cloud masking) while inheriting the rest.
- Users can also supply a simple function `(Path, AOI) -> ObservationBundle` instead of a class.
- The pipeline never touches SAFE/MTL/COG specifics — only the `ObservationBundle`.

**Override a single method**:
```python
class MySentinel2(Sentinel2Preprocessor):
    def extract_cloud_mask(self, input_path: Path) -> xr.DataArray:
        return my_ml_model.predict(input_path)  # custom ML cloud mask

siac_process(config, input_path, preprocessor=MySentinel2())
```

**Or supply a standalone function**:
```python
def my_loader(input_path: Path, aoi=None) -> ObservationBundle:
    ...  # load from COG, STAC, etc.

siac_process(config, input_path, preprocessor=my_loader)
```

**Tests** (`tests/unit/test_preprocessor.py`):

| Test | What | Assert |
|------|------|--------|
| `test_preprocessor_returns_obs_bundle` | Call `.preprocess()` with mock SAFE path | Returns `ObservationBundle` with all fields |
| `test_preprocessor_aoi_clips_toa` | Pass AOI covering half the scene | `obs.bounds` matches AOI, TOA cropped |
| `test_preprocessor_no_aoi_full_extent` | Pass `aoi=None` | `obs.bounds` equals full image extent |
| `test_preprocessor_subclass_override` | Subclass with custom `extract_cloud_mask()` | Overridden method called, other methods inherited |
| `test_preprocessor_function_injection` | Pass a plain function as preprocessor | Pipeline accepts it, receives valid bundle |
| `test_preprocessor_metadata_has_time` | Check returned bundle | `metadata["observation_time"]` is `datetime` |
| `test_preprocessor_geometry_radians` | Check returned geometry | All angles in radians (0–2π range) |

### 4.2 M2: Atmospheric Prior Provider

**Purpose**: provide atmospheric state (AOT, TCWV, TCO3 + uncertainties + elevation) over the AOI.

**Callable signature** (what the pipeline calls):
```python
AtmoPriorFn = Callable[
    [tuple[float, float, float, float], str, datetime, float],
    AtmosphericState
]
```

**Built-in providers are small classes** because they hold configured state
(data directories, auth credentials, caching):

```python
class CAMSProvider:
    """Atmospheric prior from local CAMS NetCDF archive."""
    def __init__(self, cams_dir: Path):
        self.cams_dir = cams_dir

    def get_prior(self, bounds, crs, obs_time, resolution) -> AtmosphericState:
        """Load and interpolate CAMS data over the AOI."""
        ...

class MERRA2Provider:
    """Atmospheric prior via earthaccess (NASA)."""
    def __init__(self, cache_dir: Path | None = None):
        self.cache_dir = cache_dir

    def get_prior(self, bounds, crs, obs_time, resolution) -> AtmosphericState:
        ...
```

**Why a class here?** Each provider holds real configuration — directory paths,
credentials, cached file handles. `CAMSProvider(cams_dir="/data/cams")` then
`.get_prior(...)` is cleaner than `partial(get_cams_prior, cams_dir="/data/cams")`
because the IDE shows the constructor args and the method signature separately.

**User injection example** (a function is fine when there's no state):
```python
# User provides a constant atmosphere (e.g. for testing / controlled experiments)
def constant_atmo(bounds, crs, obs_time, resolution) -> AtmosphericState:
    shape = compute_grid_shape(bounds, crs, resolution)
    return AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.2)),
        tcwv=xr.DataArray(np.full(shape, 1.5)),
        tco3=xr.DataArray(np.full(shape, 0.34)),
        aot_unc=xr.DataArray(np.full(shape, 0.5)),
        tcwv_unc=xr.DataArray(np.full(shape, 0.5)),
        tco3_unc=xr.DataArray(np.full(shape, 0.05)),
        elevation=xr.DataArray(np.zeros(shape)),
    )

siac_process(config, input_path, atmo_provider=constant_atmo)
```

**Tests** (`tests/unit/test_atmo_provider.py`):

| Test | What | Assert |
|------|------|--------|
| `test_cams_provider_returns_atmo_state` | `CAMSProvider(tmp_dir).get_prior(bounds, crs, time, res)` | Returns `AtmosphericState` with all 7 fields |
| `test_cams_provider_aoi_scoped` | Check spatial extent of returned arrays | Covers requested bounds, no over-fetch |
| `test_merra2_provider_returns_atmo_state` | `MERRA2Provider().get_prior(...)` | Returns valid `AtmosphericState` |
| `test_atmo_provider_uncertainties_nonneg` | Check all `*_unc` fields | Values ≥ 0 everywhere |
| `test_atmo_provider_function_injection` | Pass a plain function `(b,c,t,r) -> AtmosphericState` | Pipeline accepts it |
| `test_constant_atmo_function` | Call `constant_atmo()` from user example | Returns valid `AtmosphericState` with uniform values |

### 4.3 M3: Surface Prior Provider

**Purpose**: provide surface reflectance prior (BOA + uncertainty + mask) over the AOI.

**Callable signature** (what the pipeline calls):
```python
SurfacePriorFn = Callable[
    [tuple[float, float, float, float], str, datetime, SensorConfig, GeometryAngles, float],
    SurfacePrior
]
```

**Why geometry is an input**: the surface prior depends on viewing geometry when derived from BRDF kernels. Pre-built stores may ignore this argument.

**Built-in providers are small classes** (same rationale as M2 — they hold
configured paths):

```python
class BRDFDerivedPriorProvider:
    """Derive surface prior from MODIS MCD43 BRDF kernels + viewing geometry."""
    def __init__(self, mcd43_dir: Path | None = None):
        self.mcd43_dir = mcd43_dir

    def get_surface_prior(self, bounds, crs, obs_time, sensor_config, geometry, resolution) -> SurfacePrior:
        ...

class PrebuiltPriorStore:
    """Load pre-built surface prior from Zarr store."""
    def __init__(self, store_path: Path):
        self.store_path = store_path

    def get_surface_prior(self, bounds, crs, obs_time, sensor_config, geometry, resolution) -> SurfacePrior:
        ...
```

**User injection example** (function — fine when no state needed):
```python
def my_prior(bounds, crs, obs_time, sensor_config, geometry, resolution) -> SurfacePrior:
    # load from custom dataset, interpolate, convolve to sensor bands...
    return SurfacePrior(boa=my_boa, boa_unc=my_unc, mask=my_mask)

siac_process(config, input_path, surface_prior_provider=my_prior)
```

**Tests** (`tests/unit/test_surface_prior_provider.py`):

| Test | What | Assert |
|------|------|--------|
| `test_brdf_provider_returns_surface_prior` | `BRDFDerivedPriorProvider(tmp_dir).get_surface_prior(...)` | Returns `SurfacePrior` with `boa`, `boa_unc`, `mask` |
| `test_brdf_provider_geometry_used` | Compare output with different geometry angles | BRDF-derived prior varies with geometry |
| `test_prebuilt_store_returns_surface_prior` | `PrebuiltPriorStore(zarr_path).get_surface_prior(...)` | Returns valid `SurfacePrior` |
| `test_prebuilt_store_ignores_geometry` | Same store, different geometry | Output unchanged |
| `test_surface_prior_aoi_scoped` | Check spatial extent | Covers requested bounds |
| `test_surface_prior_function_injection` | Pass plain function | Pipeline accepts it, receives valid `SurfacePrior` |

### 4.4 M4: Grid Assembler

**Purpose**: resample and align all upstream outputs to the solver's working grids. Produce a `SolverInputBundle`.

**Function signature**:
```python
def assemble_grids(
    obs: ObservationBundle,
    atmo: AtmosphericState,
    surface: SurfacePrior,
    rt_model: RTModelBackend,
    aux_resolution_m: float = 500.0,
    aerosol_resolution_m: float = 1000.0,
) -> SolverInputBundle:
    """Resample and align all upstream outputs to solver grids."""
    ...

# Type alias for user-provided alternatives
GridAssemblerFn = Callable[
    [ObservationBundle, AtmosphericState, SurfacePrior, RTModelBackend, float, float],
    SolverInputBundle
]
```

**Responsibilities**:
- Resample TOA, geometry, cloud mask to aux resolution (area-weighted mean for reflectance, bilinear for angles, conservative OR for cloud mask)
- Resample atmospheric state to aux resolution
- Resample surface prior to aux resolution
- Select solver bands from `sensor_config` by wavelength
- Validate that all outputs are spatially aligned

**Design note**: Most users will use the built-in `assemble_grids()`. Custom assemblers are useful for non-standard grids or when the user wants to skip resampling (e.g. all data already at the right resolution).

**Tests** (`tests/unit/test_grid_assembler.py`):

| Test | What | Assert |
|------|------|--------|
| `test_assemble_grids_output_type` | Call with valid upstream outputs | Returns `SolverInputBundle` |
| `test_assemble_grids_spatial_alignment` | Check all raster fields | Share the same spatial grid |
| `test_assemble_grids_aux_resolution` | Check `aux_resolution_m` | Matches configured value (default 500.0) |
| `test_assemble_grids_band_selection` | Check `bands` field | Only aerosol-sensitive bands (400–520 nm) selected |
| `test_assemble_grids_cloud_mask_conservative` | Any native pixel is cloud | Aux pixel is cloud (conservative OR) |
| `test_assemble_grids_geometry_bilinear` | Angles with gradient | Smooth gradient preserved (not nearest-neighbour) |
| `test_assemble_grids_identity_same_res` | All inputs at aux resolution already | Output ≈ input (no resampling artefacts) |
| `test_assemble_grids_multi_res_input` | TOA with mixed 10/20/60 m bands | All resampled to aux without error |
| `test_assemble_grids_passes_validation` | Run `_validate_solver_input_bundle()` on output | No error |

### 4.5 M5: Aerosol Solver

**Purpose**: retrieve AOT and TCWV by minimising the cost function.

**Function signature**:
```python
def solve_aerosol(
    inputs: SolverInputBundle,
    config: SolverConfig,
) -> SolvedAtmosphere:
    """Multi-grid L-BFGS-B solver (default implementation)."""
    ...

# Type alias for user-provided alternatives
SolverFn = Callable[[SolverInputBundle, SolverConfig], SolvedAtmosphere]
```

**Key design rule**: the solver receives a **single `SolverInputBundle`** — it never fetches or aligns data itself. All data is pre-assembled and validated.

**Built-in functions**: `solve_aerosol()` (multi-grid L-BFGS-B), `solve_aerosol_pixelwise()` (per-pixel)

**User injection example**:
```python
# User bypasses the solver entirely with pre-computed AOT
def precomputed_solver(inputs: SolverInputBundle, config: SolverConfig) -> SolvedAtmosphere:
    precomputed_aot = xr.open_dataarray("my_aot.tif")
    return SolvedAtmosphere(
        atmo_state=inputs.atmo_prior.with_updated_aot_tcwv(precomputed_aot, inputs.atmo_prior.tcwv),
        aot=precomputed_aot,
        tcwv=inputs.atmo_prior.tcwv,
        aot_unc=inputs.atmo_prior.aot_unc,
        tcwv_unc=inputs.atmo_prior.tcwv_unc,
        cost_final=0.0, n_iterations=0, converged=True,
    )

siac_process(config, input_path, solver=precomputed_solver)
```

**Tests** (`tests/unit/test_solver.py` — extend existing):

| Test | What | Assert |
|------|------|--------|
| `test_solve_aerosol_output_type` | Call with mock `SolverInputBundle` | Returns `SolvedAtmosphere` |
| `test_solve_aerosol_aot_positive` | Check solved AOT | ≥ 0 everywhere |
| `test_solve_aerosol_tcwv_positive` | Check solved TCWV | ≥ 0 everywhere |
| `test_solve_aerosol_convergence` | Sufficient iterations on clean scene | `converged=True` |
| `test_solve_aerosol_state_updated` | Compare `atmo_state.aot` vs prior | Solver changed it |
| `test_solve_aerosol_known_answer` | Pre-computed forward-model scene | AOT within 20% of truth |
| `test_precomputed_solver_passthrough` | Inject pre-solved AOT function | Output uses supplied AOT directly |

### 4.6 M6: Atmospheric Corrector

**Purpose**: apply atmospheric correction (TOA → BOA) using the solved atmospheric state.

**Function signature**:
```python
def correct_atmosphere(
    obs: ObservationBundle,
    solved: SolvedAtmosphere,
    rt_model: RTModelBackend,
) -> CorrectionResult:
    """Apply atmospheric correction at native band resolutions."""
    ...

# Type alias for user-provided alternatives
CorrectorFn = Callable[[ObservationBundle, SolvedAtmosphere, RTModelBackend], CorrectionResult]
```

**Key design rule**: the corrector operates at **native band resolutions**. It upsamples the solved AOT/TCWV (at aerosol resolution) to each band's native resolution.

**Tests** (`tests/unit/test_correction.py` — extend existing):

| Test | What | Assert |
|------|------|--------|
| `test_correct_output_type` | Call with `ObservationBundle` + `SolvedAtmosphere` | Returns `CorrectionResult` |
| `test_correct_boa_bands_match_toa` | Check BOA band names | Same as input TOA |
| `test_correct_boa_range` | Check BOA values | In [–0.05, 1.5] (minor negatives allowed) |
| `test_correct_cloud_mask_preserved` | Check `result.cloud_mask` | Matches original from M1 |
| `test_correct_native_resolution` | Check BOA spatial shape | Matches M1 native TOA shape |
| `test_correct_metadata_timing` | Check `result.metadata` | Contains `"processing_time_s"` |

### 4.7 RT Model Backend (Cross-Cutting)

The RT model is not a data-sourcing module — it is a **compute service** used by both M5 and M6. It is always local (no I/O latency).

```python
class RTModelBackend(Protocol):
    def compute_coefficients(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
        compute_jacobian: bool = False,
    ) -> RTCoefficients: ...

    def supports_jacobian(self) -> bool: ...
    def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool: ...
```

**Tests** (`tests/unit/test_rt_backend.py`):

| Test | What | Assert |
|------|------|--------|
| `test_rt_backend_protocol_conformance` | `MockRTModel` satisfies `RTModelBackend` | Duck-typing check passes |
| `test_rt_backend_jacobian_shape` | `compute_coefficients(..., compute_jacobian=True)` | `d_xap` has `param` dim with `["aot", "tcwv"]` |
| `test_rt_backend_no_jacobian` | `compute_coefficients(..., compute_jacobian=False)` | `d_xap/d_xbp/d_xcp` are `None` |
| `test_rt_backend_unsupported_sensor` | `is_available_for_sensor("UNKNOWN", "SAT")` | Returns `False` |
| `test_rt_coefficients_apply_correction` | Call `apply_correction(toa)` | BOA = known analytical formula |

### 4.8 Satellite Data Access (Search + Download) — Generic Contract

**Purpose**: optionally resolve "logical inputs" (product IDs, tile/date shorthands, spatial/temporal queries) into a local path consumable by M1.

This is intentionally **not** part of the main pipeline because:
- some sensors are always local (no search/download)
- some sensors have multiple sources (S3, HTTP, public buckets)
- discovery rules differ substantially

Core contract (planned):

```python
class SatelliteDataBackend(Protocol):
    def search(self, query: SatelliteQuery) -> list[SatelliteProduct]: ...
    def download(self, product: SatelliteProduct, dest_dir: Path) -> Path: ...
```

Sensor plans define:
- query/product dataclasses (e.g. `S2Query`, `S2Product`)
- one or more concrete backends

---

## 5. Custom Provider Pattern (User-Injected Functions)

### 5.1 Design Principle

The pipeline entry point is a function that accepts optional overrides for each
module. Each override can be **a class instance** (for stateful providers) **or a
plain function** (for stateless transforms or simple user overrides):

```python
def siac_process(
    config: SIACConfig,
    input_path: Path,
    *,
    aoi: AOI | None = None,
    preprocessor: PreprocessorFn | Sentinel2Preprocessor | None = None,
    atmo_provider: AtmoPriorFn | CAMSProvider | None = None,
    surface_prior_provider: SurfacePriorFn | BRDFDerivedPriorProvider | None = None,
    grid_assembler: GridAssemblerFn | None = None,
    solver: SolverFn | None = None,
    corrector: CorrectorFn | None = None,
    rt_model: RTModelBackend | None = None,
) -> CorrectionResult:
    """Run the full SIAC pipeline.

    Stateful modules (M1, M2, M3): pass a class instance or a function.
    Stateless modules (M4, M5, M6): pass a function (the default is usually fine).
    """
    ...
```

Internally the pipeline normalises each argument to a callable:
if a class instance is passed (e.g. `CAMSProvider`), the pipeline calls its
well-known method (`.get_prior()`). If a bare function is passed, it's called
directly. Each `None` argument is resolved to the config-driven default.

### 5.2 Resolution Order

For each module, the pipeline resolves the implementation in this order:

1. **Explicit injection** (constructor argument) — highest priority
2. **Config-driven factory** (e.g. `config.atmo_prior.provider = "cams"`) — default
3. **Auto-detection** (e.g. sensor auto-detect) — fallback

### 5.3 Contract Validation

Each module output is validated at runtime before being passed downstream:

```python
def _validate_observation_bundle(obs: ObservationBundle) -> None:
    """Validate M1 output before passing to M4."""
    assert "observation_time" in obs.metadata, "metadata must include 'observation_time'"
    assert isinstance(obs.metadata["observation_time"], datetime)
    assert obs.toa.sizes, "toa must have spatial dimensions"
    assert obs.cloud_mask.shape == _spatial_shape(obs.toa)

def _validate_atmospheric_state(atmo: AtmosphericState) -> None:
    """Validate M2 output before passing to M4."""
    for field_val in [atmo.aot_unc, atmo.tcwv_unc, atmo.tco3_unc]:
        assert (field_val.values >= 0).all(), "uncertainties must be non-negative"
```

Validation is applied automatically by the pipeline orchestrator. Users get clear error messages if their custom provider returns an invalid contract.

**Tests** (`tests/integration/test_injection.py`):

| Test | What | Assert |
|------|------|--------|
| `test_inject_preprocessor_function` | Pass plain `(Path, AOI) -> ObservationBundle` function | Pipeline calls it, returns valid result |
| `test_inject_atmo_function` | Pass `constant_atmo(b,c,t,r) -> AtmosphericState` | Pipeline uses it instead of config default |
| `test_inject_surface_prior_function` | Pass custom function returning `SurfacePrior` | Pipeline uses it |
| `test_inject_solver_function` | Pass passthrough solver (fixed AOT) | Solver bypassed, correction uses fixed AOT |
| `test_inject_corrector_function` | Pass identity corrector (BOA = TOA) | Output BOA ≈ input TOA |
| `test_inject_preprocessor_class` | Pass `MySensor().preprocess` as bound method | Pipeline calls bound method |
| `test_inject_subclassed_preprocessor` | Subclass with overridden `extract_cloud_mask` | Overridden method called, rest inherited |
| `test_inject_bad_preprocessor` | Return incomplete `ObservationBundle` | Validation error with clear message |
| `test_inject_bad_atmo_wrong_type` | Return dict instead of `AtmosphericState` | `TypeError` or `AttributeError` |
| `test_inject_partial_override` | Override only M2, rest default | Only M2 custom called, M1/M3 from config |

### 5.4 Example: Full Custom Pipeline

```python
from siac import siac_process
from siac.core.types import ObservationBundle, AtmosphericState, SurfacePrior

# User provides all three data-sourcing modules as simple functions
def my_preprocessor(input_path, aoi=None) -> ObservationBundle:
    ...  # custom satellite format reader

def my_atmo(bounds, crs, obs_time, resolution) -> AtmosphericState:
    ...  # custom atmospheric model

def my_surface(bounds, crs, obs_time, sensor_config, geometry, resolution) -> SurfacePrior:
    ...  # custom surface reflectance model

result = siac_process(
    config=SIACConfig(sensor="auto"),
    input_path=Path("/path/to/data"),
    preprocessor=my_preprocessor,
    atmo_provider=my_atmo,
    surface_prior_provider=my_surface,
)
```

### 5.5 Example: Partial Override (Only Replace One Module)

```python
# Use all defaults except atmospheric prior
result = siac_process(
    config=SIACConfig.from_yaml("config.yaml"),
    input_path=Path("/path/to/S2_SAFE/"),
    atmo_provider=my_custom_atmo_provider,
)
```

### 5.6 Example: Skip the Solver (Use Pre-Solved AOT)

```python
# Provide pre-retrieved AOT — the solver is bypassed
def passthrough_solver(inputs, config) -> SolvedAtmosphere:
    aot_map = xr.open_dataarray("my_aot_retrieval.tif")
    return SolvedAtmosphere(
        atmo_state=inputs.atmo_prior.with_updated_aot_tcwv(aot_map, inputs.atmo_prior.tcwv),
        aot=aot_map, tcwv=inputs.atmo_prior.tcwv,
        aot_unc=inputs.atmo_prior.aot_unc, tcwv_unc=inputs.atmo_prior.tcwv_unc,
        cost_final=0.0, n_iterations=0, converged=True,
    )

result = siac_process(config, input_path, solver=passthrough_solver)
```

---

## 6. AOI Scoping Rules (Global)

**AOI is optional.** The following rules apply throughout the entire pipeline:

1. **If AOI is provided** (via `config.aoi` as GeoJSON path, WKT, or bounding box):
   - The input satellite TOA data is **cropped to the AOI extent** before processing.
   - All auxiliary data (atmospheric priors, BRDF/surface prior, DEM, water mask) is **retrieved/read only over the AOI extent**.
2. **If AOI is not provided**:
   - The AOI defaults to the **full extent of the input satellite image** (derived from the loaded TOA grid).
   - All auxiliary data is still scoped to this image-derived AOI — no global reads.
3. **In both cases**:
   - providers receive AOI bounds and must only fetch/crop data within that region
   - all returned rasters must be reprojectable/alignable to the pipeline grids

This ensures efficient data access (no over-fetching) and consistent spatial alignment.

### 6.1 Planned Work (earthaccess + AOI wiring)

**Status**: Planned (pre-implementation)

The target implementation introduces:

| File | Action | Description |
|------|--------|-------------|
| `python/siac/core/aoi.py` | NEW | AOI wrapper constructors + AOI → WGS84 bbox export |
| `python/siac/io/earthaccess_source.py` | NEW | earthaccess auth/search/open wrapper |
| `python/siac/priors/brdf/mcd43_earthaccess.py` | NEW | MCD43A1 BRDF provider via earthaccess |
| `python/siac/priors/atmospheric/merra2.py` | NEW | MERRA-2 atmospheric provider via earthaccess |
| `python/siac/priors/brdf/gee_stub.py` | NEW | GEE placeholder |
| `python/siac/siac.py` | MODIFY | Wire AOI + new providers; fix surface-prior provider plumbing |
| `python/siac/core/config.py` | MINOR | Ensure `merra2` provider works end-to-end |

Key design decisions:
- AOI is a first-class object passed through the pipeline; it wraps existing geometry utilities.
- AOI is optional: if explicit AOI is provided, TOA is clipped; otherwise AOI = image extent.
- All auxiliary data reads/fetches are AOI-scoped regardless of whether TOA was clipped.
- `EarthAccessSource` does lazy auth and standardises `search_granules()` / `open_granules()`.

**Tests** (`tests/unit/test_aoi.py`, `tests/unit/test_io.py` — extend existing):

| Test | What | Assert |
|------|------|--------|
| `test_aoi_from_bbox` | Construct AOI from `(xmin, ymin, xmax, ymax)` | Valid AOI with correct bounds |
| `test_aoi_from_geojson` | Construct AOI from GeoJSON file | Bounds match geometry envelope |
| `test_aoi_to_wgs84_bbox` | Export AOI to WGS84 | Coordinates in ±180/±90 |
| `test_aoi_clip_toa` | `aoi.clip(toa_dataset)` | Clipped TOA extent matches AOI |
| `test_aoi_none_uses_image_extent` | No AOI provided | `obs.bounds` equals full image extent |
| `test_earthaccess_source_lazy_auth` | Create `EarthAccessSource` | No auth until first `search_granules()` call |
| `test_merra2_provider_aoi_scoped` | Pass AOI bounds to MERRA-2 provider | Only fetches data within bounds |
| `test_mcd43_provider_aoi_scoped` | Pass AOI bounds to MCD43 provider | Only fetches tiles covering AOI |

---

## 7. Pipeline Orchestration (Sequential → Parallel)

### 7.1 Sensor-Agnostic Pipeline Flow

Expressed in terms of modules and contracts:

```
1. M1: preprocess(input_path, aoi)                                      → ObservationBundle
2. M2: get_prior(bounds, crs, time, res)                                → AtmosphericState
3. M3: get_surface_prior(bounds, crs, time, sensor_config, geometry, res) → SurfacePrior
4. M4: assemble(obs, atmo, surface, rt_model)                           → SolverInputBundle
5. M5: solve(inputs, solver_config)                                      → SolvedAtmosphere
6. M6: correct(obs, solved, rt_model)                                    → CorrectionResult
```

Steps 1–3 are independent and can execute concurrently. Steps 4–6 are sequential.

### 7.2 Parallel Fetch Groups

The pipeline supports a parallel I/O model mapped to module calls:

- **Group A (M1)**: `preprocess()` → `ObservationBundle` (may be split internally into fast metadata + slow band loading)
- **Group B (M2)**: `get_prior()` → `AtmosphericState`
- **Group C (M3)**: `get_surface_prior()` → `SurfacePrior`

**Design rule**: Groups A/B/C are independent modules returning independent contracts. They are fetched concurrently.

**Note**: M3 may depend on geometry from M1 for BRDF-based priors. In this case, M1 provides a fast metadata-only path:
```python
# M1 fast path: extract just geometry + metadata (< 1s)
obs_meta = preprocessor.get_metadata_fast(input_path)
# M3 can start with geometry from obs_meta while M1 continues loading TOA bands
```

### 7.3 Input Requirements Per Processing Stage

#### Stage 1: Aerosol & Water Vapour Retrieval (M5 Solver)

The solver minimises the cost function J = J_obs + J_prior + J_smooth to retrieve **AOT** and **TCWV**.

It receives a single `SolverInputBundle` containing:

| Field | Type | Original Source |
|-------|------|----------------|
| `toa` | `xr.DataArray` | M1 → M4 (resampled) |
| `geometry` | `GeometryAngles` | M1 → M4 (resampled) |
| `cloud_mask` | `xr.DataArray` | M1 → M4 (resampled) |
| `sensor_config` | `SensorConfig` | M1 (passthrough) |
| `bands` | `list[SensorBand]` | M4 (wavelength-selected) |
| `atmo_prior` | `AtmosphericState` | M2 → M4 (resampled) |
| `surface_prior` | `SurfacePrior` | M3 → M4 (resampled) |
| `rt_model` | `RTModelBackend` | Config (local) |

#### Stage 2: Atmospheric Correction (M6 Corrector)

| Field | Type | Source |
|-------|------|--------|
| `obs` | `ObservationBundle` | M1 (original, native resolution) |
| `solved` | `SolvedAtmosphere` | M5 output |
| `rt_model` | `RTModelBackend` | Config (local) |

The corrector upsamples solved AOT/TCWV to native band resolutions for per-pixel correction.

### 7.4 Dependency Analysis

```
  TIME ──────────────────────────────────────────────────────────────▶

  M1 ▐██████████████████████████████▌ preprocess → ObservationBundle  ~30s
  M2 ▐████████████████▌ get_prior → AtmosphericState                 ~15s
  M3 ▐████████████████████▌ get_surface_prior → SurfacePrior          ~20s
                   │
    ┌──────────────▼────────────────┐
    │ M4: assemble → SolverInputBundle │      ← starts at max(M1, M2, M3)
    └──────────────┬────────────────┘
                   │
    ┌──────────────▼────────────────┐
    │ M5: solve → SolvedAtmosphere  │
    └──────────────┬────────────────┘
                   │
    ┌──────────────▼────────────────┐
    │ M6: correct → CorrectionResult│
    └───────────────────────────────┘
```

### 7.5 Implementation: `run_pipeline()` Function

Planned implementation lives in `python/siac/pipeline.py`:

```python
from concurrent.futures import ThreadPoolExecutor

def run_pipeline(
    input_path: Path,
    aoi: AOI | None,
    config: SIACConfig,
    *,
    preprocessor: PreprocessorFn,
    atmo_provider: AtmoPriorFn,
    surface_prior_provider: SurfacePriorFn,
    grid_assembler: GridAssemblerFn,
    solver: SolverFn,
    corrector: CorrectorFn,
    rt_model: RTModelBackend,
    max_workers: int = 4,
) -> CorrectionResult:
    """Orchestrate module execution with concurrent data sourcing.

    This is a plain function — no class, no state. Each module callable
    is passed as an argument (either a bound method or a plain function).
    Easy to test, debug, and trace.
    """
    # Phase 1: Preprocess (needs to complete first to get bounds)
    obs = preprocessor(input_path, aoi)
    _validate_observation_bundle(obs)

    bounds = obs.bounds
    crs = obs.crs
    obs_time = obs.metadata["observation_time"]
    resolution = config.solver.aerosol_resolution

    # Phase 2: Concurrent data sourcing (M2, M3 are independent)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        f_m2 = executor.submit(atmo_provider, bounds, crs, obs_time, resolution)
        f_m3 = executor.submit(
            surface_prior_provider,
            bounds, crs, obs_time, obs.sensor_config, obs.geometry, resolution,
        )
        atmo = f_m2.result()
        surface = f_m3.result()

    _validate_atmospheric_state(atmo)
    _validate_surface_prior(surface)

    # Phase 3: Sequential processing
    solver_inputs = grid_assembler(obs, atmo, surface, rt_model)
    solved = solver(solver_inputs, config.solver)
    result = corrector(obs, solved, rt_model)

    return result
```

Note: no `self`, no stored state. Every dependency is an explicit argument. A stack
trace from any failure points directly to the failing function — no method-resolution
hunting.

**Convenience wrapper** (`python/siac/siac.py`):

```python
def siac_process(
    config: SIACConfig,
    input_path: Path,
    *,
    aoi: AOI | None = None,
    preprocessor: PreprocessorFn | None = None,
    atmo_provider: AtmoPriorFn | None = None,
    surface_prior_provider: SurfacePriorFn | None = None,
    grid_assembler: GridAssemblerFn | None = None,
    solver: SolverFn | None = None,
    corrector: CorrectorFn | None = None,
    rt_model: RTModelBackend | None = None,
) -> CorrectionResult:
    """Public entry point. Resolves defaults from config, then calls run_pipeline()."""
    preprocessor = preprocessor or _resolve_preprocessor(config)
    atmo_provider = atmo_provider or _resolve_atmo_provider(config)
    surface_prior_provider = surface_prior_provider or _resolve_surface_prior_provider(config)
    grid_assembler = grid_assembler or assemble_grids
    solver = solver or solve_aerosol
    corrector = corrector or correct_atmosphere
    rt_model = rt_model or _resolve_rt_model(config)

    return run_pipeline(
        input_path, aoi, config,
        preprocessor=preprocessor,
        atmo_provider=atmo_provider,
        surface_prior_provider=surface_prior_provider,
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
    )
```

The `_resolve_*` helpers are plain functions that read `config` and return
configured class instances or functions:

```python
def _resolve_atmo_provider(config: SIACConfig) -> AtmoPriorFn:
    """Return a callable (bounds, crs, time, res) -> AtmosphericState."""
    match config.atmo_prior.provider:
        case "cams":
            provider = CAMSProvider(cams_dir=config.atmo_prior.cams_dir)
            return provider.get_prior
        case "merra2":
            provider = MERRA2Provider(cache_dir=config.atmo_prior.cache_dir)
            return provider.get_prior
        case _:
            raise ValueError(f"Unknown atmo provider: {config.atmo_prior.provider}")

def _resolve_preprocessor(config: SIACConfig) -> PreprocessorFn:
    match config.sensor:
        case "sentinel2":
            return Sentinel2Preprocessor().preprocess
        case "landsat89":
            return Landsat89Preprocessor().preprocess
        case _:
            raise ValueError(f"Unknown sensor: {config.sensor}")
```

The pattern is: **construct the class** (configuring state once), then **return its method**
as a plain callable. The pipeline sees only `Callable[..., SomeContract]`.

### 7.6 Multi-Resolution Data Assembly (M4 Responsibility)

SIAC operates across three spatial scales. The Grid Assembler (M4) is solely
responsible for aligning data between them:

- **Native grid**: sensor-native per-band resolution (input to/output from M6)
- **Aux grid**: common analysis grid (e.g. 500 m) for solver inputs
- **Aerosol grid**: coarser retrieval grid (e.g. 1000 m) for AOT/TCWV

Core requirement: M4 produces a `SolverInputBundle` where all raster data is
on the same grid. M6 receives native-resolution data from M1 and upsamples solved
fields.

### 7.7 Verification Plan (Orchestration)

**Tests** (`tests/integration/test_pipeline.py` — extend existing):

**`run_pipeline()` — Happy path & call order:**

| Test | What | Assert |
|------|------|--------|
| `test_pipeline_happy_path` | All six mock modules → `CorrectionResult` | No error, valid result type |
| `test_pipeline_calls_all_modules` | Wrap each mock with call counter | Each mock called exactly once |
| `test_pipeline_module_call_order` | Record call timestamps | M1 before M4; M4 before M5; M5 before M6 |
| `test_pipeline_m2_m3_after_m1` | Record call times | M2 and M3 start after M1 returns (needs bounds) |
| `test_pipeline_returns_correction_result` | Check returned type | `isinstance(result, CorrectionResult)` |

**`run_pipeline()` — Validation integration:**

| Test | What | Assert |
|------|------|--------|
| `test_pipeline_invalid_m1_output` | M1 returns bundle missing `observation_time` | Clear error before M4 starts |
| `test_pipeline_invalid_m2_output` | M2 returns negative uncertainties | Clear error before M4 starts |
| `test_pipeline_invalid_m3_output` | M3 returns shape mismatch | Clear error before M4 starts |
| `test_pipeline_m1_exception` | M1 raises `FileNotFoundError` | Exception propagates with useful traceback |
| `test_pipeline_m2_exception` | M2 raises `ConnectionError` | Exception propagates, not swallowed by thread |

**`siac_process()` — Config resolution:**

| Test | What | Assert |
|------|------|--------|
| `test_siac_process_resolves_defaults` | Call with config only (no overrides) | `_resolve_*` helpers called, pipeline runs |
| `test_siac_process_explicit_override` | Pass custom `atmo_provider` | Config default NOT called |
| `test_siac_process_none_means_default` | All args are `None` | Resolved to config defaults |

**`tests/unit/test_resolve.py` — `_resolve_*` helpers:**

| Test | What | Assert |
|------|------|--------|
| `test_resolve_preprocessor_sentinel2` | `config.sensor = "sentinel2"` | Returns callable wrapping `Sentinel2Preprocessor` |
| `test_resolve_preprocessor_unknown` | `config.sensor = "unknown"` | `ValueError` |
| `test_resolve_atmo_cams` | `config.atmo_prior.provider = "cams"` | Returns `CAMSProvider.get_prior` bound method |
| `test_resolve_atmo_merra2` | `config.atmo_prior.provider = "merra2"` | Returns `MERRA2Provider.get_prior` bound method |
| `test_resolve_atmo_unknown` | `config.atmo_prior.provider = "xxx"` | `ValueError` |
| `test_resolve_rt_model_emulator` | `config.rt_model.backend = "emulator"` | Returns `EmulatorBackend` instance |
| `test_resolve_returns_callable` | All `_resolve_*` helpers | Returned value is `callable()` |

**Concurrency tests** (`tests/benchmarks/test_perf.py`):

| Test | What | Assert |
|------|------|--------|
| `test_m2_m3_concurrent` | M2/M3 mocks each sleep 0.5 s | Total time < 0.8 s (parallel, not 1.0 s) |
| `test_pipeline_sequential_fallback` | `max_workers=1` | Pipeline still works, M2/M3 run sequentially |
| `test_concurrent_exception_handling` | M2 raises after 0.2 s | Exception propagated, M3 not swallowed |

---

## 8. Core Requirements for Adding New Sensors

This section is the checklist that every new sensor must satisfy to plug into the core pipeline.

### 8.1 Required: Implement M1 (Return `ObservationBundle`)

You must provide something that produces an `ObservationBundle`.

**Option A — Small class** (recommended when the sensor has multiple steps):
```python
class MySensorPreprocessor:
    def load_toa(self, input_path: Path) -> xr.Dataset: ...
    def extract_geometry(self, input_path: Path) -> GeometryAngles: ...
    def extract_cloud_mask(self, input_path: Path) -> xr.DataArray: ...
    def get_metadata(self, input_path: Path) -> dict[str, Any]: ...
    def build_sensor_config(self, metadata: dict) -> SensorConfig: ...

    def preprocess(self, input_path: Path, aoi: AOI | None = None) -> ObservationBundle:
        metadata = self.get_metadata(input_path)
        toa = self.load_toa(input_path)
        geometry = self.extract_geometry(input_path)
        cloud_mask = self.extract_cloud_mask(input_path)
        sensor_config = self.build_sensor_config(metadata)
        crs = str(toa.rio.crs)
        bounds = tuple(toa.rio.bounds())
        if aoi is not None:
            toa = aoi.clip(toa)
            bounds = aoi.get_bounds(crs)
        return ObservationBundle(toa, geometry, cloud_mask, sensor_config, metadata, crs, bounds)
```

**Option B — Plain function** (fine for simple sensors or quick prototyping):
```python
def preprocess_my_sensor(input_path: Path, aoi: AOI | None = None) -> ObservationBundle:
    ...
    return ObservationBundle(toa, geometry, cloud_mask, sensor_config, metadata, crs, bounds)
```

Either way, keep methods/functions small and independently testable.
Do **not** build deep class hierarchies — one concrete class per sensor is enough.

Requirements for the `ObservationBundle`:
- `toa`: reflectance data in consistent units (typically scaled to [0,1]), with georeferencing
- `geometry`: SZA/SAA/VZA/VAA in **radians**, broadcastable to TOA spatial dims
- `cloud_mask`: boolean (True = cloud), matching TOA spatial shape
- `metadata`: must include `'observation_time'` as `datetime`
- `sensor_config`: band wavelength definitions + native resolution metadata
- `crs`, `bounds`: spatial reference for downstream AOI operations

### 8.2 Required: Provide a `SensorConfig` (wavelength-driven)

Core rule: the solver selects bands by **wavelength ranges**, not by sensor-specific names.

At minimum, the sensor must define:
- visible bands (400–700 nm) for surface prior targets
- aerosol-sensitive bands (typically 400–520 nm) used by the solver
- if available: NIR/SWIR bands (750–2500 nm) for optional surface-prior refinement

If a sensor has SRFs, they can be included for accurate spectral convolution.

### 8.3 Required: AOI Scoping Compatibility

The new sensor preprocessor must support:
- deriving bounds/CRS from TOA grid (populating `ObservationBundle.bounds` and `.crs`)
- cropping TOA when an explicit AOI is provided
- ensuring downstream modules can align data to the aux grid

### 8.4 Optional: Add a Data Access Module (Search + Download)

Only needed if you want `SIAC.process_*()` convenience APIs for remote discovery/download. This sits **before** M1 in the pipeline and is not part of the core module chain.

### 8.5 Required: Verification

For each new sensor:
- unit tests: `ObservationBundle` construction with sample data passes validation
- unit tests: `SensorConfig` band definitions, wavelength selection
- integration tests: preprocessor loads TOA + geometry + cloud mask on at least one sample input
- regression tests: the existing test suite still passes

---

## 9. Sensor-Agnostic Spectral Model & Surface Prior

**Status**: Planned (not yet implemented)

### 9.1 Motivation: Sensor-Agnostic Design

SIAC must work with any optical sensor — not just Sentinel-2 and Landsat.
Sensors range from broadband multispectral (6–13 bands) to hyperspectral (200+ bands).
The system should never reference hardcoded band names like `B8A` or `B11`; instead,
all band selection is driven by **wavelength** and (when available) **spectral response functions (SRFs)**.

### 9.2 Spectral Band Descriptor

Every input sensor describes its bands with a `SpectralBandDescriptor` that supports:

1. Gaussian approximation (multispectral): center wavelength + FWHM
2. Full SRF (hyperspectral / precision): tabulated (wavelength, response)

```python
@dataclass(frozen=True)
class SpectralBandDescriptor:
    name: str
    center_wavelength_nm: float
    fwhm_nm: float
    resolution_m: float
    srf_wavelengths_nm: np.ndarray | None = None
    srf_response: np.ndarray | None = None

    @property
    def has_srf(self) -> bool: ...

    @property
    def wavelength_um(self) -> float: ...
```

The existing `SensorBand` in `python/siac/core/types.py` is expected to be replaced/extended by `SpectralBandDescriptor`.

**Tests** (`tests/unit/test_spectral.py`):

| # | Test | Assert |
|---|------|--------|
| 1 | Gaussian-only construction (no SRF arrays) | `has_srf == False`, `center_wavelength_nm > 0` |
| 2 | Full-SRF construction (tabulated arrays) | `has_srf == True`, arrays non-empty |
| 3 | `wavelength_um` property | 550 nm → `≈ 0.55` |
| 4 | Frozen: mutating field raises | `FrozenInstanceError` |

### 9.3 SensorConfig: Wavelength-Driven Band Selection

The spectral model must support wavelength-driven selection.
Planned `SensorConfig` additions:

```python
@dataclass(frozen=True)
class SensorConfig:
    bands: tuple[SpectralBandDescriptor, ...]

    def select_bands_in_range(self, wl_min_nm: float, wl_max_nm: float) -> list[SpectralBandDescriptor]: ...
    def select_nearest_band(self, target_nm: float, tolerance_nm: float = 50.0) -> SpectralBandDescriptor | None: ...

    @property
    def vis_bands(self) -> list[SpectralBandDescriptor]: ...     # 400–700 nm
    @property
    def nir_bands(self) -> list[SpectralBandDescriptor]: ...     # 750–1000 nm
    @property
    def swir_bands(self) -> list[SpectralBandDescriptor]: ...    # 1000–2500 nm, excluding absorption windows
```

**Tests** (`tests/unit/test_spectral.py`):

| # | Test | Assert |
|---|------|--------|
| 1 | `vis_bands` returns only 400–700 nm bands | all returned bands in [400, 700] |
| 2 | `nir_bands` excludes visible | all returned bands in [750, 1000] |
| 3 | `select_nearest_band(550)` finds green | returned band centre ≈ 550 nm |
| 4 | `select_nearest_band(3000)` returns None | no match outside tolerance |
| 5 | Empty range returns empty list | `select_bands_in_range(9000, 9999) == []` |

### 9.4 Spectral Regions (Sensor-Agnostic)

All algorithms operate on spectral regions defined by wavelength boundaries:

| Region | Wavelength (nm) | Role in SIAC |
|--------|------------------|--------------|
| Coastal/Deep blue | 400–450 | Aerosol retrieval (high sensitivity) |
| Blue | 450–520 | Aerosol retrieval (primary) |
| Green | 520–600 | Surface prior target |
| Red | 600–700 | Surface prior target/query |
| Red edge | 700–750 | Optional vegetation |
| NIR | 750–1000 | Prior query (low aerosol effect) |
| Cirrus/WV | 1350–1420 | Excluded |
| SWIR-1 | 1500–1700 | Prior query (aerosol-transparent) |
| WV absorption | 1800–1950 | Excluded |
| SWIR-2 | 2000–2350 | Prior query (aerosol-transparent) |

### 9.5 Spectral Convolution: Sensor ↔ Reference Band Mapping

Surface priors are stored at a compact **reference basis** (MODIS/VIIRS land bands).
To compare against any input sensor, the mapping uses the reference sensor's **actual RSR**:

$$
\rho_{\text{ref},i} = \frac{\int R_i(\lambda)\, \rho(\lambda)\, d\lambda}{\int R_i(\lambda)\, d\lambda}
$$

Planned API — plain functions, not a class:

```python
def sensor_to_reference(
    sensor_reflectance: xr.Dataset,
    sensor_config: SensorConfig,
    reference_sensor: str = "MODIS",
) -> np.ndarray:
    """Convolve sensor-band reflectance to reference-sensor basis."""
    ...

def reference_to_sensor(
    ref_reflectance: np.ndarray,
    target_bands: list[SpectralBandDescriptor],
    sensor_config: SensorConfig,
    reference_sensor: str = "MODIS",
) -> xr.Dataset:
    """Project reference-basis reflectance back to sensor bands."""
    ...

def load_reference_rsr(reference_sensor: str = "MODIS") -> dict[str, np.ndarray]:
    """Load tabulated spectral response functions for the reference sensor."""
    ...
```

**Tests** (`tests/unit/test_spectral.py`):

| # | Test | Assert |
|---|------|--------|
| 1 | `sensor_to_reference` with identity SRF | output ≈ input (flat spectrum) |
| 2 | `reference_to_sensor` roundtrip | `ref→sensor→ref ≈ original` within tolerance |
| 3 | Output band count matches reference sensor | `len(result) == len(reference_bands)` |
| 4 | `load_reference_rsr("MODIS")` | returns dict with ≥6 bands, arrays non-empty |
| 5 | `load_reference_rsr("UNKNOWN")` | raises `ValueError` |

### 9.6 Surface Prior as an Independent Resource (Preferred)

The surface prior is **decoupled from the runtime pipeline**.
It is pre-computed offline and loaded at runtime, enabling fully parallel Group C.

#### 9.6.1 Prior Store Architecture

```
Surface Prior Store (Zarr)
└── /{tile_id}/
    ├── reflectance   (doy, band, y, x)
    ├── uncertainty   (doy, band, y, x)
    ├── n_obs         (doy, y, x)
    ├── quality       (doy, y, x)
    ├── wavelengths   (band,)
    └── .zattrs
```

Runtime loader:
1. selects tile(s) covering AOI
2. interpolates to observation DOY
3. crops to AOI
4. projects reference bands → sensor bands via `SpectralConvolver`
5. **returns `SurfacePrior`** (the M3 output contract)

#### 9.6.2 Data Sources for Building the Prior

The store can be built from multiple sources:

- **MODIS/VIIRS BRDF products** (global, long record; 500 m)
- **High-resolution composites** (sensor-specific; e.g. Sentinel-2 L2A; 10–20 m)
- **Blended** (recommended): use global BRDF as base and add high-resolution structure where available

**Tests** (`tests/unit/test_prior_store.py`):

| # | Test | Assert |
|---|------|--------|
| 1 | Tile selection covers AOI | selected tiles overlap AOI bounds |
| 2 | DOY interpolation between two snapshots | result DOY between endpoints |
| 3 | AOI crop produces smaller spatial extent | output shape ≤ input shape |
| 4 | Spectral projection (reference → sensor) | output has correct number of sensor bands |
| 5 | Loader returns valid `SurfacePrior` | passes `validate_surface_prior()` |

### 9.7 Optional Runtime Refinement (SWIR/NIR Query)

Optional refinement uses aerosol-insensitive NIR/SWIR bands to query the prior store and update the visible prediction.
This is most useful where climatology is stale or the sensor is hyperspectral.

High-level steps:
1. select query bands by wavelength (NIR+SWIR)
2. apply first-pass correction using prior atmosphere
3. convolve into reference space
4. approximate nearest-neighbour lookup (e.g. FAISS IVF)
5. project predicted visible bands back to sensor space
6. **return `SurfacePrior`** — the output contract is unchanged regardless of method

### 9.8 Comparison of Prior Approaches

| Aspect | Pre-built climatological prior | Runtime SWIR/NIR refinement |
|--------|-------------------------------|----------------------------|
| When computed | Offline | Runtime |
| Depends on TOA bands | No | Yes |
| Pipeline independence | Fully independent (Group C) | Depends on M1 (ObservationBundle) |
| Output contract | `SurfacePrior` | `SurfacePrior` (same) |
| Recommended | Default | Optional accuracy mode |

### 9.9 Planned Architecture + Config Integration

Planned files:

| File | Action | Description |
|------|--------|-------------|
| `python/siac/core/spectral.py` | NEW | `SpectralBandDescriptor`, `sensor_to_reference()`, `reference_to_sensor()`, reference RSR loading |
| `python/siac/core/types.py` | MODIFY | Add `ObservationBundle`, `SolverInputBundle`, `SolvedAtmosphere`; extend/replace `SensorBand` |
| `python/siac/core/validation.py` | NEW | Contract validation functions for all output types |
| `python/siac/pipeline.py` | NEW | `run_pipeline()` function with concurrent module execution |
| `python/siac/siac.py` | MODIFY | `siac_process()` convenience entry point + `_resolve_*` helpers that construct provider instances |
| `python/siac/priors/surface/prior_store.py` | NEW | `PrebuiltPriorStore` class → returns `SurfacePrior` |
| `python/siac/priors/surface/build_prior.py` | NEW | Offline builder tool |
| `python/siac/priors/surface/swir_refine.py` | NEW | Optional runtime refinement → returns `SurfacePrior` |
| `python/siac/grid/assembler.py` | NEW | `assemble_grids()` function → returns `SolverInputBundle` |
| `python/siac/core/config.py` | MODIFY | Surface prior store path + refinement flags |
| `tests/unit/test_spectral.py` | NEW | `SpectralBandDescriptor`, `SensorConfig` selection, spectral convolution tests |
| `tests/unit/test_contracts.py` | NEW | All contract type construction + validation function tests |
| `tests/unit/test_grid_assembler.py` | NEW | `assemble_grids()` unit tests |
| `tests/unit/test_prior_store.py` | NEW | Prior store tile selection, DOY interpolation, spectral projection |
| `tests/integration/test_injection.py` | NEW | Custom provider injection + bad-provider rejection tests |
| `tests/integration/test_orchestration.py` | NEW | Pipeline happy-path, validation-integration, concurrency tests |

---

## 10. Pluggable Data Providers (Atmosphere / BRDF / Surface Prior)

### 10.1 Atmospheric Priors (M2 Providers)

All atmospheric providers must return `AtmosphericState`. How they acquire the data is an implementation detail.

| Provider | Data Source | Notes |
|----------|-----------|-------|
| `CAMSProvider` | Local CAMS NetCDF | Fast, requires local archive |
| `MERRA2Provider` | earthaccess | Remote, NASA auth required |
| `ERA5Provider` | CDS API | Remote, ECMWF auth required |
| `UserAtmoProvider` | User-supplied rasters | Any format; user reads and constructs `AtmosphericState` |
| Custom function | Anything | `(bounds, crs, time, res) -> AtmosphericState` |

**Tests** (`tests/unit/test_atmo_providers.py`):

| # | Test | Assert |
|---|------|--------|
| 1 | Each provider returns valid `AtmosphericState` | passes `validate_atmospheric_state()` |
| 2 | Provider scopes to AOI bounds | spatial extent ≤ requested AOI |
| 3 | Custom function injection | user function called, result validated |
| 4 | Provider with missing auth raises early | `ConfigError` or skip, not silent failure |

### 10.2 Surface Prior (M3 Providers)

All surface prior providers must return `SurfacePrior`. Two internal architectures are supported,
but both produce the same output contract:

1. **BRDF-derived** (legacy / transitional):
   - Fetch `BRDFKernelWeights` from MODIS/VIIRS
   - Apply kernel model with geometry → `SurfacePrior`

2. **Pre-built prior store** (preferred):
   - Load from Zarr store, interpolate to DOY, convolve to sensor bands → `SurfacePrior`

3. **Custom function**:
   - User provides `(bounds, crs, time, sensor_config, geometry, res) -> SurfacePrior`

Core rule: whichever provider is used, it returns `SurfacePrior` with spatially-referenced `boa`, `boa_unc`, and `mask` fields that the Grid Assembler (M4) can resample.

**Tests** (`tests/unit/test_surface_providers.py`):

| # | Test | Assert |
|---|------|--------|
| 1 | BRDF-derived provider returns valid `SurfacePrior` | passes `validate_surface_prior()` |
| 2 | Pre-built store returns valid `SurfacePrior` | passes `validate_surface_prior()` |
| 3 | Custom function injection | user function called, result validated |
| 4 | Output contains `boa`, `boa_unc`, `mask` | all three fields present and non-empty |

### 10.3 RT Model Backend (Cross-Cutting)

Not a data provider — a compute module used by M5 and M6. Implementations:

| Backend | Speed | Jacobian | Sensor coverage |
|---------|-------|----------|-----------------|
| `EmulatorBackend` | Fast | Analytical | S2, L8 (pre-trained) |
| `LUTBackend` | Medium | Numerical | Any sensor (with LUT) |
| `Py6SBackend` | Slow | Numerical | Any sensor |

Users can provide a custom RT backend by implementing the `RTModelBackend` protocol
(this is the one place a multi-method protocol is justified — the backend has
genuinely coupled state: preloaded LUTs, emulator weights, etc.):
```python
class MyRTBackend:
    def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False) -> RTCoefficients: ...
    def supports_jacobian(self) -> bool: ...
    def is_available_for_sensor(self, sensor_id, satellite_id) -> bool: ...
```

**Tests** (`tests/unit/test_rt_backends.py`):

| # | Test | Assert |
|---|------|--------|
| 1 | Each backend satisfies `RTModelBackend` protocol | `isinstance` check passes |
| 2 | `EmulatorBackend` analytical jacobian shape | `(n_params, n_bands)` matches config |
| 3 | `LUTBackend` numerical jacobian finite | no NaN/Inf in output |
| 4 | `compute_coefficients` returns valid `RTCoefficients` | all fields present, physically bounded |
| 5 | Custom backend via user class | user class used, result validated |
| 6 | `is_available_for_sensor` rejects unsupported | returns `False` for unknown sensor |

---

## 11. Centralised Authentication

SIAC v2 accesses multiple remote data sources (CDSE, CAMS/CDS API, AWS S3, NASA Earthdata, GCS). Rather than handling credentials independently in each module, a single **`CredentialManager`** (`siac.core.auth`) acts as a thread-safe credential registry.

### Key types

| Type | Location | Purpose |
|---|---|---|
| `CredentialSpec` | `core/auth.py` | Frozen `(key, secret)` pair |
| `OAuthToken` | `core/auth.py` | Cached bearer token with monotonic expiry |
| `CredentialManager` | `core/auth.py` | Registry + factory + token cache |
| `CredentialConfig` | `core/config.py` | Pydantic model for config-file credentials |
| `AuthenticationError` | `core/exceptions.py` | Raised on missing/failed auth |

### Credential resolution order (in `from_config`)

1. `SIACConfig.credentials` fields (programmatic / YAML)
2. `SIAC_*` environment variables
3. External config files (`~/.cdsapirc`)

### Env-var mapping

| Provider | Key env var | Secret env var |
|---|---|---|
| CDSE | `SIAC_CDSE_USERNAME` | `SIAC_CDSE_PASSWORD` |
| CDS API | `SIAC_CDS_API_KEY` | — |
| AWS/S3 | `SIAC_AWS_ACCESS_KEY_ID` | `SIAC_AWS_SECRET_ACCESS_KEY` |
| Earthdata | `SIAC_EARTHDATA_USERNAME` | `SIAC_EARTHDATA_PASSWORD` |
| GCS | `SIAC_GCS_CREDENTIALS_FILE` | — |

AWS also falls back to standard `AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY`.

### Integration points

- **`SIAC.__init__`** creates `self._auth = CredentialManager.from_config(config)`
- **`siac_process()`** accepts optional `auth` kwarg; defaults to `from_config`
- **`CopernicusDataspaceBackend`** accepts optional `auth` for CDSE OAuth2
- **`CAMSProvider`** accepts optional `auth` for CDS API key injection
- **`_resolve_rt_model_for_pipeline`** injects AWS credentials into `storage_options` for S3 LUT paths

## 12. Earthaccess Rollout Plan (Pre-Implementation)

This section defines the implementation order for NASA Earthdata access via
`earthaccess`, with explicit module boundaries and acceptance criteria.

Detailed execution plan: see `docs/EARTHACCESS_PLAN.md`.

### 12.1 Scope

Planned Earthaccess-backed capabilities:

- Landsat scene discovery + access path in M1 data ingestion flow
- Atmospheric priors using NASA products (including MCD19 AOD path in M2)
- BRDF products for surface priors in M3:
  - MODIS MCD43 family
  - VIIRS VNP43 equivalents

### 12.2 Guardrails

- AOI-scoped queries only (no global over-fetch)
- No direct Earthaccess calls from solver/corrector modules
- Provider outputs must remain contract-pure:
  - M2 returns `AtmosphericState`
  - M3 returns `SurfacePrior`
- All product IDs/short names must be validated in code via dataset discovery
  before hard-coding defaults

### 12.3 Milestones

| Milestone | Goal | Exit Criteria |
|---|---|---|
| M0 | Earthaccess foundation + product registry | Auth/search/open wrappers tested; dataset registry documented |
| M1 | Landsat access path | Scene search/download/open works for AOI/time and feeds M1 preprocessor |
| M2 | Atmospheric priors via Earthaccess | MCD19 AOD path returns valid `AtmosphericState`; fallback strategy defined |
| M3 | BRDF providers | MCD43 + VNP43 providers return valid BRDF inputs to M3 |
| M4 | Pipeline wiring + docs/examples | Config-driven provider resolution works; docs and tests updated |

### 12.4 Definition of Done

- Unit tests for each Earthaccess wrapper/provider module
- Integration tests for at least one live Earthdata query path (network-marked)
- AOI clipping and reprojection behavior verified for M2/M3 outputs
- Config examples for Landsat + MCD19 + MCD43/VNP43 in docs
- No behavior regressions in existing S2/CAMS/local workflows
