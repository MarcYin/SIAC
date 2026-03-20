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
7. [Pipeline Orchestration (Execution Backends + Parallel Stages)](#7-pipeline-orchestration-execution-backends--parallel-stages)
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
- **Decoupled RT and AC choices**: "which radiative-transfer model supplies terms" and "which atmospheric-correction equation consumes them" are separate design axes. For example, emulator-backed coefficient correction and spectral-LUT convolution correction must both fit the architecture.
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

- `siac.algorithms.rt.lut.backend`: RT interpolation logic (`ZarrLUTBackend`)
- `siac.algorithms.rt.lut.store`: local/remote/S3/ZIP store resolution
- `siac.algorithms.rt.lut.http_zip_store`: ReadOnlyZipFileSystem-style ZIP access for local/HTTP/S3 LUT archives
- `siac.algorithms.rt.lut.create`: LUT generation utilities (`create_lut_from_py6s`)
- `siac.algorithms.rt.lut.constants`: default public LUT URL and LUT coordinate constants

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

**Planned BRDF-prior routes**:

1. **Time-local BRDF smoothing route**:
   - Forward-model surface reflectance from BRDF parameters for the target sensor geometry over a
     `±16 d` window around the sensing time.
   - Use the Rust Whittaker smoother from `whitsmooth_rust_repo` on the simulated reflectance time
     series to gap-fill missing days and produce a smooth estimate at the sensing date, with
     uncertainty from residual spread and QA.
   - This is the route closest to the current SIAC paper Appendix A workflow and is the default
     for dynamic per-scene BRDF priors.
   - Spectral mapping is still required in this route whenever the BRDF-source SRF basis differs
     from the target-sensor SRF basis.

2. **Historical monthly best-pixel composite + spectral-mapping route**:
   - Build monthly best-pixel composites of BRDF-simulated reflectance for the scene month and
     adjacent months (`m-1`, `m`, `m+1`) over the previous 5 years, i.e. 15 monthly composites.
   - Pixel selection inside each monthly composite must be driven by BRDF quality so the retained
     BRDF parameters and forward-modelled reflectance come from the best available observation for
     that month, not from a temporal median.
   - Use target-sensor NIR+SWIR observations to query a database built from those 15 monthly
     composites and retrieve the most appropriate visible-reflectance prior for the scene.
   - All MODIS/VIIRS BRDF values remain defined in the MODIS/VIIRS spectral basis; target-sensor
     SWIR/NIR matching must therefore use a spectral-mapping step rather than nearest-band lookup.
   - For Sentinel-2 in particular, SWIR mapping to MODIS/VIIRS must follow the spectral-mapping
     approach described in SIAC Appendix D, not a simple centre-wavelength substitution.

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
| `test_brdf_whittaker_route_uses_pm16d_window` | Mock BRDF time series over ±16 d | Only local window used; smoother output at sensing date |
| `test_brdf_monthly_composite_route_uses_15_candidates` | Scene month `m` over 5 years | Uses `(m-1,m,m+1)` for each year = 15 monthly candidates |
| `test_brdf_monthly_route_uses_swir_nir_query` | Provide target NIR/SWIR observations | Composite ranking/query depends on target SWIR/NIR |
| `test_brdf_monthly_route_requires_spectral_mapping` | Sentinel-2 SWIR + MODIS/VIIRS basis | Uses Appendix-D mapping path, not nearest-band substitution |
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

**Built-in atmospheric-correction family**:

1. **Coefficient-space correction** (current default):
   - Used with RT backends that return `RTCoefficients` directly.
   - Typical sources: neural-network emulators, compact per-band LUTs, Py6S wrappers that already expose `xap/xbp/xcp`, and spectral LUT backends after bandpass convolution.
   - Equation family:
     - `y = xap * toa - xbp`
     - `boa = y / (1 + xcp * y)`

2. **Spectral-LUT derivation path** (backend strategy, not a separate M6 equation):
   - Used when the RT backend starts from spectrally resolved LUT variables (for example `TOA_rho1`, `TOA_rho2`, `Eg_rho1`, `Eg_rho2`, `Eg_dir_rho1`) rather than final per-band coefficients.
   - Workflow:
     - subset LUT over atmospheric state,
     - convolve LUT wavelength slices with the sensor bandpass / SRF in radiance and irradiance space,
     - derive per-band `path_ref`, `T_total`, and `S` (spherical-albedo term),
     - convert them into standard coefficients:
       - `xap = 1 / T_total`
       - `xbp = path_ref / T_total`
       - `xcp = S`
   - Once these coefficients are derived, M6 reuses the same coefficient-space equation as every other built-in path.
   - This is the preferred backend design for hyperspectral sensors or any sensor whose band centers/FWHM vary enough that pre-baked coefficients are not the right abstraction.

3. **Custom correction**:
   - User supplies `CorrectorFn` directly.
   - This covers BRDF-aware variants, adjacency-aware variants, or mission-specific radiance-to-reflectance schemes that do not fit either built-in family.

**Dispatch rule**:
- `CorrectorFn` stays the public injection boundary.
- The built-in `correct_atmosphere()` remains a coefficient-space corrector.
- Spectral LUT logic belongs in the RT backend, which derives `RTCoefficients` before handing values to M6.
- The solver contract does not change.

**Tests** (`tests/unit/test_correction.py` — extend existing):

| Test | What | Assert |
|------|------|--------|
| `test_correct_output_type` | Call with `ObservationBundle` + `SolvedAtmosphere` | Returns `CorrectionResult` |
| `test_correct_boa_bands_match_toa` | Check BOA band names | Same as input TOA |
| `test_correct_boa_range` | Check BOA values | In [–0.05, 1.5] (minor negatives allowed) |
| `test_correct_cloud_mask_preserved` | Check `result.cloud_mask` | Matches original from M1 |
| `test_correct_native_resolution` | Check BOA spatial shape | Matches M1 native TOA shape |
| `test_correct_metadata_timing` | Check `result.diagnostics` | Contains `processing_time_s` |
| `test_correct_coefficient_mode` | Emulator-style backend | Uses `RTCoefficients` path |
| `test_correct_spectral_lut_coefficients` | Dense spectral LUT backend | Derived coefficients match direct `path_ref/T_total/S` formulation |

### 4.7 RT Model Backend (Cross-Cutting)

The RT model is a **compute-oriented service** used by both M5 and M6. A backend may lazily initialise from local files or remote stores (for example a remote ZIP-hosted Zarr LUT), but once opened it behaves like a local service from the pipeline's perspective.

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

Dense spectral LUT workflows still fit this protocol. Their internal implementation may start from wavelength-resolved terms, but after SRF/bandpass convolution they should convert to standard coefficients before returning:

- `xap = 1 / T_total`
- `xbp = path_ref / T_total`
- `xcp = S`

Design intent:
- `RTModelBackend` remains the main public contract.
- A `SpectralLUTBackend` may do substantially more internal work than an emulator backend, but it still returns `RTCoefficients` at the boundary.
- This keeps M6 and the user injection surface stable.

**Tests** (`tests/unit/test_rt_backend.py`):

| Test | What | Assert |
|------|------|--------|
| `test_rt_backend_protocol_conformance` | `MockRTModel` satisfies `RTModelBackend` | Duck-typing check passes |
| `test_rt_backend_jacobian_shape` | `compute_coefficients(..., compute_jacobian=True)` | `d_xap` has `param` dim with `["aot", "tcwv"]` |
| `test_rt_backend_no_jacobian` | `compute_coefficients(..., compute_jacobian=False)` | `d_xap/d_xbp/d_xcp` are `None` |
| `test_rt_backend_unsupported_sensor` | `is_available_for_sensor("UNKNOWN", "SAT")` | Returns `False` |
| `test_rt_coefficients_apply_correction` | Call `apply_correction(toa)` | BOA = known analytical formula |
| `test_spectral_rt_backend_terms_shape` | `prepare_band_terms(...)` | Returns aligned `(path_ref, T_total, S)` arrays |

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
from siac.runtime import ObservationBundle, AtmosphericState, SurfacePrior

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
    config=SIACConfig.from_file("config.toml"),
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

## 7. Pipeline Orchestration (Execution Backends + Parallel Stages)

### 7.1 Sensor-Agnostic Pipeline Flow

Expressed in terms of modules and contracts (target design):

```
1. M1a: metadata + TOA load                                           → partial ObservationBundle inputs
2. M1b: geometry + cloud/classes from shared TOA (parallel branches)  → complete ObservationBundle
3. M2: get_prior(bounds, crs, time, res)                              → AtmosphericState
4. M3: get_surface_prior(bounds, crs, time, sensor_config, geometry, res) → SurfacePrior
5. M4: assemble(obs, atmo, surface, rt_model)                         → SolverInputBundle
6. M5: solve(inputs, solver_config)                                   → SolvedAtmosphere
7. M6: correct(obs, solved, rt_model)                                 → CorrectionResult
```

Execution policy:

- M1b geometry and cloud tasks run concurrently and share TOA to avoid duplicate reads.
- M2 and M3 run concurrently once required metadata (`bounds`, `crs`, `time`) is available.
- M4–M6 remain ordered.

### 7.2 Execution Backends

Two orchestration backends are supported:

| Backend | Role | Notes |
|------|------|--------|
| `thread` | Compatibility fallback | Local `ThreadPoolExecutor`; minimal dependencies. |
| `dask` | Primary target backend | `dask.distributed` scheduling, retries, diagnostics dashboard, performance reports. |

Design rule:

- Module contracts do not change by backend.
- Backend selection only changes scheduling/orchestration behavior.

### 7.3 Parallel Fetch Groups

Parallel groups mapped to module responsibilities:

- **Group A (M1b-geom)**: geometry extraction from M1 inputs.
- **Group B (M1b-cloud)**: cloud/class generation from shared TOA and metadata.
- **Group C (M2)**: atmospheric prior retrieval.
- **Group D (M3)**: surface prior retrieval.

Group A and B are always independent after TOA is loaded. Group C and D are independent after M1 metadata is available.

### 7.4 Input Requirements Per Processing Stage

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

### 7.5 Dependency Analysis

```
  TIME ────────────────────────────────────────────────────────────────────▶

  M1a ▐████████████████▌ metadata + TOA load
            │
            ├──────────► M1b-geom  ▐████▌
            └──────────► M1b-cloud ▐█████▌

  M2  (after M1a metadata)  ▐████████████▌
  M3  (after M1a metadata)  ▐██████████████▌
                    │
    ┌───────────────▼─────────────────┐
    │ M4: assemble → SolverInputBundle│  ← starts at max(M1b, M2, M3)
    └───────────────┬─────────────────┘
                    │
    ┌───────────────▼─────────────────┐
    │ M5: solve → SolvedAtmosphere    │
    └───────────────┬─────────────────┘
                    │
    ┌───────────────▼─────────────────┐
    │ M6: correct → CorrectionResult  │
    └─────────────────────────────────┘
```

### 7.6 Implementation: `run_pipeline()` Backend Strategy

`python/siac/pipeline.py` should dispatch by execution backend:

```python
def run_pipeline(..., execution: ExecutionConfig) -> CorrectionResult:
    if execution.backend == "dask":
        return _run_pipeline_dask(..., execution=execution)
    return _run_pipeline_thread(..., execution=execution)
```

Dask path design:

- Create a `dask.distributed.Client`.
- Submit M1 split tasks, then dependent M2/M3 tasks.
- Use retries and per-stage timeout configuration for remote I/O tasks.
- Optionally generate `performance_report` HTML.
- Preserve existing validation checks before M4.

Thread path design:

- Maintain simple local fallback with `ThreadPoolExecutor`.
- Keep behavior parity with Dask path (same contracts, same exceptions).

Operational diagnostics:

- Dask dashboard URL emitted in logs when enabled.
- Optional run artifact: `output_dir/reports/dask-performance.html`.
- Optional stage summary JSON for CI artifacts.

### 7.7 Multi-Resolution Data Assembly (M4 Responsibility)

SIAC operates across three spatial scales. The Grid Assembler (M4) is solely
responsible for aligning data between them:

- **Native grid**: sensor-native per-band resolution (input to/output from M6)
- **Aux grid**: common analysis grid (e.g. 500 m) for solver inputs
- **Aerosol grid**: coarser retrieval grid (e.g. 1000 m) for AOT/TCWV

Core requirement: M4 produces a `SolverInputBundle` where all raster data is
on the same grid. M6 receives native-resolution data from M1 and upsamples solved
fields.

### 7.8 Verification Plan (Orchestration)

**Tests** (`tests/integration/test_pipeline.py` and `tests/unit/test_satellite.py` — extend existing):

**`run_pipeline()` — Happy path & ordering:**

| Test | What | Assert |
|------|------|--------|
| `test_pipeline_happy_path` | All modules mocked | Returns `CorrectionResult` |
| `test_pipeline_call_order` | Trace stage transitions | M1 before M4; M4 before M5; M5 before M6 |
| `test_m1_geometry_cloud_parallel` | M1 geometry/cloud each sleep 0.5 s | Combined latency < 0.8 s |
| `test_shared_toa_single_read` | Instrument TOA read counter | TOA read exactly once |

**Backend parity tests:**

| Test | What | Assert |
|------|------|--------|
| `test_dask_vs_thread_equivalence` | Same mocked inputs/backends | Output arrays equal within tolerance |
| `test_thread_backend_fallback` | `backend="thread"` | Pipeline still works without Dask |

**Resilience tests (Dask path):**

| Test | What | Assert |
|------|------|--------|
| `test_retry_transient_m2_failure` | M2 fails once then succeeds | Pipeline succeeds with retry |
| `test_stage_timeout` | Slow M3 exceeds timeout | Clear timeout error with stage context |
| `test_cancellation_propagation` | Hard failure in one branch | Dependent futures canceled/aborted |

**Observability tests:**

| Test | What | Assert |
|------|------|--------|
| `test_performance_report_output` | `performance_report` enabled | HTML file created |
| `test_progress_summary_output` | stage-summary enabled | JSON summary created |

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
    def extract_cloud_mask(self, input_path: Path, toa: xr.Dataset | None = None) -> xr.DataArray: ...
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

#### 9.6.3 BRDF Prior Usage: Two Operational Routes

The BRDF prior should be treated as two distinct operational routes, both of
which return the same `SurfacePrior` contract.

**Route A: time-local Whittaker-smoothed BRDF prior**

- Inputs:
  - BRDF kernel parameters (`f_iso`, `f_vol`, `f_geo`) from `MCD43A1`, `VNP43MA1`, or equivalent
  - target-scene illumination/view geometry
  - valid samples in a `±16 d` window around the sensing date
- Process:
  1. forward-simulate BRDF reflectance for each available day in the `±16 d` window using the
     target-scene angles
  2. mask poor-quality days/pixels with product QA
  3. apply the Rust Whittaker smoother per band/pixel to gap-fill and smooth the time series
  4. if source and target SRFs differ, run spectral mapping from the BRDF source basis to the
     target-sensor basis
  5. evaluate the smoothed signal at the sensing date to obtain the BRDF prior mean
  6. derive uncertainty from QA, residual spread, smoother diagnostics, and spectral-mapping error
- Role:
  - default dynamic BRDF prior for atmospheric retrieval
  - closest to the current SIAC observational-prior construction

**Route B: historical monthly best-pixel composite + SWIR/NIR query**

- Inputs:
  - monthly best-pixel composites of BRDF-simulated reflectance
  - months `(m-1, m, m+1)` for each of the previous 5 years (15 monthly composites total)
  - target-scene NIR+SWIR reflectance from the target sensor
- Process:
  1. precompute monthly best-pixel composites in the MODIS/VIIRS BRDF spectral basis, choosing the
     monthly pixel from the best BRDF quality rather than by median compositing
  2. for a target scene, gather the 15 monthly composites for `(m-1, m, m+1)` across the previous
     5 years
  3. if source and target SRFs differ, use spectral mapping to compare target-scene NIR+SWIR
     observations with the MODIS/VIIRS BRDF basis
  4. build or load a database whose query key contains:
     - NIR
     - SWIR-1
     - SWIR-2
     - the median summary derived from the 15 monthly composites
  5. apply a first-pass atmospheric correction with prior atmospheric parameters to the target
     scene, then use the corrected NIR/SWIR query vector to retrieve the visible reflectance
     estimate from the database
  6. return visible/NIR/SWIR prior reflectance and uncertainty for the retrieved candidate, and
     feed the visible prior into the AOD solve
- Role:
  - climatological/historical fallback route when the local time series is weak, cloudy, or sparse
  - route for stable prior stores built offline as a searchable 15-month BRDF-composite database

The planner should treat Route A and Route B as first-class provider strategies,
not as minor implementation options inside one opaque class.

**Tests** (`tests/unit/test_prior_store.py`):

| # | Test | Assert |
|---|------|--------|
| 1 | Tile selection covers AOI | selected tiles overlap AOI bounds |
| 2 | DOY interpolation between two snapshots | result DOY between endpoints |
| 3 | AOI crop produces smaller spatial extent | output shape ≤ input shape |
| 4 | Spectral projection (reference → sensor) | output has correct number of sensor bands |
| 5 | Loader returns valid `SurfacePrior` | passes `validate_surface_prior()` |
| 6 | Monthly composite builder produces 15 best-pixel composites | `(m-1,m,m+1)` over 5 years present and no temporal-median composite is used |
| 7 | Monthly composite builder respects BRDF quality | best-quality BRDF pixel is selected for each month |
| 8 | Monthly database query responds to NIR + two SWIR bands | retrieved visible prior changes with the corrected query bands |
| 9 | Monthly database feature summary is stable | median summary over the 15 monthly composites is computed consistently |

### 9.7 Spectral Mapping for BRDF Composite Usage

Spectral mapping is required for **both Route A and Route B** whenever the BRDF
source basis and the target-sensor SRFs are not the same. Route B simply makes
that requirement more obvious because the query is performed in NIR/SWIR.

The required rule is:

> **Do not map Sentinel-2/Landsat SWIR bands to MODIS/VIIRS BRDF bands by nearest wavelength alone.**
> Use the Appendix D SIAC spectral-mapping workflow whenever source and target SRFs differ.

Per SIAC Appendix D, the planned runtime/offline mapping should:

1. use the MODIS/VIIRS reflectance basis and MODIS/VIIRS SRFs as the reference basis
2. search a hyperspectral reflectance library for spectra consistent with the reference-basis reflectance
3. reconstruct a 1 nm reflectance estimate from the selected neighbours
4. convolve the reconstructed spectrum with both:
   - MODIS/VIIRS SRFs
   - target-sensor SRFs
5. carry the mapping uncertainty from the neighbour dispersion / reconstruction error

Explicit exception:

- **Hyperspectral → multispectral** projection does **not** require an external hyperspectral
  library/database when the source sensor already measures the spectrum densely enough in wavelength.
  In that case, multispectral simulation can be done directly by convolving the hyperspectral
  reflectance with the target multispectral SRFs.
- The external-library / Appendix-D reconstruction route is required when translating between
  multispectral bases, e.g. MODIS/VIIRS BRDF basis to Sentinel-2 or Landsat.

This is especially important for:
- Sentinel-2 SWIR bands (`B11`, `B12`)
- Landsat OLI SWIR bands
- any target sensor whose NIR/SWIR SRFs differ materially from MODIS/VIIRS band passes

The spectral-mapping output is not only a visible-band helper. It is also the
bridge that lets target-sensor SWIR/NIR observations query a MODIS/VIIRS BRDF
composite database consistently, including the Route-B database keyed by NIR,
two SWIR bands, and the median summary from the 15 monthly composites.

### 9.8 Optional Runtime Refinement (SWIR/NIR Query)

Optional refinement uses aerosol-insensitive NIR/SWIR bands to query the prior store and update the visible prediction.
In the BRDF-composite route, this query must happen in the MODIS/VIIRS reference basis via the Appendix-D spectral mapping.
This is most useful where climatology is stale or the sensor is hyperspectral.

High-level steps:
1. select query bands by wavelength (`NIR`, `SWIR-1`, `SWIR-2`)
2. apply first-pass correction using prior atmosphere
3. if source and target SRFs differ, map target reflectance into the BRDF source basis using
   Appendix-D spectral mapping; if the source is hyperspectral and the target is multispectral,
   use direct SRF convolution instead
4. form the Route-B database query key from:
   - corrected `NIR`
   - corrected `SWIR-1`
   - corrected `SWIR-2`
   - the median summary derived from the 15 monthly composites
5. query/select the monthly composite or prior-store candidate
6. use the retrieved candidate to estimate visible surface reflectance per pixel for the AOD solve
7. project predicted visible bands back to sensor space
8. **return `SurfacePrior`** — the output contract is unchanged regardless of method

### 9.9 Comparison of Prior Approaches

| Aspect | Pre-built climatological prior | Runtime SWIR/NIR refinement |
|--------|-------------------------------|----------------------------|
| When computed | Offline | Runtime |
| Depends on TOA bands | No | Yes |
| Pipeline independence | Fully independent (Group C) | Depends on M1 (ObservationBundle) |
| Output contract | `SurfacePrior` | `SurfacePrior` (same) |
| Recommended | Default | Optional accuracy mode |

### 9.10 Planned Architecture + Config Integration

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
| `python/siac/priors/surface/brdf_whittaker.py` | NEW | Route A: `±16 d` BRDF prior using `whitsmooth_rust_repo` |
| `python/siac/priors/surface/brdf_monthly_composite.py` | NEW | Route B: monthly best-pixel BRDF composites over 5-year history |
| `python/siac/priors/surface/brdf_monthly_database.py` | NEW | Route B database builder/query for 15 monthly composites keyed by NIR + two SWIR bands + median summary |
| `python/siac/priors/surface/swir_refine.py` | NEW | Optional runtime refinement → returns `SurfacePrior` |
| `python/siac/priors/surface/spectral_mapping.py` | NEW | SRF-dependent mapping layer: Appendix-D reconstruction for multispectral↔multispectral, direct convolution for hyperspectral→multispectral |
| `python/siac/grid/assembler.py` | NEW | `assemble_grids()` function → returns `SolverInputBundle` |
| `python/siac/core/config.py` | MODIFY | Surface prior store path + refinement flags |
| `tests/unit/test_spectral.py` | NEW | `SpectralBandDescriptor`, `SensorConfig` selection, spectral convolution tests |
| `tests/unit/test_contracts.py` | NEW | All contract type construction + validation function tests |
| `tests/unit/test_grid_assembler.py` | NEW | `assemble_grids()` unit tests |
| `tests/unit/test_prior_store.py` | NEW | Prior store tile selection, DOY interpolation, spectral projection |
| `tests/unit/test_brdf_whittaker.py` | NEW | Route A Whittaker gap-filling and sensing-date evaluation via Rust smoother |
| `tests/unit/test_brdf_monthly_composite.py` | NEW | Route B 15-month best-pixel composite logic and BRDF-quality selection |
| `tests/unit/test_brdf_monthly_database.py` | NEW | Route B database build/query with NIR + two SWIR bands + median-summary key |
| `tests/unit/test_spectral_mapping.py` | NEW | Appendix-D mapping for differing multispectral SRFs and direct hyperspectral→multispectral convolution |
| `tests/integration/test_injection.py` | NEW | Custom provider injection + bad-provider rejection tests |
| `tests/integration/test_orchestration.py` | NEW | Pipeline happy-path, validation-integration, concurrency tests |

### 9.10 Spectral Response Function (SRF) Rollout Plan

The current Gaussian centre/FWHM approximation is not sufficient once SIAC
needs to distinguish between satellite platforms such as `S2A`, `S2B`, `S2C`,
`L8`, and `L9`. The runtime needs a proper SRF plan with three separate layers:

1. **Source access**: how SIAC locates authoritative SRF publications
2. **Canonicalization**: how different vendor formats are converted into one SIAC format
3. **Runtime consumption**: how LUT, emulator, and prior code use the canonical SRF

The design rule is:

> **If an official SRF exists for a platform, SIAC should use that SRF as the primary spectral definition.**
> Gaussian centre/FWHM is a fallback only for sensors without published SRFs.

Package boundary rule:
- `python/siac/core/` should only hold cross-cutting contracts and config
- canonical SRF types should live under `python/siac/domain/spectral.py`
- external sensor response loading should go through the `RSRF` adapter in `python/siac/adapters/rsrf.py`
- LUT-specific aligned-kernel logic should live under `python/siac/rt/lut/`

#### 9.10.1 SRF Source Access Layer

SIAC should not fetch ad hoc SRF files from arbitrary URLs during runtime.
Instead, it should have a small source registry that points to authoritative
sensor-specific SRF publications.

Implemented architecture direction:
- `python/siac/adapters/rsrf.py` is the generic SRF loading entry point
- authoritative sampled curves and band specs come from the external `RSRF` package
- SIAC should not maintain its own parallel SRF source registry or remote-download layer
  - `load_sensor_config_from_local_srf_file(...)`
- metadata-driven sensors should build a `SensorConfig` from per-band
  centre-wavelength / FWHM metadata through a shared builder, rather than
  pretending they all have one mission-wide tabulated RSR workbook

Planned source families:

| Family | Typical source | Platforms |
|--------|----------------|-----------|
| ESA / Copernicus | Sentinel-2 SRF tables | `S2A`, `S2B`, `S2C` |
| USGS / NASA | Landsat RSR tables | `L8`, `L9` |
| NASA / NOAA | MODIS / VIIRS reference RSR | prior reference basis |
| User-local | CSV / TSV / NetCDF provided by user | custom sensors |

Sentinel-2 source note:
- the authoritative landing page is [SentiWiki S2 Mission](https://sentiwiki.copernicus.eu/web/s2-mission)
- the Sentinel-2 parser should resolve the linked `Sentinel-2 Spectral Response Functions (S2-SRF)` document from that page / its linked documents page, rather than pinning a brittle attachment URL in code

Cross-sensor source inventory currently identified:

| Sensor / platform | Access pattern | Spectral definition | Official source |
|------------------|----------------|---------------------|-----------------|
| `MSI / S2A` | remote official file | tabulated RSR | [SentiWiki S2 Mission](https://sentiwiki.copernicus.eu/web/s2-mission) -> [S2 documents SRF section](https://sentiwiki.copernicus.eu/web/s2-documents?inheritRedirect=true#S2Documents-SPECTRALRESPONSEFUNCTIONS) |
| `MSI / S2B` | remote official file | tabulated RSR | same as `S2A` |
| `MSI / S2C` | remote official file | tabulated RSR | same as `S2A` |
| `OLI / L8` | remote catalog / export | tabulated RSR | [USGS Landsat Spectral Characteristics Viewer](https://landsat.usgs.gov/spectral-characteristics-viewer) |
| `OLI-2 / L9` | remote catalog / export | tabulated RSR | [USGS Landsat Spectral Characteristics Viewer](https://landsat.usgs.gov/spectral-characteristics-viewer) |
| `MODIS / Terra` | remote official table | tabulated RSR | [NASA Ocean Color RSR tables](https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables/) |
| `MODIS / Aqua` | remote official table | tabulated RSR | [NASA Ocean Color RSR tables](https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables/) |
| `VIIRS / SNPP` | remote official table | tabulated RSR | [NASA Ocean Color RSR tables](https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables/) |
| `VIIRS / NOAA-20` | remote official monitoring page | tabulated / ancillary files | [NOAA STAR ICVS VIIRS N20](https://www.star.nesdis.noaa.gov/icvs/status_N20_VIIRS.php) |
| `VIIRS / NOAA-21` | remote official monitoring page | tabulated / ancillary files | [NOAA STAR ICVS VIIRS N21](https://www.star.nesdis.noaa.gov/icvs/status_N21_VIIRS.php) |
| `OLCI / S3A` | remote official document set | tabulated RSR | SentiWiki Sentinel-3 OLCI instrument performance documents |
| `OLCI / S3B` | remote official document set | tabulated RSR | SentiWiki Sentinel-3 OLCI instrument performance documents |
| `SLSTR / S3A` | remote official document set | tabulated RSR | SentiWiki Sentinel-3 SLSTR instrument performance documents |
| `PRISMA` | scene product metadata | metadata band characterization | [ASI PRISMA mission](https://www.asi.it/en/earth-science/prisma/) |
| `EnMAP` | scene product metadata | metadata band characterization | [EnMAP L1/L2 product specification](https://www.enmap.org/data/doc/EN-PCV-ICD-2009-2_HSI_Product_Specification_Level1_Level2.pdf) |
| `EMIT` | scene product metadata | metadata band characterization | [EMIT L2A reflectance ATBD](https://lpdaac.usgs.gov/documents/2147/EMIT_L2A-RFL_ATBD_V1.pdf) |
| `PlanetScope / PS2` | remote official file | tabulated RSR | [Planet RSR access article](https://support.planet.com/hc/en-us/articles/4411132050451-How-Can-I-Access-Relative-Spectral-Responses-RSRs-) |
| `PlanetScope / PS2.SD` | remote official file | tabulated RSR | [Planet RSR access article](https://support.planet.com/hc/en-us/articles/4411132050451-How-Can-I-Access-Relative-Spectral-Responses-RSRs-) |
| `PlanetScope / PSB.SD` | remote official file | tabulated RSR | [Planet RSR access article](https://support.planet.com/hc/en-us/articles/4411132050451-How-Can-I-Access-Relative-Spectral-Responses-RSRs-) |

Design rule from that inventory:
- remote tabulated SRFs and metadata-derived band characterizations are not the
  same class of source and must not share one parser path
- a source registry should tell the runtime whether to fetch a remote file,
  parse product metadata, or expect a user-supplied local file
- hyperspectral missions such as `PRISMA`, `EnMAP`, and `EMIT` should be
  treated as metadata-driven first, unless a stable official tabulated RSR file
  is verified later

Planned source manifest contract:

```python
@dataclass(frozen=True)
class SRFSourceSpec:
    source_id: str
    sensor_id: str
    satellite_id: str
    version: str
    source_url: str | None
    local_path: Path | None
    file_format: str          # "csv", "xlsx", "txt", "nc", ...
    parser_name: str
    checksum_sha256: str | None
    licence: str | None
```

Rules:
- only official or user-explicit sources are allowed
- the manifest must pin version and checksum where possible
- remote access is a build/update task, not a runtime dependency
- runtime code reads from a local SRF repository generated from the manifest

#### 9.10.2 Canonical SIAC SRF Format

Raw SRF publications come in different file layouts, wavelength units, and band
names. SIAC needs one canonical in-memory format for all sensors.

Planned canonical object:

```python
@dataclass(frozen=True)
class SpectralResponseFunction:
    sensor_id: str
    satellite_id: str
    band_name: str
    wavelengths_nm: np.ndarray
    response: np.ndarray
    response_raw: np.ndarray | None = None
    source_id: str | None = None
    source_version: str | None = None
    source_url: str | None = None
    centre_wavelength_nm: float | None = None
    effective_wavelength_nm: float | None = None
    fwhm_nm: float | None = None
```

Canonicalization rules:
- `wavelengths_nm` is strictly ascending and stored in nanometres
- `response` is dimensionless, non-negative, finite, and **area-normalized**
- `response_raw` optionally preserves the published relative-response values
- `centre_wavelength_nm`, `effective_wavelength_nm`, and `fwhm_nm` are derived diagnostics, not the authoritative definition
- band identity is always `(sensor_id, satellite_id, band_name)`, not just `band_name`

The runtime should integrate with `response`, not with `response_raw`.
Area-normalized response is the correct form for spectral convolution:

$$
\int R(\lambda)\,d\lambda = 1
$$

This avoids ambiguity around peak-normalized curves from vendor files.

#### 9.10.3 Raw Source Conversion Pipeline

Each raw SRF family should have a dedicated converter that maps its native file
format into `SpectralResponseFunction`.

Planned conversion stages:

1. load vendor file with a source-specific parser
2. map vendor band labels to SIAC canonical band names
3. convert wavelength units to nanometres
4. sort wavelengths and drop duplicates
5. clip tiny negative numerical artefacts to zero
6. preserve the published response as `response_raw`
7. compute normalized `response`
8. derive `effective_wavelength_nm`, `centre_wavelength_nm`, and `fwhm_nm`
9. validate and store in the local SRF repository

Planned converters:

| Converter | Input families | Notes |
|----------|----------------|-------|
| `parse_esa_s2_srf(...)` | Sentinel-2 ESA SRF tables | source is the SentiWiki `S2 Mission` page and its linked `S2-SRF` document; distinguish `S2A`, `S2B`, `S2C` explicitly |
| `parse_usgs_landsat_rsr(...)` | Landsat 8/9 RSR files | keep `L8` and `L9` separate |
| `parse_reference_rsr(...)` | MODIS / VIIRS RSR tables | used by surface-prior projection |
| `parse_user_srf(...)` | user CSV / TSV / NetCDF | requires explicit metadata mapping |

#### 9.10.4 Local SRF Repository

Runtime code should not know about raw vendor files. It should resolve SRFs
from a local SIAC repository.

Planned repository responsibilities:
- load SRF by `(sensor_id, satellite_id, band_name)`
- expose all bands for a platform
- return provenance metadata with each SRF
- cache interpolation of SRFs onto the active LUT wavelength grid

Planned API:

```python
class SRFRepository:
    def get_band_srf(
        self,
        sensor_id: str,
        satellite_id: str,
        band_name: str,
    ) -> SpectralResponseFunction: ...

    def get_sensor_srfs(
        self,
        sensor_id: str,
        satellite_id: str,
    ) -> dict[str, SpectralResponseFunction]: ...
```

Storage choice:
- authoritative local store should be versioned and testable
- actual on-disk format can be `zarr`, `netcdf`, or `npz`
- the repository API is the stable boundary; storage format is not

#### 9.10.5 SRF Storage Policy for LUT Usage

The SRF is used by the dense spectral LUT backend, so the storage policy must
be designed around LUT convolution rather than around plotting convenience.

The decision is:

> **Do not store the authoritative SRF on the LUT wavelength grid.**
> Store the canonical SRF in a sensor-native form, then derive a LUT-aligned
> kernel for the active LUT grid and cache that derived kernel.

This avoids coupling SRF content to one LUT version or wavelength spacing.
Different LUT products may use different wavelength axes, so storing the
authoritative SRF on the LUT grid would duplicate data and make SRFs depend on
the current LUT implementation.

##### Canonical SRF storage

Canonical SRFs should be stored only over their real spectral support:

- start at the first wavelength where the published response becomes non-zero
- end at the last wavelength where the published response is non-zero
- preserve one explicit zero-valued boundary sample on each side if the source
  format provides it, or synthesize equivalent boundary zeros during conversion
- do **not** pad the canonical SRF over the full LUT wavelength domain

This means the canonical object is compact and physically meaningful. The
authoritative SRF is the band support, not a large mostly-zero vector.

##### Derived LUT-aligned kernel

For actual convolution, the runtime should build a derived SRF kernel on the
active LUT wavelength axis:

```python
@dataclass(frozen=True)
class AlignedSRFKernel:
    sensor_id: str
    satellite_id: str
    band_name: str
    lut_id: str
    wavelength_axis_hash: str
    start_index: int
    end_index: int
    wavelengths_nm: np.ndarray
    response_on_lut: np.ndarray
    solar_weighted_response_on_lut: np.ndarray | None = None
```

Design rules:
- `wavelengths_nm` must be a slice of the LUT wavelength coordinate
- `response_on_lut` is the canonical SRF interpolated onto the LUT axis
- only the support slice should be stored, not the full zero-padded LUT axis
- `start_index:end_index` maps the support slice back into the parent LUT axis
- `solar_weighted_response_on_lut` is optional and only valid for LUT terms that
  require irradiance-weighted convolution

This gives the backend a structured SRF representation without making the SRF
repository depend on LUT layout.

##### Wavelength spacing policy

The canonical SRF should keep the source publication's effective sampling after
cleanup. SIAC should **not** force all sensors onto a universal SRF spacing such
as 1 nm, and should **not** force canonical SRFs onto the LUT spacing.

Instead:
- the LUT remains authoritative for convolution spacing during runtime
- the SRF remains authoritative for the band shape in the repository
- interpolation bridges the two only when a specific LUT is active

If a vendor SRF is extremely coarse or irregular, the converter may densify it
as a controlled preprocessing step, but that is a source-family parser detail,
not a global storage rule.

##### Support window for LUT subsetting

LUT wavelength subsetting for a band should be driven by the tabulated SRF
support, not by Gaussian sigma rules.

Planned rule:
- find the SRF support bounds from the canonical SRF
- map them onto the LUT wavelength axis
- expand the selected LUT slice by one wavelength sample on each side when
  available, so numerical integration at the support edge remains stable
- interpolate the SRF only on that slice

This replaces the current "centre wavelength ± sigma window" approximation for
platforms that have a real SRF.

##### Persisted artifacts

Only the canonical SRF repository is authoritative.
Derived LUT-aligned kernels are cache artifacts.

Allowed persistence options for derived kernels:
- in-memory cache keyed by `(lut_id, wavelength_axis_hash, sensor_id, satellite_id, band_name)`
- optional on-disk cache for expensive public LUTs

If on-disk caching is added, it must be versioned by LUT identity. A cached
`S2A B08` kernel for one LUT wavelength axis must not be reused for another LUT
with different spacing or bounds.

#### 9.10.6 Runtime Integration Rules

`SensorConfig` and `ObservationBundle` must carry the actual platform identity
used to resolve SRFs. This means Sentinel-2 processing must stop collapsing
unknown platforms to `S2A`.

Required runtime changes:
- `SensorConfig` should reference canonical SRFs by platform-specific key
- `Sentinel2Preprocessor` must resolve `S2A`, `S2B`, and `S2C` explicitly
- Landsat paths must keep `L8` and `L9` distinct
- `SpectralBandDescriptor` should treat tabulated SRF as first-class, not optional decoration
- LUT spectral convolution must use the tabulated SRF when available
- emulator and prior-projection code should use derived spectral metadata from the SRF repository, not hand-written constants

Interpolation policy:
- canonical SRFs remain in their native tabulated form
- at runtime they are interpolated onto the active LUT wavelength axis
- interpolation is cached per `(lut_id, platform, band, wavelength_axis)`
- Gaussian fallback is only allowed when the repository has no SRF for that platform

#### 9.10.7 Expected SRF Design by Sensor Family

Platform specificity must be represented at the right level:

| Sensor family | Required SRF identity | Notes |
|--------------|------------------------|-------|
| Sentinel-2 MSI | `(MSI, S2A, band)`, `(MSI, S2B, band)`, `(MSI, S2C, band)` | platform-specific SRFs are required |
| Landsat OLI | `(OLI, L8, band)`, `(OLI, L9, band)` | platform-specific RSRs are required |
| MODIS / VIIRS reference | reference-sensor band id | used for prior-space mapping |
| Custom sensors | `(sensor_id, satellite_id, band)` | supplied by user manifest |

Detector-level SRFs are out of scope for the first implementation.
The first stable unit is **platform-level band SRF**. If detector-specific SRFs
are needed later, they should extend the key with an optional detector id
without changing the platform-level API.

#### 9.10.8 Validation and Test Plan

SRF handling needs both parser tests and physics tests.

Required tests:

| Layer | Test | Assert |
|------|------|--------|
| unit | source manifest validation | bad parser / checksum config is rejected |
| unit | raw parser per source family | correct band count and provenance extracted |
| unit | normalization | `response >= 0`, finite, area integrates to `≈ 1` |
| unit | wavelength cleanup | sorted ascending, duplicates removed |
| unit | derived diagnostics | effective wavelength and FWHM are finite |
| unit | repository lookup | correct platform-specific SRF returned |
| unit | LUT-kernel build | support slice and aligned weights are correct |
| integration | Sentinel-2 platform split | `S2A`, `S2B`, `S2C` resolve different SRFs |
| integration | LUT convolution | tabulated SRF path is used when SRF exists |
| regression | fallback control | Gaussian path only used for sensors without SRF |

#### 9.10.9 Planned Files

| File | Action | Description |
|------|--------|-------------|
| `python/siac/domain/spectral.py` | NEW | `SpectralResponseFunction` dataclass + validation helpers |
| `python/siac/adapters/rsrf.py` | NEW | Thin adapter over the external `RSRF` package |
| `python/siac/rt/lut/srf_kernel.py` | NEW | LUT-aligned SRF kernel builder and cache-key helpers |
| external `RSRF` data root | REQUIRED | authoritative sampled curves and band specifications |
| `tools/build_srf_repository.py` | NEW | build/update tool from official raw SRF sources |
| `python/siac/core/spectral.py` | MODIFY | consume `SpectralResponseFunction` directly |
| `python/siac/core/types.py` | MODIFY | make band identity platform-specific and SRF-aware |
| `python/siac/satellite/sentinel2.py` | MODIFY | resolve `S2A` / `S2B` / `S2C` explicitly |
| `python/siac/rt/lut/backend.py` | MODIFY | interpolate and use tabulated SRFs in convolution |
| `tests/unit/test_srf.py` | NEW | SRF parsing, normalization, repository tests |
| `tests/unit/test_srf_kernel.py` | NEW | LUT-grid alignment and support-window tests |
| `tests/integration/test_srf_runtime.py` | NEW | end-to-end platform-specific SRF usage |

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

| Backend | Speed | Jacobian | Output contract | Sensor coverage |
|---------|-------|----------|-----------------|-----------------|
| `EmulatorBackend` | Fast | Analytical | `RTCoefficients` | S2, L8 (pre-trained) |
| `CoefficientLUTBackend` | Medium | Numerical | `RTCoefficients` | Multispectral sensors with pre-banded LUT support |
| `SpectralLUTBackend` | Medium-Slow | Numerical | `RTCoefficients` after SRF convolution | Hyperspectral or sensors needing on-the-fly SRF convolution |
| `Py6SBackend` | Slow | Numerical | Usually `RTCoefficients` | Any sensor |

Two distinct planning axes must be kept explicit:

1. **RT model family**: emulator, compact LUT, dense spectral LUT, line-by-line / Py6S.
2. **AC method**: built-in coefficient-space correction, or custom user correction.

These axes are related but not identical:
- an emulator almost always feeds coefficient correction;
- a dense spectral LUT may compute `path_ref/T_total/S` internally, then collapse back to coefficients before M6;
- Py6S may compute coefficients directly or derive them from more detailed intermediate terms, but the public boundary should still be `RTCoefficients`.

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
| 3 | `CoefficientLUTBackend` numerical jacobian finite | no NaN/Inf in output |
| 4 | `SpectralLUTBackend` bandpass convolution path | derived `xap/xbp/xcp` are stable and finite |
| 5 | `compute_coefficients` returns valid `RTCoefficients` | all fields present, physically bounded |
| 6 | Custom backend via user class | user class used, result validated |
| 7 | `is_available_for_sensor` rejects unsupported | returns `False` for unknown sensor |

---

## 11. Centralised Authentication

SIAC v2 accesses multiple remote data sources (CDSE, CAMS/CDS API, AWS S3, NASA Earthdata, GCS). The authentication layer is split into two responsibilities:

1. **Credential store / precedence resolution**
   - a single **`CredentialManager`** (`siac.adapters.auth`) owns credential loading from resolved config
   - it stores raw credentials only and acts as a factory for provider-specific auth adapters
2. **Provider-specific authentication adapters**
   - small source-specific helpers translate raw credentials into the exact shape needed by that provider:
     - OAuth token exchange for CDSE
     - `cdsapi.Client(...)` kwargs and external-config detection for CDS
     - `storage_options` for AWS/GCS
     - temporary environment activation and `EarthAccessSource` construction for Earthdata

This keeps secret discovery centralised while preventing provider-specific auth behavior from leaking into orchestration and provider modules.

### Key types

| Type | Location | Purpose |
|---|---|---|
| `CredentialSpec` | `core/auth.py` | Frozen `(key, secret)` pair |
| `OAuthToken` | `core/auth.py` | Cached bearer token with monotonic expiry |
| `CredentialManager` | `core/auth.py` | Secret registry + factory for auth adapters |
| `CDSEAuth` | `core/auth.py` | CDSE OAuth2 token cache and bearer-header helper |
| `CDSAuth` | `core/auth.py` | CDS client kwargs + external credential detection |
| `AWSAuth` | `core/auth.py` | S3 / fsspec storage option builder |
| `GCSAuth` | `core/auth.py` | GCS / fsspec storage option builder |
| `EarthdataAuth` | `core/auth.py` | Earthdata env activation + `EarthAccessSource` factory |
| `AuthConfig` | `config/schema.py` | Pydantic model for centralized auth config |
| `AuthenticationError` | `core/exceptions.py` | Raised on missing/failed auth |

### Credential resolution order (in `from_config`)

1. `SIACConfig.auth` fields
2. Environment overlays applied centrally by `siac.config.load.overlay_env_secrets`

The store remains the single place that resolves precedence. Adapters must not reimplement config/env parsing or provider-native rc-file discovery.

### Env-var mapping

| Provider | Key env var | Secret env var |
|---|---|---|
| CDSE | `SIAC_CDSE_USERNAME` | `SIAC_CDSE_PASSWORD` |
| CDS API | `CDSAPI_KEY` | — |
| AWS/S3 | `AWS_ACCESS_KEY_ID` | `AWS_SECRET_ACCESS_KEY` |
| Earthdata | `EARTHDATA_USERNAME` | `EARTHDATA_PASSWORD` |
| GCS | `GOOGLE_APPLICATION_CREDENTIALS` | — |

`CDSE` keeps SIAC-scoped variables because it does not have a stable de facto standard comparable to the other providers.
The same adapter should also mint temporary S3 credentials for `s3://eodata/...`
access so CDSE catalogue and object-store usage stay under one auth boundary.

### Adapter contracts

Adapters expose typed operations instead of generic `(key, secret)` plumbing:

- `CDSEAuth.get_token()` / `authorization_header()`
- `CDSEAuth.create_temporary_s3_credentials()` / `temporary_s3_credentials()`
- `CDSAuth.client_kwargs()`
- `AWSAuth.storage_options()`
- `GCSAuth.storage_options()`
- `EarthdataAuth.activate_environment()` / `build_earthaccess_source()`

Provider modules and backends should depend on those typed operations, not on direct `get_credentials()` branching.

### Integration points

- **`SIAC.__init__`** creates `self._auth = CredentialManager.from_config(config)`
- **`siac_process()`** accepts optional `auth` kwarg; defaults to `from_config`
- **`CopernicusDataspaceBackend`** uses `CDSEAuth` for OAuth2 bearer tokens
- **`CAMSProvider`** uses `CDSAuth` for `cdsapi` client kwargs
- **`EarthAccessSource`** is built through `EarthdataAuth`
- **`_resolve_rt_model_for_pipeline`** uses `AWSAuth` to inject S3 `storage_options` for LUT paths

### Guardrails

- `CredentialManager` may cache raw credentials and provider adapters, but only provider adapters own provider-specific behavior
- token exchange must not be embedded in search/download backends directly
- provider modules should not parse SIAC environment variables themselves
- provider modules should not honor provider-native rc files outside the centralized config/env overlay path

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
