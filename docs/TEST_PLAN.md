# SIAC v2 — Core Test Plan

This document defines the test strategy for the sensor-agnostic core architecture
described in [PLANS.md](PLANS.md). It covers what to test, how to test it, and
what test infrastructure is needed.

**Scope**: core pipeline only (M1–M6 contracts, orchestration, validation,
injection, spectral model). Sensor-specific tests (S2 SAFE parsing, CDSE download)
belong in separate sensor test plans.

---

## Table of Contents

1. [Test Organization](#1-test-organization)
2. [Test Fixtures & Mock Factory](#2-test-fixtures--mock-factory)
3. [Layer 1 — Contract Types (Unit)](#3-layer-1--contract-types-unit)
4. [Layer 2 — Validation Functions (Unit)](#4-layer-2--validation-functions-unit)
5. [Layer 3 — Module Implementations (Unit)](#5-layer-3--module-implementations-unit)
6. [Layer 4 — Pipeline Orchestration (Integration)](#6-layer-4--pipeline-orchestration-integration)
7. [Layer 5 — Custom Injection (Integration)](#7-layer-5--custom-injection-integration)
8. [Layer 6 — Spectral Model (Unit)](#8-layer-6--spectral-model-unit)
9. [Layer 7 — Config Resolution (Unit)](#9-layer-7--config-resolution-unit)
10. [Layer 8 — End-to-End with Synthetic Data (Regression)](#10-layer-8--end-to-end-with-synthetic-data-regression)
11. [Layer 9 — Performance & Concurrency](#11-layer-9--performance--concurrency)
12. [Test Matrix Summary](#12-test-matrix-summary)
13. [Fixture Dependency Graph](#13-fixture-dependency-graph)
14. [Running the Tests](#14-running-the-tests)

---

## 1. Test Organization

### Directory layout

```
tests/
├── conftest.py                 ← shared fixtures (mock_observation_bundle, etc.)
├── fixtures/
│   ├── sample_data/            ← tiny real data snippets (< 5 MB total)
│   └── synthetic/              ← generated NetCDF / GeoTIFF for reproducible tests
├── unit/
│   ├── test_contracts.py       ← Layer 1: contract type construction + properties
│   ├── test_validation.py      ← Layer 2: validation functions
│   ├── test_grid_assembler.py  ← Layer 3: M4
│   ├── test_solver.py          ← Layer 3: M5 (already exists)
│   ├── test_correction.py      ← Layer 3: M6 (already exists)
│   ├── test_spectral.py        ← Layer 6: spectral model
│   ├── test_config.py          ← Layer 7: config resolution (already exists)
│   ├── test_resolve.py         ← Layer 7: _resolve_* helpers
│   ├── test_types.py           ← existing, extend for new types
│   └── ...                     ← existing files kept as-is
├── integration/
│   ├── test_pipeline.py        ← Layer 4: orchestration (extend existing)
│   ├── test_injection.py       ← Layer 5: custom provider injection
│   └── test_e2e_synthetic.py   ← Layer 8: full pipeline with synthetic data
├── regression/
│   └── test_known_outputs.py   ← Layer 8: known-answer tests
└── benchmarks/
    └── test_perf.py            ← Layer 9: timing + concurrency
```

### Markers

| Marker | Meaning | CI gate? |
|--------|---------|----------|
| *(none)* | Fast unit test (< 1 s) | Yes |
| `@pytest.mark.slow` | Slower unit/integration (< 30 s) | Yes |
| `@pytest.mark.integration` | Requires multiple modules interacting | Yes |
| `@pytest.mark.regression` | Known-answer comparisons | Yes |
| `@pytest.mark.benchmark` | Performance timing — not a pass/fail gate | No |
| `@pytest.mark.real_data` | Requires real satellite data on disk | No |

### Naming convention

```
test_{module}_{what}_{expected_outcome}
```

Examples:
- `test_observation_bundle_missing_time_raises`
- `test_validate_atmo_state_negative_unc_fails`
- `test_pipeline_custom_preprocessor_called`
- `test_assemble_grids_spatial_alignment`

---

## 2. Test Fixtures & Mock Factory

All mock data objects should be built from a centralised factory in `conftest.py`
so that every test uses consistent, valid baseline data. Tests that need invalid
data start from a valid fixture and corrupt one field.

### 2.1 Core Fixtures to Add

```python
# ── conftest.py additions ──────────────────────────────────────────

SHAPE = (32, 32)      # fast but large enough for multi-grid solver
N_BANDS = 3           # B02, B03, B04

@pytest.fixture
def mock_sensor_config() -> SensorConfig:
    """3-band sensor config (Blue, Green, Red)."""
    ...

@pytest.fixture
def mock_observation_bundle(mock_sensor_config) -> ObservationBundle:
    """Complete valid ObservationBundle with 32×32 synthetic TOA."""
    ...

@pytest.fixture
def mock_atmospheric_state() -> AtmosphericState:
    """Spatially uniform AtmosphericState at 32×32."""
    ...

@pytest.fixture
def mock_surface_prior() -> SurfacePrior:
    """Uniform SurfacePrior at 32×32."""
    ...

@pytest.fixture
def mock_solver_input_bundle(
    mock_observation_bundle,
    mock_atmospheric_state,
    mock_surface_prior,
    mock_rt_model,
) -> SolverInputBundle:
    """Pre-assembled SolverInputBundle (skips M4)."""
    ...

@pytest.fixture
def mock_solved_atmosphere(mock_atmospheric_state) -> SolvedAtmosphere:
    """Dummy solved output (converged, 5 iterations)."""
    ...
```

### 2.2 Mock Module Callables

Each module gets a mock that returns a valid contract (immediate, no I/O):

```python
@pytest.fixture
def mock_preprocessor(mock_observation_bundle) -> Callable:
    """(Path, AOI) -> ObservationBundle. Returns fixed bundle."""
    def _preprocess(input_path, aoi=None):
        return mock_observation_bundle
    return _preprocess

@pytest.fixture
def mock_atmo_provider(mock_atmospheric_state) -> Callable:
    """(bounds, crs, time, res) -> AtmosphericState."""
    def _get_prior(bounds, crs, obs_time, resolution):
        return mock_atmospheric_state
    return _get_prior

@pytest.fixture
def mock_surface_prior_provider(mock_surface_prior) -> Callable:
    """(bounds, crs, time, sensor_config, geometry, res) -> SurfacePrior."""
    def _get_surface_prior(bounds, crs, obs_time, sensor_config, geometry, resolution):
        return mock_surface_prior
    return _get_surface_prior
```

### 2.3 Fixture Reuse Principle

> Never hard-code array shapes or values inside individual test functions.
> Always start from a fixture and mutate only the field under test.

This makes tests resilient to shape changes (e.g. switching from 32×32 to 64×64).

---

## 3. Layer 1 — Contract Types (Unit)

**File**: `tests/unit/test_contracts.py`

**Purpose**: verify that each frozen dataclass can be constructed with valid data,
rejects obviously invalid data, and exposes the expected properties.

### 3.1 `ObservationBundle`

| Test ID | What | Assert |
|---------|------|--------|
| `test_obs_bundle_construction` | Build from mock fixtures | All fields accessible, types correct |
| `test_obs_bundle_frozen` | Attempt attribute assignment | `FrozenInstanceError` raised |
| `test_obs_bundle_metadata_time` | Check `metadata["observation_time"]` | Is `datetime` |
| `test_obs_bundle_bounds_tuple` | Check `bounds` | len 4, floats |
| `test_obs_bundle_crs_string` | Check `crs` | Non-empty string starting with `"EPSG:"` |

### 3.2 `AtmosphericState`

| Test ID | What | Assert |
|---------|------|--------|
| `test_atmo_state_construction` | Build from arrays | All 7 fields populated |
| `test_atmo_state_with_updated_aot_tcwv` | Call `with_updated_aot_tcwv()` | New object, old unchanged, tco3 preserved |
| `test_atmo_state_shape_consistency` | All fields | Same spatial shape |

### 3.3 `SurfacePrior`

| Test ID | What | Assert |
|---------|------|--------|
| `test_surface_prior_construction` | Build with boa, boa_unc, mask | Fields match |
| `test_surface_prior_mask_dtype` | Check mask type | Boolean |
| `test_surface_prior_boa_unc_nonneg` | Check boa_unc ≥ 0 | True |

### 3.4 `SolverInputBundle`

| Test ID | What | Assert |
|---------|------|--------|
| `test_solver_input_bundle_construction` | Build from upstream mocks | Complete, all fields present |
| `test_solver_input_bundle_bands_subset` | `bands` ⊂ `sensor_config.bands` | True |
| `test_solver_input_bundle_resolution_metadata` | `aux_resolution_m`, `aerosol_resolution_m` | Positive floats |

### 3.5 `SolvedAtmosphere`

| Test ID | What | Assert |
|---------|------|--------|
| `test_solved_atmo_construction` | Build with solved fields | Converged flag, cost, n_iterations |
| `test_solved_atmo_state_has_updated_aot` | `atmo_state.aot` == `aot` | Values match |
| `test_solved_atmo_diagnostics_types` | `cost_final` float, `n_iterations` int, `converged` bool | Type checks |

### 3.6 `CorrectionResult`

| Test ID | What | Assert |
|---------|------|--------|
| `test_correction_result_construction` | Build from boa Dataset + metadata | All fields present |
| `test_correction_result_boa_bands` | `boa` dataset has expected band vars | Band names match |
| `test_correction_result_optional_unc` | `boa_unc=None` is valid | No error |

---

## 4. Layer 2 — Validation Functions (Unit)

**File**: `tests/unit/test_validation.py`

**Purpose**: verify that the contract validation functions (`_validate_*`) catch
every documented violation and accept all valid inputs.

Each validator has a **happy path** (valid fixture passes) and a set of
**negative cases** (one broken field per test).

### 4.1 `_validate_observation_bundle`

| Test ID | Corrupted field | Expected error |
|---------|----------------|----------------|
| `test_validate_obs_happy` | None (valid fixture) | No error |
| `test_validate_obs_missing_time` | Remove `metadata["observation_time"]` | `"metadata must include 'observation_time'"` |
| `test_validate_obs_wrong_time_type` | `metadata["observation_time"] = "2024-01-01"` (str) | Type assertion |
| `test_validate_obs_empty_toa` | `toa = xr.Dataset()` | `"toa must have spatial dimensions"` |
| `test_validate_obs_cloud_shape_mismatch` | `cloud_mask` shape ≠ TOA spatial shape | Shape assertion |
| `test_validate_obs_geometry_not_broadcastable` | Geometry shape incompatible with TOA | Broadcast error |

### 4.2 `_validate_atmospheric_state`

| Test ID | Corrupted field | Expected error |
|---------|----------------|----------------|
| `test_validate_atmo_happy` | None | No error |
| `test_validate_atmo_negative_aot_unc` | `aot_unc` contains –0.1 | `"uncertainties must be non-negative"` |
| `test_validate_atmo_negative_tcwv_unc` | `tcwv_unc` contains –0.5 | Same |
| `test_validate_atmo_negative_tco3_unc` | `tco3_unc` contains –0.01 | Same |
| `test_validate_atmo_nan_in_aot` | `aot` contains NaN | Warning or error (decide policy) |

### 4.3 `_validate_surface_prior`

| Test ID | Corrupted field | Expected error |
|---------|----------------|----------------|
| `test_validate_prior_happy` | None | No error |
| `test_validate_prior_shape_mismatch` | `boa.shape ≠ boa_unc.shape` | Shape assertion |
| `test_validate_prior_mask_broadcast` | `mask` not broadcastable to `boa` | Broadcast error |

### 4.4 `_validate_solver_input_bundle`

| Test ID | Corrupted field | Expected error |
|---------|----------------|----------------|
| `test_validate_sib_happy` | None | No error |
| `test_validate_sib_misaligned_grids` | `toa` on different spatial grid than `atmo_prior.aot` | Grid alignment error |
| `test_validate_sib_bands_not_in_config` | `bands` contains band not in `sensor_config` | Subset validation |
| `test_validate_sib_rt_unsupported_band` | `rt_model.is_available_for_sensor()` returns False | RT availability |

---

## 5. Layer 3 — Module Implementations (Unit)

### 5.1 M4: Grid Assembler — `assemble_grids()`

**File**: `tests/unit/test_grid_assembler.py`

| Test ID | What | Assert |
|---------|------|--------|
| `test_assemble_grids_output_type` | Call with valid upstream outputs | Returns `SolverInputBundle` |
| `test_assemble_grids_spatial_alignment` | Check that all raster fields share the same spatial grid | `toa`, `atmo_prior.aot`, `surface_prior.boa` grids match |
| `test_assemble_grids_aux_resolution` | `aux_resolution_m` matches config | 500.0 (default) |
| `test_assemble_grids_band_selection` | Bands selected by wavelength 400–520 nm | Only B01/B02-like bands |
| `test_assemble_grids_cloud_mask_conservative` | If any native pixel in footprint is cloud → aux pixel is cloud | Conservative OR |
| `test_assemble_grids_geometry_bilinear` | Angles interpolated (not nearest) | Smooth gradient preserved |
| `test_assemble_grids_identity_same_res` | All inputs already at aux resolution | Output ≈ input (no resampling artefacts) |
| `test_assemble_grids_multi_res_input` | TOA with mixed 10/20/60 m bands | All resampled to aux without error |
| `test_assemble_grids_passes_validation` | Output passes `_validate_solver_input_bundle()` | No error |

### 5.2 M5: Aerosol Solver — `solve_aerosol()`

Existing tests in `test_solver.py` cover the multi-grid solver and cost function.
Add contract-level tests:

| Test ID | What | Assert |
|---------|------|--------|
| `test_solve_aerosol_output_type` | Call with mock `SolverInputBundle` | Returns `SolvedAtmosphere` |
| `test_solve_aerosol_aot_positive` | Solved AOT | ≥ 0 everywhere |
| `test_solve_aerosol_tcwv_positive` | Solved TCWV | ≥ 0 everywhere |
| `test_solve_aerosol_convergence_flag` | Sufficient iterations | `converged=True` |
| `test_solve_aerosol_atmo_state_updated` | `solved.atmo_state.aot` | ≠ prior (solver changed it) |
| `test_solve_aerosol_known_answer` | Pre-computed forward-model scene | AOT within 20% of truth |

### 5.3 M6: Atmospheric Corrector — `correct_atmosphere()`

Existing tests in `test_correction.py` cover BOA physics. Add contract-level tests:

| Test ID | What | Assert |
|---------|------|--------|
| `test_correct_output_type` | Call with `ObservationBundle` + `SolvedAtmosphere` | Returns `CorrectionResult` |
| `test_correct_boa_bands_match_toa` | BOA dataset band names | Same as input TOA |
| `test_correct_boa_range` | BOA values | In [–0.05, 1.5] (allowing minor negatives) |
| `test_correct_cloud_mask_preserved` | Cloud mask in result | Matches original |
| `test_correct_metadata_has_timing` | `result.diagnostics` | Contains `processing_time_s` |
| `test_correct_native_resolution` | BOA spatial shape | Matches M1 native TOA shape |

### 5.4 RT Model Backend — Protocol Conformance

| Test ID | What | Assert |
|---------|------|--------|
| `test_rt_backend_protocol_mock` | `MockRTModel` satisfies `RTModelBackend` | `isinstance` or duck typing check |
| `test_rt_backend_jacobian_shape` | `compute_coefficients(..., compute_jacobian=True)` | d_xap has `param` dim with `["aot", "tcwv"]` |
| `test_rt_backend_no_jacobian` | `compute_coefficients(..., compute_jacobian=False)` | d_xap/d_xbp/d_xcp are `None` |
| `test_rt_backend_unsupported_sensor` | `is_available_for_sensor("UNKNOWN", "SAT")` | Returns `False` |

---

## 6. Layer 4 — Pipeline Orchestration (Integration)

**File**: `tests/integration/test_pipeline.py` (extend existing)

All tests in this layer use **mock module callables** — no real I/O. The purpose
is to verify that `run_pipeline()` and `siac_process()` wire modules correctly.

### 6.1 `run_pipeline()` — Happy Path

| Test ID | What | Assert |
|---------|------|--------|
| `test_pipeline_happy_path` | All six mock modules → `CorrectionResult` | No error, valid result |
| `test_pipeline_calls_all_modules` | Wrap each mock with a call counter | Each mock called exactly once |
| `test_pipeline_module_call_order` | Record call timestamps | M1 before M4; M4 before M5; M5 before M6 |
| `test_pipeline_m2_m3_after_m1` | Record call times | M2 and M3 start after M1 returns (need bounds) |
| `test_pipeline_returns_correction_result` | Check returned type | `isinstance(result, CorrectionResult)` |

### 6.2 `run_pipeline()` — Validation Integration

| Test ID | What | Assert |
|---------|------|--------|
| `test_pipeline_invalid_m1_output` | M1 mock returns bundle missing `observation_time` | Clear error before M4 |
| `test_pipeline_invalid_m2_output` | M2 mock returns negative uncertainties | Clear error before M4 |
| `test_pipeline_invalid_m3_output` | M3 mock returns shape mismatch | Clear error before M4 |
| `test_pipeline_m1_exception` | M1 mock raises `FileNotFoundError` | Exception propagates with useful traceback |
| `test_pipeline_m2_exception` | M2 mock raises `ConnectionError` | Exception propagates |

### 6.3 `siac_process()` — Defaults Resolution

| Test ID | What | Assert |
|---------|------|--------|
| `test_siac_process_resolves_defaults` | Call with `config` only (no overrides) | `_resolve_*` helpers called, pipeline runs |
| `test_siac_process_explicit_overrides` | Pass custom `atmo_provider` | Config-driven default NOT called |
| `test_siac_process_partial_override` | Override only M2, rest default | Only M2 custom called, M1/M3 default |
| `test_siac_process_none_defaults` | All explicit args are `None` | Treated as "use config default" |

---

## 7. Layer 5 — Custom Injection (Integration)

**File**: `tests/integration/test_injection.py`

**Purpose**: prove that users can replace any module with a custom callable (function
or class instance) and the pipeline still works end-to-end.

### 7.1 Function Injection

| Test ID | What | Assert |
|---------|------|--------|
| `test_inject_preprocessor_function` | Pass `my_loader(path, aoi) -> ObservationBundle` | Pipeline uses it, returns valid result |
| `test_inject_atmo_function` | Pass `constant_atmo(b,c,t,r) -> AtmosphericState` | Pipeline uses it |
| `test_inject_surface_prior_function` | Pass custom function returning `SurfacePrior` | Pipeline uses it |
| `test_inject_solver_function` | Pass passthrough solver (fixed AOT) | Solver skipped, correction uses fixed AOT |
| `test_inject_corrector_function` | Pass identity corrector (BOA = TOA) | Output BOA ≈ input TOA |
| `test_inject_grid_assembler_function` | Pass custom assembler | Resampling logic replaced |

### 7.2 Class Instance Injection

| Test ID | What | Assert |
|---------|------|--------|
| `test_inject_preprocessor_class` | Pass `MySensor().preprocess` as bound method | Pipeline calls it |
| `test_inject_atmo_provider_class` | Pass `CAMSProvider(tmp_cams_dir).get_prior` | Pipeline calls bound method |
| `test_inject_subclassed_preprocessor` | Subclass `Sentinel2Preprocessor`, override `extract_cloud_mask` | Overridden method called, rest inherited |

### 7.3 Contract Violation by Custom Provider

| Test ID | What | Assert |
|---------|------|--------|
| `test_inject_bad_preprocessor_missing_field` | Custom preprocessor returns incomplete `ObservationBundle` | Validation error with clear message |
| `test_inject_bad_atmo_wrong_type` | Custom atmo returns a dict instead of `AtmosphericState` | `TypeError` or `AttributeError` |
| `test_inject_bad_solver_no_converged` | Custom solver returns `SolvedAtmosphere` without `converged` | `AttributeError` |

---

## 8. Layer 6 — Spectral Model (Unit)

**File**: `tests/unit/test_spectral.py`

### 8.1 `SensorBand`

| Test ID | What | Assert |
|---------|------|--------|
| `test_gaussian_only_construction` | Create with center + FWHM only | `has_rsrf == False` |
| `test_with_rsrf` | Create with tabulated RSRF | `has_rsrf == True` |
| `test_band_descriptor_wavelength_um` | 550 nm → 0.55 µm | Property correct |

### 8.2 `SensorConfig` Band Selection

| Test ID | What | Assert |
|---------|------|--------|
| `test_select_bands_visible` | `vis_bands` (400–700 nm) on S2 config | Returns B01, B02, B03, B04 |
| `test_select_bands_nir` | `nir_bands` (750–1000 nm) | Returns B07, B08, B8A, B09 |
| `test_select_bands_swir` | `swir_bands` (1000–2500 nm) | Returns B11, B12 (excludes B10 cirrus) |
| `test_select_nearest_band` | Target 660 nm, tol 20 nm | Returns B04 (665 nm) |
| `test_select_nearest_band_no_match` | Target 1240 nm, tol 20 nm | Returns `None` |
| `test_select_bands_empty_range` | Range 2500–3000 nm | Empty list |

### 8.3 Spectral Convolution Functions

| Test ID | What | Assert |
|---------|------|--------|
| `test_sensor_to_reference_identity` | Input sensor ≈ reference bands | Output ≈ input |
| `test_reference_to_sensor_roundtrip` | `sensor_to_reference` → `reference_to_sensor` | Approximates original (within tolerance) |
| `test_sensor_to_reference_flat_spectrum` | Flat reflectance = 0.5 across all λ | All reference bands ≈ 0.5 |
| `test_sensor_to_reference_shape` | 3-band input | Output has 7 MODIS reference bands |
| `test_load_modis` | Load MODIS RSRF | 7 bands, wavelengths in [0.4, 2.5] µm |
| `test_load_unknown_raises` | Unknown sensor name | `ValueError` |

---

## 9. Layer 7 — Config Resolution (Unit)

**File**: `tests/unit/test_resolve.py`

**Purpose**: test the `_resolve_*` helper functions that translate config values
into module callables.

| Test ID | What | Assert |
|---------|------|--------|
| `test_resolve_preprocessor_sentinel2` | `config.sensor = "sentinel2"` | Returns `Sentinel2Preprocessor().preprocess` (callable) |
| `test_resolve_preprocessor_unknown` | `config.sensor = "unknown"` | `ValueError` |
| `test_resolve_atmo_cams` | `config.providers.atmo.kind = "cams"` | Returns callable that wraps `CAMSProvider` |
| `test_resolve_atmo_merra2` | `config.providers.atmo.kind = "merra2"` | Returns callable that wraps `MERRA2Provider` |
| `test_resolve_atmo_unknown` | `config.providers.atmo.kind = "xxx"` | `ValueError` |
| `test_resolve_surface_prior_brdf` | `config.algorithms.surface_prior.method = "kernel_model"` | Returns callable wrapping kernel-model surface priors |
| `test_resolve_surface_prior_store` | `config.algorithms.surface_prior.method = "monthly_database"` | Returns callable wrapping monthly database priors |
| `test_resolve_rt_model_emulator` | `config.algorithms.rt.backend = "emulator"` | Returns `EmulatorBackend` instance |
| `test_resolve_rt_model_lut` | `config.algorithms.rt.backend = "lut"` | Returns `LUTBackend` instance |
| `test_resolve_returns_callable` | All `_resolve_*` helpers | Return value is `callable()` |

---

## 10. Layer 8 — End-to-End with Synthetic Data (Regression)

**File**: `tests/integration/test_e2e_synthetic.py`  
**File**: `tests/regression/test_known_outputs.py`

### 10.1 Synthetic Scene Pipeline

Build a fully synthetic scene with **known** physical parameters (flat surface,
uniform atmosphere, known geometry) so that the correct BOA is analytically
computable.

| Test ID | What | Assert |
|---------|------|--------|
| `test_e2e_synthetic_flat_surface` | Uniform ρ=0.2, AOT=0.1, SZA=30° | BOA within 5% of truth |
| `test_e2e_synthetic_high_aot` | Same but AOT=0.8 | BOA differs from low-AOT case |
| `test_e2e_synthetic_cloudy_pixels` | 25% cloud mask | Cloudy pixels flagged, clear pixels corrected |
| `test_e2e_synthetic_aoi_crop` | Provide AOI covering half the scene | Result spatial extent matches AOI |
| `test_e2e_synthetic_all_bands` | 13-band S2-like sensor config | All bands corrected, no NaN in clear pixels |
| `test_e2e_synthetic_deterministic` | Run twice with same seed | Bit-identical outputs |

### 10.2 Known-Answer Regression

Store a small reference output (AOT, TCWV, BOA) as `.npz` in `tests/fixtures/`.
Each regression test loads the same input, runs the pipeline, and compares.

| Test ID | What | Tolerance |
|---------|------|-----------|
| `test_regression_aot_reference` | Solved AOT matches stored reference | RMSE < 0.01 |
| `test_regression_tcwv_reference` | Solved TCWV matches stored reference | RMSE < 0.05 |
| `test_regression_boa_reference` | Corrected BOA matches stored reference | RMSE < 0.005 per band |

**Maintenance rule**: when the algorithm intentionally changes (e.g. new RT
emulator weights), regenerate references and commit with a changelog entry
explaining the expected delta.

---

## 11. Layer 9 — Performance & Concurrency

**File**: `tests/benchmarks/test_perf.py`

### 11.1 Concurrency Tests

| Test ID | What | Assert |
|---------|------|--------|
| `test_m2_m3_run_concurrently` | Instrument M2/M3 mocks with 0.5 s sleep each | Total time < 0.8 s (not 1.0 s) |
| `test_pipeline_sequential_fallback` | Set `max_workers=1` | Pipeline still works, M2/M3 sequential |
| `test_concurrent_exception_handling` | M2 mock raises after 0.2 s | Exception propagated, M3 not swallowed |

### 11.2 Performance Baselines

| Test ID | What | Baseline |
|---------|------|----------|
| `test_perf_assemble_grids_32x32` | Time `assemble_grids()` on 32×32 | < 0.5 s |
| `test_perf_solver_32x32_3band` | Time `solve_aerosol()` on 32×32, 3 bands | < 5 s |
| `test_perf_corrector_32x32` | Time `correct_atmosphere()` on 32×32 | < 1 s |
| `test_perf_full_pipeline_32x32` | Time `siac_process()` end-to-end | < 10 s |

Performance tests record timing but do not fail CI — they alert on significant
regressions (> 2× slower).

---

## 12. Test Matrix Summary

| Layer | File(s) | # Tests (est.) | I/O? | Marker |
|-------|---------|---------------|------|--------|
| 1. Contract types | `test_contracts.py` | ~20 | No | *(none)* |
| 2. Validation | `test_validation.py` | ~18 | No | *(none)* |
| 3. Module impls | `test_grid_assembler.py`, existing solver/correction | ~25 | No | *(none)* |
| 4. Orchestration | `test_pipeline.py` | ~12 | No | `integration` |
| 5. Injection | `test_injection.py` | ~12 | No | `integration` |
| 6. Spectral | `test_spectral.py` | ~12 | No | *(none)* |
| 7. Config resolve | `test_resolve.py` | ~10 | No | *(none)* |
| 8. E2E + regression | `test_e2e_synthetic.py`, `test_known_outputs.py` | ~9 | Minimal | `regression` |
| 9. Performance | `test_perf.py` | ~7 | No | `benchmark` |
| **Total** | | **~125** | | |

All Layer 1–7 tests should run in < 30 s total on a single core.

---

## 13. Fixture Dependency Graph

```
mock_sensor_config
    │
    ├── mock_observation_bundle
    │       │
    │       ├── mock_preprocessor
    │       │
    │       ├── mock_solver_input_bundle ─── mock_rt_model
    │       │       │
    │       │       ├── mock_solved_atmosphere
    │       │       │
    │       │       └── (used in: test_grid_assembler, test_solver, test_pipeline)
    │       │
    │       └── (used in: test_contracts, test_validation, test_injection)
    │
    ├── mock_atmospheric_state
    │       │
    │       ├── mock_atmo_provider
    │       │
    │       └── (used in: test_contracts, test_validation, test_pipeline)
    │
    └── mock_surface_prior
            │
            ├── mock_surface_prior_provider
            │
            └── (used in: test_contracts, test_validation, test_pipeline)
```

---

## 14. Running the Tests

```bash
# Build the native module first in a fresh environment.
pixi run build-rust

# All fast tests (Layers 1–7)
pixi run test-fast

# Quick smoke test (< 10 s)
pixi run pytest tests/unit/test_contracts.py tests/unit/test_validation.py -x -q

# Integration only
pixi run pytest tests/integration/ -m integration

# Regression only
pixi run pytest tests/regression/ -m regression

# Performance benchmarks (informational)
pixi run pytest tests/benchmarks/ -m benchmark --benchmark-only

# Full suite
pixi run test

# With strict coverage gates (global >=95%, each file >90%)
pixi run coverage
```

### CI Configuration

```yaml
# .github/workflows/test.yml  (excerpt)
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: prefix-dev/setup-pixi@v0
      - run: pixi install
      - run: pixi run build-rust
      - run: pixi run test-fast
```
