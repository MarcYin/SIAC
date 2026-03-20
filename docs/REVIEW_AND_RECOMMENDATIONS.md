# SIAC v2 — Code Review & Recommendations

**Date:** 2025-07-10  
**Scope:** Full review of implementation against `PLANS.md` architecture document  
**Test status after changes:** 526/526 passing (499 unit + 27 integration)

---

## 1. Changes Already Made

The following concrete improvements were implemented and verified:

### 1.1 Dead code in `compute_boa_jacobian` (types.py)

Line `1.0 / (denom**2)` was computed but never assigned. Replaced with a `# TODO` comment explaining the derivative derivation. This was an unfinished Jacobian pathway needed for uncertainty propagation.

### 1.2 Cloud mask contract fix (atmospheric.py + test)

The corrector was computing a combined quality mask (`True = valid`) and storing it in `CorrectionResult.cloud_mask`. The plan and `ObservationBundle` contract say `True = cloudy/invalid`. The corrector now preserves the input cloud mask directly, and creates a fallback `True = cloudy` mask from invalid BOA pixels when no cloud mask is provided. Updated corresponding test.

### 1.3 BOA range relaxation (atmospheric.py)

Changed BOA validity threshold from `boa > 0` to `boa > -0.05` to allow minor negative reflectance per plan §4.6 physics (adjacency correction, very dark targets).

### 1.4 Empty-bands guard (atmospheric.py)

Added a `ValueError` when no TOA bands match the sensor config, providing diagnostic band names for debugging.

### 1.5 Pipeline validation for M4–M6 (pipeline.py + validation.py)

Added `_validate_solved_atmosphere()` and `_validate_correction_result()` validators. Pipeline now validates outputs at every stage (M1–M6), matching the plan's contract enforcement strategy.

### 1.6 `SurfacePrior.kernels` made optional (types.py)

Changed from `BRDFKernelWeights` to `BRDFKernelWeights | None` so custom surface-prior providers that don't use BRDF can return `SurfacePrior(kernels=None)`.

### 1.7 Geometry broadcastability check (validation.py)

Added a validation that SZA, VZA, RAA shapes are mutually broadcastable in `_validate_observation_bundle()`.

### 1.8 Provider recreation fix (siac.py)

`_resolve_surface_prior_provider()` was creating a new `BRDFProductProvider` and `KernelModelDeriver` on every call. Fixed to create once and close over them.

---

## 2. Gaps Between Plan and Implementation

### 2.1 `AerosolSolver` Protocol vs Plan §4.5

| Aspect | Plan | Implementation |
|---|---|---|
| Input | `SolverInputBundle` (frozen dataclass) | Individual args: `toa, surface_prior, geometry, atmo_prior, rt_model, cloud_mask, solver_config` |
| Output | `SolvedAtmosphere` | `AtmosphericState` |

**Impact:** Medium. The protocol accepts loose arguments rather than the documented bundle contract. This means the pipeline must destructure the bundle when calling the solver, reducing type safety. The return type should be `SolvedAtmosphere` (which includes `converged`, `n_iterations`, `diagnostics` beyond just atmospheric state fields).

**Recommendation:** Align the protocol to accept `SolverInputBundle` and return `SolvedAtmosphere`. This is a breaking change to the solver interface but affects only the `MultigridSolver` implementor and test mocks.

### 2.2 `BaseSatellitePreprocessor.preprocess()` returns `dict`, not `ObservationBundle`

The plan says M1 returns an `ObservationBundle`. The implementation returns `dict[str, Any]` with keys `toa`, `geometry`, `cloud_mask`, `metadata`. The pipeline then constructs the `ObservationBundle` itself.

**Impact:** Low-medium. The dict return means preprocessors don't get type-checked at the protocol boundary. However, the pipeline bridge works.

**Recommendation:** Have preprocessors return `ObservationBundle` directly, or add a `_to_observation_bundle()` method in `BaseSatellitePreprocessor`. Deprecate the dict return.

### 2.3 `SIAC.process()` bypasses `run_pipeline()`

`SIAC.process()` (the main user-facing method) runs its own orchestration loop that doesn't use the modular `run_pipeline()` function. This means:
- Pipeline validations don't apply
- The thread/dask concurrency backends aren't used
- The retry/timeout/LUT-preloading logic is unused

**Impact:** High. Two divergent execution paths increases maintenance burden and means pipeline improvements don't propagate.

**Recommendation:** Refactor `SIAC.process()` to delegate to `run_pipeline()` internally, bridging the config-oriented API to the functional pipeline interface.

### 2.4 Hardcoded `sensor_id="MSI", satellite_id="S2A"` in emulator path

`_resolve_rt_model_for_pipeline()` hardcodes Sentinel-2A when using the emulator backend. Landsat 8/9, Sentinel-2B/2C users would get wrong spectral response.

**Impact:** High for multi-sensor support.

**Recommendation:** Pass `sensor_config.sensor_id` and `sensor_config.satellite_id` (or similar) through the RT config to the emulator factory.

### 2.5 No Sentinel-2C or Landsat 9 sensor configs

`types.py` has `SENTINEL2A_CONFIG`, `SENTINEL2B_CONFIG`, `LANDSAT8_CONFIG`, but not Sentinel-2C (launched 2024) or Landsat 9. Plan §9.3 indicates extensibility.

**Recommendation:** Add `SENTINEL2C_CONFIG` and `LANDSAT9_CONFIG`. S2C band specs are nearly identical to S2A/B but with updated SRFs. Landsat 9 OLI-2 is nearly identical to Landsat 8 OLI.

### 2.6 Historical spectral-band duplication

This issue has since been resolved. `SensorBand` is now the canonical spectral
band type, and the duplicate spectral-band proposal was removed from the active
implementation.

**Impact:** Historical only. Keep future work on the single `SensorBand` model.

**Recommendation:** Do not reintroduce a parallel spectral-band type.

### 2.7 Missing dask integration tests

The pipeline supports a dask backend but no integration test exercises it. Only thread-pool concurrency is tested.

**Recommendation:** Add a test with `dask.distributed.LocalCluster` (or mock) that verifies the dask code path, including `_run_pipeline_dask()`.

---

## 3. Architecture Recommendations

### 3.1 Error hierarchy

The codebase uses bare `ValueError`/`TypeError`/`AssertionError`. The plan §7 mentions structured error handling.

**Recommendation:** Create a small exception hierarchy:
```python
class SIACError(Exception): ...
class PreprocessingError(SIACError): ...
class SolverConvergenceError(SIACError): ...
class ValidationError(SIACError): ...
class ConfigurationError(SIACError): ...
```
The `exceptions.py` file already exists but only has `SIACError` and `ConfigurationError`. Extend it.

### 3.2 Validation strategy: assertions vs exceptions

Pipeline validators currently use `assert`, which can be silenced with `python -O`. For production safety:

**Recommendation:** Replace `assert` with explicit `raise ValidationError(...)` in all validators. Keep `assert` only for development-time invariants.

### 3.3 Logging

The codebase has minimal structured logging. Pipeline operations (M1–M6 timing, retry attempts, dask task submission) should log at INFO/DEBUG level.

**Recommendation:** Add `structlog` or standard `logging` with module-level loggers. The plan mentions observability; logging is the first step.

### 3.4 Type narrowing for `metadata: dict[str, Any]`

Several contracts use `metadata: dict[str, Any]` which is opaque. Consider typed metadata:
```python
@dataclass(frozen=True)
class CorrectionMetadata:
    diagnostics.processing_time_s: float
    solver_iterations: int
    bands_corrected: list[str]
    # extensible via extra: dict[str, Any]
```

### 3.5 Configuration validation at construction time

`SIACConfig` (via pydantic) validates at init, but `SensorConfig` is a plain frozen dataclass. Band overlap, negative wavelengths, or duplicate band names aren't caught.

**Recommendation:** Add a `__post_init__` to `SensorConfig` that validates:
- No duplicate band names
- All wavelengths positive
- All resolutions positive
- Band indices unique

---

## 4. Testing Recommendations

### 4.1 Current coverage

499 unit tests + 27 integration tests pass. The project targets 95% coverage (per `pyproject.toml`). Key gaps:

- **No orchestration tests** — `test_orchestration.py` doesn't exist despite plan references
- **No dask backend test** — only thread pool is exercised
- **No property-based tests** — the dataclass contracts are ideal for Hypothesis
- **No fuzz tests for edge cases** — NaN-heavy arrays, single-pixel images, 1-band sensors

### 4.2 Recommended new tests

1. **`test_orchestration.py`**: Mock all 6 modules, verify call order, verify validation gates fire between each module
2. **`test_dask_pipeline.py`**: Use `dask.distributed.LocalCluster(n_workers=1)` to exercise the dask code path
3. **`test_sensor_configs.py`**: Parametrize over all predefined configs, verify invariants (positive wavelengths, no duplicates, ordered band indices)
4. **Property tests for contracts**: Use Hypothesis to generate random `GeometryAngles`, `AtmosphericState`, validate they round-trip through validators
5. **Edge case corrector tests**: All-NaN TOA, single-pixel image, image with only cloud pixels

---

## 5. Performance Recommendations

### 5.1 Grid assembler resampling

`assembler.py` uses `scipy.ndimage.zoom` for resampling. For large tiles, this is single-threaded.

**Recommendation:** Consider `xarray.DataArray.interp()` or `rioxarray.reproject()` for dask-aware resampling, enabling chunk-parallel reprojection.

### 5.2 RT coefficient computation

Each band computes coefficients independently. For sensors with 10+ bands, this is embarrassingly parallel.

**Recommendation:** Vectorize the RT model to accept a batch of bands, or use `concurrent.futures` within the corrector to parallelize band-level correction.

### 5.3 LUT preloading

The pipeline already has `_preload_lut()` running in parallel with M2/M3. Verify this actually achieves overlap by adding timing logs.

---

## 6. Summary Priority Matrix

| Priority | Item | Effort | Impact |
|---|---|---|---|
| **P0** | Fix `SIAC.process()` to use `run_pipeline()` (#2.3) | High | Unifies execution paths |
| **P0** | Fix hardcoded S2A in emulator path (#2.4) | Low | Blocks multi-sensor |
| **P1** | Align `AerosolSolver` protocol (#2.1) | Medium | Type safety |
| **P1** | Replace `assert` with exceptions (#3.2) | Low | Production safety |
| **P1** | Add S2C + L9 configs (#2.5) | Low | Sensor coverage |
| **P1** | Keep one canonical `SensorBand` type (#2.6) | Medium | Reduces confusion |
| **P2** | Preprocessor → `ObservationBundle` return (#2.2) | Medium | Type consistency |
| **P2** | Exception hierarchy (#3.1) | Low | Better error handling |
| **P2** | Structured logging (#3.3) | Medium | Observability |
| **P2** | Typed metadata (#3.4) | Low | Clarity |
| **P2** | `SensorConfig.__post_init__` validation (#3.5) | Low | Catches bad configs |
| **P3** | Dask integration test (#2.7, #4.2) | Medium | Test coverage |
| **P3** | Property-based tests (#4.2) | Medium | Robustness |
| **P3** | Vectorized RT computation (#5.2) | High | Performance |
