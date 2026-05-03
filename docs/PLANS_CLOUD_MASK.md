# Cloud And Cloud-Shadow Module Plan

Status: implemented baseline (validated on 2026-02-15), concurrency redesign planned

## 1. Goal

Add a reusable cloud/cloud-shadow masking module that can:

1. Use a built-in default detector (OmniCloudMask).
2. Use user-provided detector output.
3. Use existing cloud/shadow files.

Standard module output must be an `xr.DataArray` with class values:

- `0`: missing / invalid / no data
- `1`: clear
- `2`: cloud
- `3`: cloud shadow

## 2. Current SIAC Constraints

Current pipeline contracts (`ObservationBundle.cloud_mask`) are boolean and interpreted as:

- `True` = cloudy/invalid (excluded)
- `False` = valid clear pixel

Solver and corrector currently apply boolean operations (`~cloud_mask`), so direct integer class masks cannot be passed through unchanged.

## 3. Contract Strategy

Introduce a dual-mask strategy:

1. **Categorical cloud classes** (new): `cloud_classes: xr.DataArray[uint8]` with values in `{0,1,2,3}`.
2. **Boolean cloud mask**: `cloud_mask: xr.DataArray[bool]` derived as:
   - `True` for classes in `{0,2,3}`
   - `False` for class `1`

This keeps a single boolean M4/M5/M6 contract while exposing richer cloud semantics for outputs and QA.

## 4. Module API

Create module package:

- `python/siac/cloud/__init__.py`
- `python/siac/cloud/mask.py`
- `python/siac/cloud/providers/omnicloudmask.py`
- `python/siac/cloud/mapping.py`

Proposed API:

```python
# primary public entrypoint
build_cloud_classes(
    toa: xr.Dataset,
    geometry: GeometryAngles | None,
    sensor_config: SensorConfig,
    *,
    config: CloudMaskAlgorithmConfig,
    input_path: Path | None = None,
) -> xr.DataArray
```

```python
# cloud-class conversion
classes_to_bool_mask(cloud_classes: xr.DataArray) -> xr.DataArray
```

## 5. Configuration Model

Add `CloudMaskAlgorithmConfig` to `AlgorithmsConfig`:

- `mode`: `"auto" | "external_file" | "user_callable" | "none"`
- `provider`: `"omnicloudmask" | "custom"`
- `external_mask_path: Path | None`
- `target_resolution_m: float = 10.0`
- `resolution_policy: "auto" | "force"`
- `allow_upsample_to_target: bool = False`
- `class_mapping: dict[int, list[int]] | None`
  - keys are target classes `{0,1,2,3}`
  - values are source classes mapped into each target class
- `unmapped_to_missing: bool = True`
- `output_name: str = "cloud_classes"`

## 6. Input Modes

### 6.1 Built-in detector (`mode="auto"`, `provider="omnicloudmask"`)

Flow:

1. Select/derive red, green, NIR reflectance from TOA.
2. Harmonize resolution to `target_resolution_m`.
3. Run OmniCloudMask.
4. Convert detector output to target classes using mapping rules.

### 6.2 Existing cloud/shadow file (`mode="external_file"`)

Flow:

1. Read raster (`read_raster`).
2. Reproject/resample to TOA reference grid (`reproject_match` with nearest).
3. Apply class mapping to target `{0,1,2,3}`.

### 6.3 User callable (`mode="user_callable"`)

Flow:

1. Invoke user function with TOA + metadata.
2. Accept either:
   - pre-standardized class array (`0..3`), or
   - custom classes requiring mapping.
3. Validate and standardize.

## 7. Mapping Rules (Required)

### 7.1 Allowed behavior

- Multiple source classes may map to one target class.
  - Example: thin cloud + thick cloud -> target class `2`.

### 7.2 Forbidden behavior

- A source class must not map to multiple target classes.

Validation algorithm:

1. Build reverse lookup `source_class -> target_class`.
2. Raise `ValueError` if a source class appears under more than one target class.
3. Ensure target keys are subset of `{0,1,2,3}`.
4. If `unmapped_to_missing=True`, unmapped source values become `0`; else raise.

## 8. Spectral Band Selection For OmniCloudMask

OmniCloudMask needs red/green/NIR reflectance.

### 8.1 Multispectral sensors

- Use canonical bands directly when available.

### 8.2 Hyperspectral or multi-band-per-color case

- Select candidate bands by wavelength windows:
  - Green: 530–590 nm
  - Red: 630–690 nm
  - NIR: 760–900 nm
- If multiple candidates exist for one color group, average them.
- If no candidate exists for a required group, raise clear error.

## 9. Resolution Policy

Default `target_resolution_m = 10`.

`resolution_policy="auto"`:

1. If native resolution is finer than target (e.g., 2 m, 5 m), downsample to target.
2. If native resolution is coarser than target (e.g., 20 m), keep native unless
   `allow_upsample_to_target=True`.

`resolution_policy="force"`:

- Always produce mask at `target_resolution_m` (upsampling allowed).

Resampling methods:

- Reflectance inputs to detector: bilinear/area (continuous).
- Class masks: nearest (categorical).

## 10. Integration Points

### 10.1 Preprocessor integration

In each sensor preprocessor (`extract_cloud_mask` stage):

1. Build `cloud_classes` via new module.
2. Convert to boolean `cloud_mask` for existing pipeline.
3. Return both in preprocess result (`cloud_mask` required, `cloud_classes` optional extra).

### 10.2 ObservationBundle extension (recommended)

Option A (minimal-risk):

- Keep `ObservationBundle` unchanged.
- Store `cloud_classes` in metadata or preprocess result only.

Option B (cleaner long-term):

- Add optional `cloud_classes: xr.DataArray | None` to `ObservationBundle`.
- Keep `cloud_mask` as the canonical boolean processing mask.

Recommended: Option B with default `None`.

### 10.3 Output writing

Write standardized cloud classes as QA raster with values `0/1/2/3`.

## 11. Error Handling

Must fail fast with explicit messages for:

- Missing required spectral groups (red/green/NIR).
- Invalid class mapping (one source mapped to multiple targets).
- Output class values outside `{0,1,2,3}` after mapping.
- Shape/CRS mismatch that cannot be reprojected.

## 12. Test Plan

Add unit tests for:

1. Mapping validation success/failure.
2. Multi-to-one mapping (thin+thick -> cloud).
3. Forbidden one-to-multiple mapping detection.
4. Unmapped class behavior (`unmapped_to_missing=True/False`).
5. Hyperspectral averaging for red/green/NIR groups.
6. Resolution policy behavior (`auto` vs `force`).
7. External mask reproject/mapping path.
8. User callable path.
9. Boolean cloud-mask conversion (`{0,2,3}` => `True`).

Integration tests:

1. Sentinel-2 preprocess with OmniCloudMask default.
2. External mask file path end-to-end through M1 -> M6.
3. Ensure solver/corrector unchanged behavior with boolean mask.

## 13. Phased Delivery

Phase 1 (infrastructure):

- Add config model, mapping utilities, class validation, bool conversion.

Phase 2 (providers):

- Add OmniCloudMask provider wrapper and external-file/user-callable modes.

Phase 3 (sensor integration):

- Wire into Sentinel-2 preprocessor; add optional ObservationBundle class field.

Phase 4 (outputs + docs):

- Save QA cloud-class output and document usage examples.

Phase 5 (current baseline complete):

- Run geometry + cloud-mask extraction concurrently inside preprocessing.
- Reuse preloaded TOA in cloud extraction to avoid duplicate reads.

Phase 6 (next redesign target):

- Move orchestration concurrency to a Dask-based execution engine.
- Add robust retries + progress/dashboard + performance report artifacts.

## 14. Open Decisions (to lock before coding)

1. Keep `ObservationBundle` unchanged vs adding `cloud_classes` field.
2. Whether default policy should upsample coarser-than-10m masks.
3. Exact expected output format from OmniCloudMask wrapper (labels vs probabilities).

## 15. Concurrency Redesign (Dask-First)

### 15.1 Decision

Use `dask.distributed` as the primary execution backend for SIAC orchestration.

Reasoning:

- Native interoperability with `xarray` and chunked raster workflows.
- Better fault tolerance and task lifecycle controls than ad-hoc threads.
- Built-in observability (dashboard, task stream, memory profiling, timelines).
- Scales from local laptop to cluster without changing core module contracts.

Ray remains an optional future backend, but is not the default design target.

### 15.2 Execution Model

Keep module contracts unchanged; only scheduler/orchestration changes.

Planned DAG:

1. M1a: metadata + TOA load.
2. M1b: geometry and cloud classes from shared TOA (parallel branches).
3. M2 and M3 in parallel once M1 metadata/bounds are available.
4. M4 -> M5 -> M6 remain ordered, with optional distributed execution for heavy internals.

### 15.3 Robustness Requirements

- Configurable task retries per stage (default >= 2 for remote I/O stages).
- Per-stage timeout and cancellation propagation.
- Structured exception reporting with stage name + input context.
- Optional worker resource tags for memory-heavy tasks.

### 15.4 Progress and Visualization

Required user-facing progress options:

1. Dask dashboard URL (local default: `http://127.0.0.1:8787`).
2. CLI progress bars (`rich`) for stage-level status when dashboard is unavailable.
3. Optional `performance_report` HTML artifact per run.

Suggested report outputs:

- `output_dir/reports/dask-performance.html`
- `output_dir/reports/pipeline-summary.json`

### 15.5 Configuration Additions (planned)

Add execution config section (name tentative: `execution`):

- `backend: "thread" | "dask"` (default `"dask"` for new flows, `"thread"` fallback)
- `max_workers: int`
- `retries: int`
- `stage_timeout_s: int | None`
- `dashboard: bool`
- `dashboard_address: str | None`
- `performance_report_path: Path | None`
- `show_progress: bool`

### 15.6 Validation Plan for Redesign

Unit tests:

1. DAG ordering and dependency correctness.
2. Retry behavior for transient failures.
3. Cancellation propagation on hard failures.
4. Shared TOA reuse (no duplicate read) under Dask scheduling.

Integration tests:

1. End-to-end run with Dask local cluster.
2. Dashboard/report artifact generation.
3. Equivalence checks vs thread backend outputs.
