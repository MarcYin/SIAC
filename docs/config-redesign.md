# SIAC Config Redesign

## Goal

Replace the current mixed config implementation with a simpler model built around:

1. one canonical project or machine config file in TOML
2. one explicit per-run request object
3. one resolved execution plan consumed by the pipeline

This redesign intentionally drops the current dual-schema compatibility layer.

## Design Principles

- TOML is the only human-edited config file format.
- The package config stores stable settings, not scene-specific inputs.
- AOI, input paths, tile IDs, product IDs, and output paths are per-run inputs.
- Environment variables are used only as an overlay for secrets or machine overrides.
- The pipeline does not read raw config files directly.
- The pipeline consumes a resolved execution plan with no fallback logic left inside it.
- Serialization, file I/O, env overlay, resolution, and runtime snapshotting are separate concerns.

## Config Layers

### 1. `SystemConfig`

Persistent package or project settings.

Contents:

- paths to assets and caches
- provider defaults
- algorithm defaults
- runtime defaults
- output defaults
- auth references or auth values

This is the only thing loaded from `config.toml`.

### 2. `RunRequest`

Per-run inputs.

Contents:

- input path
- sensor or sensor hint
- AOI
- output path
- S2 query inputs like tile and date
- optional one-off overrides

This should be built by the CLI, Python API, or higher-level workflow code.

### 3. `ResolvedConfig`

Fully merged runtime state.

Contents:

- effective provider selection
- effective asset paths
- effective auth values
- effective algorithm settings
- effective runtime settings
- effective run inputs

This is the only config object the pipeline and resolver functions should consume.

## Canonical TOML Schema

```toml
[paths]
dem = "/data/siac/dem/copernicus_dem.vrt"
water_mask = "/data/siac/water/global_water.tif"
emulator_dir = "/data/siac/emulators"
lut_path = "s3://bucket/siac/lut.zarr"
rsrf_root = "/data/siac/RSRF"
cache_root = "/data/siac/cache"

[paths.caches]
atmo = "/data/siac/cache/atmo"
brdf = "/data/siac/cache/brdf"
s2 = "/data/siac/cache/s2"

[auth.cdse]
username = ""
password = ""
username_env = "SIAC_CDSE_USERNAME"
password_env = "SIAC_CDSE_PASSWORD"

[auth.earthdata]
username = ""
password = ""
username_env = "EARTHDATA_USERNAME"
password_env = "EARTHDATA_PASSWORD"

[auth.cds]
api_key = ""
api_key_env = "CDSAPI_KEY"

[auth.aws]
access_key_id = ""
secret_access_key = ""
access_key_id_env = "AWS_ACCESS_KEY_ID"
secret_access_key_env = "AWS_SECRET_ACCESS_KEY"

[providers.atmo]
kind = "cams"
download_missing = true
temporal_interpolation = "nearest"
data_path = "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/"

[providers.brdf]
kind = "mcd43"
temporal_window = 16
use_cache = true

[providers.s2]
backend = "cdse"
processing_level = "L1C"
max_cloud_cover = 80.0
prefer_newest_baseline = true

[algorithms.surface_prior]
method = "whittaker"
psf_sigma_x = 29.75
psf_sigma_y = 39.0
apply_psf = true
whittaker_lambda = 10.0

[algorithms.surface_prior.spectral_mapping]
enabled = true
k_neighbors = 5
neighbor_estimator = "distance_weighted_mean"
knn_backend = "scipy_ckdtree"
knn_eps = 0.0
min_valid_bands = 1

[algorithms.rt]
backend = "emulator"
lut_interpolation = "linear"
fallback_to_lut = true
fallback_to_py6s = true

[algorithms.solver]
aot_gamma = 10.0
tcwv_gamma = 5.0
alpha = -1.6
max_iterations = 300
gtol = 0.01
ftol = 1e-7
aerosol_resolution = 1000.0
use_multigrid = true
min_grid_size = 4

[algorithms.solver.bounds]
aot = [0.001, 2.5]
tcwv = [0.0, 7.0]

[algorithms.cloud_mask]
mode = "auto"
provider = "omnicloudmask"
target_resolution_m = 10.0
resolution_policy = "auto"
allow_upsample_to_target = false

[runtime]
log_level = "INFO"
n_jobs = -1
chunk_size = 2048

[runtime.execution]
backend = "thread"
max_workers = 4
retries = 2
stage_timeout_s = 300.0
dashboard = false
show_progress = false

[output.defaults]
format = "cog"
compression = "deflate"
include_rgb = true
include_uncertainty = true
include_auxiliary = true
boa_dtype = "float32"
boa_scale = 10000.0
boa_nodata = 0.0
```

## Explicitly Not In `config.toml`

These belong in `RunRequest`, not `SystemConfig`:

- `aoi`
- `input_path`
- `output_path`
- `tile`
- `product_id`
- `date`
- `start_date`
- `end_date`
- one-off query filters

Example Python usage:

```python
system = load_system_config("config.toml")
request = RunRequest(
    input_path="/path/to/S2_SAFE",
    sensor="s2",
    aoi="/path/to/aoi.geojson",
    output_path="/tmp/out",
)
resolved = resolve_config(system, request)
result = run_siac(resolved)
```

## Proposed Python Layout

```text
python/siac/config/
  __init__.py
  schema.py
  load.py
  resolve.py
  snapshot.py
```

### `schema.py`

Owns only typed models.

Recommended models:

- `SystemConfig`
- `PathsConfig`
- `CachePathsConfig`
- `AuthConfig`
- `EarthdataAuthConfig`
- `CDSAuthConfig`
- `AWSAuthConfig`
- `ProviderConfig`
- `AtmoProviderConfig`
- `BRDFProviderConfig`
- `S2ProviderConfig`
- `AlgorithmConfig`
- `SurfacePriorAlgorithmConfig`
- `SpectralMappingConfig`
- `RTAlgorithmConfig`
- `SolverAlgorithmConfig`
- `SolverBoundsConfig`
- `CloudMaskAlgorithmConfig`
- `RuntimeConfig`
- `ExecutionRuntimeConfig`
- `OutputDefaultsConfig`
- `RunRequest`
- `ResolvedConfig`

Rules:

- These are plain `BaseModel` classes.
- No file loading.
- No env reading.
- No redaction logic.
- No path existence probing.
- No backward-compat flattening.

### `load.py`

Owns file and env loading.

Recommended public functions:

- `load_system_config(path: Path | str) -> SystemConfig`
- `load_system_config_from_default() -> SystemConfig`
- `overlay_env_secrets(config: SystemConfig, env: Mapping[str, str] | None = None) -> SystemConfig`
- `write_system_config(config: SystemConfig, path: Path | str) -> None`
- `write_default_system_config(path: Path | str) -> None`

Rules:

- TOML only
- no YAML support
- env overlay is explicit
- env precedence is documented and deterministic

### `resolve.py`

Owns merging and fallback rules.

Recommended public functions:

- `resolve_config(system: SystemConfig, request: RunRequest) -> ResolvedConfig`
- `resolve_auth(system: SystemConfig) -> ResolvedAuthConfig`
- `resolve_paths(system: SystemConfig) -> ResolvedPathsConfig`
- `resolve_provider_plan(system: SystemConfig, request: RunRequest) -> ResolvedProviderConfig`

Rules:

- all fallback behavior lives here
- all defaults are frozen into `ResolvedConfig`
- all validation that depends on multiple sections lives here

Examples:

- if RT backend is `lut`, ensure `lut_path` is resolved
- if spectral mapping is enabled and needed, ensure spectral library exists
- if S2 backend is `cdse`, ensure credentials are resolved
- if AOI is missing, leave it unset and allow the preprocessor to infer it

### `snapshot.py`

Owns audit output.

Recommended public functions:

- `snapshot_system_config(config: SystemConfig, redact_secrets: bool = True) -> dict`
- `snapshot_resolved_config(config: ResolvedConfig, redact_secrets: bool = True) -> dict`
- `write_runtime_snapshot(config: ResolvedConfig, path: Path | str, redact_secrets: bool = True) -> None`

Rules:

- snapshots are diagnostics, not loadable config
- snapshots should default to JSON or YAML
- snapshots should show source path, env keys used, resolved paths, and redacted secrets

## Pipeline Integration

## Current Problem

`python/siac/siac.py` currently mixes:

- config loading
- request handling
- AOI coercion
- provider construction
- fallback logic
- runtime resolution

and duplicates resolver logic between class methods and module-level functions.

## Target

Move to a single construction path:

```python
system = load_system_config(path)
system = overlay_env_secrets(system)
request = RunRequest(...)
resolved = resolve_config(system, request)
pipeline = build_pipeline(resolved)
result = pipeline.run()
```

Recommended new modules:

```text
python/siac/runtime/
  plan.py
  builders.py
  run.py
```

### `runtime/plan.py`

- functions that turn `ResolvedConfig` into a small execution plan

### `runtime/builders.py`

- `build_preprocessor(plan)`
- `build_atmo_provider(plan)`
- `build_surface_prior_provider(plan)`
- `build_rt_model(plan)`
- `build_solver(plan)`
- `build_corrector(plan)`

### `runtime/run.py`

- `run_siac(plan) -> CorrectionResult`

Rules:

- no builder should read raw `SystemConfig`
- no builder should read TOML
- no builder should inspect env vars

## API Shape

### Keep

- `SIAC.from_config(path)`
- `SIAC.process(...)`
- `siac_process(...)`
- `siac_process_s2(...)`

### Change

- `SIAC.from_config(path)` should load `SystemConfig`, not a general-purpose mutable config blob
- `SIAC.process(...)` should build a `RunRequest` and pass it through the resolver
- `siac_process(...)` and class-based execution should share the same single implementation path

### Remove

- old flat config fields such as top-level `global_dem`, `global_water`, `log_level`, `n_jobs`, `chunk_size`
- flattening and synchronization logic
- YAML file support
- config-state files being treated as loadable config

## Migration Sequence

### Phase 1. Add new package without wiring it in

Create:

- `python/siac/config/schema.py`
- `python/siac/config/load.py`
- `python/siac/config/resolve.py`
- `python/siac/config/snapshot.py`

Add tests for:

- TOML read and write
- env overlay
- resolved config validation
- snapshot redaction

### Phase 2. Define the new canonical TOML and examples

Add:

- `docs/config.example.toml`
- `docs/config-redesign.md`

Decide the final section names and stop changing them after this phase.

### Phase 3. Introduce `RunRequest`

Move all per-run fields out of the persistent config path:

- AOI
- input path
- output path
- S2 query values

Update:

- `SIAC.process(...)`
- `siac_process(...)`
- `siac_process_s2(...)`

so they all build a `RunRequest`.

### Phase 4. Introduce `ResolvedConfig`

Move all merge and fallback logic out of the config model and into `resolve.py`.

Remove from the config model:

- path fallback methods
- synchronization validators
- redaction logic
- source-path tracking

### Phase 5. Unify runtime builders

Replace duplicated resolver and builder logic in `python/siac/siac.py` with a single runtime plan path.

No class path and module path divergence should remain after this phase.

### Phase 6. Switch public API to the new system

Update public imports and docs so all normal usage goes through:

- TOML system config
- per-run request arguments
- resolved execution plan

### Phase 7. Delete old config implementation

Remove:

- legacy flat fields
- YAML loaders
- categorized flattening logic
- handwritten config synchronization

This should be a hard cut, not a soft compatibility layer.

## Recommended Immediate Refactor Order

If this starts now, the first implementation sequence should be:

1. create new config package modules
2. move schema classes into `schema.py`
3. create `RunRequest` and `ResolvedConfig`
4. implement TOML-only loader
5. implement explicit resolver
6. make `siac_process(...)` the primary execution path
7. make `SIAC.process(...)` a thin wrapper over that path
8. delete duplicated builder logic from `python/siac/siac.py`
9. remove old flat config fields
10. remove YAML support

## Success Criteria

The redesign is complete when:

- there is one canonical config file format: TOML
- there is one canonical persistent schema: `SystemConfig`
- there is one canonical per-run schema: `RunRequest`
- there is one canonical resolved runtime schema: `ResolvedConfig`
- `siac.py` contains no config parsing logic
- the pipeline reads only `ResolvedConfig`
- snapshots are diagnostic only and cannot be mistaken for config input
- there is no config field duplication or synchronization logic left
