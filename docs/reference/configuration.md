# Configuration Reference

This page documents the public configuration surface driven by `SIACConfig`, `SystemConfig`, and the typed schema in `python/siac/config/schema.py`.

## Precedence model

SIAC resolves configuration in this order:

| Step | Effect |
| --- | --- |
| Built-in defaults | Seed the config with safe schema defaults |
| TOML file | Apply persistent machine or project settings |
| Environment overlay | Fill missing secret/auth fields from environment variables |
| Per-run request fields | Apply `sensor`, `aoi`, `input_path`, `output_path`, and `s2_query` |
| Resolution step | Produce a final `ResolvedConfig` used by assembly and workflows |

```mermaid
flowchart TD
    Defaults["Schema defaults"] --> File["TOML config"]
    File --> Env["Environment secret overlay"]
    Env --> Run["Per-run request fields"]
    Run --> Resolved["ResolvedConfig"]
```

## Top-level sections

## `paths`

Controls shared filesystem and remote resource locations.

Key fields:

- `dem`
- `water_mask`
- `emulator_dir`
- `lut_path`
- `rsrf_root`
- `cache_root`
- `caches.atmo`
- `caches.brdf`
- `caches.s2`

Use `cache_root` when you want SIAC to derive per-subsystem cache directories automatically.

## `auth`

Credentials for:

- `cdse`
- `cds`
- `aws`
- `earthdata`
- `gcs`

Important environment variables:

- `SIAC_CDSE_USERNAME`
- `SIAC_CDSE_PASSWORD`
- `CDSAPI_KEY`
- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`
- `EARTHDATA_USERNAME`
- `EARTHDATA_PASSWORD`
- `GOOGLE_APPLICATION_CREDENTIALS`

## `providers`

### `providers.atmo`

Main fields:

- `kind`
- `data_path`
- `cache_dir`
- `download_missing`
- `temporal_interpolation`
- `user_aot`
- `user_tcwv`
- `user_tco3`

Supported kinds include `cams`, `merra2`, `mcd19`, `vnp19`, `era5`, and `user`.

### `providers.brdf`

Main fields:

- `kind`
- `data_path`
- `cache_dir`
- `temporal_window`
- `use_cache`
- `gee_project`
- `zarr_url`

Supported kinds include `mcd43`, `vnp43`, `mcd19`, `gee`, `zarr`, and `user`.

### `providers.s2`

Main fields:

- `backend`
- `cache_dir`
- `max_cloud_cover`
- `prefer_newest_baseline`
- `processing_level`

Supported backends:

- `local`
- `cdse`
- `gcs`

## `algorithms`

### `algorithms.surface_prior`

Main fields:

- `method`
- `psf_sigma_x`
- `psf_sigma_y`
- `apply_psf`
- `whittaker_lambda`
- `monthly_database_resolution_policy` — for `method = "monthly_database"`, choose `provider_or_coarser` to build the Route-B database at the prepared/provider grid when it is coarser than the AOT grid, or `aerosol` to force the Route-B database/query grid to the configured solver `aerosol_resolution`
- `spectral_mapping.*`

Supported methods:

- `kernel_model`
- `whittaker`
- `monthly_database`
- `neural`
- `direct`

### `algorithms.rt`

Main fields:

- `backend`
- `lut_interpolation`
- `lut_storage_options`
- `py6s_aero_profile`
- `sixs.*`
- `fallback_to_lut`
- `fallback_to_py6s`

Supported backends:

- `emulator`
- `lut`
- `py6s`
- `sixs`

`py6s` and `sixs` are separate backends. `.[py6s]` installs the optional Py6S
dependency for the `py6s` backend only; the native `sixs` backend uses the
compiled 6SV2.1 build path documented in
[Native 6SV2.1 Backend](../user-guide/native-sixs-backend.md).

### `algorithms.rt.sixs`

Main execution and build fields:

- `source_url`
- `source_dir`
- `build_dir`
- `module_path`
- `library_path` — compatibility alias for native module override
- `auto_build`
- `compiler`
- `build_profile` — `release` or `parity`
- `mode` — `direct`, `scene_lut`, or `auto`
- `parallel_backend` — `openmp` or `worker_libraries`
- `native_threads`
- `worker_libraries`
- `chunk_size`
- `scene_lut_min_pixels`
- `scene_lut_max_nodes_per_axis`
- `scene_lut_max_cases`
- `scene_lut_required_speedup`
- `month`
- `day`
- `reference_reflectance`
- `output_variables`

Atmospheric-profile fields:

- `atmospheric_profile`
- `atmospheric_profile_latitude`
- `radiosonde_profile.*`

Supported atmospheric profiles:

- `no_gas`
- `tropical`
- `midlatitude_summer`
- `midlatitude_winter`
- `subarctic_summer`
- `subarctic_winter`
- `us_standard_62`
- `auto_latitude_date`
- `user_water_ozone`
- `user_profile`

Aerosol-profile fields:

- `aerosol_profile`
- `aerosol_mixture`
- `aerosol_distribution.*`
- `sun_photometer_aerosol.*`
- `aerosol_layer_profile.*`
- `aerosol_model_path`

Supported aerosol profiles:

- `none`
- `continental`
- `maritime`
- `urban`
- `desert`
- `biomass_burning`
- `stratospheric`
- `user_mixture`
- `multimodal_log_normal`
- `modified_gamma`
- `junge_power_law`
- `sun_photometer`
- `layered_profile`
- `user_model`

Surface and correction fields:

- `surface.mode`
- `surface.target.*`
- `surface.environment.*`
- `surface.radius_km`
- `surface.brdf.*`
- `atmospheric_correction.mode`
- `atmospheric_correction.value`

Supported surface modes:

- `homogeneous_lambertian`
- `heterogeneous_lambertian`
- `homogeneous_brdf`

Supported surface target kinds:

- `constant`
- `built_in`
- `spectrum`

Supported built-in target/environment reflectances:

- `green_vegetation`
- `clear_water`
- `sand`
- `lake_water`

Supported `surface.brdf.model` values:

- `user_defined`
- `hapke`
- `verstraete`
- `roujean`
- `walthall`
- `minnaert`
- `ocean`
- `iaquinta_pinty`
- `rahman`
- `kuusk`
- `modis`
- `ross_li_maignan`

Supported atmospheric-correction modes:

- `none`
- `lambertian_reflectance`
- `lambertian_radiance`
- `brdf_reflectance`
- `brdf_radiance`

Important validation rules:

- `atmospheric_profile = "auto_latitude_date"` requires
  `atmospheric_profile_latitude`
- `atmospheric_profile = "user_profile"` requires `radiosonde_profile`
- `aerosol_profile = "user_mixture"` requires `aerosol_mixture`, which must sum
  to `1.0`
- `aerosol_profile = "user_model"` requires `aerosol_model_path`
- BRDF atmospheric-correction modes require `surface.mode = "homogeneous_brdf"`
- `output_variables` accepts the native 6S output names documented in
  [Native 6SV2.1 Backend](../user-guide/native-sixs-backend.md)

For examples and route-selection guidance, see
[Native 6SV2.1 Backend](../user-guide/native-sixs-backend.md).

### `algorithms.solver`

Main fields:

- `aot_gamma` — regularization strength for AOT prior term (default `0.05`)
- `tcwv_gamma` — regularization strength for TCWV prior term (default `0.05`)
- `alpha` — band weight power exponent for the observation cost (default `-1.6`). Controls wavelength-dependent weighting: bands are weighted proportionally to $\lambda^{\alpha}$, so negative values give more weight to shorter (aerosol-sensitive) wavelengths. Propagated internally as `MultiGridConfig.band_weight_power` and `CostFunctionConfig.band_weight_power`.
- `max_iterations` — maximum L-BFGS-B iterations per grid level
- `gtol` — gradient tolerance for L-BFGS-B convergence
- `ftol` — function tolerance for L-BFGS-B convergence
- `aerosol_resolution` — target spatial resolution (metres) for the aerosol retrieval grid
- `grid_search_aot_points` — number of AOT candidates in the no-Jacobian grid-search path (default `11`)
- `grid_search_tcwv_points` — number of TCWV candidates in the no-Jacobian grid-search path (default `11`)
- `fixed_atmospheric_parameter` — set to `aot` or `tcwv` to hold that field at the atmospheric prior while solving only the other field (default `none`)
- `stages` — optional staged solver chain. Each stage declares `solve`, `fixed`, optional `bands`, and `initial_state = "previous"` or `"prior"`. Current production execution supports staged AOT/TCWV combinations; `tco3` is carried through the atmospheric state and rejected as a solved parameter until RT ozone Jacobian/grid-search support is added.
- `quadratic_block_size` — solve one shared AOT/TCWV pair for each `NxN` block in the no-Jacobian grid-search path, compute RT coefficients on the same block grid, then broadcast the block solution back to full resolution
- `quadratic_block_min_valid_fraction` — minimum fraction of pixels in each `quadratic_block_size` block that must have valid observations and surface-prior support before the block is solved (default `0.5`)
- `water_mask_buffer_pixels` — dilate the native water mask by this many native mask pixels before it is reprojected to the aerosol solver grid
- `sharp_transition_filter.*` — optional native-resolution edge/singularity detector used to exclude unstable transitions before the solver-grid aggregation step
- `use_multigrid` — enable the coarse-to-fine multi-grid solver strategy (default `true`)
- `min_grid_size` — minimum grid dimension (pixels) for multi-grid levels
- `bounds.aot` — `[min, max]` bounds for AOT during optimization
- `bounds.tcwv` — `[min, max]` bounds for TCWV during optimization

Example staged solve:

```toml
[[algorithms.solver.stages]]
name = "aot_pass"
solve = ["aot"]
fixed = ["tcwv", "tco3"]
bands = ["B01", "B02", "B04"]

[[algorithms.solver.stages]]
name = "tcwv_pass"
solve = ["tcwv"]
fixed = ["aot", "tco3"]
initial_state = "previous"
```

### `algorithms.cloud_mask`

Main fields:

- `mode`
- `provider`
- `external_mask_path`
- `class_mapping`
- `unmapped_to_missing`
- `target_resolution_m`
- `resolution_policy`
- `allow_upsample_to_target`

## `runtime.execution`

Execution knobs:

- `backend`
- `max_workers`
- `correction_max_workers` — Optional worker count for M6 atmospheric correction band-level parallelism. If omitted, uses `max_workers`.
- `retries`
- `stage_timeout_s`
- `dashboard`
- `dashboard_address`
- `performance_report_path`
- `show_progress`

Supported execution backends:

- `thread`
- `dask`

## `output.defaults`

Controls output persistence.

Main fields:

- `output_dir`
- `format`
- `compression`
- `include_rgb`
- `include_uncertainty`
- `include_auxiliary`
- `skip_correction` — When `true`, the aerosol retrieval solver (M5) still runs but the atmospheric correction (M6) is skipped, so no BOA reflectance is produced. The retrieved AOT, TCWV, solver QA, cloud mask, surface prior, and monthly composites are still written. Default: `false`.
- `boa_dtype`
- `boa_scale`
- `boa_nodata`
- `reopen_streamed_boa` — When `true` (default), BOA bands streamed during M6 are reopened from disk after write so downstream steps consume on-disk values. Set `false` for faster runs when write/read parity is not required.

Supported formats:

- `geotiff`
- `cog`
- `zarr`
- `netcdf`

## Minimal example

```toml
[providers.s2]
backend = "cdse"
max_cloud_cover = 40.0

[algorithms.rt]
backend = "emulator"

[output.defaults]
format = "cog"
include_uncertainty = true
```

Native 6S example:

```toml
[algorithms.rt]
backend = "sixs"

[algorithms.rt.sixs]
mode = "auto"
parallel_backend = "openmp"
build_profile = "release"
output_variables = ["xap", "xbp", "xcp", "tgasm", "sutott", "sast"]
```

## Practical examples

| Use case | Important settings |
| --- | --- |
| local-only run | `providers.s2.backend = "local"` |
| CDSE-backed Sentinel-2 | `providers.s2.backend = "cdse"` plus CDSE credentials |
| GCS-backed Sentinel-2 | `providers.s2.backend = "gcs"` plus cloud access setup |
| emulator RT | `algorithms.rt.backend = "emulator"` |
| LUT RT | `algorithms.rt.backend = "lut"` plus `paths.lut_path` |
| native 6S RT | `algorithms.rt.backend = "sixs"` plus `algorithms.rt.sixs.*` |
| Whittaker surface prior | `algorithms.surface_prior.method = "whittaker"` |

## Public configuration APIs

- `SIACConfig.from_file(path)`
- `SIACConfig.load(path=None, **overrides)`
- `SIACConfig.with_overrides(...)`
- `SIACConfig.with_env_overlay()`
- `SIACConfig.resolve(RunRequest(...))`

For user-facing narrative guidance, read [Configuration Basics](../user-guide/configuration-basics.md).
