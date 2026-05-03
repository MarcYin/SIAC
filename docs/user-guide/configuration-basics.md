# Configuration Basics

SIAC uses a typed TOML configuration model backed by `SIACConfig`, `SystemConfig`, and a runtime resolution step that merges persistent settings with per-run inputs.

## Configuration sources

| Source | Role |
| --- | --- |
| Built-in defaults | Safe baseline values in the schema |
| TOML file | Persistent project or machine settings |
| Environment overlay | Secrets and credential fields |
| Per-run request inputs | Input path, output path, AOI, sensor, and S2 query |
| Resolved runtime config | Final merged configuration used by the pipeline |

For scene processing, the per-run `sensor` value is limited to `auto` or `s2`.
Landsat sensor metadata exists in the catalog for spectral and RT components,
but Landsat scenes are not accepted by the end-to-end preprocessing pipeline.

## Default config location

SIAC looks for configuration in the default path exported by `siac.config.DEFAULT_CONFIG_PATH`, unless overridden by the `SIAC_CONFIG_FILE` environment variable.

## Main sections

### `paths`

Controls global paths such as:

- DEM
- emulator directory
- LUT path
- spectral library root
- RSRF root
- cache root and per-subsystem caches

### `auth`

Holds credentials for:

- CDSE
- CDS
- AWS
- Earthdata
- GCS

In most cases, secrets are better supplied through environment variables than committed to TOML.

### `providers`

Chooses and configures the data providers used for:

- atmospheric priors
- BRDF priors
- Sentinel-2 scene access

### `algorithms`

Controls algorithm choices and parameters for:

- surface priors
- radiative transfer backend
- solver behavior
- cloud masking

`algorithms.rt.backend` supports `emulator`, `lut`, and the native `sixs`
backend. The native backend has its own nested
`algorithms.rt.sixs.*` configuration surface for build, routing, and output
selection. The shared RT semantic layer lives under `algorithms.rt.setup.*`;
backend-specific blocks such as `algorithms.rt.sixs.*` now carry native
execution/build controls only. See
[Native 6SV2.1 Backend](native-sixs-backend.md).
The packaged remote `lut` route is a fixed libRadtran preset within that shared
framework, not a fully configurable RT family.

The native `sixs` backend can also be used outside the full end-to-end SIAC
pipeline. If you already have `GeometryAngles`, `AtmosphericState`, and TOA
reflectance, you can start directly from the RT or correction stage with the
same `algorithms.rt.setup.*` semantics plus the native `algorithms.rt.sixs.*`
runtime controls. That direct-start workflow is documented in
[Native 6SV2.1 Backend](native-sixs-backend.md).

### `runtime`

Controls execution options such as log level, chunk size, and the thread or Dask backend.

### `output`

Controls output directory, format, compression, uncertainty inclusion, auxiliary outputs, and BOA encoding details.

Execution tuning includes `runtime.execution.max_workers` and optional
`runtime.execution.correction_max_workers` for M6 band-level correction parallelism.

## Environment overlay

The config loader fills missing auth fields from environment variables. Common examples include:

- `SIAC_CDSE_USERNAME`
- `SIAC_CDSE_PASSWORD`
- `EARTHDATA_USERNAME`
- `EARTHDATA_PASSWORD`
- `CDSAPI_KEY`
- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`
- `GOOGLE_APPLICATION_CREDENTIALS`

## Practical rule of thumb

- Put stable machine or project settings in TOML.
- Put secrets in environment variables.
- Put scene-specific values such as input path, AOI, query, and output path in the actual API or CLI call.

For the full field-by-field reference, continue to
[Configuration Reference](../reference/configuration.md). For native 6S setup
and route-selection guidance, continue to
[Native 6SV2.1 Backend](native-sixs-backend.md).
