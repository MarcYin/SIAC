# Troubleshooting

## Installation and environment

### `siac._rust` cannot be imported

The Rust extension has not been built in the active environment. Rebuild it in
the Pixi workspace with:

```bash
pixi run build-rust
```

### Geospatial dependencies are difficult to install

Prefer the Pixi environment. It is the easiest way to match the repo's tested
setup for SIAC development, validation, and docs.

## Credentials and backends

### Sentinel-2 query fails with a local backend error

If the input does not exist as a local SAFE path and `providers.s2.backend` is `local`, SIAC cannot resolve the query remotely. Either provide an existing path or switch to a remote backend such as `cdse` or `gcs`.

### Remote data access fails

Check the matching credentials for your chosen providers and backend:

- CDSE credentials for CDSE-backed scene access
- Earthdata credentials for Earthdata-backed providers
- CDS API key for CDS-backed atmospheric data
- cloud credentials for AWS or GCS-backed access

### Native `sixs` backend cannot build

The native backend requires a Fortran compiler, `meson`, and `ninja`. If those
tools are missing, either install them directly or use the Pixi `rt6s`
environment and rerun:

```bash
pixi run -e rt6s build-6s-native
```

### Installed `.[py6s]` but native `sixs` still is not available

`.[py6s]` only enables the separate `py6s` backend. The native `sixs` backend
still needs the native toolchain and compiled extension.

## Runtime and configuration

### Wrong cache or output location

Check:

- `paths.cache_root`
- per-provider cache directories
- `output.defaults.output_dir`
- per-run `output_path` overrides

### Dask backend fails

If you configure `execution.backend = "dask"` without the required runtime support, switch back to `thread` or install the missing Dask dependencies in your environment.

## Data and workflow behavior

### AOI is rejected

Valid AOI forms include:

- GeoJSON path
- dict-like GeoJSON object
- bounds tuple or list with four values

### LUT or emulator issues

Check:

- `paths.lut_path`
- `paths.emulator_dir`
- `algorithms.rt.backend`
- backend fallbacks such as `fallback_to_lut`

### Native `sixs` route is unexpectedly slow

If `algorithms.rt.sixs.mode = "scene_lut"` is slower than expected, the scene
likely has too much geometry or atmospheric variability to compress well into a
reduced interpolation grid. Prefer:

- `mode = "auto"` for normal runs
- `mode = "direct"` for small or geometry-rich scenes

Relevant tuning fields:

- `scene_lut_min_pixels`
- `scene_lut_max_nodes_per_axis`
- `scene_lut_max_cases`
- `scene_lut_required_speedup`

### OpenMP-backed native execution is unstable in the target environment

Switch the native backend to the fallback worker-library path:

```toml
[algorithms.rt]
backend = "sixs"

[algorithms.rt.sixs]
parallel_backend = "worker_libraries"
worker_libraries = 4
chunk_size = 4096
```

### Native `sixs` config is rejected

Common validation failures:

- `atmospheric_profile = "auto_latitude_date"` without
  `atmospheric_profile_latitude`
- `atmospheric_profile = "user_profile"` without `radiosonde_profile`
- `aerosol_profile = "user_mixture"` without `aerosol_mixture`
- `aerosol_profile = "user_model"` without `aerosol_model_path`
- BRDF atmospheric-correction mode without `surface.mode = "homogeneous_brdf"`

### Unexpected remote fetch behavior

Remember that Sentinel-2 processing may first resolve a query through a backend and local cache before the normal scene-processing pipeline starts.
