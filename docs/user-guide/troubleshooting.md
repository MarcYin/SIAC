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

The native backend requires a Fortran compiler and a supported F2PY backend.
The repo-managed `rt6s` environment pins Python 3.11 and `setuptools < 60`,
which keeps the legacy `distutils` backend available for this extension while
still preferring the Meson backend first. If the native build fails in a
general-purpose environment, switch to `rt6s` and rerun:

```bash
pixi run -e rt6s build-6s-native
```

To mirror the hosted Linux smoke workflow exactly:

```bash
SIAC_SIXS_F2PY_BACKEND=meson,distutils PYTHONPATH=python pixi run -e rt6s python tools/rt6s_smoke.py
```

## Runtime and configuration

### Wrong cache or output location

Check:

- `paths.cache_root`
- per-provider cache directories
- `output.defaults.output_dir`
- per-run `output_path` overrides

### Dask backend fails

If you configure `execution.backend = "dask"` without the required runtime
support, switch back to `thread` or install `siac[dask]` in your environment.

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

### Remote ZIP/Zarr LUT rejects the requested RT setup

That is expected when `algorithms.rt.backend = "lut"` and the config tries to
customize atmosphere, aerosol, surface, or correction semantics. The packaged
remote LUT is a fixed libRadtran preset:

- atmosphere profile: `us_standard_62`
- aerosol family: `continental_average`
- surface mode: `homogeneous_lambertian`

Use `backend = "sixs"` if you need configurable atmosphere profiles, aerosol
families, BRDFs, or surface reflectance targets.

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

- `setup.atmosphere.profile = "auto_latitude_date"` without
  `setup.atmosphere.profile_latitude`
- `setup.atmosphere.profile = "user_profile"` without
  `setup.atmosphere.radiosonde_profile`
- `setup.aerosol.profile = "user_mixture"` without `setup.aerosol.mixture`
- `setup.aerosol.profile = "user_model"` without `setup.aerosol.model_path`
- BRDF atmospheric-correction mode without
  `setup.surface.mode = "homogeneous_brdf"`

### Unexpected remote fetch behavior

Remember that Sentinel-2 processing may first resolve a query through a backend and local cache before the normal scene-processing pipeline starts.
