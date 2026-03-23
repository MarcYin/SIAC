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

### Unexpected remote fetch behavior

Remember that Sentinel-2 processing may first resolve a query through a backend and local cache before the normal scene-processing pipeline starts.
