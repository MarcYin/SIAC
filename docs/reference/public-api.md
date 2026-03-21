# Public API Reference

This page documents the stable public entry points exposed by the package and CLI. Everything here should be treated as the supported user-facing surface; lower-level modules are implementation details unless explicitly documented otherwise.

## Python API

## `siac.SIAC`

User-facing session facade for local scene processing.

Typical use:

```python
from siac import SIAC
from siac.config import SIACConfig

config = SIACConfig(sensor="s2")
result = SIAC(config).process("/path/to/S2_SAFE/")
```

Main method:

- `process(input_path, output_path=None, *, aoi=None) -> CorrectionResult`

## `siac.SIACConfig`

Public configuration wrapper around the system schema plus lightweight Python-side defaults for `sensor` and `aoi`.

Useful constructors:

- `SIACConfig()`
- `SIACConfig.from_file(path)`
- `SIACConfig.load(path=None, **overrides)`

Useful methods:

- `with_overrides(...)`
- `with_env_overlay()`
- `resolve(RunRequest(...))`
- `to_file(path)`

## `siac.process_sentinel2`

Convenience helper for local scene processing with a default Sentinel-2 configuration.

Signature shape:

```python
process_sentinel2(input_path, output_path=None, **kwargs)
```

## `siac.process_landsat8`

Convenience helper for Landsat-8 processing with a default Landsat configuration.

Signature shape:

```python
process_landsat8(input_path, output_path=None, **kwargs)
```

## `siac.resolve_s2_input`

Resolve a Sentinel-2 query or path into a local path.

Signature shape:

```python
resolve_s2_input(query, config, *, auth=None)
```

Use this when you want to inspect or manage query resolution separately from full processing.

## `siac.search_sentinel2`

Search for Sentinel-2 products using the configured backend.

Signature shape:

```python
search_sentinel2(
    *,
    tile=None,
    date=None,
    date_value=None,
    start_date=None,
    end_date=None,
    bbox=None,
    max_cloud_cover=80.0,
    backend="cdse",
    config=None,
    auth=None,
)
```

Returns a list of typed product objects.

## `siac.siac_process_s2`

Primary query-based Sentinel-2 processing helper.

Signature shape:

```python
siac_process_s2(config, query, *, output_path=None, aoi=None, auth=None)
```

This is the clearest Python entry point when your input is not already a local SAFE path.

## CLI

## `siac process-s2`

The main CLI command currently exposed by the package.

Usage pattern:

```bash
siac process-s2 QUERY [--config path] [--output-path dir] [--aoi path | --aoi-bbox ...]
```

Accepted query forms:

- local SAFE path
- Sentinel-2 product ID
- tile/date shorthand such as `T31UDQ_20240101`

## Entry point summary

| Entry point | Audience | Input form | Return/output |
| --- | --- | --- | --- |
| `SIAC` | Python user | config plus local scene path | `CorrectionResult` |
| `SIACConfig` | Python user | defaults, TOML, overrides | typed config |
| `process_sentinel2` | Python user | local scene path | `CorrectionResult` |
| `process_landsat8` | Python user | local scene path | `CorrectionResult` |
| `resolve_s2_input` | Python user | query or path | local path |
| `search_sentinel2` | Python user | search fields | list of products |
| `siac_process_s2` | Python user | config plus query | `CorrectionResult` |
| `siac process-s2` | CLI user | query/path and options | disk outputs plus CLI status |

## Notes on compatibility

- The API surface above is the supported public layer.
- Lower-level packages such as `app`, `workflows`, `adapters`, `algorithms`, and `runtime` are documented for developers but should not be treated as stable end-user APIs.
