# Installation

SIAC can be installed in a minimal Python environment, in the repository development layout, or through the Pixi workspace used by CI and local development.

## Choose an installation track

| Audience | Recommended path | Notes |
| --- | --- | --- |
| Runtime user | `pip install -e ".[docs]"` in a checked-out repo | Good if you mainly want to run SIAC and read docs locally |
| Scientific power user | Editable install plus Rust build | Best when testing algorithms and reproducing workflows |
| Developer | Pixi workspace plus Rust build | Closest match to CI and the repo's local tooling |

## Prerequisites

- Python 3.10 to 3.12
- Rust toolchain for the native extension
- GDAL/raster stack as provided through the package environment

## Editable install

```bash
python -m pip install -e ".[dev,docs]"
maturin develop --release --manifest-path src/siac_rs/Cargo.toml
```

This gives you:

- the `siac` package
- development tools such as `pytest`, `ruff`, and `mypy`
- docs tooling including Sphinx, MyST, and Mermaid rendering
- the Rust extension compiled into the active environment

## Pixi workspace

```bash
pixi install
pixi run build-rust
```

The Pixi workspace is the easiest way to reproduce CI locally because it installs the same general toolchain used for tests and documentation builds.

## Optional extras

- `.[dev]` for development tools
- `.[docs]` for the documentation stack
- `.[gee]` for Google Earth Engine support
- `.[py6s]` for the Py6S radiative transfer backend

## Common installation issues

### Rust extension not built

If imports work but `siac._rust` is missing, rebuild the extension:

```bash
maturin develop --release --manifest-path src/siac_rs/Cargo.toml
```

### Native geospatial stack problems

Prefer the Pixi environment when GDAL, rasterio, or JP2 support is difficult to install manually.

### Optional backend not available

Some workflows require extra dependencies or credentials. Install the relevant optional extras and configure the matching auth variables before running remote-data flows.
