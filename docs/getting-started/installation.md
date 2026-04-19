# Installation

SIAC can be installed in a minimal Python environment, but the canonical
development and validation workflow uses the Pixi workspace in the repository
root.

## Choose an installation track

| Audience | Recommended path | Notes |
| --- | --- | --- |
| Runtime user | `pip install -e ".[docs]"` in a checked-out repo | Good if you mainly want to run SIAC and read docs locally |
| Scientific power user | Editable install plus Rust build | Best when testing algorithms and reproducing workflows |
| Developer | Pixi workspace plus Rust build | Matches the repo tasks used for lint, type checking, tests, coverage, and docs |

## Prerequisites

- Python 3.10 to 3.12
- Rust toolchain for the native extension
- GDAL/raster stack as provided through the package environment

## Pixi workspace

```bash
pixi install
pixi run bootstrap
pixi run build-rust
```

This gives you:

- the `siac` package
- development tools such as `pytest`, `ruff`, and `mypy`
- docs tooling including MkDocs Material and Mermaid rendering
- the Rust extension compiled into the active environment

The standard validation commands are:

```bash
pixi run build-rust
pixi run lint
pixi run typecheck-scoped
pixi run test-fast
pixi run test
pixi run coverage
```

The Pixi workspace is the easiest way to reproduce CI locally because it
installs the same general toolchain used for tests, coverage, and documentation
builds.

Run `pixi run build-rust` first in a fresh workspace so the test suite can
import `siac._rust` during collection.

`pixi run bootstrap` uses the Rust crate declared in
`src/siac_rs/Cargo.toml`. If you move the Rust crate, update the
`[tool.maturin]` `manifest-path` in `pyproject.toml` or editable installs will
stop working.

## Editable install

```bash
python -m pip install -e ".[dev,docs]"
maturin develop --release --manifest-path src/siac_rs/Cargo.toml
```

Use this only when you specifically want a non-Pixi environment.

## Optional extras

- `.[dev]` for development tools
- `.[docs]` for the documentation stack
- `.[gee]` for Google Earth Engine support
- `.[py6s]` for the Py6S radiative transfer backend

## Native 6S Backend Setup

The native `sixs` backend is separate from the optional Py6S dependency.

Use the native backend when you set:

```toml
[algorithms.rt]
backend = "sixs"
```

### Requirements

Native 6S requires:

- a Fortran compiler
- `meson` and `ninja` for the preferred Meson-backed F2PY path

The repository provides these through the Pixi `rt6s` environment/feature. If
Meson compilation fails on a platform, the SIAC builder falls back to the F2PY
`distutils` backend automatically.

### Recommended Path

```bash
pixi install
pixi run -e rt6s build-6s-native
```

This builds the compiled Python extension used by the native backend.

For the closest CI reproduction, run the native smoke script in the same
environment:

```bash
PYTHONPATH=python pixi run -e rt6s python tools/rt6s_smoke.py
```

### Manual Path

If you already have the compiler toolchain in a non-Pixi environment:

```bash
python tools/build_6s_native.py
```

`tools/build_6s_native.py` supports `--source-dir`, `--build-dir`,
`--compiler`, and `--build-profile`.

### Important Distinction

`pip install -e ".[py6s]"` enables the separate `py6s` backend only. It does
not install or build the native 6SV2.1 backend.

## Common installation issues

### Rust extension not built

If imports work but `siac._rust` is missing, rebuild the extension in the Pixi
environment:

```bash
pixi run build-rust
```

### Native geospatial stack problems

Prefer the Pixi environment when GDAL, rasterio, or JP2 support is difficult to install manually.

### Optional backend not available

Some workflows require extra dependencies or credentials. Install the relevant optional extras and configure the matching auth variables before running remote-data flows.

### Native 6S build fails

Check that the active environment has:

- a Fortran compiler on `PATH`
- `meson` and `ninja` if you want the preferred Meson path

If you want the repo-managed toolchain, switch to the Pixi `rt6s` environment
and rerun `pixi run -e rt6s build-6s-native`.
