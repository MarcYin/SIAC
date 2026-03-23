# SIAC

Sensor-Invariant Atmospheric Correction (SIAC) is a modular framework for atmospheric correction of medium-resolution optical satellite imagery. It is designed to turn top-of-atmosphere observations into bottom-of-atmosphere reflectance products with atmospheric diagnostics and uncertainty layers, while keeping the processing flow explicit enough for both operational and scientific use.

SIAC v2 in this repository is organized around:

- a small public Python API and CLI
- a typed TOML-based configuration system
- a staged runtime pipeline with explicit data contracts
- pluggable providers for satellite access, atmospheric priors, BRDF priors, and radiative transfer backends
- Rust-accelerated kernels, PSF, emulator, optimization, and Whittaker smoothing

The scientific reference for the method is:

- Yin, F., Lewis, P. E., and Gomez-Dans, J. L. (2022), [Bayesian atmospheric correction over land: Sentinel-2/MSI and Landsat 8/OLI](https://gmd.copernicus.org/articles/15/7933/2022/)

## Who This Is For

- **Basic users** who need a fast path to install SIAC and run a scene.
- **Scientific users** who need the theory, assumptions, and paper-to-code mapping.
- **Operators and engineers** who need configuration, credentials, cache, runtime, and output details.
- **Developers** who need package boundaries, runtime data flow, and extension points.

## Supported Workflows

- Sentinel-2 as the primary user-facing workflow
- Landsat entry points in the Python API
- Local SAFE input and remote Sentinel-2 resolution/search via configured backends
- Output writing as GeoTIFF, COG, NetCDF, or Zarr

## Core Outputs

- BOA reflectance
- BOA uncertainty, when enabled
- Aerosol optical thickness (AOT)
- Total column water vapour (TCWV)
- Cloud mask
- Quicklook and auxiliary artifacts, depending on output settings

## Fast Install

### Pixi workspace

```bash
pixi install
pixi run bootstrap
pixi run build-rust
```

Pixi is the canonical local setup for SIAC development and validation. It
matches the repository tasks used for linting, type checking, tests, coverage,
and the Rust-backed code paths.

Build the Rust extension before running the full test or coverage tasks in a
fresh environment; the suite imports `siac._rust` during collection.

### Python package

If you are not using Pixi, install the editable package and build the Rust
extension manually:

```bash
python -m pip install -e ".[dev,docs]"
maturin develop --release --manifest-path src/siac_rs/Cargo.toml
```

## Fast Run

### Canonical developer workflow

```bash
pixi run build-rust
pixi run lint
pixi run typecheck-scoped
pixi run test-fast
pixi run test
pixi run coverage
```

Run `pixi run build-rust` after bootstrapping a fresh environment, and rerun it
before the full test or coverage commands when you change the native extension.

### Python API

```python
from siac import SIAC
from siac.config import SIACConfig

config = SIACConfig(sensor="s2")
result = SIAC(config).process("/path/to/S2_SAFE/")

print(result.boa)
print(float(result.aot.mean()))
```

### CLI

For a direct scene run, the CLI entry point is:

```bash
siac process-s2 T31UDQ_20240101 --output-path ./outputs/example
```

## Runtime Flow

```mermaid
flowchart LR
    Input["User input or S2 query"] --> Resolve["Config load + resolve"]
    Resolve --> M1["M1: preprocess scene"]
    M1 --> M2["M2: atmospheric prior"]
    M1 --> M3["M3: surface prior"]
    M2 --> M4["M4: grid assembly"]
    M3 --> M4
    M4 --> M5["M5: solver"]
    M5 --> M6["M6: correction"]
    M6 --> Output["Outputs"]
```

## Documentation

- Hosted docs: [marcyin.github.io/SIAC](https://marcyin.github.io/SIAC/)
- Documentation portal source: [docs/index.md](docs/index.md)
- Installation guide: [docs/getting-started/installation.md](docs/getting-started/installation.md)
- First run: [docs/getting-started/first-run.md](docs/getting-started/first-run.md)
- Configuration basics: [docs/user-guide/configuration-basics.md](docs/user-guide/configuration-basics.md)
- Running SIAC: [docs/user-guide/running-siac.md](docs/user-guide/running-siac.md)
- Theory: [docs/science/theory.md](docs/science/theory.md)
- Architecture: [docs/architecture/package-map.md](docs/architecture/package-map.md)

## GitHub Pages

The repository is set up to build a MkDocs Material site in GitHub Actions and deploy it to GitHub Pages. Repository Pages must be configured to publish from GitHub Actions.

### Local docs workflow

```bash
pixi run docs-serve
pixi run docs-check
```
