# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Branch Context

Active development happens on `feat/refactor`, which contains SIAC v2 — a complete rewrite (Python package in `python/siac/` + Rust extension in `src/siac_rs/`). The `master` branch holds the legacy v1 package (flat `SIAC/` directory) and is not where this documentation applies. The leftover `SIAC/` and `SIAC.egg-info/` directories in the working tree are untracked v1 artifacts, not part of this branch.

## What SIAC Is

Sensor-Invariant Atmospheric Correction: converts top-of-atmosphere (TOA) satellite reflectance to bottom-of-atmosphere (BOA) surface reflectance with uncertainties, via Bayesian inversion of atmospheric state (AOT, TCWV, TCO3) using BRDF-derived surface priors and CAMS/MERRA-2 atmospheric priors. Sentinel-2 is the end-to-end supported workflow; Landsat exists at the catalog/spectral level only. Reference paper: Yin et al. 2022, https://gmd.copernicus.org/articles/15/7933/2022/

## Commands

Pixi is the canonical environment (`pixi.toml` defines all tasks; CI runs these exact tasks):

```bash
pixi install                 # create environment
pixi run bootstrap           # pip install -e . --no-deps
pixi run build-rust          # maturin build of siac._rust — REQUIRED before tests
```

The test suite imports `siac._rust` during collection, so `build-rust` must run in a fresh environment, and must be re-run after any change under `src/siac_rs/`.

```bash
pixi run test                # full suite with coverage
pixi run test-no-cov         # full suite, no coverage
pixi run test-fast           # excludes slow/integration/regression markers
pixi run coverage            # tests + enforce ≥95% total / ≥90% per-file
pixi run lint                # ruff check python/siac tests
pixi run format              # ruff format (format-check for CI mode)
pixi run typecheck-scoped    # mypy on the curated file list — this is the CI gate
pixi run rust-test           # cargo test
pixi run rust-clippy         # cargo clippy -D warnings (CI gate)
pixi run docs-build          # mkdocs build --strict
```

Single test (pytest `addopts` enables coverage by default — pass `--no-cov` for quick runs):

```bash
pixi run pytest tests/unit/test_validation.py --no-cov
pixi run pytest tests/unit/test_validation.py::test_name --no-cov
```

Optional native RT backends are built in dedicated pixi environments:

```bash
pixi run -e rt6s build-6s-native        # compiles 6SV2.1 (Python 3.11 + setuptools<60 pins)
pixi run -e libradtran build-libradtran # compiles libRadtran/uvspec
```

CI (`.github/workflows/ci.yml`) runs: build-rust, lint, format-check, typecheck-scoped, rust-test, clippy, unit tests, an integration slice (`test_e2e_synthetic.py`, `test_pipeline.py`, `test_injection.py`), import/wheel smoke checks. Changes under `python/siac/algorithms/rt/direct/`, `python/siac/config/`, or `python/siac/sixs_outputs.py` should also be validated against the `rt6s-smoke.yml` workflow.

Credentials come from env vars: `SIAC_CDSE_USERNAME`/`SIAC_CDSE_PASSWORD` (Copernicus Dataspace) and `SIAC_EARTHDATA_USERNAME`/`SIAC_EARTHDATA_PASSWORD` (NASA Earthdata). Live-network integration tests need these; unit tests do not.

## Architecture

Layered package under `python/siac/` (note `pythonpath = ["python"]`; maturin `python-source = "python"`). Authoritative docs: `docs/code_structure.md`, `docs/architecture/execution-flow.md`, `docs/naming-conventions.md`, `docs/developer/dev-setup.md`.

### Layers (keep responsibilities separate)

- `api/` — canonical public surface: `SIAC`, `siac_process_s2(...)`, `search_sentinel2(...)`, `prepare_monthly_composites(...)`. Everything else is internal.
- `config/` — typed TOML schema (pydantic), env-secret overlay, and resolution: `SystemConfig` + `RunRequest` → `ResolvedConfig`. Stable machine settings live in TOML (see `docs/siac-config.example.toml`); per-run inputs (AOI, scene, output path) do not.
- `app/` — turns resolved config into runtime callables: `NamedRegistry` registries in `app/registry.py` (atmo provider, BRDF provider, monthly composites, S2 backend, surface-prior method), `_assembly_*.py` modules per component, and `planning.build_execution_plan(...)` producing the `ExecutionPlan` — the boundary between setup and execution.
- `workflows/` — orchestration: `pipeline.run_pipeline(...)` executes the staged pipeline; `scene.process_scene(...)` adds output dispatch (writes only when an output path AND writer are configured); `sentinel2.py` resolves S2 queries (SAFE path, product ID, or `T31UDQ_20240101` shorthand) before delegating to the generic scene flow. Execution backends: `thread` or `dask`, with retries/timeouts resolved in `workflows/_pipeline_config.py`.
- `algorithms/` — the science: `surface/` (BRDF-derived priors, spectral mapping, SWIR refinement), `solver/` (multigrid aerosol inversion), `correction/` (TOA→BOA), `cloud/`, `grid/` (solver-input assembly), `brdf/` (kernels), `rt/` (radiative-transfer backends).
- `adapters/` — external systems: `data/` (CDSE, GCS Sentinel-2, earthaccess catalog), `atmo/` (CAMS, MERRA-2, MCD19), `brdf/` (MCD43 via earthaccess), `satellite/` (sensor preprocessors), auth.
- `domain/` (pure types + protocols), `runtime/` (xarray execution payloads + validation), `catalog/` (built-in sensor definitions), `geo/` (canonical resampling/reprojection), `storage/` (readers/writers: GeoTIFF, COG, NetCDF, Zarr).

### Pipeline stages M1–M6 (data contracts between stages)

M1 preprocessing → `ObservationBundle`; M2 atmospheric prior → `AtmosphericState`; M3 surface prior → `SurfacePrior`; M4 grid assembly → `SolverInputBundle`; M5 solver → `SolvedAtmosphere`; M6 correction → `CorrectionResult`. These payload types live in `runtime/models.py` with validators in `runtime/validation.py`.

### RT backends

`algorithms.rt` supports four backends selected by `algorithms.rt.backend` in config: `emulator` (default; NN forward pass in Rust), `lut` (Zarr LUT, remote-capable via HTTP zip store), `sixs` (native compiled 6SV2.1), `libradtran` (compiled uvspec). Backend dispatch is registry-driven in `adapters/rt.py` / `app/_assembly_rt.py`.

### Rust extension

`src/siac_rs/` compiles to `siac._rust` via maturin/PyO3: BRDF kernels (`kernels.rs`), PSF convolution (`psf.rs`), emulator forward pass (`emulator.rs`), solver grid ops (`optimization.rs`, `optimization_grid.rs`), Whittaker smoothing (`whittaker.rs`).

### CLI

`siac` entry point (`python/siac/cli.py`): `siac process-s2 <query> [--config ... --output-path ... --aoi-file ...]` and `siac prepare-monthly-composites`.

## Conventions

`docs/naming-conventions.md` is enforced in review: one canonical term per concept (`RSRF`, never `srf`/`rsr`), American English, explicit names over `cfg`-style abbreviations in shared code, and standard suffixes — `*Config`, `*Request`, `*Bundle`, `*Result`, `*Diagnostics`, `*Provider`/`*Backend`/`*Preprocessor`/`*Writer`. Only domain-standard acronyms (AOT, BOA, BRDF, LUT, RSRF, TCWV, TOA, SZA/VZA/RAA, ...) may appear in public names.

Ruff (line length 100, isort, pyupgrade, pathlib enforcement) and strict mypy settings (`disallow_untyped_defs`) are configured in `pyproject.toml` — new code should be fully typed. Tests live in `tests/unit/`, `tests/integration/`, `tests/regression/`, `tests/benchmarks/` with markers `slow`, `integration`, `regression`; match the layer of what you change, and keep the 95%/90% coverage gates in mind.
