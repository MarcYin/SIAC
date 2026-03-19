# SIAC Code Structure

Status note (2026-03-19): the package has been cut over to the new layered
layout. The old `siac.core.*` modules and the old top-level `siac.pipeline`
module have been removed.

## Package Layout

```text
siac/
├── __init__.py
├── siac.py                      # Compatibility shim to the public API
├── errors.py                    # Shared exception types
│
├── api/                         # Public API layer and request models
│   ├── __init__.py
│   ├── public.py                # Canonical SIAC facade and entrypoints
│   └── requests.py              # Typed workflow/public request objects
│
├── config/                      # System config, loading, resolution, snapshots
│   ├── __init__.py
│   ├── public.py                # SIACConfig public wrapper
│   ├── schema.py                # Typed TOML schema
│   ├── load.py                  # TOML I/O and env overlays
│   ├── resolve.py               # SystemConfig + RunRequest -> ResolvedConfig
│   └── snapshot.py              # Resolved/system snapshots
│
├── catalog/                     # Static built-in catalog data
│   ├── __init__.py
│   └── sensors.py               # Built-in sensor definitions and registry
│
├── domain/                      # Pure domain models and contracts
│   ├── __init__.py
│   ├── aoi.py                   # AOI container
│   ├── contracts.py             # Geometry, atmosphere, solver/correction contracts
│   ├── sensors.py               # SensorBand and SensorConfig types only
│   ├── protocols.py             # Structural interfaces for pluggable modules
│   ├── validation.py            # Contract validators
│   └── spectral.py              # Spectral helper utilities
│
├── adapters/                    # External systems and backend adapters
│   ├── __init__.py
│   ├── atmo/                    # Atmospheric-prior data-source adapters
│   ├── brdf/                    # BRDF product data-source adapters
│   ├── auth.py                  # CredentialManager and provider auth helpers
│   ├── earthdata.py             # Earthaccess source wiring
│   ├── earthdata_common.py      # Shared MODLAND/Earthdata grid helpers
│   ├── rt.py                    # RT backend assembly
│   ├── satellite/               # Sensor preprocessors
│   └── sentinel2.py             # Sentinel-2 backend assembly
│
├── algorithms/                  # Numerical and retrieval algorithms
│   ├── __init__.py
│   ├── brdf/                    # BRDF kernel computations
│   ├── cloud/                   # Cloud/cloud-shadow detection
│   ├── correction/              # TOA -> BOA atmospheric correction
│   ├── grid/                    # Grid assembly for solver inputs
│   ├── rt/                      # Radiative-transfer backends
│   ├── solver/                  # Aerosol inversion
│   └── surface/                 # Surface-prior derivation and spectral mapping
│
├── app/                         # Runtime planning and component assembly
│   ├── __init__.py
│   ├── registry.py              # Kind/backend registries
│   ├── assembly.py              # Config -> runtime callable assembly
│   └── planning.py              # ExecutionPlan construction
│
├── workflows/                   # Orchestration workflows
│   ├── __init__.py
│   ├── pipeline.py              # Core pipeline orchestrator
│   ├── scene.py                 # Generic scene execution flow
│   └── sentinel2.py             # Sentinel-2 search/download/process flow
│
├── geo/                         # Geometry and reprojection utilities
├── storage/                     # Raster/product read-write helpers
├── io/                          # Thin convenience facade over geo/storage/data
└── srf/                         # Spectral response function domain
```

## Layer Boundaries

- `catalog/` owns built-in static lookup data such as bundled sensor
  definitions.
- `domain/` contains contracts and types only. It does not own auth, config
  loading, filesystem discovery, or remote service access.
- `config/` owns the canonical TOML schema and config resolution.
- `adapters/` owns external-system integration concerns.
- `algorithms/` owns retrieval math, surface priors, correction, cloud masking,
  RT backends, and grid prep.
- `app/` resolves configuration into concrete runtime components.
- `workflows/` executes end-to-end processing plans.
- `api/` is the canonical public surface.
- `siac.py` is intentionally thin and only forwards to `api/`.
- `io/` is intentionally narrow and exists only as a convenience facade over
  `geo/`, `storage/`, and `adapters/data/`.

## Canonical Runtime Flow

1. Load `SIACConfig` from `siac.config`.
2. Build a `RunRequest` / resolve a `ResolvedConfig`.
3. Assemble runtime dependencies in `siac.app`.
4. Build an `ExecutionPlan`.
5. Execute via `siac.workflows.scene` or `siac.workflows.pipeline`.

## Public API

- `siac.SIAC`
- `siac.SIACConfig`
- `siac.siac_process_s2(...)`
- `siac.search_sentinel2(...)`

Everything else should generally be treated as internal package structure.
