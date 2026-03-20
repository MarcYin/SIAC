# SIAC Naming Conventions

This document defines the default naming rules for SIAC code, tests, and docs.
The goal is simple: one concept should have one name across the package.

## Core Rules

1. Use one canonical term per concept.
   Relative spectral response is `RSRF` everywhere in SIAC-owned identifiers.
   Do not mix `srf`, `rsr`, and `rsrf` for the same concept.

2. Use explicit names in production code.
   Prefer `config`, `auth_config`, `cloud_mask_config`, and `solver_config`
   over ambiguous names like `cfg` in shared implementation code.

3. Use American English in code and docs.
   Write `center`, `catalog`, `normalize`, and `organization`.

4. Use full type names for domain and runtime objects.
   Examples: `RelativeSpectralResponse`, `SensorConfig`,
   `SceneProcessRequest`, `CorrectionDiagnostics`.

5. Use standard suffixes consistently.
   - `*Config` for persistent configuration models
   - `*Request` for workflow or API requests
   - `*Bundle` for grouped runtime payloads
   - `*Result` for outputs
   - `*Diagnostics` for timing and metrics
   - `*Provider`, `*Backend`, `*Preprocessor`, `*Writer` for implementations

## Acronyms

Use only domain-standard acronyms in public names:

- `AOI`
- `AOT`
- `BOA`
- `BRDF`
- `CDS`
- `CDSE`
- `GCS`
- `LUT`
- `RAA`
- `RSRF`
- `RT`
- `S2`
- `SZA`
- `TCWV`
- `TOA`
- `VZA`

If an abbreviation is not standard for this domain, spell it out.

## Spectral Naming

- Canonical domain type: `RelativeSpectralResponse`
- Sensor-band fields: `rsrf_wavelengths_nm`, `rsrf_response`
- Capability check: `has_rsrf`
- Reference-loader functions: `load_reference_rsrf`
- LUT helpers: `rsrf_kernel`

Use `SRF` only when quoting an external source title or file name that SIAC does
not control, such as `Sentinel-2 Spectral Response Functions (S2-SRF)`.

## Module and File Names

- Use lowercase `snake_case` module names.
- Name each module for one concept when practical.
- Keep package-local adapter names aligned with the external system they wrap.
  Examples: `rsrf.py`, `s2_backend.py`, `copernicus_dataspace.py`.

## Variable Names

- Shared implementation code should prefer descriptive locals.
- Short names like `i`, `x`, and `ax` are acceptable only for tiny scopes where
  the meaning is obvious.
- Test-local helper names may be shorter, but public fixtures, helper
  functions, and reusable utilities should still be descriptive.

## Exceptions

- Preserve external names when SIAC is mirroring a third-party schema or API.
  Example: a library payload field named `rsr`.
- Preserve quoted document titles and URLs as published by the source.

## Review Checklist

Before adding a new name, check:

- Does this concept already have a canonical SIAC name?
- Is the name explicit enough outside a tiny local scope?
- Does it match American spelling?
- Does it use an approved suffix when applicable?
- Is an acronym truly standard for this domain?
