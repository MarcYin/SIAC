# SIAC Codebase Review

**Branch reviewed:** `feat/refactor` @ `d0baf45` (worktree HEAD is identical — 0 commits ahead, 0 behind)
**Compared against:** `master` @ `fc6e592` — 569 files changed, +124,568 / −31,491
**Scope:** every Python file under `python/siac/`, `tools/`, `tests/`, plus the Rust crate `src/siac_rs/`
**Methodology:** eight parallel reviewers, each reading its assigned files line-by-line; findings consolidated and de-duplicated below
**Out of scope:** `docs/**`, `pixi.lock`, `.github/**`, `.gitignore`, the lockfile-style TOMLs

Each finding is `path:line — short description`. They are organised by package layer; a "Highest-Impact Issues" curation comes first.

---

## 1. Highest-Impact Issues

These are the findings most likely to produce silent wrong outputs or broken developer flows. Address these first.

### 1.1 Broken / outright incorrect

1. **Several tools import modules that no longer exist** — they will `ModuleNotFoundError` on first run:
   - `tools/build_6s_native.py:13` — `from siac.config.schema import SixSAlgorithmConfig` (module gone, lives in `siac.config.algorithms`)
   - `tools/compare_6s_route_coefficients.py:20`, `tools/compare_native_6s_to_remote_lut.py:20` — same `siac.config.schema` import
   - `tools/compare_surface_prior_approaches.py:24`, `tools/compare_surface_prior_experiment.py:18` — `from siac.app.assembly import …` (module renamed to `_assembly_*`)
   - `tools/profile_m3.py:85` — `from siac.app.assembly import build_pipeline_runtime` and `from siac.config import load_config` (neither exists)
2. **CLI mis-routes WKT input** — `python/siac/cli.py:183-186` passes `--aoi-wkt` to `AOI.from_geojson(...)`. WKT is not GeoJSON; this silently corrupts the AOI.
3. **Public re-export drift** — `siac.api.__init__` re-exports `siac_process`, but `siac/__init__.py:_LAZY_IMPORTS` does not list it; `from siac import siac_process` fails with AttributeError.
4. **L-BFGS-B QA fields are silently zero/None** — `python/siac/algorithms/solver/multigrid.py:526, 540-547`: `final_zero_obs_mask`, `final_invalid_mask`, `final_insufficient_support_mask`, `final_fitting_cost` are only populated by the grid-search branch. The L-BFGS-B path always reports clean QA even when retrievals failed.
5. **Sentinel-2 view-angle averaging is a running-average bug** — `python/siac/adapters/satellite/sentinel2.py:705-731`: `mean_vza = (mean_vza + …) / 2` per band biases the result toward later entries; not the arithmetic mean.
6. **Solver `ftol` is effectively zero** — `python/siac/algorithms/solver/multigrid.py:113`: `1e-7 * np.finfo(float).eps ≈ 2.22e-23`. L-BFGS-B never stops on function tolerance; convergence depends on `gtol` only.
7. **Fixed-parameter objective: gradient/cost mismatch** — `python/siac/algorithms/solver/multigrid.py:1648-1666`: when `fixed_parameter == "aot"`, only `grad[n:]` is returned, but `cost` still includes prior/smoothness contributions from the fixed half — breaks Wolfe line-search heuristics.
8. **Sparse Laplacian boundary indexing likely wrong on non-square grids** — `python/siac/algorithms/solver/cost.py:601-620`: rows/cols derived as `idx // ny`, `idx % ny`, but boundary masks reference `nx-1`/`ny-1`. Off-by-one on non-square grids (and the function may be dead code — verify before fixing).
9. **`prior_store._crop_to_bounds` assumes y-up origin** — `python/siac/algorithms/surface/prior_store.py:147-158`: standard rasters with y-down are mis-cropped.
10. **Leap-year DOY math wrong** — `python/siac/algorithms/surface/prior_store.py:75-84`: uses `365` not `366`; off-by-one DOY interpolation in leap years.

### 1.2 Silent numerical / scientific risk

1. **Spectral-math sign flip near zero** — `python/siac/algorithms/rt/lut/_spectral_math.py:121`: `np.where(|d|<1e-10, 1e-10, d)` replaces small negatives with `+1e-10`, flipping sign.
2. **Single-wavelength sampling for compact-coefficient LUT** — `python/siac/algorithms/rt/lut/backend.py:357-359`: argmin-nearest instead of RSRF integration. Bands spanning 30+ nm sample only the centre.
3. **Cost weights divide by zero** — `python/siac/algorithms/solver/cost.py:174`: `weights / weights.sum()` produces NaN when all wavelengths sum to zero. The same hazard at `multigrid.py:854` IS guarded; cost.py is not.
4. **Index-clamp defeats out-of-range NaN fill** — `python/siac/algorithms/rt/lut/backend.py:501-511`: `_sanitize_point_indexers` clips to LUT bounds; combined with `bounds_error=False, fill_value=NaN` in the interpolator, out-of-range pixels silently get the boundary value rather than NaN.
5. **Heuristic ozone-units conversion** — `python/siac/algorithms/rt/lut/backend.py:482-499`: detects atm-cm vs DU by magnitude (`> 20` and `< 10`). False positives possible.
6. **Magic atmosphere-profile selection by month bucket** — `python/siac/algorithms/rt/direct/sixs_native.py:325-435`: April mapped to "January-style" winter profile.
7. **Heuristic posterior uncertainty** — `python/siac/algorithms/solver/multigrid.py:1769-1774`: scaled by spatial-mean of retrieved AOT — same scene at different crop sizes yields different uncertainty. Comment claims "posterior uncertainty" but this isn't Fisher info.
8. **Historical monthly composites are computed at *current* observation geometry** — `python/siac/algorithms/surface/swir_refine.py:842`: passes the obs geometry into a 5-year history. Confirm with science whether this is by design.
9. **Float32→Float64→Float32 round-trip per multigrid level** — `python/siac/algorithms/solver/multigrid.py:1538`. Memory and conversion cost on full scenes.
10. **Mean-over-altitude collapses RT non-linearity** — `python/siac/algorithms/rt/lut/backend.py:1163`: `subset.mean(dim="altitude")` after a multi-node slice — RT's altitude dependence is typically non-linear.

### 1.3 Resource / write-safety

1. **HDF4 dataset handles are leaked on every read** — `python/siac/adapters/earthdata_common.py:483, 507, 528`: `read_hdf4_dataset` opens an `SD` and never calls `sd.end()`. On long batches this exhausts file descriptors.
2. **No atomic-write semantics anywhere in `storage/writers.py`** — GeoTIFF / COG / NetCDF / Zarr all write directly to the final path. A crash mid-write leaves a partial file that downstream code consumes. Same issue for STAC JSON in `adapters/output.py:209`.
3. **CAMS cache writes are non-atomic** — `python/siac/adapters/atmo/cams.py:763-781`: partial cache files survive crashes and pollute later runs; no atomic rename, no `.tmp`, no lockfile.
4. **Monthly composite store overwrites before atomic swap** — `python/siac/algorithms/surface/monthly_composite_store.py:265-270`: deletes existing data before writing new; crash mid-write leaves an inconsistent period.
5. **HTTP zip store caches arbitrarily large bodies** — `python/siac/algorithms/rt/lut/http_zip_store.py:111-122, 151-157`: `_full_body_cache` is unbounded; a 4 GB LUT zip can balloon RAM.
6. **PIL-based resampling does not handle NaN** — `python/siac/algorithms/grid/assembler.py:106, 163-185`: NaN is silently propagated through `Image.fromarray(mode="F")`; `uniform_filter` window uses integer division → aliasing for non-integer ratios; failure-padding uses `0`, corrupting nodata semantics.

### 1.4 Layering violations

1. **`siac.geo.resample` imports `siac.runtime.models`** — `python/siac/geo/resample.py:18,20-21`. Bidirectional dependency between low-level geo and runtime; circular-risk.
2. **`siac.storage.writers` and `siac.storage.stac` import from `runtime` / `algorithms` / `workflows`** — same direction problem.
3. **`siac.sixs_upstream_parity` imports `_*` private helpers from `siac.algorithms.rt.direct.sixs_native`** — `python/siac/sixs_upstream_parity.py:21-29`. Top-level module reaches into a sibling subpackage's private names.
4. **`siac.adapters.brdf.mcd43_earthaccess` imports `_build_virtual_stack_vrt`, `_target_bounds_from_template`** from `adapters.earthdata` — `python/siac/adapters/brdf/mcd43_earthaccess.py:21-31`. Adapters reaching into siblings' private surfaces.
5. **`siac.config.public.SIACConfig` inherits from `SystemConfig`** — `python/siac/config/public.py:42`. Public model depends on a "system" (private) one; the wrapper-only `sensor`/`aoi` extensions are filtered out by `_system_only()` for serialization, so round-trip silently loses run-input fields.
6. **`siac.domain.aoi` imports from `siac.geo.geometry` / `siac.geo.reprojection`** — `python/siac/domain/aoi.py:14`. The `domain` package's docstring claims "pure domain-layer types and protocols."

---

## 2. Cross-Cutting Themes

These appear in many files — fix the pattern, not just the instance.

### 2.1 Bare/wide `except` swallowing meaningful errors

The codebase is dotted with `except Exception:` returning a default or fallback silently. Sample sites:
- `python/siac/_rust_compat.py:14-21` swallows everything during native load
- `python/siac/adapters/earthdata.py:113, 494-509, 587-592` (multiple)
- `python/siac/adapters/atmo/cams.py:285-290, 1058-1063` swallow auth/quota errors as "data not available"
- `python/siac/adapters/atmo/mcd19_earthaccess.py:108, 248`
- `python/siac/adapters/atmo/merra2.py:101`
- `python/siac/adapters/output.py:494-549` demote preview failures to debug
- `python/siac/adapters/rt.py:60-67` silently degrades emulator → LUT
- `python/siac/adapters/satellite/sentinel2.py:122, 325-337` silently switch reproject paths on any error
- `python/siac/algorithms/cloud/providers/omnicloudmask.py:27` masks partial-install failures
- `python/siac/algorithms/rt/lut/backend.py:165, 200-208, 251-256` etc.
- `python/siac/algorithms/surface/kernel_model.py:189`
- `python/siac/algorithms/surface/spectral_mapping.py:402`
- `python/siac/algorithms/surface/prior_store.py:127`
- `python/siac/runtime/models.py:62-109` (CRS write silenced)
- `python/siac/runtime/models.py` mixes silencing on three rioxarray operations

Recommendation: replace with narrow exception types, and when a fallback is the right semantic, *log at WARNING with original `exc_info=True`* and tag the event in the observability stream.

### 2.2 Inconsistent error taxonomy

`python/siac/errors.py` defines `SIACError`, `ConfigurationError`, `RTModelError`, etc., but most callers raise plain `ValueError` / `RuntimeError` / `TypeError`:
- `python/siac/cli.py:200-215`
- `python/siac/rt_setup.py:93-99, 184-189`
- `python/siac/sixs_upstream_parity.py:768-810`
- `python/siac/algorithms/correction/atmospheric.py:46`
- `python/siac/algorithms/solver/multigrid.py:264`
- many config-level validators

Also `ValidationError(SIACError)` does NOT subclass `ValueError`, so `except ValueError` users miss it.

### 2.3 `runtime_checkable` Protocols that don't actually check

`python/siac/domain/protocols.py:31-153` — every protocol is `@runtime_checkable`, but several have `@property` members (`source_name`, `sensor_config`, `backend_name`, `source_bands`). `runtime_checkable` only checks attribute presence, not property semantics; on Python ≥3.12 this also warns. `isinstance(...)` against these protocols will lie.

Also: `protocols.py:38-40` declares `load_toa(input_path: str)` while the rest of the codebase uses `Path | str`. The catalog/registry split (`is_available_for_sensor(sensor_id, satellite_id)`) is unclear because `SensorName` enum only has `S2/AUTO`.

### 2.4 Magic numbers without citation

Pervasive — a non-exhaustive list:
- `python/siac/algorithms/correction/atmospheric.py:98` — BOA acceptance window `(-0.05, 1.5)`
- `python/siac/algorithms/cloud/mask.py:163` — color windows hard-coded
- `python/siac/algorithms/surface/brdf_whittaker.py:120-121` — fallback `prior=0.20`, `prior_unc=0.08`
- `python/siac/algorithms/surface/kernel_model.py:146, 296, 411` — `1.5` BOA ceiling, `1e-10 + 0.01` threshold (likely a copy-paste artefact)
- `python/siac/algorithms/rt/lut/backend.py:1319-1357` — scale-height = 8.5 km, applied as Beer-Lambert correction to spherical albedo and path reflectance with no physical justification documented
- `python/siac/algorithms/rt/lut/backend.py:1380` — Jacobian step sizes `delta_aot=0.01`, `delta_tcwv=0.1`
- `python/siac/algorithms/solver/cost.py:387` — Pseudo-Huber `delta`-handling diverges at `delta == 0`
- `python/siac/adapters/satellite/sentinel2.py:734-755` — mean angles default to magic `5.0 / 100.0 / 30.0` if XML is missing
- `python/siac/algorithms/grid/assembler.py:454-464` — NaN/Inf coerced to `0/1`

These should at minimum be named constants with comments and ideally promoted to config.

### 2.5 Pydantic v1/v2 mixing & enum-vs-str comparisons

`python/siac/config/algorithms.py` is the worst offender, with `self.kind == "constant"`, `self.mode == "user_defined"`, `self.atmospheric_correction.mode.startswith("brdf_")` etc. Works only because `SIACEnum` is `str, Enum`; a small enum-value rename is silently breaking. Hits at lines 217, 249, 279, 283, 294, 347, 457, 461, 472, 489-491. Same pattern at `_assembly_correction.py` and many user-facing branches.

Also: `model_config` on `_base.py` is a plain dict, not `ConfigDict(...)`. No `frozen=True`, no `validate_assignment=True`, so direct attribute writes (used in `config/load.py:56-86 overlay_env_secrets`) skip validators silently.

### 2.6 `requests` / network usage is uneven

Three transport stacks coexist — `requests` (auth, copernicus_dataspace, water_mask), `urllib.request` (`gcs_sentinel2.py:71-76`), and `fsspec` (cams). All have:
- timeouts but no retry/backoff (`copernicus_dataspace.py:230-238, 374-389`, `water_mask.py:182-205`)
- `_retry.py:13-17`'s `_TRANSIENT_EXCEPTIONS` includes plain `OSError`, which on Linux covers many non-retryable failures (ENOSPC, EROFS, ENOENT)
- backoff at `_retry.py:43-65` has no jitter; concurrent retriers collide deterministically
- token caching missing in `copernicus_dataspace.py:262-274`: every download re-exchanges the OAuth token

### 2.7 Secret redaction is incomplete

`python/siac/adapters/_log_filter.py` is only attached to the `auth` logger (line 29). Other modules (`cams.py:1062`, `copernicus_dataspace.py`, `mcd43_earthaccess.py:586`) log raw URLs that may carry tokens. Bearer regex matches only contiguous `[A-Za-z0-9._\-]` so JWT `+` / `=` / `/` characters leak through.

### 2.8 Resource lifecycle on storage

- `siac.storage.readers.read_raster, read_netcdf_variable, read_multiband` — all open datasets and never close them (`readers.py:95, 282-303, 412-421`). With S2's 13 bands, file handles leak per scene.
- Casts to float32 destroy integer raster precision (`readers.py:175-178`).
- CRS comparison is string-equality (`readers.py:495-527`): `"EPSG:4326"` ≠ verbose WKT for the same authority.

### 2.9 Tools/Tests

- 4 of 17 tools have outright broken imports (see §1.1).
- `tools/profile_m3.py:33-80` monkey-patches dozens of internal functions; if any is renamed the patch silently no-ops.
- `tools/compare_6s_route_coefficients.py:73-78` and `tools/compare_native_6s_to_remote_lut.py:73-78` duplicate `_band_catalog()` / `_select_bands()` verbatim. Several of the `compare_*` tools fork ~250 LOC of plotting from each other.
- Test fixtures in `tests/conftest.py:18-26` reference `tests/fixtures/` which does not exist.
- `tests/regression/__init__.py` is aspirational — zero regression tests live under it.
- `tests/conftest.py:107, 130` use unseeded `default_rng()` for "sample" fixtures, making them non-deterministic.
- `tests/integration/test_e2e_synthetic.py` is not marked `@pytest.mark.integration`, so it runs in the fast suite.
- `tests/integration/test_gcs_sentinel2_real.py:22` hard-codes a public GCS product ID; if GCS ages it out, the network catch-all silently turns the test into a skip.
- 55 of ~85 unit-test files import private symbols; refactor cost on private internals is amplified by tests.

### 2.10 Rust crate

- `src/siac_rs/Cargo.toml:14-22` declares `ndarray-stats` and `num-traits` deps that are unused (no source matches).
- `src/siac_rs/benches/kernels.rs` is a no-op stub: `cargo bench` succeeds without measuring anything.
- The parallel hot loops in `optimization.rs:128-189, 196-226, 608-656, 767-799`, `whittaker.rs:96-137`, `kernels.rs:64-73` do NOT release the GIL via `py.allow_threads`. Rayon parallelism still works inside Rust, but Python threads sit blocked.
- `optimization.rs:914-1163` `refine_grid_search_with_qa` runs serially on `iy in 0..ny` despite no shared mutable state — obvious rayon target.
- `psf.rs:56-62` docstring says "DCT-based Gaussian" but the implementation is direct spatial-domain weighted averaging.
- `optimization.rs:14-19` `enum FixedParameter` has variant `None` shadowing `Option::None`.

---

## 3. Per-Module Findings

### 3.1 Entry points

#### `python/siac/__init__.py`
- `__init__.py:8-15` — Docstring example imports `from siac.config import SIACConfig`, but the lazy `__getattr__` already exposes `SIACConfig` from `siac` directly. `_LAZY_IMPORTS` does not list `siac_process` even though `siac/api/__init__.py` re-exports it.
- `__init__.py:28-37` — No `__all__`; `dir(siac)` won't list lazy names; `from siac import *` silently fails to load `SIAC` / `SIACConfig`.

#### `python/siac/_rust_compat.py`
- `_rust_compat.py:14-21` — `_load_native_rust` swallows everything via bare `Exception`, with import side effects.
- `_rust_compat.py:39-45` — `_require_native` raises `ImportError`, not the `errors.py` taxonomy; fallback path returns `Any`.
- `_rust_compat.py:48-87` — `__getattr__` recursion can opaquely raise the same error.

#### `python/siac/cli.py`
- `cli.py:183-186` — **Bug.** `--aoi-wkt` is passed to `AOI.from_geojson(...)`. WKT ≠ GeoJSON.
- `cli.py:218-225` — `_run_process_s2` hands `args.query` raw to `siac_process_s2` with no shorthand parser despite help text claiming SAFE path / product ID / tile/date are accepted.
- `cli.py:228-244` — Dead branch: `_resolve_cli_aoi` already returns `AOI` so the `if not isinstance(aoi, AOI):` guard is unreachable.
- `cli.py:262-264` — `parser.error` later raises `SystemExit` not caught here; `main()` exits process directly.
- `cli.py:200-215` — Bare `ValueError` for "either/or" + duplicate detection; not aligned with `errors.py`.
- `cli.py:224, 247-254` — `print(...)` mixed with logger output for CLI feedback.

#### `python/siac/api/__init__.py`, `python/siac/api/public.py`, `python/siac/api/requests.py`
- `api/__init__.py` — re-exports `siac_process` not registered in top-level lazy-imports (see §1.1).
- `api/public.py:97-122` — `search_sentinel2` accepts both `date` and `date_value` with no precedence; both pass through to `Sentinel2SearchRequest`.
- `api/public.py:125-130` — `process_sentinel2(**kwargs: Any)` silently breaks with `TypeError` for any kwarg outside `output_path`/`aoi`.
- `api/public.py:151-164` — `siac_process` cannot set `output_path`; `siac_process_s2` can. Latent missing parameter.
- `api/public.py:18-31` — Unused re-exports `apply_s2_query_defaults`, `coerce_date`, `coerce_s2_query`; leak of internal helpers.
- `api/requests.py:1-26` — Pure shim with no `__all__`; partial-migration pattern.

#### `python/siac/public_models.py`
- `public_models.py:1-19` — frozen dataclass with `tuple[str, ...]` fields; missing `__all__`.

#### `python/siac/runtime/models.py`
- `runtime/models.py:62-109` — three bare `except Exception` silencing rioxarray ops.
- `runtime/models.py:131-137` — `__post_init__` validates types; `RTCoefficients.__post_init__` (line 324) has to use `object.__setattr__` because of `frozen=True` — inconsistency with `GeometryAngles` which doesn't normalize.
- `runtime/models.py:142` — `raa = abs(vaa - saa) % (2π)` discards sign of relative azimuth; downstream consumers expecting signed RAA see folded values.
- `runtime/models.py:310-324` — `extras` annotated `Mapping` but replaced with `MappingProxyType` in `__post_init__`; type lies about runtime mutability.
- `runtime/models.py:355-359` vs `:362` — different stability gates between forward and inverse paths (`apply_correction` thresholds `denom`; `simulate_toa` checks `xap`).
- `runtime/models.py:383-392` — `compute_boa_jacobian` indexes via `.sel(param="aot")` without validating dims.
- `runtime/models.py:419-422` — `compute_reflectance_uncertainty` returns cached value regardless of `k_vol`/`k_geo` arguments — silent wrongness.
- `runtime/models.py:432, 466` — `monthly_composites: tuple[Any, ...]` is `Any`-typed; `bands: list[SensorBand]` inside a `frozen=True` dataclass — false sense of immutability.

#### `python/siac/runtime/validation.py`
- `validation.py:25-32` — `spatial_shape` falls back to `list(ds.dims)` whose ordering is implementation-defined for older xarray.
- `validation.py:73-78, 153-159` — `.values` materialises lazy dask arrays inside a validator. Unexpected side effect.
- `validation.py:135-140` — Shadows builtin `field`.

#### `python/siac/errors.py`
- `errors.py:17-30` — `stage` and `is_recoverable` are class attributes; `stage` is never set per-class so it's always `None`.
- `errors.py:57-58` — `ValidationError(SIACError)` does not inherit from `ValueError`; `except ValueError` users miss it.

#### `python/siac/observability.py`
- `observability.py:36, 48-53` — mixes `isinstance(value, str | int | float | bool)` PEP 604 union with tuple form elsewhere; `getattr(value, "item", None)` then call-and-swallow.
- `observability.py:217-246` — `finish()` joins sampler thread with timeout; on timeout the daemon thread keeps writing to `summary_path`, racing later atomic writes.
- `observability.py:436-445` — `_append_event_locked` opens a new file handle per event (no buffering) inside the lock — disk latency now blocks all observer methods.
- `observability.py:443-445` vs `record_error` — `record_error` directly appends to `_errors` without going through `emit`, so it skips events file persistence.
- `observability.py:519-524` — `read_bytes` / `write_bytes` summary fields are *overwritten* with current sample rather than tracked as deltas.

#### `python/siac/rt_setup.py`
- `rt_setup.py:43-52` — `_coerce_rt_setup` accepts non-pydantic via `vars()`, leaking private attributes.
- `rt_setup.py:55-60` — `_merge_model_payload` non-pydantic branch is unreachable. Dead.
- `rt_setup.py:93-99, 184-189` — `raise ValueError` instead of `ConfigurationError`.

#### `python/siac/sixs_outputs.py`
- `sixs_outputs.py:69, 72, 75` — three entries declare arithmetic-expression strings (`"sdtotr*sutotr"`) where every other tuple element is a Fortran name. Downstream parsers expecting tuple[1] as a name will break.
- `sixs_outputs.py:133-138` — `_all_names` mutated at module-level import; `_name` leaks as a module attribute.

#### `python/siac/sixs_upstream_parity.py`
- `:11, :437-470` — `pickle` over subprocess stdin; child errors re-raised as `RuntimeError`, swallowing the original traceback.
- `:21-29` — imports `_*` private helpers from a sibling subpackage.
- `:393-431` — mutates `case.sixs_config` via `model_copy(update={...})`; raises if config drops the field.
- `:446-462, 467-487` — Inline Python script as a string literal passed to subprocess; assumes a fixed package layout (`__file__.resolve().parents[2]`); wrong if installed.
- `:633-637` — `Path("")` truthy check is fragile.
- `:768-810, 1000-1014` — subprocess errors re-raised as `RuntimeError`, taxonomy lost.

### 3.2 Config + Domain

#### `python/siac/config/_base.py`
- `_base.py:9` — `model_config` is plain dict; should be `ConfigDict(...)`. No `frozen` / `validate_assignment` policy.

#### `python/siac/config/algorithms.py`
- `algorithms.py:217-491` (multiple) — enum-vs-string equality across `kind`, `mode`, `parallel_backend`, `atmospheric_correction.mode`. Works only because `SIACEnum` inherits `str`.
- `algorithms.py:255-256, 312, 316-317, 521-522` — magic numbers (`13`, `10`, `512`) without docstring/citation; cross-field validation missing for `(month, day)`; AOT defaults duplicate constants likely defined elsewhere.
- `algorithms.py:325-343` — `normalize_output_variables` silently drops empty/whitespace items, dedupes, falls back to defaults when input is `[""]`.
- `algorithms.py:347-351` — Validator `validate_custom_aerosol` actually validates worker/chunk/scene_lut bounds (misnamed); the bounds it checks are already enforced by `Field(ge=1)` — dead code.
- `algorithms.py:498-508` — `has_overrides()` `is not None` semantics rely on RT setup tree always being constructed; comment-only.
- `algorithms.py:528, 543-549, 589-590, 605, 634, 206-207` — `aot[0] >= 0` not enforced; `normalize_parameters` allows empty `solve=()` no-op stage; `normalize_blur_kernel` silently mutates user even values; `stages: tuple[…] = Field(default_factory=tuple)` allows zero-stage solver; `user_callable: Any | None` defeats schema and can't round-trip TOML; `kind="spectrum"` doesn't enforce `wavelengths_um`.

#### `python/siac/config/load.py`
- `load.py:25` — `tomli_w = _TomliWFallback()` instance shadows module name; fallback raises unconditionally.
- `load.py:42` — `isinstance(loaded, dict)` is dead defensive (tomllib always returns dict).
- `load.py:48-50` — Empty `CONFIG_PATH_ENV` becomes `Path("")` which is `Path('.')` — silent cwd redirect.
- `load.py:56-86` — `overlay_env_secrets` mutates the deep-copied model via direct attribute assignment, bypassing validators.
- `load.py:93` — `model_dump(mode="json", exclude_none=True)` makes user-cleared fields indistinguishable from un-set.

#### `python/siac/config/providers.py`
- `providers.py:39` — `lut_path` defaults to a remote https URL silently.
- `providers.py:81-83` — `GCSAuthConfig.credentials_file_env` collides with Google's own SDK env var.
- `providers.py:121, 138, 150` — undocumented constraints; `max_cloud_cover=100.0` ambiguous; `user_callable: Any | None` opaque.

#### `python/siac/config/public.py`
- `public.py:32-39` — `_deep_merge` only handles dict-on-dict; pydantic models or list-of-dict overwrite silently.
- `public.py:42` — Public `SIACConfig` inherits `SystemConfig` (private). Round-trip via `_system_only()` loses `sensor`/`aoi`.
- `public.py:52-54` — `aoi: dict | Path | str | tuple[float,…] | list[float] | None` overlapping `tuple` and `list[float]`; first match wins.
- `public.py:67-82, 90-91` — `from_yaml`/`to_yaml` raise `ValueError` rather than `NotImplementedError`; stale stubs.
- `public.py:115` — `write_state_snapshot` ignores `input_path`, `output_path`, `s2_query`.
- `public.py:139-155` — `get_jasmin_config()` hard-codes a single environment; library-vs-fixture issue.

#### `python/siac/config/request.py`, `resolve.py`, `resolved.py`, `snapshot.py`, `system.py`, `types.py`
- `request.py:14-19` — no validator forcing at least one of `input_path`/`s2_query`.
- `request.py:21-29` — `s2_query: str | Path | None` allows typo'd local path to pass through as a "remote URL".
- `resolve.py:92` — `or "auto"` literal fallback is `str`, not `SensorName.AUTO`; works only by `str, Enum`. Also `getattr(system, "sensor", …)` is a runtime-typing escape because `SystemConfig` lacks `sensor`.
- `resolve.py:97-128` — `**model.model_dump(...)` then manually setting `cache_dir` is fragile if base gains new fields. `solver=system.algorithms.solver.model_copy()` returns the base class while everything else uses `Resolved*` — naming inconsistency.
- `resolved.py:28, 84, 94-99, 122-129` — Multiple `Resolved*` classes that are pure aliases or empty subclasses. The `*_env` fields are still present in resolved output, contradicting apparent intent.
- `snapshot.py:13` — Lazy `import_module("yaml")` without `try/except`; if PyYAML is missing, every snapshot import fails.
- `snapshot.py:19, 57-70, 99-113` — Hard-coded URI scheme set, redaction list and 9 path fields; new fields silently omitted.
- `system.py:48-66` — `n_jobs` no constraint on `0`; `chunk_size` overlaps with `SixSAlgorithmConfig.chunk_size`; `stage_timeout_s` vs `stage_timeouts` precedence undocumented.
- `types.py:108-118` — `user_water_ozone` enumerated profile but no validator demanding payload.
- `types.py:121-153, 225-228` — Duplicate parallel enums (`SixSAerosolProfile` vs `RTAerosolProfile`, `FixedAtmosphericParameter` vs `AtmosphericParameterName`).
- `types.py:286-296, 309-316` — `_coerce_path_or_url` inconsistent for `Path` inputs; `_coerce_int_tuple` declared but unused (dead).

#### `python/siac/domain/aoi.py`
- `aoi.py:54-55, 80-85, 90-99, 14, 115-128` — `hasattr(raster, "rio")` masks rioxarray import errors; deprecated GeoJSON `crs` member silently honoured; `Path(str(aoi))` wasteful when already a Path; layering violation (domain → geo); `xmin == xmax` rejected only on one path.

#### `python/siac/domain/protocols.py`
- See §2.3.

#### `python/siac/domain/sensors.py`
- `sensors.py:12-25, 75, 96-131, 133-164` — `frozen=True` with numpy arrays + manual `object.__setattr__` defeats hash/pickle; `np.ndarray` typed without dtype/shape; O(N²) duplicate detection; `get_band_by_wavelength` and `select_nearest_band` near-duplicates with different defaults; `default_aerosol_solver_bands` hard-codes `"MSI"` inside generic dataclass; `swir_bands` excludes water-vapor windows by hard-coded magic numbers.

#### `python/siac/domain/spectral.py`
- `spectral.py:14-104, 130-145` — `_trapezoid` numpy 1.x/2.x compat with hand-rolled formula (undocumented); equal-y interpolation falls through silently; `from_tabulated` aggressively trims with `+2` slice convention; `__post_init__` mutates frozen fields.

### 3.3 Adapters

#### `python/siac/adapters/_log_filter.py`
- `_log_filter.py:9-15, 22-23` — Bearer regex misses JWT chars; mutates `record.msg` in-place (other handlers see redacted).

#### `python/siac/adapters/_retry.py`
- `_retry.py:13-17, 43-65` — Plain `OSError` retried; backoff has no jitter, no cap.

#### `python/siac/adapters/auth.py`
- `auth.py:107-132, 120-124, 138-145` — Token-acquisition double-check pattern leaks two concurrent OAuth round-trips; no skew margin on token expiry; POST has timeout but no retry.
- `auth.py:200, 216-226` — Catches `(OSError, PermissionError, ConnectionError)` during S3 activation polling, burying root cause; `revoke()` in `finally` re-raises and masks original auth failure.
- `auth.py:305-336` — `activate_environment` writes credentials into process env (race with concurrent calls; leaks to subprocesses); `source_kwargs` sets `persist=False` twice.

#### `python/siac/adapters/earthdata.py`
- `earthdata.py:113, 217-242, 386-391, 494-509, 587-592, 991-1013` — Multiple silenced exception paths; identity-comparison against module-level `_DEFAULT_REPROJECT_NATIVE_TO_TARGET` is brittle; `gdal.Unlink` failures leak; `_dataset_bounds`/`_dataset_resolution` swallow CRS misconfigurations and silently warp the whole tile; rasterio MemoryFile writer not in the outer `ExitStack` — leak on exception.

#### `python/siac/adapters/earthdata_common.py`
- `earthdata_common.py:42-51, 149-158, 158-187, 205-221, 616` — `_GRID_METADATA_RE` greedy across grid sections; `lru_cache(maxsize=128)` keyed on path string with no invalidation on file replacement.
- `earthdata_common.py:483, 507, 528` — **HDF4 SD handles never closed** (see §1.3).
- `earthdata_common.py:262-269, 218-220` — `except Exception: return None` swallows h5py errors as "no native grid".
- `earthdata_common.py:660-663` — `attr_scalar(...) or 1.0` corrupts when scale_factor is legitimately `0.0`.

#### `python/siac/adapters/output.py`
- `output.py:494-549` — Preview failures silently demoted to debug.
- `output.py:209-211` — Non-atomic STAC JSON write; no schema validation.

#### `python/siac/adapters/rsrf.py`
- `rsrf.py:23-30` — Reaches into `rsrf.registry._runtime_release_root` (private API of an external package).

#### `python/siac/adapters/rt.py`
- `rt.py:15-67` — Sentinel placeholders that look monkey-patchable; `from … import TwoLayerNNEmulator as emulator_cls` shadows module-globals on every call; `except Exception` silently downgrades emulator → LUT without `exc_info`.

#### `python/siac/adapters/s2_backend.py`
- `s2_backend.py:21` — `GCSSentinel2Backend()` ignores `auth=auth` and the constructor takes none.

#### `python/siac/adapters/atmo/cams.py`
- `cams.py:285-301` — Corrupt local cache silently triggers re-download (or worse: `local_candidates_found` was set during a failed open, blocking re-download forever).
- `cams.py:763-781, 769` — Non-atomic cache writes; no lockfile; concurrent partial-file race.
- `cams.py:719-722` — NetCDF dataset handle not closed.
- `cams.py:850-871, 1058-1063, 1050` — Long-lived S3 creds with one-shot poll; `cdsapi` errors swallowed; throwaway `CredentialManager()` allocation per download.
- `cams.py:909` — Naive UTC datetime parsing.

#### `python/siac/adapters/atmo/mcd19_earthaccess.py`
- `mcd19_earthaccess.py:108-114, 248-249` — Silenced parse failures; HDF4 file-handle leak via `read_hdf4_dataset`.

#### `python/siac/adapters/atmo/merra2.py`
- `merra2.py:78` — Stale "Phase M0 fallback" comment; provider always returns climatological constants regardless of probe results. Effectively dead code path.
- `merra2.py:101-102` — Silenced probe.

#### `python/siac/adapters/brdf/mcd43_earthaccess.py`
- `mcd43_earthaccess.py:21-31` — Imports `_*` private helpers from `adapters.earthdata` (see §1.4).
- `mcd43_earthaccess.py:62-69, 91, 255-267, 330-336, 478-485` — `_DATA_READ_ERRORS` is overly wide (`KeyError`, `RuntimeError`, etc.) → masks programming errors; magic `_MAX_ONE_SHOT_TEMPORAL_VRT_OUTPUT_BYTES = 1 GiB`; HDF4Error fallback to `OSError` inflates swallow surface.
- `mcd43_earthaccess.py:706-724` — `retry_transient` on a non-idempotent download retries from scratch on every failure; `Path(p).resolve().is_relative_to(cache_root)` warning only — should refuse out-of-cache paths.
- `mcd43_earthaccess.py:817-823, 874-881` — Per-band-mismatch failures hidden as "no finite values" warnings; warning text omits `short_name`.
- `mcd43_earthaccess.py:1187-1213, 1349-1358` — `_load_temporal_payload_vrt[_from_daily_payload_vrts]` are no-ops (`del …; return None`); `if … any(): return temporal` then unconditionally returns `temporal` again — both branches return the same object.
- `mcd43_earthaccess.py:1893, 2026` — Reads dataset attrs only from the first granule and reuses for every other day; collection mixing silently mis-scales.
- `mcd43_earthaccess.py:2106-2112, 2601-2607` — `gdal.Unlink` in `finally` may raise and mask the original exception; lazy import inside `finally`.

#### `python/siac/adapters/data/copernicus_dataspace.py`
- `copernicus_dataspace.py:67-80` — `_to_datetime_range` has tangled control flow; brittle `if`/`return` ordering.
- `copernicus_dataspace.py:230-238, 262-274` — No Session reuse; no 5xx retry; `_token_exchange` discards `expires_in` (see §2.6).
- `copernicus_dataspace.py:374-389, 386` — `requests.get(stream=True)` not in `with`; on `iter_content` raise the connection isn't released; `_safe_extract_zip` re-iterates extraction without applying its own validation results.
- `copernicus_dataspace.py:328` — Hard-coded 50-page pagination cap; silently truncates.
- `copernicus_dataspace.py:300, 417-422` — Stale signature accepting `access_key`/`secret_key` it never uses; mixing explicit keys + `auth_manager` silently ignores manager.

#### `python/siac/adapters/data/earthaccess_catalog.py`, `earthaccess_source.py`
- `earthaccess_catalog.py:188-194` — If `short_names` discovery returns empty, the configured short_name is silently accepted.
- `earthaccess_catalog.py:140-152` — `set_override` mutates cached instances; cross-catalog state leak.
- `earthaccess_source.py:91-99, 103-124` — `ea.login(**kwargs)` falls back to `ea.login()` on TypeError (silent kwarg drop); `_temporary_environment_credentials` is process-wide → race.
- `earthaccess_source.py:262-276` — `download_granules` returns `[]` on unknown shapes; generator outputs discarded silently.
- `earthaccess_source.py:264-268` — Hard-coded `_EARTHACCESS_DOWNLOAD_THREADS = 8`; not configurable.

#### `python/siac/adapters/data/gcs_sentinel2.py`
- `gcs_sentinel2.py:71-76` — `urllib.request.urlopen(req, timeout=timeout)`; empty/HTML body raises ValueError uncaught; no 5xx retry.
- `gcs_sentinel2.py:285-323, 357-359` — `except Exception` loses class info; `pool.map` with worker exception cancels in-flight downloads, partial files leak through `KeyboardInterrupt` path that bypasses cleanup.
- `gcs_sentinel2.py:296-323` — `attempts = retries + 1` then off-by-one in user-facing log.

#### `python/siac/adapters/data/water_mask.py`
- `water_mask.py:118-205` — Concurrent first-run race (mitigated by atomic move); no retry on transient HTTP. VRT companion `.tif` tiles re-create Session per call (perf).

#### `python/siac/adapters/satellite/base.py`
- `base.py:135-140` — `ThreadPoolExecutor(max_workers=2)` running `extract_geometry` + `extract_cloud_mask` while both implementations may mutate `self._reference_grid`, `_satellite_id` etc. → race.
- `base.py:380-413` — `read_text()` on Landsat MTL files with no size cap or encoding.

#### `python/siac/adapters/satellite/observation.py`
- `observation.py:38-41` — `try: rio.clip_box; except Exception: return field` returns un-clipped array silently.

#### `python/siac/adapters/satellite/sentinel2.py`
- `sentinel2.py:194-204` — `attrs[_TOA_BAND_LOADER_ATTR] = closure`; pickle/netcdf serialize will fail.
- `sentinel2.py:498, 507` — `ET.parse(granule_xml)` no DTD/size protection (XXE — moot for trusted ESA SAFE).
- `sentinel2.py:683-745` — XML namespace handling duplicated across helpers; **`_parse_view_angles` running-average bug** (see §1.1).
- `sentinel2.py:734-755` — Magic angle defaults (`5.0` / `100.0` / `30.0`) silently fabricated.
- `sentinel2.py:739-740` — `np.full((23,23), …)` defaults dtype to float64; later cast to float32.
- `sentinel2.py:770-828` — `_angles_to_grid` and `_angles_to_grid_batch` are unused. **Dead code**.
- `sentinel2.py:780-783, 850-851` — `np.linspace(...)` places coords at pixel edges, not centers — half-pixel bias on resample.
- `sentinel2.py:251-253` — `try: get_band(name); except: float("inf")` silently sorts missing-from-RSRF bands last.

### 3.4 Algorithms — RT / LUT / Emulator

#### `python/siac/algorithms/rt/direct/sixs.py`
- `sixs.py:47, 96` — Default `("xap","xbp","xcp")` quietly unioned with config; `preload_scene_subset` `*args, **kwargs` masks signature mismatches.

#### `python/siac/algorithms/rt/direct/sixs_build.py`
- `sixs_build.py:208-291` — Profile selection has no validation: misspelled profile silently selects O3 release. Stale cached extension may be returned for a different profile (`find_built_extension` runs first; profile path ignored when `module_path` is set).
- `sixs_build.py:237-269, 272-291` — Meson → distutils silent fallback; the two backends produce different ABI behaviour. Setuptools major-version regex sloppy on `60a1`.
- `sixs_build.py:594-602, 780-813, 866-876, 921, 1100-1192, 1408-1409, 1452-1499, 1539, 1576-1632, 1718-1748, 1822-1932, 1989-1990` — Many small Fortran patch fragility points: text-replace patches that depend on upstream ordering (modis vs ross-li-maignan); silent clamps (`taer55`, `xps`, `iwave=1` hard-set); textually fragile `\bwrite\s*\(\s*6\s*,\b` regex; non-atomic archive download with no checksum; `os.chdir(build_root)` is process-global (not thread-safe); `setup.py` invocation registers atexit handlers that linger; numpy.distutils import will break on numpy 2.x / Python 3.12.

#### `python/siac/algorithms/rt/direct/sixs_native.py`
- `sixs_native.py:43-45, 248-263` — RSRF support outside 6S grid silently falls back to delta; delta is not normalized when `np.trapezoid` returns 0.
- `sixs_native.py:325-435, 483-487` — Magic month-bucket atmosphere selection; `atmospheric_columns_mode` only honoured for `input_columns`.
- `sixs_native.py:621-748` — BRDF table layout reversed without justification; mixed `[]`/`.get()` access → late KeyError.
- `sixs_native.py:909-929` — `RegularGridInterpolator(... fill_value=NaN)` with `bounds_error=False`: out-of-LUT pixels silently NaN; user cannot tell.
- `sixs_native.py:817-893` — `_build_scene_lut_axis` may oscillate without convergence; `_build_scene_lut_plan`'s shrink-loop has no max-iter safeguard.
- `sixs_native.py:963-985` — Module cache stale-key cleanup not locked; first-time loads can double-import; module-name collisions in `sys.modules` evict previous extensions.
- `sixs_native.py:1086-1166, 1199-1206, 1290-1291, 1532-1563, 1571-1635` — Native output matrix in Fortran order assumed throughout; `__del__` during interpreter shutdown is hazardous; full-invalid-scene path skips RT silently with no warning; `_ensure_openmp_session` falls back to `isolate=False` with shared COMMON-block state across runners; no validation of `parallel_backend`; round-robin worker assignment uneven; `np.empty` outputs un-initialised if assignments mis-cover (currently safe).

#### `python/siac/algorithms/rt/emulator/two_nn.py`
- `two_nn.py:140` — `B8a` (lower-case) misses pattern.
- `two_nn.py:188-200` — `np.load(allow_pickle=True)` with only `is_relative_to` for path safety; symlinks under the dir bypass.
- `two_nn.py:233-235, 244-254, 267-273, 276-306` — Radians/degrees mix-up silently produces garbage; float32 cast hard-coded; output index assumptions on `[xap,xbp,xcp]` not asserted; Jacobian indices `[3]` / `[4]` for AOT/TCWV are hard-coded — should derive from `INPUT_FEATURES.index(...)`.
- `two_nn.py:407-438, 540-560, 572-576` — `load(satellite_id, _band_name, …)` ignores `_band_name`; `EmulatorRegistry.get_emulator` swallows everything to debug; `list_supported_sensors` hard-coded list incomplete.

#### `python/siac/algorithms/rt/lut/_spectral_math.py`
- See §1.2 finding 1.
- `_spectral_math.py:140-155` — `safe_*` clamps applied independently to denoms with opposite sign; internal inconsistency. Cast to float32 after division causes precision loss.
- `_spectral_math.py:165-202` — `finite_range(... fallback=(0,0))` collapses axis to a single point silently when no finite values.

#### `python/siac/algorithms/rt/lut/backend.py`
- `backend.py:91, 155-208` — `chunk_cache_size` accepted but never used (dead parameter). Zarr v3 detection re-tries v2 on failure; bare `except` masks config errors; `consolidated=True` silent fallback.
- `backend.py:188-189, 235-279` — Eager materialisation of all coords; `trans_total = np.maximum(td*tu, 1e-10)` clamps zero-transmittance pixels silently.
- `backend.py:357-359, 482-499, 501-511, 546-547` — Compact-coefficient LUT path uses argmin-nearest wavelength (not RSRF); ozone-units heuristic; clip + NaN-fill mismatch; altitude unit-coupling assumed.
- `backend.py:649-650, 705-737, 870-875, 947-1019` — `raa` abs-fold collapses 270° → 90° if input already `[0,360]`; `RegularGridInterpolator` rebuilt per call; `derive_standard_rt_coefficients` uses possibly-mismatched `rho1/rho2`; spectral cache key uses only `name/center/bandwidth` — RSRF differences collide; double-checked-locking has a subtle race in `_get_or_build_spectral_scene_subset`.
- `backend.py:1108-1124, 1140-1166, 1199-1203, 1264-1303, 1313-1357, 1380-1398, 1450-1453` — Linear backoff capped at 5s, no jitter; mean-over-altitude collapses non-linearity (see §1.2 finding 10); float32 casts on coords; `_reload_lut` clears caches without lock; transient-error detector substring-matches `"timeout"` (false positives); altitude correction is heuristic (scale_height=8.5 km hard-coded, sqrt for transmittance) with no citation; finite-difference Jacobian forward only; `is_available_for_sensor` always True.

#### `python/siac/algorithms/rt/lut/store.py`
- `store.py:127-153, 174-178, 194, 251-256, 269-298` — Reaches into `zip_fs._files`; atomic `.tmp` write good; `timeout` silently coerced to float; bare `except` rebuilds the cache silently.

#### `python/siac/algorithms/rt/lut/http_zip_store.py`
- `http_zip_store.py:80-122, 111-122, 151-157` — Unbounded `_size_cache` and `_full_body_cache`; multi-GB zips can blow out RAM. Whole-file caching should be opt-in.
- `http_zip_store.py:189-210` — `asyncio.Lock` created lazily and per event loop; multiple loops bypass serialization; race on lock-creation.
- `http_zip_store.py:219-291, 355-356, 419-433, 507-509, 541-547, 601-647` — EOCD/Zip64 parsing has subtle off-by-N risks with sentinel-valued fields and large local-header extra fields; close() doesn't reset `_files`; redundant filesystem-options validation; remote compressed zips fall back to fsspec which may materialise the whole file.

#### `python/siac/algorithms/rt/lut/rsrf_kernel.py`
- `rsrf_kernel.py:53-71, 91-94, 96` — When RSRF support is entirely outside the LUT, picks `[0,1]` silently (should raise); zero solar-irradiance silently falls back; sha1 hash with `usedforsecurity=False` (correct).

### 3.5 Algorithms — Surface / BRDF / Cloud / Grid / Correction / Solver

#### `python/siac/algorithms/solver/multigrid.py`
- See §1.1 (4, 6, 7) and §1.2 (7, 9) for highest-impact issues.
- `multigrid.py:236, 285, 290-293, 478, 497, 526-547, 782, 851, 932, 1066, 1183-1209, 1268, 1416-1473, 1538, 1571-1622, 1648-1666, 1680, 1769-1774, 1949-1972` — divisions by `mask.size` without zero-guard; magic uncertainty floors that disagree with `cost_config.min_*`; gap-fill seeds taken from too-loose mask; `ftol` near zero; spatial-mean uncertainty normalization; many float32→float64 round-trips; integer-arithmetic threshold off-by-one; `_resample_mask_to_grid` Python double-loop; non-monotonic level sizes for tiny scenes; lossy stride slicing drops geo metadata; redundant elevation resample per level; kernel resolution mismatch; objective gradient/cost mismatch; TCO3 in `fixed=("tco3",…)` silently ignored.

#### `python/siac/algorithms/solver/cost.py`
- See §1.2 finding 3 and §1.1 finding 8.
- `cost.py:120-174, 201-234, 282, 294-299, 387-405, 500-622` — `n_pixels` numpy scalar in slice; `_aot_da`/`_tcwv_da` mutated in-place per call (thread-unsafe); NaN priors recovered by floor-then-replace two-pass; `boa_unc` 3D vs 2D broadcast intent unclear; weights/sum NaN; clipping inside cost contradicts L-BFGS-B bounds (small Wolfe-search degradation); `zip` without `strict=True` truncates silently; once-only "2D BOA reused" warning may quietly degrade retrievals; Pseudo-Huber `delta == 0` divides by zero; `compute_laplacian_eigenvalues` and `create_sparse_laplacian` likely dead code — boundary-mask off-by-one.

#### `python/siac/algorithms/solver/aod_smoothing.py`
- `aod_smoothing.py:139, 198, 267, 328-344, 393-402` — `distance_transform_edt` failure mode for all-False mask; `last_surface` is "last" not "best"; Python `weights/sum` divide-by-zero impossible by construction but brittle; `_whittaker_smooth_axis` reshape couples to Rust binding shape; transient-NaN `surface[bad] = fallback[bad]` writes the *initial* nearest-seed value, breaking monotone convergence.

#### `python/siac/algorithms/surface/swir_refine.py`
- `swir_refine.py:447-453, 498-525` — Single-worker pool; `try/except Exception` then `raise`; `predict_visible_with_diagnostics` branch may leak the pool (no explicit shutdown).
- `swir_refine.py:516-554, 573, 632, 672, 790-800, 842, 864, 937, 1014, 1162-1212` — `predicted_source_fit` defaulted to zero short-circuits filtering; redundant finite-checks; hard-coded MODIS BRDF kernel hyperparams (`hb=2.0, br=1.0`); strict `np.allclose(rtol=0,atol=0)` may miss legitimate epsilon noise; nearest-resampled int16 is no-op rounding; **historical composites at *current* observation geometry**; `time` coordinate as integer index; missing-month silently falls through to full stack; native-resolution fallback `10.0` Sentinel-2-biased; `except (KeyError, RuntimeError): return None` swallows disk-read failures; AND-combine across bands is over-conservative.

#### `python/siac/algorithms/surface/spectral_mapping.py`
- `spectral_mapping.py:159-204, 288, 359, 543, 649-691, 765-770` — Name-mismatch alone triggers heavy mapping; common-window early-return for all-zero RSRF; tolerance `5e-3` magic; NaN unmasked through float32 matmul; `_cache_distance_metrics` is a no-op (dead method); float32 sentinel `-9999.0` for dedup loses precision and collides with real values; `np.unique(axis=0)` allocates a full copy.

#### `python/siac/algorithms/surface/monthly_composite_store.py`
- `monthly_composite_store.py:62, 115-142, 265-270, 323-324, 418, 572-575, 598-604, 937` — Bounds silently shifted to multiples of resolution; non-transactional manifest update; non-atomic period overwrite; magic block sizes; nodata `-1` colliding with real sentinel; `np.rint(NaN)` → undefined int16; transform-rotation rejected only at very tight tolerance; no log when grid inference fails.

#### `python/siac/algorithms/grid/assembler.py`
- `assembler.py:85-95, 106-185, 209-273, 303-338, 425-465, 540, 765, 828, 909` — Pixel grid origin may not match the observation grid (causing systematic offset); PIL bilinear/uniform_filter NaN-unsafe and aliasing for non-integer ratios; `(ImportError, ModuleNotFoundError)` duplicate; dynamic `__import__("contextlib")`; bare-ish `except` falling back to scipy on geo errors; conservative-True padding for masks; `_coerce_window_size` even-bumping; OpenCV `MORPH_ELLIPSE` discrete approximation differs from scipy; `_assume_aligned_native_grid=True` triggers the shape-only fast path even when template is coarser; `default_aerosol_solver_bands` may be empty (no guard); template stored as `(y,x)` is positional; stringly-typed `_siac_toa_band_loader` attribute lost on netCDF round-trip; post-assembly only validates TOA shape.

#### `python/siac/algorithms/cloud/mask.py` / `mapping.py` / `providers/omnicloudmask.py`
- `mask.py:23, 62-68, 163-166, 199, 401, 407-409` — `_EXPECTED_CLASS_VALUES` duplicated across files; non-dict cache silently replaced; hard-coded color windows; bilinear pre-resample of TOA before cloud detection couples reflectance scale; `Unsupported cloud provider` thrown for the only supported provider; `np.unique` O(N) for fast-path validation.
- `omnicloudmask.py:27-103` — `except Exception` masks partial-install issues; `NaN→0` then re-mask leaks fake "dark" context into convolutional predictions; magic class mapping dictionary undocumented.
- `mapping.py:75-81` — Default `0` (missing) used implicitly even when `unmapped_to_missing=False`; `np.isin` casts to `src.dtype` (lossy in pathological cases).

#### `python/siac/algorithms/correction/atmospheric.py`
- `atmospheric.py:46, 98, 117, 133` — `ValueError` instead of domain `siac.errors`; magic acceptance window; serial fallback duplicates iteration; coord misalignment masked by xarray-broadcast.

#### `python/siac/algorithms/surface/brdf_monthly_database.py`
- `brdf_monthly_database.py:36, 103, 230-270, 323-368` — `_feature_names_for_query_bands` mislabels non-NIR/SWIR triplets; `cKDTree(... copy_data=False)` lifetime coupling; IDW absolute eps in feature units; Python-loop exact-match path; `source_fit_rmse` defaulting to zeros silently disables filtering; trailing-name slicing fragile.

#### `python/siac/algorithms/surface/brdf_whittaker.py`
- `brdf_whittaker.py:120, 172-176, 227, 255` — Magic fallback `0.20`/`0.08`; `obs_time` consumed via `kwargs.pop` then unexpected-kwarg raise; `np.where(weights>0, ...)` NaN flow OK; days-only date math.

#### `python/siac/algorithms/surface/kernel_model.py`
- `kernel_model.py:67, 146, 189, 296, 411` — PSF FWHM comment misstates derivation; `1.5` BOA ceiling; bare `except`; `1.0e-10 + 0.01` threshold copy-paste artefact; threshold inconsistency between spatial and DCT paths.

#### `python/siac/algorithms/surface/prior_store.py`
- See §1.1 (9, 10).
- `prior_store.py:127, 266, 277-301` — Corrupt `.zattrs` silently included; MODIS-band fallback wavelengths leak into non-MODIS stores; collapses N bands → single-band BOA via `nanmean`; dummy zero-filled BRDF kernels propagate to uncertainty pipeline.

#### `python/siac/algorithms/surface/_swir_refine_resample.py` / `_swir_refine_months.py`
- `_swir_refine_resample.py:50, 110` — `np.array_equal` no tolerance; `assign_coords` silently overrides slight grid differences.
- `_swir_refine_months.py:31-43, 78-85` — Year-boundary correctness in the while-loops; convention "day 16" off by half-day from "day 15"; weekly-step heuristic counter-intuitive.

#### `python/siac/algorithms/brdf/kernels.py`
- `kernels.py:89, 184-227` — Vague shape-mismatch error; magic albedo coefficients with silent clamping (logger.warning at INFO suppressed); cubic terms zeroed via len-4 templates that should be simplified.

#### `python/siac/algorithms/__init__.py` / `surface/__init__.py`
- Re-exports of private helpers via public init paths.

### 3.6 App / Workflows

#### `python/siac/app/_assembly_correction.py`
- `_assembly_correction.py:23, 42-56, 75-89` — `_callable_supports_kwarg` returns implicit `None` swallowed; corrector re-instantiated each call; resample-every-field even when grids match; `corrected.aot/tcwv` silently overwritten from `solved.atmo_state`.

#### `python/siac/app/_assembly_preprocessor.py`, `_assembly_providers.py`
- `_assembly_preprocessor.py:42-58, 99-109` — Hard-coded color windows again; `get_metadata(...)` discarded for side effects; `setdefault` mutation of shared singleton.
- `_assembly_providers.py:38-355, 191-288` — Atmo factories duplicate `cast(...)` shape; `_CallableMonthlyCompositeProvider.source_bands` reassigned per-call (race-prone); `resolve_atmo_provider` returns a bound method, opaque ownership of provider state.

#### `python/siac/app/_assembly_solver.py`
- `_assembly_solver.py:35-66` — `isinstance(solver_config, dict)` branch dead; `inspect.signature` introspection silently drops `sharp_transition_mask`/`water_mask` if renamed; missing `result.atmo_state` falls back to prior with no log.

#### `python/siac/app/_assembly_surface.py`
- `_assembly_surface.py:31-66, 117-121, 241-501` — Hard-coded MSI/OLI band names belong in catalog; `mark_surface_prior_metadata` sets attrs on the function (lost when wrapped); database-vs-query geometry shared by reference; full monthly composite tuple held on `runtime.database` for the rest of the pipeline; `diagnostic_cache_dir` reaches into private attr; fallback factory may produce a second-tier failure with a less-informative error; `_build_monthly_surface_prior` registered in `SURFACE_PRIOR_METHOD_REGISTRY` but never retrieved (registration is dead code).

#### `python/siac/app/planning.py`, `registry.py`, `requests.py`, `s2_backend.py`, `sentinel2.py`
- `planning.py:138-176` — `or` short-circuit constructs every provider eagerly even when overrides are passed.
- `requests.py:51-52` — `Sentinel2SearchRequest` carries both `date` and `date_value` (partial migration).
- `s2_backend.py:33-35` — `_build_local_s2_backend` ignores `auth`.
- `sentinel2.py:42-96` — `S2Query(**query.__dict__)` rebuild loses dataclass safety; hard-coded `MSIL2A`/`MSIL1C` substring detection; `SIACConfig(sensor="s2")` synthesised without AOI may pick wrong env credentials.

#### `python/siac/workflows/__init__.py`, `_pipeline_*.py`, `pipeline.py`, `scene.py`, `sentinel2.py`
- `__init__.py:8-13` — `run_pipeline` re-exported as if public; ten parameters including private callable types.
- `_pipeline_config.py:38-138, 191-236` — `PipelineExecutionSettings` mixes dict + model access patterns; seeded defaults can shadow new required fields; duplicated `_select_solver_bands_for_preload` vs `_assembly_surface.select_surface_prior_bands`; expensive scatter diagnostics enabled by default when output config is absent.
- `_pipeline_diagnostics.py:48-146` — `build_aot_scatter_diagnostics` runs even when caller will not use the result; recomputes RT per band (~2× RT cost per scene); `_finite_diagnostic_field` fabricates infill in a *diagnostic* field; sort then linspace-sample logic unclear.
- `_pipeline_executors.py:113-167` — `LocalCluster(processes=False, threads_per_worker=1)` is functionally a ThreadPool with Dask overhead; LUT-preload `ThreadPoolExecutor` blocks shutdown on slow LUT reads; `LocalCluster` exits before `run_tail` (M4/M5/M6 run on calling thread despite `backend='dask'`).
- `_pipeline_outputs.py:42, 73` — Integer-band fallback names `band_01` mismatch sensor naming; `sample_index.astype(np.int16)` overflows silently for long series.
- `_pipeline_priors.py:96-204, 124-167` — Dask `Future.cancel()` without `force=True` doesn't actually stop work; LUT-preload timeout misses `cancel`; LUT preload competes for the same dask pool with M3.
- `monthly_composites.py:61-120` — Missing-period reporting INFO-only; deprecated GeoJSON CRS member used in AOI spec.
- `pipeline.py:46-758` — `as _X` re-exports indicate partial migration; `_stage_timeout` cannot disable a global timeout per stage; LUT-preload skip not observable; `_open_correction_output_stream` returns None silently when streaming is not supported; corrector signature introspection breaks if wrapped; `paths_config.water_mask = ""` falls through to a *remote URL* default; M4 stage validation runs *after* `observe_stage`; skip-correction branch lacks `validate_correction_result`; `solver_qa` divergence between paths; diagnostics failures invisible to monitoring; `output_stream.finish` runs inside `observe_stage("M6")`, conflating compute and write; `has_written` attribute is not a Protocol; mixed `execution`/`max_workers` args without warning; futures-timeout from M5 may escape as `FuturesTimeoutError` (logged as "error", not "timeout"); type-alias info lost across `_pipeline_*` boundaries; `_run_pipeline_thread` and `_run_pipeline_dask` near-duplicates; `_OUTPUTS_WRITTEN_METADATA_KEY` private symbol shared cross-file; M6 surface-prior re-emitted from M3 (`corrected.surface_prior` discarded); `_call_with_retries` cannot distinguish transient vs permanent OSError; `_set_rt_observation_time` silent no-op; monthly composite collection held in memory through M6 (large heap retention).
- `scene.py:23-108` — `derive_execution_report_path` runs even when `output_path is None`; reads private `_OUTPUTS_WRITTEN_METADATA_KEY`; `write_output` is a one-line passthrough not in module exports.
- `sentinel2.py:29-53` — Resolves config + auth twice (duplicate side effects); duplicates `request.config.sensor` coercion logic with `app/sentinel2.py:75`.

### 3.7 Storage / Geo / Catalog

#### `python/siac/storage/writers.py`
- See §1.3 finding 2.
- `writers.py:146-249, 312-427, 553-573, 957-968` — No atomic-write; missing GeoTIFF predictor / num_threads / `BIGTIFF=IF_SAFER`; COG kwargs leak to GTiff-only options; Zarr lacks consolidated/region/safe_chunks; NetCDF lacks per-variable chunking; coord encoding overwrite drops time `units`/`calendar`; RGB scale-factor branch logic looks inverted vs docstring; nodata + zero-reflectance collapse to same uint8 value; `_compute_overview_levels` dead.
- `writers.py:65-107, 294` — `ensure_writable_directory` chmods existing user dirs (privilege escalation risk on multi-user); skip rule misnamed.

#### `python/siac/storage/stac.py`
- `stac.py:208, 259, 295, 300, 316-318, 341, 373, 377, 383, 389-390, 540-554, 145-167, 152, 471-481` — FWHM ≠ bandwidth (semantic mislabel); redundant `.zarr` suffix check; `metadata.get("bounds", (0,0,1,1))` fabricates a default; bounds/geometry break across antimeridian; `view:sun_elevation` assumes radians input (no unit assertion); `view:off_nadir` no range validation; `processing:software` version hard-coded; `float(NaN)` fails json.dump; self link `href: "./"` invalid; eo extension always added even when no `eo:bands`; no root/parent links; group split by `key.split(".")[0]` collides on dotted keys.

#### `python/siac/storage/readers.py`
- See §2.8.
- `readers.py:95-107, 132, 175-178, 217-238, 282-303, 412-421, 443, 495-527` — Dataset leaks; only logs missing CRS; integer rasters cast to float32; same-CRS reproject for resample is heavy; remote-scheme detection misses `gs://`/`abfs://`; CRS string-equality compare; `crs_wkt` may be None.

#### `python/siac/storage/__init__.py`, `raster_writers.py`, `product_writers.py`
- Re-export shims; three import paths for the same function.

#### `python/siac/geo/geometry.py`
- `geometry.py:71-489, 511-570` — `Path(aoi_str).exists()` accepts directories; WKT detection misses MULTILINESTRING/GEOMETRYCOLLECTION/MULTIPOINT; `_normalize_geojson` silently drops all but first feature; deprecated GeoJSON `crs` member coupled to user intent; `_parse_wkt` discards `_target_crs`; `_simple_wkt_parse` lacks MULTIPOLYGON / Z/M / holes; `geometry.copy()` shallow; `_flatten_coordinates` IndexError on empty; `bounds_area` returns degrees² for geographic CRS; pixel-coords helpers assume axis-aligned transforms via `np.round`; return order `(row, col)` may surprise rasterio users.

#### `python/siac/geo/reprojection.py`
- `reprojection.py:88-179, 217-289, 484` — Unknown method silently downgrades to bilinear; integer rasters resampled with bilinear by default; `int(...)` truncation loses fractional border; same-CRS reproject for resample is heavy; `clip_box(*, crs=...)` rioxarray version coupling; `isinstance(geometry, list)` doesn't cover tuples; CRS object equality is not authority-based.

#### `python/siac/geo/resample.py`
- See §1.4 finding 1.
- `resample.py:113-285` — Implicit 2D assumption returns 0; bare `except`; integer-zoom over-shoot truncates without warning; mask interp swallows everything; non-2D mask returns ones silently; `xr.concat([], dim='param')` raises if `param` size is 0.

#### `python/siac/catalog/sensors/landsat.py`, `sentinel2.py`, `_common.py`, `registry.py`
- `landsat.py:9-20` — Only 7 bands; no thermal, pan, QA — feature-set mismatch with Sentinel-2 (13 bands). Default scale/offset only valid for Collection-2. RSRF naming mismatch (`B01` vs `B1`) latent.
- `sentinel2.py:1-71` — Three near-identical configs duplicated; FWHM values are approximations without source citation.
- `_common.py:6-26` — Empty `rsrf_band_id=""` mapped to `name`; `sensor_unit_id` not validated (empty allowed).
- `registry.py:11-27` — Closed registry; no plugin path; verbose error message.

### 3.8 Tools / Tests / Rust

See §2.9 and §2.10 for the consolidated lists. Additional per-file detail:

- `tools/aod_spatial_smoothing_experiment.py:18, :175` — `sys.path` mutation; default writes alongside input.
- `tools/bench_upsample.py:1-56` — No argparse, no `__main__` guard, unseeded RNG, uneven warm-up.
- `tools/benchmark_6s_routes.py:189-287` — Reaches into `backend._runner._openmp_session`, `_runner._prepare_scene_inputs`, `_runner._build_scene_lut_plan`, `_runner.close()` — tight coupling to private internals.
- `tools/run_full_s2.py:62` — Default cache under `$HOME` without confirmation.
- `tests/conftest.py` — Missing fixtures dir; `default_rng()` unseeded for "sample" fixtures.
- `tests/regression/__init__.py` — Empty marker package.
- `src/siac_rs/Cargo.toml:14-22, 32, 42-44` — Dead deps; PyO3 0.21+ migration is a tracked-but-unticketed blocker.
- `src/siac_rs/src/whittaker.rs:139-148` — 4-D index arithmetic per element; could batch.
- `src/siac_rs/src/optimization.rs:914-1163` — `refine_grid_search_with_qa` serial loop; rayon target.
- `src/siac_rs/src/psf.rs:56-62` — Function name and docstring claim DCT; implementation is direct spatial.

---

## 4. Suggested Roll-up Action Plan

A practical ordering of work items (rough effort tiers in parentheses):

**Tier A — Quick wins / bug fixes (≤ 1 day each)**
1. Repair broken tool imports (§1.1 #1) — every `from siac.config.schema import …` and every `from siac.app.assembly import …`.
2. Fix `cli.py:183-186` WKT route.
3. Register `siac_process` in top-level `_LAZY_IMPORTS` or remove its api re-export.
4. Fix `_parse_view_angles` running-average bug in `sentinel2.py:705-731`.
5. Set `multigrid.py:113` `ftol` to a usable value (or document why it's intentionally zero and remove the user-facing knob).
6. Close HDF4 SD handles in `earthdata_common.py:483, 507, 528`.
7. Wrap GeoTIFF/COG/NetCDF/Zarr/STAC writes in atomic-rename helpers.
8. Remove dead helpers identified above (`_compute_overview_levels`, `_cache_distance_metrics`, `_coerce_int_tuple`, dead `_angles_to_grid*`, `_load_temporal_payload_vrt[_from_*]`, `_build_monthly_surface_prior` registration, etc.).

**Tier B — Numerical correctness audits (1–3 days each)**
1. `cost.py:601-620` Laplacian boundary indexing — verify nx/ny convention or remove the dead path.
2. `multigrid.py:1648-1666` fixed-parameter cost/grad consistency — recompute cost as a *function of free vars only*.
3. `prior_store._crop_to_bounds` y-axis convention; DOY leap-year math.
4. RT path: backend.py `argmin`-nearest sampling vs RSRF integration; ozone heuristic; altitude correction citation; LUT subset mean-over-altitude.
5. `_spectral_math.py:121` sign-preserving denominator clamp.
6. `surface/swir_refine.py:842` historical composites' geometry assumption (science check).

**Tier C — Layering / typing / API hygiene (sustained refactor)**
1. Break `geo ↔ runtime`, `storage ↔ runtime`, `domain → geo` cycles via `siac.protocols` or by moving shared types into `domain/`.
2. Settle the `Resolved*` config story — either give every resolved class a meaningful subclass diff or alias them all and drop the indirection.
3. Replace enum-vs-string equality across `config/algorithms.py` with explicit enum members.
4. Eliminate `as _X` re-exports in `workflows/pipeline.py`; promote intentionally-public symbols to a real public module.
5. Tighten exception taxonomy: replace `ValueError` / `RuntimeError` with `siac.errors` types where appropriate; widen `siac.errors.ValidationError` to subclass `ValueError`.

**Tier D — Robustness and observability**
1. Single HTTP transport (probably `requests` with a shared `Session` + `urllib3.Retry`) used everywhere except where `fsspec` is a hard requirement.
2. Attach `_log_filter.SecretRedactionFilter` to the root logger, not just `auth`.
3. Replace `except Exception:` warnings with narrow types + `exc_info=True` and observability events.
4. Add timeouts/jitter to `_retry.py`.
5. Make remote LUT default opt-in (or make the default URL clearly point at a CDN with a pinned version).

**Tier E — Tests / Tooling**
1. Mark every test in `tests/integration/` with `@pytest.mark.integration`; ensure CI fast-suite excludes it.
2. Seed all `default_rng()` fixtures.
3. Either populate or remove `tests/regression/`.
4. Add a `test_tool_imports.py` that imports every `tools/*.py` to catch the broken imports listed in §1.1.
5. Audit `tests/unit/*` for direct imports of `_*` private helpers and add a contract layer for the helpers worth stabilising.

**Tier F — Rust crate**
1. Wrap parallel hot loops with `py.allow_threads`.
2. Parallelise `refine_grid_search_with_qa`.
3. Remove the no-op bench stub or implement real benches.
4. Drop unused deps; rename `FixedParameter::None`; align `psf::dct_convolve` name with the actual implementation (or implement DCT properly).

---

*Generated by an automated, parallel review across 8 module groups. Each finding cites a path:line reference. Many findings are "smell" rather than "smoking-gun" — verify before acting, especially the numerical ones.*
