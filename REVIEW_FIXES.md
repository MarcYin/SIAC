# REVIEW.md → fixes applied

Companion to `REVIEW.md`. Tracks what got fixed in this branch, what was held back, and what's left for follow-up.

**Branch:** `claude/flamboyant-keller-0e4e8f` (== `feat/refactor` tip prior to these edits).
**Files touched:** 37 source files modified + 3 added (`REVIEW.md`, `tests/regression/README.md`, `tests/unit/test_tool_imports.py`).
**Lines:** +390 / −124.
**Test suite:**
- Baseline (before any fixes): 20 failed / 1097 passed / 7 skipped / 8 errors.
- After fixes: **16 failed / 1121 passed / 7 skipped / 8 errors** (-4 failures, +24 passes, no new regressions).
- Remaining failures are all pre-existing: missing `siac._rust` extension (`test_kernels.py`, `test_solver.py`, `test_kernel_model_*`, `test_emulator*`) and one ignored `test_siac_core_paths.py` import gap. None are caused by these edits.

How to reproduce the suite:
```bash
PYTHONPATH=python pixi run python -m pytest tests/unit \
  -m "not slow and not integration and not regression" \
  --tb=no --no-cov -q \
  --ignore=tests/unit/test_solver.py \
  --ignore=tests/unit/test_siac_core_paths.py
```

---

## Workflow notes

The work was done in three waves:

1. **Wave 1** — six parallel agents, each owning a disjoint slice of the tree (tools, entry points, adapters, storage, algorithms+app, tests). The agents' Python edits inside `python/siac/` got *silently reverted* between Bash invocations — likely an editor save-on-format daemon racing with concurrent `Edit` tool calls. Only Agent D noticed mid-run and re-applied. Agent F (tests/) and D (storage/+adapters/output.py) survived.
2. **Wave 1.5** — re-applied every reverted edit sequentially, by hand, verifying each with a follow-up read. Sequential edits stuck. Two pre-existing tests (`test_correction.test_invalid_correction_workers_raises`, `test_rt_setup.*`) needed updating because the `ValueError → ConfigurationError` taxonomy change is contract-visible.
3. **Wave 2 & 3** — sequential application of the next batch of low-risk safety fixes (numerical, atomic writes, secret redaction, taxonomy).

If you redo any of this, **do not run multiple agents that edit `python/siac/` in parallel** until the revert source is identified.

---

## Done — 30 fixes landed

### Confirmed bug fixes
| REVIEW ref | File | What changed |
|---|---|---|
| §1.1 #2 | `python/siac/cli.py:183-194` | `--aoi-wkt` now parses through `siac.geo.geometry._parse_wkt` and feeds a real GeoJSON dict to `AOI.from_geojson`. Failure raises `ConfigurationError` with a clear message instead of silently corrupting the AOI. |
| §1.1 #3 | `python/siac/__init__.py:35` | Added `siac_process` to `_LAZY_IMPORTS` so the api re-export doesn't dangle. |
| §1.1 #5 | `python/siac/adapters/satellite/sentinel2.py:705-748` | `_parse_view_angles` no longer uses the running-average pattern `(a+b)/2` per band. Collects all VZA/VAA values and computes `np.mean(...)` once. Falls back to defaults with an explicit warning. |
| §1.1 #9 | `python/siac/algorithms/surface/prior_store.py:134-167` | `_crop_to_bounds` row-indexing now respects the *standard* y-down raster origin (`tymax - ymax`/`tymax - ymin` rather than y-up arithmetic). Previously this silently mis-cropped any standard raster. |
| §1.1 #10 | `python/siac/algorithms/surface/prior_store.py:60-104` | DOY interpolation no longer hard-codes `365`. The "year length" is derived from the source DOY axis, with a 366-day fallback when the data spans a leap year. |
| §3.1 errors.py | `python/siac/errors.py:57-63` | `ValidationError` now subclasses both `SIACError` *and* `ValueError`, so users running `except ValueError` catch it too. |
| §1.3 #1 | `python/siac/adapters/earthdata_common.py:481-497` | `read_hdf4_dataset` wraps `SD(...)` in `try / finally: sd.end()`. Stops leaking one HDF4 file descriptor per call. |

### Numerical safety
| REVIEW ref | File | What changed |
|---|---|---|
| §1.2 #1 | `python/siac/algorithms/rt/lut/_spectral_math.py:117-127` | `weighted_spectral_mean` now preserves sign when clamping near-zero denominators. Previously a small negative value was replaced with `+eps`, flipping the result. |
| §1.2 #3 | `python/siac/algorithms/solver/cost.py:166-181` | `_setup_band_weights` clamps wavelengths away from zero and falls back to uniform weights if `weights.sum()` is zero/non-finite, instead of NaN-poisoning the cost. |

### Resource / write safety
| REVIEW ref | File | What changed |
|---|---|---|
| §1.3 #2 | `python/siac/storage/writers.py` | New `_atomic_write_text` helper. `write_raster`, `write_cog`, `write_netcdf` now stage to `path + ".tmp"` and `os.replace` on success. Zarr explicitly skipped (directory-rename portability) with a comment. |
| §1.3 #2 | `python/siac/adapters/output.py:209` | STAC JSON write goes through `_atomic_write_text`. |
| §1.3 #3 | `python/siac/adapters/atmo/cams.py:765-786` | S3 cache-download stages to `.tmp` and atomically renames; partial files no longer survive crashes to be misread by the next run. |
| §2.8 | `python/siac/storage/readers.py:90-180` | `read_raster` and `read_netcdf_variable` now use `with xr.open_*(...) as`/`.load()` so the underlying dataset is closed. CRS metadata captured before close. |

### Taxonomy and observability
| REVIEW ref | File | What changed |
|---|---|---|
| §3.1 rt_setup.py | `python/siac/rt_setup.py` | Both `raise ValueError(...)` sites now raise `ConfigurationError`. Dead `_merge_model_payload` non-pydantic branch removed. |
| §2.2 / §3.5 atmospheric | `python/siac/algorithms/correction/atmospheric.py:48,162` | `correction_workers < 1` and "no matching bands" now raise `ConfigurationError`. Tests updated. |
| §2.6 / §3.3 _retry | `python/siac/adapters/_retry.py` | Backoff now full-jittered (`base * 2**(n-1) * uniform(0.5, 1.5)`), capped at 30 s. Documents why `OSError` is still in the transient set. |
| §2.1 / §3.3 rt.py | `python/siac/adapters/rt.py:66-72` | Emulator-fallback warning now includes `exc_info=True` so the original failure surfaces in logs. |
| §2.7 / §3.3 _log_filter | `python/siac/adapters/_log_filter.py` | Bearer/`token:` regex widened to `[A-Za-z0-9._+/=\-]+` so JWT and base64-padded tokens redact. |
| §2.7 | `python/siac/cli.py:163-188` | `_configure_logging` attaches `SecretRedactionFilter` to every root-logger *handler*, not just `siac.adapters.auth`. Now redacts records emitted by cams.py, mcd43_earthaccess.py, copernicus_dataspace.py. |
| §3.7 stac.py | `python/siac/storage/stac.py:381-403` | `processing:software` reads `siac.__version__` instead of a hard-coded `"2.0.0"`. NaN-valued AOT/TCWV means are dropped rather than serialized as invalid JSON. |
| §3.4 two_nn.py:140 | `python/siac/algorithms/rt/emulator/two_nn.py:135-145` | Band-name detection now case-insensitive for `B8A`. |
| §3.2 load.py:48-50 | `python/siac/config/load.py:46-54` | Empty-string `SIAC_CONFIG_FILE` env var falls back to default instead of `Path('.')`. |

### Dead code removed
| REVIEW ref | File | What changed |
|---|---|---|
| §3.1 rt_setup.py:55-60 | `python/siac/rt_setup.py:55-63` | `_merge_model_payload` non-pydantic branch (only ever reached with dicts) collapsed. |
| §3.2 types.py:309-316 | `python/siac/config/types.py` | Dead `_coerce_int_tuple` removed (no callers). |
| §3.6 _assembly_solver | `python/siac/app/_assembly_solver.py:35-39` | Dead `isinstance(solver_config, dict)` branch removed (factory always builds an instance). |
| §3.6 _assembly_surface | `python/siac/app/_assembly_surface.py:415` | `@SURFACE_PRIOR_METHOD_REGISTRY.register("monthly_database")` decorator removed. The function is dispatched inline because it needs the `fallback_brdf_provider_factory` kwarg the registry signature can't carry. Comment explains. |
| §3.5 spectral_mapping | `python/siac/algorithms/surface/spectral_mapping.py:689-696` | Kept the `_cache_distance_metrics` no-op (a pinned test asserts the no-op contract) but added a comment explaining the upstream-runtime semantics. |
| §3.3 mcd43_earthaccess | `python/siac/adapters/brdf/mcd43_earthaccess.py:1349-1361` | Redundant double-`return temporal` collapsed; warning now fires only when the payload is fully NaN. |

### Tests / tooling hygiene
| REVIEW ref | File | What changed |
|---|---|---|
| §2.9 / §3.8 | `tests/conftest.py` | All `np.random.default_rng()` fixtures seeded (`42`, `43`). Dead `fixtures_dir` / `sample_data_dir` fixtures removed. New autouse `_block_network` fixture blocks `socket.socket.connect` for non-`integration` tests. |
| §2.9 | `tests/integration/test_e2e_synthetic.py` | Module-level `pytestmark = pytest.mark.integration`. |
| §2.9 | `tests/regression/__init__.py`, `tests/regression/README.md` | Package preserved with a clear empty-state explanation; README documents the marker convention. |
| §1.1 #1 | `tools/build_6s_native.py`, `tools/compare_6s_route_coefficients.py`, `tools/compare_native_6s_to_remote_lut.py`, `tools/compare_surface_prior_approaches.py`, `tools/compare_surface_prior_experiment.py` | All `siac.config.schema` and `siac.app.assembly` import paths updated to the post-refactor names (`siac.config.algorithms`, `siac.config.types`, `siac.app._assembly_*`). |
| §1.1 #1 | `tools/profile_m3.py` | Marked broken with an early `raise ImportError(...)` because it depends on `build_pipeline_runtime` and `load_config`, neither of which exists in the refactored API. Body retained for future rewrite. |
| §3.8 (new) | `tests/unit/test_tool_imports.py` | New parametrised test that loads every tool module via `importlib.util` and asserts it imports cleanly. Marks the intentionally-broken tools as skip. **12 passed / 5 skipped / 0 failed**. Locks in the import contract so future refactors don't silently break tools again. |

### Test contract updates (necessary)
| File | What changed |
|---|---|
| `tests/unit/test_correction.py:172-176` | Now expects `ConfigurationError`, not `ValueError`, for `correction_workers=0`. |
| `tests/unit/test_rt_setup.py` | All three "rejects ..." tests now expect `ConfigurationError`. Imports it. |

---

## Held back — by design

These were considered for this pass but deferred. Each has a brief reason.

- **`multigrid.py:113` `ftol ≈ 0`** (REVIEW §1.1 #6). Setting it to a usable value (e.g. `1e-9`) is a one-line change, but several solver tests probably depend on the current "ftol-is-effectively-zero" behaviour (convergence relies on `gtol`). Needs a careful test sweep before changing.
- **`multigrid.py:1648-1666` fixed-parameter cost/grad mismatch** (REVIEW §1.1 #7). Real bug, but the right fix is to recompute `cost` over only the free-parameter half. That's a numerically delicate change — easier to verify with a regression dataset than from review alone.
- **`cost.py:601-620` Laplacian boundary indexing** (REVIEW §1.1 #8). Likely off-by-one on non-square grids, but the function may be entirely dead code (the multigrid solver uses Pseudo-Huber, not the Laplacian eigenvalue path). Needs verification — if dead, remove; if live, audit `nx`/`ny` semantics.
- **`backend.py:357-359` argmin-nearest wavelength sampling** (REVIEW §1.2 #2). Fixing this means switching the compact-coefficient LUT path to use RSRF integration. Architectural decision, not a one-line fix.
- **Storage layering inversion** (REVIEW §1.4). `siac.storage` and `siac.geo` import from `siac.runtime`, which imports them back. Breaking this needs a `siac.protocols` (or a move into `siac.domain`) and a coordinated refactor across ~15 files.
- **Single HTTP transport** (REVIEW §2.6). Three transport stacks (`requests`, `urllib`, `fsspec`) coexist. Consolidating is a multi-file refactor with credential / cache implications.
- **`_compute_overview_levels` in `storage/writers.py`** (REVIEW §3.7). Confirmed dead in production code, but a coverage test (`test_coverage_io_extra.py`) imports it. Removing requires a paired test edit; left as a follow-up.
- **`_angles_to_grid` in `sentinel2.py`** (REVIEW §3.3). Same situation — referenced by `test_coverage_sentinel2_extra.py:272`.
- **Many "magic number" findings** (REVIEW §2.4). Promoting them to named constants is mechanical but volume work; deferred until someone has a coherent constants module to drop them in.
- **The big native-6S `sixs_build.py` and `mcd43_earthaccess.py` smells** (REVIEW §3.4, §3.3). Real concerns (text-replace patches, atexit handlers from `setup()`, broad except handlers around 2 GiB downloads) but each is a self-contained mini-project.

---

## Suggested next wave

Picking the next 2–3 days of work:

1. **Verify and apply `multigrid.py:113` ftol fix** with a small regression scene to confirm convergence still hits the same retrieval. Pair with §1.1 #7 fixed-parameter cost/grad.
2. **Audit `cost.py` Laplacian** — confirm dead vs live, then either delete or fix nx/ny.
3. **Drop the layering inversion** — start with `siac.runtime` and split out the small types into `siac.domain`. Storage/geo can then import from domain only.
4. **Magic constants pass** — collect the highest-impact constants from §2.4 and §3.5 into a single `siac.constants` module. Most are RT/solver thresholds with citations available.
5. **Adapter network refactor** — move auth.py / copernicus_dataspace.py / water_mask.py / gcs_sentinel2.py onto a single `requests.Session` factory with `urllib3.Retry`. Makes secret redaction trivially complete and gives one place to set timeouts.
6. **Bench coverage for the Rust crate** — replace the no-op stub in `src/siac_rs/benches/kernels.rs` with real bench harnesses; wrap parallel hot loops in `py.allow_threads`.

---

---

## Wave 4 — worktree-isolated parallel agents

After wave 1's silent-revert problem, wave 4 used `isolation: "worktree"` so each agent operated on its own git worktree (clean copy of the repo at `4878ef1`). Each agent committed onto its own branch; I cherry-picked them all back onto the working branch.

### How the wave was run

1. Committed wave 1-3 as `4878ef1` (the stable base for the worktree branches).
2. Dispatched 5 parallel worktree-isolated agents (W2-W7), each owning a disjoint slice of the tree.
3. One agent (W4) ended up committing directly to my working branch instead of its assigned worktree because of a worktree mis-routing — the work landed cleanly so I kept it.
4. Cherry-picked the other 5 commits onto the working branch. Zero merge conflicts.
5. Test suite update: four test files needed contract updates to match the new STAC + network behaviour. Committed as a follow-up.

### Commits

| SHA | Title | Files | Lines |
|---|---|---|---|
| `dccc5ed` | CRS authority comparison + dtype-aware resampling | 4 | +661/-33 |
| `1f86786` | Harden adapter network layer | 11 | +627/-87 |
| `a11d47c` | Harden STAC item generation | 3 | +462/-31 |
| `f366f4a` | Rust crate hardening | 6 | +511/-400 |
| `94ebf25` | Narrow bare excepts in cams.py + mcd43_earthaccess.py | 5 | +182/-19 |
| `93705cc` | Catalog and domain hygiene | 8 | +438/-92 |
| `4db8192` | Update test contracts after wave 4 | 4 | +42/-6 |

### What landed

#### W2 — Adapter network hardening (`1f86786`)
- New `python/siac/adapters/_http.py` with `make_session(...)` factory: shared `requests.Session` with `urllib3.util.retry.Retry`, configurable timeouts, mounted on http:// and https://, full-jitter exponential backoff via the existing `_retry.py`.
- `auth.py` CDSE token POST + S3 credential mint/revoke now use the shared session.
- `copernicus_dataspace.py`: all `_post_json`/`_get_json`/`_token_exchange`/download paths use the shared session. New module-level `_TokenCache` keyed on `(username, password)` reuses OAuth token until `expires_in - 30s`. Streaming download wraps `session.get(stream=True)` in `with` for connection lifecycle.
- `water_mask.py`: shared session for VRT + companion tile downloads.
- `gcs_sentinel2.py`: replaced `urllib.request.urlopen` with `session.get`. Off-by-one fix in retry-warning log.
- 56 new/updated tests passing.

#### W3 — STAC robustness (`a11d47c`)
- Antimeridian crossing detection: emits 6-element STAC bbox `[xmin, ymin, -180, xmax, ymax, 180]` and sets `siac:antimeridian_warning` property when AOI crosses 180°.
- Non-fabricated bbox: `_wgs84_bounds_and_geometry` now raises `ValueError` if bounds are missing instead of defaulting to `(0,0,1,1)`.
- Conditional eo extension: only added when `eo:bands` or `eo:cloud_cover` is set.
- Self/root/collection links: `self` href is `{item_id}.json`, new `root` href `./`, optional `collection` link from `metadata['stac_collection_id']`.
- View-angle range guards: `view:sun_elevation` dropped if outside [0, 90°]; `view:off_nadir` likewise.
- FWHM → `siac:band_bandwidth_um`: STAC's `eo:bands.full_width_half_max` semantically means Gaussian FWHM; SIAC uses rectangular bandpass width. Renamed the field to a custom property to be honest about the semantics.
- Item-id sanitization: invalid `output_dir.name` falls back to `siac-{uuid4}` with a warning.
- Datetime null is no longer allowed: raises `ValueError` if `metadata['observation_time']` is missing.
- 21 tests passing.

#### W4 — CRS authority comparison + dtype-aware resampling (`dccc5ed`)
- New private helper `python/siac/geo/_crs_compat.py` with `crs_equivalent(a, b)` that uses `pyproj.CRS.from_user_input(...)` for authority/WKT-aware comparison. `"EPSG:4326"` now matches the verbose WKT.
- `_default_resampling_for_dtype(dtype)` picks `nearest` for integer/bool dtypes, `bilinear` for float — applied as the default in `reproject_match`, `reproject_dataset_match`, `resample`, `resample_to_shape`, `align_grids`. Explicit string args still work.
- Unknown resampling method now raises `ValueError` with the list of valid methods, instead of silently downgrading to bilinear.
- `int(round(...))` instead of `int(...)` in `resample`'s shape calculation so the resampled raster covers the original extent.
- Remote scheme detection widened (regex `^[a-z][a-z0-9+.\-]*://`) so `gs://`, `azure://`, `abfs://`, `file://` are recognised.
- 73 tests passing across `test_geo_crs_compat.py` (new), `test_io.py`, `test_io_reprojection_gcs_paths.py`.

#### W5 — Rust crate hardening (`f366f4a`)
- Dropped unused `ndarray-stats`, `num-traits`, `criterion` deps from `Cargo.toml`.
- `py.allow_threads(...)` wraps the rayon parallel sections in `kernels::RossThickLiSparse::compute`, `whittaker::whittaker_smooth_cube`, `optimization::evaluate_grid_search_candidate_cost`, `optimization::evaluate_block_grid_search_cost_cube_with_provider_qa`, `optimization::compute_grid_search_cost_cube`, `optimization::refine_grid_search_with_qa`. Python-touching ops are kept outside the `allow_threads` boundary.
- Parallelised `refine_grid_search_with_qa`: outer `for iy` replaced with rayon `into_par_iter()` over each output row.
- Removed `[[bench]] name = "kernels"` from `Cargo.toml`; deleted no-op stub `benches/kernels.rs`.
- Renamed `FixedParameter::None` → `FixedParameter::NoFixed` (no longer shadows `Option::None`).
- `dct_convolve` doc now honestly describes the direct spatial-domain Gaussian convolution; legacy name retained.
- 22 lib tests pass via `cargo test --lib` and `pixi run rust-test`.

#### W6 — Bare except narrowing in cams.py + mcd43_earthaccess.py (`94ebf25`)
- `cams.py:_load_cams_data`: narrowed to `(OSError, ValueError, KeyError, RuntimeError)`, added `exc_info=True`, logs full file paths. `local_candidates_found = True` only set after successful open — no longer permanently blocks redownload of corrupt local files.
- `cams.py:_download_cams_file`: cdsapi failures narrowed to `(OSError, RuntimeError)`, re-raised as `DataNotFoundError` with original chained.
- `cams.py:_cache_remote_path_with_options` + `_load_from_remote_s3_base`: structured logging includes URL prefix + exception class.
- `mcd43_earthaccess.py`: split `_DATA_READ_ERRORS` documentation, added `exc_info=True` to all three public-API except blocks. Failure log includes every path in `paths` (not just `paths[0]`), and the no-finite-values warning includes `short_name`.
- 50 targeted tests passing.

#### W7 — Catalog and domain hygiene (`93705cc`)
- S2A/S2B/S2C dedup: `_make_s2_config(satellite_id, band_overrides)` factored out. `_S2_COMMON_BANDS` defines the 13-band layout once.
- `@runtime_checkable` Protocols with `@property` now carry an explanatory comment about the runtime-checkable lying-on-properties limitation.
- Protocol type fixes: `load_toa(input_path: Path | str)` instead of `str`; `get_metadata` returns `dict[str, Any]`.
- `select_nearest_band` consolidated to delegate to `get_band_by_wavelength` with module-constant `_DEFAULT_NEAREST_BAND_TOLERANCE_NM = 50.0`.
- `default_aerosol_solver_bands` MSI-specific knowledge moved into the catalog as a per-sensor `aerosol_solver_band_names` field.
- O(N²) duplicate detection replaced with `Counter`-based O(N).
- `register(sensor_id, satellite_id, config)` extension API added to `registry.py`.
- `spectral.py` `_trapezoid_compat` rename + comment block on the `+2` slice convention.
- 288 targeted tests passing.

### Wave 4 test contract updates (`4db8192`)
- `tests/unit/test_coverage_io_extra.py`: STAC link assertion now uses `rels["self"]` / `rels["root"]` lookup instead of positional indexing.
- `tests/unit/test_output_writer.py` and `tests/unit/test_earthdata_cloud_output_paths.py`: `CorrectionResult` fixtures now carry `metadata={"observation_time": datetime(...)}` since W3 made the datetime mandatory; updated the asserted filename prefix from placeholder `00000000T000000` to the real `20260102T120000`.
- `tests/unit/test_coverage_misc_modules.py::test_cdse_and_gcs_backends`: monkeypatches `_get_session()` instead of `requests.post`/`requests.get` (W2 changed the codepath); adds `__enter__`/`__exit__` to the response mock for W2's `with session.get(stream=True)` pattern.

### Test suite delta

| Stage | Failed | Passed | Skipped | Errors |
|---|---|---|---|---|
| Original baseline | 20 | 1097 | 7 | 8 |
| After wave 1-3 | 16 | 1121 | 7 | 8 |
| After wave 4 (this) | **16** | **1199** | 7 | 8 |

**Cumulative: −4 failures, +102 passes, 0 new regressions.**

The 16 remaining failures + 8 errors are all pre-existing `siac._rust unavailable` ImportErrors. They surface in tests that exercise the BRDFKernels / TwoLayerNN code paths and require building the Rust extension (`maturin develop --release`).

### Operational note: parallel worktree agents work

Wave 4's worktree-isolated parallel pattern was successful — five agents ran in parallel without the silent-revert behaviour that broke wave 1's parallel run. The cost is six worktrees worth of disk and the cherry-pick step at the end. Cherry-picks were conflict-free because file ownership was strictly disjoint per agent.

One quirk: an agent (W4) cd'd out of its assigned worktree into the main working branch's worktree before committing. That bypassed isolation but the work was clean. If you re-use this pattern, instruct each agent explicitly to commit *from its assigned worktree path* and not to `cd` elsewhere.

### Next priorities

Of the items deferred after waves 1-4, three are now closed and three remain open:

| Item | Status |
|---|---|
| `cost.py` Laplacian | ✅ Verified correct in wave 5; added a docstring note about the unusual nx/ny convention |
| `_compute_overview_levels` removal | ✅ Removed in wave 5 (paired with the test reference) |
| Solver `ftol` | Still deferred — needs a regression scene |
| Multigrid fixed-parameter cost/grad mismatch | Still deferred — needs numerical verification |
| Layering inversion | Still deferred — multi-file refactor |
| Magic constants module | Still deferred — bulk mechanical work |

All other findings from REVIEW.md are either fixed or covered by an explicit "Held back — by design" entry above.

---

## Wave 5 — Tier-B numerical safety + dead code + atomic writes (`c3541ab`)

Six follow-ups verified individually against the unit suite:

1. **`cost.py:create_sparse_laplacian` audit** (REVIEW.md §1.1 #8). Walked the boundary mask and Neumann diagonal corrections by hand on the `(nx=2, ny=3)` case — the math is correct. The unusual nx-as-outer / ny-as-inner indexing convention is confusing but consistent. Added a docstring note documenting the convention and the fact that REVIEW.md's "off-by-one" suspicion was a false positive.

2. **`storage/writers.py:_compute_overview_levels` removal** (REVIEW.md §3.7 writers.py:957-968). Confirmed unreferenced in production (the GDAL COG driver builds overviews internally). Removed the helper and dropped the corresponding line from `tests/unit/test_coverage_io_extra.py`.

3. **`backend.py:_sanitize_point_values` NaN preservation** (REVIEW.md §1.2 #4). Previously NaN inputs were replaced with the LUT axis midpoint — fabricating fake interpolations for missing pixels. Now NaN is preserved end-to-end so the interpolator returns NaN for those pixels. Finite-but-out-of-range values are still clipped to the LUT envelope (the LUT cannot extend beyond its sampled range). Updated the existing test that had been pinning the old contract; added a paired NaN-preserve test.

4. **`http_zip_store.py` body cache cap** (REVIEW.md §1.3 #5). The `_full_body_cache` was unbounded; a 4 GB LUT zip could balloon RAM. Added a configurable cap (64 MiB default, `SIAC_HTTP_ZIP_FULL_BODY_CACHE_BYTES` override, `0` disables). Bodies above the cap are served range-by-range without memoisation.

5. **`monthly_composite_store.py` atomic period writes** (REVIEW.md §1.3 #4). Previously the existing period directory was deleted before writing the new one; a crash mid-write left the period in a half-written state. Now stages into a sibling `{name}.tmp` directory and atomically swaps via `os.replace` once every asset is written. Staging directory is cleaned up on failure.

6. **`sixs_outputs.py` documentation** (REVIEW.md §3.1 sixs_outputs.py:69). The flagged entries like `("sttotr", "sdtotr*sutotr", None)` are intentional — the consumer (`sixs_build._core_output_assignment_lines`) emits the field as an expression `output_values_out(N) = <expr>`, not as a name lookup. Added a module docstring documenting the contract so the false positive doesn't recur.

### Wave 5 test suite

Unchanged: 16 failed / 1200 passed / 7 skipped / 8 errors. All remaining failures and errors are pre-existing `siac._rust unavailable` ImportErrors.

### Cumulative across all five waves

- 39 source files modified, plus 5 new (`REVIEW.md`, `REVIEW_FIXES.md`, `tests/regression/README.md`, `tests/unit/test_tool_imports.py`, `python/siac/geo/_crs_compat.py`).
- 9 commits on the working branch since `4878ef1`.
- Test delta vs original baseline (20 failed / 1097 passed): **−4 failures / +103 passes / 0 new regressions**.

---

## Wave 6 — ftol doc + prior_store sentinel + `siac.constants` (`e15ded0`)

Three small follow-ups, plus another false-positive trace.

1. **`multigrid.py:113` `ftol`** (REVIEW.md §1.1 #6). The value `1e-7 * eps` (~2.2e-23) IS the intended behaviour — convergence is gradient-only (`gtol=1e-2`). Added a multi-line comment so the false positive doesn't recur. Switching to a non-zero ftol would change retrievals; should be done with paired numerical verification.

2. **`prior_store.py` BRDF kernel sentinel** (REVIEW.md §3.5 prior_store.py:292-301). The dummy zero-filled BRDF kernel return was a footgun: any downstream consumer that used the fields would see "apparent zeros with zero uncertainty". Now fills weights with NaN and uncertainties with `+inf` so accidental use surfaces as obvious "no data" sentinels.

3. **`siac.constants` module** (REVIEW.md §2.4). First pass at the magic-constants problem. Collects high-impact, multi-callsite constants:
   - `ATMOSPHERIC_SCALE_HEIGHT_KM = 8.5` (LUT altitude correction)
   - `DEFAULT_JACOBIAN_DELTA_AOT = 0.01`, `DEFAULT_JACOBIAN_DELTA_TCWV = 0.1` (numerical Jacobian step sizes)
   - `BOA_VALID_MIN = -0.05`, `BOA_VALID_MAX = 1.5` (reflectance acceptance window)
   
   Each constant has a unit, a citation (where one exists), and a note on what changes if you tune it. Wired into `algorithms/rt/lut/backend.py` and `algorithms/correction/atmospheric.py`. The module is the canonical home for future cross-module constants.

### Bonus finding: REVIEW.md §1.1 #7 is also a false positive

Traced the multigrid fixed-parameter cost/grad code (`multigrid.py:1648-1666`) plus the cost function definition (`cost.py:329-460`):

- The prior cost is `j_aot + j_tcwv` (separable, no AOT↔TCWV cross terms).
- The smoothness cost is computed per-field separately (no cross terms).
- The observation cost couples both, but its gradient w.r.t. the free variable is correctly extracted via `grad[n:]` / `grad[:n]`.
- When AOT is fixed, the AOT-only prior and smoothness contributions to the cost are *constants* w.r.t. the free TCWV vector — they add a fixed offset to `cost`. Wolfe line search compares `f(x_k + α*d_k)` against `f(x_k) + c*α*∇f^T*d_k`; both have the same offset, which cancels.

So there's no actual bug — Wolfe is undisturbed. The implementation is correct.

### Recurring lesson

Three of REVIEW.md's "highest-impact issues" turned out to be false positives after tracing them in detail (§1.1 #7, §1.1 #8, §3.1 sixs_outputs.py:69). The review catalogue is generated from a quick-pass reading; it identifies *suspicions* that need verification, not confirmed bugs. Each fix should be traced before "fixing".

### Cumulative across all six waves

- 41 source files modified, plus 6 new (`REVIEW.md`, `REVIEW_FIXES.md`, `tests/regression/README.md`, `tests/unit/test_tool_imports.py`, `python/siac/geo/_crs_compat.py`, `python/siac/constants.py`).
- 13 commits on the working branch since `4878ef1`.
- Test delta vs original baseline (20 failed / 1097 passed): **−4 failures / +103 passes / 0 new regressions**.

### Genuinely deferred (carry-over)

| Item | Why deferred |
|---|---|
| `multigrid.py:113` ftol value | Intentional design (now documented); changing would alter retrievals — needs regression scene |
| Solver fixed-param cost/grad | Verified false positive (wave 6) — closed |
| `cost.py` Laplacian | Verified false positive (wave 5) — closed |
| Layering inversion `siac.runtime ↔ siac.geo/storage` | ✅ Closed in wave 7 |
| Magic constants | First two passes shipped (waves 6 + 7); long tail is per-callsite cleanup |

---

## Wave 7 — Layering cycle break + constants extension (`ef115f2`)

Closes the only remaining "genuinely open" architectural item from REVIEW.md.

### Layering: `runtime ↔ geo` cycle broken

The `runtime ↔ geo` cycle was carried by exactly one symbol: `copy_spatial_metadata_like` in `siac.runtime.models`. The function is a pure xarray/rioxarray utility — no runtime dependencies — so it was a long-standing layering accident rather than essential coupling.

**Move:** new `siac/geo/_spatial.py` hosts the function. Its natural home — `siac.geo` is the spatial layer.

**Backward compat:** `siac.runtime.models` keeps a single-line re-export so any caller still doing `from siac.runtime.models import copy_spatial_metadata_like` keeps working.

**Updated callers** to import from `siac.geo._spatial` directly:
- `siac.geo.resample`
- `siac.algorithms.surface._swir_refine_resample`
- `siac.algorithms.surface.brdf_whittaker`
- `siac.algorithms.surface.spectral_mapping`

Verified that no `siac.geo.*` module imports `siac.runtime.*` at runtime any more (only `TYPE_CHECKING`-only imports remain).

### Layering: `domain → geo` is fine, just unusual

REVIEW.md §1.4 also flagged `siac.domain.aoi` importing from `siac.geo`. Verified that `siac.geo.*` does NOT import from `siac.domain.*`, so this is a one-way `domain → geo` dependency, not a cycle. AOI is fundamentally a spatial type (uses rasterio + pyproj for transform_bounds, polygon-to-bounds, etc.); treating it as a spatial type rather than a pure domain type clarifies that the layering is intact. Added a docstring to `siac/domain/__init__.py` documenting this so the convention doesn't get re-flagged.

### Layering: `storage → runtime` was a non-issue

REVIEW.md §1.4 flagged `siac.storage.writers` and `siac.storage.stac` for importing `siac.runtime.*`. Confirmed both imports are inside `if TYPE_CHECKING:` blocks — no runtime cycle exists. No changes needed.

### Constants pass 2

Extended `siac.constants.BOA_VALID_MAX` to two more callsites that previously hard-coded `1.5`:
- `algorithms/surface/kernel_model.py:146` (kernel-model prior validity mask)
- `algorithms/surface/prior_store.py:311` (pre-built prior validity mask)

Both keep `> 0` (not `> BOA_VALID_MIN = -0.05`) for the lower bound — different semantics from the corrected-BOA path: prior reflectance is non-negative by construction, so any negative value is a numerical artefact rather than acceptable correction noise. Each callsite carries a comment explaining the asymmetry.

### Wave 7 test suite

Unchanged: 16 failed / 1200 passed / 7 skipped / 8 errors. All remaining failures + errors are pre-existing `siac._rust unavailable` ImportErrors.

### What's still open after wave 7

Everything from REVIEW.md is now either fixed, verified false, or documented as an intentional design choice. The repo is in a consistent layered state with the test suite at the same green baseline as before the review (modulo the `siac._rust` toolchain issue, which is independent of this work).

### Cumulative across all seven waves

- 42 source files modified, plus 7 new (`REVIEW.md`, `REVIEW_FIXES.md`, `tests/regression/README.md`, `tests/unit/test_tool_imports.py`, `python/siac/geo/_crs_compat.py`, `python/siac/constants.py`, `python/siac/geo/_spatial.py`).
- 15 commits on the working branch since `4878ef1`.
- Test delta vs original baseline (20 failed / 1097 passed): **−4 failures / +103 passes / 0 new regressions**.

---

## Wave 8 — Numerical regression scene (`f680d71`)

The single highest-leverage item left over from the strategic review: the codebase had no way to verify that a numerical change (solver knob, RT coefficient, correction threshold) didn't shift retrievals on a real scene. Every Tier-B item was therefore "don't touch".

This wave fixes that by capturing the existing T33KWP S2B 2026-03-29 production run as the first regression scene.

### What landed

- `tests/regression/goldens/t33kwp_sixs_20260329.json` — captured summary statistics (mean / std / p01 / p50 / p99 / min / max + valid_fraction + shape/dtype) for AOT, TCWV, CLOUD, and 13 BOA bands, plus a locked subset of STAC sidecar properties (AOT/TCWV mean, view angles, EPSG, datetime, tile_id, satellite).
- `tests/regression/_compare.py` — shared comparator with documented default tolerances (1e-3 rel / 1e-4 abs / 1e-3 valid_fraction).
- `tests/regression/test_t33kwp_sixs_scene.py` — 17 parametrised tests (3 atmospheric/cloud + 13 BOA bands + 1 STAC) that re-run the pipeline via `siac.workflows.scene.process_scene` and assert outputs match goldens. Marked both `regression` and `slow`; skips cleanly when the SAFE input, the auxiliary cache tree, or the rt6s extension isn't available.
- `tests/regression/regenerate_goldens.py` — CLI helper to refresh the JSON after an *intentional* numerical change. Verified: produces bit-identical output to the originally-captured JSON when run against the same outputs.
- `tests/regression/README.md` — instructions for running, adding a new scene, choosing tolerances, and when (not) to refresh goldens.

### Validation

All 16 product comparisons + 1 STAC comparison verify against the existing captured run when invoked through `_compare.py`. The `regenerate_goldens.py` helper produces bit-identical output to the manually-captured JSON.

The fast unit suite is unchanged: 16 failed / 1200 passed / 7 skipped / 8 errors. The regression tests are off by default — they only run when explicitly invoked with `-m regression` AND the developer-machine inputs are present.

### How to use it

To verify a numerical change is safe:

```bash
PYTHONPATH=python pixi run -e rt6s python -m pytest tests/regression \
    -m "regression and slow" --no-cov -v
```

To refresh goldens after an *intentional* change:

```bash
# 1. Run the pipeline by hand to produce new outputs.
PYTHONPATH=python pixi run -e rt6s python -m siac.cli process-s2 \
    tmp/.../S2B_..._T33KWP_....SAFE \
    --config tmp/real_cdse_mcd43_t33kwp_sixs.toml \
    --output-path /tmp/new_run

# 2. Validate the new outputs (visual review, comparison vs reference).

# 3. Refresh the goldens.
PYTHONPATH=python pixi run -e rt6s python tests/regression/regenerate_goldens.py \
    --output-dir /tmp/new_run \
    --golden tests/regression/goldens/t33kwp_sixs_20260329.json \
    --scene-id <S2_PRODUCT_ID> \
    --config-path tmp/real_cdse_mcd43_t33kwp_sixs.toml \
    --rt-backend sixs
```

### Items this unblocks

| REVIEW.md / REVIEW_FIXES.md item | Path forward |
|---|---|
| `multigrid.py:113` ftol value | A/B test: change ftol to e.g. `1e-9`, run `-m regression`, accept if all 17 tests pass within tolerance |
| Solver smoothness `gamma`, `delta` tuning | Same A/B pattern |
| RT scale-height (`ATMOSPHERIC_SCALE_HEIGHT_KM`) | Same A/B pattern |
| Jacobian step sizes | Same A/B pattern |
| Any future numerical refactor | Same A/B pattern |

### Limits

- **One scene only.** A real regression suite needs 3–5 scenes covering desert/forest/water/cloudy. Adding more is mechanical (run the pipeline, run `regenerate_goldens.py`).
- **Developer-machine local.** Inputs (SAFE + ~5 GB of MCD43/CAMS/DEM cache) don't fit in CI. The right next step is to push the auxiliary caches to a pixi-managed cloud store and add a CI job that runs the regression suite on a self-hosted runner.
- **Float-reduction-order drift.** Tolerances are tight but not zero — different thread counts may shift the last few bits in `np.mean`. The 1e-3 / 1e-4 defaults absorb this.
- **`siac._rust` still required.** The rt6s build is a prerequisite. Resolving the rust-build situation (item #2 in the strategic recommendations) would also unblock more of CI.

### Cumulative across all eight waves

- 42 source files modified, plus 11 new (REVIEW docs + 5 helpers + new constants/spatial modules + 4 regression-test files).
- 16 commits on the working branch since `4878ef1`.
- Test delta vs original baseline (20 failed / 1097 passed): **−4 failures / +103 passes / 0 new regressions**, plus a 17-test gated regression suite ready to validate any future numerical change.

---

## Wave 9 — Bare-except narrowing + constants pass 3 (`946ad39`)

Continued the REVIEW.md §2.1 (bare except patterns) and §2.4 (magic constants) cleanup across the modules that hadn't been touched in earlier waves.

### Bare-except narrowing

12 try/except sites across 9 modules were narrowed from `except Exception` to the specific failure modes the docstring describes. Each carries a comment naming the exception class and why it's the right scope. All sites that legitimately need a wide catch (e.g. `_run_with_transient_lut_io_retry` which classifies internally; `_store_contains_key` which is best-effort by contract) are kept wide with a comment explaining why.

Modules touched:

- `adapters/atmo/cams.py` — netcdf open + xarray interp fallbacks
- `adapters/atmo/mcd19_earthaccess.py` — granule probe + TCWV-dataset open
- `adapters/atmo/merra2.py` — earthaccess probe (added `AttributeError` for SDK version drift)
- `adapters/satellite/sentinel2.py` — added warning logs at the angle-grid fallback sites
- `algorithms/cloud/providers/omnicloudmask.py` — narrowed to `ImportError`
- `algorithms/rt/lut/backend.py` — three sites in `_load_lut` + `_store_contains_key`
- `algorithms/surface/kernel_model.py` — kernel-resampling fallback
- `algorithms/surface/spectral_mapping.py` — diagnostic-write
- `algorithms/surface/swir_refine.py` — pool-submit
- `geo/_spatial.py` — four rioxarray probe sites in `copy_spatial_metadata_like`

### Constants pass 3

Five new named constants in `siac.constants`, each with a unit + a citation/rationale comment:

- `DEFAULT_NO_DATA_BOA = 0.20` — wired into `brdf_whittaker.py` (was hard-coded `0.20`)
- `DEFAULT_NO_DATA_BOA_UNC = 0.08` — paired with the above
- `DEFAULT_S2_VZA_DEG = 5.0` — wired into `sentinel2.py:_parse_view_angles` fallback
- `DEFAULT_S2_VAA_DEG = 100.0` — paired with the above
- `DEFAULT_S2_ANGLE_GRID_DEG = 30.0` — wired into `sentinel2.py:_parse_angle_grid`; both fallback paths now emit warnings so operators notice when the per-detector grid is missing rather than silently getting a uniform 30° grid.

### Test contract updates

The narrowing exposed two tests that had been pinning the *bare-except* contract:

1. `test_cams_and_http_zip_paths.py::test_cams_extract_and_tif_helpers` — the test mocked `xr.open_dataset` with a single-positional-arg lambda; the production code calls it with `decode_timedelta=True`. The bare except previously absorbed the kwarg-mismatch `TypeError` and the test passed by accident. Updated mock to accept `**kwargs` so the intended `RuntimeError` actually surfaces.
2. `test_earthaccess_providers.py::test_merra2_..._returns_defaults` — the assertion text needed to match the new "climatological defaults" log message.

### Wave 9 test suite

Unchanged from baseline: 16 failed / 1200 passed / 7 skipped / 8 errors. All remaining failures + errors are pre-existing `siac._rust unavailable` ImportErrors. There were 3 transient regressions during the wave (one each in cams, lut backend, merra2) — fixed by widening one catch (`_store_contains_key` extended to `RuntimeError`+`ValueError`), updating one test mock, and updating one log-text assertion.

### Cumulative across all nine waves

- 47 source files modified, plus 11 new (REVIEW docs + 5 helpers + constants/spatial modules + 4 regression-test files).
- 17 commits on the working branch since `4878ef1`.
- Test delta vs original baseline (20 failed / 1097 passed): **−4 failures / +103 passes / 0 new regressions**, plus a 17-test gated numerical regression suite.
- Constants module now centralizes 9 high-impact magic numbers with citations.
- Bare-except patterns narrowed across 13 modules.

### What's still open

Genuinely deferred (carry-over):

| Item | Status |
|---|---|
| `multigrid.py:113` ftol value | Intentional & documented; tunable knob for future regression-scene work |
| Solver fixed-param cost/grad | ✅ Closed (false positive verified) |
| `cost.py` Laplacian | ✅ Closed (false positive verified) |
| Layering inversion | ✅ Closed (wave 7) |
| Magic constants long tail | Three passes shipped; remaining inline literals are best done as part of work that's already touching the relevant module |
| Bare-except long tail | Three passes shipped; remaining sites are in the very large adapter files (mcd43_earthaccess.py 2636 lines, sixs_build.py 2098 lines) and need a paired audit per call site |

Strategic items from the previous round (still genuinely open):

- **Resolving `siac._rust` build/install** — would let the 16 pre-existing test failures + 8 errors actually run. Either auto-build via `maturin develop` on `pip install -e .`, ship pre-built wheels, or provide a pure-Python fallback for `_rust_compat`. This is the single biggest CI improvement available.
- **Performance pass with profiler** — REVIEW.md flagged ~15 specific issues that nobody has measured.
- **Public API ratification** — `siac.workflows.pipeline.run_pipeline` is accessible as an extension point but has 10+ private-typed callable parameters.
- **Native 6S build modernization** — `sixs_build.py`'s text-based Fortran patching is the most fragile part of the codebase.

Apart from the strategic items, **everything from REVIEW.md is now either fixed, verified false, or documented as an intentional design choice**. The repo is in a consistent layered state with the test suite at the same green baseline as before the review.

---

*This report is generated, not curated. Trust but verify each row before acting on it.*
