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

*This report is generated, not curated. Trust but verify each row before acting on it.*
