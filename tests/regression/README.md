# SIAC numerical regression tests

This directory hosts end-to-end pipeline tests that re-run SIAC on a real
Sentinel-2 scene and assert that AOT, TCWV, BOA, and STAC sidecar outputs
match a captured set of summary statistics within tight tolerances.

Their purpose is to make solver / RT / correction *numerical* changes
safe to apply. Without a regression scene every numerical knob is
"don't touch" — with one, every change becomes a small, verifiable diff
against a golden snapshot.

## Markers

Every test in this directory carries both:

- `@pytest.mark.regression` — opt-in via `-m regression`. The default
  test invocation in `pyproject.toml` excludes regression tests.
- `@pytest.mark.slow` — paired with `regression` because the pipeline
  takes ~3 minutes on a workstation.

## Running

These tests are not run in CI by default. Locally:

    pixi run -e rt6s python -m pytest tests/regression \
        -m "regression and slow" --no-cov -v

To run a single scene:

    pixi run -e rt6s python -m pytest \
        tests/regression/test_t33kwp_sixs_scene.py \
        -m "regression and slow" --no-cov -v

## What's needed locally

The tests skip cleanly (rather than failing) when their inputs aren't
present. To make a scene runnable you need:

1. **The SAFE input** at the path encoded in the config TOML, or
   pointed at by `SIAC_REGRESSION_SAFE`.
2. **The auxiliary caches** (MCD43, CAMS, DEM, BRDF) the config
   references — typically a multi-GB tree under `tmp/`.
3. **The native 6S extension** (`rt6s` pixi env, built via
   `pixi run -e rt6s build-6s-native`).

If any is missing the relevant fixture skips with a message explaining
which knob to set.

## Adding a new scene

1. Run the pipeline by hand against the new scene + config and save the
   outputs somewhere stable.
2. Generate goldens with the helper next to this README:

       PYTHONPATH=python pixi run -e rt6s python \
           tests/regression/regenerate_goldens.py \
           --output-dir /path/to/run/outputs \
           --golden tests/regression/goldens/<scene_slug>.json \
           --scene-id <S2_PRODUCT_ID> \
           --config-path tmp/<config>.toml \
           --rt-backend sixs

3. Copy `test_t33kwp_sixs_scene.py` to `test_<scene_slug>.py` and
   update `GOLDEN_PATH`, `DEFAULT_CONFIG`, `DEFAULT_SAFE`, and the
   parametrised band list to match your scene.

## Tolerances

`_compare.py` defines defaults:

- `DEFAULT_REL_TOL = 1e-3` — relative tolerance on means / stds /
  percentiles.
- `DEFAULT_ABS_TOL = 1e-4` — absolute floor (so values near zero
  don't trigger spurious relative-tolerance failures).
- `DEFAULT_VALID_FRACTION_ABS_TOL = 1e-3` — allow up to ~0.1% of
  pixels to drift in/out of the valid region (cloud-mask micro-shifts).

These are tight on purpose: the pipeline is deterministic on the same
inputs in principle, so a real numerical regression should exceed them
by orders of magnitude. Float-reduction-order drift between thread
counts is the only expected reason for sub-tolerance differences.

If a future change legitimately moves a stat past tolerance — for
example, a new aerosol model — refresh the goldens via
`regenerate_goldens.py` *after* a separate, principled validation
that the new outputs are correct. Don't refresh as a way to silence a
real regression.

## What's captured

For each output product:

- `shape`, `dtype` — exact-match
- `valid_fraction` — fraction of finite (non-masked) pixels
- `mean`, `std` — first two moments over the valid region
- `p01`, `p50`, `p99` — distribution percentiles (catches bimodal
  regressions that mean/std miss)
- `min`, `max` — clamps

Plus a subset of STAC properties (AOT/TCWV mean, view angles, EPSG)
and the WGS84 bbox.

The full rasters are **not** checked into the repo — only the summary
JSON. A 5-byte change in the pipeline will surface as a tolerance
breach long before any per-pixel diff would matter.

## Initial scene

The first regression scene captured here is **T33KWP / S2B / 2026-03-29**
over Namibia. It was produced by:

    PYTHONPATH=python pixi run -e rt6s python -m siac.cli process-s2 \
        tmp/real_cdse_mcd43_t33kwp/cache/s2/S2B_MSIL1C_..._T33KWP_....SAFE \
        --config tmp/real_cdse_mcd43_t33kwp_sixs.toml \
        --output-path /Users/fengyin/siac_runs/output_t33kwp_sixs

with the full chain: 6S RT backend, monthly-database surface prior,
MCD43 BRDF, CAMS atmospheric prior, omnicloudmask cloud detection.
Goldens at `goldens/t33kwp_sixs_20260329.json`.
