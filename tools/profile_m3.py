#!/usr/bin/env python
"""Profile the M3 surface-prior query to find the 30-second bottleneck.

Monkey-patches key functions with timers to show where time goes.
Run from the repo root with the .venv activated.

.. warning::
    This script depends on ``siac.app.assembly.build_pipeline_runtime`` and
    ``siac.config.load_config`` — neither symbol exists in the current package
    (the assembly module was split into ``siac.app._assembly_*`` and
    ``load_config`` was renamed to ``siac.config.load.load_system_config``).
    The script is therefore marked broken pending a manual rewrite against the
    current API. (REVIEW.md §1.1 #1)
"""
from __future__ import annotations

raise ImportError(
    "tools/profile_m3.py is broken: it depends on siac.app.assembly."
    "build_pipeline_runtime and siac.config.load_config, both removed in the "
    "current refactor. Rewrite against siac.app._assembly_* and "
    "siac.config.load.load_system_config before using."
)

import functools  # noqa: E402,F401  - dead code retained for future rewrite
import logging  # noqa: E402,F401
import time  # noqa: E402,F401
from pathlib import Path  # noqa: E402,F401

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger("m3_profile")


def _timer(label: str):
    """Decorator that logs wall-clock time for a function call."""
    def decorator(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            t0 = time.perf_counter()
            result = fn(*args, **kwargs)
            dt = time.perf_counter() - t0
            logger.info("%s: %.3f s", label, dt)
            return result
        return wrapper
    return decorator


# ── Patch before import ──────────────────────────────────────────────────────
import siac.algorithms.surface.swir_refine as swir

_orig_query_obs_valid = swir._query_observation_valid_mask
swir._query_observation_valid_mask = _timer("_query_observation_valid_mask")(_orig_query_obs_valid)

_orig_resample_with_valid = swir._resample_dataset_with_validity
swir._resample_dataset_with_validity = _timer("_resample_dataset_with_validity")(_orig_resample_with_valid)

_orig_resample_geom = swir._resample_geometry_to_target_shape
swir._resample_geometry_to_target_shape = _timer("_resample_geometry_to_target_shape")(_orig_resample_geom)

_orig_resample_atmo = swir._resample_atmo_to_target_shape
swir._resample_atmo_to_target_shape = _timer("_resample_atmo_to_target_shape")(_orig_resample_atmo)

_orig_resample_cloud = swir._resample_cloud_mask_to_target_shape
swir._resample_cloud_mask_to_target_shape = _timer("_resample_cloud_mask_to_target_shape")(_orig_resample_cloud)

import siac.algorithms.correction.atmospheric as atmo_mod

_orig_correct = atmo_mod.AtmosphericCorrector.correct
atmo_mod.AtmosphericCorrector.correct = _timer("AtmosphericCorrector.correct")(_orig_correct)

import siac.algorithms.rt.lut.backend as lut_backend

_orig_compute_coeffs = lut_backend.ZarrLUTBackend.compute_coefficients
lut_backend.ZarrLUTBackend.compute_coefficients = _timer("compute_coefficients")(_orig_compute_coeffs)

_orig_scene_subset = lut_backend.ZarrLUTBackend._subset_spectral_lut_for_scene
lut_backend.ZarrLUTBackend._subset_spectral_lut_for_scene = _timer("_subset_spectral_lut_for_scene")(_orig_scene_subset)

_orig_band_grids = lut_backend.ZarrLUTBackend._get_or_build_spectral_band_grids
lut_backend.ZarrLUTBackend._get_or_build_spectral_band_grids = _timer("_get_or_build_spectral_band_grids")(_orig_band_grids)

_orig_interp_terms = lut_backend.ZarrLUTBackend._interpolate_spectral_band_terms
lut_backend.ZarrLUTBackend._interpolate_spectral_band_terms = _timer("_interpolate_spectral_band_terms")(_orig_interp_terms)

_orig_load_lut = lut_backend.ZarrLUTBackend._load_lut
lut_backend.ZarrLUTBackend._load_lut = _timer("_load_lut")(_orig_load_lut)

import siac.algorithms.surface.brdf_monthly_database as brdf_db

_orig_predict = brdf_db.MonthlyCompositeDatabase.predict_visible_with_diagnostics
brdf_db.MonthlyCompositeDatabase.predict_visible_with_diagnostics = _timer("predict_visible_with_diagnostics")(_orig_predict)

import siac.geo.resample as geo_resample

_orig_resample_coefficients = geo_resample.resample_coefficients_to_template
geo_resample.resample_coefficients_to_template = _timer("resample_coefficients_to_template")(_orig_resample_coefficients)

# ── Now import and run the pipeline ──────────────────────────────────────────
import toml

from siac.app.assembly import build_pipeline_runtime
from siac.config import load_config
from siac.workflows.pipeline import run_pipeline

CONFIG_PATH = Path("tmp/real_gcs_mcd43_cornbelt.toml")
INPUT_PATH = Path("tmp/real_gcs_mcd43_cornbelt/cache/s2")

def find_safe_dir(root: Path) -> Path:
    """Find the first .SAFE directory under the cache."""
    for p in root.rglob("*.SAFE"):
        if p.is_dir():
            return p
    raise FileNotFoundError(f"No .SAFE directory found under {root}")

def main():
    config = load_config(str(CONFIG_PATH))
    safe_dir = find_safe_dir(INPUT_PATH)
    logger.info("SAFE dir: %s", safe_dir)
    logger.info("Config: %s", CONFIG_PATH)

    runtime = build_pipeline_runtime(config, input_path=str(safe_dir))

    logger.info("=" * 60)
    logger.info("Starting pipeline run with profiling patches")
    logger.info("=" * 60)

    t_total = time.perf_counter()
    try:
        result = run_pipeline(
            input_path=str(safe_dir),
            config=config,
            **runtime,
        )
    except Exception as e:
        logger.error("Pipeline failed: %s", e, exc_info=True)
    logger.info("Total pipeline: %.3f s", time.perf_counter() - t_total)


if __name__ == "__main__":
    main()
