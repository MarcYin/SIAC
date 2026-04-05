#!/usr/bin/env python3
"""Benchmark and export quality-aware spatial AOD smoothing variants."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import rasterio
from scipy.ndimage import gaussian_filter

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOT = PROJECT_ROOT / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from siac.algorithms.solver.aod_smoothing import (  # noqa: E402
    SEED_EXCLUSION_QA_NAMES,
    build_trusted_aod_seed_mask,
    distance_to_seed_pixels,
    harmonic_surface,
    idw_knn_surface,
    median_from_seed_surface,
    multiscale_normalized_gaussian_surface,
    normalized_gaussian_surface,
    preserve_seed_values,
    sample_holdout_mask,
    score_holdout,
    whittaker_xy_surface,
)

LOGGER = logging.getLogger("aod_spatial_smoothing_experiment")

if TYPE_CHECKING:
    from collections.abc import Callable

QA_FILENAME_MAP = {
    "invalid_retrieval": "qa_invalid_retrieval.tif",
    "low_quality": "qa_low_quality.tif",
    "parameter_boundary": "qa_parameter_boundary.tif",
    "zero_obs_support": "qa_zero_obs_support.tif",
    "sharp_transition_excluded": "qa_sharp_transition_excluded.tif",
    "aot_lower_boundary": "qa_aot_lower_boundary.tif",
    "aot_upper_boundary": "qa_aot_upper_boundary.tif",
    "tcwv_lower_boundary": "qa_tcwv_lower_boundary.tif",
    "tcwv_upper_boundary": "qa_tcwv_upper_boundary.tif",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create quality-aware AOD smoothing and gap-filling experiments "
            "from an existing SIAC auxiliary directory."
        )
    )
    parser.add_argument(
        "--aot",
        type=Path,
        required=True,
        help="Path to the retrieved AOT GeoTIFF (for example auxiliary/aot.tif).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where smoothed rasters and the summary JSON will be written.",
    )
    parser.add_argument(
        "--border-pixels",
        type=int,
        default=1,
        help="Exclude this many outer image pixels from the trusted-seed mask.",
    )
    parser.add_argument(
        "--cloud-buffer-pixels",
        type=int,
        default=1,
        help="Dilate the cloud mask by this many pixels before selecting seeds.",
    )
    parser.add_argument(
        "--sharp-transition-buffer-pixels",
        type=int,
        default=1,
        help="Dilate sharp-transition exclusions by this many pixels before selecting seeds.",
    )
    parser.add_argument(
        "--holdout-fraction",
        type=float,
        default=0.1,
        help="Fraction of trusted seeds to hold out for validation.",
    )
    parser.add_argument(
        "--holdout-seed",
        type=int,
        default=42,
        help="Random seed used when sampling holdout pixels.",
    )
    parser.add_argument(
        "--surface-prior-dir",
        type=Path,
        default=None,
        help=(
            "Optional directory of surface-prior rasters used to reject bad spectral-mapping pixels. "
            "Defaults to <scene>/surface_prior when present."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="Logging verbosity.",
    )
    return parser.parse_args()


def _setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def _read_single_band(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    with rasterio.open(path) as src:
        data = src.read(1)
        profile = src.profile.copy()
    return data, profile


def _load_qa_masks(aux_dir: Path, *, shape: tuple[int, int]) -> dict[str, np.ndarray]:
    qa_masks: dict[str, np.ndarray] = {}
    for name, filename in QA_FILENAME_MAP.items():
        path = aux_dir / filename
        if not path.exists():
            continue
        mask, _profile = _read_single_band(path)
        array = np.asarray(mask, dtype=bool)
        if array.shape != shape:
            raise ValueError(f"QA mask {path} has shape {array.shape}, expected {shape}")
        qa_masks[name] = array
    return qa_masks


def _load_bad_spectral_mapping_mask(
    surface_prior_dir: Path | None,
    *,
    shape: tuple[int, int],
) -> np.ndarray | None:
    if surface_prior_dir is None or not surface_prior_dir.exists():
        return None

    band_paths = sorted(surface_prior_dir.glob("*.tif"))
    if not band_paths:
        return None

    valid_mask = np.ones(shape, dtype=bool)
    for path in band_paths:
        values, _profile = _read_single_band(path)
        band_valid = np.isfinite(np.asarray(values))
        if band_valid.shape != shape:
            raise ValueError(f"Surface-prior raster {path} has shape {band_valid.shape}, expected {shape}")
        valid_mask &= band_valid
    return ~valid_mask


def _write_float_raster(path: Path, data: np.ndarray, profile: dict[str, Any]) -> Path:
    updated = profile.copy()
    updated.update(dtype="float32", count=1, nodata=np.nan, compress="deflate")
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", **updated) as dst:
        dst.write(np.asarray(data, dtype=np.float32), 1)
    return path


def _write_mask_raster(path: Path, data: np.ndarray, profile: dict[str, Any]) -> Path:
    updated = profile.copy()
    updated.update(dtype="uint8", count=1, nodata=255, compress="lzw")
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", **updated) as dst:
        dst.write(np.asarray(data, dtype=np.uint8), 1)
    return path


def _summarize_mask(mask: np.ndarray) -> dict[str, float | int]:
    return {
        "pixel_count": int(np.count_nonzero(mask)),
        "fraction": float(np.mean(mask)),
    }


def _clip_surface(surface: np.ndarray, aot_min: float, aot_max: float) -> np.ndarray:
    return np.clip(np.asarray(surface, dtype=np.float32), a_min=aot_min, a_max=aot_max).astype(np.float32)


def main() -> int:
    args = _parse_args()
    _setup_logging(args.log_level)

    aot_path = args.aot.resolve()
    aux_dir = aot_path.parent
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else (aux_dir / "aot_smoothing").resolve()
    )

    aot, profile = _read_single_band(aot_path)
    aot = np.asarray(aot, dtype=np.float32)
    cloud_mask, _ = _read_single_band(aux_dir / "cloud_mask.tif")
    qa_masks = _load_qa_masks(aux_dir, shape=aot.shape)
    surface_prior_dir = (
        args.surface_prior_dir.resolve()
        if args.surface_prior_dir is not None
        else (aot_path.parents[1] / "surface_prior").resolve()
    )
    bad_spectral_mapping_mask = _load_bad_spectral_mapping_mask(
        surface_prior_dir,
        shape=aot.shape,
    )

    trusted_seed_mask = build_trusted_aod_seed_mask(
        aot,
        cloud_mask=np.asarray(cloud_mask, dtype=bool),
        qa_masks=qa_masks,
        spectral_mapping_mask=bad_spectral_mapping_mask,
        border_pixels=int(args.border_pixels),
        cloud_buffer_pixels=int(args.cloud_buffer_pixels),
        sharp_transition_buffer_pixels=int(args.sharp_transition_buffer_pixels),
    )
    holdout_mask = sample_holdout_mask(
        trusted_seed_mask,
        rng_seed=int(args.holdout_seed),
        holdout_fraction=float(args.holdout_fraction),
    )
    calibration_mask = trusted_seed_mask & ~holdout_mask

    finite_aot = aot[np.isfinite(aot)]
    if finite_aot.size == 0:
        raise ValueError(f"No finite AOT values found in {aot_path}")
    aot_min = float(np.min(finite_aot))
    aot_max = float(np.max(finite_aot))

    def _gaussian_surface(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        return normalized_gaussian_surface(values, mask, sigma=2.0)

    def _multiscale_gaussian_surface(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        return multiscale_normalized_gaussian_surface(
            values,
            mask,
            sigmas=(1.5, 3.0, 6.0),
            min_support_weight=1.0e-2,
        )

    def _idw_surface(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        return idw_knn_surface(values, mask, k=16, power=2.0)

    def _harmonic_full(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        return harmonic_surface(values, mask, iterations=250)

    def _harmonic_gapfill_preserve(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        return preserve_seed_values(harmonic_surface(values, mask, iterations=250), values, mask)

    def _harmonic_then_gaussian(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        harmonic = harmonic_surface(values, mask, iterations=250)
        return np.asarray(gaussian_filter(harmonic, sigma=1.0, mode="nearest"), dtype=np.float32)

    def _whittaker_surface(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        init = _multiscale_gaussian_surface(values, mask)
        return whittaker_xy_surface(
            values,
            mask,
            lambda_x=20.0,
            lambda_y=20.0,
            seed_weight=50.0,
            carry_weight=1.0,
            passes=3,
            init_surface=init,
        )

    def _whittaker_then_gaussian(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        whittaker = _whittaker_surface(values, mask)
        return np.asarray(gaussian_filter(whittaker, sigma=1.0, mode="nearest"), dtype=np.float32)

    def _median20_surface(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        init = _multiscale_gaussian_surface(values, mask)
        return median_from_seed_surface(values, mask, window_size=(20, 20), init_surface=init)

    methods: list[tuple[str, Callable[[np.ndarray, np.ndarray], np.ndarray], str]] = [
        ("gaussian_sigma2_seed_only", _gaussian_surface, "Normalized Gaussian interpolation from trusted seeds."),
        (
            "multiscale_gaussian_sigmas1p5_3_6_seed_only",
            _multiscale_gaussian_surface,
            "Multiscale normalized Gaussian interpolation using the smallest well-supported scale.",
        ),
        ("idw_k16_power2_seed_only", _idw_surface, "Inverse-distance weighted interpolation from trusted seeds."),
        (
            "harmonic_iter250_gapfill_preserve_seeds",
            _harmonic_gapfill_preserve,
            "Harmonic gap-fill that preserves the original trusted seed values exactly.",
        ),
        (
            "harmonic_iter250_then_gaussian_sigma1",
            _harmonic_then_gaussian,
            "Harmonic interpolation from trusted seeds followed by a light full-field Gaussian smooth.",
        ),
        (
            "whittaker_xy_lambda20_pass3_preserve_seeds",
            _whittaker_surface,
            "Weighted Whittaker smoothing along x then y with multiscale Gaussian initialization.",
        ),
        (
            "whittaker_xy_lambda20_pass3_then_gaussian_sigma1",
            _whittaker_then_gaussian,
            "Weighted Whittaker x/y smoothing followed by a light full-field Gaussian smooth.",
        ),
        (
            "median_20x20_after_multiscale_gaussian",
            _median20_surface,
            "20x20 spatial median filter applied to a trusted-seed surface initialized from multiscale Gaussian interpolation.",
        ),
    ]

    LOGGER.info(
        "Trusted seed pixels: %d (%.2f%% of scene); holdout=%d calibration=%d",
        int(np.count_nonzero(trusted_seed_mask)),
        100.0 * float(np.mean(trusted_seed_mask)),
        int(np.count_nonzero(holdout_mask)),
        int(np.count_nonzero(calibration_mask)),
    )

    distance = distance_to_seed_pixels(trusted_seed_mask)
    outputs: dict[str, str] = {}
    metrics: dict[str, dict[str, float | int | str]] = {}

    _write_mask_raster(output_dir / "trusted_seed_mask.tif", trusted_seed_mask, profile)
    _write_mask_raster(output_dir / "holdout_seed_mask.tif", holdout_mask, profile)
    _write_float_raster(output_dir / "distance_to_trusted_seed_px.tif", distance, profile)

    for name, method, description in methods:
        LOGGER.info("Running %s", name)
        validated_surface = _clip_surface(method(aot, calibration_mask), aot_min, aot_max)
        holdout_metrics = score_holdout(validated_surface, aot, holdout_mask)

        final_surface = _clip_surface(method(aot, trusted_seed_mask), aot_min, aot_max)
        output_path = output_dir / f"{name}.tif"
        _write_float_raster(output_path, final_surface, profile)
        outputs[name] = str(output_path)
        metrics[name] = {
            "description": description,
            **holdout_metrics.to_dict(),
            "scene_mean": float(np.nanmean(final_surface)),
            "scene_std": float(np.nanstd(final_surface)),
            "trusted_mean": float(np.nanmean(final_surface[trusted_seed_mask])),
            "gap_mean": float(np.nanmean(final_surface[~trusted_seed_mask])),
        }

    ranking = sorted(metrics.items(), key=lambda item: (item[1]["rmse"], item[1]["mae"]))
    summary = {
        "source_aot": str(aot_path),
        "auxiliary_dir": str(aux_dir),
        "output_dir": str(output_dir),
        "seed_selection": {
            "excluded_qa_layers": list(SEED_EXCLUSION_QA_NAMES),
            "uses_cloud_mask": True,
            "uses_sharp_transition_excluded": "sharp_transition_excluded" in qa_masks,
            "uses_bad_spectral_mapping_mask": bad_spectral_mapping_mask is not None,
            "surface_prior_dir": str(surface_prior_dir) if bad_spectral_mapping_mask is not None else None,
            "border_pixels": int(args.border_pixels),
            "cloud_buffer_pixels": int(args.cloud_buffer_pixels),
            "sharp_transition_buffer_pixels": int(args.sharp_transition_buffer_pixels),
        },
        "scene_summary": {
            "shape": [int(aot.shape[0]), int(aot.shape[1])],
            "finite_aot_pixels": int(np.count_nonzero(np.isfinite(aot))),
            "aot_min": aot_min,
            "aot_max": aot_max,
            "aot_mean": float(np.nanmean(aot)),
            "trusted_seed_mask": _summarize_mask(trusted_seed_mask),
            "holdout_seed_mask": _summarize_mask(holdout_mask),
            "calibration_seed_mask": _summarize_mask(calibration_mask),
            "distance_to_seed_px_mean": float(np.mean(distance)),
            "distance_to_seed_px_p95": float(np.quantile(distance, 0.95)),
            "bad_spectral_mapping": (
                _summarize_mask(bad_spectral_mapping_mask)
                if bad_spectral_mapping_mask is not None
                else None
            ),
        },
        "outputs": outputs,
        "metrics": metrics,
        "recommended_by_holdout_rmse": ranking[0][0] if ranking else None,
    }

    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    LOGGER.info("Wrote summary to %s", summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
