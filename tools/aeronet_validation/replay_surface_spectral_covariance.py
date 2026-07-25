"""Replay SIAC cost cubes with tile-estimated surface spectral covariance.

The operational surface likelihood uses one uncertainty per band and therefore
counts correlated B02/B03/B04 surface errors as independent evidence.  This
diagnostic estimates a regularized correlation matrix from temporal anomalies
in the same tile's historical S2 surface stack, then evaluates the saved RT
residuals with that covariance.  AERONET truth is copied only for downstream
scoring and is never used in the replay.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.covariance import LedoitWolf

from siac._rust_compat import surface_driven_pool_argmin

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_BASELINE = ROOT / "phaseD_results_lowcloud20_native_maiac_adaptive_b03_chi2_20260713"
DEFAULT_CUBES = ROOT / "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
DEFAULT_DUMPS = ROOT / "calib_dumps_c250"
DEFAULT_OUTPUT = ROOT / "analysis/medium_aod_surface_spectral_covariance_development_20260713"
DEFAULT_PRODUCT_UNC_OUTPUT = ROOT / (
    "analysis/medium_aod_surface_spectral_covariance_product_unc_development_20260713"
)
DEFAULT_ADAPTIVE_OUTPUT = ROOT / (
    "analysis/medium_aod_surface_spectral_covariance_adaptive_z2576_development_20260713"
)

COMPOSITE_COLUMNS = {"B01": 0, "B02": 1, "B03": 2, "B04": 3}
MAD_TO_STD = 1.4826


def _finite_mean(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.mean(finite)) if finite.size else None


def _within_ee(retrieved: float | None, truth: float) -> bool | None:
    if retrieved is None or not math.isfinite(retrieved):
        return None
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


def estimate_surface_correlation(
    composites: np.ndarray,
    band_names: list[str],
    *,
    spatial_stride: int,
) -> tuple[np.ndarray, int, float]:
    columns = [COMPOSITE_COLUMNS[name] for name in band_names]
    stack = np.asarray(composites, dtype=np.float64)[
        :, columns, ::spatial_stride, ::spatial_stride
    ]
    with np.errstate(invalid="ignore"):
        temporal_center = np.nanmedian(stack, axis=0)
    anomalies = np.moveaxis(stack - temporal_center[np.newaxis], 1, -1).reshape(
        -1, len(columns)
    )
    valid = np.all(np.isfinite(anomalies), axis=1)
    anomalies = anomalies[valid]
    if anomalies.shape[0] < 100:
        return np.eye(len(columns), dtype=np.float64), int(anomalies.shape[0]), 1.0

    center = np.median(anomalies, axis=0)
    spread = np.median(np.abs(anomalies - center), axis=0) * MAD_TO_STD
    positive = spread[np.isfinite(spread) & (spread > 1.0e-8)]
    fallback = float(np.median(positive)) if positive.size else 1.0
    spread = np.where(np.isfinite(spread) & (spread > 1.0e-8), spread, fallback)
    standardized = (anomalies - center) / spread
    robust = np.all(np.abs(standardized) <= 5.0, axis=1)
    standardized = standardized[robust]
    if standardized.shape[0] < 100:
        return np.eye(len(columns), dtype=np.float64), int(standardized.shape[0]), 1.0

    estimator = LedoitWolf(assume_centered=True).fit(standardized)
    covariance = np.asarray(estimator.covariance_, dtype=np.float64)
    diagonal = np.sqrt(np.maximum(np.diag(covariance), 1.0e-12))
    correlation = covariance / np.outer(diagonal, diagonal)
    correlation = 0.5 * (correlation + correlation.T)
    np.fill_diagonal(correlation, 1.0)
    return correlation, int(standardized.shape[0]), float(estimator.shrinkage_)


def covariance_cost_cube(
    band_cost: np.ndarray,
    band_residual: np.ndarray,
    correlation: np.ndarray,
) -> np.ndarray:
    costs = np.asarray(band_cost, dtype=np.float64)
    residuals = np.asarray(band_residual, dtype=np.float64)
    if costs.shape != residuals.shape or costs.ndim != 4:
        raise ValueError(
            "band cost and residual cubes must both have shape (band, aot, y, x)"
        )
    if correlation.shape != (costs.shape[0], costs.shape[0]):
        raise ValueError("surface correlation matrix does not match the solve bands")
    precision = np.linalg.pinv(correlation, rcond=1.0e-8, hermitian=True)
    normalized = np.sign(residuals) * np.sqrt(np.maximum(costs, 0.0))
    return np.einsum(
        "bayx,bc,cayx->ayx",
        normalized,
        precision,
        normalized,
        optimize=True,
    )


def _surface_curve_min_aot(cube: np.ndarray, axis: np.ndarray) -> float | None:
    flat = np.asarray(cube, dtype=np.float64).reshape(cube.shape[0], -1)
    node_cost = np.full(flat.shape[0], np.nan, dtype=np.float64)
    for index, values in enumerate(flat):
        finite = values[np.isfinite(values)]
        if finite.size:
            node_cost[index] = float(np.median(finite))
    valid = np.flatnonzero(np.isfinite(node_cost))
    if valid.size == 0:
        return None
    best = int(valid[np.argmin(node_cost[valid])])
    return float(axis[best])


def _solve(
    cube: np.ndarray,
    axis: np.ndarray,
    prior: np.ndarray,
    prior_unc: np.ndarray,
    valid: np.ndarray,
    pool_window: int,
    min_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    aod, aod_unc, observation_cost = surface_driven_pool_argmin(
        np.ascontiguousarray(cube, dtype=np.float64),
        np.ascontiguousarray(axis, dtype=np.float64),
        np.ascontiguousarray(prior, dtype=np.float64),
        np.ascontiguousarray(prior_unc, dtype=np.float64),
        np.ascontiguousarray(valid, dtype=bool),
        pool_window,
        min_count,
    )
    return (
        np.asarray(aod, dtype=np.float64),
        np.asarray(aod_unc, dtype=np.float64),
        np.asarray(observation_cost, dtype=np.float64),
    )


def _build_record(
    baseline: dict[str, Any],
    *,
    aod: np.ndarray,
    aod_unc: np.ndarray,
    observation_cost: np.ndarray,
    metadata: dict[str, Any],
) -> dict[str, Any]:
    record = dict(baseline)
    retrieved = _finite_mean(aod)
    truth = float(record["truth"])
    record["status"] = "OK" if retrieved is not None else "NO_RETRIEVAL"
    record["retrieved"] = retrieved
    record["scene_mean"] = retrieved
    record["err"] = retrieved - truth if retrieved is not None else None
    record["within_ee"] = _within_ee(retrieved, truth)
    record["retrieval_extraction"] = {
        "mode": "scene_mean",
        "selected": "scene_mean",
        "scene_mean": retrieved,
    }
    record["surface_spectral_covariance"] = {
        **metadata,
        "aot_unc_mean": _finite_mean(aod_unc),
        "observation_cost_mean": _finite_mean(observation_cost),
    }
    return record


def replay_case(
    cube_path: Path,
    dump_path: Path,
    baseline_path: Path,
    output_dir: Path,
    product_unc_output_dir: Path,
    adaptive_output_dir: Path,
    *,
    spatial_stride: int,
    z_threshold: float,
) -> None:
    baseline: dict[str, Any] = json.loads(baseline_path.read_text())
    with np.load(cube_path) as data:
        band_cost = np.asarray(data["band_cost_cube"], dtype=np.float64)
        if "band_signed_residual_cube" not in data.files:
            raise ValueError(
                f"{cube_path} predates signed residual dumps; regenerate the cost cube"
            )
        band_residual = np.asarray(data["band_signed_residual_cube"], dtype=np.float64)
        axis = np.asarray(data["aot_axis"], dtype=np.float64)
        prior = np.asarray(data["aot_prior"], dtype=np.float64)
        prior_unc = np.asarray(data["aot_prior_unc"], dtype=np.float64)
        valid = np.asarray(data["solve_valid"], dtype=bool)
        band_names = [str(value) for value in np.asarray(data["band_names"]).tolist()]
        pool_window = int(np.asarray(data["pool_window"]).reshape(-1)[0])
        min_count = int(np.asarray(data["min_count"]).reshape(-1)[0])
    with np.load(dump_path, allow_pickle=False) as data:
        composites = np.asarray(data["comp"], dtype=np.float64)

    correlation, sample_count, shrinkage = estimate_surface_correlation(
        composites,
        band_names,
        spatial_stride=spatial_stride,
    )
    cube = covariance_cost_cube(band_cost, band_residual, correlation)
    product_unc = np.maximum(0.5 * prior, 0.02)
    calibrated_solution = _solve(
        cube,
        axis,
        prior,
        prior_unc,
        valid,
        pool_window,
        min_count,
    )
    product_unc_solution = _solve(
        cube,
        axis,
        prior,
        product_unc,
        valid,
        pool_window,
        min_count,
    )
    surface_min = _surface_curve_min_aot(cube, axis)
    atmo_aot = _finite_mean(prior)
    calibrated_sigma = _finite_mean(prior_unc)
    standardized_conflict = None
    if (
        surface_min is not None
        and atmo_aot is not None
        and calibrated_sigma is not None
        and calibrated_sigma > 0.0
    ):
        standardized_conflict = (surface_min - atmo_aot) / calibrated_sigma
    use_product_unc = bool(
        standardized_conflict is not None and standardized_conflict > z_threshold
    )

    common_metadata = {
        "uses_aeronet_in_retrieval": False,
        "estimator": "LedoitWolf temporal-anomaly correlation",
        "band_names": band_names,
        "correlation": correlation.tolist(),
        "sample_count": sample_count,
        "shrinkage": shrinkage,
        "surface_curve_min_aot": surface_min,
        "atmospheric_prior_aot_mean": atmo_aot,
        "calibrated_prior_sigma_mean": calibrated_sigma,
        "standardized_positive_conflict": standardized_conflict,
        "z_threshold": z_threshold,
    }
    calibrated_record = _build_record(
        baseline,
        aod=calibrated_solution[0],
        aod_unc=calibrated_solution[1],
        observation_cost=calibrated_solution[2],
        metadata={**common_metadata, "prior_uncertainty": "calibrated"},
    )
    product_unc_record = _build_record(
        baseline,
        aod=product_unc_solution[0],
        aod_unc=product_unc_solution[1],
        observation_cost=product_unc_solution[2],
        metadata={**common_metadata, "prior_uncertainty": "native_product"},
    )
    adaptive_record = dict(product_unc_record if use_product_unc else calibrated_record)
    adaptive_record["surface_spectral_covariance"] = {
        **adaptive_record["surface_spectral_covariance"],
        "prior_uncertainty": "adaptive_conflict",
        "decision": "native_product" if use_product_unc else "calibrated",
    }

    for destination, record in (
        (output_dir, calibrated_record),
        (product_unc_output_dir, product_unc_record),
        (adaptive_output_dir, adaptive_record),
    ):
        destination.mkdir(parents=True, exist_ok=True)
        (destination / baseline_path.name).write_text(json.dumps(record, indent=2) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--cube-dir", type=Path, default=DEFAULT_CUBES)
    parser.add_argument("--dump-dir", type=Path, default=DEFAULT_DUMPS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--product-unc-output-dir", type=Path, default=DEFAULT_PRODUCT_UNC_OUTPUT
    )
    parser.add_argument("--adaptive-output-dir", type=Path, default=DEFAULT_ADAPTIVE_OUTPUT)
    parser.add_argument("--spatial-stride", type=int, default=3)
    parser.add_argument("--z-threshold", type=float, default=2.576)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.spatial_stride < 1:
        raise SystemExit("--spatial-stride must be positive")
    if args.z_threshold <= 0.0:
        raise SystemExit("--z-threshold must be positive")
    completed = 0
    skipped = 0
    for cube_path in sorted(args.cube_dir.glob("*.npz")):
        baseline_path = args.baseline_dir / f"{cube_path.stem}.json"
        dump_path = args.dump_dir / f"{cube_path.stem}.npz"
        if not baseline_path.exists() or not dump_path.exists():
            skipped += 1
            continue
        try:
            with np.load(cube_path, allow_pickle=False) as data:
                has_signed_residual = "band_signed_residual_cube" in data.files
        except FileNotFoundError:
            skipped += 1
            continue
        if not has_signed_residual:
            skipped += 1
            continue
        try:
            replay_case(
                cube_path,
                dump_path,
                baseline_path,
                args.output_dir,
                args.product_unc_output_dir,
                args.adaptive_output_dir,
                spatial_stride=args.spatial_stride,
                z_threshold=args.z_threshold,
            )
        except FileNotFoundError:
            skipped += 1
            continue
        completed += 1
    print(
        f"replayed={completed} skipped={skipped} output={args.output_dir} "
        f"product_unc_output={args.product_unc_output_dir} "
        f"adaptive_output={args.adaptive_output_dir}"
    )


if __name__ == "__main__":
    main()
