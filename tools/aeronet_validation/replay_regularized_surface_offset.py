"""Replay saved SIAC cost cubes with a constrained common surface offset.

The standard surface likelihood treats every visible-band reflectance error as
independent and absolute.  This replay adds one zero-mean, spectrally flat BOA
offset per pixel and AOT node, with a fixed Gaussian uncertainty.  It therefore
absorbs only small common-mode surface-normalization errors while preserving the
spectral residual that constrains aerosol.  AERONET truth is copied only for
downstream scoring and never enters the solve or adaptive prior decision.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from siac._rust_compat import surface_driven_pool_argmin

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_BASELINE = ROOT / "phaseD_results_lowcloud20_native_maiac_adaptive_b03_chi2_20260713"
DEFAULT_CUBES = ROOT / "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
DEFAULT_OUTPUT = ROOT / "analysis/medium_aod_regularized_offset_s006_development_20260713"
DEFAULT_PRODUCT_UNC_OUTPUT = ROOT / (
    "analysis/medium_aod_regularized_offset_s006_product_unc_development_20260713"
)
DEFAULT_ADAPTIVE_OUTPUT = ROOT / (
    "analysis/medium_aod_regularized_offset_s006_adaptive_z2576_development_20260713"
)


def _finite_mean(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.mean(finite)) if finite.size else None


def _finite_median(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.median(finite)) if finite.size else None


def _within_ee(retrieved: float | None, truth: float) -> bool | None:
    if retrieved is None or not math.isfinite(retrieved):
        return None
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


def regularized_offset_cost_cube(
    signed_residual: np.ndarray,
    band_uncertainty: np.ndarray,
    *,
    offset_sigma: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Profile a Gaussian-constrained common offset at every AOT/pixel node."""
    residual = np.asarray(signed_residual, dtype=np.float64)
    uncertainty = np.asarray(band_uncertainty, dtype=np.float64)
    if residual.ndim != 4:
        raise ValueError("signed residual must have shape (band, aot, y, x)")
    if uncertainty.shape != (residual.shape[0], residual.shape[2], residual.shape[3]):
        raise ValueError("band uncertainty must have shape (band, y, x)")
    if offset_sigma <= 0.0:
        raise ValueError("offset_sigma must be positive")

    inv_var = 1.0 / np.square(np.maximum(uncertainty, 1.0e-6))
    offset_precision = 1.0 / float(offset_sigma) ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        numerator = np.sum(residual * inv_var[:, np.newaxis], axis=0)
        denominator = np.sum(inv_var, axis=0)[np.newaxis] + offset_precision
        offset = numerator / denominator
        centered = residual - offset[np.newaxis]
        cost = np.sum(centered * centered * inv_var[:, np.newaxis], axis=0)
        cost += offset_precision * offset * offset
    return np.asarray(cost, dtype=np.float64), np.asarray(offset, dtype=np.float64)


def _surface_curve_min_index(cube: np.ndarray) -> int | None:
    flat = np.asarray(cube, dtype=np.float64).reshape(cube.shape[0], -1)
    node_cost = np.full(flat.shape[0], np.nan, dtype=np.float64)
    for index, values in enumerate(flat):
        finite = values[np.isfinite(values)]
        if finite.size:
            node_cost[index] = float(np.median(finite))
    valid = np.flatnonzero(np.isfinite(node_cost))
    if valid.size == 0:
        return None
    return int(valid[np.argmin(node_cost[valid])])


def _solve(
    cube: np.ndarray,
    axis: np.ndarray,
    prior: np.ndarray,
    prior_unc: np.ndarray,
    valid: np.ndarray,
    pool_window: int,
    min_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    result = surface_driven_pool_argmin(
        np.ascontiguousarray(cube, dtype=np.float64),
        np.ascontiguousarray(axis, dtype=np.float64),
        np.ascontiguousarray(prior, dtype=np.float64),
        np.ascontiguousarray(prior_unc, dtype=np.float64),
        np.ascontiguousarray(valid, dtype=bool),
        pool_window,
        min_count,
    )
    return tuple(np.asarray(value, dtype=np.float64) for value in result)  # type: ignore[return-value]


def _build_record(
    baseline: dict[str, Any],
    *,
    solution: tuple[np.ndarray, np.ndarray, np.ndarray],
    metadata: dict[str, Any],
) -> dict[str, Any]:
    record = dict(baseline)
    retrieved = _finite_mean(solution[0])
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
    record["regularized_surface_offset"] = {
        **metadata,
        "aot_unc_mean": _finite_mean(solution[1]),
        "observation_cost_mean": _finite_mean(solution[2]),
    }
    return record


def replay_case(
    cube_path: Path,
    baseline_path: Path,
    output_dir: Path,
    product_unc_output_dir: Path,
    adaptive_output_dir: Path,
    *,
    offset_sigma: float,
    z_threshold: float,
) -> None:
    baseline: dict[str, Any] = json.loads(baseline_path.read_text())
    with np.load(cube_path, allow_pickle=False) as data:
        if "band_signed_residual_cube" not in data.files:
            raise ValueError(
                f"{cube_path} predates signed residual dumps; regenerate the cost cube"
            )
        signed_residual = np.asarray(data["band_signed_residual_cube"], dtype=np.float64)
        band_uncertainty = np.asarray(data["boa_unc"], dtype=np.float64)
        axis = np.asarray(data["aot_axis"], dtype=np.float64)
        prior = np.asarray(data["aot_prior"], dtype=np.float64)
        calibrated_unc = np.asarray(data["aot_prior_unc"], dtype=np.float64)
        valid = np.asarray(data["solve_valid"], dtype=bool)
        band_names = [str(value) for value in np.asarray(data["band_names"]).tolist()]
        pool_window = int(np.asarray(data["pool_window"]).reshape(-1)[0])
        min_count = int(np.asarray(data["min_count"]).reshape(-1)[0])

    cube, offset = regularized_offset_cost_cube(
        signed_residual,
        band_uncertainty,
        offset_sigma=offset_sigma,
    )
    native_unc = np.maximum(0.5 * prior, 0.02)
    calibrated_solution = _solve(
        cube, axis, prior, calibrated_unc, valid, pool_window, min_count
    )
    native_solution = _solve(cube, axis, prior, native_unc, valid, pool_window, min_count)

    surface_index = _surface_curve_min_index(cube)
    surface_min = float(axis[surface_index]) if surface_index is not None else None
    offset_at_min = offset[surface_index] if surface_index is not None else np.asarray([])
    atmo_aot = _finite_mean(prior)
    calibrated_sigma = _finite_mean(calibrated_unc)
    standardized_conflict = None
    if (
        surface_min is not None
        and atmo_aot is not None
        and calibrated_sigma is not None
        and calibrated_sigma > 0.0
    ):
        standardized_conflict = (surface_min - atmo_aot) / calibrated_sigma
    use_native = bool(
        standardized_conflict is not None and standardized_conflict > z_threshold
    )

    common_metadata = {
        "uses_aeronet_in_retrieval": False,
        "model": "zero-mean Gaussian common BOA offset",
        "band_names": band_names,
        "offset_sigma": offset_sigma,
        "surface_curve_min_aot": surface_min,
        "offset_at_surface_min_mean": _finite_mean(offset_at_min),
        "offset_at_surface_min_median": _finite_median(offset_at_min),
        "atmospheric_prior_aot_mean": atmo_aot,
        "calibrated_prior_sigma_mean": calibrated_sigma,
        "standardized_positive_conflict": standardized_conflict,
        "z_threshold": z_threshold,
    }
    calibrated_record = _build_record(
        baseline,
        solution=calibrated_solution,
        metadata={**common_metadata, "prior_uncertainty": "calibrated"},
    )
    native_record = _build_record(
        baseline,
        solution=native_solution,
        metadata={**common_metadata, "prior_uncertainty": "native_product"},
    )
    adaptive_record = dict(native_record if use_native else calibrated_record)
    adaptive_record["regularized_surface_offset"] = {
        **adaptive_record["regularized_surface_offset"],
        "prior_uncertainty": "adaptive_conflict",
        "decision": "native_product" if use_native else "calibrated",
    }

    for destination, record in (
        (output_dir, calibrated_record),
        (product_unc_output_dir, native_record),
        (adaptive_output_dir, adaptive_record),
    ):
        destination.mkdir(parents=True, exist_ok=True)
        (destination / baseline_path.name).write_text(json.dumps(record, indent=2) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--cube-dir", type=Path, default=DEFAULT_CUBES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--product-unc-output-dir", type=Path, default=DEFAULT_PRODUCT_UNC_OUTPUT
    )
    parser.add_argument("--adaptive-output-dir", type=Path, default=DEFAULT_ADAPTIVE_OUTPUT)
    parser.add_argument("--offset-sigma", type=float, default=0.006)
    parser.add_argument("--z-threshold", type=float, default=2.576)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.offset_sigma <= 0.0:
        raise SystemExit("--offset-sigma must be positive")
    if args.z_threshold <= 0.0:
        raise SystemExit("--z-threshold must be positive")
    completed = 0
    skipped = 0
    for cube_path in sorted(args.cube_dir.glob("*.npz")):
        baseline_path = args.baseline_dir / f"{cube_path.stem}.json"
        if not baseline_path.exists():
            skipped += 1
            continue
        try:
            with np.load(cube_path, allow_pickle=False) as data:
                has_signed_residual = "band_signed_residual_cube" in data.files
        except FileNotFoundError:
            # Slurm producers replace an older dump before atomically writing
            # the regenerated signed cube. A concurrent replay may observe that
            # short gap; leave the case for the next pass.
            skipped += 1
            continue
        if not has_signed_residual:
            skipped += 1
            continue
        try:
            replay_case(
                cube_path,
                baseline_path,
                args.output_dir,
                args.product_unc_output_dir,
                args.adaptive_output_dir,
                offset_sigma=args.offset_sigma,
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
