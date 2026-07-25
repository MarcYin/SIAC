"""Replay one saved SIAC surface likelihood with two atmospheric-prior widths.

The calibrated and native-product uncertainty solves share the exact same
surface prior, RT coefficients, masks, and cost cube.  The adaptive result uses
the native uncertainty only when the surface-only tile optimum is significantly
above the atmospheric prior.  AERONET truth is copied for downstream scoring
but never participates in either solve or the decision.
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
DEFAULT_BASELINE = ROOT / (
    "phaseD_results_lowcloud20_mediumphysics_baseline_costcube_mediumdev_20260713"
)
DEFAULT_CUBES = ROOT / (
    "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
)
DEFAULT_CALIBRATED_OUTPUT = ROOT / (
    "analysis/medium_aod_same_cube_calibrated_development_20260713"
)
DEFAULT_NATIVE_OUTPUT = ROOT / (
    "analysis/medium_aod_same_cube_product_unc_development_20260713"
)
DEFAULT_ADAPTIVE_OUTPUT = ROOT / (
    "analysis/medium_aod_same_cube_adaptive_z2576_development_20260713"
)


def _finite_mean(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.mean(finite)) if finite.size else None


def _within_ee(retrieved: float | None, truth: float) -> bool | None:
    if retrieved is None or not math.isfinite(retrieved):
        return None
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


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
    record["same_cube_prior_uncertainty"] = {
        **metadata,
        "aot_unc_mean": _finite_mean(solution[1]),
        "observation_cost_mean": _finite_mean(solution[2]),
    }
    return record


def replay_case(
    cube_path: Path,
    baseline_path: Path,
    calibrated_output_dir: Path,
    native_output_dir: Path,
    adaptive_output_dir: Path,
    *,
    z_threshold: float,
) -> str:
    baseline: dict[str, Any] = json.loads(baseline_path.read_text())
    with np.load(cube_path, allow_pickle=False) as data:
        cube = np.asarray(data["cube"], dtype=np.float64)
        axis = np.asarray(data["aot_axis"], dtype=np.float64)
        prior = np.asarray(data["aot_prior"], dtype=np.float64)
        calibrated_unc = np.asarray(data["aot_prior_unc"], dtype=np.float64)
        valid = np.asarray(data["solve_valid"], dtype=bool)
        pool_window = int(np.asarray(data["pool_window"]).reshape(-1)[0])
        min_count = int(np.asarray(data["min_count"]).reshape(-1)[0])

    native_unc = np.maximum(0.5 * prior, 0.02)
    calibrated_solution = _solve(
        cube, axis, prior, calibrated_unc, valid, pool_window, min_count
    )
    native_solution = _solve(cube, axis, prior, native_unc, valid, pool_window, min_count)
    surface_index = _surface_curve_min_index(cube)
    surface_min = float(axis[surface_index]) if surface_index is not None else None
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
        "shared_surface_prior_and_rt_cube": True,
        "surface_curve_min_aot": surface_min,
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
    decision = "native_product" if use_native else "calibrated"
    adaptive_record["same_cube_prior_uncertainty"] = {
        **adaptive_record["same_cube_prior_uncertainty"],
        "prior_uncertainty": "adaptive_conflict",
        "decision": decision,
    }

    for destination, record in (
        (calibrated_output_dir, calibrated_record),
        (native_output_dir, native_record),
        (adaptive_output_dir, adaptive_record),
    ):
        destination.mkdir(parents=True, exist_ok=True)
        (destination / baseline_path.name).write_text(json.dumps(record, indent=2) + "\n")
    return decision


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--cube-dir", type=Path, default=DEFAULT_CUBES)
    parser.add_argument(
        "--calibrated-output-dir", type=Path, default=DEFAULT_CALIBRATED_OUTPUT
    )
    parser.add_argument("--native-output-dir", type=Path, default=DEFAULT_NATIVE_OUTPUT)
    parser.add_argument("--adaptive-output-dir", type=Path, default=DEFAULT_ADAPTIVE_OUTPUT)
    parser.add_argument("--z-threshold", type=float, default=2.576)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.z_threshold <= 0.0:
        raise SystemExit("--z-threshold must be positive")
    completed = 0
    selected_native = 0
    skipped = 0
    for cube_path in sorted(args.cube_dir.glob("*.npz")):
        baseline_path = args.baseline_dir / f"{cube_path.stem}.json"
        if not baseline_path.exists():
            skipped += 1
            continue
        decision = replay_case(
            cube_path,
            baseline_path,
            args.calibrated_output_dir,
            args.native_output_dir,
            args.adaptive_output_dir,
            z_threshold=args.z_threshold,
        )
        selected_native += int(decision == "native_product")
        completed += 1
    print(
        f"replayed={completed} selected_native={selected_native} skipped={skipped} "
        f"adaptive_output={args.adaptive_output_dir}"
    )


if __name__ == "__main__":
    main()
