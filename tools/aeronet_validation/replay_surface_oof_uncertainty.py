"""Replay SIAC with tile-derived temporal surface-prediction uncertainty.

For every historical S2 realization, this diagnostic predicts its visible
bands using ExtraTrees trained on the other realizations.  The robust temporal
prediction error is added in quadrature to the existing ensemble spread before
rebuilding the saved RT residual cost.  The procedure uses one S2 surface-prior
type and no AERONET values, site coefficients, or alternate surface source.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.aeronet_validation.replay_cost_cube_prior_uncertainty import (  # noqa: E402
    _build_record,
    _finite_mean,
    _solve,
    _surface_curve_min_index,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_BASELINE = ROOT / (
    "phaseD_results_lowcloud20_mediumphysics_baseline_costcube_mediumdev_20260713"
)
DEFAULT_CUBES = ROOT / (
    "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
)
DEFAULT_DUMPS = ROOT / "calib_dumps_c250"
DEFAULT_CALIBRATED_OUTPUT = ROOT / (
    "analysis/medium_aod_surface_oof_unc_calibrated_development_20260713"
)
DEFAULT_NATIVE_OUTPUT = ROOT / (
    "analysis/medium_aod_surface_oof_unc_product_development_20260713"
)
DEFAULT_ADAPTIVE_OUTPUT = ROOT / (
    "analysis/medium_aod_surface_oof_unc_adaptive_z2576_development_20260713"
)

ANCHOR_COLUMNS = (4, 5, 6)
VISIBLE_COLUMNS = (0, 1, 2, 3)
TARGET_COLUMNS = {"B01": 0, "B02": 1, "B03": 2, "B04": 3}
MAD_TO_STD = 1.4826


def estimate_temporal_oof_uncertainty(
    composites: np.ndarray,
    band_names: list[str],
    *,
    min_samples_leaf: int = 5,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Estimate one robust cross-realization prediction sigma per solve band."""
    comp = np.asarray(composites, dtype=np.float64)
    if comp.ndim != 4 or comp.shape[1] < 7:
        raise ValueError("composites must have shape (realization, >=7, y, x)")
    target_columns = [TARGET_COLUMNS[name] for name in band_names]
    mean_visible = np.nanmean(comp[:, list(VISIBLE_COLUMNS)], axis=0)
    localizer = mean_visible.reshape(len(VISIBLE_COLUMNS), -1).T

    trees: list[ExtraTreesRegressor] = []
    used_indices: list[int] = []
    flattened: dict[int, np.ndarray] = {}
    for index in range(comp.shape[0]):
        realization = comp[index].reshape(comp.shape[1], -1).T
        flattened[index] = realization
        good = np.all(np.isfinite(realization), axis=1) & np.all(
            np.isfinite(localizer), axis=1
        )
        if int(np.count_nonzero(good)) < 200:
            continue
        train_x = np.column_stack(
            [realization[good][:, list(ANCHOR_COLUMNS)], localizer[good]]
        )
        train_y = realization[good][:, target_columns]
        tree = ExtraTreesRegressor(
            n_estimators=20,
            random_state=0,
            min_samples_leaf=int(min_samples_leaf),
            n_jobs=1,
        ).fit(train_x, train_y)
        trees.append(tree)
        used_indices.append(index)
    if len(trees) < 3:
        raise ValueError("at least three usable realizations are required for temporal OOF")

    realization_errors: list[np.ndarray] = []
    held_sample_counts: list[int] = []
    for held_index in used_indices:
        held = flattened[held_index]
        good = np.all(np.isfinite(held), axis=1) & np.all(np.isfinite(localizer), axis=1)
        if int(np.count_nonzero(good)) < 200:
            continue
        features = np.column_stack([held[good][:, list(ANCHOR_COLUMNS)], localizer[good]])
        predictions = []
        for train_index, tree in zip(used_indices, trees, strict=True):
            if train_index == held_index:
                continue
            prediction = np.asarray(tree.predict(features), dtype=np.float64)
            if prediction.ndim == 1:
                prediction = prediction[:, np.newaxis]
            predictions.append(prediction)
        ensemble = np.median(np.stack(predictions, axis=0), axis=0)
        residual = ensemble - held[good][:, target_columns]
        realization_errors.append(np.median(np.abs(residual), axis=0) * MAD_TO_STD)
        held_sample_counts.append(int(np.count_nonzero(good)))

    error_stack = np.stack(realization_errors, axis=0)
    sigma = np.median(error_stack, axis=0)
    metadata = {
        "estimator": "leave-one-realization-out ET20 median absolute error",
        "used_realizations": len(used_indices),
        "held_sample_counts": held_sample_counts,
        "realization_sigma": error_stack.tolist(),
        "band_sigma": {name: float(value) for name, value in zip(band_names, sigma, strict=True)},
    }
    return np.asarray(sigma, dtype=np.float64), metadata


def estimate_temporal_oof_error_maps(
    composites: np.ndarray,
    band_names: list[str],
    *,
    min_samples_leaf: int = 5,
    min_support: int = 3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Estimate local signed bias and residual spread from historical S2 only.

    Each historical realization is predicted by the ExtraTrees ensemble built
    from the other realizations.  The median signed error estimates persistent
    local predictor bias; the MAD about that median estimates unresolved
    temporal error.  Sparse pixels fall back to the tile-wide robust band
    spread so this diagnostic never invents zero uncertainty from missing data.
    """
    comp = np.asarray(composites, dtype=np.float64)
    if comp.ndim != 4 or comp.shape[1] < 7:
        raise ValueError("composites must have shape (realization, >=7, y, x)")
    if min_support < 1:
        raise ValueError("min_support must be positive")

    target_columns = [TARGET_COLUMNS[name] for name in band_names]
    mean_visible = np.nanmean(comp[:, list(VISIBLE_COLUMNS)], axis=0)
    localizer = mean_visible.reshape(len(VISIBLE_COLUMNS), -1).T
    height, width = comp.shape[2:]

    trees: list[ExtraTreesRegressor] = []
    used_indices: list[int] = []
    flattened: dict[int, np.ndarray] = {}
    for index in range(comp.shape[0]):
        realization = comp[index].reshape(comp.shape[1], -1).T
        flattened[index] = realization
        good = np.all(np.isfinite(realization), axis=1) & np.all(
            np.isfinite(localizer), axis=1
        )
        if int(np.count_nonzero(good)) < 200:
            continue
        train_x = np.column_stack(
            [realization[good][:, list(ANCHOR_COLUMNS)], localizer[good]]
        )
        train_y = realization[good][:, target_columns]
        trees.append(
            ExtraTreesRegressor(
                n_estimators=20,
                random_state=0,
                min_samples_leaf=int(min_samples_leaf),
                n_jobs=1,
            ).fit(train_x, train_y)
        )
        used_indices.append(index)
    if len(trees) < 3:
        raise ValueError("at least three usable realizations are required for temporal OOF")

    residual_maps: list[np.ndarray] = []
    held_sample_counts: list[int] = []
    for held_index in used_indices:
        held = flattened[held_index]
        good = np.all(np.isfinite(held), axis=1) & np.all(
            np.isfinite(localizer), axis=1
        )
        if int(np.count_nonzero(good)) < 200:
            continue
        features = np.column_stack(
            [held[good][:, list(ANCHOR_COLUMNS)], localizer[good]]
        )
        predictions = []
        for train_index, tree in zip(used_indices, trees, strict=True):
            if train_index == held_index:
                continue
            prediction = np.asarray(tree.predict(features), dtype=np.float64)
            if prediction.ndim == 1:
                prediction = prediction[:, np.newaxis]
            predictions.append(prediction)
        ensemble = np.median(np.stack(predictions, axis=0), axis=0)
        residual = ensemble - held[good][:, target_columns]
        residual_map = np.full(
            (len(band_names), height * width), np.nan, dtype=np.float64
        )
        residual_map[:, good] = residual.T
        residual_maps.append(residual_map.reshape(len(band_names), height, width))
        held_sample_counts.append(int(np.count_nonzero(good)))

    error_stack = np.stack(residual_maps, axis=0)
    support = np.sum(np.isfinite(error_stack), axis=0)
    with np.errstate(invalid="ignore"):
        bias = np.nanmedian(error_stack, axis=0)
        absolute_centered = np.abs(error_stack - bias[np.newaxis])
        sigma = np.nanmedian(absolute_centered, axis=0) * MAD_TO_STD

    global_sigma = np.empty(len(band_names), dtype=np.float64)
    global_bias = np.empty(len(band_names), dtype=np.float64)
    for band_index in range(len(band_names)):
        values = error_stack[:, band_index]
        global_bias[band_index] = float(np.nanmedian(values))
        global_sigma[band_index] = float(
            np.nanmedian(np.abs(values - global_bias[band_index])) * MAD_TO_STD
        )
        sparse = (support[band_index] < min_support) | ~np.isfinite(bias[band_index])
        bias[band_index, sparse] = global_bias[band_index]
        invalid_sigma = (support[band_index] < min_support) | ~np.isfinite(
            sigma[band_index]
        )
        sigma[band_index, invalid_sigma] = global_sigma[band_index]

    metadata = {
        "estimator": "leave-one-realization-out ET20 local signed error and MAD",
        "used_realizations": len(used_indices),
        "held_sample_counts": held_sample_counts,
        "min_support": int(min_support),
        "band_bias_median": {
            name: float(np.nanmedian(bias[index]))
            for index, name in enumerate(band_names)
        },
        "band_sigma_median": {
            name: float(np.nanmedian(sigma[index]))
            for index, name in enumerate(band_names)
        },
        "band_sigma_p90": {
            name: float(np.nanpercentile(sigma[index], 90.0))
            for index, name in enumerate(band_names)
        },
        "band_support_median": {
            name: float(np.nanmedian(support[index]))
            for index, name in enumerate(band_names)
        },
    }
    return bias, sigma, support, metadata


def replay_case(
    cube_path: Path,
    dump_path: Path,
    baseline_path: Path,
    calibrated_output_dir: Path,
    native_output_dir: Path,
    adaptive_output_dir: Path,
    *,
    z_threshold: float,
) -> str:
    baseline: dict[str, Any] = json.loads(baseline_path.read_text())
    with np.load(cube_path, allow_pickle=False) as data:
        signed_residual = np.asarray(data["band_signed_residual_cube"], dtype=np.float64)
        base_uncertainty = np.asarray(data["boa_unc"], dtype=np.float64)
        axis = np.asarray(data["aot_axis"], dtype=np.float64)
        prior = np.asarray(data["aot_prior"], dtype=np.float64)
        calibrated_unc = np.asarray(data["aot_prior_unc"], dtype=np.float64)
        valid = np.asarray(data["solve_valid"], dtype=bool)
        band_names = [str(value) for value in np.asarray(data["band_names"]).tolist()]
        pool_window = int(np.asarray(data["pool_window"]).reshape(-1)[0])
        min_count = int(np.asarray(data["min_count"]).reshape(-1)[0])
    with np.load(dump_path, allow_pickle=False) as data:
        composites = np.asarray(data["comp"], dtype=np.float64)

    oof_sigma, oof_metadata = estimate_temporal_oof_uncertainty(composites, band_names)
    effective_uncertainty = np.sqrt(
        np.square(base_uncertainty) + np.square(oof_sigma[:, np.newaxis, np.newaxis])
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        cube = np.sum(
            np.square(signed_residual)
            / np.square(effective_uncertainty[:, np.newaxis]),
            axis=0,
        )

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
        "surface_source": "historical_s2_bestpixel_only",
        "shared_surface_prior_and_rt_cube": True,
        "oof": oof_metadata,
        "base_uncertainty_mean": {
            name: _finite_mean(base_uncertainty[index])
            for index, name in enumerate(band_names)
        },
        "effective_uncertainty_mean": {
            name: _finite_mean(effective_uncertainty[index])
            for index, name in enumerate(band_names)
        },
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
    calibrated_record["surface_oof_uncertainty"] = calibrated_record.pop(
        "same_cube_prior_uncertainty"
    )
    native_record["surface_oof_uncertainty"] = native_record.pop(
        "same_cube_prior_uncertainty"
    )
    adaptive_record = dict(native_record if use_native else calibrated_record)
    decision = "native_product" if use_native else "calibrated"
    adaptive_record["surface_oof_uncertainty"] = {
        **adaptive_record["surface_oof_uncertainty"],
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
    parser.add_argument("--dump-dir", type=Path, default=DEFAULT_DUMPS)
    parser.add_argument(
        "--calibrated-output-dir", type=Path, default=DEFAULT_CALIBRATED_OUTPUT
    )
    parser.add_argument("--native-output-dir", type=Path, default=DEFAULT_NATIVE_OUTPUT)
    parser.add_argument("--adaptive-output-dir", type=Path, default=DEFAULT_ADAPTIVE_OUTPUT)
    parser.add_argument("--z-threshold", type=float, default=2.576)
    parser.add_argument("--matchup-id", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.z_threshold <= 0.0:
        raise SystemExit("--z-threshold must be positive")
    selected_ids = set(args.matchup_id)
    completed = 0
    selected_native = 0
    skipped = 0
    for cube_path in sorted(args.cube_dir.glob("*.npz")):
        if selected_ids and cube_path.stem not in selected_ids:
            continue
        baseline_path = args.baseline_dir / f"{cube_path.stem}.json"
        dump_path = args.dump_dir / f"{cube_path.stem}.npz"
        if not baseline_path.exists() or not dump_path.exists():
            skipped += 1
            continue
        decision = replay_case(
            cube_path,
            dump_path,
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
