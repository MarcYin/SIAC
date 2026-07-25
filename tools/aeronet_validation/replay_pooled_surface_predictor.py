"""Replay a coefficient-free pooled historical-S2 surface predictor.

The operational experiment fits one 20-tree forest per historical realization
and aggregates their target-scene predictions.  This replay instead fits one
20-tree forest to all usable historical S2 realization pixels.  It uses the
same visible/anchor bands, localizer planes, target anchor correction, saved RT
residuals, masks, and atmospheric prior.  No AERONET value enters prediction.
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

from tools.aeronet_validation.replay_cost_cube_prior_uncertainty import (
    _build_record,
    _finite_mean,
)
from tools.aeronet_validation.replay_surface_oof_map import _solutions

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_BASELINE = ROOT / (
    "phaseD_results_lowcloud20_mediumphysics_baseline_costcube_mediumdev_20260713"
)
DEFAULT_CUBES = ROOT / (
    "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
)
DEFAULT_DUMPS = ROOT / "calib_dumps_c250"
DEFAULT_FIELDS = ROOT / (
    "analysis/medium_aod_pooled_et20_surface_fields_development_20260713"
)

ANCHOR_COLUMNS = (4, 5, 6)
VISIBLE_COLUMNS = (0, 1, 2, 3)
TARGET_COLUMNS = {"B01": 0, "B02": 1, "B03": 2, "B04": 3}
ANCHOR_BANDS = ("B8A", "B11", "B12")
MAD_TO_STD = 1.4826
VISIBLE_DEBIAS = {
    "B02": (-0.0003, 0.0243),
    "B03": (-0.0006, 0.0235),
    "B04": (-0.0011, 0.0223),
}


def _output_dir(uncertainty_mode: str, prior_mode: str) -> Path:
    return ROOT / (
        "analysis/medium_aod_pooled_et20_surface_"
        f"{uncertainty_mode}_{prior_mode}_development_20260713"
    )


def _correct_target_anchors(
    dump: np.lib.npyio.NpzFile,
    *,
    anchor_aot: float,
) -> np.ndarray:
    aot_grid = np.asarray(dump["aot_grid"], dtype=np.float64)
    corrected = []
    for band_name in ANCHOR_BANDS:
        toa = np.asarray(dump[f"toa_{band_name}"], dtype=np.float64)
        xap = float(np.interp(anchor_aot, aot_grid, dump[f"xap_{band_name}"]))
        xbp = float(np.interp(anchor_aot, aot_grid, dump[f"xbp_{band_name}"]))
        xcp = float(np.interp(anchor_aot, aot_grid, dump[f"xcp_{band_name}"]))
        y = xap * toa - xbp
        denominator = 1.0 + xcp * y
        corrected.append(
            np.where(
                np.isfinite(denominator) & (np.abs(denominator) > 1.0e-10),
                y / denominator,
                np.nan,
            )
        )
    return np.stack(corrected, axis=-1)


def fit_pooled_surface(
    composites: np.ndarray,
    target_anchors: np.ndarray,
    band_names: list[str],
    *,
    min_samples_leaf: int = 5,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Fit one ET20 model to all historical realization pixels."""
    comp = np.asarray(composites, dtype=np.float64)
    if comp.ndim != 4 or comp.shape[1] < 7:
        raise ValueError("composites must have shape (realization, >=7, y, x)")
    if target_anchors.shape != (*comp.shape[2:], len(ANCHOR_COLUMNS)):
        raise ValueError("target anchor grid does not match historical composites")

    target_columns = [TARGET_COLUMNS[name] for name in band_names]
    mean_visible = np.nanmean(comp[:, list(VISIBLE_COLUMNS)], axis=0)
    localizer = mean_visible.reshape(len(VISIBLE_COLUMNS), -1).T
    train_x_parts: list[np.ndarray] = []
    train_y_parts: list[np.ndarray] = []
    realization_counts: list[int] = []
    for index in range(comp.shape[0]):
        realization = comp[index].reshape(comp.shape[1], -1).T
        good = np.all(np.isfinite(realization), axis=1) & np.all(
            np.isfinite(localizer), axis=1
        )
        count = int(np.count_nonzero(good))
        if count < 200:
            continue
        train_x_parts.append(
            np.column_stack(
                [realization[good][:, list(ANCHOR_COLUMNS)], localizer[good]]
            )
        )
        train_y_parts.append(realization[good][:, target_columns])
        realization_counts.append(count)
    if len(train_x_parts) < 3:
        raise ValueError("at least three historical realizations are required")

    train_x = np.concatenate(train_x_parts, axis=0)
    train_y = np.concatenate(train_y_parts, axis=0)
    model = ExtraTreesRegressor(
        n_estimators=20,
        min_samples_leaf=int(min_samples_leaf),
        random_state=0,
        n_jobs=1,
    ).fit(train_x, train_y)

    flat_anchor = target_anchors.reshape(-1, target_anchors.shape[-1])
    target_good = np.all(
        np.isfinite(flat_anchor) & (flat_anchor > 0.0) & (flat_anchor < 1.2), axis=1
    ) & np.all(np.isfinite(localizer), axis=1)
    prediction = np.full((flat_anchor.shape[0], len(band_names)), np.nan)
    tree_stack: list[np.ndarray] = []
    if np.any(target_good):
        features = np.column_stack([flat_anchor[target_good], localizer[target_good]])
        prediction[target_good] = np.asarray(model.predict(features), dtype=np.float64)
        for estimator in model.estimators_:
            tree_prediction = np.asarray(estimator.predict(features), dtype=np.float64)
            if tree_prediction.ndim == 1:
                tree_prediction = tree_prediction[:, np.newaxis]
            tree_stack.append(tree_prediction)

    tree_sigma = np.full_like(prediction, np.nan)
    if tree_stack:
        stack = np.stack(tree_stack, axis=0)
        centre = np.median(stack, axis=0)
        tree_sigma[target_good] = (
            np.median(np.abs(stack - centre[np.newaxis]), axis=0) * MAD_TO_STD
        )
    height, width = comp.shape[2:]
    surface = prediction.T.reshape(len(band_names), height, width)
    sigma = tree_sigma.T.reshape(len(band_names), height, width)
    metadata = {
        "estimator": "pooled historical-S2 ExtraTreesRegressor",
        "n_estimators": 20,
        "min_samples_leaf": int(min_samples_leaf),
        "used_realizations": len(realization_counts),
        "realization_sample_counts": realization_counts,
        "training_samples": int(train_x.shape[0]),
        "target_valid_pixels": int(np.count_nonzero(target_good)),
        "tree_sigma_median": {
            name: _finite_mean(sigma[index])
            for index, name in enumerate(band_names)
        },
    }
    return surface, sigma, metadata


def replay_case(
    cube_path: Path,
    dump_path: Path,
    baseline_path: Path,
    *,
    z_threshold: float,
) -> None:
    baseline = json.loads(baseline_path.read_text())
    with np.load(cube_path, allow_pickle=False) as cube_data:
        signed_residual = np.asarray(
            cube_data["band_signed_residual_cube"], dtype=np.float64
        )
        old_surface = np.asarray(cube_data["boa_prior"], dtype=np.float64)
        base_uncertainty = np.asarray(cube_data["boa_unc"], dtype=np.float64)
        axis = np.asarray(cube_data["aot_axis"], dtype=np.float64)
        prior = np.asarray(cube_data["aot_prior"], dtype=np.float64)
        calibrated_unc = np.asarray(cube_data["aot_prior_unc"], dtype=np.float64)
        valid = np.asarray(cube_data["solve_valid"], dtype=bool)
        band_names = [
            str(value) for value in np.asarray(cube_data["band_names"]).tolist()
        ]
        pool_window = int(np.asarray(cube_data["pool_window"]).reshape(-1)[0])
        min_count = int(np.asarray(cube_data["min_count"]).reshape(-1)[0])
    anchor_iteration = baseline.get("anchor_iterate") or {}
    anchor_aot = anchor_iteration.get("pass1_scene_mean")
    if anchor_aot is None:
        anchor_aot = baseline.get("atmo_prior", {}).get("aot_median")
    anchor_aot = float(anchor_aot)

    field_path = DEFAULT_FIELDS / cube_path.name
    if field_path.exists():
        with np.load(field_path, allow_pickle=False) as fields:
            cached_names = [
                str(value) for value in np.asarray(fields["band_names"]).tolist()
            ]
            if cached_names != band_names:
                raise ValueError(
                    f"cached pooled bands {cached_names} do not match {band_names}"
                )
            pooled_surface = np.asarray(fields["surface"], dtype=np.float64)
            tree_sigma = np.asarray(fields["tree_sigma"], dtype=np.float64)
            usable = np.asarray(fields["usable"], dtype=bool)
        predictor_metadata = {
            "estimator": "cached pooled historical-S2 ExtraTreesRegressor",
            "n_estimators": 20,
            "tree_sigma_median": {
                name: _finite_mean(tree_sigma[index])
                for index, name in enumerate(band_names)
            },
        }
    else:
        with np.load(dump_path, allow_pickle=False) as dump:
            composites = np.asarray(dump["comp"], dtype=np.float64)
            target_anchors = _correct_target_anchors(dump, anchor_aot=anchor_aot)
        pooled_surface, tree_sigma, predictor_metadata = fit_pooled_surface(
            composites, target_anchors, band_names
        )
        usable = (
            np.isfinite(pooled_surface)
            & (pooled_surface > 0.001)
            & (pooled_surface < 0.6)
        )
    if pooled_surface.shape != old_surface.shape:
        raise ValueError(
            f"pooled surface shape {pooled_surface.shape} != cost surface {old_surface.shape}"
        )

    debias_offsets = np.asarray(
        [
            VISIBLE_DEBIAS.get(name, (0.0, 0.0))[0]
            + VISIBLE_DEBIAS.get(name, (0.0, 0.0))[1] * anchor_aot
            for name in band_names
        ],
        dtype=np.float64,
    )[:, np.newaxis, np.newaxis]
    raw_ensemble_surface = old_surface - debias_offsets
    raw_valid = (
        np.isfinite(raw_ensemble_surface)
        & (raw_ensemble_surface > 0.001)
        & (raw_ensemble_surface < 0.6)
    )
    new_surface = np.where(
        usable,
        pooled_surface,
        np.where(raw_valid, raw_ensemble_surface, old_surface),
    )
    blend_valid = usable & raw_valid
    blend_surface = np.where(
        blend_valid,
        0.5 * (new_surface + raw_ensemble_surface),
        np.where(usable, new_surface, raw_ensemble_surface),
    )
    current_blend_valid = usable & np.isfinite(old_surface)
    current_blend_surface = np.where(
        current_blend_valid,
        0.5 * (new_surface + old_surface),
        np.where(usable, new_surface, old_surface),
    )
    finite_tree_sigma = np.where(np.isfinite(tree_sigma), tree_sigma, base_uncertainty)
    pooled_uncertainty = np.sqrt(
        np.square(base_uncertainty) + np.square(finite_tree_sigma)
    )
    mixture_variance = 0.5 * (
        np.square(base_uncertainty) + np.square(raw_ensemble_surface - blend_surface)
    ) + 0.5 * (
        np.square(finite_tree_sigma) + np.square(new_surface - blend_surface)
    )
    mixture_uncertainty = np.maximum(np.sqrt(mixture_variance), 0.006)
    current_mixture_variance = 0.5 * (
        np.square(base_uncertainty)
        + np.square(old_surface - current_blend_surface)
    ) + 0.5 * (
        np.square(finite_tree_sigma)
        + np.square(new_surface - current_blend_surface)
    )
    current_mixture_uncertainty = np.maximum(
        np.sqrt(current_mixture_variance), 0.006
    )
    surface_recipes = {
        "baseunc": (new_surface, base_uncertainty),
        "treeunc": (new_surface, pooled_uncertainty),
        "blend50_mixunc": (blend_surface, mixture_uncertainty),
        "blend50_current_mixunc": (
            current_blend_surface,
            current_mixture_uncertainty,
        ),
    }

    DEFAULT_FIELDS.mkdir(parents=True, exist_ok=True)
    if not field_path.exists():
        np.savez_compressed(
            field_path,
            surface=new_surface.astype(np.float32),
            tree_sigma=tree_sigma.astype(np.float32),
            usable=usable,
            band_names=np.asarray(band_names),
        )
    for uncertainty_mode, (candidate_surface, uncertainty) in surface_recipes.items():
        candidate_signed_residual = (
            signed_residual
            + old_surface[:, np.newaxis]
            - candidate_surface[:, np.newaxis]
        )
        with np.errstate(invalid="ignore", divide="ignore"):
            cube = np.sum(
                np.square(candidate_signed_residual)
                / np.square(uncertainty[:, np.newaxis]),
                axis=0,
            )
        solutions, solve_metadata = _solutions(
            cube,
            axis,
            prior,
            calibrated_unc,
            valid,
            pool_window,
            min_count,
            z_threshold=z_threshold,
        )
        for prior_mode, solution in solutions.items():
            record = _build_record(
                baseline,
                solution=solution,
                metadata={
                    "uses_aeronet_in_retrieval": False,
                    "surface_source": "historical_s2_bestpixel_only",
                    "surface_predictor": "pooled_et20",
                    "removed_aeronet_surface_debias": True,
                    "surface_anchor_aot": anchor_aot,
                    "uncertainty_mode": uncertainty_mode,
                    "surface_blend": (
                        "equal coefficient-free model average"
                        if uncertainty_mode == "blend50_mixunc"
                        else (
                            "equal current-and-pooled model average"
                            if uncertainty_mode == "blend50_current_mixunc"
                            else "pooled_et20"
                        )
                    ),
                    "prior_mode": prior_mode,
                    "predictor": predictor_metadata,
                    "surface_replacement_fraction": float(np.mean(usable)),
                    "base_uncertainty_mean": {
                        name: _finite_mean(base_uncertainty[index])
                        for index, name in enumerate(band_names)
                    },
                    "effective_uncertainty_mean": {
                        name: _finite_mean(uncertainty[index])
                        for index, name in enumerate(band_names)
                    },
                    **solve_metadata,
                },
            )
            record["pooled_surface_replay"] = record.pop(
                "same_cube_prior_uncertainty"
            )
            destination = _output_dir(uncertainty_mode, prior_mode)
            destination.mkdir(parents=True, exist_ok=True)
            (destination / baseline_path.name).write_text(
                json.dumps(record, indent=2) + "\n"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--cube-dir", type=Path, default=DEFAULT_CUBES)
    parser.add_argument("--dump-dir", type=Path, default=DEFAULT_DUMPS)
    parser.add_argument("--z-threshold", type=float, default=2.576)
    parser.add_argument("--matchup-id", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    selected = set(args.matchup_id)
    completed = skipped = 0
    for cube_path in sorted(args.cube_dir.glob("*.npz")):
        if selected and cube_path.stem not in selected:
            continue
        baseline_path = args.baseline_dir / f"{cube_path.stem}.json"
        dump_path = args.dump_dir / f"{cube_path.stem}.npz"
        if not baseline_path.exists() or not dump_path.exists():
            skipped += 1
            continue
        replay_case(
            cube_path,
            dump_path,
            baseline_path,
            z_threshold=args.z_threshold,
        )
        completed += 1
    print(f"replayed={completed} skipped={skipped}")


if __name__ == "__main__":
    main()
