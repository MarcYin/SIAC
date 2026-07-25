"""Check seed and tree-count robustness for one selected ExtraTrees AOD recipe."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from tools.aeronet_validation.aod_residual_calibration import (
    CalibrationSample,
    expected_error,
    feature_schema,
    load_samples,
    metrics,
    site_fold,
)
from tools.aeronet_validation.select_generic_aod_calibrator import (
    DEFAULT_MAX_SCENE_CLOUD_COVER,
    DEFAULT_ROOT,
    HOLDOUT_FOLD,
    TARGET_HITS,
    PredictionRecipe,
    _apply_recipe,
    _fit_predict,
    _ids,
    candidate_specs,
    filter_scene_cloud_cover,
)

SEEDS = (20260713, 20260731, 20260817, 20260903, 20261001)
TREE_COUNTS = (300, 600, 1200, 1500)
PRIMARY_VARIANT = "seed20260713_n1500"


def evaluate_slice(
    samples: list[CalibrationSample], prediction: np.ndarray, indices: np.ndarray
) -> dict[str, Any]:
    truth = np.asarray([sample.truth for sample in samples], dtype=float)[indices]
    baseline = np.asarray([sample.retrieved for sample in samples], dtype=float)[indices]
    candidate = prediction[indices]
    baseline_hit = np.abs(baseline - truth) <= expected_error(truth)
    candidate_hit = np.abs(candidate - truth) <= expected_error(truth)
    return {
        "count": len(indices),
        "baseline": metrics(truth, baseline),
        "candidate": metrics(truth, candidate),
        "gains": int((candidate_hit & ~baseline_hit).sum()),
        "losses": int((~candidate_hit & baseline_hit).sum()),
    }


def _selected_recipe(path: Path) -> PredictionRecipe:
    record = json.loads(path.read_text(encoding="utf-8"))["selected_recipe"]
    spec = next(spec for spec in candidate_specs() if spec.name == record["candidate"])
    if spec.family != "extra_trees_l1":
        raise ValueError(f"Expected ExtraTrees-L1 recipe, got {spec.family}")
    return PredictionRecipe(spec, float(record["shrinkage"]))


def _predict(
    recipe: PredictionRecipe,
    train: list[CalibrationSample],
    test: list[CalibrationSample],
    *,
    seed: int,
    tree_count: int,
) -> np.ndarray:
    names = feature_schema(train)
    baseline = np.asarray([sample.retrieved for sample in test], dtype=float)
    raw = _fit_predict(
        recipe.candidate,
        train,
        test,
        names,
        random_state=seed,
        tree_estimators=tree_count,
    )
    return _apply_recipe(baseline, raw, recipe.shrinkage)


def _variant_summary(
    target: list[CalibrationSample],
    target_prediction: np.ndarray,
    holdout: list[CalibrationSample],
    holdout_prediction: np.ndarray,
    external_sites: set[str],
) -> dict[str, Any]:
    folds = np.asarray([site_fold(sample.site) for sample in target], dtype=int)
    seen = np.asarray([sample.site in external_sites for sample in target], dtype=bool)
    return {
        "target_all": evaluate_slice(target, target_prediction, np.arange(len(target))),
        "target_development": evaluate_slice(target, target_prediction, np.flatnonzero(folds != 4)),
        "target_confirmation": evaluate_slice(target, target_prediction, np.flatnonzero(folds == 4)),
        "target_seen_sites": evaluate_slice(target, target_prediction, np.flatnonzero(seen)),
        "target_unseen_sites": evaluate_slice(target, target_prediction, np.flatnonzero(~seen)),
        "external_holdout": evaluate_slice(
            holdout, holdout_prediction, np.arange(len(holdout))
        ),
        "target_prediction": [float(value) for value in target_prediction],
        "external_holdout_prediction": [float(value) for value in holdout_prediction],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--selection-analysis", type=Path, required=True)
    parser.add_argument("--train-result-dir", type=Path, required=True)
    parser.add_argument("--train-context-dir", type=Path, required=True)
    parser.add_argument("--train-atmo-context-dir", type=Path, required=True)
    parser.add_argument("--train-mids", type=Path, required=True)
    parser.add_argument("--target-result-dir", type=Path, required=True)
    parser.add_argument("--target-context-dir", type=Path, required=True)
    parser.add_argument("--target-atmo-context-dir", type=Path, required=True)
    parser.add_argument("--target-mids", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    matchup_path = args.root / "matchups" / "matchups.csv"
    common = {"include_geography": False}
    train = filter_scene_cloud_cover(
        load_samples(
            args.train_result_dir,
            args.train_context_dir,
            matchup_path,
            _ids(args.train_mids),
            atmo_context_dir=args.train_atmo_context_dir,
            require_complete=True,
            **common,
        ),
        DEFAULT_MAX_SCENE_CLOUD_COVER,
    )
    target = load_samples(
        args.target_result_dir,
        args.target_context_dir,
        matchup_path,
        _ids(args.target_mids),
        atmo_context_dir=args.target_atmo_context_dir,
        **common,
    )
    recipe = _selected_recipe(args.selection_analysis)
    development = [sample for sample in train if site_fold(sample.site) != HOLDOUT_FOLD]
    holdout = [sample for sample in train if site_fold(sample.site) == HOLDOUT_FOLD]
    external_sites = {sample.site for sample in train}

    target_predictions: dict[str, np.ndarray] = {}
    holdout_predictions: dict[str, np.ndarray] = {}
    variants: dict[str, dict[str, Any]] = {}
    for seed in SEEDS:
        name = f"seed{seed}_n300"
        target_predictions[name] = _predict(
            recipe, train, target, seed=seed, tree_count=300
        )
        holdout_predictions[name] = _predict(
            recipe, development, holdout, seed=seed, tree_count=300
        )
    for tree_count in TREE_COUNTS[1:]:
        name = f"seed{SEEDS[0]}_n{tree_count}"
        target_predictions[name] = _predict(
            recipe, train, target, seed=SEEDS[0], tree_count=tree_count
        )
        holdout_predictions[name] = _predict(
            recipe, development, holdout, seed=SEEDS[0], tree_count=tree_count
        )

    seed_names = [f"seed{seed}_n300" for seed in SEEDS]
    target_stack = np.stack([target_predictions[name] for name in seed_names])
    holdout_stack = np.stack([holdout_predictions[name] for name in seed_names])
    target_predictions["seed_ensemble_mean_5x300"] = np.mean(target_stack, axis=0)
    target_predictions["seed_ensemble_median_5x300"] = np.median(target_stack, axis=0)
    holdout_predictions["seed_ensemble_mean_5x300"] = np.mean(holdout_stack, axis=0)
    holdout_predictions["seed_ensemble_median_5x300"] = np.median(holdout_stack, axis=0)

    for name, prediction in target_predictions.items():
        variants[name] = _variant_summary(
            target,
            prediction,
            holdout,
            holdout_predictions[name],
            external_sites,
        )
    analysis = {
        "schema_version": "siac-selected-aod-seed-robustness-v1",
        "recipe": recipe.name,
        "primary_stability_variant": PRIMARY_VARIANT,
        "primary_variant_predeclared_reason": (
            "Increase tree count to reduce finite-ensemble variance without changing features, "
            "loss, target transform, or decision rules."
        ),
        "target_hits_required": TARGET_HITS,
        "target_matchup_ids": [sample.matchup_id for sample in target],
        "external_holdout_matchup_ids": [sample.matchup_id for sample in holdout],
        "seed_300_target_hits": {
            name: variants[name]["target_all"]["candidate"]["hits"] for name in seed_names
        },
        "variants": variants,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "recipe": analysis["recipe"],
                "primary_stability_variant": PRIMARY_VARIANT,
                "seed_300_target_hits": analysis["seed_300_target_hits"],
                "primary": variants[PRIMARY_VARIANT],
            },
            indent=2,
        )
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
