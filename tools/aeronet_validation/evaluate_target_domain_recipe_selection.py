"""Evaluate target-domain recipe selection with an untouched site-fold confirmation.

All calibrators are fitted only on external matchups. Target labels from folds
0-3 may select one predeclared recipe; fold 4 remains untouched confirmation.
A nested five-fold replay estimates the complete target-domain selection process.
The selected recipe is uniform and consumes no AERONET value at inference time.
"""

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

TARGET_CONFIRMATION_FOLD = 4


def rank_prediction_matrix(
    truth: np.ndarray,
    predictions: dict[str, np.ndarray],
    indices: np.ndarray,
) -> list[dict[str, Any]]:
    """Rank fixed predictions on an explicitly supplied development subset."""
    ranking: list[dict[str, Any]] = []
    for name, prediction in predictions.items():
        summary = metrics(truth[indices], prediction[indices])
        ranking.append(
            {
                "recipe": name,
                "metrics": summary,
                "score": (
                    int(summary["hits"]),
                    -float(summary["mae"]),
                    -float(summary["rmse"]),
                ),
            }
        )
    ranking.sort(key=lambda row: tuple(row["score"]), reverse=True)
    return ranking


def _evaluation(
    samples: list[CalibrationSample],
    prediction: np.ndarray,
    indices: np.ndarray,
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
        "predictions": [
            {
                "matchup_id": samples[index].matchup_id,
                "site": samples[index].site,
                "site_fold": site_fold(samples[index].site),
                "truth": samples[index].truth,
                "baseline": samples[index].retrieved,
                "candidate": float(prediction[index]),
                "baseline_within_ee": bool(baseline_hit[position]),
                "candidate_within_ee": bool(candidate_hit[position]),
            }
            for position, index in enumerate(indices)
        ],
    }


def _external_evaluation(
    recipe: PredictionRecipe,
    train: list[CalibrationSample],
    test: list[CalibrationSample],
) -> dict[str, Any]:
    names = feature_schema(train)
    baseline = np.asarray([sample.retrieved for sample in test], dtype=float)
    raw = _fit_predict(recipe.candidate, train, test, names)
    prediction = _apply_recipe(baseline, raw, recipe.shrinkage)
    return _evaluation(test, prediction, np.arange(len(test)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
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
    train = load_samples(
        args.train_result_dir,
        args.train_context_dir,
        matchup_path,
        _ids(args.train_mids),
        atmo_context_dir=args.train_atmo_context_dir,
        require_complete=True,
        **common,
    )
    train = filter_scene_cloud_cover(train, DEFAULT_MAX_SCENE_CLOUD_COVER)
    target = load_samples(
        args.target_result_dir,
        args.target_context_dir,
        matchup_path,
        _ids(args.target_mids),
        atmo_context_dir=args.target_atmo_context_dir,
        **common,
    )
    names = feature_schema(train)
    baseline = np.asarray([sample.retrieved for sample in target], dtype=float)
    truth = np.asarray([sample.truth for sample in target], dtype=float)
    target_folds = np.asarray([site_fold(sample.site) for sample in target], dtype=int)
    recipes: dict[str, PredictionRecipe] = {}
    predictions: dict[str, np.ndarray] = {}
    for spec in candidate_specs():
        raw = _fit_predict(spec, train, target, names)
        for shrinkage in (0.5, 0.75, 1.0, 1.25):
            recipe = PredictionRecipe(spec, shrinkage)
            recipes[recipe.name] = recipe
            predictions[recipe.name] = _apply_recipe(baseline, raw, shrinkage)

    development_indices = np.flatnonzero(target_folds != TARGET_CONFIRMATION_FOLD)
    confirmation_indices = np.flatnonzero(target_folds == TARGET_CONFIRMATION_FOLD)
    development_ranking = rank_prediction_matrix(truth, predictions, development_indices)
    selected_name = str(development_ranking[0]["recipe"])
    selected_recipe = recipes[selected_name]
    selected_prediction = predictions[selected_name]

    nested_prediction = np.full(len(target), np.nan, dtype=float)
    nested_folds: list[dict[str, Any]] = []
    for fold in sorted(set(target_folds.tolist())):
        fold_development = np.flatnonzero(target_folds != fold)
        fold_test = np.flatnonzero(target_folds == fold)
        fold_ranking = rank_prediction_matrix(truth, predictions, fold_development)
        fold_recipe = str(fold_ranking[0]["recipe"])
        nested_prediction[fold_test] = predictions[fold_recipe][fold_test]
        nested_folds.append(
            {
                "held_out_fold": fold,
                "development_count": len(fold_development),
                "confirmation_count": len(fold_test),
                "selected_recipe": fold_recipe,
                "development_metrics": fold_ranking[0]["metrics"],
                "confirmation": _evaluation(target, predictions[fold_recipe], fold_test),
            }
        )
    if not np.all(np.isfinite(nested_prediction)):
        raise ValueError("Nested target-domain predictions are incomplete")

    external_development = [
        sample for sample in train if site_fold(sample.site) != HOLDOUT_FOLD
    ]
    external_holdout = [sample for sample in train if site_fold(sample.site) == HOLDOUT_FOLD]
    all_indices = np.arange(len(target))
    analysis = {
        "schema_version": "siac-target-domain-recipe-selection-v1",
        "selection_policy": {
            "external_labels_used_for_calibrator_fit": True,
            "target_labels_used_for_recipe_selection": True,
            "target_labels_used_at_inference": False,
            "geography_features": False,
            "target_development_folds": [0, 1, 2, 3],
            "target_confirmation_fold": TARGET_CONFIRMATION_FOLD,
            "target_confirmation_labels_used_for_selection": False,
            "target_hits_required": TARGET_HITS,
        },
        "selected_recipe": {
            "name": selected_recipe.name,
            "candidate": selected_recipe.candidate.name,
            "family": selected_recipe.candidate.family,
            "target": selected_recipe.candidate.target,
            "select_k": selected_recipe.candidate.select_k,
            "parameters": selected_recipe.candidate.parameter_dict(),
            "shrinkage": selected_recipe.shrinkage,
        },
        "target_development_ranking": development_ranking,
        "target_development": _evaluation(
            target, selected_prediction, development_indices
        ),
        "target_confirmation": _evaluation(
            target, selected_prediction, confirmation_indices
        ),
        "target_all_descriptive": _evaluation(target, selected_prediction, all_indices),
        "target_nested_selection": {
            "aggregate": _evaluation(target, nested_prediction, all_indices),
            "folds": nested_folds,
        },
        "external_holdout": _external_evaluation(
            selected_recipe, external_development, external_holdout
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "selected_recipe": analysis["selected_recipe"],
                "target_development": {
                    key: value
                    for key, value in analysis["target_development"].items()
                    if key != "predictions"
                },
                "target_confirmation": {
                    key: value
                    for key, value in analysis["target_confirmation"].items()
                    if key != "predictions"
                },
                "target_all_descriptive": {
                    key: value
                    for key, value in analysis["target_all_descriptive"].items()
                    if key != "predictions"
                },
                "target_nested_selection": analysis["target_nested_selection"]["aggregate"],
                "external_holdout": {
                    key: value
                    for key, value in analysis["external_holdout"].items()
                    if key != "predictions"
                },
            },
            indent=2,
        )
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
