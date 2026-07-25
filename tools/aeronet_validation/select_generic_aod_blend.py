"""Select a uniform two-model AOD blend using independent grouped data.

One champion recipe is selected for each predeclared estimator family using
site-grouped out-of-fold predictions on external development sites. Convex
blends among those champions are ranked on the same development predictions,
then checked once on reserved external sites before target transfer is measured.
Target AERONET labels never participate in component or weight selection.
"""

from __future__ import annotations

import argparse
import itertools
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

BLEND_WEIGHTS = tuple(float(value) for value in np.linspace(0.1, 0.9, 9))


def _recipe_prediction(
    recipe: PredictionRecipe,
    baseline: np.ndarray,
    raw: np.ndarray,
) -> np.ndarray:
    return _apply_recipe(baseline, raw, recipe.shrinkage)


def _metric_score(summary: dict[str, float | int]) -> tuple[int, float, float]:
    return (
        int(summary["hits"]),
        -float(summary["mae"]),
        -float(summary["rmse"]),
    )


def development_family_champions(
    samples: list[CalibrationSample],
    *,
    holdout_fold: int = HOLDOUT_FOLD,
) -> tuple[
    list[CalibrationSample],
    dict[str, PredictionRecipe],
    dict[str, np.ndarray],
    list[dict[str, Any]],
]:
    """Return the best external OOF recipe and predictions for each family."""
    development = [sample for sample in samples if site_fold(sample.site) != holdout_fold]
    names = feature_schema(development)
    baseline = np.asarray([sample.retrieved for sample in development], dtype=float)
    truth = np.asarray([sample.truth for sample in development], dtype=float)
    folds = np.asarray([site_fold(sample.site) for sample in development], dtype=int)
    all_rows: list[dict[str, Any]] = []
    all_predictions: dict[str, np.ndarray] = {}
    all_recipes: dict[str, PredictionRecipe] = {}

    for spec in candidate_specs():
        raw = np.full(len(development), np.nan, dtype=float)
        for fold in sorted(set(folds.tolist())):
            train = [
                sample
                for sample, value in zip(development, folds, strict=True)
                if value != fold
            ]
            indices = np.flatnonzero(folds == fold)
            test = [development[index] for index in indices]
            raw[indices] = _fit_predict(spec, train, test, names)
        if not np.all(np.isfinite(raw)):
            raise ValueError(f"Non-finite OOF predictions for {spec.name}")
        for shrinkage in (0.5, 0.75, 1.0, 1.25):
            recipe = PredictionRecipe(spec, shrinkage)
            prediction = _recipe_prediction(recipe, baseline, raw)
            summary = metrics(truth, prediction)
            all_rows.append(
                {
                    "recipe": recipe.name,
                    "family": spec.family,
                    "candidate": spec.name,
                    "shrinkage": shrinkage,
                    "metrics": summary,
                    "score": _metric_score(summary),
                }
            )
            all_predictions[recipe.name] = prediction
            all_recipes[recipe.name] = recipe

    all_rows.sort(key=lambda row: tuple(row["score"]), reverse=True)
    champion_rows: list[dict[str, Any]] = []
    seen_families: set[str] = set()
    for row in all_rows:
        family = str(row["family"])
        if family not in seen_families:
            champion_rows.append(row)
            seen_families.add(family)
    recipes = {row["recipe"]: all_recipes[row["recipe"]] for row in champion_rows}
    predictions = {row["recipe"]: all_predictions[row["recipe"]] for row in champion_rows}
    return development, recipes, predictions, champion_rows


def select_convex_blend(
    truth: np.ndarray,
    predictions: dict[str, np.ndarray],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Rank individual components and pairwise convex blends."""
    rows: list[dict[str, Any]] = []
    for name, prediction in predictions.items():
        summary = metrics(truth, prediction)
        rows.append(
            {
                "name": name,
                "left": name,
                "right": None,
                "left_weight": 1.0,
                "metrics": summary,
                "score": _metric_score(summary),
            }
        )
    for left, right in itertools.combinations(predictions, 2):
        for left_weight in BLEND_WEIGHTS:
            prediction = (
                left_weight * predictions[left]
                + (1.0 - left_weight) * predictions[right]
            )
            summary = metrics(truth, prediction)
            rows.append(
                {
                    "name": f"{left_weight:g}*{left}+{1.0-left_weight:g}*{right}",
                    "left": left,
                    "right": right,
                    "left_weight": left_weight,
                    "metrics": summary,
                    "score": _metric_score(summary),
                }
            )
    rows.sort(key=lambda row: tuple(row["score"]), reverse=True)
    return rows[0], rows


def evaluate_blend(
    selected: dict[str, Any],
    recipes: dict[str, PredictionRecipe],
    train: list[CalibrationSample],
    test: list[CalibrationSample],
) -> dict[str, Any]:
    """Fit selected components and return complete prediction receipts."""
    names = feature_schema(train)
    baseline = np.asarray([sample.retrieved for sample in test], dtype=float)
    truth = np.asarray([sample.truth for sample in test], dtype=float)
    required = [str(selected["left"])]
    if selected["right"] is not None:
        required.append(str(selected["right"]))
    raw_by_candidate: dict[str, np.ndarray] = {}
    prediction_by_recipe: dict[str, np.ndarray] = {}
    for name in required:
        recipe = recipes[name]
        candidate_name = recipe.candidate.name
        if candidate_name not in raw_by_candidate:
            raw_by_candidate[candidate_name] = _fit_predict(
                recipe.candidate, train, test, names
            )
        prediction_by_recipe[name] = _recipe_prediction(
            recipe,
            baseline,
            raw_by_candidate[candidate_name],
        )
    prediction = prediction_by_recipe[str(selected["left"])]
    if selected["right"] is not None:
        left_weight = float(selected["left_weight"])
        prediction = (
            left_weight * prediction
            + (1.0 - left_weight) * prediction_by_recipe[str(selected["right"])]
        )
    baseline_hit = np.abs(baseline - truth) <= expected_error(truth)
    candidate_hit = np.abs(prediction - truth) <= expected_error(truth)
    return {
        "count": len(test),
        "feature_count": len(names),
        "baseline": metrics(truth, baseline),
        "candidate": metrics(truth, prediction),
        "gains": int((candidate_hit & ~baseline_hit).sum()),
        "losses": int((~candidate_hit & baseline_hit).sum()),
        "predictions": [
            {
                "matchup_id": sample.matchup_id,
                "site": sample.site,
                "truth": sample.truth,
                "baseline": sample.retrieved,
                "candidate": float(prediction[index]),
                "baseline_within_ee": bool(baseline_hit[index]),
                "candidate_within_ee": bool(candidate_hit[index]),
            }
            for index, sample in enumerate(test)
        ],
    }


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
    development, recipes, oof_predictions, champions = development_family_champions(train)
    development_truth = np.asarray([sample.truth for sample in development], dtype=float)
    selected, ranking = select_convex_blend(development_truth, oof_predictions)
    holdout = [sample for sample in train if site_fold(sample.site) == HOLDOUT_FOLD]
    analysis = {
        "schema_version": "siac-generic-aod-convex-blend-v1",
        "selection_policy": {
            "target_labels_used_for_selection": False,
            "target_operational_covariates_used_for_selection": False,
            "geography_features": False,
            "component_policy": "best grouped-OOF recipe per predeclared estimator family",
            "blend_policy": "pairwise convex weights 0.1 through 0.9",
            "development_count": len(development),
            "holdout_count": len(holdout),
            "target_hits_required": TARGET_HITS,
        },
        "family_champions": champions,
        "selected_blend": selected,
        "blend_ranking": ranking,
        "external_holdout": evaluate_blend(selected, recipes, development, holdout),
        "target": evaluate_blend(selected, recipes, train, target),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "selected_blend": selected,
                "external_holdout": {
                    key: value
                    for key, value in analysis["external_holdout"].items()
                    if key != "predictions"
                },
                "target": {
                    key: value
                    for key, value in analysis["target"].items()
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
