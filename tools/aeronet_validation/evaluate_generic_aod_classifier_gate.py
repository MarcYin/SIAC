"""Evaluate a nested external-only classifier gate for AOD calibration."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
from sklearn.ensemble import (
    ExtraTreesClassifier,
    HistGradientBoostingClassifier,
    RandomForestClassifier,
)
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils.class_weight import compute_sample_weight
from tools.aeronet_validation.aod_residual_calibration import (
    CalibrationSample,
    expected_error,
    feature_schema,
    load_samples,
    metrics,
    site_fold,
)
from tools.aeronet_validation.evaluate_generic_aod_gate import (
    cross_fitted_candidate,
    gate_matrix,
)
from tools.aeronet_validation.select_generic_aod_calibrator import (
    DEFAULT_ROOT,
    HOLDOUT_FOLD,
    PredictionRecipe,
    _apply_recipe,
    _fit_predict,
    candidate_specs,
)

RANDOM_STATE = 20260713
GATE_THRESHOLDS = tuple(float(value) for value in np.linspace(0.1, 0.9, 9))


@dataclass(frozen=True)
class ClassifierGateSpec:
    name: str
    family: Literal["logistic", "hist", "extra_trees", "random_forest"]
    select_k: int
    min_samples_leaf: int = 10


def classifier_gate_specs() -> tuple[ClassifierGateSpec, ...]:
    return (
        ClassifierGateSpec("logistic_gate25", "logistic", 25),
        ClassifierGateSpec("logistic_gate50", "logistic", 50),
        ClassifierGateSpec("hist_gate25_leaf10", "hist", 25, 10),
        ClassifierGateSpec("hist_gate50_leaf20", "hist", 50, 20),
        ClassifierGateSpec("et_gate25_leaf5", "extra_trees", 25, 5),
        ClassifierGateSpec("et_gate50_leaf10", "extra_trees", 50, 10),
        ClassifierGateSpec("rf_gate25_leaf5", "random_forest", 25, 5),
        ClassifierGateSpec("rf_gate50_leaf10", "random_forest", 50, 10),
    )


def _classifier_model(spec: ClassifierGateSpec, n_features: int) -> Pipeline:
    if spec.family == "logistic":
        estimator: Any = LogisticRegression(
            C=0.1,
            class_weight="balanced",
            max_iter=1000,
            random_state=RANDOM_STATE,
        )
    elif spec.family == "hist":
        estimator = HistGradientBoostingClassifier(
            max_iter=150,
            learning_rate=0.05,
            min_samples_leaf=spec.min_samples_leaf,
            l2_regularization=1.0,
            class_weight="balanced",
            random_state=RANDOM_STATE,
        )
    elif spec.family == "extra_trees":
        estimator = ExtraTreesClassifier(
            n_estimators=300,
            min_samples_leaf=spec.min_samples_leaf,
            max_features=0.75,
            class_weight="balanced",
            n_jobs=4,
            random_state=RANDOM_STATE,
        )
    else:
        estimator = RandomForestClassifier(
            n_estimators=300,
            min_samples_leaf=spec.min_samples_leaf,
            max_features=0.75,
            class_weight="balanced",
            n_jobs=4,
            random_state=RANDOM_STATE,
        )
    return Pipeline(
        (
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("select", SelectKBest(f_classif, k=min(spec.select_k, n_features))),
            ("model", estimator),
        )
    )


def gate_preference_target(
    samples: list[CalibrationSample], candidate: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return preferred output and extra weight for hit-changing decisions."""
    truth = np.asarray([sample.truth for sample in samples], dtype=float)
    baseline = np.asarray([sample.retrieved for sample in samples], dtype=float)
    scale = expected_error(truth)
    baseline_error = np.abs(baseline - truth) / scale
    candidate_error = np.abs(candidate - truth) / scale
    baseline_hit = baseline_error <= 1.0
    candidate_hit = candidate_error <= 1.0
    preference = (candidate_error < baseline_error).astype(int)
    decisive = baseline_hit != candidate_hit
    sample_weight = compute_sample_weight("balanced", preference).astype(float)
    sample_weight *= np.where(decisive, 4.0, 1.0)
    return preference, sample_weight


def _fit_classifier(
    spec: ClassifierGateSpec,
    samples: list[CalibrationSample],
    feature_names: tuple[str, ...],
    candidate: np.ndarray,
) -> Pipeline:
    target, sample_weight = gate_preference_target(samples, candidate)
    if np.unique(target).size != 2:
        raise ValueError(f"Classifier gate training has one class for {spec.name}")
    model = _classifier_model(spec, len(feature_names) + 8)
    model.fit(
        gate_matrix(samples, feature_names, candidate),
        target,
        model__sample_weight=sample_weight,
    )
    return model


def nested_classifier_predictions(
    samples: list[CalibrationSample],
    recipe: PredictionRecipe,
    spec: ClassifierGateSpec,
    feature_names: tuple[str, ...],
) -> tuple[np.ndarray, np.ndarray]:
    """Return outer-fold candidate AOD and candidate-use probabilities."""
    folds = np.asarray([site_fold(sample.site) for sample in samples], dtype=int)
    candidate = np.full(len(samples), np.nan, dtype=float)
    probability = np.full(len(samples), np.nan, dtype=float)
    for fold in sorted(set(folds.tolist())):
        train_indices = np.flatnonzero(folds != fold)
        test_indices = np.flatnonzero(folds == fold)
        train = [samples[index] for index in train_indices]
        test = [samples[index] for index in test_indices]
        train_candidate = cross_fitted_candidate(train, recipe, feature_names)
        test_baseline = np.asarray([sample.retrieved for sample in test], dtype=float)
        raw_test = _fit_predict(recipe.candidate, train, test, feature_names)
        test_candidate = _apply_recipe(test_baseline, raw_test, recipe.shrinkage)
        model = _fit_classifier(spec, train, feature_names, train_candidate)
        candidate[test_indices] = test_candidate
        probability[test_indices] = model.predict_proba(
            gate_matrix(test, feature_names, test_candidate)
        )[:, 1]
    if not np.all(np.isfinite(candidate)) or not np.all(np.isfinite(probability)):
        raise ValueError(f"Non-finite nested classifier result for {spec.name}")
    return candidate, probability


def _gated_metrics(
    samples: list[CalibrationSample],
    candidate: np.ndarray,
    probability: np.ndarray,
    threshold: float,
) -> dict[str, Any]:
    truth = np.asarray([sample.truth for sample in samples], dtype=float)
    baseline = np.asarray([sample.retrieved for sample in samples], dtype=float)
    use_candidate = probability >= threshold
    prediction = np.where(use_candidate, candidate, baseline)
    baseline_hit = np.abs(baseline - truth) <= expected_error(truth)
    prediction_hit = np.abs(prediction - truth) <= expected_error(truth)
    return {
        "threshold": threshold,
        "candidate_used": int(use_candidate.sum()),
        "metrics": metrics(truth, prediction),
        "gains": int((prediction_hit & ~baseline_hit).sum()),
        "losses": int((~prediction_hit & baseline_hit).sum()),
        "prediction": prediction,
        "probability": probability,
        "use_candidate": use_candidate,
    }


def select_classifier_gate(
    samples: list[CalibrationSample], recipe: PredictionRecipe
) -> tuple[ClassifierGateSpec, float, list[dict[str, Any]]]:
    names = feature_schema(samples)
    ranking: list[dict[str, Any]] = []
    for spec in classifier_gate_specs():
        candidate, probability = nested_classifier_predictions(samples, recipe, spec, names)
        for threshold in GATE_THRESHOLDS:
            result = _gated_metrics(samples, candidate, probability, threshold)
            summary = result["metrics"]
            ranking.append(
                {
                    "gate": spec.name,
                    "threshold": threshold,
                    "candidate_used": result["candidate_used"],
                    "metrics": summary,
                    "gains": result["gains"],
                    "losses": result["losses"],
                    "score": (
                        int(summary["hits"]),
                        -int(result["losses"]),
                        -float(summary["mae"]),
                    ),
                }
            )
    ranking.sort(key=lambda item: tuple(item["score"]), reverse=True)
    winner = ranking[0]
    spec = next(spec for spec in classifier_gate_specs() if spec.name == winner["gate"])
    return spec, float(winner["threshold"]), ranking


def fit_apply_classifier_gate(
    train: list[CalibrationSample],
    test: list[CalibrationSample],
    recipe: PredictionRecipe,
    spec: ClassifierGateSpec,
    threshold: float,
) -> dict[str, Any]:
    names = feature_schema(train)
    train_candidate = cross_fitted_candidate(train, recipe, names)
    test_baseline = np.asarray([sample.retrieved for sample in test], dtype=float)
    test_truth = np.asarray([sample.truth for sample in test], dtype=float)
    raw_test = _fit_predict(recipe.candidate, train, test, names)
    test_candidate = _apply_recipe(test_baseline, raw_test, recipe.shrinkage)
    model = _fit_classifier(spec, train, names, train_candidate)
    probability = model.predict_proba(gate_matrix(test, names, test_candidate))[:, 1]
    result = _gated_metrics(test, test_candidate, probability, threshold)
    baseline_hit = np.abs(test_baseline - test_truth) <= expected_error(test_truth)
    candidate_hit = np.abs(test_candidate - test_truth) <= expected_error(test_truth)
    return {
        "count": len(test),
        "baseline": metrics(test_truth, test_baseline),
        "ungated_candidate": metrics(test_truth, test_candidate),
        "gated_candidate": result["metrics"],
        "candidate_used": result["candidate_used"],
        "gains": result["gains"],
        "losses": result["losses"],
        "ungated_gains": int((candidate_hit & ~baseline_hit).sum()),
        "ungated_losses": int((~candidate_hit & baseline_hit).sum()),
        "predictions": [
            {
                "matchup_id": sample.matchup_id,
                "site": sample.site,
                "truth": sample.truth,
                "baseline": sample.retrieved,
                "candidate": float(test_candidate[index]),
                "gated": float(result["prediction"][index]),
                "candidate_probability": float(probability[index]),
                "candidate_used": bool(result["use_candidate"][index]),
            }
            for index, sample in enumerate(test)
        ],
    }


def _recipe(path: Path) -> PredictionRecipe:
    record = json.loads(path.read_text(encoding="utf-8"))["selected_recipe"]
    candidate = next(spec for spec in candidate_specs() if spec.name == record["candidate"])
    return PredictionRecipe(candidate, float(record["shrinkage"]))


def _ids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--calibration-analysis", type=Path, required=True)
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
    target = load_samples(
        args.target_result_dir,
        args.target_context_dir,
        matchup_path,
        _ids(args.target_mids),
        atmo_context_dir=args.target_atmo_context_dir,
        **common,
    )
    recipe = _recipe(args.calibration_analysis)
    development = [sample for sample in train if site_fold(sample.site) != HOLDOUT_FOLD]
    holdout = [sample for sample in train if site_fold(sample.site) == HOLDOUT_FOLD]
    spec, threshold, ranking = select_classifier_gate(development, recipe)
    analysis = {
        "schema_version": "siac-generic-aod-classifier-gate-v1",
        "selection_policy": {
            "target_labels_used_for_selection": False,
            "target_operational_covariates_used_for_selection": False,
            "geography_features": False,
            "nested_site_grouping": True,
            "development_count": len(development),
            "holdout_count": len(holdout),
        },
        "base_recipe": recipe.name,
        "selected_gate": {"name": spec.name, "threshold": threshold},
        "development_ranking": ranking,
        "external_holdout": fit_apply_classifier_gate(
            development, holdout, recipe, spec, threshold
        ),
        "target": fit_apply_classifier_gate(train, target, recipe, spec, threshold),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "selected_gate": analysis["selected_gate"],
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
