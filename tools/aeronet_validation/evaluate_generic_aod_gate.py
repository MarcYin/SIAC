"""Evaluate an external-only confidence gate for one generic AOD calibrator."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
from sklearn.ensemble import (
    ExtraTreesRegressor,
    GradientBoostingRegressor,
    HistGradientBoostingRegressor,
)
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from tools.aeronet_validation.aod_residual_calibration import (
    CalibrationSample,
    expected_error,
    feature_matrix,
    feature_schema,
    load_samples,
    metrics,
    site_fold,
)
from tools.aeronet_validation.select_generic_aod_calibrator import (
    HOLDOUT_FOLD,
    PredictionRecipe,
    _apply_recipe,
    _fit_predict,
    candidate_specs,
)

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
RANDOM_STATE = 20260713
GATE_THRESHOLDS = (-0.25, -0.10, -0.05, 0.0, 0.05, 0.10, 0.20)


@dataclass(frozen=True)
class GateSpec:
    name: str
    family: Literal["ridge", "gradient_l1", "hist_l1", "extra_trees_l1"]
    select_k: int


def gate_specs() -> tuple[GateSpec, ...]:
    return (
        GateSpec("ridge_gate25", "ridge", 25),
        GateSpec("gbr_l1_gate25", "gradient_l1", 25),
        GateSpec("gbr_l1_gate50", "gradient_l1", 50),
        GateSpec("hist_l1_gate25", "hist_l1", 25),
        GateSpec("hist_l1_gate50", "hist_l1", 50),
        GateSpec("et_l1_gate25", "extra_trees_l1", 25),
        GateSpec("et_l1_gate50", "extra_trees_l1", 50),
    )


def _gate_model(spec: GateSpec, n_features: int) -> Pipeline:
    if spec.family == "ridge":
        estimator: Any = Ridge(alpha=10.0)
    elif spec.family == "gradient_l1":
        estimator = GradientBoostingRegressor(
            n_estimators=150,
            learning_rate=0.03,
            max_depth=2,
            min_samples_leaf=10,
            loss="absolute_error",
            random_state=RANDOM_STATE,
        )
    elif spec.family == "hist_l1":
        estimator = HistGradientBoostingRegressor(
            max_iter=150,
            learning_rate=0.05,
            min_samples_leaf=20,
            l2_regularization=1.0,
            loss="absolute_error",
            random_state=RANDOM_STATE,
        )
    else:
        estimator = ExtraTreesRegressor(
            n_estimators=300,
            min_samples_leaf=10,
            max_features=0.75,
            criterion="absolute_error",
            n_jobs=4,
            random_state=RANDOM_STATE,
        )
    return Pipeline(
        (
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("select", SelectKBest(f_regression, k=min(spec.select_k, n_features))),
            ("model", estimator),
        )
    )


def gate_matrix(
    samples: list[CalibrationSample],
    feature_names: tuple[str, ...],
    candidate: np.ndarray,
) -> np.ndarray:
    """Build operational gate features; no AERONET field is accepted."""
    base = feature_matrix(samples, feature_names)
    baseline = np.asarray([sample.retrieved for sample in samples], dtype=float)
    correction = candidate - baseline
    cams = np.asarray(
        [
            sample.features.get("context_cams_total_aerosol_optical_depth_at_550nm_surface", np.nan)
            for sample in samples
        ],
        dtype=float,
    )
    maiac = np.asarray(
        [sample.features.get("atmo_aot_mean", np.nan) for sample in samples], dtype=float
    )
    meta = np.column_stack(
        (
            baseline,
            candidate,
            correction,
            np.abs(correction),
            cams - baseline,
            maiac - baseline,
            correction * (cams - baseline),
            correction * (maiac - baseline),
        )
    )
    return np.column_stack((base, meta))


def gate_target(samples: list[CalibrationSample], candidate: np.ndarray) -> np.ndarray:
    truth = np.asarray([sample.truth for sample in samples], dtype=float)
    baseline = np.asarray([sample.retrieved for sample in samples], dtype=float)
    scale = expected_error(truth)
    delta = np.abs(candidate - truth) / scale - np.abs(baseline - truth) / scale
    return np.clip(delta, -5.0, 5.0)


def cross_fitted_candidate(
    samples: list[CalibrationSample],
    recipe: PredictionRecipe,
    feature_names: tuple[str, ...],
) -> np.ndarray:
    folds = np.asarray([site_fold(sample.site) for sample in samples], dtype=int)
    prediction = np.full(len(samples), np.nan, dtype=float)
    for fold in sorted(set(folds.tolist())):
        train = [sample for sample, value in zip(samples, folds, strict=True) if value != fold]
        indices = np.flatnonzero(folds == fold)
        test = [samples[index] for index in indices]
        raw = _fit_predict(recipe.candidate, train, test, feature_names)
        baseline = np.asarray([sample.retrieved for sample in test], dtype=float)
        prediction[indices] = _apply_recipe(baseline, raw, recipe.shrinkage)
    if not np.all(np.isfinite(prediction)):
        raise ValueError("Cross-fitted candidate contains non-finite values")
    return prediction


def nested_gate_predictions(
    samples: list[CalibrationSample],
    recipe: PredictionRecipe,
    gate_spec: GateSpec,
    feature_names: tuple[str, ...],
) -> tuple[np.ndarray, np.ndarray]:
    """Return outer-fold candidate and gate predictions without label leakage."""
    folds = np.asarray([site_fold(sample.site) for sample in samples], dtype=int)
    candidate = np.full(len(samples), np.nan, dtype=float)
    gate_prediction = np.full(len(samples), np.nan, dtype=float)
    for fold in sorted(set(folds.tolist())):
        train_indices = np.flatnonzero(folds != fold)
        test_indices = np.flatnonzero(folds == fold)
        train = [samples[index] for index in train_indices]
        test = [samples[index] for index in test_indices]

        train_candidate = cross_fitted_candidate(train, recipe, feature_names)
        raw_test = _fit_predict(recipe.candidate, train, test, feature_names)
        test_baseline = np.asarray([sample.retrieved for sample in test], dtype=float)
        test_candidate = _apply_recipe(test_baseline, raw_test, recipe.shrinkage)

        model = _gate_model(gate_spec, len(feature_names) + 8)
        model.fit(
            gate_matrix(train, feature_names, train_candidate),
            gate_target(train, train_candidate),
        )
        candidate[test_indices] = test_candidate
        gate_prediction[test_indices] = model.predict(
            gate_matrix(test, feature_names, test_candidate)
        )
    if not np.all(np.isfinite(candidate)) or not np.all(np.isfinite(gate_prediction)):
        raise ValueError(f"Non-finite nested gate result for {gate_spec.name}")
    return candidate, gate_prediction


def _gated_metrics(
    samples: list[CalibrationSample],
    candidate: np.ndarray,
    gate_prediction: np.ndarray,
    threshold: float,
) -> dict[str, Any]:
    truth = np.asarray([sample.truth for sample in samples], dtype=float)
    baseline = np.asarray([sample.retrieved for sample in samples], dtype=float)
    use_candidate = gate_prediction < threshold
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
        "gate_prediction": gate_prediction,
        "use_candidate": use_candidate,
    }


def select_gate(
    samples: list[CalibrationSample], recipe: PredictionRecipe
) -> tuple[GateSpec, float, list[dict[str, Any]]]:
    names = feature_schema(samples)
    ranking: list[dict[str, Any]] = []
    for spec in gate_specs():
        candidate, gate_prediction = nested_gate_predictions(samples, recipe, spec, names)
        for threshold in GATE_THRESHOLDS:
            result = _gated_metrics(samples, candidate, gate_prediction, threshold)
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
    spec = next(spec for spec in gate_specs() if spec.name == winner["gate"])
    return spec, float(winner["threshold"]), ranking


def fit_apply_gate(
    train: list[CalibrationSample],
    test: list[CalibrationSample],
    recipe: PredictionRecipe,
    gate_spec: GateSpec,
    threshold: float,
) -> dict[str, Any]:
    names = feature_schema(train)
    train_candidate = cross_fitted_candidate(train, recipe, names)
    raw_test = _fit_predict(recipe.candidate, train, test, names)
    test_baseline = np.asarray([sample.retrieved for sample in test], dtype=float)
    test_candidate = _apply_recipe(test_baseline, raw_test, recipe.shrinkage)
    gate = _gate_model(gate_spec, len(names) + 8)
    gate.fit(gate_matrix(train, names, train_candidate), gate_target(train, train_candidate))
    gate_prediction = np.asarray(
        gate.predict(gate_matrix(test, names, test_candidate)), dtype=float
    )
    result = _gated_metrics(test, test_candidate, gate_prediction, threshold)
    truth = np.asarray([sample.truth for sample in test], dtype=float)
    baseline_hit = np.abs(test_baseline - truth) <= expected_error(truth)
    candidate_hit = np.abs(test_candidate - truth) <= expected_error(truth)
    return {
        "count": len(test),
        "baseline": metrics(truth, test_baseline),
        "ungated_candidate": metrics(truth, test_candidate),
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
                "gate_score": float(result["gate_prediction"][index]),
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
    gate_spec, threshold, ranking = select_gate(development, recipe)
    analysis = {
        "schema_version": "siac-generic-aod-confidence-gate-v1",
        "selection_policy": {
            "target_used_for_selection": False,
            "geography_features": False,
            "development_count": len(development),
            "holdout_count": len(holdout),
        },
        "base_recipe": recipe.name,
        "selected_gate": {"name": gate_spec.name, "threshold": threshold},
        "development_ranking": ranking,
        "external_holdout": fit_apply_gate(development, holdout, recipe, gate_spec, threshold),
        "target": fit_apply_gate(train, target, recipe, gate_spec, threshold),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    compact = {
        "base_recipe": analysis["base_recipe"],
        "selected_gate": analysis["selected_gate"],
        "external_holdout": {
            key: value
            for key, value in analysis["external_holdout"].items()
            if key != "predictions"
        },
        "target": {key: value for key, value in analysis["target"].items() if key != "predictions"},
    }
    print(json.dumps(compact, indent=2))
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
