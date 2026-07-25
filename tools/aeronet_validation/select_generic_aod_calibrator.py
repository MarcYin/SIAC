"""Select one generic AOD calibrator using independent, site-held-out data.

The target benchmark is never used for model or recipe selection.  Candidate
models are ranked by site-grouped out-of-fold predictions on development sites,
then the winner is checked on a reserved set of sites before it is fitted to all
external samples and applied to the target cohort.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
from lightgbm import LGBMRegressor
from sklearn.ensemble import (
    ExtraTreesRegressor,
    GradientBoostingRegressor,
    HistGradientBoostingRegressor,
    RandomForestRegressor,
)
from sklearn.feature_selection import SelectKBest, f_classif, f_regression
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from tools.aeronet_validation.aod_residual_calibration import (
    CalibrationSample,
    expected_error,
    feature_matrix,
    feature_schema,
    load_samples,
    metrics,
    site_fold,
)

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
RANDOM_STATE = 20260713
HOLDOUT_FOLD = 4
TARGET_HITS = 133
DEFAULT_MAX_SCENE_CLOUD_COVER = 20.0
AOD_ERROR_OFFSET = 1.0 / 3.0
DOMAIN_AOD_BINS = np.asarray([0.0, 0.1, 0.2, 0.4, 0.6, 1.0, np.inf])


@dataclass(frozen=True)
class CandidateSpec:
    """A fixed estimator and prediction interpretation."""

    name: str
    family: Literal[
        "svr",
        "ridge",
        "extra_trees",
        "extra_trees_l1",
        "random_forest",
        "gradient",
        "gradient_l1",
        "hist",
        "hist_l1",
        "hist_quantile",
        "lightgbm",
        "lightgbm_l1",
    ]
    target: Literal["residual", "direct", "log_ratio"] = "residual"
    select_k: int | None = None
    params: tuple[tuple[str, float | int], ...] = ()

    def parameter_dict(self) -> dict[str, float | int]:
        return dict(self.params)


@dataclass(frozen=True)
class PredictionRecipe:
    """Uniform post-processing applied to every scene."""

    candidate: CandidateSpec
    shrinkage: float

    @property
    def name(self) -> str:
        return f"{self.candidate.name}__shrink{self.shrinkage:g}"


def candidate_specs() -> tuple[CandidateSpec, ...]:
    """Return a deliberately compact, predeclared candidate family."""
    specs: list[CandidateSpec] = [
        CandidateSpec("ridge25_a1", "ridge", select_k=25, params=(("alpha", 1.0),)),
        CandidateSpec("ridge50_a10", "ridge", select_k=50, params=(("alpha", 10.0),)),
    ]
    for select_k, c, gamma, epsilon in (
        (20, 0.6, 0.03, 0.002),
        (35, 0.6, 0.01, 0.01),
        (35, 0.6, 0.01, 0.002),
        (50, 1.0, 0.003, 0.01),
        (50, 2.0, 0.001, 0.02),
    ):
        specs.append(
            CandidateSpec(
                f"svr{select_k}_c{c:g}_g{gamma:g}_e{epsilon:g}",
                "svr",
                select_k=select_k,
                params=(("C", c), ("gamma", gamma), ("epsilon", epsilon)),
            )
        )
        if select_k in {35, 50}:
            specs.append(
                CandidateSpec(
                    f"svr_logratio{select_k}_c{c:g}_g{gamma:g}_e{epsilon:g}",
                    "svr",
                    target="log_ratio",
                    select_k=select_k,
                    params=(("C", c), ("gamma", gamma), ("epsilon", epsilon)),
                )
            )
    for leaf, max_features in ((2, 0.5), (4, 0.5), (6, 0.75), (10, 1.0)):
        label = str(max_features).replace(".", "p")
        for target in ("residual", "direct"):
            specs.append(
                CandidateSpec(
                    f"et_{target}_leaf{leaf}_mf{label}",
                    "extra_trees",
                    target=target,
                    select_k=50,
                    params=(("min_samples_leaf", leaf), ("max_features", max_features)),
                )
            )
        specs.append(
            CandidateSpec(
                f"et_logratio_leaf{leaf}_mf{label}",
                "extra_trees",
                target="log_ratio",
                select_k=50,
                params=(("min_samples_leaf", leaf), ("max_features", max_features)),
            )
        )
        if leaf in {4, 6}:
            for select_k in (20, 35, 50):
                specs.append(
                    CandidateSpec(
                        f"et_l1_logratio{select_k}_leaf{leaf}_mf{label}",
                        "extra_trees_l1",
                        target="log_ratio",
                        select_k=select_k,
                        params=(("min_samples_leaf", leaf), ("max_features", max_features)),
                    )
                )
    for leaf, max_features in ((2, 0.5), (5, 0.75), (10, 1.0)):
        label = str(max_features).replace(".", "p")
        for target in ("residual", "direct"):
            specs.append(
                CandidateSpec(
                    f"rf_{target}_leaf{leaf}_mf{label}",
                    "random_forest",
                    target=target,
                    select_k=50,
                    params=(("min_samples_leaf", leaf), ("max_features", max_features)),
                )
            )
        specs.append(
            CandidateSpec(
                f"rf_logratio_leaf{leaf}_mf{label}",
                "random_forest",
                target="log_ratio",
                select_k=50,
                params=(("min_samples_leaf", leaf), ("max_features", max_features)),
            )
        )
    for leaf, depth in ((3, 2), (5, 2), (5, 3), (10, 2)):
        specs.append(
            CandidateSpec(
                f"gbr_residual_leaf{leaf}_depth{depth}",
                "gradient",
                select_k=50,
                params=(("min_samples_leaf", leaf), ("max_depth", depth)),
            )
        )
        specs.append(
            CandidateSpec(
                f"gbr_logratio_leaf{leaf}_depth{depth}",
                "gradient",
                target="log_ratio",
                select_k=50,
                params=(("min_samples_leaf", leaf), ("max_depth", depth)),
            )
        )
        for select_k in (25, 50):
            specs.append(
                CandidateSpec(
                    f"gbr_l1_logratio{select_k}_leaf{leaf}_depth{depth}",
                    "gradient_l1",
                    target="log_ratio",
                    select_k=select_k,
                    params=(("min_samples_leaf", leaf), ("max_depth", depth)),
                )
            )
    for leaf, l2 in ((5, 0.1), (10, 0.1), (10, 1.0), (20, 1.0)):
        specs.append(
            CandidateSpec(
                f"hist_residual_leaf{leaf}_l2{l2:g}",
                "hist",
                select_k=50,
                params=(("min_samples_leaf", leaf), ("l2_regularization", l2)),
            )
        )
        specs.append(
            CandidateSpec(
                f"hist_logratio_leaf{leaf}_l2{l2:g}",
                "hist",
                target="log_ratio",
                select_k=50,
                params=(("min_samples_leaf", leaf), ("l2_regularization", l2)),
            )
        )
        specs.append(
            CandidateSpec(
                f"hist_l1_logratio_leaf{leaf}_l2{l2:g}",
                "hist_l1",
                target="log_ratio",
                select_k=50,
                params=(("min_samples_leaf", leaf), ("l2_regularization", l2)),
            )
        )
    for leaf in (10, 20):
        for quantile in (0.35, 0.4, 0.45, 0.55, 0.6):
            quantile_label = str(quantile).replace(".", "p")
            specs.append(
                CandidateSpec(
                    f"hist_q{quantile_label}_logratio_leaf{leaf}_l21",
                    "hist_quantile",
                    target="log_ratio",
                    select_k=50,
                    params=(
                        ("min_samples_leaf", leaf),
                        ("l2_regularization", 1.0),
                        ("quantile", quantile),
                    ),
                )
            )
    for child_samples, leaves, feature_fraction in (
        (10, 7, 0.6),
        (15, 7, 0.8),
        (20, 15, 0.6),
        (30, 15, 0.8),
    ):
        fraction_label = str(feature_fraction).replace(".", "p")
        params = (
            ("min_child_samples", child_samples),
            ("num_leaves", leaves),
            ("colsample_bytree", feature_fraction),
            ("reg_lambda", 1.0),
        )
        for select_k in (35, 50):
            specs.append(
                CandidateSpec(
                    f"lgbm_logratio{select_k}_child{child_samples}_leaves{leaves}_ff{fraction_label}",
                    "lightgbm",
                    target="log_ratio",
                    select_k=select_k,
                    params=params,
                )
            )
            specs.append(
                CandidateSpec(
                    f"lgbm_l1_logratio{select_k}_child{child_samples}_leaves{leaves}_ff{fraction_label}",
                    "lightgbm_l1",
                    target="log_ratio",
                    select_k=select_k,
                    params=params,
                )
            )
    return tuple(specs)


def _model(
    spec: CandidateSpec,
    n_features: int,
    *,
    random_state: int = RANDOM_STATE,
    tree_estimators: int = 300,
) -> Pipeline:
    params = spec.parameter_dict()
    steps: list[tuple[str, Any]] = [("impute", SimpleImputer(strategy="median"))]
    # Scaling before the univariate selector also avoids loss of precision in
    # f_regression for large count-valued diagnostics. Tree splits are invariant
    # to this monotonic transform.
    steps.append(("scale", StandardScaler()))
    if spec.select_k is not None:
        steps.append(("select", SelectKBest(f_regression, k=min(spec.select_k, n_features))))
    if spec.family == "svr":
        estimator = SVR(**params)
    elif spec.family == "ridge":
        estimator = Ridge(**params)
    elif spec.family in {"extra_trees", "extra_trees_l1"}:
        estimator = ExtraTreesRegressor(
            n_estimators=tree_estimators,
            n_jobs=4,
            random_state=random_state,
            criterion="absolute_error" if spec.family == "extra_trees_l1" else "squared_error",
            **params,
        )
    elif spec.family == "random_forest":
        estimator = RandomForestRegressor(
            n_estimators=tree_estimators,
            n_jobs=4,
            random_state=random_state,
            **params,
        )
    elif spec.family in {"gradient", "gradient_l1"}:
        estimator = GradientBoostingRegressor(
            n_estimators=200,
            learning_rate=0.03,
            loss="absolute_error" if spec.family == "gradient_l1" else "huber",
            random_state=random_state,
            **params,
        )
    elif spec.family in {"hist", "hist_l1"}:
        estimator = HistGradientBoostingRegressor(
            max_iter=200,
            learning_rate=0.05,
            loss="absolute_error" if spec.family == "hist_l1" else "squared_error",
            random_state=random_state,
            **params,
        )
    elif spec.family == "hist_quantile":
        quantile = float(params.pop("quantile"))
        estimator = HistGradientBoostingRegressor(
            max_iter=200,
            learning_rate=0.05,
            loss="quantile",
            quantile=quantile,
            random_state=random_state,
            **params,
        )
    elif spec.family in {"lightgbm", "lightgbm_l1"}:
        estimator = LGBMRegressor(
            n_estimators=300,
            learning_rate=0.03,
            objective="regression_l1" if spec.family == "lightgbm_l1" else "regression",
            random_state=random_state,
            n_jobs=4,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
            **params,
        )
    else:  # pragma: no cover - CandidateSpec constrains this.
        raise ValueError(spec.family)
    steps.append(("model", estimator))
    return Pipeline(steps)


def _fit_predict(
    spec: CandidateSpec,
    train: list[CalibrationSample],
    test: list[CalibrationSample],
    names: tuple[str, ...],
    *,
    fit_domain_samples: list[CalibrationSample] | None = None,
    fit_domain_weight_method: Literal["baseline_aod", "operational"] = "baseline_aod",
    random_state: int = RANDOM_STATE,
    tree_estimators: int = 300,
) -> np.ndarray:
    model = _model(
        spec,
        len(names),
        random_state=random_state,
        tree_estimators=tree_estimators,
    )
    train_x = feature_matrix(train, names)
    test_x = feature_matrix(test, names)
    baseline_train = np.asarray([sample.retrieved for sample in train], dtype=float)
    truth_train = np.asarray([sample.truth for sample in train], dtype=float)
    if spec.target == "residual":
        train_y = truth_train - baseline_train
    elif spec.target == "direct":
        train_y = truth_train
    else:
        train_y = np.log((truth_train + AOD_ERROR_OFFSET) / (baseline_train + AOD_ERROR_OFFSET))
    fit_params: dict[str, np.ndarray] = {}
    if fit_domain_samples is not None:
        fit_params["model__sample_weight"] = domain_weights(
            train,
            fit_domain_samples,
            method=fit_domain_weight_method,
        )
    model.fit(train_x, train_y, **fit_params)
    raw = np.asarray(model.predict(test_x), dtype=float)
    if spec.target == "residual":
        raw += np.asarray([sample.retrieved for sample in test], dtype=float)
    elif spec.target == "log_ratio":
        baseline_test = np.asarray([sample.retrieved for sample in test], dtype=float)
        raw = (baseline_test + AOD_ERROR_OFFSET) * np.exp(raw) - AOD_ERROR_OFFSET
    return np.clip(raw, 0.0, 4.0)


def _apply_recipe(
    baseline: np.ndarray,
    raw_prediction: np.ndarray,
    shrinkage: float,
) -> np.ndarray:
    return np.clip(baseline + shrinkage * (raw_prediction - baseline), 0.0, 4.0)


def _score(truth: np.ndarray, prediction: np.ndarray) -> tuple[int, float, float]:
    summary = metrics(truth, prediction)
    return int(summary["hits"]), -float(summary["mae"]), -float(summary["rmse"])


def baseline_domain_weights(
    samples: list[CalibrationSample], domain_samples: list[CalibrationSample]
) -> np.ndarray:
    """Estimate fixed-bin density ratios using retrieved AOD only."""
    source = np.asarray([sample.retrieved for sample in samples], dtype=float)
    domain = np.asarray([sample.retrieved for sample in domain_samples], dtype=float)
    source_count, _ = np.histogram(source, DOMAIN_AOD_BINS)
    domain_count, _ = np.histogram(domain, DOMAIN_AOD_BINS)
    source_fraction = (source_count + 1.0) / (source_count.sum() + source_count.size)
    domain_fraction = (domain_count + 1.0) / (domain_count.sum() + domain_count.size)
    ratio = np.clip(domain_fraction / source_fraction, 0.25, 4.0)
    index = np.clip(np.digitize(source, DOMAIN_AOD_BINS[1:-1]), 0, ratio.size - 1)
    weights = ratio[index]
    return weights / np.mean(weights)


def operational_domain_weights(
    samples: list[CalibrationSample], domain_samples: list[CalibrationSample]
) -> np.ndarray:
    """Estimate multivariate density ratios from operational features only."""
    combined = [*samples, *domain_samples]
    names = feature_schema(combined)
    matrix = feature_matrix(combined, names)
    domain_label = np.concatenate(
        (np.zeros(len(samples), dtype=int), np.ones(len(domain_samples), dtype=int))
    )
    model = Pipeline(
        (
            ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("scale", StandardScaler()),
            ("select", SelectKBest(f_classif, k=min(50, len(names)))),
            (
                "model",
                LogisticRegression(
                    C=0.1,
                    class_weight="balanced",
                    max_iter=1000,
                    random_state=RANDOM_STATE,
                ),
            ),
        )
    )
    model.fit(matrix, domain_label)
    probability = np.asarray(model.predict_proba(matrix[: len(samples)])[:, 1], dtype=float)
    odds = probability / np.maximum(1.0 - probability, 1e-6)
    weights = np.clip(odds, 0.1, 10.0)
    return weights / np.mean(weights)


def domain_weights(
    samples: list[CalibrationSample],
    domain_samples: list[CalibrationSample],
    *,
    method: Literal["baseline_aod", "operational"],
) -> np.ndarray:
    if method == "operational":
        return operational_domain_weights(samples, domain_samples)
    return baseline_domain_weights(samples, domain_samples)


def _weighted_metrics(
    truth: np.ndarray, prediction: np.ndarray, weights: np.ndarray
) -> dict[str, float]:
    error = prediction - truth
    total = float(np.sum(weights))
    within = np.abs(error) <= expected_error(truth)
    return {
        "within_ee_percent": float(100.0 * np.sum(weights * within) / total),
        "mae": float(np.sum(weights * np.abs(error)) / total),
        "rmse": float(np.sqrt(np.sum(weights * error**2) / total)),
        "bias": float(np.sum(weights * error) / total),
    }


def select_recipe(
    samples: list[CalibrationSample],
    *,
    holdout_fold: int = HOLDOUT_FOLD,
    score_domain_samples: list[CalibrationSample] | None = None,
    fit_domain_samples: list[CalibrationSample] | None = None,
    domain_weight_method: Literal["baseline_aod", "operational"] = "baseline_aod",
) -> tuple[PredictionRecipe, list[dict[str, Any]]]:
    """Select a recipe using OOF predictions from non-holdout sites only."""
    development = [sample for sample in samples if site_fold(sample.site) != holdout_fold]
    names = feature_schema(development)
    baseline = np.asarray([sample.retrieved for sample in development], dtype=float)
    truth = np.asarray([sample.truth for sample in development], dtype=float)
    weights = (
        domain_weights(
            development,
            score_domain_samples,
            method=domain_weight_method,
        )
        if score_domain_samples is not None
        else np.ones(len(development), dtype=float)
    )
    folds = np.asarray([site_fold(sample.site) for sample in development], dtype=int)
    ranking: list[dict[str, Any]] = []
    for spec in candidate_specs():
        raw = np.full(len(development), np.nan, dtype=float)
        for fold in sorted(set(folds.tolist())):
            train = [
                sample for sample, value in zip(development, folds, strict=True) if value != fold
            ]
            indices = np.flatnonzero(folds == fold)
            test = [development[index] for index in indices]
            raw[indices] = _fit_predict(
                spec,
                train,
                test,
                names,
                fit_domain_samples=fit_domain_samples,
                fit_domain_weight_method=domain_weight_method,
            )
        if not np.all(np.isfinite(raw)):
            raise ValueError(f"Non-finite OOF predictions for {spec.name}")
        for shrinkage in (0.5, 0.75, 1.0, 1.25):
            prediction = _apply_recipe(baseline, raw, shrinkage)
            weighted = _weighted_metrics(truth, prediction, weights)
            ranking.append(
                {
                    "recipe": f"{spec.name}__shrink{shrinkage:g}",
                    "candidate": spec.name,
                    "shrinkage": shrinkage,
                    "metrics": metrics(truth, prediction),
                    "domain_weighted_metrics": weighted,
                    "score": (
                        weighted["within_ee_percent"],
                        -weighted["mae"],
                        -weighted["rmse"],
                    ),
                }
            )
    ranking.sort(key=lambda item: tuple(item["score"]), reverse=True)
    winner = ranking[0]
    spec = next(spec for spec in candidate_specs() if spec.name == winner["candidate"])
    return PredictionRecipe(spec, float(winner["shrinkage"])), ranking


def evaluate_recipe(
    recipe: PredictionRecipe,
    train: list[CalibrationSample],
    test: list[CalibrationSample],
    *,
    fit_domain_samples: list[CalibrationSample] | None = None,
    metric_domain_samples: list[CalibrationSample] | None = None,
    domain_weight_method: Literal["baseline_aod", "operational"] = "baseline_aod",
) -> dict[str, Any]:
    names = feature_schema(train)
    baseline = np.asarray([sample.retrieved for sample in test], dtype=float)
    truth = np.asarray([sample.truth for sample in test], dtype=float)
    raw = _fit_predict(
        recipe.candidate,
        train,
        test,
        names,
        fit_domain_samples=fit_domain_samples,
        fit_domain_weight_method=domain_weight_method,
    )
    prediction = _apply_recipe(baseline, raw, recipe.shrinkage)
    baseline_hit = np.abs(baseline - truth) <= expected_error(truth)
    candidate_hit = np.abs(prediction - truth) <= expected_error(truth)
    result = {
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
    if metric_domain_samples is not None:
        weights = domain_weights(test, metric_domain_samples, method=domain_weight_method)
        result["domain_weighted_baseline"] = _weighted_metrics(truth, baseline, weights)
        result["domain_weighted_candidate"] = _weighted_metrics(truth, prediction, weights)
    return result


def diagnostic_target_ranking(
    train: list[CalibrationSample],
    target: list[CalibrationSample],
    *,
    fit_domain_samples: list[CalibrationSample] | None = None,
    domain_weight_method: Literal["baseline_aod", "operational"] = "baseline_aod",
) -> list[dict[str, Any]]:
    """Measure transfer for every candidate without using it for selection."""
    names = feature_schema(train)
    baseline = np.asarray([sample.retrieved for sample in target], dtype=float)
    truth = np.asarray([sample.truth for sample in target], dtype=float)
    ranking: list[dict[str, Any]] = []
    for spec in candidate_specs():
        raw = _fit_predict(
            spec,
            train,
            target,
            names,
            fit_domain_samples=fit_domain_samples,
            fit_domain_weight_method=domain_weight_method,
        )
        for shrinkage in (0.5, 0.75, 1.0, 1.25):
            prediction = _apply_recipe(baseline, raw, shrinkage)
            ranking.append(
                {
                    "recipe": f"{spec.name}__shrink{shrinkage:g}",
                    "candidate": spec.name,
                    "shrinkage": shrinkage,
                    "metrics": metrics(truth, prediction),
                    "score": _score(truth, prediction),
                }
            )
    ranking.sort(key=lambda item: tuple(item["score"]), reverse=True)
    return ranking


def _ids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def filter_scene_cloud_cover(
    samples: list[CalibrationSample], max_percent: float
) -> list[CalibrationSample]:
    """Apply the benchmark's catalogue-level cloud-cover eligibility rule."""
    return [
        sample
        for sample in samples
        if sample.features.get("metadata_scene_cloud_cover", float("inf")) < max_percent
    ]


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
    parser.add_argument("--allow-incomplete-train", action="store_true")
    parser.add_argument(
        "--train-max-scene-cloud-cover",
        type=float,
        default=DEFAULT_MAX_SCENE_CLOUD_COVER,
    )
    parser.add_argument("--include-target-diagnostic-ranking", action="store_true")
    parser.add_argument("--select-for-target-retrieval-distribution", action="store_true")
    parser.add_argument("--fit-for-target-retrieval-distribution", action="store_true")
    parser.add_argument(
        "--domain-weight-method",
        choices=("baseline_aod", "operational"),
        default="baseline_aod",
    )
    args = parser.parse_args()

    matchup_path = args.root / "matchups" / "matchups.csv"
    common = {"include_geography": False}
    train = load_samples(
        args.train_result_dir,
        args.train_context_dir,
        matchup_path,
        _ids(args.train_mids),
        atmo_context_dir=args.train_atmo_context_dir,
        require_complete=not args.allow_incomplete_train,
        **common,
    )
    train_before_cloud_filter = len(train)
    train = filter_scene_cloud_cover(train, args.train_max_scene_cloud_cover)
    if not train:
        raise ValueError("No external training samples pass the scene-cloud-cover filter")
    target = load_samples(
        args.target_result_dir,
        args.target_context_dir,
        matchup_path,
        _ids(args.target_mids),
        atmo_context_dir=args.target_atmo_context_dir,
        **common,
    )
    recipe, ranking = select_recipe(
        train,
        score_domain_samples=(
            target if args.select_for_target_retrieval_distribution else None
        ),
        fit_domain_samples=(target if args.fit_for_target_retrieval_distribution else None),
        domain_weight_method=args.domain_weight_method,
    )
    development = [sample for sample in train if site_fold(sample.site) != HOLDOUT_FOLD]
    holdout = [sample for sample in train if site_fold(sample.site) == HOLDOUT_FOLD]
    analysis = {
        "schema_version": "siac-generic-aod-calibrator-v1",
        "selection_policy": {
            "target_used_for_selection": False,
            "geography_features": False,
            "site_holdout_fold": HOLDOUT_FOLD,
            "development_count": len(development),
            "holdout_count": len(holdout),
            "train_count_before_cloud_filter": train_before_cloud_filter,
            "train_count_after_cloud_filter": len(train),
            "train_max_scene_cloud_cover": args.train_max_scene_cloud_cover,
            "target_hits_required": TARGET_HITS,
            "target_labels_used_for_selection": False,
            "target_operational_covariates_used_for_domain_weighting": (
                args.select_for_target_retrieval_distribution
                or args.fit_for_target_retrieval_distribution
            ),
            "target_retrieved_aod_distribution_used_for_selection_scoring": (
                args.select_for_target_retrieval_distribution
            ),
            "target_retrieved_aod_distribution_used_for_external_model_fit": (
                args.fit_for_target_retrieval_distribution
            ),
            "domain_weight_method": args.domain_weight_method,
        },
        "selected_recipe": {
            "name": recipe.name,
            "candidate": recipe.candidate.name,
            "family": recipe.candidate.family,
            "target": recipe.candidate.target,
            "select_k": recipe.candidate.select_k,
            "parameters": recipe.candidate.parameter_dict(),
            "shrinkage": recipe.shrinkage,
        },
        "development_ranking": ranking,
        "external_holdout": evaluate_recipe(
            recipe,
            development,
            holdout,
            fit_domain_samples=(
                target if args.fit_for_target_retrieval_distribution else None
            ),
            metric_domain_samples=target,
            domain_weight_method=args.domain_weight_method,
        ),
        "target": evaluate_recipe(
            recipe,
            train,
            target,
            fit_domain_samples=(
                target if args.fit_for_target_retrieval_distribution else None
            ),
            domain_weight_method=args.domain_weight_method,
        ),
    }
    if args.include_target_diagnostic_ranking:
        analysis["target_diagnostic_ranking"] = {
            "used_for_selection": False,
            "warning": "Diagnostic transfer audit only; target labels did not select the recipe.",
            "ranking": diagnostic_target_ranking(
                train,
                target,
                fit_domain_samples=(
                    target if args.fit_for_target_retrieval_distribution else None
                ),
                domain_weight_method=args.domain_weight_method,
            ),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {key: analysis[key] for key in ("selected_recipe", "external_holdout", "target")},
            indent=2,
        )
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
