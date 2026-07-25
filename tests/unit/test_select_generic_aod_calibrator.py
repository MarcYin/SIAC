from __future__ import annotations

import numpy as np
from tools.aeronet_validation.aod_residual_calibration import CalibrationSample
from tools.aeronet_validation.select_generic_aod_calibrator import (
    CandidateSpec,
    PredictionRecipe,
    _apply_recipe,
    baseline_domain_weights,
    candidate_specs,
    evaluate_recipe,
    filter_scene_cloud_cover,
    operational_domain_weights,
)


def _samples(offset: float) -> list[CalibrationSample]:
    return [
        CalibrationSample(
            matchup_id=f"m{index}",
            site=f"site{index}",
            retrieved=0.1 + 0.05 * index,
            truth=0.1 + 0.05 * index + offset,
            features={"x": float(index), "cams": 0.11 + 0.05 * index},
        )
        for index in range(8)
    ]


def test_apply_recipe_is_uniform_shrinkage() -> None:
    baseline = np.array([0.1, 0.4])
    raw = np.array([0.3, 0.2])
    np.testing.assert_allclose(_apply_recipe(baseline, raw, 0.5), [0.2, 0.3])


def test_evaluate_recipe_has_complete_prediction_receipts() -> None:
    recipe = PredictionRecipe(
        CandidateSpec("ridge", "ridge", select_k=2, params=(("alpha", 1.0),)),
        shrinkage=0.75,
    )
    train = _samples(0.03)
    test = _samples(0.02)
    result = evaluate_recipe(recipe, train, test)
    assert result["count"] == len(test)
    assert len(result["predictions"]) == len(test)
    assert result["candidate"]["count"] == len(test)


def test_filter_scene_cloud_cover_matches_strict_benchmark_rule() -> None:
    samples = _samples(0.0)
    samples[0].features["metadata_scene_cloud_cover"] = 19.99
    samples[1].features["metadata_scene_cloud_cover"] = 20.0
    filtered = filter_scene_cloud_cover(samples, 20.0)
    assert [sample.matchup_id for sample in filtered] == ["m0"]


def test_log_ratio_recipe_returns_finite_nonnegative_aod() -> None:
    recipe = PredictionRecipe(
        CandidateSpec(
            "ridge_log_ratio",
            "ridge",
            target="log_ratio",
            select_k=2,
            params=(("alpha", 1.0),),
        ),
        shrinkage=1.0,
    )
    result = evaluate_recipe(recipe, _samples(0.03), _samples(0.02))
    predictions = np.asarray([row["candidate"] for row in result["predictions"]])
    assert np.all(np.isfinite(predictions))
    assert np.all(predictions >= 0.0)


def test_baseline_domain_weights_are_finite_and_normalized() -> None:
    source = _samples(0.0)
    domain = _samples(0.0)[:4]
    weights = baseline_domain_weights(source, domain)
    assert np.all(np.isfinite(weights))
    assert np.all(weights > 0.0)
    np.testing.assert_allclose(np.mean(weights), 1.0)


def test_evaluate_recipe_reports_domain_weighted_fit_and_metrics() -> None:
    recipe = PredictionRecipe(
        CandidateSpec("ridge", "ridge", select_k=2, params=(("alpha", 1.0),)),
        shrinkage=0.75,
    )
    train = _samples(0.03)
    test = _samples(0.02)
    domain = _samples(0.01)[4:]
    result = evaluate_recipe(
        recipe,
        train,
        test,
        fit_domain_samples=domain,
        metric_domain_samples=domain,
    )
    assert result["domain_weighted_baseline"]["within_ee_percent"] >= 0.0
    assert result["domain_weighted_candidate"]["within_ee_percent"] >= 0.0
    assert len(result["predictions"]) == len(test)


def test_candidate_grid_includes_regularized_lightgbm_log_ratio_models() -> None:
    lightgbm = [spec for spec in candidate_specs() if spec.family.startswith("lightgbm")]
    assert lightgbm
    assert {spec.target for spec in lightgbm} == {"log_ratio"}
    assert {spec.select_k for spec in lightgbm} == {35, 50}


def test_candidate_grid_includes_asymmetric_log_ratio_quantiles() -> None:
    quantiles = [spec for spec in candidate_specs() if spec.family == "hist_quantile"]
    assert quantiles
    assert {spec.target for spec in quantiles} == {"log_ratio"}
    assert {spec.parameter_dict()["quantile"] for spec in quantiles} == {
        0.35,
        0.4,
        0.45,
        0.55,
        0.6,
    }


def test_operational_domain_weights_ignore_aeronet_targets() -> None:
    source = _samples(0.0)
    domain = _samples(0.0)[2:]
    weights = operational_domain_weights(source, domain)
    changed_truth = [
        CalibrationSample(
            matchup_id=sample.matchup_id,
            site=sample.site,
            retrieved=sample.retrieved,
            truth=sample.truth + 10.0,
            features=sample.features,
        )
        for sample in domain
    ]
    changed_weights = operational_domain_weights(source, changed_truth)
    np.testing.assert_allclose(weights, changed_weights)
    np.testing.assert_allclose(np.mean(weights), 1.0)
    assert np.all(weights > 0.0)
