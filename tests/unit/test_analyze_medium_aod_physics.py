from __future__ import annotations

import pytest
from tools.aeronet_validation.analyze_medium_aod_physics import (
    _paired_transitions,
    _subset_ids,
    compute_metrics,
    ee_threshold,
    within_ee,
)


def test_expected_error_definition_is_inclusive() -> None:
    truth = 0.4
    threshold = ee_threshold(truth)

    assert threshold == pytest.approx(0.11)
    assert within_ee(truth - threshold, truth)
    assert not within_ee(truth - threshold - 1e-6, truth)


def test_compute_metrics_counts_under_and_over_failures() -> None:
    metrics = compute_metrics(
        [
            {"truth": 0.4, "retrieved": 0.28},
            {"truth": 0.6, "retrieved": 0.75},
        ]
    )

    assert metrics.n == 2
    assert metrics.within_ee == 0
    assert metrics.within_ee_rate == 0.0
    assert metrics.bias == pytest.approx(0.015)
    assert metrics.mae == pytest.approx(0.135)
    assert metrics.rmse == pytest.approx((0.12**2 / 2 + 0.15**2 / 2) ** 0.5)
    assert metrics.median_error == pytest.approx(0.015)
    assert metrics.under_ee == 1
    assert metrics.over_ee == 1


def test_medium_subset_includes_requested_boundaries() -> None:
    truth = {"low": 0.249, "lower": 0.25, "upper": 0.85, "high": 0.851}
    subsets = _subset_ids(truth, {"upper"})

    assert subsets["low"] == ["low"]
    assert subsets["medium"] == ["lower", "upper"]
    assert subsets["high"] == ["high"]
    assert "upper" not in subsets["non_extreme"]


def test_paired_transitions_report_gains_and_losses() -> None:
    baseline = {
        "gain": {"truth": 0.4, "retrieved": 0.28},
        "loss": {"truth": 0.6, "retrieved": 0.6},
        "hit": {"truth": 0.2, "retrieved": 0.2},
        "miss": {"truth": 0.3, "retrieved": 0.1},
    }
    candidate = {
        "gain": {"truth": 0.4, "retrieved": 0.35},
        "loss": {"truth": 0.6, "retrieved": 0.8},
        "hit": {"truth": 0.2, "retrieved": 0.21},
        "miss": {"truth": 0.3, "retrieved": 0.1},
    }

    result = _paired_transitions(baseline, candidate, baseline)

    assert result["n_paired"] == 4
    assert result["gains"] == 1
    assert result["losses"] == 1
    assert result["net_gains"] == 0
    assert result["stable_hits"] == 1
    assert result["stable_misses"] == 1
