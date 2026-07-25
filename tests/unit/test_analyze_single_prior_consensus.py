from __future__ import annotations

import numpy as np
from tools.aeronet_validation.analyze_single_prior_consensus import (
    _apply_gate,
    _feature_values,
    _hit,
    _site_fold,
)


def test_hit_uses_expected_error_envelope() -> None:
    assert _hit(0.20, 0.20)
    assert not _hit(0.40, 0.20)


def test_site_fold_is_stable() -> None:
    assert _site_fold("example") == _site_fold("example")


def test_apply_gate_switches_only_matching_rows() -> None:
    rows = [
        {"baseline": 0.1, "cams": 0.2, "band_mean": 0.6, "cost": 1.0},
        {"baseline": 0.2, "cams": 0.3, "band_mean": 0.7, "cost": 3.0},
    ]
    gate = {
        "aggregation": "mean",
        "weight": 0.25,
        "feature": "cost",
        "direction": "high",
        "threshold": 2.0,
    }

    prediction = _apply_gate(rows, np.array([0, 1]), gate)

    np.testing.assert_allclose(prediction, [0.1, 0.4])


def test_feature_values_computes_posterior_disagreement() -> None:
    rows = [{"baseline": 0.1}, {"baseline": 0.4}]

    values = _feature_values(rows, np.array([0, 1]), "posterior_delta", np.array([0.3, 0.2]))

    np.testing.assert_allclose(values, [0.2, 0.2])
