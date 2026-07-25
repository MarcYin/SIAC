from __future__ import annotations

import numpy as np
from tools.aeronet_validation.aod_residual_calibration import CalibrationSample
from tools.aeronet_validation.evaluate_generic_aod_classifier_gate import (
    gate_preference_target,
)


def _samples() -> list[CalibrationSample]:
    return [
        CalibrationSample(
            matchup_id=f"m{index}",
            site=f"s{index}",
            retrieved=retrieved,
            truth=truth,
            features={"x": float(index)},
        )
        for index, (retrieved, truth) in enumerate(((0.1, 0.2), (0.4, 0.3), (0.5, 0.5)))
    ]


def test_gate_preference_target_marks_more_accurate_candidate() -> None:
    samples = _samples()
    candidate = np.asarray([0.2, 0.6, 0.55])
    preference, weights = gate_preference_target(samples, candidate)
    np.testing.assert_array_equal(preference, [1, 0, 0])
    assert np.all(np.isfinite(weights))
    assert np.all(weights > 0.0)


def test_hit_changing_preferences_receive_more_weight() -> None:
    samples = _samples()
    candidate = np.asarray([0.2, 0.6, 0.5])
    _, weights = gate_preference_target(samples, candidate)
    assert weights[0] > weights[2]
