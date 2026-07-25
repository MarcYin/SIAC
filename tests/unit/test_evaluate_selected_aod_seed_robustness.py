from __future__ import annotations

import numpy as np
from tools.aeronet_validation.aod_residual_calibration import CalibrationSample
from tools.aeronet_validation.evaluate_selected_aod_seed_robustness import evaluate_slice


def test_evaluate_slice_respects_requested_indices() -> None:
    samples = [
        CalibrationSample(f"m{i}", f"s{i}", baseline, truth, {"x": float(i)})
        for i, (baseline, truth) in enumerate(((0.1, 0.1), (0.1, 0.5), (0.8, 0.8)))
    ]
    prediction = np.asarray([0.1, 0.5, 0.0])
    result = evaluate_slice(samples, prediction, np.asarray([0, 1]))
    assert result["count"] == 2
    assert result["candidate"]["hits"] == 2
