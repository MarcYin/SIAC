from __future__ import annotations

import numpy as np
import pytest

from siac.algorithms.aod_calibration import GenericAODCalibrator


class _FixedEstimator:
    def __init__(self, values: list[float]) -> None:
        self.values = np.asarray(values, dtype=float)

    def predict(self, matrix: np.ndarray) -> np.ndarray:
        assert matrix.shape == (len(self.values), 2)
        return self.values


def test_generic_aod_calibrator_applies_log_ratio_and_global_offset() -> None:
    model = GenericAODCalibrator(
        feature_names=("a", "b"),
        estimator=_FixedEstimator([0.0, np.log(2.0)]),
        global_log_offset=np.log(0.5),
    )
    prediction = model.predict([0.2, 0.2], [{"a": 1.0}, {"a": 2.0, "b": 3.0}])
    expected = np.asarray([0.5 * (0.2 + 1 / 3) - 1 / 3, 0.2])
    assert np.allclose(prediction, np.clip(expected, 0.0, 4.0))
    assert np.isnan(model.feature_matrix([{"a": 1.0}])[0, 1])


def test_generic_aod_calibrator_rejects_mismatched_rows() -> None:
    model = GenericAODCalibrator(("a", "b"), _FixedEstimator([0.0]))
    with pytest.raises(ValueError, match="counts do not agree"):
        model.predict([0.1], [])


def test_generic_aod_calibrator_clips_model_before_global_offset() -> None:
    model = GenericAODCalibrator(("a", "b"), _FixedEstimator([-10.0]), global_log_offset=0.1)
    prediction = model.predict([0.0], [{"a": 1.0}])
    assert prediction[0] == pytest.approx((np.exp(0.1) - 1.0) / 3.0)
