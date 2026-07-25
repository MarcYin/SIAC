from __future__ import annotations

import numpy as np
from tools.aeronet_validation.aod_residual_calibration import CalibrationSample
from tools.aeronet_validation.evaluate_generic_aod_gate import gate_matrix, gate_target


def _samples() -> list[CalibrationSample]:
    return [
        CalibrationSample(
            matchup_id=f"m{index}",
            site=f"s{index}",
            retrieved=0.1 + index * 0.1,
            truth=0.15 + index * 0.1,
            features={
                "x": float(index),
                "context_cams_total_aerosol_optical_depth_at_550nm_surface": 0.2,
                "atmo_aot_mean": 0.18,
            },
        )
        for index in range(3)
    ]


def test_gate_matrix_uses_only_operational_inputs() -> None:
    samples = _samples()
    matrix = gate_matrix(samples, ("x",), np.array([0.14, 0.25, 0.32]))
    assert matrix.shape == (3, 9)
    assert np.all(np.isfinite(matrix))


def test_gate_target_is_normalized_candidate_error_delta() -> None:
    samples = _samples()
    perfect = np.asarray([sample.truth for sample in samples])
    target = gate_target(samples, perfect)
    assert np.all(target < 0.0)
