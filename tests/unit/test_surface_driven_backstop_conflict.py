"""Unit tests for the operational surface/prior conflict gate."""

from __future__ import annotations

import numpy as np
import pytest

from siac.algorithms.solver.surface_driven import _backstop_conflict_decision


def test_positive_surface_prior_conflict_relaxes_to_flat_width() -> None:
    use_flat, diagnostics = _backstop_conflict_decision(
        surface_min_aot=0.60,
        aot_prior=np.full((2, 2), 0.30, dtype=np.float32),
        calibrated_unc=np.full((2, 2), 0.07, dtype=np.float32),
        z_threshold=2.576,
    )

    assert use_flat is True
    assert diagnostics["surface_backstop_conflict_decision"] == "flat50"
    assert diagnostics["surface_backstop_conflict_standardized_positive"] == pytest.approx(
        0.30 / 0.07
    )


def test_small_or_unconfigured_conflict_keeps_calibrated_width() -> None:
    use_flat, diagnostics = _backstop_conflict_decision(
        surface_min_aot=0.36,
        aot_prior=np.full((2, 2), 0.30, dtype=np.float32),
        calibrated_unc=np.full((2, 2), 0.07, dtype=np.float32),
        z_threshold=2.576,
    )

    assert use_flat is False
    assert diagnostics["surface_backstop_conflict_decision"] == "calibrated"

    unconfigured, unconfigured_diagnostics = _backstop_conflict_decision(
        surface_min_aot=1.0,
        aot_prior=np.full((2, 2), 0.30, dtype=np.float32),
        calibrated_unc=np.full((2, 2), 0.07, dtype=np.float32),
        z_threshold=None,
    )
    assert unconfigured is False
    assert unconfigured_diagnostics["surface_backstop_conflict_gate_configured"] is False
