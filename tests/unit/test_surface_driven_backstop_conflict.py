"""Unit tests for the surface/prior conflict diagnostic."""

from __future__ import annotations

import numpy as np
import pytest

from siac.algorithms.solver.surface_driven import _backstop_conflict_diagnostics


def test_conflict_is_reported_in_units_of_the_calibrated_sigma() -> None:
    diagnostics = _backstop_conflict_diagnostics(
        surface_min_aot=0.60,
        aot_prior=np.full((2, 2), 0.30, dtype=np.float32),
        calibrated_unc=np.full((2, 2), 0.07, dtype=np.float32),
    )

    assert diagnostics["surface_backstop_conflict_standardized_positive"] == pytest.approx(
        0.30 / 0.07
    )
    assert diagnostics["surface_backstop_conflict_surface_min_aot"] == pytest.approx(0.60)
    assert diagnostics["surface_backstop_conflict_prior_aot_median"] == pytest.approx(0.30)
    assert diagnostics["surface_backstop_conflict_calibrated_sigma_median"] == pytest.approx(0.07)


def test_conflict_statistic_is_none_when_inputs_are_unusable() -> None:
    diagnostics = _backstop_conflict_diagnostics(
        surface_min_aot=float("nan"),
        aot_prior=np.full((2, 2), np.nan, dtype=np.float32),
        calibrated_unc=np.zeros((2, 2), dtype=np.float32),
    )

    assert diagnostics["surface_backstop_conflict_standardized_positive"] is None
    assert diagnostics["surface_backstop_conflict_prior_aot_median"] is None
    assert diagnostics["surface_backstop_conflict_calibrated_sigma_median"] is None
