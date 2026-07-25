"""Tests for native-grid MAIAC staging helpers."""

from __future__ import annotations

from datetime import date

import numpy as np
import pytest
import xarray as xr
from tools.aeronet_validation.stage_maiac_prior import (
    calibrated_uncertainty,
    native_window_values,
    select_day,
)


def _array(values: list[list[float]]) -> xr.DataArray:
    array = xr.DataArray(
        np.asarray(values, dtype=np.float32),
        dims=("y", "x"),
        coords={"y": [1.0, 0.0, -1.0], "x": [-1.0, 0.0, 1.0]},
    )
    return array.rio.write_crs("EPSG:4326")


def test_native_window_values_aggregates_before_sparse_resampling() -> None:
    aot = _array([[np.nan, 0.1, np.nan], [0.2, 0.3, 0.4], [np.nan, 0.5, np.nan]])
    unc = _array([[np.nan, 0.01, np.nan], [0.02, 0.03, np.nan], [np.nan, 0.05, np.nan]])

    aot_values, unc_values = native_window_values(aot, unc, (-0.5, -0.5, 0.5, 0.5))

    np.testing.assert_allclose(aot_values, [0.3])
    np.testing.assert_allclose(unc_values, [0.03])


def test_select_day_prefers_nearest_then_larger_support() -> None:
    values = {
        date(2024, 1, 9): ([0.1, 0.2], [0.01, 0.01]),
        date(2024, 1, 11): ([0.3, 0.4, 0.5], [0.02, 0.02, 0.02]),
        date(2024, 1, 12): ([0.9] * 20, [0.1] * 20),
    }

    selected = select_day(values, date(2024, 1, 10))

    assert selected is not None
    day, aot, uncertainty = selected
    assert day == date(2024, 1, 11)
    np.testing.assert_allclose(aot, [0.3, 0.4, 0.5])
    np.testing.assert_allclose(uncertainty, [0.02, 0.02, 0.02])


def test_calibrated_uncertainty_has_native_high_aod_and_sparse_terms() -> None:
    assert calibrated_uncertainty(0.1, 0.01, 20) == 0.04
    assert calibrated_uncertainty(0.5, 0.08, 20) == pytest.approx(0.116)
    assert calibrated_uncertainty(0.1, 0.01, 3) == 0.09
