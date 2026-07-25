from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from tools.aeronet_validation.build_bright_surface_bias_report import (
    _jsonable,
    _metrics,
    _pair_case_metrics,
    _screen_solution,
)

if TYPE_CHECKING:
    from pathlib import Path


def test_independent_pixel_screen_respects_the_gaussian_aod_prior() -> None:
    axis = np.asarray([0.1, 0.3])
    cost = np.asarray([[[0.4, 0.0]], [[0.0, 0.0]]])
    aot_prior = np.asarray([[0.1, 0.1]])
    aot_unc = np.asarray([[0.05, 0.05]])
    valid = np.asarray([[True, False]])

    result = _screen_solution(cost, axis, aot_prior, aot_unc, valid)

    assert result == 0.1


def test_metrics_and_jsonable_preserve_missing_values() -> None:
    cases = [
        {"truth": 0.2, "retrieved": 0.25},
        {"truth": 0.5, "retrieved": None},
        {"truth": 0.5, "retrieved": 0.7},
    ]

    metrics = _metrics(cases, "retrieved")

    assert metrics["n"] == 2
    assert metrics["within_ee"] == 1
    assert metrics["within_ee_percent"] == 50.0
    assert _jsonable({"values": np.asarray([1.0, np.nan])}) == {"values": [1.0, None]}


def test_pair_metric_uses_mean_surface_brightness_as_the_primary_axis(tmp_path: Path) -> None:
    l2a_brightness = np.linspace(0.1, 0.4, 30)
    l2a = np.repeat(l2a_brightness[:, np.newaxis], 3, axis=1)
    current_rt = 0.8 * l2a
    path = tmp_path / "pair.npz"
    np.savez_compressed(
        path,
        band_names=np.asarray(["blue", "green", "red"]),
        l2a=l2a,
        siac=current_rt,
    )

    metrics, sample = _pair_case_metrics(path)

    assert sample is not None
    np.testing.assert_allclose(sample[0], 0.9 * l2a_brightness)
    assert metrics["l2a_minus_current_rt_vs_mean_slope"] == pytest.approx(2.0 / 9.0)
    assert metrics["l2a_minus_current_rt_slope"] == pytest.approx(0.2)
