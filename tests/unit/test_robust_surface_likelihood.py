from __future__ import annotations

import numpy as np
import pytest
from tools.aeronet_validation.evaluate_robust_surface_likelihood import (
    aggregate_curve,
    interpolate_curve_minimum,
    robust_loss,
)


def test_robust_losses_downweight_large_standardized_residuals():
    values = np.asarray([1.0, 10.0])
    chi2 = robust_loss(values, "chi2")
    huber = robust_loss(values, "huber1p5")
    student = robust_loss(values, "student3")
    assert chi2[0] == pytest.approx(0.5)
    assert huber[1] < chi2[1]
    assert student[1] < huber[1]


def test_aggregate_curve_supports_mean_and_median():
    loss = np.asarray([[1.0, 2.0, 100.0], [2.0, 3.0, 4.0]])
    assert aggregate_curve(loss, "mean").tolist() == pytest.approx([103 / 3, 3.0])
    assert aggregate_curve(loss, "median").tolist() == pytest.approx([2.0, 3.0])


def test_interpolate_curve_minimum_on_nonuniform_axis():
    axis = np.asarray([0.1, 0.2, 0.4])
    curve = np.square(axis - 0.25) + 1.0
    assert interpolate_curve_minimum(axis, curve, 1) == pytest.approx(0.25)


def test_interpolate_curve_minimum_keeps_edge_node():
    axis = np.asarray([0.1, 0.2, 0.4])
    curve = np.asarray([1.0, 2.0, 3.0])
    assert interpolate_curve_minimum(axis, curve, 0) == pytest.approx(0.1)
