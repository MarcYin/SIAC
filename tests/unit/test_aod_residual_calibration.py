"""Tests for leakage-safe AOD residual calibration helpers."""

from __future__ import annotations

import numpy as np
from tools.aeronet_validation.aod_residual_calibration import (
    CalibrationSample,
    UniformResidualCalibrator,
    apply_guard,
    extract_operational_features,
    metrics,
    site_fold,
)


def _result() -> dict[str, object]:
    return {
        "retrieved": 0.2,
        "scene_mean": 0.2,
        "truth": 0.9,
        "err": -0.7,
        "within_ee": False,
        "ee_threshold": 0.185,
        "lat": 51.0,
        "lon": -1.0,
        "prior_boa": {"B02": 0.04, "B03": 0.06, "B04": 0.08},
        "prior_unc": {"B02": 0.01, "B03": 0.02, "B04": 0.02},
        "anchor_iterate": {"pass1_scene_mean": 0.21, "pass2_scene_mean": 0.2},
        "solver": {"cost_final": 1.5, "converged": True},
    }


def _context() -> dict[str, object]:
    return {
        "truth": 0.9,
        "sensing_time_utc": "2024-06-01T10:30:00",
        "values": {
            "cams_total_aerosol_optical_depth_at_550nm_surface": 0.25,
            "cams_dust_aerosol_optical_depth_at_550nm_surface": 0.02,
            "merra_TOTEXTTAU": 0.8,
        },
    }


def _matchup() -> dict[str, object]:
    return {
        "aeronet_aod550_mean": 0.9,
        "sensing_time_utc": "2024-06-01T10:30:00",
        "satellite": "S2A",
        "scene_cloud_cover": 4.0,
        "elevation_m": 100.0,
    }


def test_feature_extraction_ignores_all_target_derived_fields() -> None:
    spectral = {
        "truth": 0.9,
        "ref": {"B02": 0.5},
        "monthly_err": "target diagnostic",
        "sza": 30.0,
        "vza": 4.0,
        "raa": 80.0,
        "toa": {"B02": 0.1, "B03": 0.12, "B04": 0.09, "B8A": 0.3},
    }
    first = extract_operational_features(
        _result(), _context(), _matchup(), spectral_context=spectral
    )

    changed_result = _result()
    changed_result.update({"truth": 0.01, "err": 99.0, "within_ee": True, "ee_threshold": 9.0})
    changed_context = _context()
    changed_context["truth"] = 0.01
    changed_matchup = _matchup()
    changed_matchup["aeronet_aod550_mean"] = 0.01
    changed_spectral = dict(spectral)
    changed_spectral.update({"truth": 0.01, "ref": {"B02": 0.001}})
    second = extract_operational_features(
        changed_result,
        changed_context,
        changed_matchup,
        spectral_context=changed_spectral,
    )

    assert first == second
    assert not any("truth" in name or "aeronet" in name or "ref" in name for name in first)
    assert "context_merra_TOTEXTTAU" not in first
    assert first["geometry_sza"] == 30.0
    assert first["toa_B02"] == 0.1


def test_site_fold_is_stable() -> None:
    assert [site_fold(site) for site in ("A", "B", "C", "D", "E")] == [1, 0, 2, 3, 2]


def test_high_aod_guard_retains_baseline() -> None:
    baseline = np.asarray([0.2, 0.8, 1.0])
    candidate = np.asarray([0.3, 0.75, 0.76])
    np.testing.assert_allclose(apply_guard(baseline, candidate), [0.3, 0.75, 1.0])


def test_metrics_use_expected_error_envelope() -> None:
    summary = metrics(np.asarray([0.0, 1.0]), np.asarray([0.05, 1.21]))
    assert summary["hits"] == 1
    assert summary["count"] == 2


def test_uniform_calibrator_applies_high_aod_disagreement_safeguard() -> None:
    class ZeroResidualModel:
        def predict(self, matrix):  # noqa: ANN001
            return np.zeros(matrix.shape[0])

    sample = CalibrationSample(
        matchup_id="case",
        site="site",
        retrieved=0.2,
        truth=1.0,
        features={
            "cams": 0.6,
            "context_cams_total_aerosol_optical_depth_at_550nm_surface": 0.6,
            "solver_surface_band_argmin_spread": 0.5,
        },
    )
    calibrator = UniformResidualCalibrator(
        feature_names=("cams",),
        models=(ZeroResidualModel(), ZeroResidualModel(), ZeroResidualModel()),
    )

    prediction, high_disagreement = calibrator.predict([sample])

    np.testing.assert_allclose(prediction, [1.08])
    np.testing.assert_array_equal(high_disagreement, [True])
