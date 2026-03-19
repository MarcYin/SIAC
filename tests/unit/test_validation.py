"""
Layer 2 — Validation functions for pipeline contract outputs.
"""

import dataclasses

import numpy as np
import pytest
import xarray as xr

from siac.domain.validation import (
    _validate_atmospheric_state,
    _validate_correction_result,
    _validate_observation_bundle,
    _validate_solved_atmosphere,
    _validate_solver_input_bundle,
    _validate_surface_prior,
)
from siac.errors import ValidationError

# ── _validate_observation_bundle ──────────────────────────────────────

class TestValidateObservationBundle:
    def test_validate_obs_happy(self, mock_observation_bundle):
        _validate_observation_bundle(mock_observation_bundle)  # no error

    def test_validate_obs_missing_time(self, mock_observation_bundle):
        obs = mock_observation_bundle
        # Create a copy with observation_time removed
        bad_meta = {k: v for k, v in obs.metadata.items() if k != "observation_time"}
        bad_obs = dataclasses.replace(obs, metadata=bad_meta)
        with pytest.raises(ValidationError, match="observation_time"):
            _validate_observation_bundle(bad_obs)

    def test_validate_obs_wrong_time_type(self, mock_observation_bundle):
        obs = mock_observation_bundle
        bad_meta = dict(obs.metadata)
        bad_meta["observation_time"] = "2024-01-01"
        bad_obs = dataclasses.replace(obs, metadata=bad_meta)
        with pytest.raises(ValidationError, match="datetime"):
            _validate_observation_bundle(bad_obs)

    def test_validate_obs_empty_toa(self, mock_observation_bundle):
        obs = mock_observation_bundle
        bad_obs = dataclasses.replace(obs, toa=xr.Dataset())
        with pytest.raises(ValidationError, match="at least one band"):
            _validate_observation_bundle(bad_obs)

    def test_validate_obs_cloud_shape_mismatch(self, mock_observation_bundle):
        obs = mock_observation_bundle
        bad_mask = xr.DataArray(np.zeros((10, 10), dtype=bool), dims=["y", "x"])
        bad_obs = dataclasses.replace(obs, cloud_mask=bad_mask)
        with pytest.raises(ValidationError, match="cloud_mask shape"):
            _validate_observation_bundle(bad_obs)


# ── _validate_atmospheric_state ───────────────────────────────────────

class TestValidateAtmosphericState:
    def test_validate_atmo_happy(self, mock_atmospheric_state):
        _validate_atmospheric_state(mock_atmospheric_state)  # no error

    def test_validate_atmo_negative_aot_unc(self, mock_atmospheric_state):
        atmo = mock_atmospheric_state
        bad_unc = xr.DataArray(np.full(atmo.aot_unc.shape, -0.1), dims=["y", "x"])
        bad_atmo = dataclasses.replace(atmo, aot_unc=bad_unc)
        with pytest.raises(ValidationError, match="non-negative"):
            _validate_atmospheric_state(bad_atmo)

    def test_validate_atmo_negative_tcwv_unc(self, mock_atmospheric_state):
        atmo = mock_atmospheric_state
        bad_unc = xr.DataArray(np.full(atmo.tcwv_unc.shape, -0.5), dims=["y", "x"])
        bad_atmo = dataclasses.replace(atmo, tcwv_unc=bad_unc)
        with pytest.raises(ValidationError, match="non-negative"):
            _validate_atmospheric_state(bad_atmo)

    def test_validate_atmo_negative_tco3_unc(self, mock_atmospheric_state):
        atmo = mock_atmospheric_state
        bad_unc = xr.DataArray(np.full(atmo.tco3_unc.shape, -0.01), dims=["y", "x"])
        bad_atmo = dataclasses.replace(atmo, tco3_unc=bad_unc)
        with pytest.raises(ValidationError, match="non-negative"):
            _validate_atmospheric_state(bad_atmo)


# ── _validate_surface_prior ───────────────────────────────────────────

class TestValidateSurfacePrior:
    def test_validate_prior_happy(self, mock_surface_prior):
        _validate_surface_prior(mock_surface_prior)  # no error

    def test_validate_prior_shape_mismatch(self, mock_surface_prior):
        sp = mock_surface_prior
        bad_unc = xr.DataArray(np.ones((10, 10)), dims=["y", "x"])
        bad_sp = dataclasses.replace(sp, boa_unc=bad_unc)
        with pytest.raises(ValidationError, match="boa shape"):
            _validate_surface_prior(bad_sp)

    def test_validate_prior_mask_broadcast(self, mock_surface_prior):
        sp = mock_surface_prior
        # mask shape that can't broadcast to boa shape
        bad_mask = xr.DataArray(np.ones((5, 7), dtype=bool), dims=["y", "x"])
        bad_sp = dataclasses.replace(sp, mask=bad_mask)
        with pytest.raises(ValidationError, match="broadcastable"):
            _validate_surface_prior(bad_sp)


# ── _validate_solver_input_bundle ─────────────────────────────────────

class TestValidateSolverInputBundle:
    def test_validate_sib_happy(self, mock_solver_input_bundle):
        _validate_solver_input_bundle(mock_solver_input_bundle)  # no error

    def test_validate_sib_bands_not_in_config(self, mock_solver_input_bundle):
        from siac.domain import SensorBand
        sib = mock_solver_input_bundle
        bad_band = SensorBand("FAKE", 9999.0, 10.0, 10.0, 99)
        bad_sib = dataclasses.replace(sib, bands=[bad_band])
        with pytest.raises(ValidationError, match="not in sensor_config"):
            _validate_solver_input_bundle(bad_sib)

    def test_validate_sib_positive_resolution(self, mock_solver_input_bundle):
        sib = mock_solver_input_bundle
        assert sib.aux_resolution_m > 0
        assert sib.aerosol_resolution_m > 0


# ── _validate_solved_atmosphere ───────────────────────────────────────

class TestValidateSolvedAtmosphere:
    def test_validate_solved_happy(self, mock_solved_atmosphere):
        _validate_solved_atmosphere(mock_solved_atmosphere)  # no error

    def test_validate_solved_negative_aot(self, mock_solved_atmosphere):
        solved = mock_solved_atmosphere
        bad_aot = xr.DataArray(np.full(solved.aot.shape, -0.1), dims=["y", "x"])
        bad_solved = dataclasses.replace(solved, aot=bad_aot)
        with pytest.raises(ValidationError, match="non-negative"):
            _validate_solved_atmosphere(bad_solved)

    def test_validate_solved_negative_tcwv(self, mock_solved_atmosphere):
        solved = mock_solved_atmosphere
        bad_tcwv = xr.DataArray(np.full(solved.tcwv.shape, -0.5), dims=["y", "x"])
        bad_solved = dataclasses.replace(solved, tcwv=bad_tcwv)
        with pytest.raises(ValidationError, match="non-negative"):
            _validate_solved_atmosphere(bad_solved)

    def test_validate_solved_converged_type(self, mock_solved_atmosphere):
        """converged must be a boolean."""
        solved = mock_solved_atmosphere
        assert isinstance(solved.converged, bool)

    def test_validate_solved_diagnostics_types(self, mock_solved_atmosphere):
        """n_iterations must be int, cost_final must be numeric."""
        solved = mock_solved_atmosphere
        assert isinstance(solved.n_iterations, int)
        assert isinstance(solved.cost_final, (int, float))


# ── _validate_correction_result ───────────────────────────────────────

class TestValidateCorrectionResult:
    def test_validate_result_happy(self, mock_observation_bundle, mock_solved_atmosphere):
        from siac.domain import CorrectionResult
        result = CorrectionResult(
            boa=mock_observation_bundle.toa,
            boa_unc=None,
            aot=mock_solved_atmosphere.aot,
            tcwv=mock_solved_atmosphere.tcwv,
            cloud_mask=mock_observation_bundle.cloud_mask,
            metadata={"processing_time_s": 0.01},
        )
        _validate_correction_result(result)  # no error

    def test_validate_result_empty_boa(self, mock_solved_atmosphere, mock_observation_bundle):
        from siac.domain import CorrectionResult
        result = CorrectionResult(
            boa=xr.Dataset(),
            boa_unc=None,
            aot=mock_solved_atmosphere.aot,
            tcwv=mock_solved_atmosphere.tcwv,
            cloud_mask=mock_observation_bundle.cloud_mask,
            metadata={"processing_time_s": 0.01},
        )
        with pytest.raises(ValidationError, match="at least one band"):
            _validate_correction_result(result)

    def test_validate_result_missing_timing(self, mock_observation_bundle, mock_solved_atmosphere):
        """Missing processing_time_s should warn, not error."""
        import warnings

        from siac.domain import CorrectionResult
        result = CorrectionResult(
            boa=mock_observation_bundle.toa,
            boa_unc=None,
            aot=mock_solved_atmosphere.aot,
            tcwv=mock_solved_atmosphere.tcwv,
            cloud_mask=mock_observation_bundle.cloud_mask,
            metadata={},
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _validate_correction_result(result)
        assert any("processing_time_s" in str(x.message) for x in w)

    def test_validate_result_optional_unc(self, mock_observation_bundle, mock_solved_atmosphere):
        """boa_unc=None should be valid."""
        from siac.domain import CorrectionResult
        result = CorrectionResult(
            boa=mock_observation_bundle.toa,
            boa_unc=None,
            aot=mock_solved_atmosphere.aot,
            tcwv=mock_solved_atmosphere.tcwv,
            cloud_mask=mock_observation_bundle.cloud_mask,
            metadata={"processing_time_s": 0.01},
        )
        _validate_correction_result(result)  # no error
