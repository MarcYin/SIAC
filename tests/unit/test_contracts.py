"""
Layer 1 — Contract type construction, immutability, and properties.
"""

import dataclasses
from datetime import datetime

import numpy as np
import pytest
import xarray as xr

# ── ObservationBundle ──────────────────────────────────────────────────

class TestObservationBundle:
    def test_obs_bundle_construction(self, mock_observation_bundle):
        obs = mock_observation_bundle
        assert isinstance(obs.toa, xr.Dataset)
        assert isinstance(obs.crs, str)
        assert isinstance(obs.metadata, dict)

    def test_obs_bundle_frozen(self, mock_observation_bundle):
        with pytest.raises(dataclasses.FrozenInstanceError):
            mock_observation_bundle.crs = "EPSG:0000"

    def test_obs_bundle_metadata_time(self, mock_observation_bundle):
        assert "observation_time" in mock_observation_bundle.metadata
        assert isinstance(mock_observation_bundle.metadata["observation_time"], datetime)

    def test_obs_bundle_bounds_tuple(self, mock_observation_bundle):
        b = mock_observation_bundle.bounds
        assert len(b) == 4
        assert all(isinstance(v, float) for v in b)

    def test_obs_bundle_crs_string(self, mock_observation_bundle):
        assert mock_observation_bundle.crs.startswith("EPSG:")


# ── AtmosphericState ──────────────────────────────────────────────────

class TestAtmosphericStateContract:
    def test_atmo_state_construction(self, mock_atmospheric_state):
        atmo = mock_atmospheric_state
        for attr in ("aot", "tcwv", "tco3", "aot_unc", "tcwv_unc", "tco3_unc", "elevation"):
            assert isinstance(getattr(atmo, attr), xr.DataArray)

    def test_atmo_state_with_updated_aot_tcwv(self, mock_atmospheric_state):
        atmo = mock_atmospheric_state
        new_aot = atmo.aot * 2.0
        new_tcwv = atmo.tcwv * 1.5
        updated = atmo.with_updated_aot_tcwv(new_aot, new_tcwv)
        # New object returned
        assert updated is not atmo
        # Old is unchanged
        np.testing.assert_allclose(atmo.aot.values, 0.15)
        # tco3 preserved
        xr.testing.assert_equal(updated.tco3, atmo.tco3)

    def test_atmo_state_shape_consistency(self, mock_atmospheric_state):
        atmo = mock_atmospheric_state
        shape = atmo.aot.shape
        for attr in ("tcwv", "tco3", "aot_unc", "tcwv_unc", "tco3_unc", "elevation"):
            assert getattr(atmo, attr).shape == shape


# ── SurfacePrior ──────────────────────────────────────────────────────

class TestSurfacePriorContract:
    def test_surface_prior_construction(self, mock_surface_prior):
        sp = mock_surface_prior
        assert sp.boa.shape == sp.boa_unc.shape

    def test_surface_prior_mask_dtype(self, mock_surface_prior):
        assert mock_surface_prior.mask.dtype == bool

    def test_surface_prior_boa_unc_nonneg(self, mock_surface_prior):
        assert (mock_surface_prior.boa_unc.values >= 0).all()


# ── SolverInputBundle ─────────────────────────────────────────────────

class TestSolverInputBundle:
    def test_solver_input_bundle_construction(self, mock_solver_input_bundle):
        sib = mock_solver_input_bundle
        assert isinstance(sib.toa, xr.DataArray)
        assert isinstance(sib.atmo_prior.aot, xr.DataArray)
        assert isinstance(sib.surface_prior.boa, xr.DataArray)

    def test_solver_input_bundle_bands_subset(self, mock_solver_input_bundle):
        sib = mock_solver_input_bundle
        config_names = {b.name for b in sib.sensor_config.bands}
        solver_names = {b.name for b in sib.bands}
        assert solver_names <= config_names

    def test_solver_input_bundle_resolution_metadata(self, mock_solver_input_bundle):
        sib = mock_solver_input_bundle
        assert sib.aux_resolution_m > 0
        assert sib.aerosol_resolution_m > 0


# ── SolvedAtmosphere ─────────────────────────────────────────────────

class TestSolvedAtmosphere:
    def test_solved_atmo_construction(self, mock_solved_atmosphere):
        sa = mock_solved_atmosphere
        assert isinstance(sa.cost_final, float)
        assert isinstance(sa.n_iterations, int)
        assert isinstance(sa.converged, bool)

    def test_solved_atmo_state_has_updated_aot(self, mock_solved_atmosphere):
        sa = mock_solved_atmosphere
        xr.testing.assert_equal(sa.atmo_state.aot, sa.aot)

    def test_solved_atmo_diagnostics_types(self, mock_solved_atmosphere):
        sa = mock_solved_atmosphere
        assert sa.converged is True
        assert sa.n_iterations == 5
        assert sa.cost_final == pytest.approx(0.001)


# ── CorrectionResult ─────────────────────────────────────────────────

class TestCorrectionResult:
    def test_correction_result_construction(self):
        from siac.runtime import CorrectionDiagnostics, CorrectionResult
        shape = (4, 4)
        boa = xr.Dataset({"B02": xr.DataArray(np.ones(shape), dims=["y", "x"])})
        result = CorrectionResult(
            boa=boa,
            boa_unc=None,
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
            cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
            diagnostics=CorrectionDiagnostics(processing_time_s=1.0),
        )
        assert "B02" in result.boa.data_vars
        assert result.diagnostics.processing_time_s == pytest.approx(1.0)

    def test_correction_result_optional_unc(self):
        from siac.runtime import CorrectionResult
        shape = (4, 4)
        result = CorrectionResult(
            boa=xr.Dataset({"B02": xr.DataArray(np.ones(shape), dims=["y", "x"])}),
            boa_unc=None,
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
            cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
            metadata={},
        )
        assert result.boa_unc is None

    def test_correction_result_boa_bands(self):
        from siac.runtime import CorrectionResult
        shape = (4, 4)
        result = CorrectionResult(
            boa=xr.Dataset({
                "B02": xr.DataArray(np.ones(shape), dims=["y", "x"]),
                "B03": xr.DataArray(np.ones(shape), dims=["y", "x"]),
            }),
            boa_unc=None,
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.0), dims=["y", "x"]),
            cloud_mask=xr.DataArray(np.zeros(shape, dtype=bool), dims=["y", "x"]),
            metadata={},
        )
        assert set(result.boa.data_vars) == {"B02", "B03"}


# ── SensorConfig wavelength selection ─────────────────────────────────

class TestSensorConfigSelection:
    def test_vis_bands(self):
        from siac.catalog import SENTINEL2A_CONFIG
        vis = SENTINEL2A_CONFIG.vis_bands
        for b in vis:
            assert 400.0 <= b.center_wavelength <= 700.0

    def test_nir_bands(self):
        from siac.catalog import SENTINEL2A_CONFIG
        nir = SENTINEL2A_CONFIG.nir_bands
        for b in nir:
            assert 750.0 <= b.center_wavelength <= 1000.0

    def test_swir_bands_exclude_cirrus(self):
        from siac.catalog import SENTINEL2A_CONFIG
        swir = SENTINEL2A_CONFIG.swir_bands
        for b in swir:
            assert not (1350.0 <= b.center_wavelength <= 1420.0)

    def test_select_bands_in_range(self):
        from siac.catalog import SENTINEL2A_CONFIG
        aerosol = SENTINEL2A_CONFIG.select_bands_in_range(400.0, 520.0)
        names = [b.name for b in aerosol]
        assert "B01" in names
        assert "B02" in names

    def test_select_nearest_band(self):
        from siac.catalog import SENTINEL2A_CONFIG
        b = SENTINEL2A_CONFIG.select_nearest_band(660.0, tolerance_nm=20.0)
        assert b is not None
        assert b.name == "B04"

    def test_select_nearest_band_no_match(self):
        from siac.catalog import SENTINEL2A_CONFIG
        b = SENTINEL2A_CONFIG.select_nearest_band(3000.0, tolerance_nm=20.0)
        assert b is None

    def test_select_bands_empty_range(self):
        from siac.catalog import SENTINEL2A_CONFIG
        result = SENTINEL2A_CONFIG.select_bands_in_range(9000.0, 9999.0)
        assert result == []
