"""
Unit tests for atmospheric correction module.
"""

import numpy as np
import pytest
import xarray as xr

from siac.core.types import (
    SENTINEL2A_CONFIG,
    AtmosphericState,
    GeometryAngles,
    RTCoefficients,
)
from siac.correction.atmospheric import AtmosphericCorrector, CorrectionResult


class TestAtmosphericCorrector:
    """Tests for AtmosphericCorrector class."""

    @pytest.fixture
    def sample_inputs(self):
        """Create sample inputs for correction."""
        shape = (50, 50)

        # TOA dataset
        toa = xr.Dataset({
            "B02": xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            "B03": xr.DataArray(np.full(shape, 0.12), dims=["y", "x"]),
            "B04": xr.DataArray(np.full(shape, 0.10), dims=["y", "x"]),
        })

        # Geometry
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )

        # Atmospheric state
        atmo_state = AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )

        return toa, geometry, atmo_state

    def test_creation(self, mock_rt_model):
        """Corrector should be creatable."""
        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)

        assert corrector.sensor_config == SENTINEL2A_CONFIG

    def test_correct_basic(self, sample_inputs, mock_rt_model):
        """Basic correction should work."""
        toa, geometry, atmo_state = sample_inputs

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)

        result = corrector.correct(toa, geometry, atmo_state)

        assert isinstance(result, CorrectionResult)
        assert "B02" in result.boa.data_vars
        assert result.boa["B02"].shape == (50, 50)

    def test_correct_values(self, sample_inputs, mock_rt_model):
        """Correction should produce expected values."""
        toa, geometry, atmo_state = sample_inputs

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)

        result = corrector.correct(toa, geometry, atmo_state)

        # With mock coefficients: xap=0.95, xbp=0.02, xcp=0.1
        # For TOA=0.15: y = 0.95 * 0.15 - 0.02 = 0.1225
        # BOA = 0.1225 / (1 + 0.1 * 0.1225) = 0.1225 / 1.01225 ≈ 0.121
        expected_boa = 0.1225 / (1 + 0.1 * 0.1225)

        np.testing.assert_allclose(
            result.boa["B02"].values.mean(), expected_boa, rtol=1e-3
        )

    def test_correct_with_cloud_mask(self, sample_inputs, mock_rt_model):
        """Correction should respect cloud mask."""
        toa, geometry, atmo_state = sample_inputs

        cloud_mask = xr.DataArray(np.zeros((50, 50), dtype=bool), dims=["y", "x"])
        cloud_mask.values[10:20, 10:20] = True  # Cloud region

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)

        result = corrector.correct(toa, geometry, atmo_state, cloud_mask=cloud_mask)

        # Mask should exclude cloudy pixels
        assert (~result.cloud_mask.values[10:20, 10:20]).all()

    def test_result_has_aot_tcwv(self, sample_inputs, mock_rt_model):
        """Result should include AOT and TCWV."""
        toa, geometry, atmo_state = sample_inputs

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)

        result = corrector.correct(toa, geometry, atmo_state)

        assert result.aot is not None
        assert result.tcwv is not None
        np.testing.assert_allclose(result.aot.values, 0.15)

    def test_invalid_rt_model_raises(self):
        """Passing non-RTModelBackend should raise TypeError."""
        with pytest.raises(TypeError, match="rt_model must implement RTModelBackend"):
            AtmosphericCorrector(object(), SENTINEL2A_CONFIG)

    def test_apply_correction_consistency(self, sample_inputs, mock_rt_model):
        """Corrector should produce same result as RTCoefficients.apply_correction()."""
        toa, geometry, atmo_state = sample_inputs

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)
        result = corrector.correct(toa, geometry, atmo_state)

        # Manually compute via apply_correction
        band_spec = SENTINEL2A_CONFIG.get_band("B02")
        coeffs = mock_rt_model.compute_coefficients(geometry, atmo_state, band_spec, False)
        expected = coeffs.apply_correction(toa["B02"])

        # Filter to valid range like the corrector does
        expected = expected.where((expected > 0) & (expected < 1.5))

        np.testing.assert_allclose(
            result.boa["B02"].values, expected.values, rtol=1e-6
        )


class TestCorrectionPhysics:
    """Tests for physical correctness of correction."""

    def test_boa_less_than_toa(self):
        """BOA should typically be less than TOA (atmospheric scattering adds light)."""
        # For most land surfaces at moderate AOT
        coeffs = RTCoefficients(
            xap=xr.DataArray(np.array([[0.90]])),
            xbp=xr.DataArray(np.array([[0.05]])),
            xcp=xr.DataArray(np.array([[0.15]])),
        )

        toa = xr.DataArray(np.array([[0.20]]))
        boa = coeffs.apply_correction(toa)

        # y = 0.90 * 0.20 - 0.05 = 0.13
        # boa = 0.13 / (1 + 0.15 * 0.13) = 0.13 / 1.0195 ≈ 0.127
        assert boa.values[0, 0] < toa.values[0, 0]

    def test_boa_bounded(self):
        """BOA should be bounded between 0 and 1 for valid inputs."""
        shape = (10, 10)

        # Range of realistic coefficients
        xap = xr.DataArray(np.random.default_rng().uniform(0.8, 1.0, shape), dims=["y", "x"])
        xbp = xr.DataArray(np.random.default_rng().uniform(0.01, 0.1, shape), dims=["y", "x"])
        xcp = xr.DataArray(np.random.default_rng().uniform(0.05, 0.2, shape), dims=["y", "x"])

        coeffs = RTCoefficients(xap=xap, xbp=xbp, xcp=xcp)

        # Range of realistic TOA
        toa = xr.DataArray(np.random.default_rng().uniform(0.05, 0.4, shape), dims=["y", "x"])

        boa = coeffs.apply_correction(toa)

        assert np.all(boa.values > -0.1)  # Allow small negative from noise
        assert np.all(boa.values < 1.5)  # Allow slightly over 1 for bright targets
