"""
Unit tests for SIAC core data types.
"""

import numpy as np
import pytest
import xarray as xr

from siac.domain import (
    LANDSAT8_OLI_CONFIG,
    SENTINEL2A_CONFIG,
    AtmosphericState,
    GeometryAngles,
    RTCoefficients,
    get_sensor_config,
)


class TestGeometryAngles:
    """Tests for GeometryAngles dataclass."""

    @pytest.fixture
    def sample_angles(self):
        """Create sample angle arrays."""
        shape = (10, 10)
        return {
            "sza": xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            "saa": xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            "vza": xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            "vaa": xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        }

    def test_creation(self, sample_angles):
        """GeometryAngles should be creatable."""
        geom = GeometryAngles(**sample_angles)

        assert geom.sza.shape == (10, 10)
        np.testing.assert_allclose(geom.sza.values, 0.5)

    def test_raa_property(self, sample_angles):
        """RAA should be computed from vaa - saa."""
        geom = GeometryAngles(**sample_angles)

        expected_raa = 1.5 - 2.5  # vaa - saa
        np.testing.assert_allclose(geom.raa.values, expected_raa)

    def test_cos_properties(self, sample_angles):
        """Cosine properties should work."""
        geom = GeometryAngles(**sample_angles)

        np.testing.assert_allclose(geom.cos_sza.values, np.cos(0.5))
        np.testing.assert_allclose(geom.cos_vza.values, np.cos(0.1))

    def test_from_degrees(self):
        """from_degrees should convert to radians."""
        shape = (5, 5)
        sza_deg = xr.DataArray(np.full(shape, 30.0), dims=["y", "x"])
        saa_deg = xr.DataArray(np.full(shape, 150.0), dims=["y", "x"])
        vza_deg = xr.DataArray(np.full(shape, 5.0), dims=["y", "x"])
        vaa_deg = xr.DataArray(np.full(shape, 100.0), dims=["y", "x"])

        geom = GeometryAngles.from_degrees(sza_deg, saa_deg, vza_deg, vaa_deg)

        np.testing.assert_allclose(geom.sza.values, np.deg2rad(30.0))
        np.testing.assert_allclose(geom.vza.values, np.deg2rad(5.0))


class TestAtmosphericState:
    """Tests for AtmosphericState dataclass."""

    @pytest.fixture
    def sample_state(self):
        """Create sample atmospheric state."""
        shape = (10, 10)
        return AtmosphericState(
            aot=xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
            tcwv=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            tco3=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            aot_unc=xr.DataArray(np.full(shape, 0.05), dims=["y", "x"]),
            tcwv_unc=xr.DataArray(np.full(shape, 0.3), dims=["y", "x"]),
            tco3_unc=xr.DataArray(np.full(shape, 0.01), dims=["y", "x"]),
            elevation=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )

    def test_creation(self, sample_state):
        """AtmosphericState should be creatable."""
        assert sample_state.aot.shape == (10, 10)
        np.testing.assert_allclose(sample_state.aot.values, 0.15)

    def test_with_updated_aot_tcwv(self, sample_state):
        """with_updated_aot_tcwv should create new state."""
        new_aot = xr.DataArray(np.full((10, 10), 0.25), dims=["y", "x"])
        new_tcwv = xr.DataArray(np.full((10, 10), 3.0), dims=["y", "x"])

        updated = sample_state.with_updated_aot_tcwv(new_aot, new_tcwv)

        # Original unchanged
        np.testing.assert_allclose(sample_state.aot.values, 0.15)

        # Updated has new values
        np.testing.assert_allclose(updated.aot.values, 0.25)
        np.testing.assert_allclose(updated.tcwv.values, 3.0)

        # TCO3 preserved
        np.testing.assert_allclose(updated.tco3.values, 0.3)


class TestRTCoefficients:
    """Tests for RTCoefficients dataclass."""

    @pytest.fixture
    def sample_coeffs(self):
        """Create sample RT coefficients."""
        shape = (10, 10)
        return RTCoefficients(
            xap=xr.DataArray(np.full(shape, 0.95), dims=["y", "x"]),
            xbp=xr.DataArray(np.full(shape, 0.02), dims=["y", "x"]),
            xcp=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
        )

    def test_apply_correction(self, sample_coeffs):
        """apply_correction should compute BOA from TOA."""
        toa = xr.DataArray(np.full((10, 10), 0.2), dims=["y", "x"])

        boa = sample_coeffs.apply_correction(toa)

        # y = xap * toa - xbp = 0.95 * 0.2 - 0.02 = 0.17
        # boa = y / (1 + xcp * y) = 0.17 / (1 + 0.1 * 0.17) = 0.17 / 1.017
        expected = 0.17 / 1.017
        np.testing.assert_allclose(boa.values, expected, rtol=1e-5)

    def test_has_jacobian(self, sample_coeffs):
        """has_jacobian should return False when no Jacobians."""
        assert not sample_coeffs.has_jacobian


class TestSensorConfig:
    """Tests for sensor configuration."""

    def test_sentinel2a_config(self):
        """Sentinel-2A config should have 13 bands."""
        assert len(SENTINEL2A_CONFIG.bands) == 13
        assert SENTINEL2A_CONFIG.sensor_id == "MSI"
        assert SENTINEL2A_CONFIG.satellite_id == "S2A"

    def test_landsat8_config(self):
        """Landsat 8 config should have 7 bands."""
        assert len(LANDSAT8_OLI_CONFIG.bands) == 7
        assert LANDSAT8_OLI_CONFIG.sensor_id == "OLI"
        assert LANDSAT8_OLI_CONFIG.satellite_id == "L8"

    def test_get_band(self):
        """get_band should return correct band."""
        band = SENTINEL2A_CONFIG.get_band("B04")

        assert band.name == "B04"
        assert band.center_wavelength == 665.0

    def test_get_band_not_found(self):
        """get_band should raise KeyError for unknown band."""
        with pytest.raises(KeyError):
            SENTINEL2A_CONFIG.get_band("B99")

    def test_get_band_by_wavelength(self):
        """get_band_by_wavelength should find closest band."""
        band = SENTINEL2A_CONFIG.get_band_by_wavelength(660.0, tolerance_nm=20.0)

        assert band is not None
        assert band.name == "B04"

    def test_get_sensor_config(self):
        """get_sensor_config should return correct config."""
        config = get_sensor_config("MSI", "S2A")
        assert config == SENTINEL2A_CONFIG

    def test_get_sensor_config_not_found(self):
        """get_sensor_config should raise for unknown sensor."""
        with pytest.raises(KeyError):
            get_sensor_config("UNKNOWN", "SAT")
