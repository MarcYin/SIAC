"""
Unit tests for SIAC core data types.
"""

import numpy as np
import pytest
import xarray as xr

from siac.catalog import LANDSAT8_OLI_CONFIG, SENTINEL2A_CONFIG, get_sensor_config
from siac.geo.resample import resample_coefficients_to_template
from siac.runtime import (
    AtmosphericState,
    GeometryAngles,
    RTCoefficients,
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
        """RAA should be |vaa - saa| wrapped to [0, 2*pi)."""
        geom = GeometryAngles(**sample_angles)

        expected_raa = np.abs(1.5 - 2.5) % (2.0 * np.pi)  # |vaa - saa| mod 2pi
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

    def test_with_updated_parameters_can_update_tco3(self, sample_state):
        """Generic parameter updates should support staged atmospheric solving."""
        new_tco3 = xr.DataArray(np.full((10, 10), 0.34), dims=["y", "x"])
        new_tco3_unc = xr.DataArray(np.full((10, 10), 0.02), dims=["y", "x"])

        updated = sample_state.with_updated_parameters(
            {"tco3": new_tco3},
            {"tco3": new_tco3_unc},
        )

        xr.testing.assert_equal(updated.get_parameter("tco3"), new_tco3)
        xr.testing.assert_equal(updated.get_uncertainty("tco3"), new_tco3_unc)
        np.testing.assert_allclose(updated.aot.values, 0.15)
        np.testing.assert_allclose(updated.tcwv.values, 2.5)

    def test_to_emulator_input_uses_shared_rt_units(self, sample_state):
        """to_emulator_input should expose the solver/LUT/emulator unit contract."""
        shape = (10, 10)
        geom = GeometryAngles(
            sza=xr.DataArray(np.full(shape, 0.5), dims=["y", "x"]),
            saa=xr.DataArray(np.full(shape, 2.5), dims=["y", "x"]),
            vza=xr.DataArray(np.full(shape, 0.1), dims=["y", "x"]),
            vaa=xr.DataArray(np.full(shape, 1.5), dims=["y", "x"]),
        )

        ds = sample_state.to_emulator_input(geom)

        np.testing.assert_allclose(ds["cos_sza"].values, np.cos(0.5))
        np.testing.assert_allclose(ds["cos_vza"].values, np.cos(0.1))
        np.testing.assert_allclose(ds["cos_raa"].values, np.cos(1.5 - 2.5))
        np.testing.assert_allclose(ds["aot"].values, 0.15)
        np.testing.assert_allclose(ds["tcwv"].values, 2.5)
        np.testing.assert_allclose(ds["tco3"].values, 0.3)
        np.testing.assert_allclose(ds["elevation"].values, 0.1)


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

    def test_simulate_toa_roundtrips_apply_correction(self, sample_coeffs):
        toa = xr.DataArray(np.full((10, 10), 0.2), dims=["y", "x"])

        simulated = sample_coeffs.simulate_toa(sample_coeffs.apply_correction(toa))

        np.testing.assert_allclose(simulated.values, toa.values, rtol=1e-5)

    def test_has_jacobian(self, sample_coeffs):
        """has_jacobian should return False when no Jacobians."""
        assert not sample_coeffs.has_jacobian

    def test_compute_boa_jacobian_matches_finite_difference(self):
        """BOA Jacobian should include both xcp and y derivative contributions."""
        dims = ("y", "x")
        coords = {"y": [0], "x": [0]}
        toa = xr.DataArray([[0.21]], dims=dims, coords=coords)
        xap = xr.DataArray([[0.92]], dims=dims, coords=coords)
        xbp = xr.DataArray([[0.04]], dims=dims, coords=coords)
        xcp = xr.DataArray([[0.35]], dims=dims, coords=coords)
        params = ["aot", "tcwv"]
        d_xap = xr.DataArray(
            np.array([[[0.025]], [[-0.015]]]),
            dims=("param", *dims),
            coords={"param": params, **coords},
        )
        d_xbp = xr.DataArray(
            np.array([[[0.006]], [[0.002]]]),
            dims=("param", *dims),
            coords={"param": params, **coords},
        )
        d_xcp = xr.DataArray(
            np.array([[[0.55]], [[-0.25]]]),
            dims=("param", *dims),
            coords={"param": params, **coords},
        )
        coeffs = RTCoefficients(
            xap=xap,
            xbp=xbp,
            xcp=xcp,
            d_xap=d_xap,
            d_xbp=d_xbp,
            d_xcp=d_xcp,
        )

        d_boa_aot, d_boa_tcwv = coeffs.compute_boa_jacobian(toa)

        def _boa_at(delta: float, param: str) -> xr.DataArray:
            shifted = RTCoefficients(
                xap=xap + delta * d_xap.sel(param=param),
                xbp=xbp + delta * d_xbp.sel(param=param),
                xcp=xcp + delta * d_xcp.sel(param=param),
            )
            return shifted.apply_correction(toa)

        eps = 1.0e-6
        expected_aot = (_boa_at(eps, "aot") - _boa_at(-eps, "aot")) / (2.0 * eps)
        expected_tcwv = (_boa_at(eps, "tcwv") - _boa_at(-eps, "tcwv")) / (2.0 * eps)
        xr.testing.assert_allclose(d_boa_aot, expected_aot, rtol=1.0e-6, atol=1.0e-9)
        xr.testing.assert_allclose(d_boa_tcwv, expected_tcwv, rtol=1.0e-6, atol=1.0e-9)

    def test_output_selection_returns_core_and_requested_extras(self, sample_coeffs):
        trans = xr.DataArray(np.full((10, 10), 0.8), dims=["y", "x"])
        coeffs = RTCoefficients(
            xap=sample_coeffs.xap,
            xbp=sample_coeffs.xbp,
            xcp=sample_coeffs.xcp,
            extras={"tgasm": trans, "sutott": trans + 0.1},
        )

        outputs = coeffs.to_output_arrays(("xap", "tgasm"))

        assert tuple(outputs) == ("xap", "tgasm")
        xr.testing.assert_equal(outputs["xap"], coeffs.xap)
        xr.testing.assert_equal(outputs["tgasm"], trans)

    def test_get_output_rejects_unknown_name(self, sample_coeffs):
        with pytest.raises(KeyError, match="Unknown RT output"):
            sample_coeffs.get_output("not_a_real_output")

    def test_extras_reject_reserved_names(self, sample_coeffs):
        with pytest.raises(ValueError, match="conflicts with a reserved field"):
            RTCoefficients(
                xap=sample_coeffs.xap,
                xbp=sample_coeffs.xbp,
                xcp=sample_coeffs.xcp,
                extras={"xap": sample_coeffs.xap},
            )

    def test_resample_coefficients_resamples_extra_outputs(self, sample_coeffs):
        extra = xr.DataArray(
            np.arange(4, dtype=np.float32).reshape(2, 2),
            dims=["y", "x"],
            coords={"y": [0.0, 1.0], "x": [0.0, 1.0]},
        )
        coeffs = RTCoefficients(
            xap=xr.DataArray(np.full((2, 2), 0.9), dims=["y", "x"], coords=extra.coords),
            xbp=xr.DataArray(np.full((2, 2), 0.1), dims=["y", "x"], coords=extra.coords),
            xcp=xr.DataArray(np.full((2, 2), 0.05), dims=["y", "x"], coords=extra.coords),
            extras={"tgasm": extra},
        )
        template = xr.DataArray(
            np.zeros((4, 4), dtype=np.float32),
            dims=["y", "x"],
            coords={"y": np.linspace(0.0, 1.0, 4), "x": np.linspace(0.0, 1.0, 4)},
        )

        resampled = resample_coefficients_to_template(coeffs, template)

        assert resampled.extras["tgasm"].shape == (4, 4)
        assert tuple(resampled.extras["tgasm"].dims) == ("y", "x")


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

    def test_default_aerosol_solver_bands(self):
        assert [band.name for band in SENTINEL2A_CONFIG.default_aerosol_solver_bands()] == [
            "B02",
            "B04",
        ]

    def test_get_sensor_config(self):
        """get_sensor_config should return correct config."""
        config = get_sensor_config("MSI", "S2A")
        assert config == SENTINEL2A_CONFIG

    def test_get_sensor_config_not_found(self):
        """get_sensor_config should raise for unknown sensor."""
        with pytest.raises(KeyError):
            get_sensor_config("UNKNOWN", "SAT")
