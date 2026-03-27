"""
Unit tests for atmospheric correction module.
"""

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.correction.atmospheric import AtmosphericCorrector, CorrectionResult
from siac.catalog import SENTINEL2A_CONFIG
from siac.domain import SensorConfig
from siac.runtime import (
    AtmosphericState,
    GeometryAngles,
    RTCoefficients,
)


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

        # cloud_mask contract: True = cloudy (consistent with ObservationBundle)
        assert result.cloud_mask.values[10:20, 10:20].all()

    def test_correct_merges_invalid_boa_into_supplied_cloud_mask(self, sample_inputs, mock_rt_model):
        """Invalid corrected BOA pixels should be marked even when a mask is already supplied."""
        toa, geometry, atmo_state = sample_inputs
        toa["B02"].values[0, 0] = np.nan

        cloud_mask = xr.DataArray(np.zeros((50, 50), dtype=bool), dims=["y", "x"])
        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)

        result = corrector.correct(toa, geometry, atmo_state, cloud_mask=cloud_mask)

        assert bool(result.cloud_mask.values[0, 0])
        assert np.isnan(result.boa["B02"].values[0, 0])

    def test_correct_resamples_supplied_cloud_mask_before_merge(self, sample_inputs, mock_rt_model):
        """Supplied masks on shifted coordinates should be remapped before merging invalid BOA pixels."""
        toa, geometry, atmo_state = sample_inputs
        toa["B02"].values[0, 0] = np.nan

        shifted_coords = {
            "y": np.arange(50, dtype=np.float32) + 0.25,
            "x": np.arange(50, dtype=np.float32) + 0.25,
        }
        cloud_mask = xr.DataArray(np.zeros((50, 50), dtype=bool), dims=["y", "x"], coords=shifted_coords)
        cloud_mask.values[10, 10] = True

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)
        result = corrector.correct(toa, geometry, atmo_state, cloud_mask=cloud_mask)

        assert result.cloud_mask.sizes == {"y": 50, "x": 50}
        assert bool(result.cloud_mask.values[0, 0])
        assert bool(result.cloud_mask.values[10, 10])

    def test_correct_keeps_valid_sibling_bands_when_one_band_is_invalid(self, sample_inputs, mock_rt_model):
        """A one-band failure should not erase valid BOA values from other bands."""
        toa, geometry, atmo_state = sample_inputs
        toa["B02"].values[0, 0] = np.nan

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)
        result = corrector.correct(toa, geometry, atmo_state)

        assert np.isnan(result.boa["B02"].values[0, 0])
        assert np.isfinite(result.boa["B03"].values[0, 0])
        assert np.isfinite(result.boa["B04"].values[0, 0])
        assert bool(result.cloud_mask.values[0, 0])

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

    def test_correct_late_loads_missing_bands_from_toa_attrs(self, sample_inputs, mock_rt_model):
        """Correction should request missing bands on demand instead of requiring full TOA upfront."""
        toa, geometry, atmo_state = sample_inputs
        subset = xr.Dataset({"B02": toa["B02"]})
        late_calls: list[str] = []

        def _load_band(name: str) -> xr.DataArray:
            late_calls.append(name)
            if name == "B03":
                return toa["B03"]
            if name == "B04":
                return toa["B04"]
            raise KeyError(name)

        subset.attrs["_siac_toa_band_loader"] = _load_band
        sensor_config = SensorConfig(
            sensor_id="MSI",
            satellite_id="S2A",
            bands=tuple(
                SENTINEL2A_CONFIG.get_band(name)
                for name in ("B02", "B03", "B04")
            ),
        )
        corrector = AtmosphericCorrector(mock_rt_model, sensor_config)

        result = corrector.correct(subset, geometry, atmo_state)

        assert set(result.boa.data_vars) == {"B02", "B03", "B04"}
        assert late_calls == ["B03", "B04"]

    def test_correct_resamples_coefficients_to_band_grid(self):
        """Coefficient computation may stay on the atmospheric grid and upsample afterwards."""
        toa = xr.Dataset(
            {
                "B02": xr.DataArray(
                    np.full((4, 4), 0.15, dtype=np.float32),
                    dims=["y", "x"],
                    coords={"y": [3.0, 2.0, 1.0, 0.0], "x": [0.0, 1.0, 2.0, 3.0]},
                )
            }
        )
        coarse_coords = {"y": [3.0, 1.0], "x": [0.0, 2.0]}
        geometry = GeometryAngles(
            sza=xr.DataArray(np.full((2, 2), 0.5, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
            saa=xr.DataArray(np.full((2, 2), 2.5, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
            vza=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
            vaa=xr.DataArray(np.full((2, 2), 1.5, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
        )
        atmo_state = AtmosphericState(
            aot=xr.DataArray(np.full((2, 2), 0.15, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
            tcwv=xr.DataArray(np.full((2, 2), 2.5, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
            tco3=xr.DataArray(np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
            aot_unc=xr.DataArray(np.full((2, 2), 0.05, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
            tcwv_unc=xr.DataArray(np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
            tco3_unc=xr.DataArray(np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
            elevation=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"], coords=coarse_coords),
        )

        class _CoarseRTModel:
            backend_name = "coarse"

            def supports_jacobian(self) -> bool:
                return False

            def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
                _ = (sensor_id, satellite_id)
                return True

            def compute_coefficients(self, geometry, atmo_state, band, compute_jacobian=False):  # noqa: ANN001
                _ = (geometry, atmo_state, band, compute_jacobian)
                coords = coarse_coords
                return RTCoefficients(
                    xap=xr.DataArray(np.full((2, 2), 0.95, dtype=np.float32), dims=["y", "x"], coords=coords),
                    xbp=xr.DataArray(np.full((2, 2), 0.02, dtype=np.float32), dims=["y", "x"], coords=coords),
                    xcp=xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"], coords=coords),
                )

        corrector = AtmosphericCorrector(_CoarseRTModel(), SENTINEL2A_CONFIG)

        result = corrector.correct(toa, geometry, atmo_state)

        assert result.boa["B02"].shape == (4, 4)
        assert result.boa["B02"].coords["x"].identical(toa["B02"].coords["x"])
        assert result.boa["B02"].coords["y"].identical(toa["B02"].coords["y"])
        assert np.isfinite(result.boa["B02"].values).all()


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
