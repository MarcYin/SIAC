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
        toa = xr.Dataset(
            {
                "B02": xr.DataArray(np.full(shape, 0.15), dims=["y", "x"]),
                "B03": xr.DataArray(np.full(shape, 0.12), dims=["y", "x"]),
                "B04": xr.DataArray(np.full(shape, 0.10), dims=["y", "x"]),
            }
        )

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

        np.testing.assert_allclose(result.boa["B02"].values.mean(), expected_boa, rtol=1e-3)

    def test_correct_with_cloud_mask(self, sample_inputs, mock_rt_model):
        """Correction should respect cloud mask."""
        toa, geometry, atmo_state = sample_inputs

        cloud_mask = xr.DataArray(np.zeros((50, 50), dtype=bool), dims=["y", "x"])
        cloud_mask.values[10:20, 10:20] = True  # Cloud region

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)

        result = corrector.correct(toa, geometry, atmo_state, cloud_mask=cloud_mask)

        # cloud_mask contract: True = cloudy (consistent with ObservationBundle)
        assert result.cloud_mask.values[10:20, 10:20].all()

    def test_correct_merges_invalid_boa_into_supplied_cloud_mask(
        self, sample_inputs, mock_rt_model
    ):
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
        cloud_mask = xr.DataArray(
            np.zeros((50, 50), dtype=bool), dims=["y", "x"], coords=shifted_coords
        )
        cloud_mask.values[10, 10] = True

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)
        result = corrector.correct(toa, geometry, atmo_state, cloud_mask=cloud_mask)

        assert result.cloud_mask.sizes == {"y": 50, "x": 50}
        assert bool(result.cloud_mask.values[0, 0])
        assert bool(result.cloud_mask.values[10, 10])

    def test_correct_keeps_valid_sibling_bands_when_one_band_is_invalid(
        self, sample_inputs, mock_rt_model
    ):
        """A one-band failure should not erase valid BOA values from other bands."""
        toa, geometry, atmo_state = sample_inputs
        toa["B02"].values[0, 0] = np.nan

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)
        result = corrector.correct(toa, geometry, atmo_state)

        assert np.isnan(result.boa["B02"].values[0, 0])
        assert np.isfinite(result.boa["B03"].values[0, 0])
        assert np.isfinite(result.boa["B04"].values[0, 0])
        assert bool(result.cloud_mask.values[0, 0])

    def test_correct_skips_fully_invalid_band_from_scene_mask(
        self, sample_inputs, mock_rt_model
    ):
        """A band with no valid pixel anywhere (e.g. the Sentinel-2 cirrus band
        B10, which has no surface product and corrects to NaN everywhere) carries
        no per-pixel validity signal and must not flag the whole scene invalid."""
        toa, geometry, atmo_state = sample_inputs
        toa["B02"].values[:] = np.nan  # entirely invalid, like the cirrus band

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)
        result = corrector.correct(toa, geometry, atmo_state)

        # The fully-invalid band is NaN in its own output, its valid siblings are
        # preserved, and the scene is NOT marked entirely invalid because of it.
        assert np.all(np.isnan(result.boa["B02"].values))
        assert np.isfinite(result.boa["B03"].values).all()
        assert np.isfinite(result.boa["B04"].values).all()
        assert not bool(result.cloud_mask.values.all())
        assert float(result.cloud_mask.values.mean()) < 0.5

    def test_result_has_aot_tcwv(self, sample_inputs, mock_rt_model):
        """Result should include AOT and TCWV."""
        toa, geometry, atmo_state = sample_inputs

        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG)

        result = corrector.correct(toa, geometry, atmo_state)

        assert result.aot is not None
        assert result.tcwv is not None
        np.testing.assert_allclose(result.aot.values, 0.15)

    def test_correct_overlaps_writes_with_compute(self, sample_inputs, mock_rt_model):
        """Wave 18: the per-band COG write should run inside the same executor
        task as the per-band compute so writes overlap with continued compute.

        We verify this by recording the wall-clock interval each band's writer
        callback spends and checking that the *total of those intervals*
        exceeds the *clock span* of when they're observed — i.e. multiple
        writes were in flight concurrently.
        """
        import threading
        import time

        toa, geometry, atmo_state = sample_inputs
        corrector = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG, correction_workers=4)

        write_intervals: list[tuple[float, float]] = []
        write_lock = threading.Lock()
        first_call_time: list[float] = []

        def _slow_writer(band_name: str, boa: xr.DataArray) -> xr.DataArray:
            _ = band_name
            t0 = time.perf_counter()
            if not first_call_time:
                first_call_time.append(t0)
            # 50ms simulated GIL-releasing I/O (rasterio.to_raster releases
            # the GIL inside the encoder so multiple bands can write in
            # parallel). Sleep does too.
            time.sleep(0.05)
            t1 = time.perf_counter()
            with write_lock:
                write_intervals.append((t0, t1))
            return boa

        t_start = time.perf_counter()
        corrector.correct(toa, geometry, atmo_state, boa_band_writer=_slow_writer)
        t_end = time.perf_counter()

        assert len(write_intervals) >= 3, "all 3 bands should have been written"

        # Sum of write durations vs wall-clock span: if writes truly overlap
        # with compute / each other, sum(durations) > span * 1.5.
        sum_durations = sum(end - start for start, end in write_intervals)
        clock_span = t_end - t_start
        # Heuristic: 3 bands × 50ms = 150ms of write work. If purely serial
        # the clock span ≈ sum_durations. With overlap the clock span is
        # significantly less. Use a conservative threshold of 0.8× to avoid
        # flake on slow CI.
        assert clock_span < sum_durations * 0.8, (
            f"writes did not overlap: sum_durations={sum_durations:.3f}s "
            f"vs clock_span={clock_span:.3f}s"
        )

    def test_invalid_rt_model_raises(self):
        """Passing non-RTModelBackend should raise TypeError."""
        with pytest.raises(TypeError, match="rt_model must implement RTModelBackend"):
            AtmosphericCorrector(object(), SENTINEL2A_CONFIG)

    def test_invalid_correction_workers_raises(self, mock_rt_model):
        """Correction workers must be at least one."""
        from siac.errors import ConfigurationError

        with pytest.raises(ConfigurationError, match="correction_workers must be >= 1"):
            AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG, correction_workers=0)

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

        np.testing.assert_allclose(result.boa["B02"].values, expected.values, rtol=1e-6)

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
            bands=tuple(SENTINEL2A_CONFIG.get_band(name) for name in ("B02", "B03", "B04")),
        )
        corrector = AtmosphericCorrector(mock_rt_model, sensor_config)

        result = corrector.correct(subset, geometry, atmo_state)

        assert set(result.boa.data_vars) == {"B02", "B03", "B04"}
        assert late_calls == ["B03", "B04"]

    def test_parallel_correction_matches_serial_and_writes_each_band_once(
        self, sample_inputs, mock_rt_model
    ):
        """Parallel band correction should remain numerically stable; the writer
        callback should fire exactly once per band.

        Wave 18e moved the per-band COG write inside the same executor task
        as the per-band compute, so writes happen as bands finish — the
        observable order is no longer deterministic (it depends on which
        thread finishes its compute first, which varies with system load).
        The previous test asserted strict list order ``["B02", "B03", "B04"]``
        which passed in isolation by luck of thread scheduling but failed
        intermittently under suite-level load. The new contract — every
        band is written exactly once — is what actually matters for the
        downstream COG outputs (each file is independent and named after
        the band).
        """
        import threading

        toa, geometry, atmo_state = sample_inputs

        writer_calls: list[str] = []
        writer_lock = threading.Lock()

        def _writer(band_name: str, data: xr.DataArray) -> xr.DataArray:
            with writer_lock:
                writer_calls.append(band_name)
            return data

        serial = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG, correction_workers=1)
        parallel = AtmosphericCorrector(mock_rt_model, SENTINEL2A_CONFIG, correction_workers=3)

        serial_result = serial.correct(toa, geometry, atmo_state)
        parallel_result = parallel.correct(toa, geometry, atmo_state, boa_band_writer=_writer)

        assert set(serial_result.boa.data_vars) == set(parallel_result.boa.data_vars)
        for name in serial_result.boa.data_vars:
            np.testing.assert_allclose(
                serial_result.boa[name].values, parallel_result.boa[name].values, rtol=1e-6
            )
        # Each input band must be written exactly once — order isn't constrained.
        assert sorted(writer_calls) == ["B02", "B03", "B04"]
        assert len(writer_calls) == len(set(writer_calls)), (
            f"writer was called more than once for some band: {writer_calls}"
        )

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
            sza=xr.DataArray(
                np.full((2, 2), 0.5, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
            saa=xr.DataArray(
                np.full((2, 2), 2.5, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
            vza=xr.DataArray(
                np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
            vaa=xr.DataArray(
                np.full((2, 2), 1.5, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
        )
        atmo_state = AtmosphericState(
            aot=xr.DataArray(
                np.full((2, 2), 0.15, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
            tcwv=xr.DataArray(
                np.full((2, 2), 2.5, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
            tco3=xr.DataArray(
                np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
            aot_unc=xr.DataArray(
                np.full((2, 2), 0.05, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
            tcwv_unc=xr.DataArray(
                np.full((2, 2), 0.3, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
            tco3_unc=xr.DataArray(
                np.full((2, 2), 0.01, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
            elevation=xr.DataArray(
                np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"], coords=coarse_coords
            ),
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
                    xap=xr.DataArray(
                        np.full((2, 2), 0.95, dtype=np.float32), dims=["y", "x"], coords=coords
                    ),
                    xbp=xr.DataArray(
                        np.full((2, 2), 0.02, dtype=np.float32), dims=["y", "x"], coords=coords
                    ),
                    xcp=xr.DataArray(
                        np.full((2, 2), 0.1, dtype=np.float32), dims=["y", "x"], coords=coords
                    ),
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
