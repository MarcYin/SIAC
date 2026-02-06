"""
Unit tests for satellite preprocessors.
"""

import numpy as np
import pytest
import xarray as xr
from pathlib import Path

from siac.satellite import (
    get_preprocessor,
    detect_sensor,
    register_preprocessor,
    list_available_sensors,
)
from siac.satellite.base import (
    BaseSatellitePreprocessor,
    degrees_to_radians,
    radians_to_degrees,
    apply_scale_offset,
    create_valid_mask,
)


class TestUtilityFunctions:
    """Tests for satellite utility functions."""

    def test_degrees_to_radians(self):
        """Should convert degrees to radians."""
        degrees = xr.DataArray([0, 90, 180, 360])
        radians = degrees_to_radians(degrees)

        expected = np.array([0, np.pi/2, np.pi, 2*np.pi])
        np.testing.assert_allclose(radians.values, expected)

    def test_radians_to_degrees(self):
        """Should convert radians to degrees."""
        radians = xr.DataArray([0, np.pi/2, np.pi, 2*np.pi])
        degrees = radians_to_degrees(radians)

        expected = np.array([0, 90, 180, 360])
        np.testing.assert_allclose(degrees.values, expected)

    def test_apply_scale_offset(self):
        """Should apply scale and offset correctly."""
        data = xr.DataArray([1000, 2000, 3000])

        result = apply_scale_offset(data, scale=0.0001, offset=-0.1)

        expected = np.array([0.0, 0.1, 0.2])
        np.testing.assert_allclose(result.values, expected)

    def test_create_valid_mask(self):
        """Should create correct validity mask."""
        data = xr.DataArray([0.1, 0.5, 1.2, -0.1, np.nan, 1.6])

        mask = create_valid_mask(data, min_val=0.0, max_val=1.5)

        expected = np.array([True, True, True, False, False, False])
        np.testing.assert_array_equal(mask.values, expected)


class TestSensorRegistry:
    """Tests for sensor registry functions."""

    def test_list_available_sensors(self):
        """Should list registered sensors."""
        sensors = list_available_sensors()

        assert isinstance(sensors, list)
        assert "s2" in sensors

    def test_get_preprocessor_s2(self):
        """Should get Sentinel-2 preprocessor."""
        preprocessor = get_preprocessor("s2")

        assert preprocessor is not None
        assert preprocessor.sensor_config.sensor_id == "MSI"

    def test_get_preprocessor_unknown(self):
        """Should raise for unknown sensor."""
        with pytest.raises(KeyError):
            get_preprocessor("unknown_sensor")

    def test_register_preprocessor(self):
        """Should register custom preprocessor."""

        @register_preprocessor("test_sensor")
        class TestPreprocessor(BaseSatellitePreprocessor):
            @property
            def sensor_config(self):
                return None

            def load_toa(self, input_path):
                return None

            def extract_geometry(self, input_path):
                return None

            def extract_cloud_mask(self, input_path):
                return None

            def get_metadata(self, input_path):
                return {}

        sensors = list_available_sensors()
        assert "test_sensor" in sensors


class TestSentinel2Preprocessor:
    """Tests for Sentinel-2 preprocessor."""

    def test_sensor_config(self):
        """Should return correct sensor config."""
        preprocessor = get_preprocessor("s2")

        config = preprocessor.sensor_config
        assert config.sensor_id == "MSI"
        assert config.satellite_id in ["S2A", "S2B"]
        assert len(config.bands) == 13

    def test_band_patterns(self):
        """Should have correct band file patterns."""
        from siac.satellite.sentinel2 import Sentinel2Preprocessor

        preprocessor = Sentinel2Preprocessor()

        assert "B02" in preprocessor.BAND_PATTERNS
        assert "B8A" in preprocessor.BAND_PATTERNS


class TestDetectSensor:
    """Tests for sensor detection."""

    def test_detect_safe_suffix(self, tmp_path):
        """Should detect S2 from .SAFE suffix."""
        safe_path = tmp_path / "S2A_MSIL1C_20230715.SAFE"
        safe_path.mkdir()

        sensor = detect_sensor(safe_path)
        assert sensor == "s2"

    def test_detect_mtd_file(self, tmp_path):
        """Should detect S2 from MTD file."""
        data_path = tmp_path / "data"
        data_path.mkdir()
        (data_path / "MTD_MSIL1C.xml").touch()

        sensor = detect_sensor(data_path)
        assert sensor == "s2"

    def test_detect_unknown(self, tmp_path):
        """Should raise for unknown format."""
        data_path = tmp_path / "unknown_data"
        data_path.mkdir()

        with pytest.raises(ValueError):
            detect_sensor(data_path)


class TestPreprocessorBase:
    """Tests for BaseSatellitePreprocessor."""

    def test_preprocess_flow(self):
        """Preprocess should call all abstract methods."""

        class TestPreprocessor(BaseSatellitePreprocessor):
            def __init__(self):
                super().__init__()
                self.calls = []

            @property
            def sensor_config(self):
                from siac.core.types import SENTINEL2A_CONFIG
                return SENTINEL2A_CONFIG

            def load_toa(self, input_path):
                self.calls.append("load_toa")
                return xr.Dataset({"B02": xr.DataArray(np.ones((10, 10)))})

            def extract_geometry(self, input_path):
                self.calls.append("extract_geometry")
                from siac.core.types import GeometryAngles
                arr = xr.DataArray(np.ones((10, 10)), dims=["y", "x"])
                return GeometryAngles(sza=arr, saa=arr, vza=arr, vaa=arr)

            def extract_cloud_mask(self, input_path):
                self.calls.append("extract_cloud_mask")
                return xr.DataArray(np.zeros((10, 10), dtype=bool), dims=["y", "x"])

            def get_metadata(self, input_path):
                self.calls.append("get_metadata")
                return {"sensor": "TEST"}

        preprocessor = TestPreprocessor()
        result = preprocessor.preprocess("/fake/path")

        assert "load_toa" in preprocessor.calls
        assert "extract_geometry" in preprocessor.calls
        assert "extract_cloud_mask" in preprocessor.calls
        assert "get_metadata" in preprocessor.calls
        assert "toa" in result
        assert "geometry" in result
