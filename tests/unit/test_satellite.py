"""
Unit tests for satellite preprocessors.
"""

import threading
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import rioxarray  # noqa: F401
import xarray as xr

from siac.adapters.satellite import (
    detect_sensor,
    get_preprocessor,
    list_available_sensors,
    register_preprocessor,
)
from siac.adapters.satellite.base import (
    BaseSatellitePreprocessor,
    apply_scale_offset,
    create_valid_mask,
    degrees_to_radians,
    radians_to_degrees,
)


class TestUtilityFunctions:
    """Tests for satellite utility functions."""

    def test_degrees_to_radians(self):
        """Should convert degrees to radians."""
        degrees = xr.DataArray([0, 90, 180, 360])
        radians = degrees_to_radians(degrees)

        expected = np.array([0, np.pi / 2, np.pi, 2 * np.pi])
        np.testing.assert_allclose(radians.values, expected)

    def test_radians_to_degrees(self):
        """Should convert radians to degrees."""
        radians = xr.DataArray([0, np.pi / 2, np.pi, 2 * np.pi])
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
        from siac.adapters.satellite.sentinel2 import Sentinel2Preprocessor

        preprocessor = Sentinel2Preprocessor()

        assert "B02" in preprocessor.BAND_PATTERNS
        assert "B8A" in preprocessor.BAND_PATTERNS

    def test_parse_observation_time_from_safe_name(self):
        """Should infer sensing datetime from SAFE product names."""
        from siac.adapters.satellite.sentinel2 import Sentinel2Preprocessor

        product = "S2C_MSIL1C_20260102T024121_N0511_R089_T50QLD_20260102T035433.SAFE"
        out = Sentinel2Preprocessor._parse_observation_time_from_name(product)
        assert out == datetime(2026, 1, 2, 2, 41, 21)
        assert Sentinel2Preprocessor._parse_observation_time_from_name("invalid.SAFE") is None

    def test_load_toa_aligns_bands_to_reference_grid(self, monkeypatch, tmp_path):
        """Mixed-resolution bands should be aligned before Dataset assembly."""
        from siac.adapters.satellite import sentinel2 as s2mod

        pre = s2mod.Sentinel2Preprocessor()
        granule_dir = tmp_path / "GRANULE" / "L1C_TEST"
        img_data = granule_dir / "IMG_DATA"
        img_data.mkdir(parents=True)
        (img_data / "T50QLD_TEST_B02.jp2").touch()
        (img_data / "T50QLD_TEST_B04.jp2").touch()
        (img_data / "T50QLD_TEST_B11.jp2").touch()

        pre._granule_path = granule_dir
        pre._satellite_id = "S2A"
        monkeypatch.setattr(pre, "_resolve_paths", lambda _: None)
        monkeypatch.setattr(
            pre,
            "get_metadata",
            lambda _: {
                "quantification_value": 10000.0,
                "radiometric_offsets": {},
                "observation_time": datetime(2026, 1, 2, 2, 41, 21),
            },
        )

        def _da(values: np.ndarray, x: np.ndarray, y: np.ndarray) -> xr.DataArray:
            da = xr.DataArray(values.astype(np.float32), dims=("y", "x"), coords={"x": x, "y": y})
            return da.rio.write_crs("EPSG:32650")

        def fake_read_raster(path):
            name = str(path)
            if "B11" in name:
                return _da(
                    np.array([[500.0, 600.0], [700.0, 800.0]]),
                    np.array([0.0, 30.0]),
                    np.array([30.0, 0.0]),
                )
            if "B02" in name:
                return _da(
                    np.full((4, 4), 1000.0),
                    np.array([0.0, 10.0, 20.0, 30.0]),
                    np.array([30.0, 20.0, 10.0, 0.0]),
                )
            return _da(
                np.full((4, 4), 2000.0),
                np.array([0.0, 10.0, 20.0, 30.0]),
                np.array([30.0, 20.0, 10.0, 0.0]),
            )

        def fake_reproject_match(source, target, resampling="bilinear"):
            del resampling
            out = source.interp(y=target.coords["y"], x=target.coords["x"], method="linear")
            return out.rio.write_crs(target.rio.crs)

        monkeypatch.setattr(s2mod, "read_raster", fake_read_raster)
        monkeypatch.setattr(s2mod, "reproject_match", fake_reproject_match)

        ds = pre.load_toa(tmp_path)

        assert ds["B02"].shape == (4, 4)
        assert ds["B04"].shape == (4, 4)
        assert ds["B11"].shape == (4, 4)
        assert np.isfinite(ds["B11"].values).all()

    def test_load_toa_uses_windowed_reads_when_subset_is_active(self, monkeypatch, tmp_path):
        """Explicit AOI subsets should read only the requested raster window."""
        from siac.adapters.satellite import sentinel2 as s2mod

        pre = s2mod.Sentinel2Preprocessor()
        granule_dir = tmp_path / "GRANULE" / "L1C_TEST"
        img_data = granule_dir / "IMG_DATA"
        img_data.mkdir(parents=True)
        (img_data / "T50QLD_TEST_B02.jp2").touch()
        (img_data / "T50QLD_TEST_B04.jp2").touch()

        pre._granule_path = granule_dir
        pre._satellite_id = "S2A"
        monkeypatch.setattr(pre, "_resolve_paths", lambda _: None)
        monkeypatch.setattr(
            pre,
            "get_metadata",
            lambda _: {
                "quantification_value": 10000.0,
                "radiometric_offsets": {},
                "observation_time": datetime(2026, 1, 2, 2, 41, 21),
            },
        )

        window_calls: list[tuple[str, tuple[float, float, float, float], str]] = []
        full_calls: list[str] = []

        def _da(values: np.ndarray) -> xr.DataArray:
            da = xr.DataArray(
                values.astype(np.float32),
                dims=("y", "x"),
                coords={"x": np.array([0.0, 10.0]), "y": np.array([10.0, 0.0])},
            )
            return da.rio.write_crs("EPSG:32650")

        def fake_read_raster_window(path, *, bounds, bounds_crs):  # noqa: ANN001
            window_calls.append((str(path), bounds, bounds_crs))
            value = 1000.0 if "B02" in str(path) else 2000.0
            return _da(np.full((2, 2), value, dtype=np.float32))

        def fake_read_raster(path):  # noqa: ANN001
            full_calls.append(str(path))
            return _da(np.full((2, 2), 9999.0, dtype=np.float32))

        monkeypatch.setattr(s2mod, "read_raster_window", fake_read_raster_window)
        monkeypatch.setattr(s2mod, "read_raster", fake_read_raster)

        pre.set_spatial_subset((100.0, 200.0, 300.0, 400.0), crs="EPSG:32650")
        ds = pre.load_toa(tmp_path)
        pre.clear_spatial_subset()

        assert full_calls == []
        assert len(window_calls) == 2
        assert all(call[1] == (100.0, 200.0, 300.0, 400.0) for call in window_calls)
        assert all(call[2] == "EPSG:32650" for call in window_calls)
        assert ds["B02"].shape == (2, 2)
        assert ds["B04"].shape == (2, 2)

    def test_load_toa_supports_selective_preload_and_late_band_loader(self, monkeypatch, tmp_path):
        """Selective TOA loads should defer missing bands until explicitly requested."""
        from siac.adapters.satellite import sentinel2 as s2mod

        pre = s2mod.Sentinel2Preprocessor()
        granule_dir = tmp_path / "GRANULE" / "L1C_TEST"
        img_data = granule_dir / "IMG_DATA"
        img_data.mkdir(parents=True)
        for band_name in ("B02", "B04", "B11"):
            (img_data / f"T50QLD_TEST_{band_name}.jp2").touch()

        pre._granule_path = granule_dir
        pre._satellite_id = "S2A"
        monkeypatch.setattr(pre, "_resolve_paths", lambda _: None)
        monkeypatch.setattr(
            pre,
            "get_metadata",
            lambda _: {
                "quantification_value": 10000.0,
                "radiometric_offsets": {},
                "observation_time": datetime(2026, 1, 2, 2, 41, 21),
            },
        )

        read_calls: list[str] = []

        def _da(values: np.ndarray, x: np.ndarray, y: np.ndarray) -> xr.DataArray:
            da = xr.DataArray(values.astype(np.float32), dims=("y", "x"), coords={"x": x, "y": y})
            return da.rio.write_crs("EPSG:32650")

        def fake_read_raster(path):  # noqa: ANN001
            read_calls.append(Path(path).name)
            name = str(path)
            if "B11" in name:
                return _da(
                    np.array([[500.0, 600.0], [700.0, 800.0]]),
                    np.array([0.0, 30.0]),
                    np.array([30.0, 0.0]),
                )
            if "B02" in name:
                return _da(
                    np.full((4, 4), 1000.0),
                    np.array([0.0, 10.0, 20.0, 30.0]),
                    np.array([30.0, 20.0, 10.0, 0.0]),
                )
            return _da(
                np.full((4, 4), 2000.0),
                np.array([0.0, 10.0, 20.0, 30.0]),
                np.array([30.0, 20.0, 10.0, 0.0]),
            )

        def fake_reproject_match(source, target, resampling="bilinear"):  # noqa: ANN001
            del resampling
            out = source.interp(y=target.coords["y"], x=target.coords["x"], method="linear")
            return out.rio.write_crs(target.rio.crs)

        monkeypatch.setattr(s2mod, "read_raster", fake_read_raster)
        monkeypatch.setattr(s2mod, "reproject_match", fake_reproject_match)

        ds = pre.load_toa(tmp_path, band_names=("B02", "B04"))

        assert set(ds.data_vars) == {"B02", "B04"}
        assert read_calls == ["T50QLD_TEST_B04.jp2", "T50QLD_TEST_B02.jp2"]

        loader = ds.attrs["_siac_toa_band_loader"]
        late_b11 = loader("B11")

        assert late_b11.shape == (4, 4)
        assert np.isfinite(late_b11.values).all()
        assert read_calls == [
            "T50QLD_TEST_B04.jp2",
            "T50QLD_TEST_B02.jp2",
            "T50QLD_TEST_B11.jp2",
        ]

    def test_load_toa_uses_finest_available_reference_grid_for_late_loaded_bands(
        self, monkeypatch, tmp_path
    ):
        """A coarse-first subset should still keep the loader anchored to the finest SAFE grid."""
        from siac.adapters.satellite import sentinel2 as s2mod

        pre = s2mod.Sentinel2Preprocessor()
        granule_dir = tmp_path / "GRANULE" / "L1C_TEST"
        img_data = granule_dir / "IMG_DATA"
        img_data.mkdir(parents=True)
        for band_name in ("B02", "B04", "B11"):
            (img_data / f"T50QLD_TEST_{band_name}.jp2").touch()

        pre._granule_path = granule_dir
        pre._satellite_id = "S2A"
        monkeypatch.setattr(pre, "_resolve_paths", lambda _: None)
        monkeypatch.setattr(
            pre,
            "get_metadata",
            lambda _: {
                "quantification_value": 10000.0,
                "radiometric_offsets": {},
                "observation_time": datetime(2026, 1, 2, 2, 41, 21),
            },
        )

        read_calls: list[str] = []

        def _da(values: np.ndarray, x: np.ndarray, y: np.ndarray) -> xr.DataArray:
            da = xr.DataArray(values.astype(np.float32), dims=("y", "x"), coords={"x": x, "y": y})
            return da.rio.write_crs("EPSG:32650")

        def fake_read_raster(path):  # noqa: ANN001
            read_calls.append(Path(path).name)
            name = str(path)
            if "B11" in name:
                return _da(
                    np.array([[500.0, 600.0], [700.0, 800.0]]),
                    np.array([0.0, 30.0]),
                    np.array([30.0, 0.0]),
                )
            if "B02" in name:
                return _da(
                    np.full((4, 4), 1000.0),
                    np.array([0.0, 10.0, 20.0, 30.0]),
                    np.array([30.0, 20.0, 10.0, 0.0]),
                )
            return _da(
                np.full((4, 4), 2000.0),
                np.array([0.0, 10.0, 20.0, 30.0]),
                np.array([30.0, 20.0, 10.0, 0.0]),
            )

        def fake_reproject_match(source, target, resampling="bilinear"):  # noqa: ANN001
            del resampling
            out = source.interp(y=target.coords["y"], x=target.coords["x"], method="linear")
            return out.rio.write_crs(target.rio.crs)

        monkeypatch.setattr(s2mod, "read_raster", fake_read_raster)
        monkeypatch.setattr(s2mod, "reproject_match", fake_reproject_match)

        ds = pre.load_toa(tmp_path, band_names=("B11",))

        assert set(ds.data_vars) == {"B11"}
        assert ds["B11"].shape == (4, 4)
        assert read_calls[:2] == ["T50QLD_TEST_B04.jp2", "T50QLD_TEST_B11.jp2"]

        loader = ds.attrs["_siac_toa_band_loader"]
        late_b02 = loader("B02")

        assert late_b02.shape == (4, 4)
        assert np.isfinite(late_b02.values).all()
        assert read_calls == [
            "T50QLD_TEST_B04.jp2",
            "T50QLD_TEST_B11.jp2",
            "T50QLD_TEST_B02.jp2",
        ]

    def test_extract_cloud_mask_auto_requires_green_band(self, tmp_path):
        pre = get_preprocessor("s2")
        toa = xr.Dataset(
            {
                "B04": xr.DataArray(np.full((4, 4), 0.2, dtype=np.float32), dims=["y", "x"]),
                "B08": xr.DataArray(np.full((4, 4), 0.3, dtype=np.float32), dims=["y", "x"]),
            }
        )

        with pytest.raises(ValueError, match="Could not find any green band"):
            pre.extract_cloud_mask(tmp_path / "missing-green.SAFE", toa=toa)

    def test_discover_band_files_distinguishes_b08_and_b8a(self, tmp_path):
        from siac.adapters.satellite.sentinel2 import Sentinel2Preprocessor

        img_data = tmp_path / "IMG_DATA"
        img_data.mkdir()
        (img_data / "T50QLD_TEST_B08.jp2").touch()
        (img_data / "T50QLD_TEST_B8A.jp2").touch()
        (img_data / "T50QLD_TEST_B11.jp2").touch()

        discovered = Sentinel2Preprocessor._discover_band_files(img_data)

        assert discovered["B08"].name.endswith("B08.jp2")
        assert discovered["B8A"].name.endswith("B8A.jp2")
        assert discovered["B11"].name.endswith("B11.jp2")


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
                from siac.catalog import SENTINEL2A_CONFIG

                return SENTINEL2A_CONFIG

            def load_toa(self, input_path):
                self.calls.append("load_toa")
                return xr.Dataset({"B02": xr.DataArray(np.ones((10, 10)))})

            def extract_geometry(self, input_path):
                self.calls.append("extract_geometry")
                from siac.runtime import GeometryAngles

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

    def test_preprocess_parallelizes_geometry_and_cloud_and_reuses_toa(self):
        """Preprocess should run geometry/cloud in parallel and avoid TOA re-reads."""

        geom_started = threading.Event()
        cloud_started = threading.Event()

        class TestPreprocessor(BaseSatellitePreprocessor):
            def __init__(self):
                super().__init__()
                self.load_toa_calls = 0
                self._last_toa = None

            @property
            def sensor_config(self):
                from siac.catalog import SENTINEL2A_CONFIG

                return SENTINEL2A_CONFIG

            def load_toa(self, input_path):
                del input_path
                self.load_toa_calls += 1
                self._last_toa = xr.Dataset(
                    {
                        "B02": xr.DataArray(
                            np.ones((6, 6), dtype=np.float32),
                            dims=["y", "x"],
                        )
                    }
                )
                return self._last_toa

            def extract_geometry(self, input_path):
                del input_path
                from siac.runtime import GeometryAngles

                geom_started.set()
                assert cloud_started.wait(timeout=1.0)
                arr = xr.DataArray(np.ones((6, 6)), dims=["y", "x"])
                return GeometryAngles(sza=arr, saa=arr, vza=arr, vaa=arr)

            def extract_cloud_mask(self, input_path, toa=None):
                del input_path
                cloud_started.set()
                assert geom_started.wait(timeout=1.0)
                assert toa is self._last_toa
                return xr.DataArray(np.zeros((6, 6), dtype=bool), dims=["y", "x"])

            def get_metadata(self, input_path):
                del input_path
                return {"sensor": "TEST"}

        preprocessor = TestPreprocessor()
        result = preprocessor.preprocess("/fake/path")

        assert preprocessor.load_toa_calls == 1
        assert geom_started.is_set()
        assert cloud_started.is_set()
        assert "cloud_mask" in result

    def test_preprocess_passes_configured_toa_band_subset_when_supported(self):
        """Preprocess should forward configured preload bands to compatible loaders."""

        class TestPreprocessor(BaseSatellitePreprocessor):
            def __init__(self):
                super().__init__({"preload_toa_bands": ("B02", "B03")})
                self.seen_band_names = None

            @property
            def sensor_config(self):
                from siac.catalog import SENTINEL2A_CONFIG

                return SENTINEL2A_CONFIG

            def load_toa(self, input_path, *, band_names=None):
                del input_path
                self.seen_band_names = band_names
                return xr.Dataset(
                    {"B02": xr.DataArray(np.ones((4, 4), dtype=np.float32), dims=["y", "x"])}
                )

            def extract_geometry(self, input_path):
                del input_path
                from siac.runtime import GeometryAngles

                arr = xr.DataArray(np.ones((4, 4), dtype=np.float32), dims=["y", "x"])
                return GeometryAngles(sza=arr, saa=arr, vza=arr, vaa=arr)

            def extract_cloud_mask(self, input_path, toa=None):
                del input_path, toa
                return xr.DataArray(np.zeros((4, 4), dtype=bool), dims=["y", "x"])

            def get_metadata(self, input_path):
                del input_path
                return {"observation_time": datetime(2024, 1, 1)}

        preprocessor = TestPreprocessor()
        preprocessor.preprocess("/fake/path")

        assert preprocessor.seen_band_names == ("B02", "B03")
