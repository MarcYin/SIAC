"""Parametrized invariant checks for all predefined SensorConfig instances."""

from __future__ import annotations

import pytest

from siac.catalog import (
    LANDSAT8_OLI_CONFIG,
    LANDSAT9_OLI2_CONFIG,
    SENSOR_CONFIGS,
    SENTINEL2A_CONFIG,
    SENTINEL2B_CONFIG,
    SENTINEL2C_CONFIG,
)
from siac.domain import SensorBand, SensorConfig

ALL_CONFIGS = [
    SENTINEL2A_CONFIG,
    SENTINEL2B_CONFIG,
    SENTINEL2C_CONFIG,
    LANDSAT8_OLI_CONFIG,
    LANDSAT9_OLI2_CONFIG,
]


def _config_id(cfg: SensorConfig) -> str:
    return f"{cfg.sensor_id}_{cfg.satellite_id}"


@pytest.fixture(params=ALL_CONFIGS, ids=_config_id)
def sensor_config(request: pytest.FixtureRequest) -> SensorConfig:
    return request.param


class TestSensorConfigInvariants:
    """Every predefined config must satisfy basic structural invariants."""

    def test_has_at_least_one_band(self, sensor_config: SensorConfig) -> None:
        assert len(sensor_config.bands) >= 1

    def test_band_names_unique(self, sensor_config: SensorConfig) -> None:
        names = [b.name for b in sensor_config.bands]
        assert len(names) == len(set(names)), f"Duplicate band names: {names}"

    def test_band_indices_unique(self, sensor_config: SensorConfig) -> None:
        indices = [b.band_index for b in sensor_config.bands]
        assert len(indices) == len(set(indices))

    def test_positive_wavelengths(self, sensor_config: SensorConfig) -> None:
        for b in sensor_config.bands:
            assert b.center_wavelength > 0, f"Band {b.name} has non-positive wavelength"

    def test_positive_resolutions(self, sensor_config: SensorConfig) -> None:
        for b in sensor_config.bands:
            assert b.resolution > 0, f"Band {b.name} has non-positive resolution"

    def test_positive_bandwidth(self, sensor_config: SensorConfig) -> None:
        for b in sensor_config.bands:
            assert b.bandwidth > 0, f"Band {b.name} has non-positive bandwidth"

    def test_band_wavelengths_property(self, sensor_config: SensorConfig) -> None:
        wls = sensor_config.band_wavelengths
        assert len(wls) == len(sensor_config.bands)
        assert all(isinstance(w, float) for w in wls)

    def test_band_names_property(self, sensor_config: SensorConfig) -> None:
        names = sensor_config.band_names
        assert len(names) == len(sensor_config.bands)

    def test_vis_nir_swir_classification(self, sensor_config: SensorConfig) -> None:
        """Each sensor should have at least one VIS and one NIR band."""
        assert len(sensor_config.vis_bands) >= 1
        assert len(sensor_config.nir_bands) >= 1

    def test_get_band_roundtrip(self, sensor_config: SensorConfig) -> None:
        for b in sensor_config.bands:
            assert sensor_config.get_band(b.name) is b

    def test_gaussian_response_callable(self, sensor_config: SensorConfig) -> None:
        """SensorBand.gaussian_response is available and returns correct shape."""
        import numpy as np
        wl = np.linspace(400.0, 2500.0, 100)
        for b in sensor_config.bands:
            resp = b.gaussian_response(wl)
            assert resp.shape == wl.shape
            assert resp.max() > 0

    def test_effective_response_callable(self, sensor_config: SensorConfig) -> None:
        import numpy as np
        wl = np.linspace(400.0, 2500.0, 100)
        for b in sensor_config.bands:
            resp = b.effective_response(wl)
            assert resp.shape == wl.shape


class TestSensorConfigRegistry:
    """The SENSOR_CONFIGS registry must reflect all built-in configs."""

    def test_all_predefined_configs_in_registry(self) -> None:
        for cfg in ALL_CONFIGS:
            key = (cfg.sensor_id, cfg.satellite_id)
            assert key in SENSOR_CONFIGS, f"{key} missing from SENSOR_CONFIGS"
            assert SENSOR_CONFIGS[key] is cfg

    def test_registry_size(self) -> None:
        assert len(SENSOR_CONFIGS) == len(ALL_CONFIGS)


class TestSensorConfigPostInit:
    """__post_init__ rejects bad data."""

    def test_empty_bands_raises(self) -> None:
        with pytest.raises(ValueError, match="at least one band"):
            SensorConfig(sensor_id="X", satellite_id="Y", bands=())

    def test_duplicate_band_names_raises(self) -> None:
        b1 = SensorBand("B1", 500.0, 20.0, 10.0, 0)
        b2 = SensorBand("B1", 600.0, 20.0, 10.0, 1)
        with pytest.raises(ValueError, match="Duplicate band names"):
            SensorConfig(sensor_id="X", satellite_id="Y", bands=(b1, b2))

    def test_duplicate_indices_raises(self) -> None:
        b1 = SensorBand("B1", 500.0, 20.0, 10.0, 0)
        b2 = SensorBand("B2", 600.0, 20.0, 10.0, 0)
        with pytest.raises(ValueError, match="Duplicate band_index"):
            SensorConfig(sensor_id="X", satellite_id="Y", bands=(b1, b2))

    def test_non_positive_wavelength_raises(self) -> None:
        b = SensorBand("B1", -100.0, 20.0, 10.0, 0)
        with pytest.raises(ValueError, match="non-positive center_wavelength"):
            SensorConfig(sensor_id="X", satellite_id="Y", bands=(b,))

    def test_non_positive_resolution_raises(self) -> None:
        b = SensorBand("B1", 500.0, 20.0, 0.0, 0)
        with pytest.raises(ValueError, match="non-positive resolution"):
            SensorConfig(sensor_id="X", satellite_id="Y", bands=(b,))
