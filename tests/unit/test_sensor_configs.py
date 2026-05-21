"""Parametrized invariant checks for all predefined SensorConfig instances."""

from __future__ import annotations

import numpy as np
import pytest

from siac.adapters.rsrf import band_convolution_weights
from siac.catalog import (
    LANDSAT8_OLI_CONFIG,
    LANDSAT9_OLI2_CONFIG,
    SENSOR_CONFIGS,
    SENTINEL2A_CONFIG,
    SENTINEL2B_CONFIG,
    SENTINEL2C_CONFIG,
)
from siac.catalog.sensors import register
from siac.catalog.sensors.registry import get_sensor_config
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

    def test_band_convolution_weights_callable(self, sensor_config: SensorConfig) -> None:
        wl = np.linspace(400.0, 2500.0, 100)
        for b in sensor_config.bands:
            weights = band_convolution_weights(b, wl)
            assert weights.shape == wl.shape
            assert np.all(np.isfinite(weights))
            assert float(np.sum(weights, dtype=np.float64)) == pytest.approx(1.0, rel=1e-5)


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

    def test_duplicate_band_names_lists_each_dupe_once(self) -> None:
        # The O(N) Counter-based detector should report a single, sorted list
        # of duplicate names rather than re-emitting every occurrence.
        b1 = SensorBand("B1", 500.0, 20.0, 10.0, 0)
        b2 = SensorBand("B1", 600.0, 20.0, 10.0, 1)
        b3 = SensorBand("B1", 700.0, 20.0, 10.0, 2)
        with pytest.raises(ValueError) as ei:
            SensorConfig(sensor_id="X", satellite_id="Y", bands=(b1, b2, b3))
        msg = str(ei.value)
        # 'B1' appears once in the message even though it occurs three times.
        assert msg.count("'B1'") == 1


class TestAerosolSolverBandConfig:
    """``aerosol_solver_band_names`` drives ``default_aerosol_solver_bands``."""

    def test_msi_uses_catalog_field(self) -> None:
        # S2 catalog entries declare ("B02", "B04"); the dataclass field
        # replaces the previously hard-coded 'MSI' branch.
        assert SENTINEL2A_CONFIG.aerosol_solver_band_names == ("B02", "B04")
        assert [b.name for b in SENTINEL2A_CONFIG.default_aerosol_solver_bands()] == [
            "B02",
            "B04",
        ]

    def test_oli_uses_catalog_field(self) -> None:
        assert LANDSAT8_OLI_CONFIG.aerosol_solver_band_names == ("B1", "B2")
        assert [b.name for b in LANDSAT8_OLI_CONFIG.default_aerosol_solver_bands()] == [
            "B1",
            "B2",
        ]

    def test_generic_fallback_when_field_unset(self) -> None:
        # When no catalog override is set, the generic 400-520 nm fallback
        # (aerosol-sensitive blue) is used — matches the pre-refactor
        # non-MSI behaviour exactly.
        bands = (
            SensorBand("B1", 450.0, 20.0, 10.0, 0),
            SensorBand("B2", 500.0, 20.0, 10.0, 1),
            SensorBand("B3", 600.0, 20.0, 10.0, 2),
            SensorBand("B4", 1500.0, 20.0, 10.0, 3),
        )
        cfg = SensorConfig(sensor_id="X", satellite_id="Y", bands=bands)
        assert cfg.aerosol_solver_band_names == ()
        # B3 (600 nm) and B4 (1500 nm) are outside the 400-520 nm window.
        assert [b.name for b in cfg.default_aerosol_solver_bands()] == ["B1", "B2"]

    def test_unknown_aerosol_band_falls_back_to_generic(self) -> None:
        # If an aerosol_solver_band_names entry is not present in bands,
        # the helper falls through to the generic blue window rather than
        # emitting an empty list.
        bands = (
            SensorBand("B1", 450.0, 20.0, 10.0, 0),
            SensorBand("B2", 600.0, 20.0, 10.0, 1),
        )
        cfg = SensorConfig(
            sensor_id="X",
            satellite_id="Y",
            bands=bands,
            aerosol_solver_band_names=("BX",),
        )
        # Generic 400-520 fallback returns B1 only.
        assert [b.name for b in cfg.default_aerosol_solver_bands()] == ["B1"]


class TestSelectNearestBandConsolidation:
    """``select_nearest_band`` shares its implementation with
    ``get_band_by_wavelength``; only the default tolerance differs."""

    def test_same_match_when_tolerance_explicit(self) -> None:
        # With the same explicit tolerance, both helpers must return the
        # same band — they share the same loop.
        b1 = SENTINEL2A_CONFIG.get_band_by_wavelength(660.0, tolerance_nm=20.0)
        b2 = SENTINEL2A_CONFIG.select_nearest_band(660.0, tolerance_nm=20.0)
        assert b1 is b2

    def test_select_nearest_band_default_tolerance_wider(self) -> None:
        # 50 nm default catches B04 (665 nm) at 620 nm; 20 nm default would not.
        assert SENTINEL2A_CONFIG.get_band_by_wavelength(620.0) is None
        assert SENTINEL2A_CONFIG.select_nearest_band(620.0) is not None


class TestRegisterApi:
    """The ``register`` function lets external packages add sensors."""

    def _make_cfg(self, sensor_id: str, satellite_id: str) -> SensorConfig:
        return SensorConfig(
            sensor_id=sensor_id,
            satellite_id=satellite_id,
            bands=(SensorBand("B1", 500.0, 20.0, 10.0, 0),),
        )

    def test_register_adds_entry(self) -> None:
        try:
            cfg = self._make_cfg("FAKE", "F1")
            register("FAKE", "F1", cfg)
            assert SENSOR_CONFIGS[("FAKE", "F1")] is cfg
            assert get_sensor_config("FAKE", "F1") is cfg
        finally:
            SENSOR_CONFIGS.pop(("FAKE", "F1"), None)

    def test_register_duplicate_without_overwrite_raises(self) -> None:
        try:
            cfg = self._make_cfg("FAKE", "F2")
            register("FAKE", "F2", cfg)
            with pytest.raises(KeyError, match="already registered"):
                register("FAKE", "F2", cfg)
        finally:
            SENSOR_CONFIGS.pop(("FAKE", "F2"), None)

    def test_register_overwrite_replaces(self) -> None:
        try:
            cfg1 = self._make_cfg("FAKE", "F3")
            cfg2 = self._make_cfg("FAKE", "F3")
            register("FAKE", "F3", cfg1)
            register("FAKE", "F3", cfg2, overwrite=True)
            assert SENSOR_CONFIGS[("FAKE", "F3")] is cfg2
        finally:
            SENSOR_CONFIGS.pop(("FAKE", "F3"), None)

    def test_register_mismatched_keys_raises(self) -> None:
        cfg = self._make_cfg("FAKE", "F4")
        with pytest.raises(ValueError, match="must match keys"):
            register("OTHER", "F4", cfg)


class TestSentinel2ConfigSharedLayout:
    """The S2A/S2B/S2C de-dup helper preserves the shared band layout."""

    @pytest.mark.parametrize("cfg", [SENTINEL2A_CONFIG, SENTINEL2B_CONFIG, SENTINEL2C_CONFIG])
    def test_band_count(self, cfg: SensorConfig) -> None:
        assert len(cfg.bands) == 13

    @pytest.mark.parametrize("cfg", [SENTINEL2A_CONFIG, SENTINEL2B_CONFIG, SENTINEL2C_CONFIG])
    def test_band_names_and_indices_match(self, cfg: SensorConfig) -> None:
        # All three S2 satellites share band names and indices verbatim;
        # only band centres / FWHMs differ between satellites.
        assert cfg.band_names == (
            "B01",
            "B02",
            "B03",
            "B04",
            "B05",
            "B06",
            "B07",
            "B08",
            "B8A",
            "B09",
            "B10",
            "B11",
            "B12",
        )
        assert tuple(b.band_index for b in cfg.bands) == tuple(range(13))

    def test_s2b_centre_overrides(self) -> None:
        # Spot-check a couple of S2B-specific overrides; the rest of the
        # band layout should be identical to S2A.
        assert SENTINEL2B_CONFIG.get_band("B01").center_wavelength == 442.0
        assert SENTINEL2B_CONFIG.get_band("B12").center_wavelength == 2186.0

    def test_s2a_and_s2c_share_centres(self) -> None:
        # S2C reuses the S2A band-centre table verbatim in this catalog
        # snapshot; the helper's no-override path must preserve that.
        for band_a, band_c in zip(SENTINEL2A_CONFIG.bands, SENTINEL2C_CONFIG.bands):
            assert band_a.center_wavelength == band_c.center_wavelength
            assert band_a.bandwidth == band_c.bandwidth
