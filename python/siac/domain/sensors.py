"""Sensor definitions for SIAC."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SensorBand:
    """Specification for a single sensor band."""

    name: str
    center_wavelength: float
    bandwidth: float
    resolution: float
    band_index: int
    srf_wavelengths_nm: np.ndarray | None = None
    srf_response: np.ndarray | None = None
    rsrf_sensor_unit_id: str | None = None
    rsrf_representation_variant: str | None = None
    rsrf_band_id: str | None = None

    @property
    def wavelength_um(self) -> float:
        return self.center_wavelength / 1000.0

    @property
    def has_srf(self) -> bool:
        return self.srf_wavelengths_nm is not None and self.srf_response is not None

    def gaussian_response(self, wavelengths_nm: np.ndarray) -> np.ndarray:
        sigma = self.bandwidth / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        return np.exp(-0.5 * ((wavelengths_nm - self.center_wavelength) / sigma) ** 2)

    def effective_response(self, wavelengths_nm: np.ndarray) -> np.ndarray:
        if self.has_srf:
            return np.interp(
                wavelengths_nm,
                self.srf_wavelengths_nm,
                self.srf_response,
                left=0.0,
                right=0.0,
            )
        return self.gaussian_response(wavelengths_nm)


@dataclass(frozen=True)
class SensorConfig:
    """Complete sensor configuration for a satellite instrument."""

    sensor_id: str
    satellite_id: str
    bands: tuple[SensorBand, ...]
    default_ref_scale: float = 1.0 / 10000.0
    default_ref_offset: float = 0.0

    def __post_init__(self) -> None:
        if not self.bands:
            raise ValueError("SensorConfig must have at least one band.")
        names = [b.name for b in self.bands]
        if len(names) != len(set(names)):
            dupes = [n for n in names if names.count(n) > 1]
            raise ValueError(f"Duplicate band names in SensorConfig: {sorted(set(dupes))}")
        indices = [b.band_index for b in self.bands]
        if len(indices) != len(set(indices)):
            raise ValueError("Duplicate band_index values in SensorConfig.")
        for band in self.bands:
            if band.center_wavelength <= 0:
                raise ValueError(
                    f"Band {band.name!r} has non-positive center_wavelength: {band.center_wavelength}"
                )
            if band.resolution <= 0:
                raise ValueError(f"Band {band.name!r} has non-positive resolution: {band.resolution}")

    def get_band(self, name: str) -> SensorBand:
        for band in self.bands:
            if band.name == name:
                return band
        raise KeyError(f"Band {name!r} not found in sensor {self.sensor_id}")

    def get_band_by_wavelength(
        self, wavelength_nm: float, tolerance_nm: float = 20.0
    ) -> SensorBand | None:
        closest = None
        min_diff = float("inf")

        for band in self.bands:
            diff = abs(band.center_wavelength - wavelength_nm)
            if diff < min_diff and diff <= tolerance_nm:
                min_diff = diff
                closest = band

        return closest

    @property
    def band_names(self) -> tuple[str, ...]:
        return tuple(b.name for b in self.bands)

    @property
    def band_wavelengths(self) -> tuple[float, ...]:
        return tuple(b.center_wavelength for b in self.bands)

    def select_bands_in_range(
        self, wl_min_nm: float, wl_max_nm: float
    ) -> list[SensorBand]:
        return [b for b in self.bands if wl_min_nm <= b.center_wavelength <= wl_max_nm]

    def select_nearest_band(
        self, target_nm: float, tolerance_nm: float = 50.0
    ) -> SensorBand | None:
        closest = None
        min_diff = float("inf")
        for band in self.bands:
            diff = abs(band.center_wavelength - target_nm)
            if diff < min_diff and diff <= tolerance_nm:
                min_diff = diff
                closest = band
        return closest

    @property
    def vis_bands(self) -> list[SensorBand]:
        return self.select_bands_in_range(400.0, 700.0)

    @property
    def nir_bands(self) -> list[SensorBand]:
        return self.select_bands_in_range(750.0, 1000.0)

    @property
    def swir_bands(self) -> list[SensorBand]:
        return [
            b
            for b in self.bands
            if 1000.0 <= b.center_wavelength <= 2500.0
            and not (1350.0 <= b.center_wavelength <= 1420.0)
            and not (1800.0 <= b.center_wavelength <= 1950.0)
        ]


def _rsrf_band(
    name: str,
    center_wavelength: float,
    bandwidth: float,
    resolution: float,
    band_index: int,
    *,
    sensor_unit_id: str,
    representation_variant: str = "band_average",
    rsrf_band_id: str | None = None,
) -> SensorBand:
    """Construct a band with an attached RSRF repository identity."""
    return SensorBand(
        name=name,
        center_wavelength=center_wavelength,
        bandwidth=bandwidth,
        resolution=resolution,
        band_index=band_index,
        rsrf_sensor_unit_id=sensor_unit_id,
        rsrf_representation_variant=representation_variant,
        rsrf_band_id=rsrf_band_id or name,
    )


SENTINEL2A_CONFIG = SensorConfig(
    sensor_id="MSI",
    satellite_id="S2A",
    bands=(
        _rsrf_band("B01", 443.0, 20.0, 60.0, 0, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B02", 490.0, 65.0, 10.0, 1, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B03", 560.0, 35.0, 10.0, 2, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B04", 665.0, 30.0, 10.0, 3, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B05", 705.0, 15.0, 20.0, 4, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B06", 740.0, 15.0, 20.0, 5, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B07", 783.0, 20.0, 20.0, 6, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B08", 842.0, 115.0, 10.0, 7, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B8A", 865.0, 20.0, 20.0, 8, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B09", 945.0, 20.0, 60.0, 9, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B10", 1375.0, 30.0, 60.0, 10, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B11", 1610.0, 90.0, 20.0, 11, sensor_unit_id="sentinel-2a_msi"),
        _rsrf_band("B12", 2190.0, 180.0, 20.0, 12, sensor_unit_id="sentinel-2a_msi"),
    ),
    default_ref_scale=1.0 / 10000.0,
    default_ref_offset=0.0,
)

SENTINEL2B_CONFIG = SensorConfig(
    sensor_id="MSI",
    satellite_id="S2B",
    bands=(
        _rsrf_band("B01", 442.0, 20.0, 60.0, 0, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B02", 492.0, 65.0, 10.0, 1, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B03", 559.0, 35.0, 10.0, 2, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B04", 665.0, 30.0, 10.0, 3, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B05", 704.0, 15.0, 20.0, 4, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B06", 739.0, 15.0, 20.0, 5, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B07", 780.0, 20.0, 20.0, 6, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B08", 833.0, 115.0, 10.0, 7, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B8A", 864.0, 20.0, 20.0, 8, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B09", 943.0, 20.0, 60.0, 9, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B10", 1377.0, 30.0, 60.0, 10, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B11", 1610.0, 90.0, 20.0, 11, sensor_unit_id="sentinel-2b_msi"),
        _rsrf_band("B12", 2186.0, 180.0, 20.0, 12, sensor_unit_id="sentinel-2b_msi"),
    ),
    default_ref_scale=1.0 / 10000.0,
    default_ref_offset=0.0,
)

LANDSAT8_OLI_CONFIG = SensorConfig(
    sensor_id="OLI",
    satellite_id="L8",
    bands=(
        _rsrf_band("B1", 443.0, 16.0, 30.0, 0, sensor_unit_id="landsat-8_oli"),
        _rsrf_band("B2", 482.0, 60.0, 30.0, 1, sensor_unit_id="landsat-8_oli"),
        _rsrf_band("B3", 561.5, 57.0, 30.0, 2, sensor_unit_id="landsat-8_oli"),
        _rsrf_band("B4", 654.5, 37.0, 30.0, 3, sensor_unit_id="landsat-8_oli"),
        _rsrf_band("B5", 865.0, 28.0, 30.0, 4, sensor_unit_id="landsat-8_oli"),
        _rsrf_band("B6", 1608.5, 85.0, 30.0, 5, sensor_unit_id="landsat-8_oli"),
        _rsrf_band("B7", 2200.5, 187.0, 30.0, 6, sensor_unit_id="landsat-8_oli"),
    ),
    default_ref_scale=2.75e-5,
    default_ref_offset=-0.2,
)

SENTINEL2C_CONFIG = SensorConfig(
    sensor_id="MSI",
    satellite_id="S2C",
    bands=(
        _rsrf_band("B01", 443.0, 20.0, 60.0, 0, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B02", 490.0, 65.0, 10.0, 1, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B03", 560.0, 35.0, 10.0, 2, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B04", 665.0, 30.0, 10.0, 3, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B05", 705.0, 15.0, 20.0, 4, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B06", 740.0, 15.0, 20.0, 5, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B07", 783.0, 20.0, 20.0, 6, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B08", 842.0, 115.0, 10.0, 7, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B8A", 865.0, 20.0, 20.0, 8, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B09", 945.0, 20.0, 60.0, 9, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B10", 1375.0, 30.0, 60.0, 10, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B11", 1610.0, 90.0, 20.0, 11, sensor_unit_id="sentinel-2c_msi"),
        _rsrf_band("B12", 2190.0, 180.0, 20.0, 12, sensor_unit_id="sentinel-2c_msi"),
    ),
    default_ref_scale=1.0 / 10000.0,
    default_ref_offset=0.0,
)

LANDSAT9_OLI2_CONFIG = SensorConfig(
    sensor_id="OLI",
    satellite_id="L9",
    bands=(
        _rsrf_band("B1", 443.0, 16.0, 30.0, 0, sensor_unit_id="landsat-9_oli2"),
        _rsrf_band("B2", 482.0, 60.0, 30.0, 1, sensor_unit_id="landsat-9_oli2"),
        _rsrf_band("B3", 561.5, 57.0, 30.0, 2, sensor_unit_id="landsat-9_oli2"),
        _rsrf_band("B4", 654.5, 37.0, 30.0, 3, sensor_unit_id="landsat-9_oli2"),
        _rsrf_band("B5", 865.0, 28.0, 30.0, 4, sensor_unit_id="landsat-9_oli2"),
        _rsrf_band("B6", 1608.5, 85.0, 30.0, 5, sensor_unit_id="landsat-9_oli2"),
        _rsrf_band("B7", 2200.5, 187.0, 30.0, 6, sensor_unit_id="landsat-9_oli2"),
    ),
    default_ref_scale=2.75e-5,
    default_ref_offset=-0.2,
)

SENSOR_CONFIGS: dict[tuple[str, str], SensorConfig] = {
    ("MSI", "S2A"): SENTINEL2A_CONFIG,
    ("MSI", "S2B"): SENTINEL2B_CONFIG,
    ("MSI", "S2C"): SENTINEL2C_CONFIG,
    ("OLI", "L8"): LANDSAT8_OLI_CONFIG,
    ("OLI", "L9"): LANDSAT9_OLI2_CONFIG,
}


def get_sensor_config(sensor_id: str, satellite_id: str) -> SensorConfig:
    """Get sensor configuration by ID."""
    key = (sensor_id, satellite_id)
    if key not in SENSOR_CONFIGS:
        raise KeyError(
            f"Unknown sensor {sensor_id}/{satellite_id}. "
            f"Available: {list(SENSOR_CONFIGS.keys())}"
        )
    return SENSOR_CONFIGS[key]


__all__ = [
    "SensorBand",
    "SensorConfig",
    "SENTINEL2A_CONFIG",
    "SENTINEL2B_CONFIG",
    "SENTINEL2C_CONFIG",
    "LANDSAT8_OLI_CONFIG",
    "LANDSAT9_OLI2_CONFIG",
    "SENSOR_CONFIGS",
    "get_sensor_config",
]
