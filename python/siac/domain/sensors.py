"""Pure sensor domain types."""

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


__all__ = ["SensorBand", "SensorConfig"]
