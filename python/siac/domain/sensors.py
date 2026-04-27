"""Pure sensor domain types."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, init=False)
class SensorBand:
    """Specification for a single sensor band."""

    name: str
    center_wavelength: float
    bandwidth: float
    resolution: float
    band_index: int
    rsrf_wavelengths_nm: np.ndarray | None = None
    rsrf_response: np.ndarray | None = None
    rsrf_sensor_unit_id: str | None = None
    rsrf_representation_variant: str | None = None
    rsrf_band_id: str | None = None

    def __init__(
        self,
        name: str,
        center_wavelength: float,
        bandwidth: float,
        resolution: float,
        band_index: int,
        rsrf_wavelengths_nm: np.ndarray | None = None,
        rsrf_response: np.ndarray | None = None,
        rsrf_sensor_unit_id: str | None = None,
        rsrf_representation_variant: str | None = None,
        rsrf_band_id: str | None = None,
        *,
        srf_wavelengths_nm: np.ndarray | None = None,
        srf_response: np.ndarray | None = None,
    ) -> None:
        if rsrf_wavelengths_nm is None:
            rsrf_wavelengths_nm = srf_wavelengths_nm
        elif srf_wavelengths_nm is not None and not np.array_equal(
            np.asarray(rsrf_wavelengths_nm),
            np.asarray(srf_wavelengths_nm),
        ):
            raise ValueError("Specify either rsrf_wavelengths_nm or srf_wavelengths_nm, not both.")

        if rsrf_response is None:
            rsrf_response = srf_response
        elif srf_response is not None and not np.array_equal(
            np.asarray(rsrf_response),
            np.asarray(srf_response),
        ):
            raise ValueError("Specify either rsrf_response or srf_response, not both.")

        object.__setattr__(self, "name", name)
        object.__setattr__(self, "center_wavelength", center_wavelength)
        object.__setattr__(self, "bandwidth", bandwidth)
        object.__setattr__(self, "resolution", resolution)
        object.__setattr__(self, "band_index", band_index)
        object.__setattr__(self, "rsrf_wavelengths_nm", rsrf_wavelengths_nm)
        object.__setattr__(self, "rsrf_response", rsrf_response)
        object.__setattr__(self, "rsrf_sensor_unit_id", rsrf_sensor_unit_id)
        object.__setattr__(self, "rsrf_representation_variant", rsrf_representation_variant)
        object.__setattr__(self, "rsrf_band_id", rsrf_band_id)

    @property
    def wavelength_um(self) -> float:
        return self.center_wavelength / 1000.0

    @property
    def has_rsrf(self) -> bool:
        return self.rsrf_wavelengths_nm is not None and self.rsrf_response is not None

    @property
    def has_srf(self) -> bool:
        return self.has_rsrf

    @property
    def srf_wavelengths_nm(self) -> np.ndarray | None:
        return self.rsrf_wavelengths_nm

    @property
    def srf_response(self) -> np.ndarray | None:
        return self.rsrf_response


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
                raise ValueError(
                    f"Band {band.name!r} has non-positive resolution: {band.resolution}"
                )

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

    def select_bands_in_range(self, wl_min_nm: float, wl_max_nm: float) -> list[SensorBand]:
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

    def default_aerosol_solver_bands(self) -> list[SensorBand]:
        """Return the default bands used for aerosol retrieval on this sensor."""
        preferred_by_sensor = {
            "MSI": ("B02", "B04"),
        }
        preferred_names = preferred_by_sensor.get(self.sensor_id, ())
        preferred = [band for name in preferred_names for band in self.bands if band.name == name]
        if len(preferred) == len(preferred_names) and preferred:
            return preferred

        aerosol = self.select_bands_in_range(400.0, 520.0)
        if aerosol:
            return aerosol
        return list(self.bands[:2])

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
