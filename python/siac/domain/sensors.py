"""Pure sensor domain types."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np


# Default tolerance (nm) when callers ask for a "nearest" band without a
# specific window in mind. 50 nm covers typical S2/Landsat spacing without
# bridging adjacent VIS bands.
_DEFAULT_NEAREST_BAND_TOLERANCE_NM = 50.0


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
    ) -> None:
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


@dataclass(frozen=True)
class SensorConfig:
    """Complete sensor configuration for a satellite instrument."""

    sensor_id: str
    satellite_id: str
    bands: tuple[SensorBand, ...]
    default_ref_scale: float = 1.0 / 10000.0
    default_ref_offset: float = 0.0
    # Band names preferred for aerosol retrieval on this sensor. When set,
    # default_aerosol_solver_bands() returns these bands in declaration order.
    # Catalog entries (e.g. Sentinel-2 MSI) populate this; the generic
    # fallback is documented on default_aerosol_solver_bands().
    aerosol_solver_band_names: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if not self.bands:
            raise ValueError("SensorConfig must have at least one band.")
        names = [b.name for b in self.bands]
        # O(N) duplicate detection via Counter (replaces O(N^2) names.count loop).
        counts = Counter(names)
        dupes = sorted(n for n, c in counts.items() if c > 1)
        if dupes:
            raise ValueError(f"Duplicate band names in SensorConfig: {dupes}")
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
        """Return the band whose center is closest to ``wavelength_nm`` within
        ``tolerance_nm``. Single shared implementation for nearest-band lookup;
        :meth:`select_nearest_band` is a thin alias with a wider default tolerance.
        """
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
        self,
        target_nm: float,
        tolerance_nm: float = _DEFAULT_NEAREST_BAND_TOLERANCE_NM,
    ) -> SensorBand | None:
        """Backwards-compatible alias for :meth:`get_band_by_wavelength` with a
        wider default tolerance (``_DEFAULT_NEAREST_BAND_TOLERANCE_NM``).

        Implementation is shared with :meth:`get_band_by_wavelength`; only the
        default tolerance differs.
        """
        return self.get_band_by_wavelength(target_nm, tolerance_nm=tolerance_nm)

    def default_aerosol_solver_bands(self) -> list[SensorBand]:
        """Return the default bands used for aerosol retrieval on this sensor.

        Resolution order:

        1. ``aerosol_solver_band_names`` set on the catalog entry — returned
           in declaration order. Used by S2 MSI (``("B02", "B04")``), OLI
           (``("B1", "B2")``), and any externally registered sensors.
        2. Legacy fallback for ``sensor_id == "MSI"`` configs without the
           catalog field — still picks ``("B02", "B04")`` if both bands
           exist. This is preserved only so test mocks that build
           ``SensorConfig(sensor_id="MSI", ...)`` ad-hoc keep working;
           the canonical path is to set ``aerosol_solver_band_names``.
           TODO(REVIEW.md §3.2): drop this branch once downstream tests
           migrate to the catalog field.
        3. Generic 400-520 nm (aerosol-sensitive blue) wavelength window —
           the historical behaviour for non-catalogued sensors.
        4. First two bands as a final fallback.

        The previous implementation hard-coded ``sensor_id == "MSI"`` inside
        this dataclass; that branching now lives in the catalog
        (``aerosol_solver_band_names``) — see REVIEW.md §3.2.
        """
        if self.aerosol_solver_band_names:
            band_lookup = {b.name: b for b in self.bands}
            preferred = [
                band_lookup[name] for name in self.aerosol_solver_band_names if name in band_lookup
            ]
            if len(preferred) == len(self.aerosol_solver_band_names) and preferred:
                return preferred

        # Step 2 (legacy MSI fallback). Catalog-built S2 configs hit step 1
        # above; this branch only matters for ad-hoc inline MSI configs.
        if self.sensor_id == "MSI":
            band_lookup = {b.name: b for b in self.bands}
            legacy = [band_lookup[n] for n in ("B02", "B04") if n in band_lookup]
            if len(legacy) == 2:
                return legacy

        # Step 3: aerosol-sensitive blue window. Matches the pre-refactor
        # non-MSI behaviour exactly (no behaviour change for uncatalogued
        # sensors).
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
