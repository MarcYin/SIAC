"""Compact coefficient-LUT path of :class:`~siac.algorithms.rt.lut.backend.ZarrLUTBackend`.

Implements the legacy COMPACT representation (pre-derived ``path_reflectance``
/ ``transmittance_down`` / ``transmittance_up`` / ``spherical_albedo``
variables): per-pixel point-indexer construction, multi-dimensional
interpolation of the four coefficient variables, and the Beer-Lambert
elevation correction for LUTs without an explicit altitude axis. The mixin is
composed into ``ZarrLUTBackend`` (which owns loading, dispatch, and the
shared sanitizer/interpolation helpers declared in the ``TYPE_CHECKING``
block below).
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, cast

import numpy as np
import xarray as xr

from siac.algorithms.rt.lut.constants import normalize_compact_ozone

logger = logging.getLogger(__name__)


class _CompactLUTMixin:
    """Compact xap/xbp/xcp coefficient-LUT methods for ``ZarrLUTBackend``."""

    _COEFFICIENT_VARS = (
        "path_reflectance",
        "transmittance_down",
        "transmittance_up",
        "spherical_albedo",
    )

    if TYPE_CHECKING:
        # Attributes and shared helpers provided by the ZarrLUTBackend facade.
        interpolation_method: str
        _lut_coords: dict[str, np.ndarray]

        @property
        def lut(self) -> xr.Dataset: ...

        @staticmethod
        def _interpolate_variable_fast(
            var: xr.DataArray,
            coords: dict[str, xr.DataArray],
            method: str,
        ) -> np.ndarray: ...

        @staticmethod
        def _require_finite_values(values: np.ndarray, *, name: str) -> np.ndarray: ...

        @classmethod
        def _sanitize_point_values(cls, values: np.ndarray, axis: np.ndarray) -> np.ndarray: ...

    def _interpolate_lut(
        self,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
        wavelength: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Interpolate LUT to get RT parameters.

        Returns:
            Tuple of (path_reflectance, trans_down, trans_up, sph_albedo)
        """
        # Ensure LUT and coordinate cache are initialized.
        if not self._lut_coords:
            _ = self.lut
        wavelength_axis = np.asarray(
            self._lut_coords.get("wavelength", np.array([], dtype=np.float32)),
            dtype=np.float32,
        )
        if wavelength_axis.size == 0:
            raise ValueError("Coefficient LUT must define a non-empty wavelength coordinate")

        original_shape = sza.shape

        # Flatten inputs
        sza_flat = sza.ravel()
        vza_flat = vza.ravel()
        raa_flat = raa.ravel()
        aot_flat = aot.ravel()
        tcwv_flat = tcwv.ravel()

        # Find nearest wavelength in LUT
        wl_idx = np.argmin(np.abs(wavelength_axis - wavelength))
        wl_value = wavelength_axis[wl_idx]

        logger.debug(f"Using LUT wavelength {wl_value:.1f}nm for band at {wavelength:.1f}nm")

        # Select wavelength slice
        lut_wl = self.lut.sel(wavelength=wl_value, method="nearest")
        coords = self._build_point_coords(
            sza=sza_flat,
            vza=vza_flat,
            raa=raa_flat,
            aot=aot_flat,
            tcwv=tcwv_flat,
            tco3=tco3.ravel(),
            elevation=elevation.ravel(),
        )

        # Perform interpolation for each pixel
        # For efficiency, we batch pixels with similar coordinates
        if self.interpolation_method == "linear":
            path_ref, trans_down, trans_up, sph_alb = self._linear_interpolate(
                lut_wl,
                coords,
            )
        else:
            path_ref, trans_down, trans_up, sph_alb = self._nearest_interpolate(
                lut_wl,
                coords,
            )

        return (
            path_ref.reshape(original_shape),
            trans_down.reshape(original_shape),
            trans_up.reshape(original_shape),
            sph_alb.reshape(original_shape),
        )

    def _linear_interpolate(
        self,
        lut: xr.Dataset,
        coords: dict[str, xr.DataArray],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Multi-dimensional linear interpolation."""
        path_ref = self._interpolate_variable_fast(lut["path_reflectance"], coords, "linear")
        trans_down = self._interpolate_variable_fast(lut["transmittance_down"], coords, "linear")
        trans_up = self._interpolate_variable_fast(lut["transmittance_up"], coords, "linear")
        sph_alb = self._interpolate_variable_fast(lut["spherical_albedo"], coords, "linear")

        return (path_ref, trans_down, trans_up, sph_alb)

    def _nearest_interpolate(
        self,
        lut: xr.Dataset,
        coords: dict[str, xr.DataArray],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Nearest-neighbor interpolation."""
        path_ref = self._interpolate_variable_fast(lut["path_reflectance"], coords, "nearest")
        trans_down = self._interpolate_variable_fast(lut["transmittance_down"], coords, "nearest")
        trans_up = self._interpolate_variable_fast(lut["transmittance_up"], coords, "nearest")
        sph_alb = self._interpolate_variable_fast(lut["spherical_albedo"], coords, "nearest")

        return (path_ref, trans_down, trans_up, sph_alb)

    @staticmethod
    def _point_indexer(values: np.ndarray) -> xr.DataArray:
        """Wrap flattened interpolation values in the LUT point dimension."""
        return xr.DataArray(np.asarray(values, dtype=np.float32), dims=["point"])

    @classmethod
    def _required_point_indexer(
        cls,
        values: np.ndarray,
        *,
        name: str,
    ) -> xr.DataArray:
        """Build a point indexer for a required LUT interpolation axis."""
        return cls._point_indexer(cls._require_finite_values(values, name=name))

    def _normalize_ozone_point_values(self, tco3: np.ndarray) -> np.ndarray:
        """Normalize ozone values to the LUT unit convention when needed.

        SIAC carries ozone in atm-cm. Some LUTs index ozone in DU, so values are
        upscaled by 1000 when the LUT axis is clearly in the DU regime.
        """
        ozone_axis = np.asarray(self._lut_coords["ozone"], dtype=np.float32)
        return cast(
            "np.ndarray", normalize_compact_ozone(np.asarray(tco3, dtype=np.float32), ozone_axis)
        )

    def _sanitize_point_indexers(self, coords: dict[str, xr.DataArray]) -> dict[str, xr.DataArray]:
        """Clamp available point indexers to the finite domain of each LUT axis."""
        sanitized = dict(coords)
        for name, coord in list(sanitized.items()):
            if name not in self._lut_coords:
                continue
            axis = np.asarray(self._lut_coords[name], dtype=np.float32)
            if axis.size == 0:
                continue
            sanitized[name] = self._point_indexer(self._sanitize_point_values(coord.values, axis))
        return sanitized

    def _build_point_coords(
        self,
        *,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
    ) -> dict[str, xr.DataArray]:
        """Build point-indexer coordinates for interpolation.

        Expected units:

        - ``sza``, ``vza``, ``raa`` in degrees
        - ``aot`` unitless
        - ``tcwv`` in cm precipitable water
        - ``tco3`` in atm-cm; converted to DU automatically when the LUT ozone
          axis uses DU
        - ``elevation`` in km
        """
        coords: dict[str, xr.DataArray] = {
            "sza": self._required_point_indexer(sza, name="sza"),
            "vza": self._required_point_indexer(vza, name="vza"),
            "raa": self._required_point_indexer(raa, name="raa"),
            "aot": self._required_point_indexer(aot, name="aot"),
            "tcwv": self._required_point_indexer(tcwv, name="tcwv"),
        }

        if "ozone" in self._lut_coords:
            coords["ozone"] = self._point_indexer(self._normalize_ozone_point_values(tco3))

        if "altitude" in self._lut_coords:
            coords["altitude"] = self._point_indexer(elevation)

        return self._sanitize_point_indexers(coords)

    def _coefficient_lut_has_altitude_axis(self) -> bool:
        """Return True when coefficient LUT variables already model altitude explicitly."""
        return all(
            "altitude" in self.lut[name].dims for name in self._COEFFICIENT_VARS if name in self.lut
        )

    def _apply_elevation_correction(
        self,
        path_ref: np.ndarray,
        trans_down: np.ndarray,
        trans_up: np.ndarray,
        sph_alb: np.ndarray,
        elevation: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Apply elevation correction to RT parameters.

        Uses Beer-Lambert law approximation for pressure-altitude relationship.

        Args:
            path_ref: Path reflectance
            trans_down: Downward transmittance
            trans_up: Upward transmittance
            sph_alb: Spherical albedo
            elevation: Surface elevation in km

        Returns:
            Corrected RT parameters
        """
        # Atmospheric Rayleigh scale height — see siac.constants.
        from siac.constants import ATMOSPHERIC_SCALE_HEIGHT_KM

        # Elevation correction factor
        correction = np.exp(-elevation / ATMOSPHERIC_SCALE_HEIGHT_KM)

        # Apply correction (reduce path effects at altitude)
        path_ref_corr = path_ref * correction
        sph_alb_corr = sph_alb * correction

        # Transmittance increases with altitude
        trans_correction = np.power(correction, 0.5)  # Approximate
        trans_down_corr = 1.0 - (1.0 - trans_down) * trans_correction
        trans_up_corr = 1.0 - (1.0 - trans_up) * trans_correction

        return path_ref_corr, trans_down_corr, trans_up_corr, sph_alb_corr
