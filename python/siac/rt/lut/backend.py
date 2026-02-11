"""Zarr LUT backend implementation."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from siac.core.types import AtmosphericState, GeometryAngles, RTCoefficients, SensorBand
from siac.rt.lut.constants import LUT_COORD_DIMS
from siac.rt.lut.store import as_local_path, build_lut_store

logger = logging.getLogger(__name__)


class ZarrLUTBackend:
    """
    Look-up table backend for RT coefficients using Zarr storage.

    Implements the RTModelBackend protocol using multi-dimensional
    interpolation in a pre-computed look-up table.

    Args:
        lut_path: Path to the Zarr LUT store
        interpolation_method: Interpolation method ("linear" or "nearest")
        chunk_cache_size: Size of chunk cache in bytes
    """

    def __init__(
        self,
        lut_path: str | Path,
        interpolation_method: str = "linear",
        chunk_cache_size: int = 128 * 1024 * 1024,  # 128 MB
        storage_options: dict[str, Any] | None = None,
    ):
        self.lut_path = str(lut_path)
        self.interpolation_method = interpolation_method
        self.chunk_cache_size = chunk_cache_size
        self.storage_options = dict(storage_options or {})

        # Lazy load the LUT
        self._lut: xr.Dataset | None = None
        self._lut_coords: dict[str, np.ndarray] = {}

    @property
    def lut(self) -> xr.Dataset:
        """Lazily load the LUT dataset."""
        if self._lut is None:
            self._load_lut()
        return self._lut

    def _load_lut(self) -> None:
        """Load the Zarr LUT dataset."""
        local_path = as_local_path(self.lut_path)
        if local_path is not None and not local_path.exists():
            raise FileNotFoundError(f"LUT not found at {local_path}")

        logger.info(f"Loading LUT from {self.lut_path}")

        # Open with dask for lazy loading. Store may be local path, fsspec mapper, or zip mapper.
        store = build_lut_store(self.lut_path, self.storage_options)
        if self._store_contains_key(store, "zarr.json"):
            os.environ.setdefault("ZARR_V3_EXPERIMENTAL_API", "1")
            try:
                self._lut = xr.open_zarr(
                    store=store,
                    consolidated=False,
                    zarr_version=3,
                )
            except Exception as e:
                logger.warning(
                    "Failed to open LUT as zarr v3 (%s); retrying with zarr v2 paths",
                    e,
                )

        if self._lut is None:
            try:
                self._lut = xr.open_zarr(store=store, consolidated=True)
            except Exception as e:
                logger.warning(
                    "Failed to open LUT with consolidated metadata (%s); retrying with consolidated=False",
                    e,
                )
                self._lut = xr.open_zarr(store=store, consolidated=False)

        # Cache coordinate arrays for fast interpolation
        for dim in LUT_COORD_DIMS:
            if dim in self._lut.coords:
                self._lut_coords[dim] = self._lut.coords[dim].values

        logger.info(
            f"LUT loaded with dimensions: "
            f"SZA={len(self._lut_coords.get('sza', []))}, "
            f"VZA={len(self._lut_coords.get('vza', []))}, "
            f"AOT={len(self._lut_coords.get('aot', []))}"
        )

    @staticmethod
    def _store_contains_key(store: Any, key: str) -> bool:
        """Best-effort key existence check for mapping/path style stores."""
        if isinstance(store, str):
            return (Path(store) / key).exists()
        try:
            return key in store
        except Exception:
            return False

    def compute_coefficients(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
        compute_jacobian: bool = False,
    ) -> RTCoefficients:
        """
        Compute radiative transfer coefficients via LUT interpolation.

        Args:
            geometry: Viewing geometry (sza, vza, raa in radians)
            atmo_state: Atmospheric state (aot, tcwv, tco3, elevation)
            band: Sensor band specification
            compute_jacobian: Whether to compute d(coeff)/d(aot,tcwv)

        Returns:
            RTCoefficients with xap, xbp, xcp and optional Jacobians.
        """
        # Convert geometry to degrees for LUT
        sza_deg = np.rad2deg(geometry.sza.values)
        vza_deg = np.rad2deg(geometry.vza.values)
        raa_deg = np.abs(np.rad2deg(geometry.raa.values))

        # Normalize RAA to [0, 180]
        raa_deg = np.where(raa_deg > 180, 360 - raa_deg, raa_deg)

        # Get atmospheric parameters
        aot = atmo_state.aot.values
        tcwv = atmo_state.tcwv.values
        elevation = atmo_state.elevation.values

        # Get band wavelength for LUT selection
        wavelength = band.center_wavelength

        # Interpolate LUT
        path_ref, trans_down, trans_up, sph_alb = self._interpolate_lut(
            sza_deg, vza_deg, raa_deg, aot, tcwv, wavelength
        )

        # Apply elevation correction
        path_ref, trans_down, trans_up, sph_alb = self._apply_elevation_correction(
            path_ref, trans_down, trans_up, sph_alb, elevation
        )

        # Convert to RT coefficients
        # The correction equation is:
        #   y = xap * toa - xbp
        #   boa = y / (1 + xcp * y)
        #
        # From RT parameters:
        #   xap = 1 / (trans_down * trans_up)
        #   xbp = path_ref / (trans_down * trans_up)
        #   xcp = sph_alb / trans_up

        trans_total = trans_down * trans_up
        trans_total = np.maximum(trans_total, 1e-10)  # Avoid division by zero

        xap = 1.0 / trans_total
        xbp = path_ref / trans_total
        xcp = sph_alb / np.maximum(trans_up, 1e-10)

        # Create DataArrays
        template = geometry.sza
        xap_da = xr.DataArray(xap, dims=template.dims, coords=template.coords)
        xbp_da = xr.DataArray(xbp, dims=template.dims, coords=template.coords)
        xcp_da = xr.DataArray(xcp, dims=template.dims, coords=template.coords)

        # Compute Jacobians via finite differences if requested
        d_xap = None
        d_xbp = None
        d_xcp = None

        if compute_jacobian:
            d_xap, d_xbp, d_xcp = self._compute_jacobian_numerical(
                geometry, atmo_state, band, xap, xbp, xcp
            )

        return RTCoefficients(
            xap=xap_da,
            xbp=xbp_da,
            xcp=xcp_da,
            d_xap=d_xap,
            d_xbp=d_xbp,
            d_xcp=d_xcp,
        )

    def _interpolate_lut(
        self,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
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

        original_shape = sza.shape

        # Flatten inputs
        sza_flat = sza.ravel()
        vza_flat = vza.ravel()
        raa_flat = raa.ravel()
        aot_flat = aot.ravel()
        tcwv_flat = tcwv.ravel()

        n_pixels = len(sza_flat)

        # Find nearest wavelength in LUT
        wl_idx = np.argmin(np.abs(self._lut_coords["wavelength"] - wavelength))
        wl_value = self._lut_coords["wavelength"][wl_idx]

        logger.debug(f"Using LUT wavelength {wl_value:.1f}nm for band at {wavelength:.1f}nm")

        # Select wavelength slice
        lut_wl = self.lut.sel(wavelength=wl_value, method="nearest")

        # Initialize output arrays
        path_ref = np.zeros(n_pixels, dtype=np.float32)
        trans_down = np.zeros(n_pixels, dtype=np.float32)
        trans_up = np.zeros(n_pixels, dtype=np.float32)
        sph_alb = np.zeros(n_pixels, dtype=np.float32)

        # Perform interpolation for each pixel
        # For efficiency, we batch pixels with similar coordinates
        if self.interpolation_method == "linear":
            path_ref, trans_down, trans_up, sph_alb = self._linear_interpolate(
                lut_wl, sza_flat, vza_flat, raa_flat, aot_flat, tcwv_flat
            )
        else:
            path_ref, trans_down, trans_up, sph_alb = self._nearest_interpolate(
                lut_wl, sza_flat, vza_flat, raa_flat, aot_flat, tcwv_flat
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
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Multi-dimensional linear interpolation."""
        n_pixels = len(sza)

        # Create coordinate arrays for interpolation
        coords = {
            "sza": xr.DataArray(sza, dims=["point"]),
            "vza": xr.DataArray(vza, dims=["point"]),
            "raa": xr.DataArray(raa, dims=["point"]),
            "aot": xr.DataArray(aot, dims=["point"]),
            "tcwv": xr.DataArray(tcwv, dims=["point"]),
        }

        # Interpolate each variable
        path_ref = lut["path_reflectance"].interp(**coords).values.astype(np.float32)
        trans_down = lut["transmittance_down"].interp(**coords).values.astype(np.float32)
        trans_up = lut["transmittance_up"].interp(**coords).values.astype(np.float32)
        sph_alb = lut["spherical_albedo"].interp(**coords).values.astype(np.float32)

        return path_ref, trans_down, trans_up, sph_alb

    def _nearest_interpolate(
        self,
        lut: xr.Dataset,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Nearest-neighbor interpolation."""
        coords = {
            "sza": xr.DataArray(sza, dims=["point"]),
            "vza": xr.DataArray(vza, dims=["point"]),
            "raa": xr.DataArray(raa, dims=["point"]),
            "aot": xr.DataArray(aot, dims=["point"]),
            "tcwv": xr.DataArray(tcwv, dims=["point"]),
        }

        path_ref = lut["path_reflectance"].sel(**coords, method="nearest").values.astype(np.float32)
        trans_down = lut["transmittance_down"].sel(**coords, method="nearest").values.astype(np.float32)
        trans_up = lut["transmittance_up"].sel(**coords, method="nearest").values.astype(np.float32)
        sph_alb = lut["spherical_albedo"].sel(**coords, method="nearest").values.astype(np.float32)

        return path_ref, trans_down, trans_up, sph_alb

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
        # Scale height for Rayleigh scattering (~8.5 km)
        scale_height = 8.5

        # Elevation correction factor
        correction = np.exp(-elevation / scale_height)

        # Apply correction (reduce path effects at altitude)
        path_ref_corr = path_ref * correction
        sph_alb_corr = sph_alb * correction

        # Transmittance increases with altitude
        trans_correction = np.power(correction, 0.5)  # Approximate
        trans_down_corr = 1.0 - (1.0 - trans_down) * trans_correction
        trans_up_corr = 1.0 - (1.0 - trans_up) * trans_correction

        return path_ref_corr, trans_down_corr, trans_up_corr, sph_alb_corr

    def _compute_jacobian_numerical(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
        xap: np.ndarray,
        xbp: np.ndarray,
        xcp: np.ndarray,
    ) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
        """
        Compute Jacobians via numerical finite differences.

        Args:
            geometry: Viewing geometry
            atmo_state: Atmospheric state
            band: Sensor band
            xap, xbp, xcp: Current coefficient values

        Returns:
            Tuple of Jacobian DataArrays for (d_xap, d_xbp, d_xcp)
        """
        # Perturbation size
        delta_aot = 0.01
        delta_tcwv = 0.1

        template = geometry.sza
        original_shape = template.shape

        # Perturb AOT
        aot_plus = atmo_state.with_updated_aot_tcwv(
            aot=atmo_state.aot + delta_aot,
            tcwv=atmo_state.tcwv,
        )

        coeffs_aot_plus = self.compute_coefficients(
            geometry, aot_plus, band, compute_jacobian=False
        )

        d_xap_aot = (coeffs_aot_plus.xap.values - xap) / delta_aot
        d_xbp_aot = (coeffs_aot_plus.xbp.values - xbp) / delta_aot
        d_xcp_aot = (coeffs_aot_plus.xcp.values - xcp) / delta_aot

        # Perturb TCWV
        tcwv_plus = atmo_state.with_updated_aot_tcwv(
            aot=atmo_state.aot,
            tcwv=atmo_state.tcwv + delta_tcwv,
        )

        coeffs_tcwv_plus = self.compute_coefficients(
            geometry, tcwv_plus, band, compute_jacobian=False
        )

        d_xap_tcwv = (coeffs_tcwv_plus.xap.values - xap) / delta_tcwv
        d_xbp_tcwv = (coeffs_tcwv_plus.xbp.values - xbp) / delta_tcwv
        d_xcp_tcwv = (coeffs_tcwv_plus.xcp.values - xcp) / delta_tcwv

        # Create DataArrays with param dimension
        d_xap = xr.concat(
            [
                xr.DataArray(d_xap_aot, dims=template.dims, coords=template.coords),
                xr.DataArray(d_xap_tcwv, dims=template.dims, coords=template.coords),
            ],
            dim="param",
        ).assign_coords(param=["aot", "tcwv"])

        d_xbp = xr.concat(
            [
                xr.DataArray(d_xbp_aot, dims=template.dims, coords=template.coords),
                xr.DataArray(d_xbp_tcwv, dims=template.dims, coords=template.coords),
            ],
            dim="param",
        ).assign_coords(param=["aot", "tcwv"])

        d_xcp = xr.concat(
            [
                xr.DataArray(d_xcp_aot, dims=template.dims, coords=template.coords),
                xr.DataArray(d_xcp_tcwv, dims=template.dims, coords=template.coords),
            ],
            dim="param",
        ).assign_coords(param=["aot", "tcwv"])

        return d_xap, d_xbp, d_xcp

    def supports_jacobian(self) -> bool:
        """Check if this backend supports Jacobian computation."""
        return True  # Via numerical differentiation

    @property
    def backend_name(self) -> str:
        """Name of the RT backend."""
        return "lut"

    def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
        """Check if this backend has data for the specified sensor."""
        # LUT backend works for any sensor with known wavelengths
        return True

    def get_available_wavelengths(self) -> np.ndarray:
        """Get wavelengths available in the LUT."""
        return self._lut_coords.get("wavelength", np.array([]))
