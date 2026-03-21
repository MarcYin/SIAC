"""Zarr LUT backend implementation."""

from __future__ import annotations

import logging
import os
import threading
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeVar, cast

import numpy as np
import xarray as xr

from siac.algorithms.rt.lut.rsrf_kernel import build_aligned_rsrf_kernel
from siac.algorithms.rt.lut.store import as_local_path, build_lut_store
from siac.domain.spectral import RelativeSpectralResponse
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients

if TYPE_CHECKING:
    from collections.abc import Callable

    from siac.domain import SensorBand

logger = logging.getLogger(__name__)
_T = TypeVar("_T")


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

    _COEFFICIENT_VARS = (
        "path_reflectance",
        "transmittance_down",
        "transmittance_up",
        "spherical_albedo",
    )
    _SPECTRAL_VARS = (
        "TOA_rho1",
        "TOA_rho2",
        "Eg_rho1",
        "Eg_rho2",
    )
    _SOLAR_IRRADIANCE_NAMES = (
        "extraterrestrial_solar_irradiance",
        "solar_irradiance",
    )
    _DEFAULT_SURFACE_RHO1 = 0.15
    _DEFAULT_SURFACE_RHO2 = 0.5
    _SPECTRAL_IO_MAX_RETRIES = 8
    _SUPPORTED_INTERPOLATION_METHODS = frozenset({"linear", "nearest"})

    def __init__(
        self,
        lut_path: str | Path,
        interpolation_method: str = "linear",
        chunk_cache_size: int = 128 * 1024 * 1024,  # 128 MB
        storage_options: dict[str, Any] | None = None,
    ):
        self.lut_path = str(lut_path)
        if interpolation_method not in self._SUPPORTED_INTERPOLATION_METHODS:
            supported = ", ".join(sorted(self._SUPPORTED_INTERPOLATION_METHODS))
            raise ValueError(
                f"Unsupported LUT interpolation_method {interpolation_method!r}; "
                f"expected one of: {supported}"
            )
        self.interpolation_method = interpolation_method
        self.chunk_cache_size = chunk_cache_size
        self.storage_options = dict(storage_options or {})
        self._scene_subset_logged = False

        # Lazy load the LUT
        self._lut: xr.Dataset | None = None
        self._lut_coords: dict[str, np.ndarray] = {}
        self._cache_lock = threading.Lock()
        self._spectral_scene_key: tuple[float, ...] | None = None
        self._spectral_scene_subset: xr.Dataset | None = None
        self._spectral_band_grid_cache: dict[
            tuple[tuple[float, ...], tuple[str, float, float]],
            tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray],
        ] = {}

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
                    zarr_format=3,
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
        for dim, coord in self._lut.coords.items():
            self._lut_coords[dim] = coord.values

        logger.info(
            "LUT index loaded (available axis sizes): SZA=%d, VZA=%d, RAA=%d, AOT=%d, TCWV=%d",
            len(self._lut_coords.get("sza", [])),
            len(self._lut_coords.get("vza", [])),
            len(self._lut_coords.get("raa", [])),
            len(self._lut_coords.get("aot", [])),
            len(self._lut_coords.get("tcwv", [])),
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
        template = self._require_matching_grid_shapes(geometry, atmo_state)

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

        if self._supports_coefficient_lut():
            # Interpolate compact LUT variables, then derive RT coefficients.
            path_ref, trans_down, trans_up, sph_alb = self._interpolate_lut(
                sza_deg,
                vza_deg,
                raa_deg,
                aot,
                tcwv,
                atmo_state.tco3.values,
                elevation,
                wavelength,
            )

            # Legacy compact LUTs may not carry altitude as an explicit dimension.
            if not self._coefficient_lut_has_altitude_axis():
                path_ref, trans_down, trans_up, sph_alb = self._apply_elevation_correction(
                    path_ref,
                    trans_down,
                    trans_up,
                    sph_alb,
                    elevation,
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
        elif self._supports_spectral_lut():
            xap, xbp, xcp = self._compute_spectral_with_retry(
                sza_deg=sza_deg,
                vza_deg=vza_deg,
                raa_deg=raa_deg,
                aot=aot,
                tcwv=tcwv,
                tco3=atmo_state.tco3.values,
                elevation=elevation,
                band=band,
            )
        else:
            raise ValueError(
                "LUT does not contain a supported RT representation. "
                "Expected compact coefficient variables or dense spectral LUT variables."
            )

        # Create DataArrays
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
        # Interpolate each variable
        path_ref = self._interpolate_variable(lut["path_reflectance"], coords, "linear")
        trans_down = self._interpolate_variable(lut["transmittance_down"], coords, "linear")
        trans_up = self._interpolate_variable(lut["transmittance_up"], coords, "linear")
        sph_alb = self._interpolate_variable(lut["spherical_albedo"], coords, "linear")

        return (
            path_ref.values.astype(np.float32),
            trans_down.values.astype(np.float32),
            trans_up.values.astype(np.float32),
            sph_alb.values.astype(np.float32),
        )

    def _nearest_interpolate(
        self,
        lut: xr.Dataset,
        coords: dict[str, xr.DataArray],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Nearest-neighbor interpolation."""
        path_ref = self._interpolate_variable(lut["path_reflectance"], coords, "nearest")
        trans_down = self._interpolate_variable(lut["transmittance_down"], coords, "nearest")
        trans_up = self._interpolate_variable(lut["transmittance_up"], coords, "nearest")
        sph_alb = self._interpolate_variable(lut["spherical_albedo"], coords, "nearest")

        return (
            path_ref.values.astype(np.float32),
            trans_down.values.astype(np.float32),
            trans_up.values.astype(np.float32),
            sph_alb.values.astype(np.float32),
        )

    @staticmethod
    def _axis_bounds(axis: np.ndarray) -> tuple[float, float]:
        axis_values = np.asarray(axis, dtype=np.float32)
        return float(np.nanmin(axis_values)), float(np.nanmax(axis_values))

    @classmethod
    def _axis_midpoint(cls, axis: np.ndarray) -> float:
        lo, hi = cls._axis_bounds(axis)
        return (lo + hi) * 0.5

    @staticmethod
    def _finite_range(
        values: np.ndarray,
        *,
        fallback: tuple[float, float],
    ) -> tuple[float, float]:
        arr = np.asarray(values, dtype=np.float32)
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            return fallback
        return float(np.min(finite)), float(np.max(finite))

    @staticmethod
    def _finite_mean(values: np.ndarray, *, fallback: float) -> float:
        arr = np.asarray(values, dtype=np.float32)
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            return fallback
        return float(np.mean(finite))

    @classmethod
    def _sanitize_point_values(
        cls,
        values: np.ndarray,
        axis: np.ndarray,
    ) -> np.ndarray:
        axis_min, axis_max = cls._axis_bounds(axis)
        midpoint = cls._axis_midpoint(axis)
        arr = np.asarray(values, dtype=np.float32)
        sanitized = np.where(np.isfinite(arr), arr, midpoint)
        clipped = np.asarray(np.clip(sanitized, axis_min, axis_max), dtype=np.float32)
        return cast("np.ndarray", clipped)

    @staticmethod
    def _require_finite_values(
        values: np.ndarray,
        *,
        name: str,
    ) -> np.ndarray:
        arr = np.asarray(values, dtype=np.float32)
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} must contain only finite values for LUT interpolation")
        return cast("np.ndarray", arr)

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
        """Build point-indexer coordinates for interpolation."""
        coords: dict[str, xr.DataArray] = {
            "sza": xr.DataArray(self._require_finite_values(sza, name="sza"), dims=["point"]),
            "vza": xr.DataArray(self._require_finite_values(vza, name="vza"), dims=["point"]),
            "raa": xr.DataArray(self._require_finite_values(raa, name="raa"), dims=["point"]),
            "aot": xr.DataArray(self._require_finite_values(aot, name="aot"), dims=["point"]),
            "tcwv": xr.DataArray(self._require_finite_values(tcwv, name="tcwv"), dims=["point"]),
        }

        if "ozone" in self._lut_coords:
            ozone = np.asarray(tco3, dtype=np.float32)
            ozone_axis = np.asarray(self._lut_coords["ozone"], dtype=np.float32)
            # Atmospheric ozone often arrives in atm-cm (~0.3), while LUTs may use DU (~300).
            finite_ozone = ozone[np.isfinite(ozone)]
            if ozone_axis.size and finite_ozone.size and np.nanmax(np.abs(ozone_axis)) > 20 and np.nanmax(np.abs(finite_ozone)) < 10:
                ozone = ozone * 1000.0
            coords["ozone"] = xr.DataArray(ozone, dims=["point"])

        if "altitude" in self._lut_coords:
            coords["altitude"] = xr.DataArray(
                np.asarray(elevation, dtype=np.float32),
                dims=["point"],
            )

        for name, coord in list(coords.items()):
            if name not in self._lut_coords:
                continue
            axis = np.asarray(self._lut_coords[name], dtype=np.float32)
            if axis.size == 0:
                continue
            coords[name] = xr.DataArray(self._sanitize_point_values(coord.values, axis), dims=["point"])

        return coords

    @staticmethod
    def _require_matching_grid_shapes(
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
    ) -> xr.DataArray:
        """Validate that geometry and atmospheric fields share the same retrieval grid."""
        template = geometry.sza
        expected_shape = tuple(template.shape)
        arrays = {
            "geometry.sza": geometry.sza,
            "geometry.vza": geometry.vza,
            "geometry.raa": geometry.raa,
            "atmo_state.aot": atmo_state.aot,
            "atmo_state.tcwv": atmo_state.tcwv,
            "atmo_state.tco3": atmo_state.tco3,
            "atmo_state.elevation": atmo_state.elevation,
        }
        mismatches = [
            f"{name}={tuple(array.shape)}"
            for name, array in arrays.items()
            if tuple(array.shape) != expected_shape
        ]
        if mismatches:
            details = ", ".join(mismatches)
            raise ValueError(
                "Geometry and atmospheric fields must share the same grid shape; "
                f"expected {expected_shape} from geometry.sza but got {details}"
            )
        return template

    def _compute_spectral_with_retry(
        self,
        *,
        sza_deg: np.ndarray,
        vza_deg: np.ndarray,
        raa_deg: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
        band: SensorBand,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute spectral coefficients with retry for transient remote LUT read failures."""
        return self._run_with_transient_lut_io_retry(
            lambda: self._compute_coefficients_from_spectral_lut(
                sza_deg,
                vza_deg,
                raa_deg,
                aot,
                tcwv,
                tco3,
                elevation,
                band,
            ),
            operation="read",
        )

    @staticmethod
    def _interpolate_variable(
        var: xr.DataArray,
        coords: dict[str, xr.DataArray],
        method: str,
    ) -> xr.DataArray:
        """Interpolate one LUT variable using only coordinates it actually depends on."""
        applicable = {name: coord for name, coord in coords.items() if name in var.dims}
        if not applicable:
            return var
        if method == "linear":
            return var.interp(**applicable)
        return var.sel(**applicable, method="nearest")

    def _compute_coefficients_from_spectral_lut(
        self,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
        band: SensorBand,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Derive standard RT coefficients from dense spectral LUT terms."""
        if "wavelength" not in self._lut_coords:
            raise ValueError("Spectral LUT must define a wavelength coordinate")

        scene_key, lut_scene = self._get_or_build_spectral_scene_subset(
            sza=sza,
            vza=vza,
            raa=raa,
            tco3=tco3,
            elevation=elevation,
        )
        if not self._scene_subset_logged:
            tco3_bounds = self._finite_range(
                np.asarray(tco3, dtype=np.float32),
                fallback=(0.0, 0.0),
            )
            elevation_bounds = self._finite_range(
                np.asarray(elevation, dtype=np.float32),
                fallback=(0.0, 0.0),
            )
            logger.info(
                (
                    "Scene LUT subset: sza=%.3f deg, vza=%.3f deg, raa=%.3f deg, "
                    "ozone=[%.3f, %.3f], altitude=[%.3f, %.3f]"
                ),
                self._finite_mean(sza, fallback=0.0),
                self._finite_mean(vza, fallback=0.0),
                self._finite_mean(raa, fallback=0.0),
                tco3_bounds[0],
                tco3_bounds[1],
                elevation_bounds[0],
                elevation_bounds[1],
            )
            self._scene_subset_logged = True

        target_shape = aot.shape
        coords = self._build_aot_tcwv_point_coords(
            lut_scene,
            aot=aot.ravel(),
            tcwv=tcwv.ravel(),
        )
        toa_rho1_grid, toa_rho2_grid, eg_rho1_grid, eg_rho2_grid = self._get_or_build_spectral_band_grids(
            scene_key,
            lut_scene,
            band,
        )

        toa_rho1 = self._interpolate_variable(toa_rho1_grid, coords, self.interpolation_method).values.astype(
            np.float32
        )
        toa_rho2 = self._interpolate_variable(toa_rho2_grid, coords, self.interpolation_method).values.astype(
            np.float32
        )
        eg_rho1 = self._interpolate_variable(eg_rho1_grid, coords, self.interpolation_method).values.astype(
            np.float32
        )
        eg_rho2 = self._interpolate_variable(eg_rho2_grid, coords, self.interpolation_method).values.astype(
            np.float32
        )

        rho1, rho2 = self._spectral_reference_reflectances()
        eps = 1e-10

        denom = rho2 * eg_rho2 - rho1 * eg_rho1
        denom = np.where(np.abs(denom) < eps, eps, denom)

        s_term = (eg_rho2 - eg_rho1) / denom
        path_ref = (
            toa_rho2 * rho1 * eg_rho1 - toa_rho1 * rho2 * eg_rho2
        ) / np.where(np.abs(rho1 * eg_rho1 - rho2 * eg_rho2) < eps, eps, rho1 * eg_rho1 - rho2 * eg_rho2)
        t_up = (toa_rho2 - toa_rho1) / denom
        eg0 = eg_rho1 * (1.0 - rho1 * s_term)
        t_total = np.maximum(eg0 * t_up, eps)

        xap = 1.0 / t_total
        xbp = path_ref / t_total
        xcp = s_term

        return (
            xap.reshape(target_shape).astype(np.float32),
            xbp.reshape(target_shape).astype(np.float32),
            xcp.reshape(target_shape).astype(np.float32),
        )

    @classmethod
    def _spectral_scene_cache_key(
        cls,
        *,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
    ) -> tuple[float, ...]:
        """Build a stable scene cache key from summary geometry/atmosphere stats."""
        tco3_arr = np.asarray(tco3, dtype=np.float32)
        elevation_arr = np.asarray(elevation, dtype=np.float32)
        return (
            round(cls._finite_mean(sza, fallback=0.0), 3),
            round(cls._finite_mean(vza, fallback=0.0), 3),
            round(cls._finite_mean(raa, fallback=0.0), 3),
            round(cls._finite_range(tco3_arr, fallback=(0.0, 0.0))[0], 3),
            round(cls._finite_range(tco3_arr, fallback=(0.0, 0.0))[1], 3),
            round(cls._finite_range(elevation_arr, fallback=(0.0, 0.0))[0], 3),
            round(cls._finite_range(elevation_arr, fallback=(0.0, 0.0))[1], 3),
        )

    @staticmethod
    def _spectral_band_cache_key(band: SensorBand) -> tuple[str, float, float]:
        """Build a stable per-band cache key for spectral compression grids."""
        return (band.name, float(band.center_wavelength), float(band.bandwidth))

    def _get_or_build_spectral_scene_subset(
        self,
        *,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
    ) -> tuple[tuple[float, ...], xr.Dataset]:
        """Return cached scene LUT subset or build and cache it once."""
        scene_key = self._spectral_scene_cache_key(
            sza=sza,
            vza=vza,
            raa=raa,
            tco3=tco3,
            elevation=elevation,
        )
        with self._cache_lock:
            if (
                self._spectral_scene_key == scene_key
                and self._spectral_scene_subset is not None
            ):
                return scene_key, self._spectral_scene_subset

        subset = self._subset_spectral_lut_for_scene(
            self.lut,
            sza=sza,
            vza=vza,
            raa=raa,
            tco3=tco3,
            elevation=elevation,
        )

        with self._cache_lock:
            if self._spectral_scene_key != scene_key:
                self._spectral_scene_key = scene_key
                self._spectral_scene_subset = subset
                self._spectral_band_grid_cache.clear()
                self._scene_subset_logged = False
            if self._spectral_scene_subset is None:
                self._spectral_scene_subset = subset
            return scene_key, self._spectral_scene_subset

    def _get_or_build_spectral_band_grids(
        self,
        scene_key: tuple[float, ...],
        lut_scene: xr.Dataset,
        band: SensorBand,
    ) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
        """Return cached per-band spectral-compressed grids or build once."""
        band_key = self._spectral_band_cache_key(band)
        cache_key = (scene_key, band_key)

        with self._cache_lock:
            cached = self._spectral_band_grid_cache.get(cache_key)
        if cached is not None:
            return cached

        lut_band = self._subset_wavelength_for_band(lut_scene, band)
        weights = self._spectral_integration_weights(band, lut_band)
        grids = (
            self._weighted_spectral_mean(lut_band["TOA_rho1"], weights),
            self._weighted_spectral_mean(lut_band["TOA_rho2"], weights),
            self._weighted_spectral_mean(lut_band["Eg_rho1"], weights),
            self._weighted_spectral_mean(lut_band["Eg_rho2"], weights),
        )

        with self._cache_lock:
            self._spectral_band_grid_cache[cache_key] = grids
        return grids

    def preload_scene_subset(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand] | None = None,
    ) -> None:
        """
        Preload LUT metadata and cache scene/band spectral subsets.

        Intended for background execution while other modules are running.
        """
        self._run_with_transient_lut_io_retry(
            lambda: self._preload_scene_subset_once(
                geometry=geometry,
                atmo_state=atmo_state,
                bands=bands,
            ),
            operation="preload",
        )

    def _preload_scene_subset_once(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand] | None,
    ) -> None:
        """Populate scene-level spectral LUT caches once for the requested bands."""
        # Ensure LUT metadata index is loaded exactly once.
        _ = self.lut
        if not self._supports_spectral_lut():
            return
        self._require_matching_grid_shapes(geometry, atmo_state)

        sza = np.rad2deg(geometry.sza.values)
        vza = np.rad2deg(geometry.vza.values)
        raa = np.abs(np.rad2deg(geometry.raa.values))
        raa = np.where(raa > 180, 360 - raa, raa)
        tco3 = atmo_state.tco3.values
        elevation = atmo_state.elevation.values

        scene_key, lut_scene = self._get_or_build_spectral_scene_subset(
            sza=sza,
            vza=vza,
            raa=raa,
            tco3=tco3,
            elevation=elevation,
        )
        for band in bands or []:
            self._get_or_build_spectral_band_grids(scene_key, lut_scene, band)

    def _run_with_transient_lut_io_retry(
        self,
        fn: Callable[[], _T],
        *,
        operation: str,
    ) -> _T:
        """Retry transient remote LUT I/O for preload/read paths."""
        max_attempts = max(1, int(self._SPECTRAL_IO_MAX_RETRIES))
        for attempt in range(1, max_attempts + 1):
            try:
                return fn()
            except Exception as exc:
                if attempt >= max_attempts or not self._is_transient_lut_io_error(exc):
                    raise
                logger.warning(
                    "Transient LUT %s failure (%s); reloading and retrying (%d/%d).",
                    operation,
                    type(exc).__name__,
                    attempt,
                    max_attempts,
                )
                time.sleep(min(float(attempt), 5.0))
                self._reload_lut()

        raise RuntimeError(f"Failed to {operation} LUT data after retries")

    def _subset_spectral_lut_for_scene(
        self,
        lut: xr.Dataset,
        *,
        sza: np.ndarray,
        vza: np.ndarray,
        raa: np.ndarray,
        tco3: np.ndarray,
        elevation: np.ndarray,
    ) -> xr.Dataset:
        """Subset spectral LUT using scene-angle means and ozone/elevation ranges."""
        out = lut

        for dim, scene_values in {"sza": sza, "vza": vza, "raa": raa}.items():
            if dim in out.coords:
                axis = np.asarray(out.coords[dim].values, dtype=np.float32)
                finite_scene_values = self._require_finite_values(scene_values, name=dim)
                value = self._finite_mean(finite_scene_values, fallback=self._axis_midpoint(axis))
                out = out.sel({dim: value}, method="nearest")
            if dim in out.dims and out.sizes.get(dim, 0) == 1:
                out = out.squeeze(dim=dim, drop=True)

        for dim, scene_values in {"ozone": tco3, "altitude": elevation}.items():
            if dim not in out.coords:
                continue
            axis = np.asarray(out.coords[dim].values, dtype=np.float32)
            if axis.size == 0:
                continue
            axis_min, axis_max = self._axis_bounds(axis)
            vmin, vmax = self._finite_range(scene_values, fallback=(axis_min, axis_max))
            lo = float(np.clip(min(vmin, vmax), axis_min, axis_max))
            hi = float(np.clip(max(vmin, vmax), axis_min, axis_max))
            subset = out.sel({dim: slice(lo, hi)})
            if subset.sizes.get(dim, 0) == 0:
                subset = out.sel({dim: (lo + hi) * 0.5}, method="nearest")
            if dim in subset.dims and subset.sizes.get(dim, 1) > 1:
                subset = subset.mean(dim=dim)
            if dim in subset.dims and subset.sizes.get(dim, 0) == 1:
                subset = subset.squeeze(dim=dim, drop=True)
            out = subset

        return out

    @staticmethod
    def _build_aot_tcwv_point_coords(
        lut: xr.Dataset,
        *,
        aot: np.ndarray,
        tcwv: np.ndarray,
    ) -> dict[str, xr.DataArray]:
        """Build point interpolation coordinates for variables solved per-pixel."""
        coords = {
            "aot": xr.DataArray(ZarrLUTBackend._require_finite_values(aot, name="aot"), dims=["point"]),
            "tcwv": xr.DataArray(ZarrLUTBackend._require_finite_values(tcwv, name="tcwv"), dims=["point"]),
        }
        for name in ("aot", "tcwv"):
            if name not in lut.coords:
                continue
            axis = np.asarray(lut.coords[name].values, dtype=np.float32)
            if axis.size == 0:
                continue
            coords[name] = xr.DataArray(ZarrLUTBackend._sanitize_point_values(coords[name].values, axis), dims=["point"])
        return coords

    def _subset_wavelength_for_band(
        self,
        lut: xr.Dataset,
        band: SensorBand,
        sigma_factor: float = 4.0,
    ) -> xr.Dataset:
        """Trim spectral LUT to a narrow wavelength window around the target band."""
        if "wavelength" not in lut.coords:
            return lut

        wl_axis = np.asarray(lut.coords["wavelength"].values, dtype=np.float32)
        if wl_axis.size == 0:
            return lut

        if band.has_rsrf:
            kernel = build_aligned_rsrf_kernel(
                self._band_rsrf(band),
                lut_wavelengths_nm=wl_axis,
                lut_id=self.lut_path,
                support_padding=1,
            )
            return lut.isel(wavelength=slice(kernel.start_index, kernel.end_index))

        sigma = max(float(band.bandwidth) / (2.0 * np.sqrt(2.0 * np.log(2.0))), 1e-6)
        half_window = max(3.0, sigma_factor * sigma)
        wl_min = float(np.nanmin(wl_axis))
        wl_max = float(np.nanmax(wl_axis))
        lo = float(np.clip(float(band.center_wavelength) - half_window, wl_min, wl_max))
        hi = float(np.clip(float(band.center_wavelength) + half_window, wl_min, wl_max))

        subset = lut.sel(wavelength=slice(min(lo, hi), max(lo, hi)))
        if subset.sizes.get("wavelength", 0) == 0:
            nearest_idx = int(np.argmin(np.abs(wl_axis - float(band.center_wavelength))))
            subset = lut.isel(wavelength=slice(nearest_idx, nearest_idx + 1))
        return subset

    def _spectral_integration_weights(self, band: SensorBand, lut: xr.Dataset | None = None) -> xr.DataArray:
        """Build wavelength weights for spectral convolution (bandpass * optional solar spectrum)."""
        source = lut if lut is not None else self.lut
        if "wavelength" not in source.coords:
            raise ValueError("Spectral LUT must define a wavelength coordinate")

        wl_axis = np.asarray(source.coords["wavelength"].values, dtype=np.float32)
        wavelength = xr.DataArray(
            wl_axis,
            dims=["wavelength"],
            coords={"wavelength": wl_axis},
        )
        if band.has_rsrf:
            solar_values = None
            for name in self._SOLAR_IRRADIANCE_NAMES:
                if name not in source:
                    continue
                solar = source[name]
                if "wavelength" not in solar.dims:
                    continue
                extra_dims = [dim for dim in solar.dims if dim != "wavelength"]
                if extra_dims:
                    solar = solar.mean(dim=extra_dims)
                solar_values = np.asarray(solar.values, dtype=np.float32)
                break

            kernel = build_aligned_rsrf_kernel(
                self._band_rsrf(band),
                lut_wavelengths_nm=wl_axis,
                lut_id=self.lut_path,
                solar_irradiance=solar_values,
                support_padding=0,
            )
            rsrf_weights = kernel.solar_weighted_response_on_lut
            if rsrf_weights is None:
                rsrf_weights = kernel.response_on_lut
            full_weights = np.zeros_like(wl_axis, dtype=np.float32)
            full_weights[kernel.start_index:kernel.end_index] = rsrf_weights
            return xr.DataArray(
                full_weights,
                dims=["wavelength"],
                coords={"wavelength": wl_axis},
            )

        sigma = max(
            float(band.bandwidth) / (2.0 * np.sqrt(2.0 * np.log(2.0))),
            1e-6,
        )
        bandpass = np.exp(
            -0.5 * np.square((wavelength - float(band.center_wavelength)) / sigma)
        ).astype(np.float32)
        weights: xr.DataArray = bandpass

        for name in self._SOLAR_IRRADIANCE_NAMES:
            if name not in source:
                continue
            solar = source[name]
            if "wavelength" not in solar.dims:
                continue
            extra_dims = [dim for dim in solar.dims if dim != "wavelength"]
            if extra_dims:
                solar = solar.mean(dim=extra_dims)
            weights = bandpass * solar.astype(np.float32)
            break

        return weights

    @staticmethod
    def _band_rsrf(band: SensorBand) -> RelativeSpectralResponse:
        """Build a canonical RSRF object from a band-carried tabulated response."""
        if not band.has_rsrf:
            raise ValueError(f"Band {band.name!r} does not define a tabulated RSRF")
        return RelativeSpectralResponse.from_tabulated(
            sensor_id="UNKNOWN",
            satellite_id="UNKNOWN",
            band_name=band.name,
            wavelengths_nm=np.asarray(band.rsrf_wavelengths_nm, dtype=np.float32),
            response=np.asarray(band.rsrf_response, dtype=np.float32),
        )

    @staticmethod
    def _weighted_spectral_mean(data: xr.DataArray, weights: xr.DataArray) -> xr.DataArray:
        """Weighted mean over wavelength with coordinate-aware integration."""
        if "wavelength" not in data.dims:
            return data

        local_weights = weights.reindex(wavelength=data["wavelength"], fill_value=0.0)
        if data.sizes["wavelength"] == 1:
            numerator = (data * local_weights).isel(wavelength=0, drop=True)
            denominator = local_weights.isel(wavelength=0, drop=True)
            return numerator / xr.where(np.abs(denominator) < 1e-10, 1e-10, denominator)

        numerator = (data * local_weights).integrate("wavelength")
        denominator = local_weights.integrate("wavelength")
        denominator = xr.where(np.abs(denominator) < 1e-10, 1e-10, denominator)
        return numerator / denominator

    def _spectral_reference_reflectances(self) -> tuple[float, float]:
        """Return reference surface reflectances used by dense spectral LUT variables."""
        attrs = self.lut.attrs
        rho1 = float(attrs.get("rho1", attrs.get("reference_reflectance_1", self._DEFAULT_SURFACE_RHO1)))
        rho2 = float(attrs.get("rho2", attrs.get("reference_reflectance_2", self._DEFAULT_SURFACE_RHO2)))
        return rho1, rho2

    def _reload_lut(self) -> None:
        """Drop cached handles and reopen LUT metadata/chunks lazily."""
        self._lut = None
        self._lut_coords = {}
        self._scene_subset_logged = False
        with self._cache_lock:
            self._spectral_scene_key = None
            self._spectral_scene_subset = None
            self._spectral_band_grid_cache.clear()
        _ = self.lut

    @staticmethod
    def _is_transient_lut_io_error(exc: Exception) -> bool:
        """Return True for transient remote I/O errors that should be retried."""
        probes = []
        current: Exception | None = exc
        for _ in range(10):
            if current is None:
                break
            probes.append(current)
            nxt = current.__cause__ if isinstance(current.__cause__, Exception) else None
            if nxt is None and isinstance(current.__context__, Exception):
                nxt = current.__context__
            current = nxt

        tokens = (
            "referencenotreachable",
            "server disconnected",
            "serverdisconnectederror",
            "connection reset",
            "connection aborted",
            "timed out",
            "timeout",
            "temporarily unavailable",
        )
        for err in probes:
            name = type(err).__name__.lower()
            text = str(err).lower()
            if any(tok in name or tok in text for tok in tokens):
                return True
        return False

    def _supports_coefficient_lut(self) -> bool:
        """Check whether the loaded LUT exposes compact coefficient variables."""
        return all(name in self.lut.data_vars for name in self._COEFFICIENT_VARS)

    def _supports_spectral_lut(self) -> bool:
        """Check whether the loaded LUT exposes dense spectral radiative terms."""
        return all(name in self.lut.data_vars for name in self._SPECTRAL_VARS)

    def _coefficient_lut_has_altitude_axis(self) -> bool:
        """Return True when coefficient LUT variables already model altitude explicitly."""
        return all("altitude" in self.lut[name].dims for name in self._COEFFICIENT_VARS if name in self.lut)

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
        return False

    @property
    def backend_name(self) -> str:
        """Name of the RT backend."""
        return "lut"

    def is_available_for_sensor(self, _sensor_id: str, _satellite_id: str) -> bool:
        """Check if this backend has data for the specified sensor."""
        # LUT backend works for any sensor with known wavelengths
        return True

    def get_available_wavelengths(self) -> np.ndarray:
        """Get wavelengths available in the LUT."""
        wavelengths = self._lut_coords.get("wavelength")
        if wavelengths is None:
            return np.zeros(0, dtype=np.float32)
        result: np.ndarray = np.array(wavelengths, dtype=np.float32)
        return result
