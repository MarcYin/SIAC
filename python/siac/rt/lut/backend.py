"""Zarr LUT backend implementation."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from siac.core.types import AtmosphericState, GeometryAngles, RTCoefficients, SensorBand
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
        for dim, coord in self._lut.coords.items():
            self._lut_coords[dim] = coord.values

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
            xap, xbp, xcp = self._compute_coefficients_from_spectral_lut(
                sza_deg,
                vza_deg,
                raa_deg,
                aot,
                tcwv,
                atmo_state.tco3.values,
                elevation,
                band,
            )
        else:
            raise ValueError(
                "LUT does not contain a supported RT representation. "
                "Expected compact coefficient variables or dense spectral LUT variables."
            )

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

        original_shape = sza.shape

        # Flatten inputs
        sza_flat = sza.ravel()
        vza_flat = vza.ravel()
        raa_flat = raa.ravel()
        aot_flat = aot.ravel()
        tcwv_flat = tcwv.ravel()

        # Find nearest wavelength in LUT
        wl_idx = np.argmin(np.abs(self._lut_coords["wavelength"] - wavelength))
        wl_value = self._lut_coords["wavelength"][wl_idx]

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
            "sza": xr.DataArray(sza, dims=["point"]),
            "vza": xr.DataArray(vza, dims=["point"]),
            "raa": xr.DataArray(raa, dims=["point"]),
            "aot": xr.DataArray(aot, dims=["point"]),
            "tcwv": xr.DataArray(tcwv, dims=["point"]),
        }

        if "ozone" in self._lut_coords:
            ozone = np.asarray(tco3, dtype=np.float32)
            ozone_axis = np.asarray(self._lut_coords["ozone"], dtype=np.float32)
            # Atmospheric ozone often arrives in atm-cm (~0.3), while LUTs may use DU (~300).
            if ozone_axis.size and np.nanmax(np.abs(ozone_axis)) > 20 and np.nanmax(np.abs(ozone)) < 10:
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
            coords[name] = xr.DataArray(
                np.clip(coord.values, float(np.nanmin(axis)), float(np.nanmax(axis))),
                dims=["point"],
            )

        return coords

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

        coords = self._build_point_coords(
            sza=sza.ravel(),
            vza=vza.ravel(),
            raa=raa.ravel(),
            aot=aot.ravel(),
            tcwv=tcwv.ravel(),
            tco3=tco3.ravel(),
            elevation=elevation.ravel(),
        )
        weights = self._spectral_integration_weights(band)

        toa_rho1 = self._weighted_spectral_mean(
            self._interpolate_variable(self.lut["TOA_rho1"], coords, self.interpolation_method),
            weights,
        ).values.astype(np.float32)
        toa_rho2 = self._weighted_spectral_mean(
            self._interpolate_variable(self.lut["TOA_rho2"], coords, self.interpolation_method),
            weights,
        ).values.astype(np.float32)
        eg_rho1 = self._weighted_spectral_mean(
            self._interpolate_variable(self.lut["Eg_rho1"], coords, self.interpolation_method),
            weights,
        ).values.astype(np.float32)
        eg_rho2 = self._weighted_spectral_mean(
            self._interpolate_variable(self.lut["Eg_rho2"], coords, self.interpolation_method),
            weights,
        ).values.astype(np.float32)

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
            xap.reshape(sza.shape).astype(np.float32),
            xbp.reshape(sza.shape).astype(np.float32),
            xcp.reshape(sza.shape).astype(np.float32),
        )

    def _spectral_integration_weights(self, band: SensorBand) -> xr.DataArray:
        """Build wavelength weights for spectral convolution (bandpass * optional solar spectrum)."""
        wavelength = xr.DataArray(
            np.asarray(self._lut_coords["wavelength"], dtype=np.float32),
            dims=["wavelength"],
            coords={"wavelength": self._lut_coords["wavelength"]},
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
            if name not in self.lut:
                continue
            solar = self.lut[name]
            if "wavelength" not in solar.dims:
                continue
            extra_dims = [dim for dim in solar.dims if dim != "wavelength"]
            if extra_dims:
                solar = solar.mean(dim=extra_dims)
            weights = bandpass * solar.astype(np.float32)
            break

        return weights

    @staticmethod
    def _weighted_spectral_mean(data: xr.DataArray, weights: xr.DataArray) -> xr.DataArray:
        """Weighted mean over wavelength with coordinate-aware integration."""
        if "wavelength" not in data.dims:
            return data

        local_weights = weights.sel(wavelength=data["wavelength"])
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
