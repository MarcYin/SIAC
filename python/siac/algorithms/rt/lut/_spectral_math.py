"""Helpers for dense spectral LUT compression and coefficient recovery."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from siac.adapters.rsrf import band_convolution_weights
from siac.algorithms.rt.lut.constants import state_tcwv_to_lut
from siac.algorithms.rt.lut.rsrf_kernel import build_aligned_rsrf_kernel

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping

    from siac.domain import SensorBand
    from siac.domain.spectral import RelativeSpectralResponse


def build_point_interpolation_coords(
    lut: xr.Dataset,
    *,
    aot: np.ndarray,
    tcwv: np.ndarray,
    require_finite_values: Callable[..., np.ndarray],
    sanitize_point_values: Callable[[np.ndarray, np.ndarray], np.ndarray],
) -> dict[str, xr.DataArray]:
    """Build point interpolation coordinates for variables solved per-pixel.

    ``aot`` is unitless on both sides. ``tcwv`` arrives in SIAC state units
    (cm) and is converted to the LUT axis unit (mm) BEFORE clipping to the
    axis envelope — querying cm values against the mm axis silently evaluated
    the LUT at a ~10x too-dry atmosphere (retrieved AOT biased low).
    """
    coords = {
        "aot": xr.DataArray(require_finite_values(aot, name="aot"), dims=["point"]),
        "tcwv": xr.DataArray(
            state_tcwv_to_lut(require_finite_values(tcwv, name="tcwv")),
            dims=["point"],
        ),
    }
    for name in ("aot", "tcwv"):
        if name not in lut.coords:
            continue
        axis = np.asarray(lut.coords[name].values, dtype=np.float32)
        if axis.size == 0:
            continue
        coords[name] = xr.DataArray(
            sanitize_point_values(coords[name].values, axis),
            dims=["point"],
        )
    return coords


def build_spectral_integration_weights(
    source: xr.Dataset,
    band: SensorBand,
    *,
    lut_path: str,
    solar_irradiance_names: tuple[str, ...],
    band_rsrf: Callable[[SensorBand], RelativeSpectralResponse],
) -> xr.DataArray:
    """Build wavelength weights for spectral convolution."""
    if "wavelength" not in source.coords:
        raise ValueError("Spectral LUT must define a wavelength coordinate")

    wl_axis = np.asarray(source.coords["wavelength"].values, dtype=np.float32)
    if band.has_rsrf:
        solar_values = None
        for name in solar_irradiance_names:
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
            band_rsrf(band),
            lut_wavelengths_nm=wl_axis,
            lut_id=lut_path,
            solar_irradiance=solar_values,
            support_padding=0,
        )
        rsrf_weights = kernel.solar_weighted_response_on_lut
        if rsrf_weights is None:
            rsrf_weights = kernel.response_on_lut
        full_weights = np.zeros_like(wl_axis, dtype=np.float32)
        full_weights[kernel.start_index : kernel.end_index] = rsrf_weights
        return xr.DataArray(
            full_weights,
            dims=["wavelength"],
            coords={"wavelength": wl_axis},
        )

    weights = xr.DataArray(
        band_convolution_weights(band, wl_axis),
        dims=["wavelength"],
        coords={"wavelength": wl_axis},
    )

    for name in solar_irradiance_names:
        if name not in source:
            continue
        solar = source[name]
        if "wavelength" not in solar.dims:
            continue
        extra_dims = [dim for dim in solar.dims if dim != "wavelength"]
        if extra_dims:
            solar = solar.mean(dim=extra_dims)
        weights = weights * solar.astype(np.float32)
        break

    return weights


def weighted_spectral_mean(data: xr.DataArray, weights: xr.DataArray) -> xr.DataArray:
    """Weighted mean over wavelength with coordinate-aware integration."""
    if "wavelength" not in data.dims:
        return data

    local_weights = weights.reindex(wavelength=data["wavelength"], fill_value=0.0)
    eps = 1e-10
    if data.sizes["wavelength"] == 1:
        numerator = (data * local_weights).isel(wavelength=0, drop=True)
        denominator = local_weights.isel(wavelength=0, drop=True)
        # REVIEW.md §1.2 #1: preserve sign when clamping near-zero
        # denominators. Previously a small negative ``denominator`` was
        # replaced with ``+eps``, flipping the result's sign.
        sign = xr.where(denominator >= 0, 1.0, -1.0)
        safe_den = xr.where(np.abs(denominator) < eps, sign * eps, denominator)
        return numerator / safe_den

    numerator = (data * local_weights).integrate("wavelength")
    denominator = local_weights.integrate("wavelength")
    sign = xr.where(denominator >= 0, 1.0, -1.0)
    denominator = xr.where(np.abs(denominator) < eps, sign * eps, denominator)
    return numerator / denominator


def derive_standard_rt_coefficients(
    *,
    toa_rho1: np.ndarray,
    toa_rho2: np.ndarray,
    eg_rho1: np.ndarray,
    eg_rho2: np.ndarray,
    rho1: float,
    rho2: float,
    eps: float = 1e-10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert dense spectral radiative terms into standard RT coefficients."""
    denom = rho2 * eg_rho2 - rho1 * eg_rho1
    safe_denom = np.where(np.abs(denom) < eps, eps, denom)
    path_ref_denom = rho1 * eg_rho1 - rho2 * eg_rho2
    safe_path_ref_denom = np.where(np.abs(path_ref_denom) < eps, eps, path_ref_denom)

    s_term = (eg_rho2 - eg_rho1) / safe_denom
    path_ref = (toa_rho2 * rho1 * eg_rho1 - toa_rho1 * rho2 * eg_rho2) / safe_path_ref_denom
    t_up = (toa_rho2 - toa_rho1) / safe_denom
    eg0 = eg_rho1 * (1.0 - rho1 * s_term)
    t_total = np.maximum(eg0 * t_up, eps)

    return (
        (1.0 / t_total).astype(np.float32),
        (path_ref / t_total).astype(np.float32),
        s_term.astype(np.float32),
    )


def finite_range(
    values: np.ndarray,
    *,
    fallback: tuple[float, float],
) -> tuple[float, float]:
    """Return the min/max over finite values or a fallback when none exist."""
    arr = np.asarray(values, dtype=np.float32)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return fallback
    return float(np.min(finite)), float(np.max(finite))


def finite_mean(values: np.ndarray, *, fallback: float) -> float:
    """Return the mean over finite values or the provided fallback."""
    arr = np.asarray(values, dtype=np.float32)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return fallback
    return float(np.mean(finite))


def summarize_spectral_scene(
    *,
    sza: np.ndarray,
    vza: np.ndarray,
    raa: np.ndarray,
    tco3: np.ndarray,
    elevation: np.ndarray,
) -> tuple[float, float, float, tuple[float, float], tuple[float, float]]:
    """Build a rounded scene summary for cache keys and one-time logging."""
    tco3_arr = np.asarray(tco3, dtype=np.float32)
    elevation_arr = np.asarray(elevation, dtype=np.float32)
    return (
        round(finite_mean(sza, fallback=0.0), 3),
        round(finite_mean(vza, fallback=0.0), 3),
        round(finite_mean(raa, fallback=0.0), 3),
        (
            round(finite_range(tco3_arr, fallback=(0.0, 0.0))[0], 3),
            round(finite_range(tco3_arr, fallback=(0.0, 0.0))[1], 3),
        ),
        (
            round(finite_range(elevation_arr, fallback=(0.0, 0.0))[0], 3),
            round(finite_range(elevation_arr, fallback=(0.0, 0.0))[1], 3),
        ),
    )


def spectral_scene_cache_key(
    *,
    sza: np.ndarray,
    vza: np.ndarray,
    raa: np.ndarray,
    tco3: np.ndarray,
    elevation: np.ndarray,
    snap_axes: Mapping[str, np.ndarray] | None = None,
) -> tuple[float, ...]:
    """Build a stable cache key from the summarized scene state.

    ``snap_axes`` optionally maps ``sza``/``vza``/``raa`` to LUT coordinate
    axes; when given, each (unrounded) angle mean is snapped to its nearest
    axis node before rounding — mirroring the subset selection, so different
    spatial resolutions of the same scene (preload grid vs the coarse M3
    grid) produce the *same* key and share one download.
    """
    means = {
        "sza": finite_mean(sza, fallback=0.0),
        "vza": finite_mean(vza, fallback=0.0),
        "raa": finite_mean(raa, fallback=0.0),
    }
    if snap_axes:
        for dim, mean in means.items():
            axis = np.asarray(snap_axes.get(dim, ()), dtype=np.float64).ravel()
            if axis.size:
                means[dim] = float(axis[np.argmin(np.abs(axis - mean))])
    tco3_bounds = finite_range(np.asarray(tco3, dtype=np.float32), fallback=(0.0, 0.0))
    elevation_bounds = finite_range(np.asarray(elevation, dtype=np.float32), fallback=(0.0, 0.0))
    return (
        round(means["sza"], 3),
        round(means["vza"], 3),
        round(means["raa"], 3),
        round(tco3_bounds[0], 3),
        round(tco3_bounds[1], 3),
        round(elevation_bounds[0], 3),
        round(elevation_bounds[1], 3),
    )
