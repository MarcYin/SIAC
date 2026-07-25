"""Convolve the TOA observation with the sensor PSF and co-register it to the prior.

The MODIS-derived surface prior is a coarse, PSF-weighted view of the surface. To
compare the fine Sentinel-2 TOA against it for the aerosol solve, the TOA must be
convolved with the (fixed-width) MODIS PSF and aligned to the prior grid — exactly
as the original SIAC v1 did (``the_aerosol.py::_get_convolved_toa`` +
``psf_optimize.py``). v2 had moved the PSF onto the prior instead; this module
restores the physically-consistent observation-side convolution plus a per-scene
integer co-registration shift (accepted only when the TOA↔prior correlation clears
``min_correlation``, mirroring v1's ``1 - r >= 0.6`` gate).

The PSF width is fixed (the calibrated ``psf_sigma_x``/``psf_sigma_y`` in pixels at
``target_resolution_m``); only the integer shift is fitted per scene.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import cast

import numpy as np
import xarray as xr

from siac._rust_compat import PSFConvolver
from siac.algorithms.surface.kernel_model import KernelModelDeriver
from siac.runtime.models import copy_spatial_metadata_like

logger = logging.getLogger(__name__)

FloatArray = np.ndarray

#: Minimum overlap pixels required to score a candidate shift (guards tiny windows).
_MIN_SHIFT_PIXELS = 64


@dataclass(frozen=True)
class ToaPsfConfig:
    """Settings for observation-side PSF convolution and integer shift fitting."""

    sigma_x: float = 29.75
    sigma_y: float = 39.0
    target_resolution_m: float = 10.0
    convolve_resolution: str = "native10m"  # "native10m" (default) | "grid" (coarse/fast)
    shift_search_radius_m: float = 500.0  # ≈ v1's ±50 px @ 10 m
    min_correlation: float = 0.6
    reference_bands: tuple[str, ...] = ("B12", "B11")
    enabled: bool = True


@dataclass(frozen=True)
class ShiftFitResult:
    """Outcome of the integer co-registration shift fit."""

    dx: int = 0
    dy: int = 0
    correlation: float = field(default=float("nan"))
    accepted: bool = False


def _scaled_sigmas(cfg: ToaPsfConfig, grid_resolution_m: float) -> tuple[float, float]:
    """Rescale the fixed pixel sigmas to the working grid (reuse the prior-path helper).

    ``sigma_x``/``sigma_y`` count pixels at ``target_resolution_m``; on a grid of
    ``grid_resolution_m`` the same physical blur spans
    ``sigma * target_resolution_m / grid_resolution_m`` pixels. Delegating to
    :meth:`KernelModelDeriver._scale_psf_sigmas` keeps the footprint identical to
    the legacy prior-side convolution.
    """
    helper = KernelModelDeriver(
        psf_sigma_x=cfg.sigma_x,
        psf_sigma_y=cfg.sigma_y,
        target_resolution_m=cfg.target_resolution_m,
    )
    return helper._scale_psf_sigmas(cfg.sigma_x, cfg.sigma_y, grid_resolution_m)


def convolve_toa_band(arr2d: FloatArray, sigma_x: float, sigma_y: float) -> FloatArray:
    """NaN-aware separable Gaussian convolution of a 2-D band via the Rust PSFConvolver."""
    if sigma_x <= 0.0 or sigma_y <= 0.0:
        return cast("FloatArray", np.asarray(arr2d, dtype=np.float32))
    conv = PSFConvolver(sigma_x=float(sigma_x), sigma_y=float(sigma_y))
    out = conv.convolve(np.asarray(arr2d, dtype=np.float64))
    return cast("FloatArray", np.asarray(out, dtype=np.float32))


def convolve_toa_cube(toa_da: xr.DataArray, sigma_x: float, sigma_y: float) -> xr.DataArray:
    """Convolve every band of a ``(band, y, x)`` (or ``(y, x)``) TOA DataArray."""
    if "band" in toa_da.dims:
        stacks = [
            convolve_toa_band(toa_da.isel(band=i).values, sigma_x, sigma_y)
            for i in range(toa_da.sizes["band"])
        ]
        data: FloatArray = np.stack(stacks, axis=0)
    else:
        data = convolve_toa_band(toa_da.values, sigma_x, sigma_y)
    out = xr.DataArray(data, dims=toa_da.dims, coords=toa_da.coords, attrs=toa_da.attrs)
    return copy_spatial_metadata_like(out, toa_da)


def _overlap_slices(length: int, shift: int) -> tuple[slice, slice]:
    """Return (prior_slice, toa_slice) along one axis for an aligned shift.

    ``aligned_toa[i] == toa[i - shift]``, so ``prior[i]`` pairs with ``toa[i-shift]``
    over the in-bounds overlap.
    """
    if shift >= 0:
        return slice(shift, length), slice(0, length - shift)
    return slice(0, length + shift), slice(-shift, length)


def _pearson(a: FloatArray, b: FloatArray) -> float:
    """Pearson correlation of two flat arrays; NaN if degenerate."""
    if a.size < _MIN_SHIFT_PIXELS:
        return float("nan")
    a = a - a.mean()
    b = b - b.mean()
    denom = float(np.sqrt(np.sum(a * a) * np.sum(b * b)))
    if denom <= 0.0:
        return float("nan")
    return float(np.sum(a * b) / denom)


def fit_integer_shift(
    toa_ref: FloatArray,
    prior_ref: FloatArray,
    valid: FloatArray,
    *,
    search_radius_px: int,
    min_correlation: float,
    subsample: int = 1,
) -> ShiftFitResult:
    """Deterministic exhaustive 2-D integer co-registration shift.

    For each integer ``(dy, dx)`` in ``[-R, R]`` the convolved TOA is aligned to the
    prior and scored by the Pearson correlation over the jointly-valid overlap. The
    best shift is accepted only if its correlation clears ``min_correlation`` (v1's
    ``1 - r >= 0.6`` gate); otherwise the scene falls back to ``(0, 0)``.

    ``subsample`` strides the correlation sampling (the shift granularity is
    unchanged) so a fine-resolution fit stays tractable — like v1, which correlated
    only at the sparse MODIS-grid points rather than every native pixel.
    """
    height, width = prior_ref.shape
    finite_prior = np.isfinite(prior_ref)
    finite_toa = np.isfinite(toa_ref)
    valid_b = np.asarray(valid, dtype=bool) & finite_prior
    s = max(1, int(subsample))

    best_r = -np.inf
    best_dx = 0
    best_dy = 0
    for dy in range(-search_radius_px, search_radius_px + 1):
        py, ty = _overlap_slices(height, dy)
        for dx in range(-search_radius_px, search_radius_px + 1):
            px, tx = _overlap_slices(width, dx)
            m = valid_b[py, px][::s, ::s] & finite_toa[ty, tx][::s, ::s]
            if int(m.sum()) < _MIN_SHIFT_PIXELS:
                continue
            r = _pearson(prior_ref[py, px][::s, ::s][m], toa_ref[ty, tx][::s, ::s][m])
            if np.isfinite(r) and r > best_r:
                best_r, best_dx, best_dy = r, dx, dy

    if not np.isfinite(best_r):
        return ShiftFitResult(0, 0, float("nan"), False)
    if best_r < min_correlation:
        return ShiftFitResult(0, 0, float(best_r), False)
    return ShiftFitResult(best_dx, best_dy, float(best_r), True)


def apply_integer_shift(toa_da: xr.DataArray, dx: int, dy: int) -> xr.DataArray:
    """Shift a ``(band, y, x)`` cube by integer ``(dy, dx)``; blank wrapped edges to NaN."""
    if dx == 0 and dy == 0:
        return toa_da
    y_axis = toa_da.get_axis_num("y")
    x_axis = toa_da.get_axis_num("x")
    data = np.asarray(toa_da.values, dtype=np.float32)
    shifted = np.roll(data, shift=(dy, dx), axis=(y_axis, x_axis))

    # np.roll wraps; overwrite the wrapped border with NaN so no opposite-edge
    # data bleeds across the scene (the solver's isfinite mask then drops them).
    height = data.shape[y_axis]
    width = data.shape[x_axis]
    y_idx: list[int] = list(range(0, dy)) if dy >= 0 else list(range(height + dy, height))
    x_idx: list[int] = list(range(0, dx)) if dx >= 0 else list(range(width + dx, width))
    if y_idx:
        sl: list[slice | list[int]] = [slice(None)] * shifted.ndim
        sl[y_axis] = y_idx
        shifted[tuple(sl)] = np.nan
    if x_idx:
        sl = [slice(None)] * shifted.ndim
        sl[x_axis] = x_idx
        shifted[tuple(sl)] = np.nan

    out = xr.DataArray(shifted, dims=toa_da.dims, coords=toa_da.coords, attrs=toa_da.attrs)
    return copy_spatial_metadata_like(out, toa_da)


def psf_convolve_and_align_toa(
    toa_da: xr.DataArray,
    toa_ref_band: xr.DataArray | None,
    prior_ref_band: xr.DataArray | None,
    prior_valid: FloatArray | None,
    *,
    grid_resolution_m: float,
    cfg: ToaPsfConfig,
) -> tuple[xr.DataArray, ShiftFitResult]:
    """Convolve the gridded TOA cube and apply the fitted integer co-registration shift.

    ``toa_ref_band`` / ``prior_ref_band`` are the gridded SWIR reference (TOA and
    prior) used only for the shift fit; if either is missing the TOA is convolved
    but not shifted.
    """
    if not cfg.enabled:
        return toa_da, ShiftFitResult()

    sigma_x, sigma_y = _scaled_sigmas(cfg, grid_resolution_m)
    conv = convolve_toa_cube(toa_da, sigma_x, sigma_y)

    if toa_ref_band is None or prior_ref_band is None:
        return conv, ShiftFitResult()

    conv_ref = convolve_toa_band(np.asarray(toa_ref_band.values), sigma_x, sigma_y)
    prior_ref = np.asarray(prior_ref_band.values, dtype=np.float64)
    if conv_ref.shape != prior_ref.shape:
        logger.warning(
            "PSF shift reference shape mismatch (toa %s vs prior %s); skipping shift.",
            conv_ref.shape,
            prior_ref.shape,
        )
        return conv, ShiftFitResult()
    valid = np.isfinite(prior_ref) if prior_valid is None else np.asarray(prior_valid, dtype=bool)

    radius = max(0, int(round(cfg.shift_search_radius_m / max(grid_resolution_m, 1e-6))))
    # Subsample the correlation to keep a fine-resolution fit tractable: the
    # exhaustive (2R+1)^2 search would otherwise scan every native pixel per shift.
    subsample = max(1, int(round(np.sqrt(prior_ref.size / 40000.0))))
    fit = fit_integer_shift(
        conv_ref,
        prior_ref,
        valid,
        search_radius_px=radius,
        min_correlation=cfg.min_correlation,
        subsample=subsample,
    )
    if not fit.accepted:
        return conv, fit
    return apply_integer_shift(conv, fit.dx, fit.dy), fit
