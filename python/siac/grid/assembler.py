"""
M4 Grid Assembler — resample and align all upstream outputs to solver grids.

Produces a ``SolverInputBundle`` where every raster field lives on the same
spatial grid (the *aux resolution*, typically 500 m).

See PLANS.md §4.4 for the full specification.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import xarray as xr
from scipy.ndimage import uniform_filter, zoom

from siac.core.types import (
    AtmosphericState,
    GeometryAngles,
    ObservationBundle,
    SolverInputBundle,
    SurfacePrior,
)

logger = logging.getLogger(__name__)


# ── Internal resampling helpers ────────────────────────────────────────

def _compute_target_shape(
    source_shape: tuple[int, int],
    source_resolution: float,
    target_resolution: float,
) -> tuple[int, int]:
    """Compute target (H, W) for a resolution change."""
    factor = source_resolution / target_resolution
    h = max(1, int(round(source_shape[0] * factor)))
    w = max(1, int(round(source_shape[1] * factor)))
    return (h, w)


def _resample_da(
    da: xr.DataArray,
    target_shape: tuple[int, int],
    method: str = "bilinear",
) -> xr.DataArray:
    """Resample a 2-D DataArray to *target_shape* using scipy zoom.

    Args:
        da: Input array with dims (y, x).
        target_shape: (H, W) of the output.
        method: ``"bilinear"`` (order=1), ``"nearest"`` (order=0), or
                ``"area"`` (block mean via uniform_filter then subsampling).
    """
    if da.ndim == 3 and "band" in da.dims:
        band_vals = da.coords["band"].values if "band" in da.coords else np.arange(da.sizes["band"])
        out_bands = [
            _resample_da(da.sel(band=band, drop=True), target_shape, method)
            for band in band_vals
        ]
        out = xr.concat(out_bands, dim="band")
        return out.assign_coords(band=band_vals).transpose("band", "y", "x")

    if da.ndim != 2:
        raise ValueError(
            f"_resample_da expects 2D or banded 3D DataArray, got dims={da.dims}"
        )

    src = da.values.astype(np.float64)
    h_out, w_out = target_shape

    if src.shape == target_shape:
        return da

    if method == "area":
        # Block-mean downsampling
        factor_y = max(1, src.shape[0] // h_out)
        factor_x = max(1, src.shape[1] // w_out)
        smoothed = uniform_filter(src, size=(factor_y, factor_x), mode="nearest")
        # Then bilinear resize to exact shape
        zoom_y = h_out / smoothed.shape[0]
        zoom_x = w_out / smoothed.shape[1]
        out = zoom(smoothed, (zoom_y, zoom_x), order=1)
    elif method == "nearest":
        zoom_y = h_out / src.shape[0]
        zoom_x = w_out / src.shape[1]
        out = zoom(src, (zoom_y, zoom_x), order=0)
    else:  # bilinear
        zoom_y = h_out / src.shape[0]
        zoom_x = w_out / src.shape[1]
        out = zoom(src, (zoom_y, zoom_x), order=1)

    # Trim/pad to exact target shape (zoom can be off by 1)
    out = out[:h_out, :w_out]
    if out.shape != target_shape:
        padded = np.zeros(target_shape, dtype=out.dtype)
        h = min(out.shape[0], h_out)
        w = min(out.shape[1], w_out)
        padded[:h, :w] = out[:h, :w]
        out = padded

    return xr.DataArray(
        out.astype(np.float32),
        dims=["y", "x"],
        attrs=da.attrs,
    )


def _resample_cloud_mask(
    mask: xr.DataArray,
    target_shape: tuple[int, int],
) -> xr.DataArray:
    """Conservative OR resampling for cloud masks.

    If any source pixel within a target-pixel footprint is True (cloudy),
    the target pixel is True.
    """
    src = mask.values.astype(np.float64)
    h_out, w_out = target_shape

    if src.shape == target_shape:
        return mask

    # Use max-pooling: zoom with nearest to upscale the "any-cloud" signal
    # First apply block-max via uniform_filter > 0 trick
    factor_y = max(1, src.shape[0] // h_out)
    factor_x = max(1, src.shape[1] // w_out)
    # Block-max: dilate by box, then check > 0
    from scipy.ndimage import maximum_filter
    dilated = maximum_filter(src, size=(factor_y, factor_x), mode="nearest")
    zoom_y = h_out / dilated.shape[0]
    zoom_x = w_out / dilated.shape[1]
    out = zoom(dilated, (zoom_y, zoom_x), order=0)
    out = out[:h_out, :w_out]

    if out.shape != target_shape:
        padded = np.ones(target_shape, dtype=np.float64)
        h = min(out.shape[0], h_out)
        w = min(out.shape[1], w_out)
        padded[:h, :w] = out[:h, :w]
        out = padded

    return xr.DataArray(
        out > 0.5,
        dims=["y", "x"],
        attrs=mask.attrs,
    )


def _resample_geometry(
    geom: GeometryAngles,
    target_shape: tuple[int, int],
) -> GeometryAngles:
    """Resample all geometry angles to target shape via bilinear interpolation."""
    return GeometryAngles(
        sza=_resample_da(geom.sza, target_shape, "bilinear"),
        saa=_resample_da(geom.saa, target_shape, "bilinear"),
        vza=_resample_da(geom.vza, target_shape, "bilinear"),
        vaa=_resample_da(geom.vaa, target_shape, "bilinear"),
    )


def _resample_atmo_state(
    atmo: AtmosphericState,
    target_shape: tuple[int, int],
) -> AtmosphericState:
    """Resample all atmospheric state fields to target shape."""
    return AtmosphericState(
        aot=_resample_da(atmo.aot, target_shape, "bilinear"),
        tcwv=_resample_da(atmo.tcwv, target_shape, "bilinear"),
        tco3=_resample_da(atmo.tco3, target_shape, "bilinear"),
        aot_unc=_resample_da(atmo.aot_unc, target_shape, "bilinear"),
        tcwv_unc=_resample_da(atmo.tcwv_unc, target_shape, "bilinear"),
        tco3_unc=_resample_da(atmo.tco3_unc, target_shape, "bilinear"),
        elevation=_resample_da(atmo.elevation, target_shape, "bilinear"),
    )


def _resample_surface_prior(
    sp: SurfacePrior,
    target_shape: tuple[int, int],
) -> SurfacePrior:
    """Resample surface prior fields to target shape."""
    mask = sp.mask
    if "band" in mask.dims:
        mask = mask.any(dim="band")

    return SurfacePrior(
        boa=_resample_da(sp.boa, target_shape, "bilinear"),
        boa_unc=_resample_da(sp.boa_unc, target_shape, "bilinear"),
        kernels=sp.kernels,  # keep original (used for spectral, not spatial)
        mask=_resample_cloud_mask(mask, target_shape),
    )


# ── Public API ─────────────────────────────────────────────────────────

def assemble_grids(
    obs: ObservationBundle,
    atmo: AtmosphericState,
    surface: SurfacePrior,
    rt_model: Any,
    aux_resolution_m: float = 500.0,
    aerosol_resolution_m: float = 1000.0,
) -> SolverInputBundle:
    """Resample and align all upstream outputs to solver grids.

    This is the M4 module callable.  It:

    1. Selects aerosol-sensitive bands (400-520 nm) from the sensor config.
    2. Resamples TOA, geometry, cloud mask to the *aux* resolution.
    3. Resamples atmospheric state and surface prior to the same grid.
    4. Returns a validated ``SolverInputBundle``.

    Args:
        obs: M1 output (satellite observation at native resolution).
        atmo: M2 output (atmospheric prior, possibly at coarse resolution).
        surface: M3 output (surface prior).
        rt_model: RT model backend (passed through, not resampled).
        aux_resolution_m: Target resolution for solver inputs (default 500 m).
        aerosol_resolution_m: Target resolution for AOT/TCWV retrieval (default 1000 m).

    Returns:
        SolverInputBundle ready for the solver.
    """
    # 1. Select solver bands by wavelength (aerosol-sensitive: 400-520 nm)
    bands = obs.sensor_config.select_bands_in_range(400.0, 520.0)
    if not bands:
        # Fallback: use first two bands
        bands = list(obs.sensor_config.bands[:2])
    logger.info(f"Selected {len(bands)} solver bands: {[b.name for b in bands]}")

    # 2. Determine target shape from the TOA native shape + resolution
    #    We need a representative native resolution to compute the scale factor.
    #    Use the first solver band's resolution (or fallback to 10 m).
    native_res = bands[0].resolution if bands else 10.0
    first_var = list(obs.toa.data_vars)[0]
    native_shape = obs.toa[first_var].shape  # (y, x)
    target_shape = _compute_target_shape(native_shape, native_res, aux_resolution_m)
    logger.info(
        f"Resampling from {native_shape} @ {native_res}m "
        f"to {target_shape} @ {aux_resolution_m}m"
    )

    # 3. Resample TOA for solver bands into (band, y, x) DataArray
    band_names = [b.name for b in bands]
    toa_arrays = []
    for bn in band_names:
        if bn in obs.toa.data_vars:
            resampled = _resample_da(obs.toa[bn], target_shape, "area")
            toa_arrays.append(resampled)
    if toa_arrays:
        toa_da = xr.concat(toa_arrays, dim="band")
    else:
        # Fallback: use first available variable
        first = list(obs.toa.data_vars)[0]
        toa_da = _resample_da(obs.toa[first], target_shape, "area").expand_dims("band")

    # 4. Resample geometry, cloud mask, atmo, surface
    geometry = _resample_geometry(obs.geometry, target_shape)
    cloud_mask = _resample_cloud_mask(obs.cloud_mask, target_shape)
    atmo_resampled = _resample_atmo_state(atmo, target_shape)
    surface_resampled = _resample_surface_prior(surface, target_shape)

    return SolverInputBundle(
        toa=toa_da,
        geometry=geometry,
        cloud_mask=cloud_mask,
        sensor_config=obs.sensor_config,
        bands=bands,
        atmo_prior=atmo_resampled,
        surface_prior=surface_resampled,
        rt_model=rt_model,
        aux_resolution_m=aux_resolution_m,
        aerosol_resolution_m=aerosol_resolution_m,
    )
