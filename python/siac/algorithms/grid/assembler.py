"""
M4 Grid Assembler — resample and align all upstream outputs to solver grids.

Produces a ``SolverInputBundle`` where every raster field lives on the same
spatial grid (the aerosol retrieval grid, typically 120 m).

See PLANS.md §4.4 for the full specification.
"""

from __future__ import annotations

import logging
from contextlib import suppress
from typing import Any, TypeAlias, cast

import numpy as np
import xarray as xr
from numpy import typing as npt
from scipy.ndimage import (
    binary_dilation,
    maximum_filter,
    median_filter,
    minimum_filter,
    sobel,
    uniform_filter,
    zoom,
)

from siac.runtime import (
    AtmosphericState,
    GeometryAngles,
    ObservationBundle,
    SolverInputBundle,
    SurfacePrior,
)
from siac.runtime.models import copy_spatial_metadata_like

logger = logging.getLogger(__name__)

BoolArray: TypeAlias = npt.NDArray[np.bool_]
Float32Array: TypeAlias = npt.NDArray[np.float32]
Float64Array: TypeAlias = npt.NDArray[np.float64]


# ── Internal resampling helpers ────────────────────────────────────────


def _template_coords(template: xr.DataArray | None) -> dict[str, xr.DataArray]:
    if template is None:
        return {}
    return {
        dim: template.coords[dim]
        for dim in template.dims
        if dim in template.coords
    }


def _ensure_template_transform(template: xr.DataArray) -> xr.DataArray:
    try:
        transform = template.rio.transform(recalc=True)
    except Exception:
        return template
    return template.rio.write_transform(transform)


def _build_target_template(
    bounds: tuple[float, float, float, float],
    crs: str,
    resolution: float,
) -> xr.DataArray:
    import rioxarray  # noqa: F401

    xmin, ymin, xmax, ymax = bounds
    width = max(1, int(np.ceil((xmax - xmin) / resolution)))
    height = max(1, int(np.ceil((ymax - ymin) / resolution)))
    x = xmin + (np.arange(width, dtype=np.float64) + 0.5) * resolution
    y = ymax - (np.arange(height, dtype=np.float64) + 0.5) * resolution
    template = xr.DataArray(
        np.full((height, width), np.nan, dtype=np.float32),
        dims=("y", "x"),
        coords={"x": x, "y": y},
    )
    template = template.rio.set_spatial_dims(x_dim="x", y_dim="y")
    template = template.rio.write_crs(crs)
    return _ensure_template_transform(template)


def _resample_da(
    da: xr.DataArray,
    target_shape: tuple[int, int],
    method: str = "bilinear",
    *,
    template: xr.DataArray | None = None,
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
            _resample_da(da.sel(band=band, drop=True), target_shape, method, template=template)
            for band in band_vals
        ]
        out = xr.concat(out_bands, dim="band")
        spatial_dims = tuple(template.dims) if template is not None else ("y", "x")
        out = out.assign_coords(band=band_vals).transpose("band", *spatial_dims)
        return copy_spatial_metadata_like(out, template) if template is not None else out

    if da.ndim != 2:
        raise ValueError(
            f"_resample_da expects 2D or banded 3D DataArray, got dims={da.dims}"
        )

    src = da.values.astype(np.float64)
    h_out, w_out = target_shape

    if src.shape == target_shape:
        return copy_spatial_metadata_like(da, template) if template is not None else da

    if method == "area":
        gdal_average = _resample_da_gdal_average(da, template)
        if gdal_average is not None:
            return gdal_average
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
        padded: np.ndarray[Any, Any] = np.zeros(target_shape, dtype=out.dtype)
        h = min(out.shape[0], h_out)
        w = min(out.shape[1], w_out)
        padded[:h, :w] = out[:h, :w]
        out = padded

    out_da = xr.DataArray(
        out.astype(np.float32),
        dims=tuple(template.dims) if template is not None else ("y", "x"),
        coords=_template_coords(template),
        attrs=da.attrs,
    )
    return copy_spatial_metadata_like(out_da, template) if template is not None else out_da


def _resample_da_gdal_average(
    da: xr.DataArray,
    template: xr.DataArray | None,
) -> xr.DataArray | None:
    """Use GDAL-backed average resampling when the source is georeferenced."""
    if template is None:
        return None

    try:
        import rioxarray  # noqa: F401
        from rasterio.enums import Resampling as RasterioResampling
    except Exception:
        logger.debug("rioxarray/rasterio not available; skipping GDAL average resampling.")
        return None

    try:
        x_dim = template.rio.x_dim
        y_dim = template.rio.y_dim
    except Exception:
        if "x" in template.dims and "y" in template.dims:
            x_dim, y_dim = "x", "y"
        else:
            logger.debug("Cannot determine spatial dims for GDAL resampling.")
            return None

    if x_dim not in da.dims or y_dim not in da.dims:
        return None
    if x_dim not in da.coords or y_dim not in da.coords:
        return None

    try:
        source = da.rio.set_spatial_dims(x_dim=x_dim, y_dim=y_dim)
    except Exception:
        logger.debug("Failed to set spatial dims for GDAL resampling.", exc_info=True)
        return None

    try:
        source_crs = source.rio.crs
    except Exception:
        source_crs = None

    try:
        target = _ensure_template_transform(template)
        target_crs = target.rio.crs
    except Exception:
        logger.debug("Failed to read target CRS/transform.", exc_info=True)
        return None

    if target_crs is None:
        return None
    if source_crs is None:
        source = source.rio.write_crs(target_crs)

    with suppress(Exception):
        source = source.rio.write_transform(source.rio.transform(recalc=True))

    try:
        out = source.rio.reproject_match(target, resampling=RasterioResampling.average)
    except Exception:
        logger.debug("GDAL reproject_match failed; falling back to scipy.", exc_info=True)
        return None

    out = out.astype(np.float32).assign_attrs(da.attrs)
    return copy_spatial_metadata_like(out, target)


def _resample_cloud_mask(
    mask: xr.DataArray,
    target_shape: tuple[int, int],
    *,
    template: xr.DataArray | None = None,
) -> xr.DataArray:
    """Conservative OR resampling for cloud masks.

    If any source pixel within a target-pixel footprint is True (cloudy),
    the target pixel is True.
    """
    src = mask.values.astype(np.float64)
    h_out, w_out = target_shape

    if src.shape == target_shape:
        return copy_spatial_metadata_like(mask, template) if template is not None else mask

    if src.ndim == 2 and h_out <= src.shape[0] and w_out <= src.shape[1]:
        # Match the center-based coarse-grid assignment used by the Rust field
        # remapper so masks and resampled fields describe the same footprint.
        out: Float64Array = np.zeros(target_shape, dtype=np.float64)
        for iy in range(src.shape[0]):
            dst_y = min(((2 * iy + 1) * h_out) // (2 * src.shape[0]), h_out - 1)
            for ix in range(src.shape[1]):
                if src[iy, ix] <= 0.5:
                    continue
                dst_x = min(((2 * ix + 1) * w_out) // (2 * src.shape[1]), w_out - 1)
                out[dst_y, dst_x] = 1.0
    else:
        out = cast(
            "Float64Array",
            zoom(src, (h_out / src.shape[0], w_out / src.shape[1]), order=0),
        )
        out = out[:h_out, :w_out]

    if out.shape != target_shape:
        padded: np.ndarray[Any, Any] = np.ones(target_shape, dtype=np.float64)
        h = min(out.shape[0], h_out)
        w = min(out.shape[1], w_out)
        padded[:h, :w] = out[:h, :w]
        out = padded

    out_mask = xr.DataArray(
        out > 0.5,
        dims=tuple(template.dims) if template is not None else ("y", "x"),
        coords=_template_coords(template),
        attrs=mask.attrs,
    )
    return copy_spatial_metadata_like(out_mask, template) if template is not None else out_mask


def _native_solver_band_stack(
    obs: ObservationBundle,
    band_names: list[str],
) -> xr.DataArray:
    available_names = [band_name for band_name in band_names if band_name in obs.toa.data_vars]
    arrays = [obs.toa[band_name] for band_name in available_names]
    if not arrays:
        first = list(obs.toa.data_vars)[0]
        return obs.toa[first].expand_dims("band").assign_coords(band=[first])
    return xr.concat(arrays, dim="band").assign_coords(band=available_names)


def _coerce_window_size(window_pixels_native: int) -> int:
    window = max(3, int(window_pixels_native))
    if window % 2 == 0:
        window += 1
    return window


def _dilate_mask(mask: np.ndarray, iterations: int) -> BoolArray:
    if iterations <= 0 or not np.any(mask):
        return cast("BoolArray", np.asarray(mask, dtype=bool))
    return cast("BoolArray", binary_dilation(mask, iterations=iterations))


def _robust_local_positive_z(
    metric: np.ndarray,
    valid_mask: np.ndarray,
    window: int,
) -> Float32Array:
    if not np.any(valid_mask):
        return cast("Float32Array", np.zeros_like(metric, dtype=np.float32))

    fill_value = float(np.median(metric[valid_mask]))
    filled = np.where(valid_mask, metric, fill_value).astype(np.float32, copy=False)
    local_background = median_filter(filled, size=window, mode="reflect")
    local_scale = median_filter(np.abs(filled - local_background), size=window, mode="reflect")
    z = (filled - local_background) / np.maximum(local_scale, 1.0e-3)
    return cast(
        "Float32Array",
        np.where(valid_mask, z, 0.0).astype(np.float32, copy=False),
    )


def _detect_sharp_transition_mask_native(
    toa_native: xr.DataArray,
    cloud_mask: xr.DataArray,
    filter_config: Any,
) -> xr.DataArray:
    band_count = toa_native.sizes.get("band", 1)
    if toa_native.ndim == 2:
        values = toa_native.values[np.newaxis, ...].astype(np.float32, copy=False)
    else:
        values = toa_native.values.astype(np.float32, copy=False)

    window = _coerce_window_size(getattr(filter_config, "window_pixels_native", 7))
    context_window = _coerce_window_size(
        getattr(filter_config, "context_window_pixels_native", max(15, (2 * window) + 1))
    )
    coherence_window = _coerce_window_size(
        getattr(filter_config, "coherence_window_pixels_native", window)
    )
    road_std_z_threshold = float(getattr(filter_config, "road_std_z_threshold_native", 4.0))
    road_coherence_threshold = float(
        getattr(filter_config, "road_coherence_threshold_native", 0.75)
    )
    road_std_floor = float(getattr(filter_config, "road_std_floor_native", 0.02))
    point_range_z_threshold = float(getattr(filter_config, "point_range_z_threshold_native", 4.0))
    point_outlier_fraction_max = float(
        getattr(filter_config, "point_outlier_fraction_max_native", 0.12)
    )
    point_range_floor = float(getattr(filter_config, "point_range_floor_native", 0.08))
    outlier_sigma = float(getattr(filter_config, "outlier_sigma_native", 2.5))
    outlier_floor = float(getattr(filter_config, "outlier_floor_native", 0.01))
    dilation_pixels = int(getattr(filter_config, "dilation_pixels", 1))
    cloud_buffer_pixels = int(getattr(filter_config, "cloud_buffer_pixels", 0))

    clear_mask = ~cloud_mask.values.astype(bool)
    if cloud_buffer_pixels > 0:
        clear_mask = ~_dilate_mask(~clear_mask, cloud_buffer_pixels)

    finite_valid = np.all(np.isfinite(values), axis=0)
    toa_valid = np.all((values > 0.0) & (values < 1.0), axis=0)
    valid_mask = clear_mask & finite_valid & toa_valid

    if not np.any(valid_mask):
        empty = xr.DataArray(
            np.zeros(valid_mask.shape, dtype=bool),
            dims=cloud_mask.dims,
            coords=cloud_mask.coords,
        )
        return copy_spatial_metadata_like(empty, cloud_mask)

    range_max = np.zeros(valid_mask.shape, dtype=np.float32)
    std_max = np.zeros(valid_mask.shape, dtype=np.float32)
    outlier_frac_max = np.zeros(valid_mask.shape, dtype=np.float32)
    coherence_max = np.zeros(valid_mask.shape, dtype=np.float32)

    for band_index in range(band_count):
        band_values = values[band_index]

        valid_band_values = band_values[valid_mask]
        fill_value = float(np.median(valid_band_values)) if valid_band_values.size else 0.0
        filled = np.where(valid_mask, band_values, fill_value)

        local_max = maximum_filter(filled, size=window, mode="reflect")
        local_min = minimum_filter(filled, size=window, mode="reflect")
        local_range = local_max - local_min

        local_mean = uniform_filter(filled, size=window, mode="reflect")
        local_mean_sq = uniform_filter(np.square(filled, dtype=np.float32), size=window, mode="reflect")
        local_std = np.sqrt(np.maximum(local_mean_sq - np.square(local_mean, dtype=np.float32), 0.0))

        outlier_hits = np.abs(filled - local_mean) >= np.maximum(outlier_sigma * local_std, outlier_floor)
        outlier_frac = uniform_filter(outlier_hits.astype(np.float32), size=window, mode="reflect")

        grad_x = sobel(filled, axis=1, mode="reflect")
        grad_y = sobel(filled, axis=0, mode="reflect")
        jxx = uniform_filter(np.square(grad_x, dtype=np.float32), size=coherence_window, mode="reflect")
        jyy = uniform_filter(np.square(grad_y, dtype=np.float32), size=coherence_window, mode="reflect")
        jxy = uniform_filter((grad_x * grad_y).astype(np.float32, copy=False), size=coherence_window, mode="reflect")
        coherence = np.sqrt(np.square(jxx - jyy) + 4.0 * np.square(jxy)) / np.maximum(jxx + jyy, 1.0e-6)

        range_max = np.maximum(range_max, local_range.astype(np.float32, copy=False))
        std_max = np.maximum(std_max, local_std.astype(np.float32, copy=False))
        outlier_frac_max = np.maximum(outlier_frac_max, outlier_frac.astype(np.float32, copy=False))
        coherence_max = np.maximum(coherence_max, coherence.astype(np.float32, copy=False))

    range_z = _robust_local_positive_z(range_max, valid_mask, context_window)
    std_z = _robust_local_positive_z(std_max, valid_mask, context_window)

    road_mask = (
        valid_mask
        & (std_z >= road_std_z_threshold)
        & (coherence_max >= road_coherence_threshold)
        & (std_max >= road_std_floor)
    )
    point_mask = (
        valid_mask
        & (range_z >= point_range_z_threshold)
        & (outlier_frac_max <= point_outlier_fraction_max)
        & (range_max >= point_range_floor)
    )
    detected = road_mask | point_mask
    detected = _dilate_mask(detected, dilation_pixels) & clear_mask
    detected_da = xr.DataArray(detected, dims=cloud_mask.dims, coords=cloud_mask.coords)
    return copy_spatial_metadata_like(detected_da, cloud_mask)


def _aggregate_native_exclusion_mask(
    native_mask: xr.DataArray,
    target_shape: tuple[int, int],
    *,
    template: xr.DataArray,
    fraction_threshold: float,
) -> xr.DataArray:
    threshold = float(fraction_threshold)
    if threshold <= 0.0:
        return _resample_cloud_mask(native_mask, target_shape, template=template)

    fraction = _resample_da(native_mask.astype(np.float32), target_shape, "area", template=template)
    out_mask = xr.DataArray(
        fraction.values >= threshold,
        dims=template.dims,
        coords=_template_coords(template),
        attrs=native_mask.attrs,
    )
    return copy_spatial_metadata_like(out_mask, template)


def _resample_geometry(
    geom: GeometryAngles,
    target_shape: tuple[int, int],
    *,
    template: xr.DataArray | None = None,
) -> GeometryAngles:
    """Resample all geometry angles to target shape via bilinear interpolation."""
    return GeometryAngles(
        sza=_resample_da(geom.sza, target_shape, "bilinear", template=template),
        saa=_resample_da(geom.saa, target_shape, "bilinear", template=template),
        vza=_resample_da(geom.vza, target_shape, "bilinear", template=template),
        vaa=_resample_da(geom.vaa, target_shape, "bilinear", template=template),
    )


def _resample_atmo_state(
    atmo: AtmosphericState,
    target_shape: tuple[int, int],
    *,
    template: xr.DataArray | None = None,
) -> AtmosphericState:
    """Resample all atmospheric state fields to target shape."""
    return AtmosphericState(
        aot=_resample_da(atmo.aot, target_shape, "bilinear", template=template),
        tcwv=_resample_da(atmo.tcwv, target_shape, "bilinear", template=template),
        tco3=_resample_da(atmo.tco3, target_shape, "bilinear", template=template),
        aot_unc=_resample_da(atmo.aot_unc, target_shape, "bilinear", template=template),
        tcwv_unc=_resample_da(atmo.tcwv_unc, target_shape, "bilinear", template=template),
        tco3_unc=_resample_da(atmo.tco3_unc, target_shape, "bilinear", template=template),
        elevation=_resample_da(atmo.elevation, target_shape, "bilinear", template=template),
    )


def _resample_surface_prior(
    sp: SurfacePrior,
    target_shape: tuple[int, int],
    *,
    template: xr.DataArray | None = None,
) -> SurfacePrior:
    """Resample surface prior fields to target shape."""
    mask = sp.mask
    if "band" in mask.dims:
        mask = mask.all(dim="band")

    return SurfacePrior(
        boa=_resample_da(sp.boa, target_shape, "area", template=template),
        boa_unc=_resample_da(sp.boa_unc, target_shape, "area", template=template),
        kernels=sp.kernels,  # keep original (used for spectral, not spatial)
        mask=_resample_cloud_mask(mask, target_shape, template=template),
        monthly_composites=sp.monthly_composites,
    )


def _align_surface_prior_to_bands(
    sp: SurfacePrior,
    band_names: list[str],
) -> SurfacePrior:
    """Reorder/select a banded surface prior to match the solver-band order."""
    boa = sp.boa
    boa_unc = sp.boa_unc
    mask = sp.mask

    if "band" in boa.dims and "band" in boa.coords:
        available = [str(value) for value in np.asarray(boa.coords["band"].values).tolist()]
        # Only align if the coords look like actual band names (strings matching
        # solver band names).  When coords are numeric indices (e.g. [1, 2]),
        # the surface prior is positional — skip alignment.
        if any(name in available for name in band_names):
            missing = [name for name in band_names if name not in available]
            if missing:
                raise ValueError(f"Surface prior is missing solver bands: {missing}")
            boa = boa.sel(band=band_names)
            if "band" in boa_unc.dims and "band" in boa_unc.coords:
                boa_unc = boa_unc.sel(band=band_names)
            if "band" in mask.dims and "band" in mask.coords:
                mask = mask.sel(band=band_names)

    return SurfacePrior(
        boa=boa,
        boa_unc=boa_unc,
        kernels=sp.kernels,
        mask=mask,
        monthly_composites=sp.monthly_composites,
    )


# ── Public API ─────────────────────────────────────────────────────────

def assemble_grids(
    obs: ObservationBundle,
    atmo: AtmosphericState,
    surface: SurfacePrior,
    rt_model: Any,
    aux_resolution_m: float = 500.0,
    aerosol_resolution_m: float = 120.0,
    sharp_transition_filter: Any | None = None,
) -> SolverInputBundle:
    """Resample and align all upstream outputs to solver grids.

    This is the M4 module callable.  It:

    1. Selects default aerosol-retrieval bands from the sensor config.
    2. Resamples TOA, geometry, cloud mask to the aerosol retrieval resolution.
    3. Resamples atmospheric state and surface prior to the same grid.
    4. Returns a validated ``SolverInputBundle``.

    Args:
        obs: M1 output (satellite observation at native resolution).
        atmo: M2 output (atmospheric prior, possibly at coarse resolution).
        surface: M3 output (surface prior).
        rt_model: RT model backend (passed through, not resampled).
        aux_resolution_m: Legacy metadata field retained for compatibility.
        aerosol_resolution_m: Target resolution for AOT/TCWV retrieval (default 120 m).

    Returns:
        SolverInputBundle ready for the solver.
    """
    # 1. Select solver bands using the sensor defaults.
    bands = obs.sensor_config.default_aerosol_solver_bands()
    logger.info(f"Selected {len(bands)} solver bands: {[b.name for b in bands]}")

    # 2. Determine the solver grid from the configured aerosol retrieval
    #    resolution and scene bounds. If callers still set only the legacy
    #    aux-resolution knob, preserve that behavior.
    resolved_aerosol_resolution = float(aerosol_resolution_m)
    if aerosol_resolution_m == 120.0 and aux_resolution_m != 500.0:
        resolved_aerosol_resolution = float(aux_resolution_m)
        logger.info(
            "assemble_grids: using aux_resolution_m=%s as the solver grid for legacy compatibility",
            aux_resolution_m,
        )
    first_var = list(obs.toa.data_vars)[0]
    native_shape = obs.toa[first_var].shape  # (y, x)
    target_template = _build_target_template(obs.bounds, obs.crs, resolved_aerosol_resolution)
    target_shape = target_template.shape
    logger.info(
        f"Resampling from {native_shape} to {target_shape} @ {resolved_aerosol_resolution}m"
    )

    # 3. Detect native sharp transitions before solver-grid averaging.
    band_names = [b.name for b in bands]
    sharp_transition_mask: xr.DataArray | None = None
    if sharp_transition_filter is not None and bool(getattr(sharp_transition_filter, "enabled", False)):
        native_solver_toa = _native_solver_band_stack(obs, band_names)
        native_mask = _detect_sharp_transition_mask_native(
            native_solver_toa,
            obs.cloud_mask,
            sharp_transition_filter,
        )
        sharp_transition_mask = _aggregate_native_exclusion_mask(
            native_mask,
            target_shape,
            template=target_template,
            fraction_threshold=float(
                getattr(sharp_transition_filter, "solver_cell_fraction_threshold", 0.0)
            ),
        )

    # 4. Resample TOA for solver bands into (band, y, x) DataArray
    toa_arrays = []
    for bn in band_names:
        if bn in obs.toa.data_vars:
            resampled = _resample_da(obs.toa[bn], target_shape, "area", template=target_template)
            toa_arrays.append(resampled)
    if toa_arrays:
        toa_da = xr.concat(toa_arrays, dim="band")
        toa_da = toa_da.assign_coords(band=band_names)
    else:
        # Fallback: use first available variable
        first = list(obs.toa.data_vars)[0]
        toa_da = _resample_da(obs.toa[first], target_shape, "area", template=target_template).expand_dims("band")
        toa_da = toa_da.assign_coords(band=[first])
    toa_da = copy_spatial_metadata_like(toa_da, target_template)

    # 5. Resample geometry, cloud mask, atmo, surface
    geometry = _resample_geometry(obs.geometry, target_shape, template=target_template)
    cloud_mask = _resample_cloud_mask(obs.cloud_mask, target_shape, template=target_template)
    atmo_resampled = _resample_atmo_state(atmo, target_shape, template=target_template)
    surface_aligned = _align_surface_prior_to_bands(surface, band_names)
    surface_resampled = _resample_surface_prior(surface_aligned, target_shape, template=target_template)

    return SolverInputBundle(
        toa=toa_da,
        geometry=geometry,
        cloud_mask=cloud_mask,
        sharp_transition_mask=sharp_transition_mask,
        sensor_config=obs.sensor_config,
        bands=bands,
        atmo_prior=atmo_resampled,
        surface_prior=surface_resampled,
        rt_model=rt_model,
        aux_resolution_m=aux_resolution_m,
        aerosol_resolution_m=resolved_aerosol_resolution,
    )
