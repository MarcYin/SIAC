"""
M4 Grid Assembler — resample and align all upstream outputs to solver grids.

Produces a ``SolverInputBundle`` where every raster field lives on the same
spatial grid (the aerosol retrieval grid, typically 120 m).

See PLANS.md §4.4 for the full specification.
"""

from __future__ import annotations

import hashlib
import logging
import warnings
from typing import TYPE_CHECKING, Any, TypeAlias, cast

import cv2
import numpy as np
import xarray as xr
from numpy import typing as npt
from scipy.ndimage import (
    binary_dilation,
    uniform_filter,
    zoom,
)

from siac.adapters.data.water_mask import load_water_mask_subset
from siac.geo.resample import shares_template_grid
from siac.runtime import (
    AtmosphericState,
    GeometryAngles,
    ObservationBundle,
    SolverInputBundle,
    SurfacePrior,
)
from siac.runtime.models import copy_spatial_metadata_like

if TYPE_CHECKING:
    from pathlib import Path

    from siac.algorithms.grid.toa_psf import ToaPsfConfig

logger = logging.getLogger(__name__)

BoolArray: TypeAlias = npt.NDArray[np.bool_]


# Atmospheric providers may set this attribute when they have a compact,
# stable source identifier (for example, a MAIAC granule or matching L2A
# asset).  Fields without it are fingerprinted below so cache correctness
# does not depend on a provider remembering to supply one.
_ATMO_CACHE_ID_ATTR = "siac_atmo_cache_identity"


# ── Internal resampling helpers ────────────────────────────────────────


def _template_coords(template: xr.DataArray | None) -> dict[str, xr.DataArray]:
    if template is None:
        return {}
    return {dim: template.coords[dim] for dim in template.dims if dim in template.coords}


def _ensure_template_transform(template: xr.DataArray) -> xr.DataArray:
    try:
        transform = template.rio.transform(recalc=True)
    except (AttributeError, RuntimeError, ValueError):
        return template
    return template.rio.write_transform(transform)


def _spatial_dims(data: xr.DataArray) -> tuple[str, str] | None:
    try:
        return data.rio.x_dim, data.rio.y_dim
    except (AttributeError, RuntimeError, ValueError):
        pass
    if "x" in data.dims and "y" in data.dims:
        return "x", "y"
    if "longitude" in data.dims and "latitude" in data.dims:
        return "longitude", "latitude"
    return None


def _is_monotonic_axis(coord: xr.DataArray) -> bool:
    values = np.asarray(coord.values, dtype=np.float64)
    if values.ndim != 1 or values.size < 2:
        return True
    diffs = np.diff(values)
    if not np.all(np.isfinite(diffs)):
        return False
    return bool(np.all(diffs > 0.0) or np.all(diffs < 0.0))


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


def _pil_resize(src: np.ndarray, w_out: int, h_out: int, method: str) -> np.ndarray:
    """Fast float32 resize using PIL (7× faster than scipy.ndimage.zoom)."""
    from PIL import Image

    resampling = Image.NEAREST if method == "nearest" else Image.BILINEAR
    return np.array(
        Image.fromarray(src.astype(np.float32), mode="F").resize(
            (w_out, h_out),
            resampling,
        ),
        dtype=np.float32,
    )


def _resample_da(
    da: xr.DataArray,
    target_shape: tuple[int, int],
    method: str = "bilinear",
    *,
    template: xr.DataArray | None = None,
) -> xr.DataArray:
    """Resample a 2-D DataArray to *target_shape*.

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
        raise ValueError(f"_resample_da expects 2D or banded 3D DataArray, got dims={da.dims}")

    src = np.asarray(da.values, dtype=np.float32)
    h_out, w_out = target_shape

    if template is not None and shares_template_grid(da, template):
        return copy_spatial_metadata_like(da, template)
    if template is None and src.shape == target_shape:
        return copy_spatial_metadata_like(da, template) if template is not None else da

    if method == "area":
        gdal_average = _resample_da_gdal(da, template, method="area")
        if gdal_average is not None:
            return gdal_average
        # Block-mean downsampling
        factor_y = max(1, src.shape[0] // h_out)
        factor_x = max(1, src.shape[1] // w_out)
        smoothed = uniform_filter(src, size=(factor_y, factor_x), mode="nearest")
        # Then bilinear resize to exact shape
        out = _pil_resize(smoothed, w_out, h_out, "bilinear")
    elif method == "nearest":
        gdal_nearest = _resample_da_gdal(da, template, method="nearest")
        if gdal_nearest is not None:
            return gdal_nearest
        out = _pil_resize(src, w_out, h_out, "nearest")
    else:  # bilinear
        gdal_bilinear = _resample_da_gdal(da, template, method="bilinear")
        if gdal_bilinear is not None:
            return gdal_bilinear
        out = _pil_resize(src, w_out, h_out, "bilinear")

    # Trim/pad to exact target shape (PIL gives exact size, but guard anyway)
    out = out[:h_out, :w_out]
    if out.shape != target_shape:
        padded: np.ndarray[Any, Any] = np.zeros(target_shape, dtype=out.dtype)
        h = min(out.shape[0], h_out)
        w = min(out.shape[1], w_out)
        padded[:h, :w] = out[:h, :w]
        out = padded

    out_da = xr.DataArray(
        np.asarray(out, dtype=np.float32),
        dims=tuple(template.dims) if template is not None else ("y", "x"),
        coords=_template_coords(template),
        attrs=da.attrs,
    )
    return copy_spatial_metadata_like(out_da, template) if template is not None else out_da


def _resample_da_gdal(
    da: xr.DataArray,
    template: xr.DataArray | None,
    *,
    method: str,
) -> xr.DataArray | None:
    """Use GDAL-backed resampling when the source is georeferenced."""
    if template is None:
        return None

    try:
        import rioxarray  # noqa: F401
        from rasterio.enums import Resampling as RasterioResampling
        from rasterio.errors import NotGeoreferencedWarning
    except (ImportError, ModuleNotFoundError):
        logger.debug("rioxarray/rasterio not available; skipping GDAL resampling.")
        return None

    resampling_lookup = {
        "area": RasterioResampling.average,
        "bilinear": RasterioResampling.bilinear,
        "nearest": RasterioResampling.nearest,
    }
    resampling = resampling_lookup.get(method)
    if resampling is None:
        return None

    source_dims = _spatial_dims(da)
    target_dims = _spatial_dims(template)
    if source_dims is None or target_dims is None:
        logger.debug("Cannot determine spatial dims for GDAL resampling.")
        return None
    source_x_dim, source_y_dim = source_dims
    target_x_dim, target_y_dim = target_dims

    if source_x_dim not in da.dims or source_y_dim not in da.dims:
        return None
    if source_x_dim not in da.coords or source_y_dim not in da.coords:
        return None
    if target_x_dim not in template.dims or target_y_dim not in template.dims:
        return None
    if target_x_dim not in template.coords or target_y_dim not in template.coords:
        return None
    if not _is_monotonic_axis(da.coords[source_x_dim]) or not _is_monotonic_axis(
        da.coords[source_y_dim]
    ):
        logger.debug("Skipping GDAL resampling for non-monotonic source coordinates.")
        return None

    # Validate georeferencing in a single warnings context to avoid repeated
    # context-manager overhead on each call.
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("error", category=NotGeoreferencedWarning)

            source = da.rio.set_spatial_dims(x_dim=source_x_dim, y_dim=source_y_dim)
            source_crs = source.rio.crs

            target = _ensure_template_transform(template)
            target_crs = target.rio.crs

            if target_crs is None or source_crs is None:
                return None

            with __import__("contextlib").suppress(AttributeError, RuntimeError, ValueError):
                source = source.rio.write_transform(source.rio.transform(recalc=True))

            out = source.rio.reproject_match(target, resampling=resampling)
    except NotGeoreferencedWarning:
        logger.debug("Skipping GDAL resampling: missing georeferencing metadata.")
        return None
    except (RuntimeError, ValueError, AttributeError) as exc:
        logger.debug(
            "GDAL reproject_match failed (%s); falling back to scipy.",
            type(exc).__name__,
            exc_info=True,
        )
        return None

    out = out.astype(np.float32).assign_attrs(da.attrs)
    return copy_spatial_metadata_like(out, target)


def _resample_cloud_mask(
    mask: xr.DataArray,
    target_shape: tuple[int, int],
    *,
    template: xr.DataArray | None = None,
    assume_aligned_native_grid: bool = False,
) -> xr.DataArray:
    """Conservative OR resampling for cloud masks.

    If any source pixel within a target-pixel footprint is True (cloudy),
    the target pixel is True.
    """
    src = np.asarray(mask.values, dtype=bool)
    h_out, w_out = target_shape

    if src.shape == target_shape:
        return copy_spatial_metadata_like(mask, template) if template is not None else mask

    use_shape_only = (
        template is None
        or assume_aligned_native_grid
        or _shape_only_mask_remap_is_safe(mask, template)
    )

    if use_shape_only and src.ndim == 2 and h_out <= src.shape[0] and w_out <= src.shape[1]:
        # Match the center-based coarse-grid assignment used by the Rust field
        # remapper so masks and resampled fields describe the same footprint.
        out = np.zeros(target_shape, dtype=bool)
        src_y, src_x = np.nonzero(src)
        if src_y.size:
            dst_y = np.minimum(((2 * src_y + 1) * h_out) // (2 * src.shape[0]), h_out - 1).astype(
                np.intp,
                copy=False,
            )
            dst_x = np.minimum(((2 * src_x + 1) * w_out) // (2 * src.shape[1]), w_out - 1).astype(
                np.intp,
                copy=False,
            )
            out[dst_y, dst_x] = True
    else:
        geospatial = _resample_mask_geospatial(mask, target_shape, template=template)
        if geospatial is not None:
            return geospatial
        out = cast(
            "BoolArray",
            zoom(
                src.astype(np.uint8, copy=False),
                (h_out / src.shape[0], w_out / src.shape[1]),
                order=0,
            )
            > 0,
        )
        out = out[:h_out, :w_out]

    if out.shape != target_shape:
        padded: np.ndarray[Any, Any] = np.ones(target_shape, dtype=bool)
        h = min(out.shape[0], h_out)
        w = min(out.shape[1], w_out)
        padded[:h, :w] = out[:h, :w]
        out = padded

    out_mask = xr.DataArray(
        out,
        dims=tuple(template.dims) if template is not None else ("y", "x"),
        coords=_template_coords(template),
        attrs=mask.attrs,
    )
    return copy_spatial_metadata_like(out_mask, template) if template is not None else out_mask


def _shape_only_mask_remap_is_safe(
    mask: xr.DataArray,
    template: xr.DataArray | None,
) -> bool:
    """Return True when a shape-only remap is safe for *mask* -> *template*.

    The fast shape-only remapper is only valid when the source mask already
    lives on the same CRS and geographic footprint as the target grid. This
    helper enforces that contract for any georeferenced source/target pair.
    """
    if template is None:
        return False
    if shares_template_grid(mask, template):
        return True

    source_dims = _spatial_dims(mask)
    target_dims = _spatial_dims(template)
    if source_dims is None or target_dims is None:
        return False

    try:
        source = mask.rio.set_spatial_dims(x_dim=source_dims[0], y_dim=source_dims[1])
        target = template.rio.set_spatial_dims(x_dim=target_dims[0], y_dim=target_dims[1])
        source_crs = source.rio.crs
        target_crs = target.rio.crs
        if source_crs is None or target_crs is None or str(source_crs) != str(target_crs):
            return False
        source_bounds = tuple(float(v) for v in source.rio.bounds())
        target_bounds = tuple(float(v) for v in target.rio.bounds())
        source_res = source.rio.resolution()
        target_res = target.rio.resolution()
    except (AttributeError, RuntimeError, ValueError, TypeError):
        return False

    tolerance = max(
        abs(float(source_res[0])),
        abs(float(source_res[1])),
        abs(float(target_res[0])),
        abs(float(target_res[1])),
        1e-6,
    )
    return all(
        abs(src_bound - dst_bound) <= tolerance
        for src_bound, dst_bound in zip(source_bounds, target_bounds)
    )


def _resample_mask_geospatial(
    mask: xr.DataArray,
    target_shape: tuple[int, int],
    *,
    template: xr.DataArray | None,
    fraction_threshold: float = 0.0,
) -> xr.DataArray | None:
    """Geospatially reproject a mask to *template* when metadata allows it."""
    if template is None:
        return None
    fraction = _resample_da(mask.astype(np.float32), target_shape, "area", template=template)
    values = np.asarray(fraction.values, dtype=np.float32)
    threshold = float(fraction_threshold)
    out = values > 0.0 if threshold <= 0.0 else values >= threshold
    out_mask = xr.DataArray(
        out,
        dims=template.dims,
        coords=_template_coords(template),
        attrs=mask.attrs,
    )
    return copy_spatial_metadata_like(out_mask, template)


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


def _dilate_mask_ellipse(mask: np.ndarray, radius: int) -> BoolArray:
    if radius <= 0 or not np.any(mask):
        return cast("BoolArray", np.asarray(mask, dtype=bool))
    size = (2 * int(radius)) + 1
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (size, size))
    dilated = cv2.dilate(mask.astype(np.uint8, copy=False) * 255, kernel, iterations=1)
    return cast("BoolArray", dilated.astype(bool, copy=False))


def _sharp_transition_proxy_uint8(
    values: np.ndarray, valid_mask: np.ndarray
) -> npt.NDArray[np.uint8]:
    clipped = np.clip(
        np.nan_to_num(values, nan=0.0, posinf=1.0, neginf=0.0),
        0.0,
        1.0,
    ).astype(np.float32, copy=False)
    proxy = clipped.mean(axis=0, dtype=np.float32)
    fill_value = float(np.median(proxy[valid_mask])) if np.any(valid_mask) else 0.0
    proxy = np.where(valid_mask, proxy, fill_value)
    return np.clip(np.rint(proxy * 255.0), 0.0, 255.0).astype(np.uint8, copy=False)


def _detect_sharp_transition_mask_native(
    toa_native: xr.DataArray,
    cloud_mask: xr.DataArray,
    filter_config: Any,
) -> xr.DataArray:
    spatial_template = toa_native if toa_native.ndim == 2 else toa_native.isel(band=0, drop=True)
    if toa_native.ndim == 2:
        values = toa_native.values[np.newaxis, ...].astype(np.float32, copy=False)
    else:
        values = toa_native.values.astype(np.float32, copy=False)

    blur_kernel = _coerce_window_size(
        getattr(
            filter_config,
            "blur_kernel_pixels_native",
            getattr(
                filter_config,
                "context_window_pixels_native",
                getattr(filter_config, "window_pixels_native", 31),
            ),
        )
    )
    residual_threshold = int(getattr(filter_config, "residual_threshold_uint8", 12))
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
        return copy_spatial_metadata_like(empty, spatial_template)

    proxy_uint8 = _sharp_transition_proxy_uint8(values, valid_mask)
    background = cv2.blur(proxy_uint8, (blur_kernel, blur_kernel))
    residual = cv2.absdiff(proxy_uint8, background)
    # Thresholds are defined in uint8 residual space to keep the detector
    # fast and reproducible across scenes.
    _, thresholded = cv2.threshold(
        residual,
        residual_threshold,
        255,
        cv2.THRESH_BINARY,
    )
    detected = thresholded.astype(bool, copy=False) & valid_mask
    detected = _dilate_mask_ellipse(detected, dilation_pixels) & clear_mask
    detected_da = xr.DataArray(detected, dims=cloud_mask.dims, coords=cloud_mask.coords)
    return copy_spatial_metadata_like(detected_da, spatial_template)


def _aggregate_native_exclusion_mask(
    native_mask: xr.DataArray,
    target_shape: tuple[int, int],
    *,
    template: xr.DataArray,
    fraction_threshold: float,
    assume_aligned_native_grid: bool = True,
) -> xr.DataArray:
    threshold = float(fraction_threshold)
    if threshold <= 0.0:
        return _resample_cloud_mask(
            native_mask,
            target_shape,
            template=template,
            assume_aligned_native_grid=assume_aligned_native_grid,
        )

    fraction = _resample_da(native_mask.astype(np.float32), target_shape, "area", template=template)
    out_mask = xr.DataArray(
        fraction.values >= threshold,
        dims=template.dims,
        coords=_template_coords(template),
        attrs=native_mask.attrs,
    )
    return copy_spatial_metadata_like(out_mask, template)


def _resample_external_exclusion_mask(
    mask: xr.DataArray,
    target_shape: tuple[int, int],
    *,
    template: xr.DataArray,
    fraction_threshold: float = 0.0,
) -> xr.DataArray:
    """Reproject an externally sourced exclusion mask onto the solver grid.

    Native SIAC masks such as cloud and sharp-transition layers are already on
    the observation grid, so shape-based conservative remapping is sufficient.
    External masks such as the Zenodo land/water product can arrive in a
    different CRS/grid (e.g. EPSG:4326 versus UTM), so they must be
    reprojected geospatially before any boolean thresholding.
    """
    if shares_template_grid(mask, template):
        return _aggregate_native_exclusion_mask(
            mask,
            target_shape,
            template=template,
            fraction_threshold=fraction_threshold,
        )

    fraction = _resample_da(mask.astype(np.float32), target_shape, "area", template=template)
    values = np.asarray(fraction.values, dtype=np.float32)
    threshold = float(fraction_threshold)
    out = values > 0.0 if threshold <= 0.0 else values >= threshold
    out_mask = xr.DataArray(
        out,
        dims=template.dims,
        coords=_template_coords(template),
        attrs=mask.attrs,
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


def _update_atmo_cache_hash(digest: Any, value: object) -> None:
    """Add an array-like value to an atmospheric-field cache fingerprint."""
    array = np.ascontiguousarray(np.asarray(value))
    digest.update(str(array.dtype).encode("utf-8"))
    digest.update(repr(array.shape).encode("utf-8"))
    digest.update(array.tobytes())


def _atmo_field_cache_token(field: xr.DataArray) -> str:
    """Return a stable identity for one atmospheric input field.

    M4 can receive distinct providers or explicit atmospheric overrides for
    the same Sentinel-2 acquisition.  Scene bounds and sensing time alone
    therefore do not identify a safe reprojection-cache entry.  Prefer a
    provider-supplied provenance key; otherwise fingerprint the small source
    field together with its grid definition.
    """
    explicit = field.attrs.get(_ATMO_CACHE_ID_ATTR)
    if explicit is not None:
        value = str(explicit).strip()
        if value:
            return f"provenance:{value}"

    digest = hashlib.sha256()
    digest.update(b"siac-atmo-field-v1\0")
    for dim in field.dims:
        digest.update(dim.encode("utf-8"))
        digest.update(b"\0")
    _update_atmo_cache_hash(digest, field.values)
    for coord_name in sorted(field.coords):
        coord = field.coords[coord_name]
        digest.update(coord_name.encode("utf-8"))
        digest.update(b"\0")
        for dim in coord.dims:
            digest.update(dim.encode("utf-8"))
            digest.update(b"\0")
        _update_atmo_cache_hash(digest, coord.values)
        # CRS/transform metadata can live on the scalar spatial_ref coordinate.
        for key, value in sorted(coord.attrs.items()):
            digest.update(str(key).encode("utf-8"))
            digest.update(repr(value).encode("utf-8"))
            digest.update(b"\0")
    for key in ("crs", "transform", "grid_mapping"):
        if key in field.attrs:
            digest.update(key.encode("utf-8"))
            digest.update(repr(field.attrs[key]).encode("utf-8"))
            digest.update(b"\0")
    return f"content:{digest.hexdigest()}"


def _atmo_cache_source_identity(
    scene_identity: str,
    field_name: str,
    field: xr.DataArray,
) -> str:
    """Build the content-addressable identity for one M2-to-M4 field."""
    return f"atmo-v2:{scene_identity}:{field_name}:{_atmo_field_cache_token(field)}"


def _resample_atmo_state(
    atmo: AtmosphericState,
    target_shape: tuple[int, int],
    *,
    template: xr.DataArray | None = None,
    cache_dir: Path | str | None = None,
    scene_identity: str | None = None,
) -> AtmosphericState:
    """Resample all atmospheric state fields to target shape.

    When both ``cache_dir`` and ``scene_identity`` are provided, each field
    routes through ``cached_reproject_match``.  The key includes the scene,
    field name, and the field's source provenance (or a source-grid/content
    fingerprint when provenance is unavailable).  This keeps cache hits for
    repeated runs while preventing a scalar override or different provider
    from being reused as a spatial atmospheric field for the same scene.
    """
    if cache_dir is not None and scene_identity:
        from siac.geo.cached_reprojection import cached_reproject_match
        from siac.geo.target_grid import TargetGrid

        def _cached(field: xr.DataArray, field_name: str, resampling: str) -> xr.DataArray:
            """Cache one field's reprojection under a per-field identity."""
            if template is None:
                return _resample_da(field, target_shape, resampling, template=template)
            # Only cache the same-CRS reproject path — if the field is
            # already on the target grid we still pay one call to
            # _resample_da (which is a near-no-op), but we don't add
            # a cache layer where there's nothing to cache.
            try:
                return cached_reproject_match(
                    field,
                    TargetGrid.from_template(template),
                    source_identity=_atmo_cache_source_identity(
                        scene_identity,
                        field_name,
                        field,
                    ),
                    cache_dir=cache_dir,
                    resampling=resampling,
                )
            except Exception:
                # Any failure inside the cache wrapper (e.g. missing CRS
                # on the source) falls back to the legacy resample path.
                logger.debug(
                    "Atmo prior cache wrapper failed for %s; using legacy resample",
                    field_name,
                    exc_info=True,
                )
                return _resample_da(field, target_shape, resampling, template=template)

        return AtmosphericState(
            aot=_cached(atmo.aot, "aot", "bilinear"),
            tcwv=_cached(atmo.tcwv, "tcwv", "bilinear"),
            tco3=_cached(atmo.tco3, "tco3", "bilinear"),
            aot_unc=_cached(atmo.aot_unc, "aot_unc", "bilinear"),
            tcwv_unc=_cached(atmo.tcwv_unc, "tcwv_unc", "bilinear"),
            tco3_unc=_cached(atmo.tco3_unc, "tco3_unc", "bilinear"),
            elevation=_cached(atmo.elevation, "elevation", "bilinear"),
        )

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

    tau_predictor = sp.tau_predictor
    if tau_predictor is not None and "localizer" in tau_predictor:
        # Resample the climatological localizer planes onto the solver grid so
        # M5 can evaluate the tau-dependent prediction per candidate AOD.
        localizer_grid = _resample_da(
            tau_predictor["localizer"], target_shape, "area", template=template
        )
        tau_predictor = {
            **tau_predictor,
            "localizer_grid": np.asarray(localizer_grid.values, dtype=np.float64),
        }

    return SurfacePrior(
        boa=_resample_da(sp.boa, target_shape, "area", template=template),
        boa_unc=_resample_da(sp.boa_unc, target_shape, "area", template=template),
        kernels=sp.kernels,  # keep original (used for spectral, not spatial)
        mask=_resample_cloud_mask(
            mask,
            target_shape,
            template=template,
            assume_aligned_native_grid=True,
        ),
        monthly_composites=sp.monthly_composites,
        tau_predictor=tau_predictor,
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
        tau_predictor=sp.tau_predictor,
    )


def _resolve_solver_bands(
    sensor_config: Any,
    requested_band_names: tuple[str, ...] | None,
) -> list[Any]:
    """Resolve the full solver-band set needed for grid assembly.

    Default aerosol bands remain the baseline. When staged solving requests
    extra bands, include them here as well so M4 assembles everything M5 may
    later select.
    """
    default_bands = list(sensor_config.default_aerosol_solver_bands())
    if not requested_band_names:
        return default_bands

    requested = tuple(
        dict.fromkeys(str(name).strip() for name in requested_band_names if str(name).strip())
    )
    available_by_name = {band.name: band for band in sensor_config.bands}
    missing = [name for name in requested if name not in available_by_name]
    if missing:
        raise ValueError(
            "Requested solver band(s) are not available for sensor "
            f"{getattr(sensor_config, 'sensor_id', 'unknown')}: {', '.join(missing)}"
        )

    selected_names = {band.name for band in default_bands}
    selected_names.update(requested)
    return [band for band in sensor_config.bands if band.name in selected_names]


# ── Public API ─────────────────────────────────────────────────────────


def _resolve_psf_shift_reference(
    obs: ObservationBundle,
    surface: SurfacePrior,
    reference_bands: tuple[str, ...],
    target_shape: tuple[int, int],
    target_template: xr.DataArray,
    toa_band_loader: Any,
) -> tuple[xr.DataArray | None, xr.DataArray | None]:
    """Grid a shared AOT-insensitive (TOA, prior) reference band for the shift fit.

    Returns the first ``reference_bands`` entry present in BOTH the (pre-alignment)
    prior and the loadable TOA, each resampled to the solver grid — or ``(None,
    None)`` when no shared band exists (shift then skipped, TOA only convolved).
    """
    boa = surface.boa
    if "band" not in boa.dims or "band" not in boa.coords:
        return None, None
    prior_band_names = {str(value) for value in np.asarray(boa.coords["band"].values).tolist()}
    for name in reference_bands:
        if name not in prior_band_names:
            continue
        toa_native: xr.DataArray | None = None
        if name in obs.toa.data_vars:
            toa_native = obs.toa[name]
        elif callable(toa_band_loader):
            try:
                toa_native = toa_band_loader(name, native=True)
            except (KeyError, RuntimeError):
                toa_native = None
        if toa_native is None:
            continue
        toa_ref = _resample_da(toa_native, target_shape, "area", template=target_template)
        prior_ref = _resample_da(boa.sel(band=name), target_shape, "area", template=target_template)
        return toa_ref, prior_ref
    return None, None


def assemble_grids(
    obs: ObservationBundle,
    atmo: AtmosphericState,
    surface: SurfacePrior,
    rt_model: Any,
    aux_resolution_m: float = 500.0,
    aerosol_resolution_m: float = 120.0,
    sharp_transition_filter: Any | None = None,
    water_mask_path: str | Path | None = None,
    water_mask_cache_dir: str | Path | None = None,
    water_mask_buffer_pixels: int = 0,
    solver_band_names: tuple[str, ...] | None = None,
    reproject_cache_dir: str | Path | None = None,
    dem_path: str | Path | None = None,
    toa_psf_config: ToaPsfConfig | None = None,
    resample_workers: int = 1,
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
        aux_resolution_m: Auxiliary product resolution metadata.
        aerosol_resolution_m: Target resolution for AOT/TCWV retrieval (default 120 m).

    Returns:
        SolverInputBundle ready for the solver.
    """
    # 0. Validate inputs at the M4 boundary.
    t0 = __import__("time").monotonic()
    if not list(obs.toa.data_vars):
        raise ValueError("M4 grid assembly received an empty TOA dataset (no data variables).")
    if aerosol_resolution_m <= 0:
        raise ValueError(f"aerosol_resolution_m must be > 0, got {aerosol_resolution_m}")
    if obs.bounds[0] >= obs.bounds[2] or obs.bounds[1] >= obs.bounds[3]:
        raise ValueError(f"Invalid observation bounds: {obs.bounds}")

    # 1. Select solver bands using the sensor defaults.
    bands = _resolve_solver_bands(obs.sensor_config, solver_band_names)
    logger.info(f"Selected {len(bands)} solver bands: {[b.name for b in bands]}")

    # 2. Determine the solver grid from the configured aerosol retrieval
    #    resolution and scene bounds.
    resolved_aerosol_resolution = float(aerosol_resolution_m)
    first_var = list(obs.toa.data_vars)[0]
    native_shape = obs.toa[first_var].shape  # (y, x)
    target_template = _build_target_template(obs.bounds, obs.crs, resolved_aerosol_resolution)
    # Surface-driven "resolve on the prior's native grid": when the surface prior
    # carries its own (composite native) grid, adopt it as the solver target so
    # the prior is never resampled (the shares_template_grid no-op in
    # _resample_surface_prior) and TOA, geometry, atmo and masks are all
    # co-registered onto it. Default keeps the observation-bounds grid above
    # (byte-identical legacy behaviour when surface.solver_grid is None).
    if surface.solver_grid is not None:
        target_template = _ensure_template_transform(surface.solver_grid)
        # The adopted grid sets the true solver pixel size; align
        # resolved_aerosol_resolution to it so physical->pixel conversions (the
        # solver's spatial pooling window) stay correct even when the configured
        # aerosol resolution differs from the composite's native resolution.
        adopted_res = abs(float(target_template.rio.resolution()[0]))
        if adopted_res > 0:
            resolved_aerosol_resolution = adopted_res
        logger.info(
            "M4 adopting surface-prior native solver grid: shape=%s @ %.1fm",
            target_template.shape,
            resolved_aerosol_resolution,
        )
    target_shape = target_template.shape
    logger.info(
        f"Resampling from {native_shape} to {target_shape} @ {resolved_aerosol_resolution}m"
    )

    # 3. Detect native sharp transitions before solver-grid averaging.
    band_names = [b.name for b in bands]
    sharp_transition_mask: xr.DataArray | None = None
    water_mask: xr.DataArray | None = None
    if sharp_transition_filter is not None and bool(
        getattr(sharp_transition_filter, "enabled", False)
    ):
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
            assume_aligned_native_grid=(surface.solver_grid is None),
        )

    if water_mask_path is not None:
        native_water_mask = load_water_mask_subset(
            obs.bounds,
            obs.crs,
            source=water_mask_path,
            cache_dir=water_mask_cache_dir,
        )
        if water_mask_buffer_pixels > 0:
            native_water_reference = native_water_mask
            buffered_water = _dilate_mask(
                np.asarray(native_water_mask.values, dtype=bool),
                int(water_mask_buffer_pixels),
            )
            native_water_mask = xr.DataArray(
                buffered_water,
                dims=native_water_mask.dims,
                coords=native_water_mask.coords,
                attrs=native_water_mask.attrs,
                name=native_water_mask.name,
            )
            native_water_mask = copy_spatial_metadata_like(
                native_water_mask, native_water_reference
            )
        water_mask = _resample_external_exclusion_mask(
            native_water_mask,
            target_shape,
            template=target_template,
        )

    # 4. Resample TOA for solver bands into (band, y, x) DataArray.
    #    Bands already in obs.toa (preloaded at native resolution) are
    #    resampled from native to the solver grid.  Missing bands are loaded
    #    lazily via the TOA band loader at native resolution and then
    #    resampled – this avoids an expensive intermediate reproject (e.g.
    #    20 m → 10 m) for bands that only need to reach the aerosol grid.
    from concurrent.futures import ThreadPoolExecutor

    _toa_band_loader = obs.toa.attrs.get("_siac_toa_band_loader")

    # On the prior native grid, prefer the RAW native band (reprojected once onto
    # the composite grid) over the M1-harmonized obs.toa[bn] — which was already
    # resampled to the 10 m reference grid — avoiding a 60m->10m->grid double
    # resample on the AOD-load-bearing bands (notably the 60 m deep-blue B01).
    # Matches the harness "raw TOA reprojected once onto the composite grid".
    use_raw_toa = surface.solver_grid is not None and callable(_toa_band_loader)

    def _resample_toa_band(bn: str) -> tuple[str, xr.DataArray] | None:
        if use_raw_toa:
            try:
                band_da = _toa_band_loader(bn, native=True)
                return bn, _resample_da(band_da, target_shape, "area", template=target_template)
            except (KeyError, RuntimeError):
                logger.warning("Raw TOA load failed for %s; using the preloaded band", bn)
        if bn in obs.toa.data_vars:
            return bn, _resample_da(obs.toa[bn], target_shape, "area", template=target_template)
        if callable(_toa_band_loader):
            try:
                band_da = _toa_band_loader(bn, native=True)
                return bn, _resample_da(band_da, target_shape, "area", template=target_template)
            except (KeyError, RuntimeError):
                logger.warning("Could not load solver band %s via band loader", bn)
        return None

    # 5. Resample geometry, cloud mask, atmo, surface — all independent,
    #    run alongside the TOA band resampling.
    surface_aligned = _align_surface_prior_to_bands(surface, band_names)

    # resample_workers=1 (default) keeps M4 bit-reproducible: concurrent GDAL
    # warps in this pool return non-deterministic grids, which moves the solver
    # argmin at flat-cost sites (±2 within-EE sites run-to-run).
    with ThreadPoolExecutor(max_workers=max(1, int(resample_workers))) as pool:
        toa_futures = {bn: pool.submit(_resample_toa_band, bn) for bn in band_names}
        geom_future = pool.submit(
            _resample_geometry, obs.geometry, target_shape, template=target_template
        )
        cloud_future = pool.submit(
            _resample_cloud_mask,
            obs.cloud_mask,
            target_shape,
            template=target_template,
            # On the prior's native grid the cloud mask (observation grid) and the
            # solver target differ in footprint, so it must be reprojected
            # geospatially rather than shape-only index-remapped (which assumes a
            # shared footprint). Default (obs-bounds target) stays shape-only.
            assume_aligned_native_grid=(surface.solver_grid is None),
        )
        # Wave 19b: derive a per-scene identity for the atmo-prior cache.
        # The atmo prior depends on (bounds, sensing time, source) — those
        # uniquely pin the bytes returned by CAMSProvider.get_prior. We
        # round to whole metres / second to absorb any insignificant
        # float-roundtrip variance across runs.
        atmo_scene_identity: str | None = None
        if reproject_cache_dir is not None:
            obs_time = obs.metadata.get("observation_time")
            obs_time_key = obs_time.isoformat() if obs_time is not None else "unknown"
            atmo_scene_identity = (
                f"{obs.crs}:"
                f"{round(obs.bounds[0])}:{round(obs.bounds[1])}:"
                f"{round(obs.bounds[2])}:{round(obs.bounds[3])}:"
                f"{obs_time_key}:res{int(round(resolved_aerosol_resolution))}m"
            )
        atmo_future = pool.submit(
            _resample_atmo_state,
            atmo,
            target_shape,
            template=target_template,
            cache_dir=reproject_cache_dir,
            scene_identity=atmo_scene_identity,
        )
        surface_future = pool.submit(
            _resample_surface_prior, surface_aligned, target_shape, template=target_template
        )

        toa_arrays = []
        resolved_band_names: list[str] = []
        for bn in band_names:
            result = toa_futures[bn].result()
            if result is not None:
                resolved_band_names.append(result[0])
                toa_arrays.append(result[1])

        geometry = geom_future.result()
        cloud_mask = cloud_future.result()
        atmo_resampled = atmo_future.result()
        surface_resampled = surface_future.result()

    # Populate terrain elevation per solver-grid pixel from the DEM (the CAMS
    # prior carries none, so the RT elevation correction would otherwise be a
    # no-op — under-attributing AOD over high terrain). Reading at the solver
    # grid (vs the coarse CAMS cell) avoids the regional over-estimate in
    # valleys. A falsy / sea-level ``dem_path`` (see siac.geo.dem) leaves the
    # resampled (sea-level) elevation in place, reproducing the legacy no-op.
    from siac.geo.dem import read_elevation_km, use_sea_level_elevation

    dem_arg = None if dem_path is None else str(dem_path)
    if not use_sea_level_elevation(dem_arg):
        from dataclasses import replace

        atmo_resampled = replace(
            atmo_resampled,
            elevation=read_elevation_km(target_template, dem_arg),
        )

    # tau-dependent prior: the solver additionally needs the scene anchor-band
    # TOA on the solver grid to re-correct at each candidate AOD.
    if surface_resampled.tau_predictor is not None:
        payload = dict(surface_resampled.tau_predictor)
        anchor_planes = []
        anchor_ok = True
        for anchor_name in payload.get("anchor_bands", ()):
            resampled = _resample_toa_band(anchor_name)
            if resampled is None:
                logger.warning(
                    "tau-dependent prior: anchor band %s unavailable; disabling.", anchor_name
                )
                anchor_ok = False
                break
            anchor_planes.append(np.asarray(resampled[1].values, dtype=np.float64))
        if anchor_ok and anchor_planes:
            from dataclasses import replace as _replace

            payload["anchor_toa_grid"] = np.stack(anchor_planes, axis=0)
            payload["anchor_sensor_bands"] = [
                obs.sensor_config.get_band(name) for name in payload.get("anchor_bands", ())
            ]
            surface_resampled = _replace(surface_resampled, tau_predictor=payload)
        else:
            from dataclasses import replace as _replace

            surface_resampled = _replace(surface_resampled, tau_predictor=None)

    if toa_arrays:
        toa_da = xr.concat(toa_arrays, dim="band")
        toa_da = toa_da.assign_coords(band=resolved_band_names)
    else:
        # Fallback: use first available variable
        first = list(obs.toa.data_vars)[0]
        toa_da = _resample_da(
            obs.toa[first], target_shape, "area", template=target_template
        ).expand_dims("band")
        toa_da = toa_da.assign_coords(band=[first])
    toa_da = copy_spatial_metadata_like(toa_da, target_template)

    # Observation-side PSF: convolve the gridded TOA with the fixed-width sensor
    # PSF and co-register it to the coarse MODIS prior with a per-scene integer
    # shift (SIAC v1 methodology). The shift is fitted on an AOT-insensitive SWIR
    # reference band that the prior carries only for this purpose; that band is
    # dropped from the bundle by `_align_surface_prior_to_bands`, so the solver is
    # unchanged. No-op when `toa_psf_config` is None / disabled.
    if toa_psf_config is not None and getattr(toa_psf_config, "enabled", False):
        from siac.algorithms.grid.toa_psf import psf_convolve_and_align_toa

        fine_mode = getattr(toa_psf_config, "convolve_resolution", "grid") == "native10m"
        if fine_mode:
            # Convolve + fit the co-registration shift at the PSF calibration
            # resolution (~10 m), where the MODIS-footprint structure (~500 m ≈
            # 50 px) is resolvable and the integer shift is meaningful — then
            # downsample to the solver grid. On the coarse solver grid the PSF
            # blur washes out the ~4 px features, leaving the shift unlocalizable.
            fine_res = float(getattr(toa_psf_config, "target_resolution_m", 10.0))
            fine_template = _build_target_template(obs.bounds, obs.crs, fine_res)
            fine_shape = fine_template.shape
            # `resolved_band_names` already excludes bands that failed to load in
            # the solver-grid pass above, so each one resolves again here.
            fine_bands: list[xr.DataArray] = []
            for bn in resolved_band_names:
                src = obs.toa[bn] if bn in obs.toa.data_vars else _toa_band_loader(bn, native=True)
                fine_bands.append(_resample_da(src, fine_shape, "area", template=fine_template))
            fine_toa = xr.concat(fine_bands, dim="band").assign_coords(band=resolved_band_names)
            fine_toa = copy_spatial_metadata_like(fine_toa, fine_template)
            toa_ref, prior_ref = _resolve_psf_shift_reference(
                obs,
                surface,
                toa_psf_config.reference_bands,
                fine_shape,
                fine_template,
                _toa_band_loader,
            )
            conv_fine, shift_fit = psf_convolve_and_align_toa(
                fine_toa, toa_ref, prior_ref, None, grid_resolution_m=fine_res, cfg=toa_psf_config
            )
            toa_da = copy_spatial_metadata_like(
                _resample_da(conv_fine, target_shape, "area", template=target_template),
                target_template,
            )
        else:
            toa_ref_grid, prior_ref_grid = _resolve_psf_shift_reference(
                obs,
                surface,
                toa_psf_config.reference_bands,
                target_shape,
                target_template,
                _toa_band_loader,
            )
            # `_resample_surface_prior` already reduces a banded mask to 2-D (y, x).
            prior_valid = (
                None
                if surface_resampled.mask is None
                else np.asarray(surface_resampled.mask.values, dtype=bool)
            )
            toa_da, shift_fit = psf_convolve_and_align_toa(
                toa_da,
                toa_ref_grid,
                prior_ref_grid,
                prior_valid,
                grid_resolution_m=resolved_aerosol_resolution,
                cfg=toa_psf_config,
            )
        logger.info(
            "M4 PSF-on-TOA (%s): shift=(dx=%d, dy=%d) r=%.3f accepted=%s",
            "native10m" if fine_mode else "grid",
            shift_fit.dx,
            shift_fit.dy,
            shift_fit.correlation,
            shift_fit.accepted,
        )

    bundle = SolverInputBundle(
        toa=toa_da,
        geometry=geometry,
        cloud_mask=cloud_mask,
        sharp_transition_mask=sharp_transition_mask,
        water_mask=water_mask,
        sensor_config=obs.sensor_config,
        bands=bands,
        atmo_prior=atmo_resampled,
        surface_prior=surface_resampled,
        rt_model=rt_model,
        aux_resolution_m=aux_resolution_m,
        aerosol_resolution_m=resolved_aerosol_resolution,
    )

    # Post-assembly validation: ensure shapes are consistent.
    expected_shape = target_shape
    if toa_da.shape[-2:] != expected_shape:
        raise ValueError(
            f"M4 post-assembly shape mismatch: TOA {toa_da.shape[-2:]} != target {expected_shape}"
        )

    elapsed = __import__("time").monotonic() - t0
    logger.info(
        "M4 grid assembly complete: %d bands, shape=%s, resolution=%.1fm (%.2fs)",
        len(resolved_band_names),
        target_shape,
        resolved_aerosol_resolution,
        elapsed,
    )
    return bundle
