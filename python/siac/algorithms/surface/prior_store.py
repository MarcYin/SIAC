"""
Pre-built surface prior store.

Loads surface reflectance priors from a Zarr store, interpolates to
the observation day-of-year, crops to the AOI, and projects reference
bands to sensor bands using the spectral convolution module.

See PLANS.md §9.6 for architecture details.

Zarr store layout::

    /{tile_id}/
        reflectance   (doy, band, y, x)
        uncertainty   (doy, band, y, x)
        n_obs         (doy, y, x)
        quality       (doy, y, x)
        wavelengths   (band,)
        .zattrs       → crs, bounds, resolution_m, doy_values

"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from siac.runtime import BRDFKernelWeights, GeometryAngles, SurfacePrior

if TYPE_CHECKING:
    from datetime import datetime

    from siac.domain import SensorConfig

logger = logging.getLogger(__name__)


def _doy_from_datetime(dt: datetime) -> int:
    """Extract day-of-year (1–366) from a datetime."""
    return dt.timetuple().tm_yday


def _interpolate_doy(
    data: xr.DataArray,
    doy_values: np.ndarray,
    target_doy: int,
) -> xr.DataArray:
    """Linearly interpolate along the 'doy' dimension to target DOY.

    Handles wrap-around at year boundaries (DOY 365 → DOY 1).

    Args:
        data: Array with a 'doy' dimension.
        doy_values: The DOY coordinate values.
        target_doy: The target day-of-year.

    Returns:
        Interpolated DataArray with the 'doy' dimension removed.
    """
    doy_sorted = np.sort(doy_values)

    if target_doy in doy_sorted:
        return data.sel(doy=target_doy).drop_vars("doy", errors="ignore")

    # Find bracketing DOYs
    idx_after = np.searchsorted(doy_sorted, target_doy)

    if idx_after == 0:
        # target is before first DOY — wrap from last DOY of previous year
        doy_lo = doy_sorted[-1]
        doy_hi = doy_sorted[0]
        gap = (365 - doy_lo) + doy_hi
        dist_lo = (365 - doy_lo) + target_doy if target_doy < doy_lo else target_doy - doy_lo
        weight = dist_lo / gap if gap > 0 else 0.5
    elif idx_after >= len(doy_sorted):
        # target is after last DOY — wrap to first DOY of next year
        doy_lo = doy_sorted[-1]
        doy_hi = doy_sorted[0]
        gap = (365 - doy_lo) + doy_hi
        dist_lo = target_doy - doy_lo
        weight = dist_lo / gap if gap > 0 else 0.5
    else:
        doy_lo = doy_sorted[idx_after - 1]
        doy_hi = doy_sorted[idx_after]
        gap = doy_hi - doy_lo
        weight = (target_doy - doy_lo) / gap if gap > 0 else 0.5

    val_lo = data.sel(doy=doy_lo).drop_vars("doy", errors="ignore")
    val_hi = data.sel(doy=doy_hi).drop_vars("doy", errors="ignore")

    return val_lo * (1.0 - weight) + val_hi * weight


def _select_tiles(store_path: Path, bounds: tuple[float, float, float, float]) -> list[str]:
    """Select tile IDs whose spatial extent overlaps the given bounds.

    Looks for subdirectories in the store and checks their .zattrs for
    spatial bounds.  If no attrs are found, the tile is included by default.
    """
    import json

    xmin, ymin, xmax, ymax = bounds
    tiles: list[str] = []

    if not store_path.is_dir():
        return tiles

    for entry in sorted(store_path.iterdir()):
        if not entry.is_dir():
            continue
        tile_id = entry.name

        zattrs_path = entry / ".zattrs"
        if zattrs_path.exists():
            try:
                attrs = json.loads(zattrs_path.read_text())
                tb = attrs.get("bounds")
                if tb is not None:
                    txmin, tymin, txmax, tymax = tb
                    # Check overlap
                    if txmax < xmin or txmin > xmax or tymax < ymin or tymin > ymax:
                        continue
            except (json.JSONDecodeError, KeyError, TypeError):
                pass  # include tile if attrs are broken

        tiles.append(tile_id)

    return tiles


def _crop_to_bounds(
    data: xr.DataArray,
    tile_bounds: tuple[float, float, float, float],
    target_bounds: tuple[float, float, float, float],
    resolution_m: float,
) -> xr.DataArray:
    """Crop a raster DataArray to target bounds using pixel indices.

    Assumes data is on a regular grid with origin at tile_bounds[:2].
    """
    txmin, tymin, txmax, tymax = tile_bounds
    xmin, ymin, xmax, ymax = target_bounds

    ny, nx = data.shape[-2], data.shape[-1]

    # Compute pixel coordinates
    col_start = max(0, int((xmin - txmin) / resolution_m))
    col_end = min(nx, int(np.ceil((xmax - txmin) / resolution_m)))
    row_start = max(0, int((ymin - tymin) / resolution_m))
    row_end = min(ny, int(np.ceil((ymax - tymin) / resolution_m)))

    if data.ndim == 2:
        return data[row_start:row_end, col_start:col_end]
    elif data.ndim == 3:
        return data[:, row_start:row_end, col_start:col_end]
    return data


def _project_to_sensor(
    reflectance: np.ndarray,
    store_wavelengths: np.ndarray,
    sensor_config: SensorConfig,
) -> np.ndarray:
    """Project reference-band reflectance to sensor bands.

    Falls back to nearest-wavelength selection.

    Args:
        reflectance: (n_ref_bands, y, x) array in the reference basis.
        store_wavelengths: Center wavelengths (nm) of the store bands.
        sensor_config: Target sensor configuration.

    Returns:
        (n_sensor_bands, y, x) array with sensor-band reflectance.
    """
    n_sensor = len(sensor_config.bands)
    spatial_shape = reflectance.shape[1:]
    result = np.zeros((n_sensor, *spatial_shape), dtype=reflectance.dtype)

    for i, band in enumerate(sensor_config.bands):
        # Find nearest reference band by wavelength
        diffs = np.abs(store_wavelengths - band.center_wavelength)
        nearest_idx = int(np.argmin(diffs))
        if nearest_idx < reflectance.shape[0]:
            result[i] = reflectance[nearest_idx]

    return result


class PrebuiltPriorStore:
    """Load pre-built surface prior from a Zarr store.

    Implements the M3 surface prior provider interface::

        (bounds, crs, obs_time, sensor_config, geometry, resolution) -> SurfacePrior

    The ``geometry`` argument is accepted for interface compatibility but is
    **ignored** — the prior is climatological and does not depend on viewing
    geometry.

    Args:
        store_path: Path to the root Zarr store directory.
    """

    def __init__(self, store_path: Path):
        self.store_path = Path(store_path)

    def get_surface_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        sensor_config: SensorConfig,
        geometry: GeometryAngles,
        resolution: float,
    ) -> SurfacePrior:
        """Load, interpolate, crop, and project the prior.

        Steps:
            1. Select tile(s) covering the AOI.
            2. Load reflectance/uncertainty data.
            3. Interpolate to observation DOY.
            4. Crop to AOI bounds.
            5. Project reference bands → sensor bands.
            6. Return ``SurfacePrior``.
        """
        del crs, geometry
        target_doy = _doy_from_datetime(obs_time)
        tiles = _select_tiles(self.store_path, bounds)

        if not tiles:
            raise ValueError(f"No tiles in {self.store_path} overlap bounds {bounds}")

        # For simplicity, use the first overlapping tile.
        # Multi-tile mosaicking can be added later.
        tile_path = self.store_path / tiles[0]
        ds = xr.open_zarr(str(tile_path))

        # Extract DOY coordinate
        doy_values = (
            ds["doy"].values
            if "doy" in ds.coords
            else ds.coords.get("doy", np.array([target_doy])).values
        )

        # Interpolate to target DOY
        refl = _interpolate_doy(ds["reflectance"], doy_values, target_doy)
        unc = _interpolate_doy(ds["uncertainty"], doy_values, target_doy)

        # Crop to AOI
        tile_attrs = ds.attrs
        tile_bounds = tuple(tile_attrs.get("bounds", bounds))
        tile_res = float(tile_attrs.get("resolution_m", resolution))

        refl = _crop_to_bounds(refl, tile_bounds, bounds, tile_res)
        unc = _crop_to_bounds(unc, tile_bounds, bounds, tile_res)

        # Get store wavelengths
        if "wavelengths" in ds:
            store_wl = ds["wavelengths"].values
        else:
            # Fallback: MODIS land band centers
            store_wl = np.array([645.0, 858.5, 469.0, 555.0, 1240.0, 1640.0, 2130.0])

        # Project to sensor bands
        refl_vals = refl.values if hasattr(refl, "values") else np.asarray(refl)
        unc_vals = unc.values if hasattr(unc, "values") else np.asarray(unc)

        sensor_refl = _project_to_sensor(refl_vals, store_wl, sensor_config)
        sensor_unc = _project_to_sensor(unc_vals, store_wl, sensor_config)

        # Collapse to 2D by averaging across sensor bands for the scalar prior
        # (the solver expects a single-band boa/boa_unc)
        boa = xr.DataArray(
            np.nanmean(sensor_refl, axis=0).astype(np.float32),
            dims=["y", "x"],
        )
        boa_unc = xr.DataArray(
            np.nanmean(sensor_unc, axis=0).astype(np.float32),
            dims=["y", "x"],
        )

        # Validity mask
        mask = xr.DataArray(
            np.isfinite(boa.values) & (boa.values > 0) & (boa.values < 1.5),
            dims=["y", "x"],
        )

        # Dummy BRDF kernel weights (not meaningful for pre-built priors)
        shape = boa.shape
        kernels = BRDFKernelWeights(
            f0=xr.DataArray(np.full(shape, 0.0, dtype=np.float32), dims=["y", "x"]),
            f1=xr.DataArray(np.full(shape, 0.0, dtype=np.float32), dims=["y", "x"]),
            f2=xr.DataArray(np.full(shape, 0.0, dtype=np.float32), dims=["y", "x"]),
            f0_unc=xr.DataArray(np.full(shape, 0.0, dtype=np.float32), dims=["y", "x"]),
            f1_unc=xr.DataArray(np.full(shape, 0.0, dtype=np.float32), dims=["y", "x"]),
            f2_unc=xr.DataArray(np.full(shape, 0.0, dtype=np.float32), dims=["y", "x"]),
        )

        return SurfacePrior(
            boa=boa,
            boa_unc=boa_unc,
            kernels=kernels,
            mask=mask,
        )
