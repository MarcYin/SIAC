"""Read terrain elevation from a configured DEM onto an arbitrary target grid.

The atmospheric RT applies an elevation correction (shorter air column at
altitude → less Rayleigh/molecular path radiance), but the CAMS prior carries
no elevation. Populating it from a DEM at the *solver* grid makes that correction
physical; leaving it at sea level (0 km) reproduces the legacy no-op behaviour.

**Forcing sea level.** Pass ``dem_path`` as ``None`` or an empty/blank string
(``""``), or set ``paths.dem`` to one of the sentinels ``"none"`` / ``"sealevel"``
/ ``"0"`` (case-insensitive), to get a 0 km elevation field directly — no DEM is
read and the RT elevation correction becomes a no-op.
"""

from __future__ import annotations

import logging

import numpy as np
import xarray as xr

logger = logging.getLogger(__name__)

# Explicit "use sea level (elevation = 0)" sentinels for ``paths.dem`` — handy
# because TOML has no null literal, so a user can write dem = "none".
_SEA_LEVEL_SENTINELS = frozenset(
    {"none", "sealevel", "sea_level", "sea-level", "0", "off", "false"}
)


def use_sea_level_elevation(dem_path: str | None) -> bool:
    """Whether ``dem_path`` requests a sea-level (0 km) elevation field."""
    if not dem_path:
        return True
    stripped = str(dem_path).strip().lower()
    return stripped == "" or stripped in _SEA_LEVEL_SENTINELS


def read_elevation_km(template: xr.DataArray, dem_path: str | None) -> xr.DataArray:
    """Surface elevation (km above MSL) sampled onto ``template``'s grid.

    Reads ``dem_path`` (e.g. the Copernicus GLO-30 ``/vsicurl`` VRT) windowed to
    the template footprint and resamples it to the template grid. Returns sea
    level (0 km) when ``dem_path`` is falsy / a sea-level sentinel (see
    :func:`use_sea_level_elevation`) or when the read fails — so a failed DEM
    read degrades to the legacy behaviour rather than aborting the run.
    """
    if use_sea_level_elevation(dem_path):
        return xr.zeros_like(template, dtype=np.float32)
    try:
        import rioxarray  # noqa: F401  (registers the .rio accessor)

        from siac.geo.reprojection import reproject_match, transform_bounds

        target_crs = template.rio.crs
        if target_crs is None:
            return xr.zeros_like(template, dtype=np.float32)
        minx, miny, maxx, maxy = (float(v) for v in template.rio.bounds())

        dem = rioxarray.open_rasterio(str(dem_path), masked=True)
        if "band" in dem.dims:
            dem = dem.isel(band=0, drop=True)
        dem_crs = dem.rio.crs
        if dem_crs is not None and str(dem_crs) != str(target_crs):
            clip = transform_bounds((minx, miny, maxx, maxy), str(target_crs), str(dem_crs))
        else:
            clip = (minx, miny, maxx, maxy)
        # Pad the window so reprojection has neighbours and a tiny AOI still spans
        # at least one DEM posting.
        pad_x = max(0.02, abs(clip[2] - clip[0]) * 0.1)
        pad_y = max(0.02, abs(clip[3] - clip[1]) * 0.1)
        dem = dem.rio.clip_box(clip[0] - pad_x, clip[1] - pad_y, clip[2] + pad_x, clip[3] + pad_y)

        elev_m = reproject_match(dem.astype("float32"), template, resampling="bilinear")
        elev_km = np.asarray(elev_m.values, dtype=np.float32) / np.float32(1000.0)
        elev_km = np.where(np.isfinite(elev_km), elev_km, np.float32(0.0)).astype(np.float32)
        return xr.DataArray(elev_km, dims=template.dims, coords=template.coords)
    except Exception:
        logger.warning(
            "DEM elevation read failed (%s); falling back to sea level (0 km).",
            dem_path,
            exc_info=True,
        )
        return xr.zeros_like(template, dtype=np.float32)
