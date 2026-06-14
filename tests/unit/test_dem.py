"""Tests for siac.geo.dem (DEM elevation read + sea-level override)."""

from __future__ import annotations

import numpy as np
import pytest
import rasterio
import rioxarray  # noqa: F401  (registers the .rio accessor)
import xarray as xr
from rasterio.transform import from_origin

from siac.geo.dem import read_elevation_km, use_sea_level_elevation


def _write_dem(path: str, value_m: float = 2400.0) -> None:
    """A constant-elevation DEM over lon/lat [0, 0.12] x [0.88, 1.0] (EPSG:4326)."""
    n = 120
    res = 0.001
    transform = from_origin(0.0, 1.0, res, res)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=n,
        width=n,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
        nodata=-9999.0,
    ) as dst:
        dst.write(np.full((n, n), value_m, dtype=np.float32), 1)


def _template() -> xr.DataArray:
    """A small 4326 target grid wholly inside the synthetic DEM."""
    res = 0.002
    n = 20
    xs = 0.02 + res * (np.arange(n) + 0.5)
    ys = 0.98 - res * (np.arange(n) + 0.5)
    da = xr.DataArray(
        np.zeros((n, n), dtype=np.float32), dims=["y", "x"], coords={"y": ys, "x": xs}
    )
    return da.rio.write_crs("EPSG:4326")


def test_read_elevation_km_returns_km_on_template_grid(tmp_path) -> None:
    dem = tmp_path / "dem.tif"
    _write_dem(str(dem), value_m=2400.0)

    elev = read_elevation_km(_template(), str(dem))

    assert elev.shape == (20, 20)
    assert elev.dims == ("y", "x")
    # 2400 m -> 2.4 km, sampled onto the template grid.
    np.testing.assert_allclose(float(np.nanmean(elev.values)), 2.4, atol=0.05)


@pytest.mark.parametrize("dem_path", [None, "", "  ", "none", "NONE", "sea_level", "0", "off"])
def test_read_elevation_km_sea_level_override_returns_zeros(dem_path) -> None:
    # A falsy or sentinel dem_path forces sea level (0 km) directly.
    assert use_sea_level_elevation(dem_path) is True
    elev = read_elevation_km(_template(), dem_path)
    assert float(np.max(np.abs(elev.values))) == 0.0
    assert elev.shape == (20, 20)


def test_use_sea_level_elevation_false_for_real_path() -> None:
    assert use_sea_level_elevation("/vsicurl/https://example/dem.vrt") is False
    assert use_sea_level_elevation("dem.tif") is False


def test_read_elevation_km_unreadable_path_falls_back_to_sea_level() -> None:
    # A bad path must degrade to sea level, not abort the run.
    elev = read_elevation_km(_template(), "/no/such/dem.tif")
    assert float(np.max(np.abs(elev.values))) == 0.0
    assert elev.shape == (20, 20)
