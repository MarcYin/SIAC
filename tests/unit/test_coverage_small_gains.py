"""Small focused tests to close remaining coverage gaps."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

from siac.io.geometry import (
    geo_to_pixel_coords,
    get_raster_bounds,
    get_raster_center,
    get_raster_footprint,
    load_aoi,
    pixel_coords_to_geo,
    raster_to_geojson_feature,
)
from siac.io.s2_data_source import S2DataAccess, S2Product, S2Query, _parse_query

if TYPE_CHECKING:
    from pathlib import Path


def _da(shape: tuple[int, int] = (10, 10), crs: str = "EPSG:32632") -> xr.DataArray:
    data = np.arange(shape[0] * shape[1], dtype=np.float32).reshape(shape)
    x = np.linspace(500000.0, 500000.0 + (shape[1] - 1) * 10.0, shape[1])
    y = np.linspace(4500000.0 + (shape[0] - 1) * 10.0, 4500000.0, shape[0])
    return xr.DataArray(data, dims=["y", "x"], coords={"y": y, "x": x}).rio.write_crs(crs)


class TestGeometryCoverage:
    def test_load_aoi_errors_and_transform(self):
        with pytest.raises(ValueError, match="Could not parse AOI"):
            load_aoi("not-an-aoi")

        with pytest.raises(ValueError, match="Empty FeatureCollection"):
            load_aoi({"type": "FeatureCollection", "features": []})

        geom = {
            "type": "Polygon",
            "coordinates": [[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]],
            "crs": {"type": "name", "properties": {"name": "EPSG:4326"}},
        }
        transformed = load_aoi(geom, target_crs="EPSG:3857")
        assert transformed["type"] == "Polygon"

    def test_raster_helpers_and_pixel_geo_roundtrip(self):
        da = _da()
        b4326 = get_raster_bounds(da, crs="EPSG:4326")
        assert b4326[0] < b4326[2]
        assert b4326[1] < b4326[3]

        fp = get_raster_footprint(da, crs="EPSG:4326")
        center = get_raster_center(da, crs="EPSG:4326")
        feat = raster_to_geojson_feature(da, properties={"k": "v"})
        assert fp["type"] == "Polygon"
        assert len(center) == 2
        assert feat["properties"]["k"] == "v"

        tr = da.rio.transform()
        rows = np.array([0, 5, 9])
        cols = np.array([0, 5, 9])
        x, y = pixel_coords_to_geo(rows, cols, tr)
        rr, cc = geo_to_pixel_coords(x, y, tr)
        np.testing.assert_array_equal(rr, rows)
        np.testing.assert_array_equal(cc, cols)


@dataclass
class _FakeBackend:
    products: list[S2Product]
    downloaded: Path | None = None

    def search(self, query: S2Query) -> list[S2Product]:
        return self.products

    def download(self, product: S2Product, dest_dir: Path) -> Path:
        self.downloaded = dest_dir / f"{product.product_id}.SAFE"
        return self.downloaded


class TestS2DataAccessCoverage:
    def test_parse_query_variants(self):
        q1 = _parse_query("S2A_MSIL1C_20240101T100000_N0500_R123_T31UDQ_20240101T120000.SAFE")
        assert q1.product_id is not None

        q2 = _parse_query("T31UDQ_20240101")
        assert q2.mgrs_tile == "31UDQ"

        q3 = _parse_query("something-random")
        assert q3.product_id == "something-random"

    def test_s2_data_access_get_and_search(self, tmp_path: Path):
        p_old = S2Product(
            product_id="P_OLD",
            mgrs_tile="31UDQ",
            sensing_date=datetime(2024, 1, 1, 10),
            processing_baseline="N0400",
            cloud_cover=10.0,
            satellite="S2A",
            orbit_number=1,
            source_url="s3://x",
        )
        p_new = S2Product(
            product_id="P_NEW",
            mgrs_tile="31UDQ",
            sensing_date=datetime(2024, 1, 1, 10),
            processing_baseline="N0500",
            cloud_cover=9.0,
            satellite="S2A",
            orbit_number=1,
            source_url="s3://y",
        )
        backend = _FakeBackend([p_old, p_new])
        access = S2DataAccess(backend=backend, cache_dir=tmp_path / "cache")

        found = access.search(S2Query(mgrs_tile="31UDQ"))
        assert len(found) == 1
        assert found[0].product_id == "P_NEW"

        out = access.get(S2Query(mgrs_tile="31UDQ"), dest_dir=tmp_path)
        assert out.name.endswith(".SAFE")
