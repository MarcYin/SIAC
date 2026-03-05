"""
Unit tests for SIAC I/O module.
"""

import json
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from siac.io.geometry import (
    bounds_area,
    bounds_contains,
    bounds_intersect,
    bounds_intersection,
    bounds_to_polygon,
    bounds_union,
    create_grid_from_bounds,
    load_aoi,
    polygon_to_bounds,
)
from siac.io.reprojection import (
    compute_common_bounds,
    get_bounds,
    get_resolution,
)

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def sample_dataarray() -> xr.DataArray:
    """Create a sample DataArray with spatial coordinates."""
    data = np.random.default_rng(0).random((100, 100)).astype(np.float32)

    # Create coordinates (10m resolution, UTM-like)
    x = np.linspace(500000, 500990, 100)
    y = np.linspace(4500990, 4500000, 100)

    da = xr.DataArray(
        data,
        dims=["y", "x"],
        coords={"y": y, "x": x},
    )

    # Set CRS and transform
    da = da.rio.write_crs("EPSG:32632")
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")

    return da


@pytest.fixture
def sample_dataset(sample_dataarray: xr.DataArray) -> xr.Dataset:
    """Create a sample Dataset with multiple bands."""
    ds = xr.Dataset(
        {
            "B02": sample_dataarray,
            "B03": sample_dataarray * 1.1,
            "B04": sample_dataarray * 0.9,
        }
    )
    return ds


@pytest.fixture
def sample_geojson(tmp_path: Path) -> Path:
    """Create a sample GeoJSON file."""
    geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [500000, 4500000],
                            [500500, 4500000],
                            [500500, 4500500],
                            [500000, 4500500],
                            [500000, 4500000],
                        ]
                    ],
                },
            }
        ],
    }

    path = tmp_path / "aoi.geojson"
    with path.open("w") as f:
        json.dump(geojson, f)

    return path


# =============================================================================
# Geometry Tests
# =============================================================================


class TestBoundsOperations:
    """Tests for bounding box operations."""

    def test_bounds_to_polygon(self):
        """bounds_to_polygon should create valid GeoJSON."""
        bounds = (0.0, 0.0, 100.0, 100.0)
        polygon = bounds_to_polygon(bounds)

        assert polygon["type"] == "Polygon"
        assert len(polygon["coordinates"]) == 1
        assert len(polygon["coordinates"][0]) == 5  # Closed ring

    def test_polygon_to_bounds(self):
        """polygon_to_bounds should extract correct bounds."""
        polygon = {
            "type": "Polygon",
            "coordinates": [[[0, 0], [100, 0], [100, 100], [0, 100], [0, 0]]],
        }
        bounds = polygon_to_bounds(polygon)

        assert bounds == (0, 0, 100, 100)

    def test_bounds_roundtrip(self):
        """bounds -> polygon -> bounds should be identity."""
        original = (10.0, 20.0, 30.0, 40.0)
        polygon = bounds_to_polygon(original)
        recovered = polygon_to_bounds(polygon)

        assert recovered == original

    def test_bounds_intersect_true(self):
        """bounds_intersect should return True for overlapping bounds."""
        bounds1 = (0, 0, 100, 100)
        bounds2 = (50, 50, 150, 150)

        assert bounds_intersect(bounds1, bounds2)

    def test_bounds_intersect_false(self):
        """bounds_intersect should return False for non-overlapping bounds."""
        bounds1 = (0, 0, 100, 100)
        bounds2 = (200, 200, 300, 300)

        assert not bounds_intersect(bounds1, bounds2)

    def test_bounds_intersect_adjacent(self):
        """bounds_intersect should return False for adjacent bounds."""
        bounds1 = (0, 0, 100, 100)
        bounds2 = (100, 0, 200, 100)  # Touching at x=100

        assert not bounds_intersect(bounds1, bounds2)

    def test_bounds_contains_true(self):
        """bounds_contains should return True when outer contains inner."""
        outer = (0, 0, 100, 100)
        inner = (25, 25, 75, 75)

        assert bounds_contains(outer, inner)

    def test_bounds_contains_false(self):
        """bounds_contains should return False when not contained."""
        outer = (0, 0, 100, 100)
        inner = (50, 50, 150, 150)  # Partially outside

        assert not bounds_contains(outer, inner)

    def test_bounds_intersection(self):
        """bounds_intersection should compute correct overlap."""
        bounds1 = (0, 0, 100, 100)
        bounds2 = (50, 50, 150, 150)

        intersection = bounds_intersection(bounds1, bounds2)

        assert intersection == (50, 50, 100, 100)

    def test_bounds_intersection_none(self):
        """bounds_intersection should return None for non-overlapping."""
        bounds1 = (0, 0, 100, 100)
        bounds2 = (200, 200, 300, 300)

        assert bounds_intersection(bounds1, bounds2) is None

    def test_bounds_union(self):
        """bounds_union should compute correct combined extent."""
        bounds1 = (0, 0, 100, 100)
        bounds2 = (50, 50, 150, 150)

        union = bounds_union(bounds1, bounds2)

        assert union == (0, 0, 150, 150)

    def test_bounds_area(self):
        """bounds_area should compute correct area."""
        bounds = (0, 0, 100, 50)

        assert bounds_area(bounds) == 5000.0


class TestLoadAOI:
    """Tests for AOI loading."""

    def test_load_aoi_from_bounds(self):
        """load_aoi should accept bounding box tuple."""
        bounds = (500000, 4500000, 500500, 4500500)
        aoi = load_aoi(bounds)

        assert aoi["type"] == "Polygon"
        recovered_bounds = polygon_to_bounds(aoi)
        assert recovered_bounds == bounds

    def test_load_aoi_from_list(self):
        """load_aoi should accept bounding box list."""
        bounds = [500000, 4500000, 500500, 4500500]
        aoi = load_aoi(bounds)

        assert aoi["type"] == "Polygon"

    def test_load_aoi_from_geojson_file(self, sample_geojson: Path):
        """load_aoi should load from GeoJSON file."""
        aoi = load_aoi(sample_geojson)

        assert aoi["type"] == "Polygon"
        assert len(aoi["coordinates"]) == 1

    def test_load_aoi_from_dict(self):
        """load_aoi should accept GeoJSON dict."""
        geojson = {
            "type": "Polygon",
            "coordinates": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]],
        }
        aoi = load_aoi(geojson)

        assert aoi["type"] == "Polygon"

    def test_load_aoi_from_wkt(self):
        """load_aoi should parse WKT string."""
        wkt = "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"
        aoi = load_aoi(wkt)

        assert aoi["type"] == "Polygon"

    def test_load_aoi_invalid_type(self):
        """load_aoi should raise for invalid input type."""
        with pytest.raises(ValueError):
            load_aoi(12345)


class TestCreateGrid:
    """Tests for grid creation."""

    def test_create_grid_from_bounds(self):
        """create_grid_from_bounds should create valid DataArray."""
        bounds = (500000, 4500000, 501000, 4501000)
        resolution = 10.0
        crs = "EPSG:32632"

        da = create_grid_from_bounds(bounds, resolution, crs)

        assert da.rio.crs.to_string() == crs
        assert da.shape == (100, 100)

        # Check resolution
        res = get_resolution(da)
        np.testing.assert_allclose(res, (10.0, 10.0), rtol=0.01)


# =============================================================================
# Reprojection Tests
# =============================================================================


class TestReprojectionUtils:
    """Tests for reprojection utilities."""

    def test_get_resolution(self, sample_dataarray: xr.DataArray):
        """get_resolution should return correct values."""
        res = get_resolution(sample_dataarray)

        # Should be approximately 10m
        np.testing.assert_allclose(res, (10.0, 10.0), rtol=0.01)

    def test_get_bounds(self, sample_dataarray: xr.DataArray):
        """get_bounds should return correct extent."""
        bounds = get_bounds(sample_dataarray)

        # Check approximate values
        assert bounds[0] < bounds[2]  # xmin < xmax
        assert bounds[1] < bounds[3]  # ymin < ymax

    def test_compute_common_bounds_intersection(
        self, sample_dataarray: xr.DataArray
    ):
        """compute_common_bounds should find intersection."""
        # Create two overlapping arrays
        da1 = sample_dataarray
        da2 = sample_dataarray.copy()

        # Shift da2 by 500m
        da2 = da2.assign_coords(
            x=da2.x + 500,
            y=da2.y - 500,
        )
        da2 = da2.rio.write_crs("EPSG:32632")

        bounds = compute_common_bounds(da1, da2, method="intersection")

        # Should have valid bounds
        assert bounds[0] < bounds[2]
        assert bounds[1] < bounds[3]

    def test_compute_common_bounds_union(self, sample_dataarray: xr.DataArray):
        """compute_common_bounds union should expand extent."""
        da1 = sample_dataarray
        bounds1 = get_bounds(da1)

        # Create shifted array
        da2 = sample_dataarray.copy()
        da2 = da2.assign_coords(
            x=da2.x + 500,
        )
        da2 = da2.rio.write_crs("EPSG:32632")
        bounds2 = get_bounds(da2)

        union_bounds = compute_common_bounds(da1, da2, method="union")

        # Union should be larger than either
        assert union_bounds[0] <= min(bounds1[0], bounds2[0])
        assert union_bounds[2] >= max(bounds1[2], bounds2[2])


# =============================================================================
# Reader/Writer Integration Tests (require actual files)
# =============================================================================


class TestReadWriteRoundtrip:
    """Integration tests for read/write operations."""

    def test_write_read_geotiff(
        self, sample_dataarray: xr.DataArray, tmp_path: Path
    ):
        """Write and read GeoTIFF should preserve data."""
        from siac.io import read_raster, write_raster

        output_path = tmp_path / "test.tif"

        # Write
        write_raster(sample_dataarray, output_path, compression="deflate")

        assert output_path.exists()

        # Read back
        read_da = read_raster(output_path)

        # Compare shapes
        assert read_da.shape == sample_dataarray.shape

        # Compare values (allowing for float precision)
        np.testing.assert_allclose(
            read_da.values,
            sample_dataarray.values,
            rtol=1e-5,
        )

    def test_write_read_cog(
        self, sample_dataarray: xr.DataArray, tmp_path: Path
    ):
        """Write and read COG should preserve data."""
        from siac.io import read_raster, write_cog

        output_path = tmp_path / "test_cog.tif"

        # Write as COG
        write_cog(sample_dataarray, output_path)

        assert output_path.exists()

        # Read back
        read_da = read_raster(output_path)

        # Compare
        np.testing.assert_allclose(
            read_da.values,
            sample_dataarray.values,
            rtol=1e-5,
        )

    def test_write_dataset(self, sample_dataset: xr.Dataset, tmp_path: Path):
        """write_dataset should create files for each variable."""
        from siac.io import write_dataset

        output_dir = tmp_path / "output"

        paths = write_dataset(sample_dataset, output_dir, prefix="test_")

        # Check all bands written
        assert len(paths) == 3
        for var_name, path in paths.items():
            assert path.exists()
            assert var_name in str(path)

    def test_write_zarr_roundtrip(
        self, sample_dataarray: xr.DataArray, tmp_path: Path
    ):
        """Write and read Zarr should preserve data."""
        from siac.io import read_zarr_array, write_zarr

        output_path = tmp_path / "test.zarr"

        # Need to convert to dataset for zarr
        ds = sample_dataarray.to_dataset(name="data")

        # Write
        write_zarr(ds, output_path)

        assert output_path.exists()

        # Read back
        read_ds = read_zarr_array(output_path)

        # Compare
        np.testing.assert_allclose(
            read_ds["data"].values,
            sample_dataarray.values,
            rtol=1e-5,
        )
