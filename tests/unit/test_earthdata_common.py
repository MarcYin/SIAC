"""Tests for shared Earthdata grid helpers."""

from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING

import h5py
import numpy as np

from siac.io.reprojection import transform_bounds
from siac.priors.atmospheric.mcd19_earthaccess import VNP19AODProvider
from siac.priors.earthdata_common import (
    granule_intersects_bounds,
    granule_native_bounds,
    make_native_grid_dataarray,
)

if TYPE_CHECKING:
    from pathlib import Path


def _write_viirs_like_hdf5(path: Path) -> tuple[np.ndarray, np.ndarray]:
    x = np.array([3320000.0, 3360000.0, 3400000.0, 3440000.0, 3480000.0], dtype=np.float64)
    y = np.array([4155000.0, 4125000.0, 4095000.0, 4065000.0], dtype=np.float64)
    struct_metadata = """
GROUP=GridStructure
  GROUP=GRID_1
    GridName="grid750m"
    XDim=5
    YDim=4
    UpperLeftPointMtrs=(3300000.000000,4170000.000000)
    LowerRightMtrs=(3500000.000000,4050000.000000)
    Projection=HE5_GCTP_SNSOID
    ProjParams=(6371007.181000,0,0,0,-160000000,0,0,0,0,0,0,0,0)
  END_GROUP=GRID_1
END_GROUP=GridStructure
"""

    with h5py.File(path, "w") as handle:
        handle.attrs["WestBoundingCoord"] = -131.1325
        handle.attrs["EastBoundingCoord"] = -114.3107
        handle.attrs["SouthBoundingCoord"] = 30.0
        handle.attrs["NorthBoundingCoord"] = 40.0

        grid = handle.create_group("HDFEOS/GRIDS/grid750m")
        grid.create_dataset("XDim", data=x)
        grid.create_dataset("YDim", data=y)
        data_fields = handle.create_group("HDFEOS/GRIDS/grid750m/Data Fields")
        projection = data_fields.create_dataset("Projection", data=np.array([0], dtype=np.int32))
        projection.attrs["earth_radius"] = np.array([6371007.181], dtype=np.float64)
        projection.attrs["false_easting"] = np.array([0.0], dtype=np.float64)
        projection.attrs["false_northing"] = np.array([0.0], dtype=np.float64)
        projection.attrs["grid_mapping_name"] = "sinusoidal"
        projection.attrs["longitude_of_central_meridian"] = np.array([0.0], dtype=np.float64)

        info = handle.create_group("HDFEOS INFORMATION")
        info.create_dataset("StructMetadata.0", data=np.bytes_(struct_metadata))

    return x, y


def test_make_native_grid_dataarray_uses_hdf5_grid_metadata(tmp_path: Path):
    granule = tmp_path / "VNP19A2.A2024165.h25v05.002.fake.h5"
    x, y = _write_viirs_like_hdf5(granule)

    da = make_native_grid_dataarray(np.ones((4, 5), dtype=np.float32), granule_path=granule)

    assert np.allclose(da.x.values, x)
    assert np.allclose(da.y.values, y)

    california_bounds = (-122.0, 36.5, -121.0, 37.5)
    projected = transform_bounds(california_bounds, "EPSG:4326", da.rio.crs)
    xres = float(x[1] - x[0])
    yres = float(y[0] - y[1])
    assert projected[0] > 0.0
    assert float(x.min() - xres / 2.0) < projected[0] < float(x.max() + xres / 2.0)
    assert float(y.min() - yres / 2.0) < projected[1] < float(y.max() + yres / 2.0)


def test_granule_native_bounds_use_hdf5_grid_metadata(tmp_path: Path):
    granule = tmp_path / "VNP19A2.A2024165.h25v05.002.fake.h5"
    _write_viirs_like_hdf5(granule)

    bounds, crs = granule_native_bounds(granule)

    assert bounds == (3300000.0, 4050000.0, 3500000.0, 4170000.0)
    projected = transform_bounds((-122.0, 36.5, -121.0, 37.5), "EPSG:4326", crs)
    assert bounds[0] < projected[0] < bounds[2]
    assert bounds[1] < projected[1] < bounds[3]


def test_vnp19_select_candidate_paths_uses_hdf5_geolocation(tmp_path: Path):
    granule = tmp_path / "VNP19A2.A2024165.h25v05.002.fake.h5"
    _write_viirs_like_hdf5(granule)

    assert granule_intersects_bounds(
        granule,
        bounds=(-122.0, 36.5, -121.0, 37.5),
        crs="EPSG:4326",
    )

    selected = VNP19AODProvider._select_candidate_paths(
        [granule],
        obs_time=datetime(2024, 6, 15),
        bounds=(-122.0, 36.5, -121.0, 37.5),
        crs="EPSG:4326",
    )

    assert selected == [granule]
