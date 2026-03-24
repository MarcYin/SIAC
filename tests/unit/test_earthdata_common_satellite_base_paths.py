"""Focused coverage lifts for Earthdata common helpers and satellite base."""

from __future__ import annotations

import importlib.util
import sys
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pytest
import xarray as xr

import siac.adapters.earthdata as earth_adapter
import siac.adapters.earthdata_common as earth_mod

if TYPE_CHECKING:
    from types import ModuleType


def _load_satellite_base_module() -> ModuleType:
    path = Path("/Users/fengyin/Documents/SIAC/python/siac/adapters/satellite/base.py")
    spec = importlib.util.spec_from_file_location("siac_satellite_base_test_module", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_earthdata_scalar_and_parse_helpers_cover_fallbacks() -> None:
    assert earth_mod.decode_attr(np.array([np.bytes_(b"abc")])) == "abc"
    assert earth_mod.decode_attr(np.array([1, 2], dtype=np.int32)) == [1, 2]
    assert earth_mod.decode_attr(np.float32(1.5)) == pytest.approx(1.5)
    assert earth_mod.decode_attr(np.bytes_(b"xyz")) == "xyz"

    assert earth_mod.attr_scalar({}, "missing", default=9) == 9
    assert earth_mod.attr_scalar({"empty": np.array([], dtype=np.int32)}, "empty", default=7) == 7
    assert earth_mod.attr_scalar({"none": None}, "none", default=5) == 5
    assert earth_mod.attr_scalar({"ints": np.array([4, 5], dtype=np.int32)}, "ints") == 4
    assert earth_mod.attr_scalar({"floats": np.array([2.5, 3.5], dtype=np.float32)}, "floats") == 2.5

    assert earth_mod.parse_tile_indices("MCD43A1.A2024001.h12v04.061.hdf") == (12, 4)
    with pytest.raises(ValueError, match="tile indices"):
        earth_mod.parse_tile_indices("bad-name.hdf")

    assert earth_mod.parse_granule_date("MCD43A1.A2024366.h12v04.061.hdf") == datetime(2024, 12, 31)
    with pytest.raises(ValueError, match="acquisition date"):
        earth_mod.parse_granule_date("bad-name.hdf")

    assert earth_mod._decode_gctp_angle(120.0) == 120.0
    assert earth_mod._decode_gctp_angle(-160000000.0) == -160.0
    assert "+lon_0=-160.0" in earth_mod._build_sinusoidal_crs(radius=6371007.181, central_meridian=-160.0)


def test_earthdata_metadata_grid_and_bounds_helpers_cover_error_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    earth_mod._parse_hdf5_grid_metadata.cache_clear()
    monkeypatch.setattr(earth_mod, "_read_hdf5_struct_metadata", lambda _path: (_ for _ in ()).throw(RuntimeError("boom")))
    assert earth_mod._parse_hdf5_grid_metadata("broken.h5") == {}

    struct_metadata = """
GROUP=GridStructure
  GROUP=GRID_1
    GridName="grid750m"
    XDim=2
    YDim=3
    UpperLeftPointMtrs=(10.0,20.0)
    LowerRightMtrs=(30.0,-10.0)
    Projection=HE5_GCTP_SNSOID
    ProjParams=(6371007.181000,0,0,0,-160000000,0,0,0,0,0,0,0,0)
  END_GROUP=GRID_1
END_GROUP=GridStructure
"""
    earth_mod._parse_hdf5_grid_metadata.cache_clear()
    monkeypatch.setattr(earth_mod, "_read_hdf5_struct_metadata", lambda _path: struct_metadata)
    parsed = earth_mod._parse_hdf5_grid_metadata("grid.h5")
    assert parsed["grid750m"]["xdim"] == 2
    assert parsed["grid750m"]["ydim"] == 3
    assert "+lon_0=-160.0" in parsed["grid750m"]["crs"]

    grid = earth_mod._grid_definition_from_extent(
        width=2,
        height=3,
        upper_left=(10.0, 20.0),
        lower_right=(30.0, -10.0),
        crs="EPSG:4326",
    )
    assert grid.x == (15.0, 25.0)
    assert grid.y == (15.0, 5.0, -5.0)

    rootless = tmp_path / "rootless.h5"
    with h5py.File(rootless, "w"):
        pass
    earth_mod._read_hdf5_root_bounds.cache_clear()
    assert earth_mod._read_hdf5_root_bounds(str(rootless)) is None
    assert earth_mod._read_hdf5_root_bounds(str(tmp_path / "missing.h5")) is None

    bad_grid_file = tmp_path / "VNP19A2.A2024165.h25v05.fake.h5"
    with h5py.File(bad_grid_file, "w") as handle:
        handle.create_group("HDFEOS/GRIDS")
    assert earth_mod._read_hdf5_native_grid_definition(bad_grid_file) is None

    one_pixel = tmp_path / "VNP19A2.A2024165.h26v06.fake.h5"
    with h5py.File(one_pixel, "w") as handle:
        grid_group = handle.create_group("HDFEOS/GRIDS/grid")
        grid_group.create_dataset("XDim", data=np.array([1000.0], dtype=np.float64))
        grid_group.create_dataset("YDim", data=np.array([2000.0], dtype=np.float64))
    bounds, crs = earth_mod.granule_native_bounds(one_pixel, height=1, width=1)
    assert bounds == (500.0, 1500.0, 1500.0, 2500.0)
    assert crs == earth_mod.MODLAND_SINUSOIDAL_CRS


def test_earthdata_native_grid_reproject_and_dataset_helpers(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    assert earth_mod.granule_geographic_bounds("scene.hdf") is None

    fallback_path = tmp_path / "MCD43A1.A2024001.h10v04.061.hdf"
    fallback_bounds, fallback_crs = earth_mod.granule_native_bounds(fallback_path)
    assert fallback_bounds == earth_mod.modland_tile_bounds(10, 4)
    assert fallback_crs == earth_mod.MODLAND_SINUSOIDAL_CRS

    with pytest.raises(ValueError, match="positive"):
        earth_mod.modland_tile_coords(10, 4, 0, 2)

    with pytest.raises(ValueError, match="at least 2 dimensions"):
        earth_mod.make_native_grid_dataarray(np.array([1.0], dtype=np.float32), granule_path=fallback_path)

    monkeypatch.setattr(earth_mod, "_read_hdf5_native_grid_definition", lambda *_args, **_kwargs: None)
    da = earth_mod.make_native_grid_dataarray(
        np.arange(12, dtype=np.float32).reshape(2, 2, 3),
        granule_path=fallback_path,
        dims=("band", "y", "x"),
        coords={"band": ["a", "b"]},
    )
    assert da.dims == ("band", "y", "x")
    assert list(da.coords["band"].values) == ["a", "b"]
    assert da.rio.crs is not None

    with pytest.raises(ValueError, match="resolution must be > 0"):
        earth_mod.build_target_template((0.0, 0.0, 1.0, 1.0), earth_mod.MODLAND_SINUSOIDAL_CRS, 0.0)

    target = earth_mod.build_target_template(
        (0.0, 0.0, 3000.0, 2000.0),
        earth_mod.MODLAND_SINUSOIDAL_CRS,
        1000.0,
        fill_value=-5.0,
    )
    assert target.shape == (2, 3)
    assert float(target.values[0, 0]) == -5.0

    clip_source = earth_mod.make_native_grid_dataarray(
        np.arange(9, dtype=np.float32).reshape(3, 3),
        granule_path=fallback_path,
    )
    clip_bounds = (
        float(clip_source.x.values[1] - 10.0),
        float(clip_source.y.values[1] - 10.0),
        float(clip_source.x.values[1] + 10.0),
        float(clip_source.y.values[1] + 10.0),
    )
    monkeypatch.setattr(earth_mod, "transform_bounds", lambda *_args, **_kwargs: clip_bounds)
    clipped = earth_mod.clip_native_to_target_bounds(
        clip_source,
        target_bounds=(0.0, 0.0, 1.0, 1.0),
        target_crs="EPSG:4326",
        pad_pixels=-2,
    )
    assert clipped.shape == (1, 1)

    class _FakeClipped:
        def __init__(self) -> None:
            self.rio = _FakeRio(self)

    class _FakeRio:
        def __init__(self, parent: _FakeClipped) -> None:
            self.parent = parent
            self.nodata = None

        def write_nodata(self, nodata: float):
            self.nodata = nodata
            return self.parent

        def reproject_match(self, target_da: xr.DataArray, resampling: object):
            return {"target_shape": target_da.shape, "resampling": resampling, "nodata": self.nodata}

    fake_clipped = _FakeClipped()
    monkeypatch.setattr(earth_mod, "clip_native_to_target_bounds", lambda *_args, **_kwargs: fake_clipped)
    reprojected = earth_mod.reproject_native_to_target(
        xr.DataArray(np.zeros((2, 2), dtype=np.float32), dims=("y", "x")),
        target_bounds=(0.0, 0.0, 2000.0, 2000.0),
        target_crs=earth_mod.MODLAND_SINUSOIDAL_CRS,
        target_resolution=1000.0,
        resampling="nearest",
        nodata=-9999.0,
    )
    assert reprojected["target_shape"] == (2, 2)
    assert reprojected["resampling"] == "nearest"
    assert reprojected["nodata"] == -9999.0

    h5_path = tmp_path / "dataset.h5"
    with h5py.File(h5_path, "w") as handle:
        ds = handle.create_dataset("data", data=np.array([[1, 2], [3, 4]], dtype=np.int16))
        ds.attrs["units"] = np.bytes_(b"reflectance")
    values, attrs = earth_mod.read_hdf5_dataset(h5_path, "data")
    assert values.shape == (2, 2)
    assert attrs["units"] == "reflectance"

    scaled = earth_mod.apply_scale_and_mask(
        np.array([[0, 5, 12], [99, 6, 7]], dtype=np.float32),
        {
            "_FillValue": np.array([99], dtype=np.int16),
            "valid_range": np.array([1, 10], dtype=np.int16),
            "scale_factor": np.array([0.5], dtype=np.float32),
            "add_offset": np.array([1.0], dtype=np.float32),
        },
    )
    assert np.isnan(scaled[0, 0])
    assert scaled[0, 1] == pytest.approx(3.5)
    assert np.isnan(scaled[0, 2])
    assert np.isnan(scaled[1, 0])

    reduced = earth_mod.reduce_orbit_stack(np.array([[[1.0, np.nan], [3.0, 5.0]], [[3.0, 7.0], [np.nan, 5.0]]], dtype=np.float32))
    assert reduced.shape == (2, 2)
    assert reduced[0, 0] == pytest.approx(2.0)
    assert reduced[0, 1] == pytest.approx(7.0)
    unchanged = earth_mod.reduce_orbit_stack(np.array([[1.0, 2.0]], dtype=np.float32))
    assert unchanged.shape == (1, 2)


def test_earthdata_attribute_readers_and_gdal_subdataset_resolution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class _FakeSDS:
        def __init__(self, name: str) -> None:
            self.name = name

        def attributes(self) -> dict[str, object]:
            return {"units": np.bytes_(self.name.encode("utf-8"))}

    class _FakeSD:
        def __init__(self, path: str, mode: object) -> None:
            self.path = path
            self.mode = mode

        def select(self, dataset_name: str) -> _FakeSDS:
            return _FakeSDS(dataset_name)

    monkeypatch.setattr(earth_mod, "SD", _FakeSD)
    monkeypatch.setattr(earth_mod, "SDC", SimpleNamespace(READ="READ"))
    assert earth_mod.read_hdf4_dataset_attrs("scene.hdf", "band") == {"units": "band"}
    assert earth_mod.read_hdf4_datasets_attrs("scene.hdf", ("band", "qa")) == {
        "band": {"units": "band"},
        "qa": {"units": "qa"},
    }

    h5_path = tmp_path / "attrs.h5"
    with h5py.File(h5_path, "w") as handle:
        band = handle.create_dataset("band", data=np.array([[1]], dtype=np.int16))
        band.attrs["scale_factor"] = np.array([2.0], dtype=np.float32)
        qa = handle.create_dataset("qa", data=np.array([[0]], dtype=np.int16))
        qa.attrs["units"] = np.bytes_(b"bits")

    assert earth_mod.read_hdf5_dataset_attrs(h5_path, "band")["scale_factor"] == 2.0
    assert earth_mod.read_hdf5_datasets_attrs(h5_path, ("band", "qa")) == {
        "band": {"scale_factor": 2.0},
        "qa": {"units": "bits"},
    }

    earth_mod.resolve_gdal_subdataset_path.cache_clear()
    fake_gdal = SimpleNamespace(
        UseExceptions=lambda: None,
        OF_RASTER=1,
        OF_READONLY=2,
        OpenEx=lambda _path, _flags: SimpleNamespace(
            GetSubDatasets=lambda: [
                (
                    'HDF4_EOS:EOS_GRID:"scene.hdf":MOD_Grid_BRDF:BRDF_Albedo_Parameters_Band1',
                    "[2400x2400x3] BRDF_Albedo_Parameters_Band1",
                ),
                (
                    'HDF4_SDS:UNKNOWN:"scene.hdf":7',
                    "[2400x2400] Status_QA",
                ),
            ]
        ),
    )
    monkeypatch.setitem(sys.modules, "osgeo", SimpleNamespace(gdal=fake_gdal))
    assert earth_mod.gdal_available()
    assert earth_mod.resolve_gdal_subdataset_path("scene.hdf", "BRDF_Albedo_Parameters_Band1") == (
        'HDF4_EOS:EOS_GRID:"scene.hdf":MOD_Grid_BRDF:BRDF_Albedo_Parameters_Band1'
    )
    assert earth_mod.resolve_gdal_subdataset_path("scene.hdf", "Status_QA") == 'HDF4_SDS:UNKNOWN:"scene.hdf":7'

    earth_mod.resolve_gdal_subdataset_path.cache_clear()
    fake_empty_gdal = SimpleNamespace(
        UseExceptions=lambda: None,
        OF_RASTER=1,
        OF_READONLY=2,
        OpenEx=lambda _path, _flags: SimpleNamespace(GetSubDatasets=lambda: []),
    )
    monkeypatch.setitem(sys.modules, "osgeo", SimpleNamespace(gdal=fake_empty_gdal))
    assert earth_mod.resolve_gdal_subdataset_path(
        "scene.h5",
        "HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M1",
    ) == 'HDF5:"scene.h5"://HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/BRDF_Albedo_Parameters_M1'


def test_earthdata_intersection_and_hdf4_reader_branches(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(earth_mod, "granule_geographic_bounds", lambda _path: (-10.0, -10.0, -5.0, -5.0))
    monkeypatch.setattr(earth_mod, "transform_bounds", lambda bounds, src, dst: bounds if src == dst else (0.0, 0.0, 1.0, 1.0))
    assert not earth_mod.granule_intersects_bounds("scene.h5", bounds=(2.0, 2.0, 3.0, 3.0), crs="EPSG:4326")

    monkeypatch.setattr(earth_mod, "granule_geographic_bounds", lambda _path: None)
    monkeypatch.setattr(earth_mod, "granule_native_bounds", lambda _path: ((0.0, 0.0, 1.0, 1.0), earth_mod.MODLAND_SINUSOIDAL_CRS))
    monkeypatch.setattr(earth_mod, "transform_bounds", lambda *_args, **_kwargs: (2.0, 2.0, 3.0, 3.0))
    assert not earth_mod.granule_intersects_bounds("scene.hdf", bounds=(0.0, 0.0, 1.0, 1.0), crs="EPSG:4326")

    class _FakeSDS:
        def get(self) -> np.ndarray:
            return np.array([[1, 2], [3, 4]], dtype=np.int16)

        def attributes(self) -> dict[str, object]:
            return {"scale_factor": np.array([2.0], dtype=np.float32)}

    class _FakeSD:
        def __init__(self, path: str, mode: object) -> None:
            self.path = path
            self.mode = mode

        def select(self, dataset_name: str) -> _FakeSDS:
            assert dataset_name == "band"
            return _FakeSDS()

    monkeypatch.setattr(earth_mod, "SD", _FakeSD)
    monkeypatch.setattr(earth_mod, "SDC", SimpleNamespace(READ="READ"))
    values, attrs = earth_mod.read_hdf4_dataset("scene.hdf", "band")
    assert values.shape == (2, 2)
    assert attrs["scale_factor"] == 2.0


def test_read_virtual_stack_to_target_builds_group_vrts_then_one_stack(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    build_calls: list[tuple[str, list[str], object | None]] = []
    unlink_calls: list[str] = []

    class _FakeDataset:
        def __init__(self, values: np.ndarray) -> None:
            self._values = values

        def ReadAsArray(self) -> np.ndarray:
            return self._values

    def _build_vrt(path: str, sources: list[str], options: object | None = None) -> object:
        build_calls.append((path, list(sources), options))
        return SimpleNamespace(path=path)

    fake_gdal = SimpleNamespace(
        UseExceptions=lambda: None,
        GDT_Float32=6,
        BuildVRT=_build_vrt,
        BuildVRTOptions=lambda **kwargs: kwargs,
        Warp=lambda _dest, _src, **_kwargs: _FakeDataset(
            np.arange(24, dtype=np.float32).reshape(4, 2, 3)
        ),
        Unlink=lambda path: unlink_calls.append(path),
    )
    monkeypatch.setitem(sys.modules, "osgeo", SimpleNamespace(gdal=fake_gdal))

    out = earth_adapter.read_virtual_stack_to_target(
        [["left-a.tif", "right-a.tif"], ["left-b.tif"]],
        group_band_counts=[3, 1],
        bounds=(0.0, 0.0, 3.0, 2.0),
        crs="EPSG:4326",
        resolution=1.0,
        resampling="nearest",
        nodata=np.nan,
    )

    assert out.dims == ("layer", "y", "x")
    assert out.shape == (4, 2, 3)
    assert build_calls[0][1] == ["left-a.tif", "right-a.tif"]
    assert build_calls[1][1] == ["left-b.tif"]
    assert build_calls[2][1] == [build_calls[0][0], build_calls[1][0]]
    assert build_calls[2][2] == {"separate": True}
    assert len(unlink_calls) == 3


def test_merge_reprojected_tiles_handles_single_pixel_target_and_irregular_grid() -> None:
    ref = xr.DataArray(
        np.array([[np.nan, 1.0], [1.0, 1.0]], dtype=np.float32),
        dims=("y", "x"),
        coords={"y": [1.0, 0.0], "x": [0.0, 1.0]},
    )
    cand = xr.DataArray(
        np.array([[5.0, 6.0], [7.0, 8.0]], dtype=np.float32),
        dims=("y", "x"),
        coords={"y": [1.0, 0.0], "x": [0.0, 2.0]},
    )
    ref = ref.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs("EPSG:4326")
    cand = cand.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs("EPSG:4326")

    out = earth_adapter.merge_reprojected_tiles(
        [ref, cand],
        bounds=(-0.5, 0.5, 0.5, 1.5),
        crs="EPSG:4326",
        resolution=1.0,
        resampling="nearest",
        nodata=None,
    )

    assert tuple(out.sizes[dim] for dim in ("y", "x")) == (1, 1)
    np.testing.assert_allclose(out.coords["x"].values, np.array([0.0]))
    np.testing.assert_allclose(out.coords["y"].values, np.array([1.0]))
    np.testing.assert_allclose(out.values, np.array([[5.0]], dtype=np.float32))
    transform = out.rio.transform()
    assert transform.a == pytest.approx(1.0)
    assert transform.e == pytest.approx(-1.0)
    assert transform.c == pytest.approx(-0.5)
    assert transform.f == pytest.approx(1.5)


def test_merge_reprojected_tiles_handles_mixed_crs_sources() -> None:
    from rasterio.warp import transform

    left = xr.DataArray(
        np.array([[1.0, np.nan], [3.0, 4.0]], dtype=np.float32),
        dims=("y", "x"),
        coords={"y": [1.0, 0.0], "x": [0.0, 1.0]},
    )
    mx, my = transform("EPSG:4326", "EPSG:3857", [1.0, 2.0], [1.0, 0.0])
    right = xr.DataArray(
        np.array([[9.0, 10.0], [11.0, 12.0]], dtype=np.float32),
        dims=("y", "x"),
        coords={"y": my, "x": mx},
    )
    left = left.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs("EPSG:4326")
    right = right.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs("EPSG:3857")

    out = earth_adapter.merge_reprojected_tiles(
        [left, right],
        bounds=(-0.5, -0.5, 2.5, 1.5),
        crs="EPSG:4326",
        resolution=1.0,
        resampling="nearest",
        nodata=None,
    )

    assert out.dims == ("y", "x")
    np.testing.assert_allclose(out.coords["x"].values, np.array([0.0, 1.0, 2.0]))
    np.testing.assert_allclose(out.coords["y"].values, np.array([1.0, 0.0]))
    np.testing.assert_allclose(
        out.values,
        np.array(
            [
                [1.0, 9.0, 10.0],
                [3.0, 4.0, 12.0],
            ],
            dtype=np.float32,
        ),
    )


def test_merge_reprojected_tiles_preserves_priority_and_extra_dims() -> None:
    left = xr.DataArray(
        np.array(
            [
                [[1.0, np.nan], [3.0, 4.0]],
                [[10.0, 11.0], [12.0, np.nan]],
            ],
            dtype=np.float32,
        ),
        dims=("parameter", "y", "x"),
        coords={"parameter": ["f0", "f1"], "y": [1.0, 0.0], "x": [0.0, 1.0]},
    )
    right = xr.DataArray(
        np.array(
            [
                [[9.0, 10.0], [11.0, 12.0]],
                [[19.0, 20.0], [21.0, 22.0]],
            ],
            dtype=np.float32,
        ),
        dims=("parameter", "y", "x"),
        coords={"parameter": ["f0", "f1"], "y": [1.0, 0.0], "x": [1.0, 2.0]},
    )
    left = left.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs("EPSG:4326")
    right = right.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs("EPSG:4326")

    out = earth_adapter.merge_reprojected_tiles(
        [left, right],
        bounds=(-0.5, -0.5, 2.5, 1.5),
        crs="EPSG:4326",
        resolution=1.0,
        resampling="nearest",
        nodata=None,
    )

    assert out.dims == ("parameter", "y", "x")
    assert list(out.coords["parameter"].values) == ["f0", "f1"]
    np.testing.assert_allclose(out.coords["x"].values, np.array([0.0, 1.0, 2.0]))
    np.testing.assert_allclose(out.coords["y"].values, np.array([1.0, 0.0]))
    np.testing.assert_allclose(
        out.values,
        np.array(
            [
                [[1.0, 9.0, 10.0], [3.0, 4.0, 12.0]],
                [[10.0, 11.0, 20.0], [12.0, 21.0, 22.0]],
            ],
            dtype=np.float32,
        ),
    )


def test_merge_reprojected_tiles_passes_full_sources_to_vrt_reader(monkeypatch: pytest.MonkeyPatch) -> None:
    left = xr.DataArray(
        np.arange(100, dtype=np.float32).reshape(10, 10),
        dims=("y", "x"),
        coords={"y": np.arange(9.0, -1.0, -1.0), "x": np.arange(10.0)},
    )
    left = left.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs("EPSG:4326")

    original_writer = earth_adapter._write_dataarray_to_memory_dataset
    seen_shapes: list[tuple[int, int]] = []

    def _recording_writer(memfile: object, data: xr.DataArray, *, nodata: np.float32) -> None:
        seen_shapes.append((int(data.sizes["y"]), int(data.sizes["x"])))
        original_writer(memfile, data, nodata=nodata)

    monkeypatch.setattr(earth_adapter, "_write_dataarray_to_memory_dataset", _recording_writer)

    out = earth_adapter.merge_reprojected_tiles(
        [left],
        bounds=(4.5, 4.5, 5.5, 5.5),
        crs="EPSG:4326",
        resolution=1.0,
        resampling="nearest",
        nodata=None,
    )

    assert seen_shapes == [(10, 10)]
    np.testing.assert_allclose(out.values, np.array([[45.0]], dtype=np.float32))


def test_merge_reprojected_tiles_returns_empty_target_when_sources_miss_aoi() -> None:
    source = xr.DataArray(
        np.array([[5.0, 6.0], [7.0, 8.0]], dtype=np.float32),
        dims=("y", "x"),
        coords={"y": [11.0, 10.0], "x": [10.0, 11.0]},
    )
    source = source.rio.set_spatial_dims(x_dim="x", y_dim="y").rio.write_crs("EPSG:4326")

    out = earth_adapter.merge_reprojected_tiles(
        [source],
        bounds=(-0.5, -0.5, 1.5, 1.5),
        crs="EPSG:4326",
        resolution=1.0,
        resampling="nearest",
        nodata=np.nan,
    )

    assert out.shape == (2, 2)
    assert np.isnan(out.values).all()


def test_satellite_base_preprocessor_registry_and_detection_paths(
    tmp_path: Path,
    mock_sensor_config,
) -> None:
    base_mod = _load_satellite_base_module()

    class _NoToaPreprocessor(base_mod.BaseSatellitePreprocessor):
        def __init__(self) -> None:
            super().__init__()
            self.received_toa = None

        @property
        def sensor_config(self):
            return mock_sensor_config

        def load_toa(self, input_path):
            del input_path
            return xr.Dataset({"B02": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=("y", "x"))})

        def extract_geometry(self, input_path):
            del input_path
            arr = xr.DataArray(np.full((2, 2), 0.1, dtype=np.float32), dims=("y", "x"))
            return SimpleNamespace(sza=arr, saa=arr, vza=arr, vaa=arr)

        def extract_cloud_mask(self, input_path):
            del input_path
            self.received_toa = "missing"
            return xr.DataArray(np.zeros((2, 2), dtype=bool), dims=("y", "x"))

        def get_metadata(self, input_path):
            del input_path
            return {"observation_time": datetime(2026, 1, 2, 3, 4, 5)}

    class _KwargPreprocessor(_NoToaPreprocessor):
        def extract_cloud_mask(self, input_path, **kwargs):
            del input_path
            self.received_toa = kwargs.get("toa")
            return xr.DataArray(np.ones((2, 2), dtype=bool), dims=("y", "x"))

    no_toa = _NoToaPreprocessor()
    raw = no_toa.preprocess(tmp_path / "scene")
    assert no_toa.received_toa == "missing"
    assert "cloud_classes" not in raw

    kwarg = _KwargPreprocessor()
    kwarg._last_cloud_classes = xr.DataArray(np.ones((2, 2), dtype=np.uint8), dims=("y", "x"))
    bundle = kwarg.to_observation_bundle(tmp_path / "scene", bounds=(1.0, 2.0, 3.0, 4.0), crs="EPSG:32632")
    assert kwarg.received_toa is raw["toa"] or kwarg.received_toa is not None
    assert bundle.bounds == (1.0, 2.0, 3.0, 4.0)
    assert bundle.crs == "EPSG:32632"
    assert bundle.sensor_config is mock_sensor_config

    @base_mod.register_preprocessor("Demo")
    class _RegisteredPreprocessor(_NoToaPreprocessor):
        pass

    assert "demo" in base_mod.list_available_sensors()
    assert isinstance(base_mod.get_preprocessor("DEMO"), _RegisteredPreprocessor)
    with pytest.raises(KeyError, match="Unknown sensor"):
        base_mod.get_preprocessor("unknown")

    l8_dir = tmp_path / "landsat8"
    l8_dir.mkdir()
    (l8_dir / "scene_MTL.txt").write_text("LANDSAT_8", encoding="utf-8")
    assert base_mod.detect_sensor(l8_dir) == "l8"

    l5_dir = tmp_path / "landsat5"
    l5_dir.mkdir()
    (l5_dir / "scene_MTL.xml").write_text("LANDSAT_5", encoding="utf-8")
    assert base_mod.detect_sensor(l5_dir) == "l5"

    aws_dir = tmp_path / "aws_s2"
    aws_dir.mkdir()
    (aws_dir / "metadata.xml").write_text("<xml/>", encoding="utf-8")
    assert base_mod.detect_sensor(aws_dir) == "s2"

    unknown_dir = tmp_path / "unknown"
    unknown_dir.mkdir()
    with pytest.raises(ValueError, match="Could not detect sensor"):
        base_mod.detect_sensor(unknown_dir)


def test_satellite_base_simple_math_helpers(mock_sensor_config) -> None:
    base_mod = _load_satellite_base_module()
    saa = xr.DataArray([10.0, 20.0], dims=("x",))
    vaa = xr.DataArray([40.0, 5.0], dims=("x",))
    raa = base_mod.compute_relative_azimuth(saa, vaa)
    assert list(raa.values) == [30.0, -15.0]

    scaled = base_mod.apply_scale_offset(xr.DataArray([1.0, 2.0]), scale=0.5, offset=1.0)
    assert list(scaled.values) == [1.5, 2.0]
    mask = base_mod.create_valid_mask(xr.DataArray([0.0, np.nan, 2.0]), min_val=0.0, max_val=1.5)
    assert list(mask.values) == [True, False, False]

    def _fake_reproject_match(angles, target, resampling="linear"):
        return xr.full_like(target, float(angles.mean()) + (0.0 if resampling == "linear" else 1.0))

    import siac.geo.reprojection as reproj_mod

    original = reproj_mod.reproject_match
    reproj_mod.reproject_match = _fake_reproject_match
    try:
        result = base_mod.resample_angles_to_data(
            xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=("y", "x")),
            xr.DataArray(np.zeros((3, 3), dtype=np.float32), dims=("y", "x")),
            method="nearest",
        )
    finally:
        reproj_mod.reproject_match = original

    assert float(result.mean()) == 2.0
    assert mock_sensor_config.sensor_id == "MOCK"
