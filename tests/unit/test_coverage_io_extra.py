"""Additional coverage tests for readers/reprojection/writers modules."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from siac.catalog import SENTINEL2A_CONFIG
from siac.geo.reprojection import (
    align_grids,
    clip_to_bounds,
    clip_to_geometry,
    compute_common_bounds,
    get_crs,
    get_transform,
    reproject_dataset_match,
    reproject_match,
    reproject_to_crs,
    resample,
    resample_to_shape,
    transform_points,
)
from siac.runtime import CorrectionResult, GeometryAngles, ObservationBundle
from siac.storage.readers import (
    check_rasters_aligned,
    get_raster_info,
    read_hdf_subdataset,
    read_jp2,
    read_multiband,
    read_multiband_stack,
    read_netcdf_variable,
    read_raster,
    read_raster_at_resolution,
    read_raster_window,
    read_zarr_array,
)
from siac.storage.stac import build_stac_item, write_stac_item
from siac.storage.writers import (
    _compute_overview_levels,
    _prepare_for_write,
    write_auxiliary_products,
    write_boa_products,
    write_cog,
    write_dataset,
    write_netcdf,
    write_rgb_quicklook,
    write_zarr,
)


def _make_da(shape: tuple[int, int] = (32, 32), crs: str = "EPSG:32632") -> xr.DataArray:
    """Create a small spatial DataArray with CRS."""
    data = np.arange(shape[0] * shape[1], dtype=np.float32).reshape(shape)
    x = np.linspace(500000.0, 500000.0 + (shape[1] - 1) * 10.0, shape[1])
    y = np.linspace(4500000.0 + (shape[0] - 1) * 10.0, 4500000.0, shape[0])
    da = xr.DataArray(data, dims=["y", "x"], coords={"y": y, "x": x})
    return da.rio.write_crs(crs)


def _write_multiband(path: Path) -> None:
    """Write a 2-band test raster."""
    band1 = _make_da((16, 16))
    band2 = band1 + 1000.0
    stacked = xr.concat([band1, band2], dim="band").assign_coords(band=[1, 2])
    stacked.rio.to_raster(path)


class TestReadersExtra:
    def test_read_raster_band_select_and_squeeze(self, tmp_path: Path):
        p = tmp_path / "multi.tif"
        _write_multiband(p)
        da = read_raster(p, band=1)
        assert da.dims == ("y", "x")
        assert da.shape == (16, 16)

    def test_read_raster_window_and_resolution_paths(self, tmp_path: Path):
        p = tmp_path / "base.tif"
        _make_da((32, 32)).rio.to_raster(p)

        clipped = read_raster_window(p, bounds=(500050.0, 4500050.0, 500150.0, 4500150.0))
        assert clipped.size > 0

        # Already close to target resolution -> early return path
        same = read_raster_at_resolution(p, target_resolution=10.0)
        assert same.shape == (32, 32)

        # Coarsen path
        coarse = read_raster_at_resolution(p, target_resolution=20.0)
        assert coarse.shape[0] < 32
        assert coarse.shape[1] < 32

    def test_read_multiband_and_stack_and_mismatch(self, tmp_path: Path):
        p1 = tmp_path / "b1.tif"
        p2 = tmp_path / "b2.tif"
        _make_da((12, 12)).rio.to_raster(p1)
        (_make_da((12, 12)) + 1).rio.to_raster(p2)

        ds = read_multiband([p1, p2], band_names=["B1", "B2"])
        assert set(ds.data_vars) == {"B1", "B2"}

        stacked = read_multiband_stack([p1, p2], band_names=["B1", "B2"])
        assert stacked.shape == (2, 12, 12)

        with pytest.raises(ValueError):
            read_multiband([p1, p2], band_names=["B1"])

    def test_specialized_reader_wrappers(self, tmp_path: Path):
        p = tmp_path / "a.tif"
        _make_da((8, 8)).rio.to_raster(p)
        assert read_jp2(p).shape == (8, 8)

        # HDF wrapper just builds a GDAL URI and forwards to read_raster.
        # We validate wiring by monkeypatching read_raster.
        import siac.storage.readers as m

        seen = {}
        orig = m.read_raster
        try:
            def _fake(path, **kwargs):
                seen["path"] = path
                return _make_da((4, 4))
            m.read_raster = _fake  # type: ignore[assignment]
            out = read_hdf_subdataset("f.hdf", "SUB")
            assert out.shape == (4, 4)
            assert str(seen["path"]).startswith('HDF4_EOS:EOS_GRID:"f.hdf"')
        finally:
            m.read_raster = orig  # type: ignore[assignment]

    def test_netcdf_and_zarr_readers(self, tmp_path: Path):
        da = _make_da((10, 10))
        ds = da.to_dataset(name="var")
        ds.attrs["crs"] = "EPSG:32632"
        nc = tmp_path / "x.nc"
        ds.to_netcdf(nc)
        da_nc = read_netcdf_variable(nc, "var")
        assert da_nc.shape == (10, 10)

        zarr_path = tmp_path / "x.zarr"
        ds.to_zarr(zarr_path, mode="w")
        out_ds = read_zarr_array(zarr_path)
        out_da = read_zarr_array(zarr_path, variable="var")
        assert isinstance(out_ds, xr.Dataset)
        assert isinstance(out_da, xr.DataArray)

    def test_raster_info_and_alignment(self, tmp_path: Path):
        p1 = tmp_path / "r1.tif"
        p2 = tmp_path / "r2.tif"
        p3 = tmp_path / "r3.tif"
        _make_da((10, 10), crs="EPSG:32632").rio.to_raster(p1)
        _make_da((10, 10), crs="EPSG:32632").rio.to_raster(p2)
        _make_da((10, 10), crs="EPSG:4326").rio.to_raster(p3)

        info = get_raster_info(p1)
        assert info["driver"] in {"GTiff", "COG"}
        assert info["shape"][0] == 1
        assert check_rasters_aligned(p1, p2)
        assert not check_rasters_aligned(p1, p3)


class TestReprojectionExtra:
    def test_reproject_and_resample_variants(self):
        da = _make_da((24, 24), crs="EPSG:32632")
        out = reproject_to_crs(da, "EPSG:4326", resolution=0.0002, resampling="nearest")
        assert out.rio.crs is not None

        target = _make_da((12, 12), crs="EPSG:32632")
        matched = reproject_match(da, target, resampling="bilinear")
        assert matched.shape == target.shape

        ds = xr.Dataset({"a": da})
        out_ds = reproject_dataset_match(ds, target)
        assert out_ds["a"].shape == target.shape

        r1 = resample(da, target_resolution=20.0)
        r2 = resample_to_shape(da, target_shape=(8, 8))
        assert r1.shape[0] < da.shape[0]
        assert r2.shape == (8, 8)

    def test_clip_transform_and_grid_utils(self):
        da = _make_da((20, 20), crs="EPSG:32632")
        b = da.rio.bounds()
        clipped = clip_to_bounds(da, (b[0], b[1], b[0] + 50.0, b[1] + 50.0))
        assert clipped.size > 0

        geom = {
            "type": "Polygon",
            "coordinates": [[
                [b[0], b[1]],
                [b[0] + 70.0, b[1]],
                [b[0] + 70.0, b[1] + 70.0],
                [b[0], b[1] + 70.0],
                [b[0], b[1]],
            ]],
        }
        clipped_geom = clip_to_geometry(da, geom)
        assert clipped_geom.size > 0

        x2, y2 = transform_points([0.0, 1.0], [0.0, 1.0], "EPSG:4326", "EPSG:3857")
        assert x2.shape == (2,)
        assert y2.shape == (2,)

        assert get_transform(da) is not None
        assert get_crs(da) == "EPSG:32632"

        same = align_grids(da)
        assert len(same) == 1

        other = _make_da((18, 18), crs="EPSG:32632")
        aligned = align_grids(da, other, reference_idx=0)
        assert len(aligned) == 2
        assert aligned[1].shape == da.shape

        u = compute_common_bounds(da, other, method="union")
        i = compute_common_bounds(da, other, method="intersection")
        assert u[0] <= i[0]


class TestWritersExtra:
    def test_prepare_helpers(self):
        da = _make_da((16, 16))
        out = _prepare_for_write(da, dtype="uint16", nodata=None)
        assert out.dtype == np.uint16

        out2 = _prepare_for_write(da, dtype=None, nodata=None)
        assert out2.rio.nodata is not None

        levels = _compute_overview_levels(_make_da((2048, 2048)))
        assert levels and levels[0] == 2

    def test_build_and_write_stac_item(self, tmp_path: Path):
        template = _make_da((16, 16))
        toa = xr.Dataset({"B02": template.copy(), "B08": (template + 0.1).rename("B08")})
        geometry = GeometryAngles(
            sza=xr.full_like(template, np.deg2rad(30.0)),
            saa=xr.full_like(template, np.deg2rad(150.0)),
            vza=xr.full_like(template, np.deg2rad(5.0)),
            vaa=xr.full_like(template, np.deg2rad(60.0)),
        )
        cloud_mask = xr.DataArray(
            np.zeros((16, 16), dtype=bool),
            dims=["y", "x"],
            coords={"y": template.coords["y"], "x": template.coords["x"]},
        )
        cloud_mask.values[:4, :8] = True
        obs = ObservationBundle(
            toa=toa,
            geometry=geometry,
            cloud_mask=cloud_mask,
            sensor_config=SENTINEL2A_CONFIG,
            metadata={
                "observation_time": datetime(2026, 1, 2, 2, 41, 21),
                "satellite": "S2A",
                "sensor": "MSI",
                "tile_id": "T50QLD",
                "processing_baseline": "N0511",
            },
            crs="EPSG:32632",
            bounds=tuple(map(float, template.rio.bounds())),
        )
        result = CorrectionResult(
            boa=xr.Dataset({"B02": template.copy(), "B08": (template + 0.1).rename("B08")}),
            boa_unc=None,
            aot=xr.full_like(template, 0.12),
            tcwv=xr.full_like(template, 1.8),
            cloud_mask=xr.full_like(template, True, dtype=bool),
            metadata={"processing_time_s": 12.5},
        )

        boa_dir = tmp_path / "boa"
        boa_dir.mkdir()
        boa_assets = {"B02": boa_dir / "B02.tif", "B08": boa_dir / "B08.tif"}
        for path in boa_assets.values():
            path.write_bytes(b"boa")
        atmosphere_path = tmp_path / "atmosphere.nc"
        atmosphere_path.write_bytes(b"netcdf")
        qa_dir = tmp_path / "qa"
        qa_dir.mkdir()
        qa_assets = {"cloud_mask": qa_dir / "cloud_mask.tif"}
        qa_assets["cloud_mask"].write_bytes(b"mask")
        summary_path = tmp_path / "run_summary.json"
        summary_path.write_text("{}", encoding="utf-8")

        item = build_stac_item(
            obs,
            result,
            output_dir=tmp_path,
            boa_assets=boa_assets,
            atmosphere_asset=atmosphere_path,
            qa_assets=qa_assets,
            summary_asset=summary_path,
            input_href=tmp_path / "S2C_MSIL1C_20260102T024121_N0511_R089_T50QLD_20260102T035433.SAFE",
        )

        assert item["stac_version"] == "1.0.0"
        assert item["id"] == tmp_path.name
        assert item["properties"]["platform"] == "sentinel-2c"
        assert item["properties"]["constellation"] == "sentinel-2"
        assert item["properties"]["instruments"] == ["msi"]
        assert item["properties"]["siac:satellite"] == "S2C"
        assert item["properties"]["eo:cloud_cover"] == pytest.approx(12.5)
        assert item["properties"]["view:sun_elevation"] == pytest.approx(60.0)
        assert item["properties"]["view:sun_azimuth"] == pytest.approx(150.0)
        assert item["properties"]["view:off_nadir"] == pytest.approx(5.0)
        assert item["properties"]["proj:epsg"] == 32632
        assert item["assets"]["B02"]["href"] == "boa/B02.tif"
        assert item["assets"]["B02"]["eo:bands"][0]["common_name"] == "blue"
        assert item["assets"]["atmosphere"]["href"] == "atmosphere.nc"
        assert item["assets"]["cloud-mask"]["href"] == "qa/cloud_mask.tif"
        assert item["assets"]["summary"]["href"] == "run_summary.json"
        assert item["links"][0]["href"] == "item.json"
        assert item["links"][1]["rel"] == "derived_from"

        item_path = write_stac_item(
            obs,
            result,
            output_dir=tmp_path,
            boa_assets=boa_assets,
            atmosphere_asset=atmosphere_path,
            qa_assets=qa_assets,
            summary_asset=summary_path,
            input_href=tmp_path / "S2C_MSIL1C_20260102T024121_N0511_R089_T50QLD_20260102T035433.SAFE",
        )
        assert item_path.exists()
        written = json.loads(item_path.read_text(encoding="utf-8"))
        assert written["assets"]["B08"]["eo:bands"][0]["name"] == "B08"

    def test_write_netcdf_and_siac_product_writers(self, tmp_path: Path):
        da = _make_da((16, 16))
        ds = xr.Dataset(
            {
                "B02": da / 10000.0,
                "B03": da / 11000.0,
                "B04": da / 12000.0,
                "B04_unc": xr.full_like(da, 0.01),
            }
        )

        nc = tmp_path / "a.nc"
        write_netcdf(ds, nc, compression={})
        assert nc.exists()

        out_dir = tmp_path / "boa"
        paths = write_boa_products(
            ds,
            output_dir=out_dir,
            sensor="S2A",
            tile_id="T31UDQ",
            datetime_str="20240101T100000",
            include_uncertainty=False,
        )
        assert "B04_unc" not in paths
        assert len(paths) == 3

        aux = write_auxiliary_products(
            da / 1000.0,
            da / 2000.0,
            output_dir=tmp_path / "aux",
            sensor="S2A",
            tile_id="T31UDQ",
            datetime_str="20240101T100000",
        )
        assert aux["aot"].exists()
        assert aux["tcwv"].exists()

    def test_write_rgb_quicklook(self, tmp_path: Path):
        da = _make_da((64, 64))
        boa = xr.Dataset(
            {"B02": da / 10000.0, "B03": da / 10000.0, "B04": da / 10000.0}
        )
        out = tmp_path / "quicklook.tif"
        p = write_rgb_quicklook(boa, out, target_resolution=40.0)
        assert p.exists()

    def test_write_dataset_skip_nonspatial_and_write_zarr_chunks(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ):
        spatial = _make_da((16, 16))
        ds_skip = xr.Dataset(
            {
                "meta": xr.DataArray(np.ones((16, 16), dtype=np.int16), dims=["y", "x"]),
            }
        )
        paths = write_dataset(ds_skip, tmp_path / "out_skip", as_cog=False)
        assert paths == {}

        ds = xr.Dataset({"B02": spatial})
        paths2 = write_dataset(ds, tmp_path / "out_ok", as_cog=False)
        assert "B02" in paths2

        monkeypatch.setattr(xr.Dataset, "chunk", lambda self, _chunks: self)
        zarr_path = tmp_path / "chunked.zarr"
        write_zarr(ds[["B02"]], zarr_path, chunks={"y": 8, "x": 8})
        assert zarr_path.exists()

    def test_write_netcdf_dataarray_default_compression_and_write_cog_options(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ):
        da = _make_da((16, 16))
        captures: list[dict[str, object]] = []

        monkeypatch.setattr(
            xr.DataArray,
            "to_netcdf",
            lambda _self, path, encoding=None, **kwargs: captures.append(  # noqa: ANN001
                {"path": path, "encoding": encoding, "kwargs": kwargs}
            ),
        )
        write_netcdf(da.rename(None), tmp_path / "d.nc", compression=None)
        assert captures and "data" in captures[-1]["encoding"]
        assert captures[-1]["kwargs"]["engine"] == "h5netcdf"
        assert captures[-1]["encoding"]["x"]["_FillValue"] is None
        assert captures[-1]["encoding"]["y"]["_FillValue"] is None
        assert captures[-1]["encoding"]["spatial_ref"]["_FillValue"] is None

        cog_calls: list[dict[str, object]] = []

        def _fake_to_raster(self, path, **kwargs):  # noqa: ANN001
            cog_calls.append(kwargs)
            Path(path).write_bytes(b"x")

        monkeypatch.setattr(type(da.rio), "to_raster", _fake_to_raster, raising=False)
        out = write_cog(da, tmp_path / "x.tif", compression="deflate")
        assert out.exists()
        assert cog_calls
        assert cog_calls[-1]["driver"] == "COG"
        assert cog_calls[-1]["compress"] == "deflate"
        assert cog_calls[-1]["level"] == 6
        assert "zlevel" not in cog_calls[-1]

    def test_write_netcdf_preserves_grid_mapping_encoding(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ):
        spatial = _make_da((8, 8))
        ds = xr.Dataset({"aot": spatial.rename("aot")}).rio.write_crs("EPSG:32632")
        ds["aot"].encoding.clear()
        captures: list[dict[str, object]] = []

        monkeypatch.setattr(
            xr.Dataset,
            "to_netcdf",
            lambda _self, path, encoding=None, **kwargs: captures.append(  # noqa: ANN001
                {"path": path, "encoding": encoding, "kwargs": kwargs}
            ),
        )

        write_netcdf(ds, tmp_path / "grid.nc", compression=None)
        assert captures
        assert captures[-1]["encoding"]["aot"]["grid_mapping"] == "spatial_ref"

    def test_write_netcdf_rejects_scipy_and_missing_supported_engines(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ):
        da = _make_da((8, 8)).rename("data")

        with pytest.raises(ValueError, match="not supported"):
            write_netcdf(da, tmp_path / "scipy.nc", engine="scipy")

        def _missing(_name: str):  # noqa: ANN001
            return None

        monkeypatch.setattr("siac.storage.writers.importlib.util.find_spec", _missing)
        with pytest.raises(RuntimeError, match="requires h5netcdf or netCDF4"):
            write_netcdf(da, tmp_path / "missing.nc")

    def test_write_netcdf_sets_no_fillvalue_for_integer_data(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ):
        cloud_mask = _make_da((8, 8)).astype(np.uint8)
        cloud_mask.attrs["_FillValue"] = np.nan
        ds = xr.Dataset({"cloud_mask": cloud_mask})
        captures: list[dict[str, object]] = []

        monkeypatch.setattr(
            xr.Dataset,
            "to_netcdf",
            lambda _self, path, encoding=None, **kwargs: captures.append(  # noqa: ANN001
                {"path": path, "encoding": encoding, "kwargs": kwargs, "attrs": dict(_self["cloud_mask"].attrs)}
            ),
        )

        write_netcdf(ds, tmp_path / "mask.nc", compression=None)

        assert captures
        assert captures[-1]["kwargs"]["engine"] == "h5netcdf"
        assert captures[-1]["encoding"]["cloud_mask"]["_FillValue"] is None
        assert "_FillValue" not in captures[-1]["attrs"]
