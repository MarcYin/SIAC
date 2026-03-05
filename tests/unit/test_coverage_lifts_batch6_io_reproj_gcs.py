"""Coverage lifts for io.readers, io.reprojection, and io.gcs_sentinel2."""

from __future__ import annotations

import sys
from datetime import date, datetime
from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr
from rasterio.enums import Resampling

import siac.io.gcs_sentinel2 as gcs_mod
from siac.core.exceptions import DataNotFoundError
from siac.io import readers as readers_mod
from siac.io import reprojection as reproj_mod
from siac.io.gcs_sentinel2 import GCS_DOWNLOAD_BASE
from siac.io.s2_data_source import S2Product, S2Query

if TYPE_CHECKING:
    from pathlib import Path


def _spatial_da(shape: tuple[int, int] = (8, 8), crs: str = "EPSG:32632") -> xr.DataArray:
    data = np.arange(shape[0] * shape[1], dtype=np.float32).reshape(shape)
    x = np.linspace(500000.0, 500000.0 + (shape[1] - 1) * 10.0, shape[1])
    y = np.linspace(4500000.0 + (shape[0] - 1) * 10.0, 4500000.0, shape[0])
    return xr.DataArray(data, dims=["y", "x"], coords={"x": x, "y": y}).rio.write_crs(crs)


def _product(product_id: str, source_url: str = "") -> S2Product:
    return S2Product(
        product_id=product_id,
        mgrs_tile="31UDQ",
        sensing_date=datetime(2024, 1, 3, 10, 30, 21),
        processing_baseline="N0500",
        cloud_cover=100.0,
        satellite="S2A",
        orbit_number=51,
        source_url=source_url,
        size_mb=None,
    )


def test_readers_open_kwargs_defaults_remote_and_alignment_branches(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    tmp_path: Path,
) -> None:
    captured: dict[str, object] = {}

    def _fake_open_dataarray(path, engine=None, **kwargs):  # noqa: ANN001
        captured["kwargs"] = kwargs
        return xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])

    monkeypatch.setattr(xr, "open_dataarray", _fake_open_dataarray)

    with caplog.at_level("WARNING"):
        da = readers_mod.read_raster("dummy.tif", band=1, chunks={"x": 1}, overview_level=2)
    assert da.shape == (2, 2)
    assert "No CRS found" in caplog.text
    assert captured["kwargs"]["chunks"] == {"x": 1}
    assert captured["kwargs"]["overview_level"] == 2

    # Default band-name paths + squeeze in multiband readers.
    band1 = xr.DataArray(np.ones((1, 3, 3), dtype=np.float32), dims=["band", "y", "x"], coords={"band": [1]})
    monkeypatch.setattr(readers_mod, "read_raster", lambda *_args, **_kwargs: band1)

    ds = readers_mod.read_multiband([tmp_path / "a.tif"])  # band_names=None branch
    assert set(ds.data_vars) == {"a"}

    stacked = readers_mod.read_multiband_stack([tmp_path / "a.tif", tmp_path / "b.tif"])  # default names
    assert stacked.shape == (2, 3, 3)
    assert list(stacked.coords["band"].values) == ["a", "b"]

    # NetCDF spatial_ref branch.
    nc = tmp_path / "with_spatial_ref.nc"
    xr.Dataset(
        {
            "var": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"]),
            "spatial_ref": xr.DataArray(0, attrs={"crs_wkt": "EPSG:4326"}),
        }
    ).to_netcdf(nc)
    out = readers_mod.read_netcdf_variable(nc, "var")
    assert out.shape == (2, 2)

    # Remote zarr mapper path.
    monkeypatch.setitem(sys.modules, "fsspec", SimpleNamespace(get_mapper=lambda path: {"remote": path}))
    opened: dict[str, object] = {}

    def _fake_open_zarr(obj, chunks=None):  # noqa: ANN001
        opened["obj"] = obj
        return xr.Dataset({"v": xr.DataArray([1], dims=["x"])})

    monkeypatch.setattr(xr, "open_zarr", _fake_open_zarr)
    readers_mod.read_zarr_array("https://example.com/data.zarr")
    assert opened["obj"] == {"remote": "https://example.com/data.zarr"}

    # Alignment false-path branches: resolution mismatch then bounds mismatch.
    def _info_res(path):  # noqa: ANN001
        return {
            "crs": "EPSG:32632",
            "resolution": (10.0, 10.0) if str(path).endswith("a") else (20.0, 10.0),
            "bounds": SimpleNamespace(left=0.0, bottom=0.0, right=1.0, top=1.0),
        }

    monkeypatch.setattr(readers_mod, "get_raster_info", _info_res)
    assert not readers_mod.check_rasters_aligned("a", "b")

    def _info_bounds(path):  # noqa: ANN001
        if str(path).endswith("a"):
            b = SimpleNamespace(left=0.0, bottom=0.0, right=1.0, top=1.0)
        else:
            b = SimpleNamespace(left=0.1, bottom=0.0, right=1.1, top=1.0)
        return {"crs": "EPSG:32632", "resolution": (10.0, 10.0), "bounds": b}

    monkeypatch.setattr(readers_mod, "get_raster_info", _info_bounds)
    assert not readers_mod.check_rasters_aligned("a", "c")


def test_reprojection_missing_branches(monkeypatch: pytest.MonkeyPatch) -> None:
    da = _spatial_da((12, 12), crs="EPSG:32632")

    out = reproj_mod.reproject_to_crs(
        da,
        "EPSG:4326",
        resolution=(0.0002, 0.0002),
        resampling="nearest",
        nodata=-9999.0,
    )
    assert out.rio.crs is not None

    target = _spatial_da((6, 6), crs="EPSG:32632")
    matched = reproj_mod.reproject_match(da, target, resampling=Resampling.nearest, nodata=-9999.0)
    assert matched.shape == target.shape

    src_ds = xr.Dataset({"meta": xr.DataArray(np.array([1, 2, 3], dtype=np.int16), dims=["k"])})
    target_no_crs = xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])
    ds_out = reproj_mod.reproject_dataset_match(src_ds, target_no_crs, resampling=Resampling.nearest)
    assert "meta" in ds_out

    r_tuple = reproj_mod.resample(da, target_resolution=(20.0, 20.0), resampling=Resampling.nearest)
    assert r_tuple.shape[0] < da.shape[0]

    r_shape = reproj_mod.resample_to_shape(da, target_shape=(4, 4), resampling=Resampling.nearest)
    assert r_shape.shape == (4, 4)

    assert reproj_mod.get_crs(xr.DataArray(np.ones((2, 2)), dims=["y", "x"])) is None

    with pytest.raises(ValueError, match="At least one array required"):
        reproj_mod.compute_common_bounds()

    calls: list[tuple[tuple[float, float, float, float], object, object]] = []

    def _fake_transform_bounds(bounds, src, dst):  # noqa: ANN001
        calls.append((bounds, src, dst))
        return bounds

    monkeypatch.setattr(reproj_mod, "transform_bounds", _fake_transform_bounds)
    _ = reproj_mod.compute_common_bounds(_spatial_da((4, 4), "EPSG:32632"), _spatial_da((4, 4), "EPSG:4326"), method="union")
    assert calls

    a = _spatial_da((4, 4), "EPSG:32632")
    b = _spatial_da((4, 4), "EPSG:32632").assign_coords(x=np.linspace(600000.0, 600030.0, 4))
    with pytest.raises(ValueError, match="do not overlap"):
        reproj_mod.compute_common_bounds(a, b, method="intersection")


def test_gcs_helpers_and_error_paths(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Invalid MGRS tile"):
        gcs_mod._normalize_mgrs_tile("bad")

    with pytest.raises(ValueError, match="Cannot infer MGRS tile"):
        gcs_mod._safe_prefix_from_product_id("S2A_INVALID")

    with pytest.raises(ValueError, match="Cannot parse sensing datetime"):
        gcs_mod._parse_sensing_datetime("S2A_INVALID")

    assert gcs_mod._parse_baseline("S2A_X") == "N0000"
    assert gcs_mod._parse_orbit("S2A_X") == 0
    assert gcs_mod._parse_satellite("X_S2") == "S2"

    with pytest.raises(ValueError, match="Invalid SAFE prefix"):
        gcs_mod._product_from_safe_prefix("tiles/31/U/DQ/not_safe")

    with pytest.raises(ValueError, match="Cannot parse MGRS tile"):
        gcs_mod._product_from_safe_prefix(
            "tiles/31/U/DQ/S2A_MSIL1C_20240103T103021_N0500_R051_BAD_20240103T120000.SAFE/"
        )

    with pytest.raises(ValueError, match="Unsupported GCS source URL"):
        gcs_mod._safe_prefix_from_source_url("https://example.com/abc")

    with pytest.raises(ValueError, match="Cannot infer SAFE prefix"):
        gcs_mod._safe_prefix_from_source_url(f"gs://{gcs_mod.GCS_BUCKET}/tiles/31/U/DQ/no_safe_here")

    encoded = GCS_DOWNLOAD_BASE + "tiles/31/U/DQ/P.SAFE/GRANULE/L1C.xml"
    assert gcs_mod._safe_prefix_from_source_url(encoded).endswith("P.SAFE/")

    prod = _product("S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000", source_url="gs://invalid/path")
    fallback = gcs_mod._resolve_safe_prefix(prod)
    assert fallback.endswith(".SAFE/")

    safe_dir = tmp_path / "safe"
    safe_dir.mkdir()
    safe_prefix = "tiles/31/U/DQ/P.SAFE/"
    assert gcs_mod._target_path_for_object("other/p.txt", safe_prefix=safe_prefix, safe_dir=safe_dir) is None
    assert gcs_mod._target_path_for_object(safe_prefix, safe_prefix=safe_prefix, safe_dir=safe_dir) is None
    with pytest.raises(DataNotFoundError, match="outside SAFE root"):
        gcs_mod._target_path_for_object(f"{safe_prefix}../escape.txt", safe_prefix=safe_prefix, safe_dir=safe_dir)

    done = safe_dir / "done.bin"
    done.write_bytes(b"abc")
    assert gcs_mod._is_fully_downloaded(done, None)

    monkeypatch.setattr(gcs_mod, "_list_api", lambda **_kwargs: {"items": [], "nextPageToken": None})
    assert not gcs_mod._prefix_exists("tiles/31/U/DQ/P.SAFE/")

    class _Resp:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):  # noqa: ANN001
            return False

        def read(self, _size):  # noqa: ANN001
            raise OSError("boom")

    import urllib.request

    monkeypatch.setattr(urllib.request, "urlopen", lambda _req, **_kwargs: _Resp())
    target = tmp_path / "x.bin"
    with pytest.raises(OSError):
        gcs_mod._download_url_to_file("https://example.com/x", target)
    assert not target.with_suffix(".bin.part").exists()


def test_gcs_retry_search_and_download_edge_cases(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    target = tmp_path / "retry.bin"

    def _short_write(url, path, timeout=300):  # noqa: ANN001
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"x")

    monkeypatch.setattr(gcs_mod, "_download_url_to_file", _short_write)
    with pytest.raises(RuntimeError, match="Failed downloading"):
        gcs_mod._download_with_retry("https://example.com/x", target, expected_size=2, retries=0, backoff_sec=0.0)

    # no-jobs branch
    gcs_mod._download_jobs_parallel([])

    q_bad = S2Query(product_id="S2A_INVALID")
    assert gcs_mod.search_gcs(q_bad) == []

    q_good = S2Query(product_id="S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000")
    monkeypatch.setattr(gcs_mod, "_prefix_exists", lambda _prefix: False)
    assert gcs_mod.search_gcs(q_good) == []

    q_tile = S2Query(mgrs_tile="31UDQ", date=date(2024, 1, 3))
    monkeypatch.setattr(
        gcs_mod,
        "_list_safe_prefixes",
        lambda _prefix: [
            "bad_prefix",
            "tiles/31/U/DQ/S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000.SAFE/",
        ],
    )
    out = gcs_mod.search_gcs(q_tile)
    assert len(out) == 1

    product_id = "S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000"
    safe_prefix = f"tiles/31/U/DQ/{product_id}.SAFE/"
    product = _product(product_id, source_url=f"gs://{gcs_mod.GCS_BUCKET}/{safe_prefix}")

    monkeypatch.setattr(gcs_mod, "_list_objects_under", lambda _prefix: [{"name": 123, "size": "1"}])
    with pytest.raises(DataNotFoundError, match="SAFE directory is empty"):
        gcs_mod.download_gcs(product, tmp_path)
