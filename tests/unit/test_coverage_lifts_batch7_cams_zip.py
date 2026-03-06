"""Coverage lifts for atmospheric CAMS provider and ZIP store helpers."""

from __future__ import annotations

import asyncio
import sys
from datetime import datetime
from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

from siac.core.auth import CredentialManager
from siac.priors.atmospheric.cams import CAMSProvider
from siac.rt.lut import http_zip_store as zip_store

if TYPE_CHECKING:
    from pathlib import Path


class _Resp:
    def __init__(self, *, ok=True, status_code=200, headers=None, content=b""):
        self.ok = ok
        self.status_code = status_code
        self.headers = headers or {}
        self.content = content

    def raise_for_status(self):
        if not self.ok:
            raise RuntimeError("http error")


class _SessionProbeContentLength:
    def head(self, path, allow_redirects=True, headers=None, timeout=None):  # noqa: ARG002
        return _Resp(ok=False, headers={})

    def get(self, path, headers=None, timeout=None):  # noqa: ARG002
        if headers and "Range" in headers:
            return _Resp(ok=True, status_code=206, headers={"Content-Length": "12"}, content=b"x")
        return _Resp(ok=True, status_code=200, headers={"Content-Length": "12"}, content=b"x" * 12)

    def close(self):
        pass


class _SessionLastResortLength:
    def head(self, path, allow_redirects=True, headers=None, timeout=None):  # noqa: ARG002
        return _Resp(ok=False, headers={})

    def get(self, path, headers=None, timeout=None):  # noqa: ARG002
        if headers and "Range" in headers:
            return _Resp(ok=False, status_code=416, headers={}, content=b"")
        return _Resp(ok=True, status_code=200, headers={"Content-Length": "9"}, content=b"abcdefghi")

    def close(self):
        pass


class _SessionRangeFallbackShort:
    def __init__(self, payload: bytes):
        self.payload = payload

    def head(self, path, allow_redirects=True, headers=None, timeout=None):  # noqa: ARG002
        return _Resp(ok=True, status_code=200, headers={"Content-Length": "20"})

    def get(self, path, headers=None, timeout=None):  # noqa: ARG002
        if headers and "Range" in headers:
            return _Resp(ok=False, status_code=416, headers={}, content=b"")
        return _Resp(ok=True, status_code=200, headers={}, content=self.payload)

    def close(self):
        pass


class _BytesFS:
    def __init__(self, data: bytes):
        self.data = data

    async def _cat_file(self, path, start=None, end=None, **kwargs):  # noqa: ARG002
        s = len(self.data) + start if start is not None and start < 0 else (0 if start is None else start)
        e = len(self.data) if end is None else end
        return self.data[s:e]



def test_cams_source_and_load_branch_paths(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    p_missing = CAMSProvider(tmp_path / "missing")
    assert p_missing.source_name == "CAMS"
    assert p_missing._load_cams_data(datetime(2024, 1, 1)) is None

    p = CAMSProvider(tmp_path)
    tif_ds = xr.Dataset({"aod550": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])})
    monkeypatch.setattr(p, "_load_cams_tif_group", lambda _d, _i: tif_ds)
    loaded = p._load_cams_data(datetime(2024, 1, 1))
    assert loaded is tif_ds

    p_download = CAMSProvider(tmp_path, download_missing=True)
    downloaded = tmp_path / "CAMS_2024-01-01.nc"
    monkeypatch.setattr(p_download, "_load_cams_tif_group", lambda _d, _i: None)
    monkeypatch.setattr(p_download, "_download_cams_file", lambda _obs_time: downloaded)
    monkeypatch.setattr(p_download, "_load_from_explicit_path", lambda _path: "loaded")
    assert p_download._load_cams_data(datetime(2024, 1, 1)) == "loaded"



def test_cams_extract_and_tif_helpers(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    p = CAMSProvider(tmp_path, temporal_interp=True)

    ds = xr.Dataset(
        {
            "aod550": xr.DataArray(
                np.array([[[[0.2], [0.2], [0.2]]], [[[0.4], [0.4], [0.4]]]], dtype=np.float32),
                dims=["time", "band", "latitude", "longitude"],
                coords={
                    "time": [np.datetime64("2024-01-01T00:00:00"), np.datetime64("2024-01-01T02:00:00")],
                    "band": [1],
                    "latitude": [1.0, 0.0, -1.0],
                    "longitude": [0.0],
                },
            )
        }
    )
    out = p._extract_variable(ds, "aod550", (0.0, -1.0, 1.0, 1.0), "EPSG:4326", 1.0, datetime(2024, 1, 1, 1))
    assert "band" not in out.dims
    assert out.ndim == 2
    assert out.sizes["latitude"] > 0

    nc = tmp_path / "cams_bad.nc"
    nc.write_text("x")
    monkeypatch.setattr(xr, "open_dataset", lambda _path: (_ for _ in ()).throw(RuntimeError("bad nc")))
    assert p._load_from_explicit_path(nc) is None

    f1 = tmp_path / "cams_20240101_a.tif"
    f2 = tmp_path / "cams_20240101_b.tif"
    f1.write_text("x")
    f2.write_text("x")

    seq = [
        xr.Dataset({"aod550": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])}),
        None,
        xr.Dataset({"tcwv": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])}),
    ]

    def _fake_load_tif_dataset(path):  # noqa: ANN001
        return seq.pop(0)

    monkeypatch.setattr(p, "_load_tif_dataset", _fake_load_tif_dataset)
    merged = p._merge_tif_files([f1, f2, f1])
    assert merged is not None
    assert {"aod550", "tcwv"}.issubset(set(merged.data_vars))

    monkeypatch.setattr(
        p,
        "_merge_tif_files",
        lambda _files: xr.Dataset({"aod550": xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["y", "x"])}),
    )
    grouped = p._load_cams_tif_group("20240101", "2024-01-01")
    assert grouped is not None

    p2 = CAMSProvider(tmp_path)
    tif = tmp_path / "cams.tif"
    tif.write_text("x")
    def _raise_bad_tif(_path, engine=None):  # noqa: ANN001
        _ = engine
        raise RuntimeError("bad tif")

    monkeypatch.setattr(xr, "open_dataarray", _raise_bad_tif)
    assert p2._load_tif_dataset(tif) is None

    da_single = xr.DataArray(np.ones((1, 2, 2), dtype=np.float32), dims=["band", "y", "x"], coords={"band": [1]})
    def _return_single(_path, engine=None):  # noqa: ANN001
        _ = engine
        return da_single

    monkeypatch.setattr(xr, "open_dataarray", _return_single)
    monkeypatch.setattr(p2, "_infer_variable_name", lambda _name: None)
    assert p2._load_tif_dataset(tif) is None

    monkeypatch.setattr(p2, "_infer_variable_name", lambda _name: "tcwv")
    single = p2._load_tif_dataset(tif)
    assert single is not None and "tcwv" in single

    da_labels = xr.DataArray(
        np.ones((3, 1, 1), dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": [1, 2, 3]},
        attrs={"long_name": ["aod550", "aod550", "tcwv"]},
    )
    ds_labels = p2._dataset_from_multiband_tif(da_labels)
    assert set(ds_labels.data_vars) == {"aod550", "tcwv"}

    da_no_labels = xr.DataArray(
        np.ones((3, 1, 1), dtype=np.float32),
        dims=["band", "y", "x"],
        coords={"band": [1, 2, 3]},
    )
    ds_fallback = p2._dataset_from_multiband_tif(da_no_labels)
    assert set(ds_fallback.data_vars) == {"aod550", "tcwv", "gtco3"}

    assert p2._extract_band_labels(xr.DataArray(np.ones((1, 1, 1), dtype=np.float32), dims=["band", "y", "x"], attrs={"long_name": "aod550"})) == ["aod550"]
    assert p2._normalize_variable_name("unrelated-name") is None



def test_cams_download_auth_key_missing_branch(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    captured: dict[str, object] = {}

    class _FakeReq:
        def download(self, path):
            captured["download_path"] = path

    class _FakeClient:
        def __init__(self, **kwargs):
            captured["kwargs"] = kwargs

        def retrieve(self, dataset, request):
            captured["dataset"] = dataset
            return _FakeReq()

    monkeypatch.setitem(sys.modules, "cdsapi", SimpleNamespace(Client=_FakeClient))

    auth = CredentialManager()
    auth.set_credentials("cds", key="")
    p = CAMSProvider(tmp_path, auth=auth)
    out = p._download_cams_file(datetime(2024, 1, 2))
    assert out is not None
    assert captured["kwargs"] == {}



def test_http_zip_store_additional_branches(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    import requests

    monkeypatch.setattr(requests, "Session", lambda: _SessionProbeContentLength())
    fs_probe = zip_store._HTTPRangeFileSystem(timeout=2.0)
    assert fs_probe._discover_size("https://example.com/a.zip") == 12
    fs_probe.close()

    monkeypatch.setattr(requests, "Session", lambda: _SessionLastResortLength())
    fs_last = zip_store._HTTPRangeFileSystem(timeout=2.0)
    assert fs_last._discover_size("https://example.com/b.zip") == 9
    fs_last.close()

    monkeypatch.setattr(requests, "Session", lambda: _SessionRangeFallbackShort(b"short"))
    fs_short = zip_store._HTTPRangeFileSystem(timeout=2.0)
    fs_short._size_cache["https://example.com/c.zip"] = 20
    out = asyncio.run(fs_short._cat_file("https://example.com/c.zip", start=1, end=4))
    assert out == b"hor"
    fs_short.close()

    local = zip_store._build_local_filesystem()
    f = tmp_path / "x.bin"
    f.write_bytes(b"abcd")
    assert asyncio.run(local._cat_file(str(f), start=2, end=2)) == b""

    z_short = zip_store._ReadOnlyZipFileSystem(_BytesFS(b"abc"), "dummy.zip")
    with pytest.raises(ValueError, match="EOCD does not fit"):
        asyncio.run(z_short._initialize())

    z_noeocd = zip_store._ReadOnlyZipFileSystem(_BytesFS(b"x" * 40), "dummy.zip")
    with pytest.raises(ValueError, match="No EOCD found"):
        asyncio.run(z_noeocd._initialize())

    captured_opts: dict[str, object] = {}

    class _FakeHTTPFS:
        def __init__(self, **kwargs):
            captured_opts["http_kwargs"] = kwargs

    class _FakeZipFS:
        def __init__(self, fs, path, **kwargs):  # noqa: ANN001, ARG002
            self.fs = fs

    monkeypatch.setattr(zip_store, "_HTTPRangeFileSystem", _FakeHTTPFS)
    monkeypatch.setattr(zip_store, "_ReadOnlyZipFileSystem", _FakeZipFS)
    monkeypatch.setattr(zip_store, "_detect_zarr_prefix", lambda _zfs: "")
    def _fake_fsmap(root, fs, check=False, create=False):  # noqa: ANN001
        _ = (check, create)
        return {"root": root, "fs": fs}

    monkeypatch.setattr(zip_store, "FSMap", _fake_fsmap)
    _ = zip_store.build_readonly_zip_mapper(
        "https://example.com/lut.zip",
        {"headers": {"Authorization": "token"}},
    )
    assert captured_opts["http_kwargs"]["headers"] == {"Authorization": "token"}

    class _FakeS3:
        def __init__(self, **kwargs):
            captured_opts["s3_kwargs"] = kwargs

    monkeypatch.setitem(sys.modules, "s3fs", SimpleNamespace(S3FileSystem=_FakeS3))
    _ = zip_store._build_s3_filesystem({"anon": True, "key": "AK", "secret": "SK"})
    assert captured_opts["s3_kwargs"]["anon"] is True
    assert "key" not in captured_opts["s3_kwargs"]
