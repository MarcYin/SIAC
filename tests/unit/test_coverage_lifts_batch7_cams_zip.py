"""Coverage lifts for atmospheric CAMS provider and ZIP store helpers."""

from __future__ import annotations

import asyncio
import sys
from contextlib import contextmanager
from datetime import datetime
from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr
from pyproj import Transformer

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
    monkeypatch.setattr(p_missing, "_complete_cams_dataset", lambda dataset, _obs_time: dataset)
    assert p_missing.source_name == "CAMS"
    assert p_missing._load_cams_data(datetime(2024, 1, 1)) is None

    p = CAMSProvider(tmp_path)
    tif_ds = xr.Dataset({"aod550": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])})
    monkeypatch.setattr(p, "_complete_cams_dataset", lambda dataset, _obs_time: dataset)
    monkeypatch.setattr(p, "_load_cams_tif_group", lambda _d, _i: tif_ds)
    loaded = p._load_cams_data(datetime(2024, 1, 1))
    assert loaded is tif_ds

    p_download = CAMSProvider(tmp_path, download_missing=True)
    downloaded = tmp_path / "CAMS_2024-01-01.nc"
    monkeypatch.setattr(
        p_download,
        "_complete_cams_dataset",
        lambda dataset, _obs_time: "loaded" if dataset is None else dataset,
    )
    monkeypatch.setattr(p_download, "_load_cams_tif_group", lambda _d, _i: None)
    monkeypatch.setattr(p_download, "_download_cams_file", lambda _obs_time: downloaded)
    monkeypatch.setattr(p_download, "_load_from_explicit_path", lambda _path: "loaded")
    assert p_download._load_cams_data(datetime(2024, 1, 1)) == "loaded"


def test_cams_remote_base_and_remote_file_paths(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cached = tmp_path / "cached"
    cached.write_text("x")
    seen: list[str] = []

    remote_base = CAMSProvider(
        "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/",
        cache_dir=tmp_path / "cache",
    )
    monkeypatch.setattr(remote_base, "_complete_cams_dataset", lambda dataset, _obs_time: dataset)

    def _cache_for_base(url: str, storage_options=None) -> Path:  # noqa: ARG001
        seen.append(url)
        if url.endswith("2024-01-01.nc"):
            return cached
        raise FileNotFoundError(url)

    monkeypatch.setattr(remote_base, "_cache_remote_file", _cache_for_base)
    monkeypatch.setattr(remote_base, "_load_from_local_explicit_path", lambda _path, source_name=None: source_name)
    loaded = remote_base._load_cams_data(datetime(2024, 1, 1))
    assert loaded == "2024-01-01.nc"
    assert seen == ["https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/2024-01-01.nc"]

    remote_file = CAMSProvider(
        "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/2024-01-02.nc",
        cache_dir=tmp_path / "cache",
    )
    monkeypatch.setattr(remote_file, "_complete_cams_dataset", lambda dataset, _obs_time: dataset)
    def _cache_remote_file(_url, storage_options=None):  # noqa: ARG001
        return cached

    monkeypatch.setattr(remote_file, "_cache_remote_file", _cache_remote_file)
    monkeypatch.setattr(remote_file, "_load_from_local_explicit_path", lambda _path, source_name=None: source_name)
    assert remote_file._load_cams_data(datetime(2024, 1, 2)) == "2024-01-02.nc"
    assert remote_file._is_remote_source("s3://eodata/CAMS/GLOBAL") is True


def test_cams_select_cdse_s3_files(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    import fsspec

    provider = CAMSProvider("s3://eodata/CAMS/GLOBAL", cache_dir=tmp_path / "cache")

    class _FakeFS:
        def ls(self, path, detail=False):  # noqa: ARG002
            assert path == "eodata/CAMS/GLOBAL/2024/01/01"
            return [
                "eodata/CAMS/GLOBAL/2024/01/01/z_cams_c_ecmf_20240101000000_prod_an_sfc_000_aod550",
                "eodata/CAMS/GLOBAL/2024/01/01/z_cams_c_ecmf_20240101000000_prod_fc_sfc_003_aod550",
                "eodata/CAMS/GLOBAL/2024/01/01/z_cams_c_ecmf_20240101000000_prod_an_sfc_000_gtco3",
            ]

    def _fake_filesystem(_protocol, **_kwargs):
        return _FakeFS()

    monkeypatch.setattr(fsspec, "filesystem", _fake_filesystem)

    selected = provider._select_cdse_cams_files(
        "s3://eodata/CAMS/GLOBAL",
        datetime(2024, 1, 1, 2, 30),
        {},
    )

    assert selected == [
        "s3://eodata/CAMS/GLOBAL/2024/01/01/z_cams_c_ecmf_20240101000000_prod_fc_sfc_003_aod550/z_cams_c_ecmf_20240101000000_prod_fc_sfc_003_aod550.nc",
        "s3://eodata/CAMS/GLOBAL/2024/01/01/z_cams_c_ecmf_20240101000000_prod_an_sfc_000_gtco3/z_cams_c_ecmf_20240101000000_prod_an_sfc_000_gtco3.nc",
    ]


def test_cams_remote_s3_base_merges_selected_datasets(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    provider = CAMSProvider("s3://eodata/CAMS/GLOBAL", cache_dir=tmp_path / "cache")

    @contextmanager
    def _fake_storage_context(_url):
        yield {"key": "k", "secret": "s"}

    monkeypatch.setattr(provider, "_remote_storage_options_context", _fake_storage_context)
    monkeypatch.setattr(
        provider,
        "_select_cdse_cams_files",
        lambda _base, _time, _opts: ["s3://eodata/aod550.nc", "s3://eodata/gtco3.nc"],
    )

    def _load(url, *, missing_ok=False, storage_options=None):  # noqa: ARG001
        if url.endswith("aod550.nc"):
            return xr.Dataset({"aod550": xr.DataArray(np.ones((1, 2, 2), dtype=np.float32), dims=["time", "latitude", "longitude"])})
        if url.endswith("gtco3.nc"):
            return xr.Dataset({"gtco3": xr.DataArray(np.full((1, 2, 2), 0.5, dtype=np.float32), dims=["time", "latitude", "longitude"])})
        return None

    monkeypatch.setattr(provider, "_load_from_remote_url", _load)

    merged = provider._load_from_remote_s3_base("s3://eodata/CAMS/GLOBAL", datetime(2024, 1, 1))

    assert merged is not None
    assert set(merged.data_vars) == {"aod550", "gtco3"}


def test_cams_partial_cdse_dataset_is_supplemented_from_jasmin(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    provider = CAMSProvider("s3://eodata/CAMS/GLOBAL", download_missing=True, cache_dir=tmp_path / "cache")
    primary = xr.Dataset(
        {
            "aod550": xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["latitude", "longitude"]),
            "gtco3": xr.DataArray(np.full((1, 1), 0.5, dtype=np.float32), dims=["latitude", "longitude"]),
        }
    )
    jasmin = xr.Dataset(
        {
            "tcwv": xr.DataArray(np.full((1, 1), 2.5, dtype=np.float32), dims=["latitude", "longitude"]),
        }
    )
    calls: list[str] = []

    def _load_remote(base_url: str, obs_time: datetime):  # noqa: ARG001
        calls.append(base_url)
        if base_url == "s3://eodata/CAMS/GLOBAL":
            return primary
        if base_url == provider._JASMIN_CAMS_BASE_URL:
            return jasmin
        return None

    monkeypatch.setattr(provider, "_load_from_explicit_path", lambda _path: None)
    monkeypatch.setattr(provider, "_load_from_remote_base", _load_remote)
    monkeypatch.setattr(provider, "_load_cds_dataset", lambda _obs_time: None)

    loaded = provider._load_cams_data(datetime(2024, 1, 1))

    assert loaded is not None
    assert set(loaded.data_vars) == {"aod550", "gtco3", "tcwv"}
    assert calls == ["s3://eodata/CAMS/GLOBAL", provider._JASMIN_CAMS_BASE_URL]


def test_cams_missing_cdse_day_falls_back_to_jasmin(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    provider = CAMSProvider("s3://eodata/CAMS/GLOBAL", download_missing=True, cache_dir=tmp_path / "cache")
    jasmin = xr.Dataset(
        {
            "aod550": xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["latitude", "longitude"]),
            "tcwv": xr.DataArray(np.full((1, 1), 2.5, dtype=np.float32), dims=["latitude", "longitude"]),
            "gtco3": xr.DataArray(np.full((1, 1), 0.5, dtype=np.float32), dims=["latitude", "longitude"]),
        }
    )
    calls: list[str] = []

    def _load_remote(base_url: str, obs_time: datetime):  # noqa: ARG001
        calls.append(base_url)
        if base_url == provider._JASMIN_CAMS_BASE_URL:
            return jasmin
        return None

    monkeypatch.setattr(provider, "_load_from_explicit_path", lambda _path: None)
    monkeypatch.setattr(provider, "_load_from_remote_base", _load_remote)
    monkeypatch.setattr(provider, "_load_cds_dataset", lambda _obs_time: None)

    loaded = provider._load_cams_data(datetime(2024, 1, 1))

    assert loaded is not None
    assert loaded.equals(jasmin)
    assert calls == ["s3://eodata/CAMS/GLOBAL", provider._JASMIN_CAMS_BASE_URL]


def test_cams_missing_variables_can_fall_back_to_cds_download(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    provider = CAMSProvider("s3://eodata/CAMS/GLOBAL", download_missing=True, cache_dir=tmp_path / "cache")
    primary = xr.Dataset(
        {
            "aod550": xr.DataArray(np.ones((1, 1), dtype=np.float32), dims=["latitude", "longitude"]),
        }
    )
    cds = xr.Dataset(
        {
            "tcwv": xr.DataArray(np.full((1, 1), 2.5, dtype=np.float32), dims=["latitude", "longitude"]),
            "gtco3": xr.DataArray(np.full((1, 1), 0.5, dtype=np.float32), dims=["latitude", "longitude"]),
        }
    )

    monkeypatch.setattr(provider, "_load_from_explicit_path", lambda _path: None)
    monkeypatch.setattr(
        provider,
        "_load_from_remote_base",
        lambda base_url, obs_time: primary if base_url == "s3://eodata/CAMS/GLOBAL" else None,  # noqa: ARG005
    )
    monkeypatch.setattr(provider, "_load_cds_dataset", lambda _obs_time: cds)

    loaded = provider._load_cams_data(datetime(2024, 1, 1))

    assert loaded is not None
    assert set(loaded.data_vars) == {"aod550", "tcwv", "gtco3"}


def test_cams_missing_variable_alignment_preserves_primary_grid(tmp_path: Path) -> None:
    primary = xr.Dataset(
        {
            "aod550": xr.DataArray(
                np.ones((2, 2), dtype=np.float32),
                dims=["latitude", "longitude"],
                coords={"latitude": [1.0, 0.0], "longitude": [0.0, 1.0]},
            ),
            "gtco3": xr.DataArray(
                np.full((2, 2), 0.5, dtype=np.float32),
                dims=["latitude", "longitude"],
                coords={"latitude": [1.0, 0.0], "longitude": [0.0, 1.0]},
            ),
        }
    )
    fallback = xr.Dataset(
        {
            "tcwv": xr.DataArray(
                np.full((2, 2), 2.5, dtype=np.float32),
                dims=["latitude", "longitude"],
                coords={"latitude": [1.5, 0.5], "longitude": [0.5, 1.5]},
            ),
        }
    )

    merged, added = CAMSProvider._merge_missing_variables(primary, fallback)

    assert added == ["tcwv"]
    assert merged["tcwv"].coords["latitude"].identical(primary["aod550"].coords["latitude"])
    assert merged["tcwv"].coords["longitude"].identical(primary["aod550"].coords["longitude"])
    assert np.isfinite(merged["gtco3"].values).all()


def test_cams_remote_missing_can_fallback_to_download(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    remote = CAMSProvider(
        "https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/",
        download_missing=True,
        cache_dir=tmp_path / "cache",
    )

    def _missing_cache(_url, _storage_options=None):
        raise FileNotFoundError("missing")

    monkeypatch.setattr(remote, "_cache_remote_file", _missing_cache)
    monkeypatch.setattr(
        remote,
        "_complete_cams_dataset",
        lambda dataset, _obs_time: "downloaded" if dataset is None else dataset,
    )

    assert remote._load_cams_data(datetime(2024, 1, 1)) == "downloaded"


def test_cams_extract_supports_forecast_style_time_axes(tmp_path: Path) -> None:
    provider = CAMSProvider(tmp_path, temporal_interp=False)
    ds = xr.Dataset(
        {
            "aod550": xr.DataArray(
                np.array(
                    [
                        [[[0.2], [0.2], [0.2]]],
                        [[[0.4], [0.4], [0.4]]],
                    ],
                    dtype=np.float32,
                ),
                dims=["forecast_period", "forecast_reference_time", "latitude", "longitude"],
                coords={
                    "forecast_period": np.array([0, 3], dtype="timedelta64[h]"),
                    "forecast_reference_time": [np.datetime64("2024-01-01T00:00:00")],
                    "valid_time": (
                        ("forecast_reference_time", "forecast_period"),
                        np.array(
                            [[
                                np.datetime64("2024-01-01T00:00:00"),
                                np.datetime64("2024-01-01T03:00:00"),
                            ]],
                        ),
                    ),
                    "latitude": [1.0, 0.0, -1.0],
                    "longitude": [0.0],
                },
            )
        }
    )

    out = provider._extract_variable(
        ds,
        "aod550",
        (0.0, -1.0, 1.0, 1.0),
        "EPSG:4326",
        1.0,
        datetime(2024, 1, 1, 1),
    )

    assert out.ndim == 2
    assert out.dims == ("latitude", "longitude")
    assert float(out.mean()) == pytest.approx(0.2)


def test_cams_extract_transforms_projected_bounds(tmp_path: Path) -> None:
    provider = CAMSProvider(tmp_path)
    ds = xr.Dataset(
        {
            "aod550": xr.DataArray(
                np.array(
                    [
                        [0.1, 0.2, 0.3],
                        [0.4, 0.5, 0.6],
                        [0.7, 0.8, 0.9],
                    ],
                    dtype=np.float32,
                ),
                dims=["latitude", "longitude"],
                coords={
                    "latitude": [17.0, 16.6, 16.2],
                    "longitude": [117.0, 117.4, 117.8],
                },
            )
        }
    )

    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32650", always_xy=True)
    west, south = transformer.transform(117.3, 16.3)
    east, north = transformer.transform(117.7, 16.7)
    out = provider._extract_variable(
        ds,
        "aod550",
        (west, south, east, north),
        "EPSG:32650",
        10.0,
        datetime(2024, 1, 1),
    )

    assert out.dims == ("latitude", "longitude")
    assert out.shape == (1, 1)
    assert float(out.values[0, 0]) == pytest.approx(0.5)



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
