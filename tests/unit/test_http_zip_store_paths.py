"""Coverage lifts for rt.lut.http_zip_store branch-heavy helpers."""

from __future__ import annotations

import asyncio
import sys
from types import SimpleNamespace
from typing import TYPE_CHECKING

import pytest

import siac.algorithms.rt.lut.http_zip_store as zip_store

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


class _SessionProbeContentRange:
    def head(self, path, allow_redirects=True, headers=None, timeout=None):  # noqa: ARG002
        return _Resp(ok=False, status_code=403, headers={})

    def get(self, path, headers=None, timeout=None):  # noqa: ARG002
        hdrs = dict(headers or {})
        if "Range" in hdrs:
            return _Resp(
                ok=True, status_code=206, headers={"Content-Range": "bytes 0-0/1234"}, content=b"x"
            )
        return _Resp(
            ok=True, status_code=200, headers={"Content-Length": "1234"}, content=b"x" * 1234
        )

    def close(self):
        pass


class _SessionRangeFallback:
    def __init__(self, payload: bytes):
        self.payload = payload

    def head(self, path, allow_redirects=True, headers=None, timeout=None):  # noqa: ARG002
        return _Resp(ok=True, status_code=200, headers={"Content-Length": str(len(self.payload))})

    def get(self, path, headers=None, timeout=None):  # noqa: ARG002
        hdrs = dict(headers or {})
        if "Range" in hdrs:
            return _Resp(ok=False, status_code=416, headers={}, content=b"")
        return _Resp(ok=True, status_code=200, headers={}, content=self.payload)

    def close(self):
        pass


class _TinyBytesFS:
    async def _cat_file(self, path, start=None, end=None, **kwargs):  # noqa: ARG002
        data = b"abc"
        s = 0 if start is None else start
        e = len(data) if end is None else end
        return data[s:e]


class _FakeZipFS:
    def __init__(self, base, path):
        self.base = base
        self.path = path
        self.loop = asyncio.new_event_loop()

    async def _initialize(self):
        return


def test_http_range_fs_probe_and_fallback_branches(monkeypatch):
    import requests

    monkeypatch.setattr(requests, "Session", lambda: _SessionProbeContentRange())
    fs = zip_store._HTTPRangeFileSystem(timeout=2.0)
    size = fs._discover_size("https://example.com/x.zip")
    assert size == 1234

    # start>=end branch
    out_empty = asyncio.run(fs._cat_file("https://example.com/x.zip", start=2, end=2))
    assert out_empty == b""
    fs.close()

    payload = b"abcdefghijklmnopqrstuvwxyz"
    monkeypatch.setattr(requests, "Session", lambda: _SessionRangeFallback(payload))
    fs2 = zip_store._HTTPRangeFileSystem(timeout=2.0)
    fs2._size_cache["https://example.com/y.zip"] = len(payload)
    out = asyncio.run(fs2._cat_file("https://example.com/y.zip", start=3, end=8))
    assert out == payload[3:8]
    fs2.close()


def test_local_fs_wrapper_directory_detail_and_file_path(tmp_path: Path):
    d = tmp_path / "d"
    d.mkdir()
    f = d / "a.bin"
    f.write_bytes(b"0123")

    fs = zip_store._build_local_filesystem()
    detail = asyncio.run(fs._ls(str(d), detail=True))
    assert detail and detail[0]["type"] == "file"

    one = asyncio.run(fs._ls(str(f), detail=True))
    assert one[0]["name"].endswith("a.bin")

    with pytest.raises(FileNotFoundError):
        asyncio.run(fs._ls(str(tmp_path / "missing"), detail=True))


def test_read_only_zip_ls_info_cat_error_branches(monkeypatch):
    zfs = zip_store._ReadOnlyZipFileSystem(_TinyBytesFS(), "dummy.zip")
    zfs._files = {"": {"children": ["a"]}, "a": {"size": 3, "offset": 0}}

    async def _noop():
        return None

    monkeypatch.setattr(zfs, "_initialize", _noop)

    # listing detail=False exercises string listing path
    listing = asyncio.run(zfs._ls("", detail=False))
    assert listing == ["/a"]

    info = asyncio.run(zfs._info("a"))
    assert info["type"] == "file"

    # zero-length slice path
    assert asyncio.run(zfs._cat_file("a", start=1, end=1)) == b""

    with pytest.raises(FileNotFoundError):
        asyncio.run(zfs._ls("missing", detail=True))
    with pytest.raises(FileNotFoundError):
        asyncio.run(zfs._info("missing"))


def test_parse_local_path_and_detect_prefix_helpers():
    p = zip_store._parse_local_path("file:///tmp/demo.zip")
    assert str(p).endswith("/tmp/demo.zip")

    assert zip_store._detect_zarr_prefix_from_names(["foo/bar.txt"]) == ""
    assert zip_store._detect_zarr_prefix_from_names(["nested/.zgroup", "nested/x/0"]) == "nested"


def test_build_s3_fs_and_mapper_validation_paths(monkeypatch):
    captured = {}

    class _FakeS3:
        def __init__(self, **kwargs):
            captured.update(kwargs)

    monkeypatch.setitem(sys.modules, "s3fs", SimpleNamespace(S3FileSystem=_FakeS3))

    _ = zip_store._build_s3_filesystem(
        {
            "client_kwargs": {
                "region_name": "eu-west-1",
                "endpoint_url": "https://example.invalid",
            },
            "anon": False,
            "key": "AK",
            "secret": "SK",
        }
    )
    assert captured["client_kwargs"]["region_name"] == "eu-west-1"
    assert captured["client_kwargs"]["endpoint_url"] == "https://example.invalid"
    assert captured["anon"] is False
    assert captured["key"] == "AK"
    assert captured["secret"] == "SK"

    with pytest.raises(TypeError, match="Top-level S3 storage option"):
        zip_store._build_s3_filesystem({"region": "eu-west-1"})

    with pytest.raises(TypeError, match="Unsupported HTTP zip storage option"):
        zip_store.build_readonly_zip_mapper("https://example.com/lut.zip", {"foo": 1})

    # s3 branch + mapper build path.
    monkeypatch.setattr(zip_store, "_build_s3_filesystem", lambda _options: "fake-base")
    monkeypatch.setattr(zip_store, "_ReadOnlyZipFileSystem", _FakeZipFS)
    monkeypatch.setattr(zip_store, "_detect_zarr_prefix", lambda _zfs: "")

    def _fake_fsmap(root, fs, check=False, create=False):  # noqa: ANN001
        _ = (check, create)
        return {"root": root, "fs": fs}

    monkeypatch.setattr(zip_store, "FSMap", _fake_fsmap)
    m = zip_store.build_readonly_zip_mapper("s3://bucket/key.zip", {})
    assert m["root"] == ""

    # ValueError passthrough path (message does not match compressed marker).
    monkeypatch.setattr(
        zip_store,
        "_detect_zarr_prefix",
        lambda _zfs: (_ for _ in ()).throw(ValueError("other parse error")),
    )
    with pytest.raises(ValueError, match="other parse error"):
        zip_store.build_readonly_zip_mapper("/tmp/lut.zip", {})
