"""Coverage tests for the Zarr LUT backend."""

from __future__ import annotations

import asyncio
import shutil
import struct
import sys
from io import BytesIO
from pathlib import Path
from types import SimpleNamespace
import zipfile

import numpy as np
import pytest
import xarray as xr

from siac.core.types import AtmosphericState, GeometryAngles, SensorBand
from siac.rt.lut.zarr_lut import (
    ZarrLUTBackend,
    _HTTPRangeFileSystem,
    _HTTPZipReadOnlyStore,
    _ReadOnlyZipFileSystem,
    create_lut_from_py6s,
)


def _geometry(shape: tuple[int, int] = (3, 3)) -> GeometryAngles:
    sza = xr.DataArray(np.full(shape, np.deg2rad(30.0), dtype=np.float32), dims=["y", "x"])
    saa = xr.DataArray(np.full(shape, np.deg2rad(150.0), dtype=np.float32), dims=["y", "x"])
    vza = xr.DataArray(np.full(shape, np.deg2rad(10.0), dtype=np.float32), dims=["y", "x"])
    vaa = xr.DataArray(np.full(shape, np.deg2rad(60.0), dtype=np.float32), dims=["y", "x"])
    return GeometryAngles(sza=sza, saa=saa, vza=vza, vaa=vaa)


def _atmo(shape: tuple[int, int] = (3, 3)) -> AtmosphericState:
    return AtmosphericState(
        aot=xr.DataArray(np.full(shape, 0.2, dtype=np.float32), dims=["y", "x"]),
        tcwv=xr.DataArray(np.full(shape, 2.0, dtype=np.float32), dims=["y", "x"]),
        tco3=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
        aot_unc=xr.DataArray(np.full(shape, 0.05, dtype=np.float32), dims=["y", "x"]),
        tcwv_unc=xr.DataArray(np.full(shape, 0.3, dtype=np.float32), dims=["y", "x"]),
        tco3_unc=xr.DataArray(np.full(shape, 0.01, dtype=np.float32), dims=["y", "x"]),
        elevation=xr.DataArray(np.full(shape, 0.1, dtype=np.float32), dims=["y", "x"]),
    )


def _write_small_lut(path: Path, consolidated: bool = True) -> Path:
    sza = np.array([0.0, 30.0, 60.0], dtype=np.float32)
    vza = np.array([0.0, 10.0, 30.0], dtype=np.float32)
    raa = np.array([0.0, 90.0, 180.0], dtype=np.float32)
    aot = np.array([0.1, 0.5], dtype=np.float32)
    tcwv = np.array([1.0, 3.0], dtype=np.float32)
    wl = np.array([490.0, 560.0], dtype=np.float32)

    shape = (len(sza), len(vza), len(raa), len(aot), len(tcwv), len(wl))
    # Deterministic values with non-zero transmittance
    base = np.ones(shape, dtype=np.float32)
    path_ref = base * 0.02
    trans_down = base * 0.85
    trans_up = base * 0.90
    sph_alb = base * 0.03

    ds = xr.Dataset(
        {
            "path_reflectance": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], path_ref),
            "transmittance_down": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], trans_down),
            "transmittance_up": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], trans_up),
            "spherical_albedo": (["sza", "vza", "raa", "aot", "tcwv", "wavelength"], sph_alb),
        },
        coords={"sza": sza, "vza": vza, "raa": raa, "aot": aot, "tcwv": tcwv, "wavelength": wl},
    )
    ds.to_zarr(path, mode="w", consolidated=consolidated)
    return path


class TestZarrLUTBackend:
    def test_missing_lut_raises(self, tmp_path: Path):
        b = ZarrLUTBackend(tmp_path / "missing.zarr")
        with pytest.raises(FileNotFoundError):
            _ = b.lut

    def test_linear_and_nearest_paths(self, tmp_path: Path):
        lut_path = _write_small_lut(tmp_path / "lut.zarr")
        geom = _geometry()
        atmo = _atmo()
        band = SensorBand("B02", 490.0, 65.0, 10.0, 0)

        b_lin = ZarrLUTBackend(lut_path, interpolation_method="linear")
        coeffs_lin = b_lin.compute_coefficients(geom, atmo, band, compute_jacobian=False)
        assert coeffs_lin.xap.shape == (3, 3)
        assert float(coeffs_lin.xap.mean()) > 0.0
        assert b_lin.backend_name == "lut"
        assert b_lin.supports_jacobian()
        assert b_lin.is_available_for_sensor("ANY", "ANY")
        assert b_lin.get_available_wavelengths().size == 2

        b_near = ZarrLUTBackend(lut_path, interpolation_method="nearest")
        coeffs_near = b_near.compute_coefficients(geom, atmo, band, compute_jacobian=False)
        assert coeffs_near.xbp.shape == (3, 3)
        assert np.isfinite(coeffs_near.xbp.values).all()

    def test_jacobian_numerical_path(self, tmp_path: Path):
        lut_path = _write_small_lut(tmp_path / "lut_jac.zarr")
        geom = _geometry((2, 2))
        atmo = _atmo((2, 2))
        band = SensorBand("B03", 560.0, 35.0, 10.0, 1)

        b = ZarrLUTBackend(lut_path, interpolation_method="linear")
        coeffs = b.compute_coefficients(geom, atmo, band, compute_jacobian=True)
        assert coeffs.d_xap is not None
        assert coeffs.d_xbp is not None
        assert coeffs.d_xcp is not None
        assert coeffs.d_xap.sizes["param"] == 2

    def test_load_local_zipped_zarr(self, tmp_path: Path):
        lut_dir = _write_small_lut(tmp_path / "lut_zip.zarr")
        zip_path = tmp_path / "lut_zip.zarr.zip"
        shutil.make_archive(str(zip_path)[:-4], "zip", root_dir=lut_dir)

        b = ZarrLUTBackend(zip_path, interpolation_method="nearest")
        assert "path_reflectance" in b.lut

    def test_s3_storage_options_are_normalized(self, tmp_path: Path, monkeypatch):
        import fsspec

        lut_path = _write_small_lut(tmp_path / "lut_mapper.zarr")
        captured: dict[str, object] = {}
        original_get_mapper = fsspec.get_mapper

        def _fake_get_mapper(path: str, **kwargs):
            captured["path"] = path
            captured["kwargs"] = kwargs
            return original_get_mapper(str(lut_path))

        monkeypatch.setattr(fsspec, "get_mapper", _fake_get_mapper)

        b = ZarrLUTBackend(
            "s3://example/lut.zarr",
            storage_options={
                "region": "eu-west-2",
                "endpoint_url": "https://s3.eu-west-2.amazonaws.com",
                "anon": True,
            },
        )
        _ = b.lut

        assert captured["path"] == "s3://example/lut.zarr"
        kwargs = captured["kwargs"]
        assert isinstance(kwargs, dict)
        assert kwargs.get("anon") is True
        assert kwargs.get("client_kwargs") == {
            "region_name": "eu-west-2",
            "endpoint_url": "https://s3.eu-west-2.amazonaws.com",
        }

    def test_load_non_consolidated_zarr(self, tmp_path: Path):
        lut_path = _write_small_lut(tmp_path / "lut_no_consolidated.zarr", consolidated=False)
        b = ZarrLUTBackend(lut_path, interpolation_method="nearest")
        assert "path_reflectance" in b.lut

    def test_http_zip_store_uses_custom_zip_mapper(self, monkeypatch):
        import siac.rt.lut.store as lut_store

        sentinel_store = {"dummy": b"1"}

        monkeypatch.setattr(
            lut_store,
            "build_readonly_zip_mapper",
            lambda path, opts: sentinel_store,
        )

        store = lut_store.build_lut_store("https://example.com/lut.zarr.zip", storage_options={})
        assert store is sentinel_store

    def test_http_zip_headers_validation(self):
        import siac.rt.lut.http_zip_store as zip_store

        with pytest.raises(TypeError):
            zip_store.build_readonly_zip_mapper(
                "https://example.com/lut.zarr.zip",
                {"headers": "not-a-dict"},
            )


class _FakeResponse:
    def __init__(self, *, status_code: int, content: bytes, headers: dict[str, str]):
        self.status_code = status_code
        self.content = content
        self.headers = headers
        self.ok = 200 <= status_code < 300

    def raise_for_status(self):
        if not self.ok:
            raise RuntimeError(f"HTTP {self.status_code}")


class _FakeSession:
    def __init__(self, data: bytes, *, support_range: bool = True):
        self._data = data
        self._support_range = support_range
        self.get_calls = 0

    def head(self, url, allow_redirects=True, headers=None, timeout=None):  # noqa: ARG002
        return _FakeResponse(
            status_code=200,
            content=b"",
            headers={"Content-Length": str(len(self._data))},
        )

    def get(self, url, headers=None, timeout=None):  # noqa: ARG002
        self.get_calls += 1
        headers = headers or {}
        range_header = headers.get("Range")
        if self._support_range and range_header is not None:
            bounds = range_header.removeprefix("bytes=").split("-", 1)
            if bounds[0] == "":
                suffix = int(bounds[1])
                start = max(0, len(self._data) - suffix)
                end = len(self._data) - 1
            else:
                start = int(bounds[0])
                end = int(bounds[1]) if bounds[1] else len(self._data) - 1
            payload = self._data[start:end + 1]
            return _FakeResponse(
                status_code=206,
                content=payload,
                headers={"Content-Range": f"bytes {start}-{end}/{len(self._data)}"},
            )

        return _FakeResponse(status_code=200, content=self._data, headers={})

    def close(self):
        pass


class _BytesRangeFS:
    """In-memory async byte-range source used for ZIP filesystem tests."""

    def __init__(self, payload: bytes):
        self.payload = payload

    async def _cat_file(self, path, start=None, end=None, **kwargs):  # noqa: ARG002
        size = len(self.payload)
        if start is None:
            start_i = 0
        elif start < 0:
            start_i = max(0, size + start)
        else:
            start_i = start

        if end is None:
            end_i = size
        elif end < 0:
            end_i = max(0, size + end)
        else:
            end_i = end

        end_i = min(end_i, size)
        if start_i >= size or end_i <= start_i:
            return b""
        return self.payload[start_i:end_i]


def _build_eocd(
    *,
    disk: int = 0,
    cd_disk: int = 0,
    entries_disk: int = 0,
    entries_total: int = 0,
    cd_size: int = 0,
    cd_offset: int = 0,
    comment_len: int = 0,
) -> bytes:
    return struct.pack(
        "<IHHHHIIH",
        0x06054B50,
        disk,
        cd_disk,
        entries_disk,
        entries_total,
        cd_size,
        cd_offset,
        comment_len,
    )


def _build_cd_header(
    *,
    signature: int = 0x02014B50,
    flag: int = 0,
    compression: int = 0,
    comp_size: int = 1,
    uncomp_size: int = 1,
    fname_len: int = 0,
    extra_len: int = 0,
    comment_len: int = 0,
    disk_start: int = 0,
    rel_offset: int = 100,
) -> bytes:
    return struct.pack(
        "<IHHHHHHIIIHHHHHII",
        signature,
        20,
        20,
        flag,
        compression,
        0,
        0,
        0,
        comp_size,
        uncomp_size,
        fname_len,
        extra_len,
        comment_len,
        disk_start,
        0,
        0,
        rel_offset,
    )


class TestHTTPRangeHelpers:
    def test_http_range_filesystem_with_range_support(self, monkeypatch):
        import requests

        data = b"0123456789abcdef"
        fake = _FakeSession(data, support_range=True)
        monkeypatch.setattr(requests, "Session", lambda: fake)

        fs = _HTTPRangeFileSystem(timeout=10.0)
        out = asyncio.run(fs._cat_file("https://example.com/test.zip", start=2, end=7))
        assert out == b"23456"
        fs.close()
        assert fake.get_calls >= 1

    def test_http_range_filesystem_server_ignores_range(self, monkeypatch):
        import requests

        data = b"abcdefghijklmnop"
        fake = _FakeSession(data, support_range=False)
        monkeypatch.setattr(requests, "Session", lambda: fake)

        fs = _HTTPRangeFileSystem(timeout=10.0)
        out = asyncio.run(fs._cat_file("https://example.com/test.zip", start=1, end=11))
        assert out == data[1:11]
        fs.close()
        assert fake.get_calls == 1

    def test_read_only_zip_filesystem_reads_zarr_entries(self):
        buf = BytesIO()
        with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_STORED) as zf:
            zf.writestr(".zgroup", '{"zarr_format":2}')
            zf.writestr(".zattrs", "{}")
            zf.writestr("arr/", b"")
            zf.writestr(
                "arr/.zarray",
                '{"zarr_format":2,"shape":[1],"chunks":[1],"dtype":"<u1","compressor":null,"fill_value":0,"order":"C","filters":null}',
            )
            zf.writestr("arr/0", b"\x01")

        zip_fs = _ReadOnlyZipFileSystem(_BytesRangeFS(buf.getvalue()), "dummy.zip")
        async def _exercise():
            listing = await zip_fs._ls("", detail=False)
            payload = await zip_fs._cat_file(".zgroup")
            return listing, payload

        ls_root, payload = asyncio.run(_exercise())
        assert "/.zgroup" in ls_root
        assert "/arr" in ls_root
        assert payload.startswith(b"{")

    def test_http_zip_store_reads_zarr_entries(self, monkeypatch):
        import requests

        buf = BytesIO()
        with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_STORED) as zf:
            zf.writestr(".zgroup", '{"zarr_format":2}')
            zf.writestr(".zattrs", "{}")
            zf.writestr("arr/", b"")
            zf.writestr(
                "arr/.zarray",
                '{"zarr_format":2,"shape":[1],"chunks":[1],"dtype":"<u1","compressor":null,"fill_value":0,"order":"C","filters":null}',
            )
            zf.writestr("arr/0", b"\x01")
        fake = _FakeSession(buf.getvalue(), support_range=True)
        monkeypatch.setattr(requests, "Session", lambda: fake)

        store = _HTTPZipReadOnlyStore("https://example.com/lut.zarr.zip")
        assert ".zgroup" in store
        assert "arr/.zarray" in store
        assert store[".zgroup"].startswith(b"{")
        with pytest.raises(KeyError):
            _ = store["missing"]
        store.close()

    def test_zipfs_initialize_error_branches_zip64_and_central_directory(self):
        # ZIP64 sentinel values + too-short tail for ZIP64 EOCD locator.
        tail_short = b"A" * 10 + _build_eocd(
            disk=0xFFFF,
            cd_disk=0xFFFF,
            entries_disk=0xFFFF,
            entries_total=0xFFFF,
            cd_size=0xFFFFFFFF,
            cd_offset=0xFFFFFFFF,
            comment_len=0,
        )
        zip_fs = _ReadOnlyZipFileSystem(_BytesRangeFS(tail_short), "dummy.zip")
        with pytest.raises(ValueError, match="ZIP64 EOCD and locator do not fit"):
            asyncio.run(zip_fs._initialize())

        # ZIP64 markers present but ZIP64 EOCD signature missing in the searched window.
        tail_long = b"A" * 120 + _build_eocd(
            disk=0xFFFF,
            cd_disk=0xFFFF,
            entries_disk=0xFFFF,
            entries_total=0xFFFF,
            cd_size=0xFFFFFFFF,
            cd_offset=0xFFFFFFFF,
            comment_len=0,
        )
        zip_fs2 = _ReadOnlyZipFileSystem(_BytesRangeFS(tail_long), "dummy.zip")
        with pytest.raises(ValueError, match="No ZIP64 EOCD found"):
            asyncio.run(zip_fs2._initialize())

        # Multi-disk central directory is rejected.
        tail_multidisk = b"A" * 10 + _build_eocd(disk=1, cd_disk=0, entries_disk=1, entries_total=1)
        zip_fs3 = _ReadOnlyZipFileSystem(_BytesRangeFS(tail_multidisk), "dummy.zip")
        with pytest.raises(ValueError, match="Unsupported multi-disk"):
            asyncio.run(zip_fs3._initialize())

    def test_zipfs_initialize_error_branches_cd_entries(self):
        class _TwoPhaseFS:
            def __init__(self, tail: bytes, cd_payload: bytes):
                self.tail = tail
                self.cd_payload = cd_payload

            async def _cat_file(self, path, start=None, end=None, **kwargs):  # noqa: ANN001, ARG002
                if start is not None and start < 0:
                    return self.tail
                return self.cd_payload

        # cd_size > returned bytes
        tail = b"A" * 20 + _build_eocd(entries_disk=1, entries_total=1, cd_size=50, cd_offset=0)
        fs_short = _TwoPhaseFS(tail, b"x" * 10)
        zip_fs = _ReadOnlyZipFileSystem(fs_short, "dummy.zip")
        with pytest.raises(ValueError, match="Failed to read central directory"):
            asyncio.run(zip_fs._initialize())

        # Truncated entry (< 46 bytes)
        tail2 = b"A" * 20 + _build_eocd(entries_disk=1, entries_total=1, cd_size=20, cd_offset=0)
        fs_trunc_entry = _TwoPhaseFS(tail2, b"x" * 20)
        zip_fs2 = _ReadOnlyZipFileSystem(fs_trunc_entry, "dummy.zip")
        with pytest.raises(ValueError, match="Truncated central directory entry"):
            asyncio.run(zip_fs2._initialize())

        # Invalid central directory signature.
        bad_sig = _build_cd_header(signature=0x12345678)
        tail3 = b"A" * 20 + _build_eocd(entries_disk=1, entries_total=1, cd_size=len(bad_sig), cd_offset=0)
        fs_bad_sig = _TwoPhaseFS(tail3, bad_sig)
        zip_fs3 = _ReadOnlyZipFileSystem(fs_bad_sig, "dummy.zip")
        with pytest.raises(ValueError, match="Invalid central directory signature"):
            asyncio.run(zip_fs3._initialize())

    def test_zipfs_initialize_error_branches_filename_and_extra(self):
        class _TwoPhaseFS:
            def __init__(self, tail: bytes, cd_payload: bytes):
                self.tail = tail
                self.cd_payload = cd_payload

            async def _cat_file(self, path, start=None, end=None, **kwargs):  # noqa: ANN001, ARG002
                if start is not None and start < 0:
                    return self.tail
                return self.cd_payload

        def _run_with_cd(cd_payload: bytes, msg: str) -> None:
            tail = b"A" * 20 + _build_eocd(
                entries_disk=1,
                entries_total=1,
                cd_size=len(cd_payload),
                cd_offset=0,
            )
            zip_fs = _ReadOnlyZipFileSystem(_TwoPhaseFS(tail, cd_payload), "dummy.zip")
            with pytest.raises(ValueError, match=msg):
                asyncio.run(zip_fs._initialize())

        # filename bytes missing
        _run_with_cd(_build_cd_header(fname_len=5), "Truncated filename")

        # extra field body missing
        cd_extra_short = _build_cd_header(fname_len=1, extra_len=5) + b"a"
        _run_with_cd(cd_extra_short, "Truncated extra field")

        # extra field header missing (need at least 4 bytes)
        cd_extra_header = _build_cd_header(fname_len=1, extra_len=3) + b"a" + b"\x01\x02\x03"
        _run_with_cd(cd_extra_header, "Truncated extra field header")

        # extra field payload missing after valid tag/size
        cd_extra_payload = _build_cd_header(fname_len=1, extra_len=5) + b"a" + struct.pack("<HH", 1, 5) + b"x"
        _run_with_cd(cd_extra_payload, "Truncated extra field payload")

    def test_zipfs_parent_directory_construction_branch(self):
        buf = BytesIO()
        with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_STORED) as zf:
            # No explicit directory entries; parent dirs are inferred from file names.
            zf.writestr("nested/dir/a.bin", b"a")
            zf.writestr("nested/dir/b.bin", b"b")

        zip_fs = _ReadOnlyZipFileSystem(_BytesRangeFS(buf.getvalue()), "dummy.zip")
        asyncio.run(zip_fs._initialize())
        assert zip_fs._files is not None
        assert "nested" in zip_fs._files
        assert "nested/dir" in zip_fs._files
        assert "nested/dir/a.bin" in zip_fs._files


class _FallbackSession:
    """Session stub that forces HTTP range fallback branches."""

    def __init__(self, payload: bytes):
        self.payload = payload
        self.calls: list[tuple[str, dict[str, str]]] = []

    def head(self, url, allow_redirects=True, headers=None, timeout=None):  # noqa: ARG002
        return _FakeResponse(status_code=403, content=b"", headers={})

    def get(self, url, headers=None, timeout=None):  # noqa: ARG002
        hdrs = dict(headers or {})
        self.calls.append((url, hdrs))
        if "Range" in hdrs:
            return _FakeResponse(status_code=403, content=b"", headers={})
        return _FakeResponse(status_code=200, content=self.payload, headers={})

    def close(self):
        pass


class TestZipStoreUtilities:
    def test_helper_path_and_slice_utilities(self):
        import siac.rt.lut.http_zip_store as zip_store

        assert zip_store._normalize_fs_path("") == ""
        assert zip_store._normalize_fs_path("/") == ""
        assert zip_store._normalize_fs_path("a/b/") == "a/b"

        assert zip_store._slice_bounds(10, None, None) == (0, 10)
        assert zip_store._slice_bounds(10, -3, None) == (7, 10)
        assert zip_store._slice_bounds(10, 3, -2) == (3, 8)
        assert zip_store._slice_bounds(10, 20, None) == (10, 10)

    def test_http_range_filesystem_full_body_fallback(self, monkeypatch):
        import requests

        payload = b"abcdefghijklmnopqrstuvwxyz"
        fake = _FallbackSession(payload)
        monkeypatch.setattr(requests, "Session", lambda: fake)

        fs = _HTTPRangeFileSystem(timeout=5.0)
        out = asyncio.run(fs._cat_file("https://example.com/fallback.zip", start=2, end=9))
        assert out == payload[2:9]
        info = asyncio.run(fs._info("https://example.com/fallback.zip"))
        assert info["size"] == len(payload)
        listing = asyncio.run(fs._ls("https://example.com/fallback.zip", detail=False))
        assert listing == ["https://example.com/fallback.zip"]
        fs.close()

    def test_local_range_filesystem(self, tmp_path: Path):
        import siac.rt.lut.http_zip_store as zip_store

        file_path = tmp_path / "x.bin"
        file_path.write_bytes(b"0123456789")

        local_fs = zip_store._LocalRangeFileSystem()
        info = asyncio.run(local_fs._info(str(file_path)))
        assert info["type"] == "file"
        assert info["size"] == 10

        listing = asyncio.run(local_fs._ls(str(tmp_path), detail=False))
        assert str(file_path) in listing

        out = asyncio.run(local_fs._cat_file(str(file_path), start=-4, end=None))
        assert out == b"6789"

        with pytest.raises(FileNotFoundError):
            asyncio.run(local_fs._info(str(tmp_path / "missing.bin")))

    def test_read_only_zip_filesystem_error_paths(self):
        buf = BytesIO()
        with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_STORED) as zf:
            zf.writestr(".zgroup", '{"zarr_format":2}')
            zf.writestr("arr/", b"")
            zf.writestr(
                "arr/.zarray",
                '{"zarr_format":2,"shape":[1],"chunks":[1],"dtype":"<u1","compressor":null,"fill_value":0,"order":"C","filters":null}',
            )
        zip_fs = _ReadOnlyZipFileSystem(_BytesRangeFS(buf.getvalue()), "dummy.zip")

        with pytest.raises(FileNotFoundError):
            asyncio.run(zip_fs._info("missing"))
        with pytest.raises(FileNotFoundError):
            asyncio.run(zip_fs._cat_file("arr"))

    def test_s3_parse_and_import_error(self):
        import siac.rt.lut.http_zip_store as zip_store

        assert zip_store._parse_s3_url("s3://bucket/key/path") == ("bucket", "key/path")
        with pytest.raises(ValueError):
            zip_store._parse_s3_url("http://bucket/key")
        with pytest.raises(ValueError):
            zip_store._parse_s3_url("s3://bucket-only")

        try:
            fs = zip_store._build_s3_filesystem({})
        except ImportError:
            # Expected in lightweight test envs where s3fs is not installed.
            return
        assert fs is not None

    def test_s3_builder_defaults_to_ambient_credentials(self, monkeypatch):
        import siac.rt.lut.http_zip_store as zip_store

        captured: dict[str, object] = {}

        class _FakeS3FileSystem:
            def __init__(self, **kwargs):
                captured.update(kwargs)

        monkeypatch.setitem(sys.modules, "s3fs", SimpleNamespace(S3FileSystem=_FakeS3FileSystem))
        _ = zip_store._build_s3_filesystem({})

        assert captured.get("asynchronous") is True
        assert "anon" not in captured
        assert "key" not in captured
        assert "secret" not in captured

    def test_build_mapper_validation_and_remote_compressed_branch(self, monkeypatch):
        import siac.rt.lut.http_zip_store as zip_store

        with pytest.raises(TypeError):
            zip_store.build_readonly_zip_mapper(
                "https://example.com/lut.zip",
                {"headers": "bad"},
            )
        with pytest.raises(TypeError):
            zip_store.build_readonly_zip_mapper(
                "/tmp/lut.zip",
                {"unexpected": 1},
            )

        monkeypatch.setattr(
            zip_store,
            "_detect_zarr_prefix",
            lambda fs: (_ for _ in ()).throw(ValueError("not stored (uncompressed)")),
        )
        with pytest.raises(ValueError):
            zip_store.build_readonly_zip_mapper(
                "https://example.com/lut.zip",
                {"timeout": 1.0},
            )

    def test_local_compressed_zip_fallback_detects_nested_zarr_root(self, tmp_path: Path):
        import siac.rt.lut.http_zip_store as zip_store

        zip_path = tmp_path / "lut_nested.zarr.zip"
        with zipfile.ZipFile(zip_path, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("lut.zarr/.zgroup", '{"zarr_format":2}')
            zf.writestr("lut.zarr/.zattrs", "{}")

        mapper = zip_store.build_readonly_zip_mapper(str(zip_path), {})
        assert ".zgroup" in mapper
        assert mapper[".zgroup"].startswith(b"{")

    def test_http_zip_store_read_only_methods(self, monkeypatch):
        import siac.rt.lut.http_zip_store as zip_store

        class _FakeMapper(dict):
            fs = None

        monkeypatch.setattr(zip_store, "build_readonly_zip_mapper", lambda path, options: _FakeMapper())
        store = zip_store._HTTPZipReadOnlyStore("https://example.com/lut.zip", timeout=1.0)
        with pytest.raises(TypeError):
            store["x"] = b"1"
        with pytest.raises(TypeError):
            del store["x"]
        store.close()

    def test_read_only_zip_listing_detail_and_close(self):
        import siac.rt.lut.http_zip_store as zip_store

        buf = BytesIO()
        with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_STORED) as zf:
            zf.writestr(".zgroup", '{"zarr_format":2}')
            zf.writestr("arr/", b"")
            zf.writestr(
                "arr/.zarray",
                '{"zarr_format":2,"shape":[1],"chunks":[1],"dtype":"<u1","compressor":null,"fill_value":0,"order":"C","filters":null}',
            )
            zf.writestr("arr/0", b"\x01")

        class _ClosableBytesFS(_BytesRangeFS):
            def __init__(self, payload: bytes):
                super().__init__(payload)
                self.closed = False

            def close(self):
                self.closed = True

        fs = _ClosableBytesFS(buf.getvalue())
        zip_fs = zip_store._ReadOnlyZipFileSystem(fs, "dummy.zip")

        async def _exercise():
            root = await zip_fs._ls("", detail=True)
            one = await zip_fs._ls(".zgroup", detail=True)
            return root, one

        root_listing, file_listing = asyncio.run(_exercise())
        assert any(item["type"] == "directory" for item in root_listing)
        assert file_listing[0]["type"] == "file"
        zip_fs.close()
        assert fs.closed is True

    def test_store_helper_paths_and_local_mapper_branch(self, monkeypatch, tmp_path: Path):
        import siac.rt.lut.store as lut_store

        local_file_uri = (tmp_path / "lut.zarr").as_uri()
        local_path = lut_store.as_local_path(local_file_uri)
        assert local_path is not None
        assert str(local_path).endswith("lut.zarr")

        opts = lut_store.normalize_storage_options(
            "s3://bucket/key",
            {
                "region": "eu-west-1",
                "endpoint_url": "https://example.invalid",
                "client_kwargs": {"region_name": "keep-me"},
            },
        )
        assert opts["client_kwargs"]["region_name"] == "keep-me"
        assert opts["client_kwargs"]["endpoint_url"] == "https://example.invalid"

        captured: dict[str, object] = {}
        monkeypatch.setattr(
            "fsspec.get_mapper",
            lambda path, **kwargs: captured.update({"path": path, "kwargs": kwargs}) or {"ok": True},
        )
        out = lut_store.build_lut_store(str(tmp_path / "x.zarr"), {"anon": True})
        assert out == {"ok": True}
        assert captured["path"] == str(tmp_path / "x.zarr")
        assert captured["kwargs"] == {"anon": True}


class TestCreateLUTFromPy6S:
    def test_import_error_branch(self, tmp_path: Path):
        saved = sys.modules.pop("Py6S", None)
        try:
            with pytest.raises(ImportError):
                create_lut_from_py6s(tmp_path / "x.zarr", wavelengths=[500.0])
        finally:
            if saved is not None:
                sys.modules["Py6S"] = saved

    def test_create_lut_with_fake_py6s(self, tmp_path: Path, monkeypatch):
        class _FakeGeometry:
            @staticmethod
            def User():
                return SimpleNamespace()

        class _FakeWavelength:
            def __init__(self, value):
                self.value = value

        class _FakeAtmosProfile:
            @staticmethod
            def UserWaterAndOzone(tcwv, ozone):
                return (tcwv, ozone)

        class _FakeAeroProfile:
            Continental = "continental"
            Maritime = "maritime"

        class _FakeSixS:
            def __init__(self):
                self.outputs = None
                self.geometry = None
                self.aot550 = 0.0
                self.wavelength = None

            def run(self):
                wl = float(getattr(self.wavelength, "value", 0.5))
                val = 0.01 + 0.1 * self.aot550 + 0.001 * wl
                self.outputs = SimpleNamespace(
                    atmospheric_intrinsic_reflectance=val,
                    transmittance_total_scattering=SimpleNamespace(
                        downward=0.8,
                        upward=0.85,
                    ),
                    spherical_albedo=0.03,
                )

        fake_module = SimpleNamespace(
            SixS=_FakeSixS,
            AtmosProfile=_FakeAtmosProfile,
            AeroProfile=_FakeAeroProfile,
            Geometry=_FakeGeometry,
            Wavelength=_FakeWavelength,
        )
        monkeypatch.setitem(sys.modules, "Py6S", fake_module)

        out = tmp_path / "fake_lut.zarr"
        create_lut_from_py6s(
            output_path=out,
            wavelengths=[500.0],
            sza_range=(0, 1, 1),
            vza_range=(0, 1, 1),
            raa_range=(0, 1, 1),
            aot_values=[0.1],
            tcwv_values=[1.0],
            aerosol_type="unknown",  # hits fallback branch
            ozone=0.3,
        )

        ds = xr.open_zarr(out, consolidated=True)
        assert "path_reflectance" in ds
        assert ds["path_reflectance"].shape == (1, 1, 1, 1, 1, 1)
        assert ds.attrs["aerosol_type"] == "unknown"

    def test_create_lut_defaults_maritime_and_progress(self, tmp_path: Path, monkeypatch):
        import siac.rt.lut.create as create_mod

        class _FakeGeometry:
            @staticmethod
            def User():
                return SimpleNamespace()

        class _FakeWavelength:
            def __init__(self, value):
                self.value = value

        class _FakeAtmosProfile:
            @staticmethod
            def UserWaterAndOzone(tcwv, ozone):
                return (tcwv, ozone)

        class _FakeAeroProfile:
            Continental = "continental"
            Maritime = "maritime"

        class _FakeSixS:
            def __init__(self):
                self.outputs = None
                self.geometry = None
                self.aot550 = 0.0
                self.wavelength = None

            def run(self):
                val = 0.02 + 0.1 * self.aot550 + 0.001 * float(getattr(self.wavelength, "value", 0.5))
                self.outputs = SimpleNamespace(
                    atmospheric_intrinsic_reflectance=val,
                    transmittance_total_scattering=SimpleNamespace(downward=0.8, upward=0.85),
                    spherical_albedo=0.03,
                )

        fake_module = SimpleNamespace(
            SixS=_FakeSixS,
            AtmosProfile=_FakeAtmosProfile,
            AeroProfile=_FakeAeroProfile,
            Geometry=_FakeGeometry,
            Wavelength=_FakeWavelength,
        )
        monkeypatch.setitem(sys.modules, "Py6S", fake_module)
        monkeypatch.setattr(create_mod.np, "logspace", lambda *args, **kwargs: np.array([0.2], dtype=np.float32))

        # Defaults path (aot_values/tcwv_values None) + maritime branch.
        create_lut_from_py6s(
            output_path=tmp_path / "defaults.zarr",
            wavelengths=[500.0],
            sza_range=(0, 1, 1),
            vza_range=(0, 1, 1),
            raa_range=(0, 1, 1),
            aot_values=None,
            tcwv_values=None,
            aerosol_type="maritime",
        )

        logged: list[str] = []
        monkeypatch.setattr(
            create_mod.logger,
            "info",
            lambda msg: logged.append(str(msg)),
        )
        create_lut_from_py6s(
            output_path=tmp_path / "progress.zarr",
            wavelengths=[500.0],
            sza_range=(0, 10, 1),
            vza_range=(0, 10, 1),
            raa_range=(0, 10, 1),
            aot_values=[0.1],
            tcwv_values=[1.0],
            aerosol_type="continental",
        )
        assert any("Progress: 1000/1000" in m for m in logged)
