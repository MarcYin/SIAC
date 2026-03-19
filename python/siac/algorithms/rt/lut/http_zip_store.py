"""Read-only ZIP store helpers for local/HTTP/S3 LUT Zarr archives.

This module implements the `ReadOnlyZipFileSystem` approach used in
`archive_scripts/load_remote_rt.py` and adapts it for the LUT backend.
"""

from __future__ import annotations

import asyncio
import posixpath
import struct
from pathlib import Path
from typing import TYPE_CHECKING, Any
from urllib.parse import unquote, urlparse

from fsspec.asyn import AsyncFileSystem, sync
from fsspec.implementations.asyn_wrapper import AsyncFileSystemWrapper
from fsspec.implementations.local import LocalFileSystem
from fsspec.mapping import FSMap

if TYPE_CHECKING:
    from collections.abc import Iterable, MutableMapping


def _normalize_fs_path(path: str) -> str:
    """Normalize an FS path by removing leading/trailing slashes."""
    normalized = posixpath.normpath(path or "")
    if normalized in (".", "/"):
        return ""
    return normalized.lstrip("/").rstrip("/")


def _slice_bounds(size: int, start: int | None, end: int | None) -> tuple[int, int]:
    """Convert optional possibly-negative [start, end) bounds to clamped positives."""
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
        return size, size
    return start_i, end_i


class _HTTPRangeFileSystem(AsyncFileSystem):
    """Minimal HTTP byte-range filesystem for one-object reads."""

    protocol = "httprange"
    cachable = False

    def __init__(
        self,
        *,
        timeout: float = 30.0,
        headers: dict[str, str] | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        import requests

        self.timeout = timeout
        self.headers = {
            "User-Agent": "siac-lut-loader/1.0 (+https://github.com/MarcYin/SIAC)",
        }
        self.headers.update(dict(headers or {}))
        self._session = requests.Session()
        self._size_cache: dict[str, int] = {}
        self._full_body_cache: dict[str, bytes] = {}

    def _discover_size(self, path: str) -> int:
        if path in self._size_cache:
            return self._size_cache[path]

        response = self._session.head(
            path,
            allow_redirects=True,
            headers=self.headers,
            timeout=self.timeout,
        )
        if response.ok and response.headers.get("Content-Length"):
            size = int(response.headers["Content-Length"])
            self._size_cache[path] = size
            return size

        probe_headers = dict(self.headers)
        probe_headers["Range"] = "bytes=0-0"
        response = self._session.get(path, headers=probe_headers, timeout=self.timeout)
        if response.ok:
            content_range = response.headers.get("Content-Range")
            if content_range and "/" in content_range:
                size = int(content_range.rsplit("/", 1)[1])
                self._size_cache[path] = size
                return size
            if response.headers.get("Content-Length"):
                size = int(response.headers["Content-Length"])
                self._size_cache[path] = size
                return size

        # Last resort: regular GET. If no content-length is available we cache the
        # full body once and serve slices from memory.
        response = self._session.get(path, headers=self.headers, timeout=self.timeout)
        response.raise_for_status()
        if response.headers.get("Content-Length"):
            size = int(response.headers["Content-Length"])
            self._size_cache[path] = size
            return size

        payload = response.content
        size = len(payload)
        self._full_body_cache[path] = payload
        self._size_cache[path] = size
        return size

    async def _cat_file(
        self,
        path: str,
        start: int | None = None,
        end: int | None = None,
        **kwargs: Any,
    ) -> bytes:
        del kwargs
        size = self._discover_size(path)
        start_i, end_i = _slice_bounds(size, start, end)
        if start_i >= end_i:
            return b""

        if path in self._full_body_cache:
            return self._full_body_cache[path][start_i:end_i]

        range_headers = dict(self.headers)
        range_headers["Range"] = f"bytes={start_i}-{end_i - 1}"
        response = self._session.get(path, headers=range_headers, timeout=self.timeout)
        if not response.ok:
            # Some endpoints reject ranged GETs. Fallback to full-body fetch.
            response = self._session.get(path, headers=self.headers, timeout=self.timeout)
            response.raise_for_status()
            payload = response.content
            if len(payload) == size:
                self._full_body_cache[path] = payload
            return payload[start_i:end_i]
        payload = response.content

        # Some servers may ignore Range and send full body.
        if response.status_code == 200 and len(payload) == size:
            self._full_body_cache[path] = payload
            return payload[start_i:end_i]
        return payload

    async def _info(self, path: str, **kwargs: Any) -> dict[str, Any]:
        del kwargs
        size = self._discover_size(path)
        return {
            "name": path,
            "type": "file",
            "size": size,
            "created": None,
            "islink": False,
        }

    async def _ls(self, path: str, detail: bool = True, **kwargs: Any) -> list[Any]:
        del kwargs
        info = await self._info(path)
        return [info] if detail else [path]

    def close(self) -> None:
        self._session.close()


def _build_local_filesystem() -> AsyncFileSystem:
    """Return an async wrapper over fsspec's LocalFileSystem."""
    return AsyncFileSystemWrapper(
        fs=LocalFileSystem(auto_mkdir=False),
        asynchronous=True,
    )


class _ReadOnlyZipFileSystem(AsyncFileSystem):
    """Async read-only filesystem for uncompressed ZIP files (ZIP_STORED).

    This class is adapted from the archived implementation used for remote LUT
    access. It indexes the central directory and serves file content via
    `_cat_file` without extracting the archive.
    """

    protocol = "zipfs"
    cachable = False
    MAX_ZIP_TAIL_READ = 64 * 1024

    def __init__(self, fs: AsyncFileSystem, path: str, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.fs = fs
        self.path = path
        self._files: dict[str, dict[str, Any]] | None = None
        self._lock: asyncio.Lock | None = None

    async def _initialize(self) -> None:
        if self._lock is None:
            self._lock = asyncio.Lock()
        async with self._lock:
            if self._files is not None:
                return

            data = await self.fs._cat_file(
                self.path,
                start=-self.MAX_ZIP_TAIL_READ,
                end=None,
            )
            if len(data) < 22:
                raise ValueError(f"EOCD does not fit in {self.path}: {len(data)} bytes")

            eocd_var_lengths = {
                "end_of_central_dir_signature": 4,
                "number_of_this_disk": 2,
                "number_of_disk_with_central_directory": 2,
                "total_entries_on_this_disk": 2,
                "total_entries": 2,
                "central_directory_size": 4,
                "central_directory_offset": 4,
                "zip_file_comment_length": 2,
            }
            zip64_eocd_var_lengths = {
                "zip64_end_of_central_dir_signature": 4,
                "zip64_end_of_central_dir_size": 8,
                "version_made_by": 2,
                "version_needed_to_extract": 2,
                "number_of_this_disk": 4,
                "number_of_disk_with_central_directory": 4,
                "total_entries_on_this_disk": 8,
                "total_entries": 8,
                "central_directory_size": 8,
                "central_directory_offset": 8,
            }
            var_length_to_format = {2: "H", 4: "L", 8: "Q"}

            eocd: dict[str, int | None] = {}
            is_zip64 = False
            eocd_pos = data[:-16].rfind(b"\x50\x4b\x05\x06")
            if eocd_pos == -1:
                raise ValueError(
                    f"No EOCD found in the last {self.MAX_ZIP_TAIL_READ} bytes of {self.path}"
                )
            pos = eocd_pos
            for var, length in eocd_var_lengths.items():
                value = struct.unpack_from(f"<{var_length_to_format[length]}", data, pos)[0]
                if value == 2 ** (length * 8) - 1:
                    value = None
                    is_zip64 = True
                eocd[var] = value
                pos += length

            if is_zip64:
                if len(data) - 22 < 56 + 20:
                    raise ValueError(
                        f"ZIP64 EOCD and locator do not fit in tail read for {self.path}"
                    )
                zip64_pos = data[: eocd_pos - 52].rfind(b"\x50\x4b\x06\x06")
                if zip64_pos == -1:
                    raise ValueError(
                        f"No ZIP64 EOCD found in the last {self.MAX_ZIP_TAIL_READ} bytes of {self.path}"
                    )
                pos = zip64_pos
                for var, length in zip64_eocd_var_lengths.items():
                    eocd[var] = struct.unpack_from(f"<{var_length_to_format[length]}", data, pos)[0]
                    pos += length

            if (
                eocd["number_of_this_disk"] != 0
                or eocd["number_of_disk_with_central_directory"] != 0
                or eocd["total_entries_on_this_disk"] != eocd["total_entries"]
            ):
                raise ValueError(f"Unsupported multi-disk central directory in {self.path}")

            cd_size = eocd["central_directory_size"]
            cd_offset = eocd["central_directory_offset"]
            cd_entries = eocd["total_entries"]

            if cd_size in (None, 0) or cd_offset is None or cd_entries in (None, 0):
                self._files = {"": {"children": []}}
                return

            cd_data = await self.fs._cat_file(
                self.path,
                start=int(cd_offset),
                end=int(cd_offset + cd_size),
            )
            if len(cd_data) != int(cd_size):
                raise ValueError(
                    f"Failed to read central directory: expected {cd_size} bytes, got {len(cd_data)}"
                )

            cd_header_var_lengths = {
                "central_file_header_signature": 4,
                "version_made_by": 2,
                "version_needed_to_extract": 2,
                "general_purpose_bit_flag": 2,
                "compression_method": 2,
                "last_mod_file_time": 2,
                "last_mod_file_date": 2,
                "crc_32": 4,
                "compressed_size": 4,
                "uncompressed_size": 4,
                "file_name_length": 2,
                "extra_field_length": 2,
                "file_comment_length": 2,
                "disk_number_start": 2,
                "internal_file_attributes": 2,
                "external_file_attributes": 4,
                "relative_offset_of_local_header": 4,
            }
            zip64_header_var_lengths = {
                "uncompressed_size": 8,
                "compressed_size": 8,
                "relative_offset_of_local_header": 8,
                "disk_number_start": 4,
            }

            self._files = {"": {"children": []}}
            pos = 0
            previous_file: dict[str, Any] | None = None

            for _ in range(int(cd_entries)):
                if pos + 46 > len(cd_data):
                    raise ValueError(f"Truncated central directory entry in {self.path}")

                header: dict[str, int | None] = {}
                for var, length in cd_header_var_lengths.items():
                    value = struct.unpack_from(f"<{var_length_to_format[length]}", cd_data, pos)[0]
                    if value == 2 ** (length * 8) - 1:
                        value = None
                    header[var] = value
                    pos += length

                if header["central_file_header_signature"] != 0x02014B50:
                    raise ValueError(f"Invalid central directory signature in {self.path}")
                if (
                    header["compression_method"] != 0
                    or header["compressed_size"] != header["uncompressed_size"]
                ):
                    raise ValueError(f"File in {self.path} is not stored (uncompressed)")

                utf8 = (header["general_purpose_bit_flag"] or 0) & 0x800 != 0
                fname_len = int(header["file_name_length"] or 0)
                extra_len = int(header["extra_field_length"] or 0)
                comment_len = int(header["file_comment_length"] or 0)

                if pos + fname_len > len(cd_data):
                    raise ValueError(f"Truncated filename in {self.path}")
                raw_name = cd_data[pos : pos + fname_len]
                name = raw_name.decode("utf-8" if utf8 else "ascii")
                pos += fname_len

                extra_end = pos + extra_len
                if extra_end > len(cd_data):
                    raise ValueError(f"Truncated extra field in {self.path}")
                while pos < extra_end:
                    if pos + 4 > extra_end:
                        raise ValueError(f"Truncated extra field header in {self.path}")
                    tag, size = struct.unpack_from("<HH", cd_data, pos)
                    pos += 4
                    if pos + size > extra_end:
                        raise ValueError(f"Truncated extra field payload in {self.path}")
                    backup_pos = pos
                    if tag == 0x0001:
                        for var, length in zip64_header_var_lengths.items():
                            if header[var] is None:
                                header[var] = struct.unpack_from(
                                    f"<{var_length_to_format[length]}",
                                    cd_data,
                                    pos,
                                )[0]
                                pos += length
                    pos = backup_pos + size

                size = int(header["compressed_size"] or 0)
                offset = int(header["relative_offset_of_local_header"] or 0)
                pos += comment_len

                is_dir = name.endswith("/")
                if is_dir:
                    name = name[:-1]

                if is_dir:
                    entry: dict[str, Any] = {"children": []}
                else:
                    entry = {"size": size, "offset": offset}
                self._files[name] = entry

                if name:
                    parent = posixpath.dirname(name)
                    if parent not in self._files:
                        parts = parent.split("/")
                        cursor = ""
                        for part in parts:
                            cursor = f"{cursor}/{part}" if cursor else part
                            if cursor not in self._files:
                                self._files[cursor] = {"children": []}
                                parent_of_cursor = posixpath.dirname(cursor)
                                if parent_of_cursor in self._files:
                                    self._files[parent_of_cursor]["children"].append(cursor)
                    self._files[parent]["children"].append(name)

                if previous_file is not None and "children" not in previous_file:
                    if offset <= int(previous_file["offset"]):
                        raise ValueError(
                            f"Non-ascending local header offsets in central directory for {self.path}"
                        )
                    previous_file["offset"] = offset - int(previous_file["size"])

                previous_file = entry

            if previous_file is not None and "children" not in previous_file:
                if int(cd_offset) <= int(previous_file["offset"]):
                    raise ValueError(
                        f"Non-ascending local header offsets against central directory for {self.path}"
                    )
                previous_file["offset"] = int(cd_offset) - int(previous_file["size"])

    async def _ls(self, path: str, detail: bool = True, **kwargs: Any) -> list[Any]:
        del kwargs
        await self._initialize()
        assert self._files is not None

        normalized = _normalize_fs_path(path)
        if normalized not in self._files:
            raise FileNotFoundError(f"Path {normalized} not found")

        file = self._files[normalized]

        def to_listing(fname: str, info: dict[str, Any]) -> Any:
            if detail:
                is_dir = "children" in info
                return {
                    "name": f"/{fname}",
                    "type": "directory" if is_dir else "file",
                    "size": 0 if is_dir else int(info["size"]),
                    "created": None,
                    "islink": False,
                }
            return f"/{fname}"

        if "children" in file:
            return [to_listing(child, self._files[child]) for child in file["children"]]
        return [to_listing(normalized, file)]

    async def _info(self, path: str, **kwargs: Any) -> dict[str, Any]:
        del kwargs
        await self._initialize()
        assert self._files is not None

        normalized = _normalize_fs_path(path)
        if normalized not in self._files:
            raise FileNotFoundError(f"Path {normalized} not found")
        info = self._files[normalized]
        is_dir = "children" in info
        return {
            "name": f"/{normalized}",
            "type": "directory" if is_dir else "file",
            "size": 0 if is_dir else int(info["size"]),
            "created": None,
            "islink": False,
        }

    async def _cat_file(
        self,
        path: str,
        start: int | None = None,
        end: int | None = None,
        **kwargs: Any,
    ) -> bytes:
        del kwargs
        await self._initialize()
        assert self._files is not None

        normalized = _normalize_fs_path(path)
        if normalized not in self._files:
            raise FileNotFoundError(f"File {normalized} not found")
        if "children" in self._files[normalized]:
            raise FileNotFoundError(f"{normalized} is a directory")

        info = self._files[normalized]
        size = int(info["size"])
        start_i, end_i = _slice_bounds(size, start, end)
        if start_i >= end_i:
            return b""

        read_start = int(info["offset"]) + start_i
        read_end = read_start + (end_i - start_i)
        return await self.fs._cat_file(self.path, start=read_start, end=read_end)

    def close(self) -> None:
        if hasattr(self.fs, "close"):
            self.fs.close()


def _parse_local_path(path: str) -> Path:
    if path.startswith("file://"):
        parsed = urlparse(path)
        return Path(unquote(parsed.path))
    return Path(path)


def _parse_s3_url(url: str) -> tuple[str, str]:
    if not url.startswith("s3://"):
        raise ValueError(f"Expected an s3:// URL, got {url}")
    try:
        bucket, key = url.replace("s3://", "", 1).split("/", 1)
    except ValueError as exc:
        raise ValueError(f"S3 URL must include a key, got {url}") from exc
    return bucket, key


def _build_s3_filesystem(storage_options: dict[str, Any]) -> AsyncFileSystem:
    try:
        import s3fs
    except ImportError as exc:
        raise ImportError(
            "s3fs is required to read s3:// LUT zip paths. Install with `pip install s3fs`."
        ) from exc

    options = dict(storage_options)
    key = options.pop("key", None)
    secret = options.pop("secret", None)
    anon = options.pop("anon", None)
    region = options.pop("region", None)
    endpoint_url = options.pop("endpoint_url", None)
    client_kwargs = dict(options.pop("client_kwargs", {}))

    if region is not None and "region_name" not in client_kwargs:
        client_kwargs["region_name"] = region
    if endpoint_url is not None and "endpoint_url" not in client_kwargs:
        client_kwargs["endpoint_url"] = endpoint_url

    s3_kwargs: dict[str, Any] = {
        "asynchronous": True,
        "client_kwargs": client_kwargs,
        **options,
    }
    if anon is not None:
        s3_kwargs["anon"] = bool(anon)

    # Keep ambient AWS credential resolution as the default path.
    if not bool(anon):
        if key is not None:
            s3_kwargs["key"] = key
        if secret is not None:
            s3_kwargs["secret"] = secret

    return s3fs.S3FileSystem(**s3_kwargs)


def _detect_zarr_prefix_from_names(names: Iterable[str]) -> str:
    """Detect optional nested root containing top-level Zarr metadata markers."""
    marker_set = {".zgroup", ".zattrs", ".zmetadata", "zarr.json"}
    names_list = list(names)

    if any(name in marker_set for name in names_list):
        return ""

    candidates: list[str] = []
    for name in names_list:
        for marker in marker_set:
            suffix = f"/{marker}"
            if name.endswith(suffix):
                candidates.append(name[: -len(suffix)])
                break
    if not candidates:
        return ""
    return min(candidates, key=len)


def _detect_zarr_prefix(zip_fs: _ReadOnlyZipFileSystem) -> str:
    """Detect optional nested root that contains top-level Zarr metadata files."""
    sync(zip_fs.loop, zip_fs._initialize)
    assert zip_fs._files is not None
    return _detect_zarr_prefix_from_names(zip_fs._files.keys())


def build_readonly_zip_mapper(path: str, storage_options: dict[str, Any]) -> MutableMapping[str, bytes]:
    """Build an FSMap for a ZIP-hosted Zarr store at local/HTTP/S3 path."""
    options = dict(storage_options)

    if path.startswith(("http://", "https://")):
        timeout = float(options.pop("timeout", 30.0))
        headers = options.pop("headers", None)
        if headers is not None and not isinstance(headers, dict):
            raise TypeError("storage_options['headers'] must be a dictionary if provided")
        if options:
            # Unsupported HTTP options should be explicit to avoid silent misuse.
            unknown = ", ".join(sorted(options))
            raise TypeError(f"Unsupported HTTP zip storage option(s): {unknown}")
        http_fs_kwargs: dict[str, Any] = {"timeout": timeout}
        if headers is not None:
            http_fs_kwargs["headers"] = headers
        base_fs = _HTTPRangeFileSystem(**http_fs_kwargs)
        zip_path = path
    elif path.startswith("s3://"):
        bucket, key = _parse_s3_url(path)
        base_fs = _build_s3_filesystem(options)
        zip_path = f"{bucket}/{key}"
    else:
        if options:
            unknown = ", ".join(sorted(options))
            raise TypeError(f"Unsupported local zip storage option(s): {unknown}")
        base_fs = _build_local_filesystem()
        zip_path = str(_parse_local_path(path))

    zip_fs = _ReadOnlyZipFileSystem(base_fs, zip_path)
    try:
        root = _detect_zarr_prefix(zip_fs)
        return FSMap(root=root, fs=zip_fs, check=False, create=False)
    except ValueError as exc:
        message = str(exc).lower()
        if "not stored (uncompressed)" not in message:
            raise
        if path.startswith(("http://", "https://", "s3://")):
            raise ValueError(
                f"Compressed remote ZIP archives are not supported by ReadOnlyZipFileSystem: {path}"
            ) from exc

        # Local compressed ZIP fallback: delegate to fsspec's standard zip mapper.
        import fsspec

        fallback_zip_fs = fsspec.filesystem("zip", fo=zip_path)
        fallback_mapper = fallback_zip_fs.get_mapper("")
        root = _detect_zarr_prefix_from_names(fallback_mapper.keys())
        return fallback_zip_fs.get_mapper(root)
