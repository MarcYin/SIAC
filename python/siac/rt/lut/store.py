"""Store loading helpers for local/remote LUT Zarr backends."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

from siac.rt.lut.http_zip_store import build_readonly_zip_mapper


def is_remote_path(path: str) -> bool:
    """Return True for non-local LUT URIs."""
    return path.startswith(("http://", "https://", "s3://"))


def is_zip_store(path: str) -> bool:
    """Return True when path points to a zip archive."""
    return path.lower().endswith(".zip")


def as_local_path(path: str) -> Path | None:
    """Convert local/file URI path to Path; return None for remote URIs."""
    if is_remote_path(path):
        return None

    if path.startswith("file://"):
        parsed = urlparse(path)
        return Path(unquote(parsed.path))

    return Path(path)


def normalize_storage_options(path: str, storage_options: dict[str, Any]) -> dict[str, Any]:
    """Normalize generic storage options, including S3 region/endpoint aliases."""
    options = dict(storage_options)

    if not path.startswith("s3://"):
        return options

    region = options.pop("region", None)
    endpoint_url = options.pop("endpoint_url", None)
    client_kwargs = dict(options.get("client_kwargs", {}))

    if region is not None and "region_name" not in client_kwargs:
        client_kwargs["region_name"] = region
    if endpoint_url is not None and "endpoint_url" not in client_kwargs:
        client_kwargs["endpoint_url"] = endpoint_url

    if client_kwargs:
        options["client_kwargs"] = client_kwargs

    return options


def build_lut_store(path: str, storage_options: dict[str, Any]) -> Any:
    """Build a zarr-readable store from local/HTTP/S3 path or zip archive."""
    import fsspec

    resolved_storage_options = normalize_storage_options(path, storage_options)
    local_path = as_local_path(path)
    path_or_url = str(local_path) if local_path is not None else path

    if is_zip_store(path_or_url):
        return build_readonly_zip_mapper(path_or_url, resolved_storage_options)

    if is_remote_path(path):
        return fsspec.get_mapper(path, **resolved_storage_options)

    if resolved_storage_options:
        return fsspec.get_mapper(path_or_url, **resolved_storage_options)

    return path_or_url
