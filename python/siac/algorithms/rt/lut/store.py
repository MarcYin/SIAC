"""Store loading helpers for local/remote LUT Zarr backends."""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast
from urllib.parse import unquote, urlparse

from siac.algorithms.rt.lut.http_zip_store import build_readonly_zip_mapper

if TYPE_CHECKING:
    from collections.abc import MutableMapping

logger = logging.getLogger(__name__)

_ZARR_MARKER_KEYS = (".zgroup", ".zattrs", ".zmetadata", "zarr.json")
_DEFAULT_REFERENCE_CACHE_DIR = Path.home() / ".cache" / "siac" / "lut_refs"


@dataclass(frozen=True)
class _ReferenceOptions:
    refresh: bool
    reference_json: Path | None
    cache_dir: Path | None


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


def _coerce_bool(value: Any, *, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "yes", "on"}:
            return True
        if normalized in {"0", "false", "no", "off"}:
            return False
    return bool(value)


def _split_storage_options(
    storage_options: dict[str, Any],
) -> tuple[dict[str, Any], _ReferenceOptions]:
    """Split reference-cache options from backend storage options."""
    options = dict(storage_options)
    if "use_reference" in options:
        raise TypeError(
            "storage_options['use_reference'] is no longer supported; "
            "remote zipped LUTs always use reference JSON mapping."
        )
    reference_json = options.pop("reference_json", None)
    reference_cache_dir = options.pop("reference_cache_dir", None)
    reference_options = _ReferenceOptions(
        refresh=_coerce_bool(options.pop("reference_refresh", False), default=False),
        reference_json=Path(str(reference_json)).expanduser()
        if reference_json is not None
        else None,
        cache_dir=Path(str(reference_cache_dir)).expanduser()
        if reference_cache_dir is not None
        else None,
    )
    return options, reference_options


def normalize_storage_options(path: str, storage_options: dict[str, Any]) -> dict[str, Any]:
    """Validate and return storage options for the selected backend."""
    options = dict(storage_options)

    if not path.startswith("s3://"):
        return options

    unsupported = {"region", "endpoint_url"} & options.keys()
    if unsupported:
        names = ", ".join(sorted(unsupported))
        raise TypeError(
            f"Top-level S3 storage option(s) are not supported: {names}. "
            "Use storage_options['client_kwargs'] instead."
        )

    return options


def _normalize_mapper_root(root: str) -> str:
    normalized = (root or "").strip("/")
    return "" if normalized in ("", ".") else normalized


def _relative_reference_key(name: str, root: str) -> str | None:
    normalized_name = name.strip("/")
    if normalized_name == "":
        return None
    if not root:
        return normalized_name
    if normalized_name == root:
        return None
    prefix = f"{root}/"
    if not normalized_name.startswith(prefix):
        return None
    return normalized_name[len(prefix) :]


def _build_reference_document(path: str, zip_mapper: Any) -> dict[str, Any]:
    """Build a reference-spec document from a zip-backed FSMap."""
    zip_fs = getattr(zip_mapper, "fs", None)
    files = getattr(zip_fs, "_files", None)
    if not isinstance(files, dict):
        raise ValueError("ZIP mapper does not expose indexed file metadata")

    root = _normalize_mapper_root(str(getattr(zip_mapper, "root", "")))
    refs: dict[str, list[Any]] = {}

    for name, info in files.items():
        if not isinstance(info, dict) or "children" in info:
            continue
        key = _relative_reference_key(str(name), root)
        if key is None:
            continue
        offset = int(info.get("offset", 0))
        size = int(info.get("size", 0))
        if offset < 0 or size < 0:
            raise ValueError(f"Invalid byte range for {name}: offset={offset}, size={size}")
        refs[key] = [path, offset, size]

    if not refs:
        raise ValueError("No files available to build reference mapping")

    return {"version": 1, "refs": refs}


def _is_valid_reference_document(document: dict[str, Any]) -> bool:
    refs = document.get("refs")
    if not isinstance(refs, dict) or not refs:
        return False
    return any(marker in refs for marker in _ZARR_MARKER_KEYS)


def _reference_json_path(path: str, reference_options: _ReferenceOptions) -> Path:
    if reference_options.reference_json is not None:
        return reference_options.reference_json

    cache_dir = reference_options.cache_dir or _DEFAULT_REFERENCE_CACHE_DIR
    parsed = urlparse(path)
    stem = Path(parsed.path).name or "lut.zarr.zip"
    digest = hashlib.sha256(path.encode("utf-8")).hexdigest()[:16]
    return cache_dir / f"{stem}.{digest}.references.json"


def _write_reference_document(reference_path: Path, document: dict[str, Any]) -> None:
    reference_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = reference_path.with_suffix(reference_path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(document, separators=(",", ":")), encoding="utf-8")
    tmp_path.replace(reference_path)


def _load_reference_document(reference_path: Path) -> dict[str, Any]:
    return cast("dict[str, Any]", json.loads(reference_path.read_text(encoding="utf-8")))


def _reference_remote(
    path: str, storage_options: dict[str, Any]
) -> tuple[str | None, dict[str, Any]]:
    """Resolve remote protocol/options for reference:// mapper."""
    scheme = urlparse(path).scheme or None
    if scheme in {"http", "https"}:
        headers = storage_options.get("headers")
        if headers is not None and not isinstance(headers, dict):
            raise TypeError("storage_options['headers'] must be a dictionary if provided")
        timeout = float(storage_options.get("timeout", 30.0))
        options: dict[str, Any] = {"timeout": timeout, "asynchronous": True}
        if headers is not None:
            options["headers"] = headers
        return scheme, options
    if scheme is None:
        return None, {}
    options = dict(storage_options)
    options.setdefault("asynchronous", True)
    return scheme, options


def _open_reference_mapper(
    path: str,
    storage_options: dict[str, Any],
    reference_path: Path,
    document: dict[str, Any] | None = None,
) -> MutableMapping[str, bytes]:
    import fsspec  # type: ignore[import-untyped]

    doc = document or _load_reference_document(reference_path)
    refs = doc.get("refs")
    if not isinstance(refs, dict):
        raise ValueError(f"Invalid reference JSON at {reference_path}")

    remote_protocol, remote_options = _reference_remote(path, storage_options)
    kwargs: dict[str, Any] = {"fo": doc, "asynchronous": True}
    if remote_protocol is not None:
        kwargs["remote_protocol"] = remote_protocol
    if remote_options:
        kwargs["remote_options"] = remote_options

    return cast("MutableMapping[str, bytes]", fsspec.get_mapper("reference://", **kwargs))


def _build_remote_zip_reference_mapper(
    path: str,
    storage_options: dict[str, Any],
    reference_options: _ReferenceOptions,
) -> MutableMapping[str, bytes]:
    """Build/open a cached reference mapper for remote ZIP LUTs."""
    reference_path = _reference_json_path(path, reference_options)

    if reference_path.exists() and not reference_options.refresh:
        try:
            document = _load_reference_document(reference_path)
            if _is_valid_reference_document(document):
                return _open_reference_mapper(
                    path,
                    storage_options,
                    reference_path,
                    document=document,
                )
            logger.warning(
                "Invalid cached LUT reference JSON at %s; rebuilding.",
                reference_path,
            )
        except Exception as exc:
            logger.warning(
                "Failed to read LUT reference JSON at %s (%s); rebuilding.",
                reference_path,
                exc,
            )

    zip_mapper = build_readonly_zip_mapper(path, storage_options)
    document = _build_reference_document(path, zip_mapper)
    _write_reference_document(reference_path, document)
    return _open_reference_mapper(
        path,
        storage_options,
        reference_path,
        document=document,
    )


def build_lut_store(path: str, storage_options: dict[str, Any]) -> str | MutableMapping[str, bytes]:
    """Build a zarr-readable store from local/HTTP/S3 path or zip archive."""
    import fsspec

    backend_options, reference_options = _split_storage_options(storage_options)
    resolved_storage_options = normalize_storage_options(path, backend_options)
    local_path = as_local_path(path)
    path_or_url = str(local_path) if local_path is not None else path

    if is_zip_store(path_or_url):
        if is_remote_path(path_or_url):
            return _build_remote_zip_reference_mapper(
                path_or_url,
                resolved_storage_options,
                reference_options,
            )
        return build_readonly_zip_mapper(path_or_url, resolved_storage_options)

    if is_remote_path(path):
        return cast(
            "MutableMapping[str, bytes]", fsspec.get_mapper(path, **resolved_storage_options)
        )

    if resolved_storage_options:
        return cast(
            "MutableMapping[str, bytes]",
            fsspec.get_mapper(path_or_url, **resolved_storage_options),
        )

    return path_or_url
