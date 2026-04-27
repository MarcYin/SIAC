"""
Diagnostic snapshots for SIAC config objects.
"""

from __future__ import annotations

from datetime import datetime, timezone
from importlib import import_module
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast
from urllib.parse import urlparse

yaml = cast("Any", import_module("yaml"))

if TYPE_CHECKING:
    from siac.config.schema import ResolvedConfig, SystemConfig

_REMOTE_URI_SCHEMES = {"http", "https", "s3", "file", "gs"}
_REDACTION_TOKEN = "<redacted>"


def _normalize(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return [_normalize(item) for item in value]
    if isinstance(value, list):
        return [_normalize(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _normalize(item) for key, item in value.items()}
    return value


def _path_state(value: Any) -> dict[str, Any]:
    normalized = _normalize(value)
    if normalized is None:
        return {"value": None, "kind": "unset", "exists": None}
    if isinstance(normalized, str):
        if (
            normalized.startswith("/vsi")
            or urlparse(normalized).scheme.lower() in _REMOTE_URI_SCHEMES
        ):
            return {"value": normalized, "kind": "remote", "exists": None}
        path = Path(normalized).expanduser().resolve(strict=False)
    else:
        path = Path(normalized).expanduser().resolve(strict=False)
    if path.exists():
        return {
            "value": str(path),
            "kind": "directory" if path.is_dir() else "file",
            "exists": True,
        }
    return {"value": str(path), "kind": "missing", "exists": False}


def _redact_auth(payload: dict[str, Any]) -> dict[str, Any]:
    redacted = {str(key): _normalize(item) for key, item in payload.items()}
    for provider, secret_keys in (
        ("cdse", ("username", "password")),
        ("cds", ("api_key",)),
        ("aws", ("access_key_id", "secret_access_key")),
        ("earthdata", ("username", "password")),
        ("gcs", ("credentials_file",)),
    ):
        section = redacted.get(provider) or {}
        for key in secret_keys:
            if section.get(key) is not None:
                section[key] = _REDACTION_TOKEN
    return redacted


def snapshot_system_config(
    config: SystemConfig,
    *,
    redact_secrets: bool = True,
    source_path: Path | str | None = None,
) -> dict[str, Any]:
    resolved_source = None if source_path is None else str(Path(source_path).expanduser())
    payload: dict[str, Any] = config.model_dump(mode="json")
    if redact_secrets:
        payload["auth"] = _redact_auth(payload["auth"])
    return {
        "source_path": resolved_source,
        "config": payload,
    }


def snapshot_resolved_config(
    config: ResolvedConfig,
    *,
    redact_secrets: bool = True,
    source_path: Path | str | None = None,
) -> dict[str, Any]:
    resolved_source = None if source_path is None else str(Path(source_path).expanduser())
    payload: dict[str, Any] = config.model_dump(mode="json")
    if redact_secrets:
        payload["auth"] = _redact_auth(payload["auth"])
    return {
        "source_path": resolved_source,
        "config": payload,
        "path_states": {
            "paths.dem": _path_state(config.paths.dem),
            "paths.water_mask": _path_state(config.paths.water_mask),
            "paths.emulator_dir": _path_state(config.paths.emulator_dir),
            "paths.lut_path": _path_state(config.paths.lut_path),
            "paths.rsrf_root": _path_state(config.paths.rsrf_root),
            "paths.caches.atmo": _path_state(config.paths.caches.atmo),
            "paths.caches.brdf": _path_state(config.paths.caches.brdf),
            "paths.caches.s2": _path_state(config.paths.caches.s2),
            "output.defaults.output_dir": _path_state(config.output.defaults.output_dir),
        },
    }


def write_runtime_snapshot(
    config: ResolvedConfig,
    path: Path | str,
    *,
    redact_secrets: bool = True,
    source_path: Path | str | None = None,
) -> None:
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "state": snapshot_resolved_config(
            config,
            redact_secrets=redact_secrets,
            source_path=source_path,
        ),
    }
    resolved = Path(path).expanduser()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    with resolved.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, default_flow_style=False, sort_keys=False)
