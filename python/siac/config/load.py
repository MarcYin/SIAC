"""
TOML loading and environment overlay utilities for SIAC system config.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING, cast

try:
    import tomllib  # type: ignore[import-not-found]
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib

try:
    import tomli_w
except ModuleNotFoundError:  # pragma: no cover - optional write dependency

    class _TomliWFallback:
        @staticmethod
        def dump(_payload: object, _handle: object) -> None:
            raise ModuleNotFoundError("tomli_w is required for writing SIAC TOML config files")

    tomli_w = _TomliWFallback()

from siac.config.schema import SystemConfig

if TYPE_CHECKING:
    from collections.abc import Mapping

DEFAULT_CONFIG_PATH = Path.home() / ".config" / "siac" / "config.toml"
CONFIG_PATH_ENV = "SIAC_CONFIG_FILE"


def load_system_config(path: Path | str) -> SystemConfig:
    """Load a system config from TOML without applying environment overlays."""
    resolved = Path(path).expanduser()
    with resolved.open("rb") as handle:
        loaded = tomllib.load(handle) or {}
    if not isinstance(loaded, dict):
        raise ValueError("System config TOML must contain a table/object at the top level.")
    return cast("SystemConfig", SystemConfig.model_validate(loaded))


def load_system_config_from_default() -> SystemConfig:
    """Load from SIAC_CONFIG_FILE or the default user config path."""
    raw_path = os.getenv(CONFIG_PATH_ENV)
    path = Path(raw_path) if raw_path else DEFAULT_CONFIG_PATH
    return load_system_config(path)


def overlay_env_secrets(
    config: SystemConfig,
    env: Mapping[str, str] | None = None,
) -> SystemConfig:
    """Fill missing auth values from environment variables."""
    resolved_env = dict(os.environ if env is None else env)
    updated = cast("SystemConfig", config.model_copy(deep=True))

    if updated.auth.cdse.username is None:
        updated.auth.cdse.username = resolved_env.get(updated.auth.cdse.username_env)
    if updated.auth.cdse.password is None:
        updated.auth.cdse.password = resolved_env.get(updated.auth.cdse.password_env)

    if updated.auth.cds.api_key is None:
        updated.auth.cds.api_key = resolved_env.get(updated.auth.cds.api_key_env)

    if updated.auth.aws.access_key_id is None:
        updated.auth.aws.access_key_id = resolved_env.get(updated.auth.aws.access_key_id_env)
    if updated.auth.aws.secret_access_key is None:
        updated.auth.aws.secret_access_key = resolved_env.get(
            updated.auth.aws.secret_access_key_env
        )

    if updated.auth.earthdata.username is None:
        updated.auth.earthdata.username = resolved_env.get(updated.auth.earthdata.username_env)
    if updated.auth.earthdata.password is None:
        updated.auth.earthdata.password = resolved_env.get(updated.auth.earthdata.password_env)

    if updated.auth.gcs.credentials_file is None:
        raw_gcs = resolved_env.get(updated.auth.gcs.credentials_file_env)
        if raw_gcs:
            updated.auth.gcs.credentials_file = Path(raw_gcs).expanduser()

    return updated


def write_system_config(config: SystemConfig, path: Path | str) -> None:
    """Write a system config to TOML."""
    resolved = Path(path).expanduser()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    payload = config.model_dump(mode="json", exclude_none=True)
    with resolved.open("wb") as handle:
        tomli_w.dump(payload, handle)


def write_default_system_config(path: Path | str = DEFAULT_CONFIG_PATH) -> Path:
    """Write a default TOML config file."""
    resolved = Path(path).expanduser()
    if not resolved.suffix:
        resolved = resolved.with_suffix(".toml")
    write_system_config(SystemConfig(), resolved)
    return resolved
