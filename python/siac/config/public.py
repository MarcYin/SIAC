"""
Public config facade built on the categorized SIAC config package.
"""

from __future__ import annotations

from pathlib import Path  # noqa: TC003
from typing import Any, cast

from pydantic import Field

from siac.config.load import (
    DEFAULT_CONFIG_PATH,
    load_system_config,
    load_system_config_from_default,
    overlay_env_secrets,
    write_default_system_config,
    write_system_config,
)
from siac.config.resolve import resolve_config
from siac.config.schema import (
    DEFAULT_LUT_URL,
    AlgorithmsConfig,
    AtmoProviderConfig,
    AuthConfig,
    BRDFProviderConfig,
    CloudMaskAlgorithmConfig,
    ExecutionRuntimeConfig,
    MonthlyCompositeProviderConfig,
    OutputDefaultsConfig,
    PathsConfig,
    ProvidersConfig,
    RTAlgorithmConfig,
    RunRequest,
    RuntimeConfig,
    S2ProviderConfig,
    SensorName,
    SolverAlgorithmConfig,
    SurfacePriorAlgorithmConfig,
    SystemConfig,
)
from siac.config.snapshot import snapshot_system_config, write_runtime_snapshot


def _deep_merge(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


class SIACConfig(SystemConfig):
    """
    Public SIAC config wrapper.

    The canonical TOML schema excludes per-run inputs, but this wrapper keeps
    lightweight Python-only defaults for `sensor` and `aoi` so call sites can
    still set those without encoding them into the file format.
    """

    sensor: SensorName = "auto"
    aoi: dict[str, Any] | Path | str | tuple[float, float, float, float] | list[float] | None = (
        Field(default=None)
    )

    @classmethod
    def from_file(cls, path: Path | str) -> SIACConfig:
        system = load_system_config(path)
        payload = system.model_dump(mode="python")
        return cast("SIACConfig", cls.model_validate(payload))

    @classmethod
    def from_toml(cls, path: Path | str) -> SIACConfig:
        return cls.from_file(path)

    @classmethod
    def from_yaml(cls, _path: Path | str) -> SIACConfig:
        raise ValueError("YAML config files are no longer supported; use TOML.")

    @classmethod
    def load(cls, path: Path | str | None = None, **overrides: Any) -> SIACConfig:
        config = (
            cls.from_file(path)
            if path is not None
            else cast(
                "SIACConfig",
                cls.model_validate(load_system_config_from_default().model_dump(mode="python")),
            )
        )
        if not overrides:
            return config
        return config.with_overrides(**overrides)

    def to_file(self, path: Path | str) -> None:
        write_system_config(self._system_only(), path)

    def to_toml(self, path: Path | str) -> None:
        self.to_file(path)

    def to_yaml(self, _path: Path | str) -> None:
        raise ValueError("YAML config files are no longer supported; use TOML.")

    def to_dict(self) -> dict[str, Any]:
        return cast("dict[str, Any]", self._system_only().model_dump(mode="python"))

    def with_overrides(self, **kwargs: Any) -> SIACConfig:
        payload = self.model_dump(mode="python")
        return cast("SIACConfig", type(self).model_validate(_deep_merge(payload, kwargs)))

    def with_env_overlay(self) -> SIACConfig:
        system = overlay_env_secrets(self._system_only())
        payload = system.model_dump(mode="python")
        payload["sensor"] = self.sensor
        payload["aoi"] = self.aoi
        return cast("SIACConfig", type(self).model_validate(payload))

    def resolve(self, request: RunRequest) -> Any:
        return resolve_config(self.with_env_overlay(), request)

    def snapshot(self, *, redact_secrets: bool = True) -> dict[str, Any]:
        return snapshot_system_config(self, redact_secrets=redact_secrets)

    def write_state_snapshot(self, path: Path | str, *, redact_secrets: bool = True) -> None:
        resolved = resolve_config(
            self.with_env_overlay(), RunRequest(sensor=self.sensor, aoi=self.aoi)
        )
        write_runtime_snapshot(resolved, path, redact_secrets=redact_secrets)

    @classmethod
    def write_default_config(cls, path: Path | str = DEFAULT_CONFIG_PATH) -> Path:
        return write_default_system_config(path)

    def _system_only(self) -> SystemConfig:
        return cast(
            "SystemConfig",
            SystemConfig.model_validate(
                self.model_dump(
                    mode="python",
                    exclude={"sensor", "aoi"},
                )
            ),
        )


AtmoPriorConfig = AtmoProviderConfig
BRDFConfig = BRDFProviderConfig
S2DataAccessConfig = S2ProviderConfig
MonthlyCompositeConfig = MonthlyCompositeProviderConfig
SurfacePriorConfig = SurfacePriorAlgorithmConfig
RTModelConfig = RTAlgorithmConfig
SolverConfig = SolverAlgorithmConfig
CloudMaskConfig = CloudMaskAlgorithmConfig
ExecutionConfig = ExecutionRuntimeConfig
CredentialConfig = AuthConfig
OutputConfig = OutputDefaultsConfig


def get_default_config() -> SIACConfig:
    return SIACConfig()


def get_jasmin_config() -> SIACConfig:
    return SIACConfig(
        providers=ProvidersConfig(
            atmo=AtmoProviderConfig(
                kind="cams",
                data_path="/work/scratch-pw3/marc/CAMS/",
            ),
            brdf=BRDFProviderConfig(
                kind="mcd43",
                data_path="/gws/nopw/j04/nceo_ard/public/MCD43/",
            ),
        ),
        paths=PathsConfig(
            dem="/vsicurl/https://raw.githubusercontent.com/MarcYin/Copernicus_GLO_30_DEM_VRT/refs/heads/main/copernicus_GLO_30_dem.vrt",
        ),
        runtime=RuntimeConfig(n_jobs=8),
    )


def get_lut_config(lut_path: str | Path) -> SIACConfig:
    return SIACConfig(
        algorithms=AlgorithmsConfig(rt=RTAlgorithmConfig(backend="lut")),
        paths=PathsConfig(lut_path=lut_path),
    )


__all__ = [
    "DEFAULT_LUT_URL",
    "DEFAULT_CONFIG_PATH",
    "SIACConfig",
    "SystemConfig",
    "RunRequest",
    "PathsConfig",
    "AuthConfig",
    "AtmoPriorConfig",
    "BRDFConfig",
    "S2DataAccessConfig",
    "MonthlyCompositeConfig",
    "RTModelConfig",
    "SolverConfig",
    "CloudMaskConfig",
    "ExecutionConfig",
    "CredentialConfig",
    "OutputConfig",
    "get_default_config",
    "get_jasmin_config",
    "get_lut_config",
]
