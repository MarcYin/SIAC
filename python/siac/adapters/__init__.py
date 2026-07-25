"""Adapter-layer helpers for external systems and backends."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "CAMSProvider": ("siac.adapters.atmo", "CAMSProvider"),
    "MERRA2Provider": ("siac.adapters.atmo", "MERRA2Provider"),
    "MCD19AODProvider": ("siac.adapters.atmo", "MCD19AODProvider"),
    "CachedMCD19AODProvider": ("siac.adapters.atmo", "CachedMCD19AODProvider"),
    "VNP19AODProvider": ("siac.adapters.atmo", "VNP19AODProvider"),
    "AWSAuth": ("siac.adapters.auth", "AWSAuth"),
    "CDSAuth": ("siac.adapters.auth", "CDSAuth"),
    "CDSEAuth": ("siac.adapters.auth", "CDSEAuth"),
    "CredentialManager": ("siac.adapters.auth", "CredentialManager"),
    "CredentialSpec": ("siac.adapters.auth", "CredentialSpec"),
    "EarthdataAuth": ("siac.adapters.auth", "EarthdataAuth"),
    "GCSAuth": ("siac.adapters.auth", "GCSAuth"),
    "MCD43EarthAccessProvider": ("siac.adapters.brdf", "MCD43EarthAccessProvider"),
    "VNP43EarthAccessProvider": ("siac.adapters.brdf", "VNP43EarthAccessProvider"),
    "MCD19EarthAccessProvider": ("siac.adapters.brdf", "MCD19EarthAccessProvider"),
    "MODLAND_SINUSOIDAL_CRS": ("siac.adapters.earthdata_common", "MODLAND_SINUSOIDAL_CRS"),
    "modland_tile_coords": ("siac.adapters.earthdata_common", "modland_tile_coords"),
    "ConfiguredOutputWriter": ("siac.adapters.output", "ConfiguredOutputWriter"),
    "load_band_rsrf": ("siac.adapters.rsrf", "load_band_rsrf"),
    "load_sensor_config_with_rsrf": ("siac.adapters.rsrf", "load_sensor_config_with_rsrf"),
    "BaseSatellitePreprocessor": ("siac.adapters.satellite", "BaseSatellitePreprocessor"),
    "Sentinel2Preprocessor": ("siac.adapters.satellite", "Sentinel2Preprocessor"),
    "detect_sensor": ("siac.adapters.satellite", "detect_sensor"),
    "get_preprocessor": ("siac.adapters.satellite", "get_preprocessor"),
    "earthaccess_source_from_auth": ("siac.adapters.earthdata", "earthaccess_source_from_auth"),
    "build_rt_model": ("siac.adapters.rt", "build_rt_model"),
    "build_s2_backend": ("siac.adapters.s2_backend", "build_s2_backend"),
}
_SUBMODULES = frozenset(
    {
        "atmo",
        "auth",
        "brdf",
        "data",
        "earthdata",
        "earthdata_common",
        "output",
        "rsrf",
        "rt",
        "s2_backend",
        "satellite",
    }
)


def __getattr__(name: str) -> Any:
    if name in _SUBMODULES:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module

    try:
        module_name, attr_name = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc

    value = getattr(import_module(module_name), attr_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(list(globals()) + list(_EXPORTS) + list(_SUBMODULES))
