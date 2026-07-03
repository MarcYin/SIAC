"""Public SIAC configuration helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

from pydantic import Field

from siac.config.algorithms import AlgorithmsConfig, RTAlgorithmConfig
from siac.config.load import (
    DEFAULT_CONFIG_PATH,
    load_system_config,
    load_system_config_from_default,
    overlay_env_secrets,
    write_default_system_config,
    write_system_config,
)
from siac.config.providers import (
    AtmoProviderConfig,
    BRDFProviderConfig,
    PathsConfig,
    ProvidersConfig,
)
from siac.config.request import RunRequest
from siac.config.resolve import resolve_config
from siac.config.snapshot import snapshot_system_config, write_runtime_snapshot
from siac.config.system import RuntimeConfig, SystemConfig
from siac.config.types import SensorName


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

    sensor: SensorName = SensorName.AUTO
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


def _get_surface_driven_bestpixel_config(
    *,
    cache_root: str | Path | None = None,
    bestpixel_endpoint: str,
    bestpixel_source: str,
) -> SIACConfig:
    """Return the opt-in clean-day bestpixel surface-driven AOD recipe.

    This is the production-facing form of the clean-surface prior setup:
    bestpixel composites, MAIAC low-AOD day gating, and the surface-driven AOT
    solver. It intentionally does not change package defaults; callers must opt
    into this config.
    """

    root = Path(cache_root).expanduser() if cache_root is not None else None

    payload: dict[str, Any] = {
        "paths": {
            "dem": "none",
            "water_mask": (
                "https://zenodo.org/records/14899246/files/landWater2020.vrt?download=1"
            ),
        },
        "providers": {
            "atmo": {
                "kind": "mcd19",
                "maiac_best_quality_qa": True,
            },
            "brdf": {
                "kind": "mcd43",
                "temporal_window": 16,
            },
            "s2": {
                "processing_level": "L1C",
                "max_cloud_cover": 100.0,
            },
            "monthly_composites": {
                "kind": "bestpixel",
                "bestpixel_endpoint": bestpixel_endpoint,
                "bestpixel_lookback_years": 5,
                "bestpixel_seasonal_window_months": 1,
                "bestpixel_top_k": 3,
                "bestpixel_max_cloud_cover": 90.0,
            },
        },
        "algorithms": {
            "cloud_mask": {
                "mode": "auto",
                "provider": "omnicloudmask",
                "target_resolution_m": 10.0,
            },
            "surface_prior": {
                "method": "bestpixel",
                "bestpixel_source": bestpixel_source,
                "bestpixel_window_reduction": "daily_median",
                "bestpixel_low_aod_frac": 0.6,
                "bestpixel_robust_clip": 1.5,
                "spectral_mapping": {
                    "enabled": True,
                    "k_neighbors": 5,
                    "knn_backend": "scipy_ckdtree",
                },
            },
            "rt": {
                "backend": "lut",
                "fallback_to_lut": True,
            },
            "solver": {
                "method": "surface_driven",
                "aerosol_resolution": 60.0,
                "grid_search_aot_points": 11,
                "surface_driven_backstop_calibrated": True,
                "surface_driven_reference_tcwv": 2.0,
                "surface_driven_scene_mean_geometry": True,
                "surface_driven_resolve_on_prior_grid": False,
                "surface_driven_pool_radius_m": 600.0,
                "water_mask_buffer_pixels": 32,
                "surface_driven_solve_bands": ["B01", "B02", "B04"],
                "surface_driven_cost_mode": "auto2",
                "surface_driven_aot_axis": "acixthree",
                "surface_driven_ignore_cloud_water": True,
                "surface_driven_aod_clean": 0.15,
                "surface_driven_aod_high": 0.6,
                "surface_driven_aerosol_species": "none",
                "surface_driven_aerosol_species_candidates": 3,
            },
        },
    }

    if root is not None:
        payload["paths"]["cache_root"] = root
        payload["providers"]["atmo"]["cache_dir"] = root / "atmo"
        payload["providers"]["brdf"]["cache_dir"] = root / "brdf"
        payload["providers"]["s2"]["cache_dir"] = root / "s2"
        payload["providers"]["monthly_composites"]["bestpixel_disk_cache"] = root / "bestpixel"

    return SIACConfig.model_validate(payload)


def get_surface_driven_l2a_config(*, cache_root: str | Path | None = None) -> SIACConfig:
    """Return the opt-in S2 L2A clean-day surface-driven AOD recipe.

    Uses the Planetary Computer Sentinel-2 L2A bestpixel endpoint
    (``bestpixel_endpoint="pc"``), MAIAC low-AOD day gating, and the
    surface-driven AOT solver.
    """

    return _get_surface_driven_bestpixel_config(
        cache_root=cache_root,
        bestpixel_endpoint="pc",
        bestpixel_source="l2a",
    )


#: Global visible-predictor debias for Sen2Cor L2A dictionaries solved in the
#: packaged libRadtran LUT space: ``prediction += intercept + slope * anchor_aot``.
#: Calibrated 2026-07-03 by Theil-Sen over 61 AERONET-site scenes using ONLY
#: same-day QA-MAIAC AOD as the reference (leave-one-site-out validated:
#: 52/62 = 83.9% within-EE, AOD RMSE 0.108, deterministic). Pair-specific:
#: valid for Sen2Cor-composite priors + LUT solve only.
L2A_LUT_PREDICT_VISIBLE_DEBIAS: dict[str, tuple[float, float]] = {
    "B02": (-0.0003, 0.0243),
    "B04": (-0.0011, 0.0223),
}


def get_surface_driven_l2a_monthly_predictor_config(
    *, cache_root: str | Path | None = None
) -> SIACConfig:
    """Return the default single-surface monthly L2A predictor AOD recipe (K2).

    This preset uses Planetary Computer Sentinel-2 L2A monthly bestpixel
    composites as the dictionary for the ExtraTree visible-band surface prior
    (anchor bands corrected through the configured RT backend), applies the
    globally MAIAC-calibrated Sen2Cor->LUT debias and a 0.015 BOA uncertainty
    floor, and solves surface-driven AOD from B02/B04 with the standard 1.2 km
    cost-pooling window. Validated 2026-07-03: 52/62 = 83.9% within-EE,
    AOD RMSE 0.108 over the 62-site AERONET set (deterministic).
    """

    return _get_surface_driven_bestpixel_config(
        cache_root=cache_root,
        bestpixel_endpoint="pc",
        bestpixel_source="l2a",
    ).with_overrides(
        providers={
            "monthly_composites": {
                "bestpixel_top_k": 15,
            },
        },
        algorithms={
            "surface_prior": {
                "bestpixel_window_reduction": "window",
                "bestpixel_predict_visible": True,
                "bestpixel_predict_visible_bands": ("B02", "B04"),
                "bestpixel_predict_visible_uncertainty_floor": 0.015,
                "bestpixel_predict_visible_debias": dict(L2A_LUT_PREDICT_VISIBLE_DEBIAS),
            },
            "solver": {
                "surface_driven_solve_bands": ("B02", "B04"),
            },
        },
    )


def get_surface_driven_hls_s30_config(*, cache_root: str | Path | None = None) -> SIACConfig:
    """Return the opt-in HLS S30 clean-day surface-driven AOD recipe.

    Uses the Planetary Computer HLS S30-only bestpixel endpoint
    (``bestpixel_endpoint="hls-s30"``), MAIAC low-AOD day gating, and the same
    surface-driven AOT solver settings as the L2A helper.
    """

    return _get_surface_driven_bestpixel_config(
        cache_root=cache_root,
        bestpixel_endpoint="hls-s30",
        bestpixel_source="hls-s30",
    )
