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


#: Copernicus GLO-30 digital elevation model as a remote VRT. Terrain height
#: sets surface pressure, which scales Rayleigh scattering — the dominant
#: atmospheric signal in the blue bands the AOD solve relies on.
COPERNICUS_GLO30_DEM_VRT = (
    "/vsicurl/https://raw.githubusercontent.com/MarcYin/Copernicus_GLO_30_DEM_VRT/refs/heads/main/copernicus_GLO_30_dem.vrt"
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
            dem=COPERNICUS_GLO30_DEM_VRT,
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
                "surface_driven_backstop_escape_enabled": False,
                "surface_driven_backstop_escape_low_aot": 0.08,
                "surface_driven_backstop_escape_cost_threshold": 20.0,
                "surface_driven_backstop_escape_delta": 0.05,
                "surface_driven_backstop_escape_band_spread": 0.25,
                "surface_driven_backstop_escape_min_jump": 0.35,
                "surface_driven_backstop_escape_cost_ratio": 1.8,
                "surface_driven_reference_tcwv": 2.0,
                "surface_driven_scene_mean_geometry": True,
                "surface_driven_resolve_on_prior_grid": False,
                "surface_driven_pool_radius_m": 600.0,
                "water_mask_buffer_pixels": 32,
                "surface_driven_solve_bands": ["B01", "B02", "B04"],
                "surface_driven_cost_mode": "auto2",
                "surface_driven_aot_axis": "acixthree",
                "surface_driven_allow_cloud_retrieval": False,
                "surface_driven_ignore_cloud_water": False,
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
    # Spectral interpolation between the independently calibrated B02/B04
    # terms.  It is dormant in the shipped B02/B04 preset and used only when a
    # caller explicitly adds B03 to the same visible predictor.
    "B03": (-0.0006, 0.0235),
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


def get_surface_driven_v1_config(
    *,
    prepared_library_path: str | Path,
    cache_root: str | Path | None = None,
) -> SIACConfig:
    """Return the validated v1 surface-driven AOD recipe.

    This is the best-measured configuration on the 152-matchup AERONET campaign:
    **84.6% within-EE, AOD RMSE 0.110, bias +0.001**, against 79.2% / 0.120 /
    −0.021 for the same recipe on a single-source aerosol prior. It combines
    four results, each measured by holding everything else fixed:

    * **Prepared surface library, corrected in the solve's own RT space.** The
      library is built offline from clear-sky Sentinel-2 L1C selected per pixel
      by same-day aerosol loading, then atmospherically corrected with the same
      6S + CCI-climatology model the solver uses. Correcting it in a different
      RT space instead costs 5.4 points (74.3% vs 79.7%) on identical imagery,
      which is why ``prepared_library_path`` is required rather than optional:
      the library declares its RT space and the pipeline checks it.
    * **Exact CCI climatology aerosol with a bounded absorbing fraction.**
      Scene-resolved species cut AOD RMSE by a quarter against a fixed
      continental profile and fix the high-AOD over-retrievals; capping the
      soot-like fraction halves the mixture overshoots. Letting the solver
      *select* species from the cost instead is degenerate — it prefers a
      compensating-wrong model — so the mixture is prescribed, not searched.
    * **A 0.006 surface-uncertainty floor.** The library is accurate enough to
      be trusted tightly; at 0.015 the cost curve flattens and the retrieval
      falls back on the aerosol prior (75.9% vs 72.4% measured on this route).
    * **A fused ``max(MAIAC, CAMS)`` aerosol prior.** Both sources under-read,
      independently, so their maximum removes most of the shared bias without
      any fitting against reference data. This is the single largest term:
      +5.4 points, concentrated in the prior-limited moderate band (69% → 84%).

    ``prepared_library_path`` points at the per-scene library store; see
    :class:`siac.adapters.surface_library.PreparedSurfaceLibrary` for its layout.
    """

    return _get_surface_driven_bestpixel_config(
        cache_root=cache_root,
        bestpixel_endpoint="pc",
        bestpixel_source="l1c",
    ).with_overrides(
        paths={
            # Real terrain, not the sea-level sentinel the bestpixel presets
            # inherit. Elevation sets surface pressure and therefore Rayleigh
            # optical depth, which is the dominant atmospheric term in the blue
            # bands this recipe solves AOD from; assuming sea level over an
            # elevated site mis-attributes molecular scattering to aerosol.
            "dem": COPERNICUS_GLO30_DEM_VRT,
        },
        providers={
            "atmo": {
                "kind": "mcd19",
                "maiac_best_quality_qa": True,
                "fuse_aod_with": ("cams",),
                "fuse_aod_op": "max",
            },
        },
        algorithms={
            "surface_prior": {
                "prepared_library_path": Path(prepared_library_path).expanduser(),
                "bestpixel_window_reduction": "window",
                "bestpixel_predict_visible": True,
                "bestpixel_predict_visible_bands": ("B02", "B03", "B04"),
                "bestpixel_predict_visible_uncertainty_floor": 0.006,
                "bestpixel_predict_visible_anchor_source": "atmo_prior",
                # The library is already in the solver's RT space, so the
                # empirical cross-space debias must stay off: correcting the
                # same offset twice re-introduces it.
                "bestpixel_predict_visible_debias": None,
                # The validated predictor is the 20-tree ensemble; the single
                # tree's noisier realizations widen sigma ~15% and cost the
                # thick-aerosol sites one AOD node.
                "bestpixel_predict_visible_model": "extra_trees_20",
            },
            "rt": {
                "backend": "sixs",
                "fallback_to_lut": False,
            },
            "solver": {
                "surface_driven_solve_bands": ("B02", "B03", "B04"),
                "surface_driven_tau_dependent_prior": True,
                "surface_driven_aerosol_species": "cci_climatology_exact",
                # The validated recipe scores the plain chi-squared cost. The
                # bestpixel base preset selects "auto2", whose shape-augmented
                # curve is ~4x flatter at thick-aerosol scenes, letting the
                # backstop drag the minimum well below the surface solution
                # (Cinzana: 0.80 vs 1.10 with identical priors and sigma).
                "surface_driven_cost_mode": "chi2",
                # One AOI-mean geometry for anchor correction and cost-cube RT,
                # as validated. Sentinel-2 view angles vary little over a solve
                # AOI, and per-pixel geometry shifts retrievals by ~one AOD node
                # relative to the validated configuration.
                "surface_driven_scene_mean_geometry": True,
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
