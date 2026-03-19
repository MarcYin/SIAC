"""
Public SIAC facade.

This module keeps the user-facing API small and delegates orchestration to the
new app/workflow layers.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

from siac.adapters.atmo import CAMSProvider
from siac.adapters.auth import CredentialManager
from siac.adapters.earthdata import earthaccess_source_from_auth
from siac.adapters.satellite import detect_sensor, get_preprocessor
from siac.algorithms.correction import AtmosphericCorrector
from siac.algorithms.solver import MultiGridConfig, MultiGridSolver
from siac.algorithms.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.algorithms.surface.kernel_model import KernelModelDeriver
from siac.algorithms.surface.swir_refine import (
    build_monthly_surface_prior_database,
    query_surface_prior_from_monthly_database,
    resample_geometry_for_surface_prior,
)
from siac.app.assembly import (
    build_preprocessor_runtime,
)
from siac.app.assembly import (
    resolve_brdf_provider as _app_resolve_brdf_provider,
)
from siac.app.assembly import (
    resolve_grid_assembler as _app_resolve_grid_assembler,
)
from siac.app.assembly import (
    resolve_s2_backend as _app_resolve_s2_backend,
)
from siac.app.planning import coerce_aoi_spec as _plan_coerce_aoi_spec
from siac.app.planning import resolve_run_config as _plan_resolve_run_config
from siac.config import SIACConfig
from siac.domain import AtmosphericState, SolvedAtmosphere
from siac.domain.aoi import AOI
from siac.rt.emulator import TwoLayerNNEmulator
from siac.rt.lut import ZarrLUTBackend
from siac.workflows.pipeline import run_pipeline
from siac.workflows.scene import save_output as _workflow_save_output
from siac.workflows.sentinel2 import (
    apply_s2_query_defaults as _workflow_apply_s2_query_defaults,
)
from siac.workflows.sentinel2 import (
    coerce_date as _workflow_coerce_date,
)
from siac.workflows.sentinel2 import (
    coerce_s2_query as _workflow_coerce_s2_query,
)
from siac.workflows.sentinel2 import (
    resolve_s2_input as _workflow_resolve_s2_input,
)
from siac.workflows.sentinel2 import (
    search_sentinel2 as _workflow_search_sentinel2,
)

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    import xarray as xr

    from siac.algorithms.correction import CorrectionResult
    from siac.domain import GeometryAngles, SensorConfig
    from siac.io.s2_data_source import S2Product, S2Query


_resolved_run_config = _plan_resolve_run_config
_coerce_aoi_spec = _plan_coerce_aoi_spec
_coerce_date = _workflow_coerce_date
_apply_s2_query_defaults = _workflow_apply_s2_query_defaults
_coerce_s2_query = _workflow_coerce_s2_query


def _resolve_preprocessor(config: SIACConfig):
    sensor = config.sensor
    if sensor in ("s2", "sentinel2"):
        from inspect import signature

        from siac.adapters.satellite.sentinel2 import Sentinel2Preprocessor

        cloud_mask = getattr(config, "cloud_mask", None)
        cloud_cfg = (
            cloud_mask.model_dump(exclude={"user_callable"})
            if cloud_mask is not None and hasattr(cloud_mask, "model_dump")
            else None
        )
        params = signature(Sentinel2Preprocessor).parameters
        if "config" in params:
            pp = Sentinel2Preprocessor(config={"cloud_mask": cloud_cfg} if cloud_cfg is not None else {})
        else:
            pp = Sentinel2Preprocessor()
            if cloud_cfg is not None and hasattr(pp, "config") and isinstance(pp.config, dict):
                pp.config.setdefault("cloud_mask", cloud_cfg)
        return pp.preprocess
    raise ValueError(f"Unknown sensor: {sensor!r}")


def _select_surface_prior_bands(sensor_config: SensorConfig | None) -> list[Any]:
    if sensor_config is None:
        return list(range(1, 8))
    bands = sensor_config.select_bands_in_range(400.0, 520.0)
    if not bands:
        bands = list(sensor_config.bands[:2])
    return bands


def _select_visible_surface_prior_bands(sensor_config: SensorConfig) -> list[Any]:
    bands = [band for band in sensor_config.bands if 400.0 <= band.center_wavelength < 700.0]
    if bands:
        return bands
    return list(sensor_config.bands[: min(3, len(sensor_config.bands))])


def _select_route_b_query_bands(sensor_config: SensorConfig) -> list[Any]:
    preferred_by_sensor = {
        "MSI": ("B08", "B11", "B12"),
        "OLI": ("B5", "B6", "B7"),
    }
    preferred_names = preferred_by_sensor.get(sensor_config.sensor_id, ())
    selected = [band for name in preferred_names for band in sensor_config.bands if band.name == name]
    if len(selected) == 3:
        return selected

    targets = (860.0, 1610.0, 2200.0)
    remaining = list(sensor_config.bands)
    resolved: list[Any] = []
    for target in targets:
        candidates = [band for band in remaining if band.center_wavelength >= 700.0] or remaining
        match = min(candidates, key=lambda band: abs(band.center_wavelength - target))
        resolved.append(match)
        remaining = [band for band in remaining if band.name != match.name]
    return resolved


def _surface_spectral_mapping_runtime(
    config: SIACConfig,
    *,
    source_bands: list[Any] | tuple[Any, ...] | None = None,
    target_bands: list[Any] | tuple[Any, ...] | None = None,
    context: str = "surface priors",
) -> tuple[object | None, int]:
    settings = config.surface_prior.spectral_mapping
    paths = getattr(config, "paths", None)
    cache_paths = getattr(paths, "caches", None)
    siac_library_root = getattr(settings, "siac_library_root", None)
    if siac_library_root is None and paths is not None:
        siac_library_root = getattr(paths, "spectral_library_root", None)
    rsrf_root = getattr(settings, "rsrf_root", None)
    if rsrf_root is None and paths is not None:
        rsrf_root = getattr(paths, "rsrf_root", None)
    cache_dir = getattr(settings, "cache_dir", None)
    if cache_dir is None and cache_paths is not None:
        cache_dir = getattr(cache_paths, "spectral_mapping", None)

    if source_bands is not None and target_bands is not None and not settings.enabled:
        from siac.algorithms.surface.spectral_mapping import needs_spectral_mapping

        if needs_spectral_mapping(tuple(source_bands), tuple(target_bands)):
            raise ValueError(
                f"surface_prior.spectral_mapping.enabled=false cannot be used for {context} "
                "when BRDF/source bands differ from the requested target bands."
            )

    if not settings.enabled:
        return None, int(settings.k_neighbors)

    from siac.algorithms.surface.spectral_mapping import (
        SpectralMappingConfig as RuntimeSpectralMappingConfig,
    )

    return (
        RuntimeSpectralMappingConfig(
            siac_library_root=siac_library_root,
            rsrf_root=rsrf_root,
            cache_dir=cache_dir,
            neighbor_estimator=settings.neighbor_estimator,
            knn_backend=settings.knn_backend,
            knn_eps=settings.knn_eps,
            min_valid_bands=settings.min_valid_bands,
        ),
        int(settings.k_neighbors),
    )


def _resolve_atmo_provider(config: SIACConfig, auth: CredentialManager | None = None):
    provider_name = config.atmo_prior.provider
    if provider_name == "cams":
        provider = CAMSProvider(
            config.atmo_prior.data_path,
            temporal_interp=config.atmo_prior.temporal_interpolation == "linear",
            download_missing=config.atmo_prior.download_missing,
            auth=auth,
            cache_dir=config.atmo_prior.cache_dir,
        )
        return provider.get_prior
    if provider_name == "merra2":
        from siac.adapters.atmo.merra2 import MERRA2Provider

        provider = MERRA2Provider(
            cache_dir=config.atmo_prior.cache_dir,
            source=earthaccess_source_from_auth(auth),
        )
        return provider.get_prior
    if provider_name == "mcd19":
        from siac.adapters.atmo.mcd19_earthaccess import MCD19AODProvider

        provider = MCD19AODProvider(
            cache_dir=config.atmo_prior.cache_dir,
            source=earthaccess_source_from_auth(auth),
        )
        return provider.get_prior
    if provider_name == "vnp19":
        from siac.adapters.atmo.mcd19_earthaccess import VNP19AODProvider

        provider = VNP19AODProvider(
            cache_dir=config.atmo_prior.cache_dir,
            source=earthaccess_source_from_auth(auth),
        )
        return provider.get_prior
    raise ValueError(f"Unknown atmo provider: {provider_name!r}")


def _resolve_surface_prior_provider(config: SIACConfig, auth: CredentialManager | None = None):
    provider_name = getattr(config.brdf, "provider", "mcd43")

    if provider_name == "mcd43":
        from siac.adapters.brdf.mcd43_earthaccess import MCD43EarthAccessProvider

        provider_cls = MCD43EarthAccessProvider
    elif provider_name == "vnp43":
        from siac.adapters.brdf.vnp43_earthaccess import VNP43EarthAccessProvider

        provider_cls = VNP43EarthAccessProvider
    elif provider_name == "mcd19":
        from siac.adapters.brdf.mcd43_earthaccess import MCD19EarthAccessProvider

        provider_cls = MCD19EarthAccessProvider
    else:
        raise ValueError(f"Unknown BRDF provider for surface prior: {provider_name!r}")

    brdf_prov = provider_cls(
        cache_dir=config.brdf.cache_dir,
        source=earthaccess_source_from_auth(auth),
    )
    method = getattr(config.surface_prior, "method", "kernel_model")
    fallback_deriver = KernelModelDeriver(
        psf_sigma_x=config.surface_prior.psf_sigma_x,
        psf_sigma_y=config.surface_prior.psf_sigma_y,
        apply_psf=config.surface_prior.apply_psf,
    )
    if method == "whittaker":
        deriver = BRDFWhittakerDeriver(
            temporal_lambda=config.surface_prior.whittaker_lambda,
            psf_sigma_x=config.surface_prior.psf_sigma_x,
            psf_sigma_y=config.surface_prior.psf_sigma_y,
            apply_psf=config.surface_prior.apply_psf,
        )

        def _brdf_surface_prior(observation, atmo_prior, rt_model, resolution):
            _ = (atmo_prior, rt_model)
            target_bands = _select_surface_prior_bands(observation.sensor_config)
            spectral_library, spectral_k_neighbors = _surface_spectral_mapping_runtime(
                config,
                source_bands=brdf_prov.source_bands,
                target_bands=target_bands,
                context="Whittaker surface priors",
            )
            brdf_weights = brdf_prov.get_temporal_brdf_parameters(
                bounds=observation.bounds,
                crs=observation.crs,
                obs_time=observation.metadata["observation_time"],
                target_resolution=resolution,
                bands=brdf_prov.source_bands,
                temporal_window=config.brdf.temporal_window,
            )
            return deriver.compute_surface_prior(
                brdf_weights,
                observation.geometry,
                obs_time=observation.metadata["observation_time"],
                source_bands=brdf_prov.source_bands,
                target_bands=target_bands,
                spectral_library=spectral_library,
                spectral_k_neighbors=spectral_k_neighbors,
            )

        _brdf_surface_prior.requires_atmo_prior = False
        return _brdf_surface_prior

    if method == "monthly_database":
        def _monthly_surface_prior(observation, atmo_prior, rt_model, resolution):
            if atmo_prior is None:
                raise ValueError("Route-B monthly_database surface prior requires an atmospheric prior")

            visible_bands = _select_visible_surface_prior_bands(observation.sensor_config)
            query_bands = _select_route_b_query_bands(observation.sensor_config)
            spectral_library, spectral_k_neighbors = _surface_spectral_mapping_runtime(
                config,
                source_bands=brdf_prov.source_bands,
                target_bands=[*visible_bands, *query_bands],
                context="Route-B monthly-database surface priors",
            )
            target_geometry = resample_geometry_for_surface_prior(observation, resolution=resolution)
            database = build_monthly_surface_prior_database(
                observation=observation,
                brdf_provider=brdf_prov,
                resolution=resolution,
                geometry=target_geometry,
                visible_bands=visible_bands,
                query_bands=query_bands,
                spectral_library=spectral_library,
                spectral_k_neighbors=spectral_k_neighbors,
            )
            prior = query_surface_prior_from_monthly_database(
                observation=observation,
                atmo_prior=atmo_prior,
                rt_model=rt_model,
                database=database,
                query_band_names=tuple(band.name for band in query_bands),
                visible_band_names=tuple(band.name for band in visible_bands),
                k_neighbors=3,
            )
            if bool(getattr(prior.mask, "values", prior.mask).any()):
                return prior

            logger.warning(
                "Route-B monthly database produced no valid surface-prior pixels; "
                "falling back to kernel-model BRDF priors."
            )
            brdf_weights = brdf_prov.get_brdf_parameters(
                bounds=observation.bounds,
                crs=observation.crs,
                obs_time=observation.metadata["observation_time"],
                target_resolution=resolution,
                bands=brdf_prov.source_bands,
                temporal_window=config.brdf.temporal_window,
            )
            return fallback_deriver.compute_surface_prior(
                brdf_weights,
                target_geometry,
                source_bands=brdf_prov.source_bands,
                target_bands=visible_bands,
                spectral_library=spectral_library,
                spectral_k_neighbors=spectral_k_neighbors,
            )

        _monthly_surface_prior.requires_atmo_prior = True
        return _monthly_surface_prior

    def _brdf_surface_prior(observation, atmo_prior, rt_model, resolution):
        _ = (atmo_prior, rt_model)
        target_bands = _select_surface_prior_bands(observation.sensor_config)
        spectral_library, spectral_k_neighbors = _surface_spectral_mapping_runtime(
            config,
            source_bands=brdf_prov.source_bands,
            target_bands=target_bands,
            context="surface priors",
        )
        brdf_weights = brdf_prov.get_brdf_parameters(
            bounds=observation.bounds,
            crs=observation.crs,
            obs_time=observation.metadata["observation_time"],
            target_resolution=resolution,
            bands=brdf_prov.source_bands,
            temporal_window=config.brdf.temporal_window,
        )
        return fallback_deriver.compute_surface_prior(
            brdf_weights,
            observation.geometry,
            source_bands=brdf_prov.source_bands,
            target_bands=target_bands,
            spectral_library=spectral_library,
            spectral_k_neighbors=spectral_k_neighbors,
        )

    _brdf_surface_prior.requires_atmo_prior = False
    return _brdf_surface_prior


def _resolve_grid_assembler():
    return _app_resolve_grid_assembler()


def _resolve_solver(config: SIACConfig):
    def _default_solver(inputs, _cfg):
        solver_config = MultiGridConfig(
            aot_gamma=config.solver.aot_gamma,
            tcwv_gamma=config.solver.tcwv_gamma,
            aot_bounds=tuple(config.solver.aot_bounds),
            tcwv_bounds=tuple(config.solver.tcwv_bounds),
        )
        mg_solver = MultiGridSolver(solver_config)
        result = mg_solver.solve(
            inputs.toa,
            inputs.surface_prior,
            inputs.geometry,
            inputs.atmo_prior,
            inputs.rt_model,
            inputs.cloud_mask,
            inputs.bands,
        )
        solved_atmo = inputs.atmo_prior.with_updated_aot_tcwv(
            aot=result.aot,
            tcwv=result.tcwv,
            aot_unc=result.aot_unc,
            tcwv_unc=result.tcwv_unc,
        )
        return SolvedAtmosphere(
            atmo_state=solved_atmo,
            aot=result.aot,
            tcwv=result.tcwv,
            aot_unc=result.aot_unc,
            tcwv_unc=result.tcwv_unc,
            cost_final=float(result.final_cost),
            n_iterations=result.n_iterations,
            converged=result.success,
        )

    return _default_solver


def _resample_field_to_template(field, template):
    if field.shape == template.shape:
        return field

    if (
        all(dim in field.dims for dim in ("y", "x"))
        and all(dim in template.dims for dim in ("y", "x"))
        and all(coord in field.coords for coord in ("y", "x"))
        and all(coord in template.coords for coord in ("y", "x"))
    ):
        try:
            return field.interp(y=template.coords["y"], x=template.coords["x"], method="linear")
        except Exception:
            pass

    import numpy as np
    from scipy import ndimage

    src = np.asarray(field.values, dtype=np.float32)
    if src.ndim != 2:
        return field

    h_out = int(template.sizes["y"])
    w_out = int(template.sizes["x"])
    if src.shape[0] == 0 or src.shape[1] == 0:
        out = np.full((h_out, w_out), np.nan, dtype=np.float32)
    else:
        out = ndimage.zoom(src, (h_out / src.shape[0], w_out / src.shape[1]), order=1)
        out = out[:h_out, :w_out]
        if out.shape != (h_out, w_out):
            padded = np.full((h_out, w_out), np.nan, dtype=np.float32)
            padded[: out.shape[0], : out.shape[1]] = out
            out = padded

    return field.__class__(
        out,
        dims=template.dims,
        coords={d: template.coords[d] for d in template.dims if d in template.coords},
    )


def _resolve_corrector(_config: SIACConfig):
    def _default_corrector(obs, solved, rt_model):
        corrector_obj = AtmosphericCorrector(rt_model, obs.sensor_config)
        first_band = obs.toa[next(iter(obs.toa.data_vars))]
        atmo = solved.atmo_state
        matched_atmo = AtmosphericState(
            aot=_resample_field_to_template(atmo.aot, first_band),
            tcwv=_resample_field_to_template(atmo.tcwv, first_band),
            tco3=_resample_field_to_template(atmo.tco3, first_band),
            aot_unc=_resample_field_to_template(atmo.aot_unc, first_band),
            tcwv_unc=_resample_field_to_template(atmo.tcwv_unc, first_band),
            tco3_unc=_resample_field_to_template(atmo.tco3_unc, first_band),
            elevation=_resample_field_to_template(atmo.elevation, first_band),
        )
        return corrector_obj.correct(obs.toa, obs.geometry, matched_atmo, obs.cloud_mask)

    return _default_corrector


def _resolve_rt_model_for_pipeline(
    config: SIACConfig,
    auth: CredentialManager | None = None,
    *,
    sensor_config: SensorConfig | None = None,
):
    sensor_defaults = {
        "s2": ("MSI", "S2A"),
        "s2a": ("MSI", "S2A"),
        "s2b": ("MSI", "S2B"),
        "s2c": ("MSI", "S2C"),
        "sentinel2": ("MSI", "S2A"),
        "l8": ("OLI", "L8"),
        "l9": ("OLI", "L9"),
        "auto": ("MSI", "S2A"),
    }
    rt_config = config.rt_model
    if sensor_config is not None:
        sid, satid = sensor_config.sensor_id, sensor_config.satellite_id
    else:
        sid, satid = sensor_defaults.get(getattr(config, "sensor", None), ("MSI", "S2A"))

    paths = getattr(config, "paths", None)
    emulator_dir = getattr(rt_config, "emulator_dir", None) or getattr(paths, "emulator_dir", None)
    lut_path = getattr(rt_config, "lut_path", None) or getattr(paths, "lut_path", None)
    backend = rt_config.backend

    if backend == "emulator":
        try:
            return TwoLayerNNEmulator(
                emulator_dir=emulator_dir,
                sensor_id=sid,
                satellite_id=satid,
            )
        except Exception as exc:
            if not rt_config.fallback_to_lut or not lut_path:
                raise ValueError(
                    "Cannot resolve emulator RT model and LUT fallback is unavailable."
                ) from exc
            logger.warning("Emulator RT model unavailable; falling back to LUT backend.")
            backend = "lut"

    if backend == "lut" and lut_path:
        storage_options = dict(rt_config.lut_storage_options)
        if (
            auth is not None
            and auth.aws().has_credentials()
            and "key" not in storage_options
            and str(lut_path).startswith("s3://")
        ):
            storage_options.update(auth.aws().storage_options())
        return ZarrLUTBackend(
            lut_path,
            interpolation_method=rt_config.lut_interpolation,
            storage_options=storage_options,
        )

    raise ValueError(f"Cannot resolve RT model from config: backend={rt_config.backend!r}")


def _resolve_s2_backend(config: SIACConfig, *, auth: CredentialManager | None = None):
    return _app_resolve_s2_backend(config, auth=auth)


def resolve_s2_input(
    query: S2Query | str | Path,
    config: SIACConfig,
    *,
    auth: CredentialManager | None = None,
) -> Path:
    return _workflow_resolve_s2_input(
        query,
        config,
        auth=auth,
        resolve_s2_backend_fn=_resolve_s2_backend,
    )


def search_sentinel2(
    *,
    tile: str | None = None,
    date: datetime | str | None = None,
    date_value: datetime | str | None = None,
    start_date: datetime | str | None = None,
    end_date: datetime | str | None = None,
    bbox: tuple[float, float, float, float] | None = None,
    max_cloud_cover: float = 80.0,
    backend: str = "cdse",
    config: SIACConfig | None = None,
    auth: CredentialManager | None = None,
) -> list[S2Product]:
    return _workflow_search_sentinel2(
        tile=tile,
        date=date,
        date_value=date_value,
        start_date=start_date,
        end_date=end_date,
        bbox=bbox,
        max_cloud_cover=max_cloud_cover,
        backend=backend,
        config=config,
        auth=auth,
        resolve_s2_backend_fn=_resolve_s2_backend,
    )


class SIAC:
    """Thin session facade over the SIAC app/workflow layers."""

    def __init__(self, config: SIACConfig):
        self.config = config
        self._auth = CredentialManager.from_config(config)
        self._preprocessor = None
        self._atmo_provider = None
        self._brdf_provider = None
        self._rt_model = None
        self._solver = None

    @classmethod
    def from_config(cls, config_path: str | Path) -> SIAC:
        return cls(SIACConfig.from_file(config_path))

    @classmethod
    def from_defaults(cls, sensor: str = "auto") -> SIAC:
        return cls(SIACConfig(sensor=sensor))

    def process(
        self,
        input_path: str | Path,
        output_path: str | Path | None = None,
        *,
        aoi: AOI | Path | str | tuple[float, float, float, float] | list[float] | None = None,
    ) -> CorrectionResult:
        input_path = Path(input_path)
        logger.info("Processing: %s", input_path)
        resolved_config = _resolved_run_config(
            self.config,
            input_path=input_path,
            output_path=output_path,
            sensor=self.config.sensor,
            aoi=aoi if aoi is not None else self.config.aoi,
        )
        self._auth = CredentialManager.from_config(resolved_config)
        runtime_aoi = _coerce_aoi_spec(resolved_config.aoi)
        sensor_name = resolved_config.sensor if resolved_config.sensor != "auto" else detect_sensor(input_path)
        preprocessor_runtime = build_preprocessor_runtime(
            resolved_config,
            input_path=input_path,
            sensor=sensor_name,
            default_aoi_resolver=self._resolve_aoi,
            detect_sensor_fn=detect_sensor,
            get_preprocessor_fn=get_preprocessor,
        )
        atmo_provider = _resolve_atmo_provider(resolved_config, auth=self._auth)
        surface_prior_provider = _resolve_surface_prior_provider(resolved_config, auth=self._auth)
        grid_assembler = _resolve_grid_assembler()
        solver = _resolve_solver(resolved_config)
        corrector = _resolve_corrector(resolved_config)
        rt_model = _resolve_rt_model_for_pipeline(
            resolved_config,
            auth=self._auth,
            sensor_config=preprocessor_runtime.sensor_config,
        )
        result = run_pipeline(
            input_path,
            runtime_aoi,
            resolved_config,
            preprocessor=preprocessor_runtime.preprocessor,
            atmo_provider=atmo_provider,
            surface_prior_provider=surface_prior_provider,
            grid_assembler=grid_assembler,
            solver=solver,
            corrector=corrector,
            rt_model=rt_model,
        )
        if output_path is not None:
            self._save_output(result, output_path)
        logger.info("Complete. Mean AOT: %.3f", float(result.aot.mean()))
        return result

    def _resolve_aoi(self, toa: xr.Dataset) -> AOI:
        if self.config.aoi is not None:
            return _coerce_aoi_spec(self.config.aoi)
        first_var = list(toa.data_vars)[0]
        return AOI.from_raster(toa[first_var])

    def _get_atmospheric_prior(self, aoi: AOI, metadata: dict) -> AtmosphericState:
        provider_name = self.config.atmo_prior.provider
        if provider_name == "cams" or provider_name in {"merra2", "mcd19", "vnp19"}:
            provider_fn = _resolve_atmo_provider(self.config, auth=self._auth)
        else:
            logger.warning(
                f"Unknown atmospheric provider '{provider_name}', falling back to CAMS"
            )
            fallback_cfg = self.config.with_overrides(providers={"atmo": {"kind": "cams"}})
            provider_fn = _resolve_atmo_provider(fallback_cfg, auth=self._auth)

        bounds = aoi.get_bounds()
        crs = aoi.crs
        resolution = 10.0
        obs_time = metadata.get("observation_time", datetime.now())
        return provider_fn(bounds, crs, obs_time, resolution)

    def _get_surface_prior(self, aoi: AOI, geometry: GeometryAngles, metadata: dict):
        method = getattr(self.config.surface_prior, "method", "kernel_model")
        if method == "monthly_database":
            raise RuntimeError(
                "Route-B monthly-database surface priors require the full ObservationBundle "
                "and are only available through the main pipeline."
            )

        _ = self._get_brdf_provider()
        obs = SimpleNamespace(
            bounds=aoi.get_bounds(),
            crs=aoi.crs,
            metadata=metadata,
            geometry=geometry,
            sensor_config=metadata.get("sensor_config"),
        )
        surface_fn = _resolve_surface_prior_provider(self.config, auth=self._auth)
        return surface_fn(obs, None, "rt", 500.0)

    def _get_brdf_provider(self):
        if self._brdf_provider is not None:
            return self._brdf_provider

        provider_name = self.config.brdf.provider
        if provider_name in {"mcd43", "vnp43", "mcd19", "gee"}:
            self._brdf_provider = _app_resolve_brdf_provider(self.config, auth=self._auth)
        else:
            raise ValueError(
                f"Unknown BRDF provider '{provider_name}'. "
                "Available: 'mcd43', 'vnp43', 'mcd19', 'gee'"
            )
        return self._brdf_provider

    def _get_rt_model(self, sensor_config: SensorConfig):
        if self._rt_model is not None:
            return self._rt_model

        rt_config = self.config.rt_model
        emulator_dir = getattr(rt_config, "emulator_dir", None) or self.config.paths.emulator_dir
        lut_path = getattr(rt_config, "lut_path", None) or self.config.paths.lut_path
        if rt_config.backend == "emulator":
            try:
                self._rt_model = TwoLayerNNEmulator(
                    emulator_dir=emulator_dir,
                    sensor_id=sensor_config.sensor_id,
                    satellite_id=sensor_config.satellite_id,
                )
                return self._rt_model
            except Exception as exc:
                logger.warning("Emulator not available: %s, falling back to LUT", exc)
                rt_config.backend = "lut"

        if rt_config.backend == "lut" and lut_path:
            self._rt_model = ZarrLUTBackend(
                lut_path,
                interpolation_method=rt_config.lut_interpolation,
                storage_options=rt_config.lut_storage_options,
            )
            return self._rt_model

        return _resolve_rt_model_for_pipeline(self.config, auth=self._auth, sensor_config=sensor_config)

    def _solve_atmosphere(self, toa, surface_prior, geometry, atmo_prior, rt_model, cloud_mask, sensor_config):
        if self._solver is None:
            solver_config = MultiGridConfig(
                aot_gamma=self.config.solver.aot_gamma,
                tcwv_gamma=self.config.solver.tcwv_gamma,
                aot_bounds=tuple(self.config.solver.aot_bounds),
                tcwv_bounds=tuple(self.config.solver.tcwv_bounds),
            )
            self._solver = MultiGridSolver(solver_config)

        bands = [sensor_config.get_band(name) for name in sensor_config.band_names[:6]]
        return self._solver.solve(toa, surface_prior, geometry, atmo_prior, rt_model, cloud_mask, bands)

    def _save_output(self, result: CorrectionResult, output_path: Path):
        _workflow_save_output(result, output_path)


def process_sentinel2(input_path: str, output_path: str | None = None, **kwargs) -> CorrectionResult:
    return SIAC.from_defaults(sensor="s2").process(input_path, output_path, **kwargs)


def process_landsat8(input_path: str, output_path: str | None = None, **kwargs) -> CorrectionResult:
    return SIAC.from_defaults(sensor="l8").process(input_path, output_path, **kwargs)


def siac_process_s2(
    config: SIACConfig,
    query: S2Query | str | Path,
    *,
    output_path: str | Path | None = None,
    aoi: AOI | Path | str | tuple[float, float, float, float] | list[float] | None = None,
    auth: CredentialManager | None = None,
) -> CorrectionResult:
    resolved_config = _resolved_run_config(
        config,
        sensor=config.sensor if config.sensor != "auto" else "s2",
        aoi=aoi if aoi is not None else config.aoi,
        s2_query=query,
    )
    auth_obj = auth or CredentialManager.from_config(resolved_config)
    input_path = resolve_s2_input(query, config, auth=auth_obj)
    siac_obj = SIAC(config)
    siac_obj._auth = auth_obj
    if aoi is None:
        return siac_obj.process(input_path, output_path)
    return siac_obj.process(input_path, output_path, aoi=aoi)


def siac_process(
    config: SIACConfig,
    input_path: Path,
    *,
    aoi: AOI | Path | str | tuple[float, float, float, float] | list[float] | None = None,
    auth: CredentialManager | None = None,
    preprocessor=None,
    atmo_provider=None,
    surface_prior_provider=None,
    grid_assembler=None,
    solver=None,
    corrector=None,
    rt_model: Any | None = None,
) -> CorrectionResult:
    resolved_config = _resolved_run_config(
        config,
        input_path=input_path,
        output_path=None,
        sensor=config.sensor,
        aoi=aoi if aoi is not None else config.aoi,
    )
    auth_obj = auth or CredentialManager.from_config(resolved_config)
    preprocessor = preprocessor or _resolve_preprocessor(resolved_config)
    atmo_provider = atmo_provider or _resolve_atmo_provider(resolved_config, auth=auth_obj)
    surface_prior_provider = surface_prior_provider or _resolve_surface_prior_provider(resolved_config, auth=auth_obj)
    grid_assembler = grid_assembler or _resolve_grid_assembler()
    solver = solver or _resolve_solver(resolved_config)
    corrector = corrector or _resolve_corrector(resolved_config)
    rt_model = rt_model or _resolve_rt_model_for_pipeline(resolved_config, auth=auth_obj)
    runtime_aoi = _coerce_aoi_spec(resolved_config.aoi)
    return run_pipeline(
        Path(input_path),
        runtime_aoi,
        resolved_config,
        preprocessor=preprocessor,
        atmo_provider=atmo_provider,
        surface_prior_provider=surface_prior_provider,
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
    )


__all__ = [
    "SIAC",
    "process_sentinel2",
    "process_landsat8",
    "resolve_s2_input",
    "search_sentinel2",
    "siac_process_s2",
    "siac_process",
    "_resolved_run_config",
    "_coerce_aoi_spec",
    "_resolve_preprocessor",
    "_resolve_atmo_provider",
    "_resolve_surface_prior_provider",
    "_resolve_grid_assembler",
    "_resolve_solver",
    "_resolve_corrector",
    "_resolve_rt_model_for_pipeline",
    "_resolve_s2_backend",
    "_coerce_date",
    "_apply_s2_query_defaults",
    "_coerce_s2_query",
]
