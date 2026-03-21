"""Assembly layer: resolve configured components into runtime callables."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from siac.adapters.atmo import CAMSProvider
from siac.adapters.earthdata import earthaccess_source_from_auth
from siac.adapters.output import ConfiguredOutputWriter
from siac.adapters.rt import build_rt_model
from siac.adapters.s2_backend import build_s2_backend
from siac.adapters.satellite import detect_sensor, get_preprocessor
from siac.algorithms.correction import AtmosphericCorrector, CorrectionResult
from siac.algorithms.solver import MultiGridConfig, MultiGridSolver
from siac.algorithms.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.algorithms.surface.kernel_model import KernelModelDeriver
from siac.algorithms.surface.swir_refine import (
    build_monthly_surface_prior_database,
    query_surface_prior_from_monthly_database,
    resample_geometry_for_surface_prior,
)
from siac.app.registry import (
    ATMO_PROVIDER_REGISTRY,
    BRDF_PROVIDER_REGISTRY,
    S2_BACKEND_REGISTRY,
    SURFACE_PRIOR_METHOD_REGISTRY,
)
from siac.domain.aoi import AOI
from siac.runtime import (
    AtmosphericState,
    ObservationBundle,
    SolvedAtmosphere,
    SolverInputBundle,
    SurfacePrior,
)

if TYPE_CHECKING:
    from pathlib import Path

    import xarray as xr

    from siac.adapters.auth import CredentialManager
    from siac.domain.protocols import AtmosphericPriorProvider
    from siac.domain.sensors import SensorConfig
    from siac.workflows.pipeline import (
        AtmoPriorFn,
        CorrectorFn,
        GridAssemblerFn,
        PreprocessorFn,
        SolverFn,
        SurfacePriorFn,
    )

logger = logging.getLogger(__name__)


if TYPE_CHECKING:
    from typing import Protocol

    class SurfacePriorCallable(Protocol):
        requires_atmo_prior: bool

        def __call__(
            self,
            observation: ObservationBundle,
            atmo_prior: AtmosphericState | None,
            rt_model: Any,
            resolution: float,
        ) -> SurfacePrior:
            ...


@dataclass(frozen=True)
class PreprocessorRuntime:
    """Preprocessor callable plus its sensor configuration."""

    preprocessor: PreprocessorFn
    sensor_config: SensorConfig


def _select_surface_prior_bands(sensor_config: SensorConfig | None) -> list[Any]:
    if sensor_config is None:
        return list(range(1, 8))
    selected: list[Any] = list(sensor_config.select_bands_in_range(400.0, 520.0))
    if selected:
        return selected
    fallback: list[Any] = list(sensor_config.bands[:2])
    return fallback


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
    config: Any,
    *,
    source_bands: list[Any] | tuple[Any, ...] | None = None,
    target_bands: list[Any] | tuple[Any, ...] | None = None,
    context: str = "surface priors",
) -> tuple[Any | None, int]:
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


def _mark_surface_prior_metadata(provider: SurfacePriorFn, *, requires_atmo_prior: bool) -> SurfacePriorFn:
    cast("SurfacePriorCallable", provider).requires_atmo_prior = requires_atmo_prior
    return provider


def build_preprocessor_runtime(
    config: Any,
    *,
    input_path: Path | None = None,
    sensor: str | None = None,
    default_aoi_resolver: Any | None = None,
    detect_sensor_fn: Any | None = None,
    get_preprocessor_fn: Any | None = None,
) -> PreprocessorRuntime:
    """Build the M1 preprocessor runtime."""
    if detect_sensor_fn is None:
        detect_sensor_fn = detect_sensor
    if get_preprocessor_fn is None:
        get_preprocessor_fn = get_preprocessor

    sensor_name = sensor or config.sensor
    if sensor_name == "auto":
        if input_path is None:
            raise ValueError("Cannot resolve preprocessor for sensor='auto' without an input path.")
        sensor_name = detect_sensor_fn(input_path)

    cloud_mask_config = config.cloud_mask.model_dump(exclude={"user_callable"})
    paths = getattr(config, "paths", None)
    rsrf_root = getattr(paths, "rsrf_root", None)
    preprocessor_config = {"cloud_mask": cloud_mask_config}
    if rsrf_root is not None:
        preprocessor_config["rsrf_root"] = rsrf_root
    try:
        preprocessor_obj = get_preprocessor_fn(sensor_name, config=preprocessor_config)
    except KeyError as exc:
        raise ValueError(f"Unknown sensor: {sensor_name!r}") from exc
    except TypeError:
        try:
            preprocessor_obj = get_preprocessor_fn(sensor_name)
        except KeyError as exc:
            raise ValueError(f"Unknown sensor: {sensor_name!r}") from exc
        if hasattr(preprocessor_obj, "config") and isinstance(preprocessor_obj.config, dict):
            preprocessor_obj.config.setdefault("cloud_mask", cloud_mask_config)
            if rsrf_root is not None:
                preprocessor_obj.config.setdefault("rsrf_root", rsrf_root)

    sensor_config = preprocessor_obj.sensor_config

    def _resolve_default_aoi(toa: xr.Dataset) -> AOI:
        if callable(default_aoi_resolver):
            return cast("AOI", default_aoi_resolver(toa))
        return AOI.from_raster(toa[list(toa.data_vars)[0]])

    def _preprocessor(path: Path, aoi: AOI | None = None) -> ObservationBundle:
        raw = preprocessor_obj.preprocess(path)
        toa = raw["toa"]
        resolved_aoi = aoi or _resolve_default_aoi(toa)
        return ObservationBundle(
            toa=toa,
            geometry=raw["geometry"],
            cloud_mask=raw["cloud_mask"],
            sensor_config=sensor_config,
            metadata=raw["metadata"],
            crs=str(resolved_aoi.crs),
            bounds=resolved_aoi.get_bounds(),
        )

    return PreprocessorRuntime(preprocessor=_preprocessor, sensor_config=sensor_config)


def resolve_preprocessor(
    config: Any,
    *,
    input_path: Path | None = None,
    sensor: str | None = None,
    default_aoi_resolver: Any | None = None,
) -> PreprocessorFn:
    return build_preprocessor_runtime(
        config,
        input_path=input_path,
        sensor=sensor,
        default_aoi_resolver=default_aoi_resolver,
    ).preprocessor


def resolve_output_writer(config: Any) -> ConfiguredOutputWriter:
    return ConfiguredOutputWriter(config.output.defaults)


@ATMO_PROVIDER_REGISTRY.register("cams")
def _build_cams_provider(
    config: Any,
    auth: CredentialManager | None = None,
) -> AtmosphericPriorProvider:
    return CAMSProvider(
        config.atmo_prior.data_path,
        temporal_interp=config.atmo_prior.temporal_interpolation == "linear",
        download_missing=config.atmo_prior.download_missing,
        auth=auth,
        cache_dir=config.atmo_prior.cache_dir,
    )


@ATMO_PROVIDER_REGISTRY.register("merra2")
def _build_merra2_provider(
    config: Any,
    auth: CredentialManager | None = None,
) -> AtmosphericPriorProvider:
    from siac.adapters.atmo.merra2 import MERRA2Provider

    return MERRA2Provider(
        cache_dir=config.atmo_prior.cache_dir,
        source=earthaccess_source_from_auth(auth),
    )


@ATMO_PROVIDER_REGISTRY.register("mcd19")
def _build_mcd19_provider(
    config: Any,
    auth: CredentialManager | None = None,
) -> AtmosphericPriorProvider:
    from siac.adapters.atmo.mcd19_earthaccess import MCD19AODProvider

    return MCD19AODProvider(
        cache_dir=config.atmo_prior.cache_dir,
        source=earthaccess_source_from_auth(auth),
    )


@ATMO_PROVIDER_REGISTRY.register("vnp19")
def _build_vnp19_provider(
    config: Any,
    auth: CredentialManager | None = None,
) -> AtmosphericPriorProvider:
    from siac.adapters.atmo.mcd19_earthaccess import VNP19AODProvider

    return VNP19AODProvider(
        cache_dir=config.atmo_prior.cache_dir,
        source=earthaccess_source_from_auth(auth),
    )


@BRDF_PROVIDER_REGISTRY.register("mcd43")
def _build_mcd43_provider(config: Any, auth: CredentialManager | None = None) -> Any:
    from siac.adapters.brdf.mcd43_earthaccess import MCD43EarthAccessProvider

    return MCD43EarthAccessProvider(
        cache_dir=config.brdf.cache_dir,
        source=earthaccess_source_from_auth(auth),
    )


@BRDF_PROVIDER_REGISTRY.register("vnp43")
def _build_vnp43_provider(config: Any, auth: CredentialManager | None = None) -> Any:
    from siac.adapters.brdf.vnp43_earthaccess import VNP43EarthAccessProvider

    return VNP43EarthAccessProvider(
        cache_dir=config.brdf.cache_dir,
        source=earthaccess_source_from_auth(auth),
    )


@BRDF_PROVIDER_REGISTRY.register("mcd19")
def _build_mcd19_brdf_provider(config: Any, auth: CredentialManager | None = None) -> Any:
    from siac.adapters.brdf.mcd43_earthaccess import MCD19EarthAccessProvider

    return MCD19EarthAccessProvider(
        cache_dir=config.brdf.cache_dir,
        source=earthaccess_source_from_auth(auth),
    )


@BRDF_PROVIDER_REGISTRY.register("gee")
def _build_gee_brdf_provider(_config: Any, _auth: CredentialManager | None = None) -> Any:
    from siac.adapters.brdf.gee_stub import GEEBRDFProvider

    return GEEBRDFProvider()


def resolve_brdf_provider(config: Any, *, auth: CredentialManager | None = None) -> Any:
    provider_name = getattr(config.brdf, "provider", "mcd43")
    return BRDF_PROVIDER_REGISTRY.get(provider_name)(config, auth)


@SURFACE_PRIOR_METHOD_REGISTRY.register("kernel_model")
def _build_kernel_surface_prior(config: Any, brdf_prov: Any) -> SurfacePriorFn:
    deriver = KernelModelDeriver(
        psf_sigma_x=config.surface_prior.psf_sigma_x,
        psf_sigma_y=config.surface_prior.psf_sigma_y,
        apply_psf=config.surface_prior.apply_psf,
    )

    def _surface_prior(
        observation: ObservationBundle,
        atmo_prior: AtmosphericState | None,
        rt_model: Any,
        resolution: float,
    ) -> SurfacePrior:
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
        return deriver.compute_surface_prior(
            brdf_weights,
            observation.geometry,
            source_bands=brdf_prov.source_bands,
            target_bands=target_bands,
            spectral_library=spectral_library,
            spectral_k_neighbors=spectral_k_neighbors,
        )

    return _mark_surface_prior_metadata(_surface_prior, requires_atmo_prior=False)


@SURFACE_PRIOR_METHOD_REGISTRY.register("whittaker")
def _build_whittaker_surface_prior(config: Any, brdf_prov: Any) -> SurfacePriorFn:
    deriver = BRDFWhittakerDeriver(
        temporal_lambda=config.surface_prior.whittaker_lambda,
        psf_sigma_x=config.surface_prior.psf_sigma_x,
        psf_sigma_y=config.surface_prior.psf_sigma_y,
        apply_psf=config.surface_prior.apply_psf,
    )

    def _surface_prior(
        observation: ObservationBundle,
        atmo_prior: AtmosphericState | None,
        rt_model: Any,
        resolution: float,
    ) -> SurfacePrior:
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

    return _mark_surface_prior_metadata(_surface_prior, requires_atmo_prior=False)


@SURFACE_PRIOR_METHOD_REGISTRY.register("monthly_database")
def _build_monthly_surface_prior(config: Any, brdf_prov: Any) -> SurfacePriorFn:
    fallback_deriver = KernelModelDeriver(
        psf_sigma_x=config.surface_prior.psf_sigma_x,
        psf_sigma_y=config.surface_prior.psf_sigma_y,
        apply_psf=config.surface_prior.apply_psf,
    )

    def _surface_prior(
        observation: ObservationBundle,
        atmo_prior: AtmosphericState | None,
        rt_model: Any,
        resolution: float,
    ) -> SurfacePrior:
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
        if bool(np.asarray(prior.mask.values).any()):
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

    return _mark_surface_prior_metadata(_surface_prior, requires_atmo_prior=True)


def resolve_atmo_provider(config: Any, auth: CredentialManager | None = None) -> AtmoPriorFn:
    provider_name = config.atmo_prior.provider
    provider = ATMO_PROVIDER_REGISTRY.get(provider_name)(config, auth)
    return cast("AtmosphericPriorProvider", provider).get_prior


def resolve_surface_prior_provider(config: Any, auth: CredentialManager | None = None) -> SurfacePriorFn:
    brdf_prov = resolve_brdf_provider(config, auth=auth)
    method = getattr(config.surface_prior, "method", "kernel_model")
    return cast("SurfacePriorFn", SURFACE_PRIOR_METHOD_REGISTRY.get(method)(config, brdf_prov))


def resolve_grid_assembler() -> GridAssemblerFn:
    from siac.algorithms.grid.assembler import assemble_grids

    return assemble_grids


def resolve_solver(config: Any) -> SolverFn:
    def _default_solver(inputs: SolverInputBundle, _config: Any) -> SolvedAtmosphere:
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


def _resample_field_to_template(field: Any, template: Any) -> Any:
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

    src = np.asarray(field.values, dtype=np.float32)
    if src.ndim != 2:
        return field

    from scipy import ndimage

    h_out = int(template.sizes["y"])
    w_out = int(template.sizes["x"])
    if src.shape[0] == 0 or src.shape[1] == 0:
        out: np.ndarray[Any, Any] = np.full((h_out, w_out), np.nan, dtype=np.float32)
    else:
        out = ndimage.zoom(src, (h_out / src.shape[0], w_out / src.shape[1]), order=1)
        out = out[:h_out, :w_out]
        if out.shape != (h_out, w_out):
            padded: np.ndarray[Any, Any] = np.full((h_out, w_out), np.nan, dtype=np.float32)
            padded[: out.shape[0], : out.shape[1]] = out
            out = padded

    return field.__class__(
        out,
        dims=template.dims,
        coords={d: template.coords[d] for d in template.dims if d in template.coords},
    )


def resolve_corrector(_config: Any) -> CorrectorFn:
    def _default_corrector(
        obs: ObservationBundle,
        solved: SolvedAtmosphere,
        rt_model: Any,
    ) -> CorrectionResult:
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


@S2_BACKEND_REGISTRY.register("cdse")
def _build_cdse_s2_backend(config: Any, auth: CredentialManager | None = None) -> Any:
    return build_s2_backend(config, auth=auth)


@S2_BACKEND_REGISTRY.register("gcs")
def _build_gcs_s2_backend(config: Any, auth: CredentialManager | None = None) -> Any:
    return build_s2_backend(config, auth=auth)


@S2_BACKEND_REGISTRY.register("local")
def _build_local_s2_backend(config: Any, auth: CredentialManager | None = None) -> Any:
    _ = auth
    return build_s2_backend(config, auth=None)


def resolve_s2_backend(config: Any, *, auth: CredentialManager | None = None) -> Any:
    return S2_BACKEND_REGISTRY.get(config.s2_data.backend)(config, auth)


def resolve_rt_model_for_pipeline(
    config: Any,
    auth: CredentialManager | None = None,
    *,
    sensor_config: SensorConfig | None = None,
) -> Any:
    return build_rt_model(config, auth=auth, sensor_config=sensor_config)


__all__ = [
    "PreprocessorRuntime",
    "build_preprocessor_runtime",
    "resolve_atmo_provider",
    "resolve_brdf_provider",
    "resolve_corrector",
    "resolve_grid_assembler",
    "resolve_output_writer",
    "resolve_preprocessor",
    "resolve_rt_model_for_pipeline",
    "resolve_s2_backend",
    "resolve_solver",
    "resolve_surface_prior_provider",
]
