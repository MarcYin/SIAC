"""Surface-prior assembly helpers shared by runtime builders."""

from __future__ import annotations

import logging
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from siac.algorithms.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.algorithms.surface.kernel_model import KernelModelDeriver
from siac.app._assembly_providers import (
    resolve_brdf_provider,
    resolve_monthly_composite_provider,
)
from siac.app.registry import SURFACE_PRIOR_METHOD_REGISTRY
from siac.observability import current_execution_observer, observe_stage

if TYPE_CHECKING:
    from siac.adapters.auth import CredentialManager
    from siac.domain.sensors import SensorConfig
    from siac.runtime import AtmosphericState, ObservationBundle, SurfacePrior
    from siac.workflows.pipeline import SurfacePriorFn

logger = logging.getLogger(__name__)


def select_surface_prior_bands(sensor_config: SensorConfig | None) -> list[Any]:
    if sensor_config is None:
        return list(range(1, 8))
    selected: list[Any] = list(sensor_config.default_aerosol_solver_bands())
    if selected:
        return selected
    selected = list(sensor_config.select_bands_in_range(400.0, 520.0))
    if selected:
        return selected
    return list(sensor_config.bands[:2])


def _select_visible_surface_prior_bands(sensor_config: SensorConfig) -> list[Any]:
    preferred_by_sensor = {
        "MSI": ("B01", "B02", "B04"),
        "OLI": ("B1", "B2", "B4"),
    }
    preferred_names = preferred_by_sensor.get(sensor_config.sensor_id, ())
    preferred = [
        band for name in preferred_names for band in sensor_config.bands if band.name == name
    ]
    if preferred:
        return preferred
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
    selected = [
        band for name in preferred_names for band in sensor_config.bands if band.name == name
    ]
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
    settings = config.algorithms.surface_prior.spectral_mapping

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
            neighbor_estimator=settings.neighbor_estimator,
            knn_backend=settings.knn_backend,
            knn_eps=settings.knn_eps,
            min_valid_bands=settings.min_valid_bands,
        ),
        int(settings.k_neighbors),
    )


def mark_surface_prior_metadata(
    provider: SurfacePriorFn, *, requires_atmo_prior: bool
) -> SurfacePriorFn:
    provider.requires_atmo_prior = requires_atmo_prior
    return provider


def surface_prior_mapping_state(
    config: Any,
    *,
    source_bands: list[Any] | tuple[Any, ...],
    target_bands: list[Any] | tuple[Any, ...],
    context: str,
) -> tuple[Any | None, int]:
    """Resolve spectral-mapping resources for a surface-prior builder."""
    return _surface_spectral_mapping_runtime(
        config,
        source_bands=source_bands,
        target_bands=target_bands,
        context=context,
    )


def surface_prior_brdf_request(
    observation: ObservationBundle,
    *,
    brdf_provider: Any,
    target_resolution: float,
    temporal_window: int,
) -> dict[str, Any]:
    """Build the shared BRDF request kwargs for surface-prior providers."""
    return {
        "bounds": observation.bounds,
        "crs": observation.crs,
        "obs_time": observation.metadata["observation_time"],
        "target_resolution": target_resolution,
        "bands": brdf_provider.source_bands,
        "temporal_window": temporal_window,
    }


@dataclass(frozen=True)
class MonthlySurfacePriorRuntime:
    visible_bands: list[Any]
    query_bands: list[Any]
    geometry: Any
    database_resolution: float
    query_resolution: float
    spectral_library: Any | None
    spectral_k_neighbors: int
    max_prediction_uncertainty: float | None
    max_composite_quality: float | None
    max_source_fit_rmse: float | None
    max_knn_feature_distance: float | None
    database: Any


def _resolve_monthly_database_resolution(
    monthly_composite_provider: Any,
    requested_resolution: float,
    *,
    policy: str = "provider_or_coarser",
) -> float:
    requested = float(requested_resolution)
    normalized_policy = str(policy).strip().lower()
    if normalized_policy == "aerosol":
        return requested
    if normalized_policy != "provider_or_coarser":
        raise ValueError(
            "surface_prior.monthly_database_resolution_policy must be "
            "'provider_or_coarser' or 'aerosol'"
        )
    provider_resolution = getattr(monthly_composite_provider, "resolution_m", None)
    if provider_resolution is None:
        return requested
    try:
        resolved = float(provider_resolution)
    except (TypeError, ValueError):
        return requested
    if not np.isfinite(resolved) or resolved <= 0.0:
        return requested
    # Keep the original monthly grid when the solver is finer, but aggregate to
    # the solver grid before Route-B database build when the solver is coarser.
    return max(requested, resolved)


def prepare_monthly_surface_prior_runtime(
    config: Any,
    monthly_composite_provider: Any,
    *,
    observation: ObservationBundle,
    resolution: float,
    build_database_fn: Any | None = None,
    resample_geometry_fn: Any | None = None,
) -> MonthlySurfacePriorRuntime:
    if build_database_fn is None or resample_geometry_fn is None:
        from siac.algorithms.surface.swir_refine import (
            build_monthly_surface_prior_database,
            resample_geometry_for_surface_prior,
        )

        if build_database_fn is None:
            build_database_fn = build_monthly_surface_prior_database
        if resample_geometry_fn is None:
            resample_geometry_fn = resample_geometry_for_surface_prior

    visible_bands = _select_visible_surface_prior_bands(observation.sensor_config)
    query_bands = _select_route_b_query_bands(observation.sensor_config)
    monthly_filter = config.algorithms.surface_prior.monthly_database_filter
    filter_enabled = bool(getattr(monthly_filter, "enabled", True))
    max_prediction_uncertainty = (
        float(getattr(monthly_filter, "max_prediction_uncertainty", 0.05))
        if filter_enabled
        else None
    )
    max_composite_quality = (
        float(getattr(monthly_filter, "max_composite_quality", 0.05)) if filter_enabled else None
    )
    max_source_fit_rmse = (
        float(getattr(monthly_filter, "max_source_fit_rmse", 0.05)) if filter_enabled else None
    )
    max_knn_feature_distance = (
        float(getattr(monthly_filter, "max_knn_feature_distance", 0.05)) if filter_enabled else None
    )
    database_resolution = _resolve_monthly_database_resolution(
        monthly_composite_provider,
        resolution,
        policy=str(config.algorithms.surface_prior.monthly_database_resolution_policy),
    )
    query_resolution = float(resolution)
    database_geometry = resample_geometry_fn(observation, resolution=database_resolution)
    if np.isclose(query_resolution, database_resolution, rtol=0.0, atol=1.0e-6):
        query_geometry = database_geometry
    else:
        query_geometry = resample_geometry_fn(observation, resolution=query_resolution)
    get_monthly_composites = getattr(monthly_composite_provider, "get_monthly_composites", None)
    if not callable(get_monthly_composites):
        raise TypeError("Monthly composite providers must implement get_monthly_composites().")
    monthly_composites = get_monthly_composites(observation, database_resolution)
    spectral_library, spectral_k_neighbors = surface_prior_mapping_state(
        config,
        source_bands=tuple(monthly_composites.source_bands),
        target_bands=[*visible_bands, *query_bands],
        context="Route-B monthly-database surface priors",
    )
    median_key = str(
        getattr(config.algorithms.surface_prior, "monthly_database_median_key", "query")
    )
    database = build_database_fn(
        monthly_composites=monthly_composites,
        geometry=database_geometry,
        visible_bands=visible_bands,
        query_bands=query_bands,
        spectral_library=spectral_library,
        spectral_k_neighbors=spectral_k_neighbors,
        max_source_fit_rmse=max_source_fit_rmse,
        median_key=median_key,
    )
    return MonthlySurfacePriorRuntime(
        visible_bands=visible_bands,
        query_bands=query_bands,
        geometry=query_geometry,
        database_resolution=database_resolution,
        query_resolution=query_resolution,
        spectral_library=spectral_library,
        spectral_k_neighbors=spectral_k_neighbors,
        max_prediction_uncertainty=max_prediction_uncertainty,
        max_composite_quality=max_composite_quality,
        max_source_fit_rmse=max_source_fit_rmse,
        max_knn_feature_distance=max_knn_feature_distance,
        database=database,
    )


def query_monthly_surface_prior(
    observation: ObservationBundle,
    atmo_prior: AtmosphericState,
    rt_model: Any,
    runtime: MonthlySurfacePriorRuntime,
    query_database_fn: Any | None = None,
) -> SurfacePrior:
    if query_database_fn is None:
        from siac.algorithms.surface.swir_refine import query_surface_prior_from_monthly_database

        query_database_fn = query_surface_prior_from_monthly_database

    prior = query_database_fn(
        observation=observation,
        atmo_prior=atmo_prior,
        rt_model=rt_model,
        database=runtime.database,
        query_band_names=tuple(band.name for band in runtime.query_bands),
        visible_band_names=tuple(band.name for band in runtime.visible_bands),
        target_resolution=runtime.query_resolution,
        k_neighbors=runtime.spectral_k_neighbors,
        max_prediction_uncertainty=runtime.max_prediction_uncertainty,
        max_composite_quality=runtime.max_composite_quality,
        max_knn_feature_distance=runtime.max_knn_feature_distance,
        diagnostic_cache_dir=getattr(runtime.spectral_library, "cache_dir", None),
    )
    composites = tuple(getattr(runtime.database, "composites", ()))
    if not composites or not hasattr(prior, "monthly_composites"):
        return prior
    return replace(prior, monthly_composites=composites)


def _build_surface_prior_component(
    method: str,
    config: Any,
    brdf_provider: Any,
) -> SurfacePriorFn:
    """Resolve and invoke a registry-backed surface-prior factory."""
    return cast(
        "SurfacePriorFn",
        SURFACE_PRIOR_METHOD_REGISTRY.get(method)(config, brdf_provider),
    )


@SURFACE_PRIOR_METHOD_REGISTRY.register("kernel_model")
def _build_kernel_surface_prior(config: Any, brdf_prov: Any) -> SurfacePriorFn:
    def _surface_prior(
        observation: ObservationBundle,
        atmo_prior: AtmosphericState | None,
        rt_model: Any,
        resolution: float,
    ) -> SurfacePrior:
        _ = (atmo_prior, rt_model)
        target_bands = select_surface_prior_bands(observation.sensor_config)
        spectral_library, spectral_k_neighbors = surface_prior_mapping_state(
            config,
            source_bands=tuple(brdf_prov.source_bands),
            target_bands=target_bands,
            context="surface priors",
        )
        brdf_weights = brdf_prov.get_brdf_parameters(
            **surface_prior_brdf_request(
                observation,
                brdf_provider=brdf_prov,
                target_resolution=resolution,
                temporal_window=config.providers.brdf.temporal_window,
            ),
        )
        deriver = KernelModelDeriver(
            psf_sigma_x=config.algorithms.surface_prior.psf_sigma_x,
            psf_sigma_y=config.algorithms.surface_prior.psf_sigma_y,
            apply_psf=config.algorithms.surface_prior.apply_psf,
        )
        return deriver.compute_surface_prior(
            brdf_weights,
            observation.geometry,
            source_bands=brdf_prov.source_bands,
            target_bands=target_bands,
            spectral_library=spectral_library,
            spectral_k_neighbors=spectral_k_neighbors,
        )

    return mark_surface_prior_metadata(_surface_prior, requires_atmo_prior=False)


@SURFACE_PRIOR_METHOD_REGISTRY.register("whittaker")
def _build_whittaker_surface_prior(config: Any, brdf_prov: Any) -> SurfacePriorFn:
    def _surface_prior(
        observation: ObservationBundle,
        atmo_prior: AtmosphericState | None,
        rt_model: Any,
        resolution: float,
    ) -> SurfacePrior:
        _ = (atmo_prior, rt_model)
        target_bands = select_surface_prior_bands(observation.sensor_config)
        spectral_library, spectral_k_neighbors = surface_prior_mapping_state(
            config,
            source_bands=tuple(brdf_prov.source_bands),
            target_bands=target_bands,
            context="Whittaker surface priors",
        )
        brdf_weights = brdf_prov.get_temporal_brdf_parameters(
            **surface_prior_brdf_request(
                observation,
                brdf_provider=brdf_prov,
                target_resolution=resolution,
                temporal_window=config.providers.brdf.temporal_window,
            ),
        )
        deriver = BRDFWhittakerDeriver(
            temporal_lambda=config.algorithms.surface_prior.whittaker_lambda,
            psf_sigma_x=config.algorithms.surface_prior.psf_sigma_x,
            psf_sigma_y=config.algorithms.surface_prior.psf_sigma_y,
            apply_psf=config.algorithms.surface_prior.apply_psf,
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

    return mark_surface_prior_metadata(_surface_prior, requires_atmo_prior=False)


# NOTE: This function is intentionally NOT registered in
# SURFACE_PRIOR_METHOD_REGISTRY. Its construction needs the
# ``fallback_brdf_provider_factory`` keyword argument, which the registry's
# uniform ``(config, brdf_provider)`` signature cannot supply. ``resolve_surface
# _prior_provider`` dispatches the ``monthly_database`` method to this builder
# inline (REVIEW.md §3.6 _assembly_surface.py).
def _build_monthly_surface_prior(
    config: Any,
    monthly_composite_provider: Any,
    *,
    fallback_brdf_provider_factory: Any | None = None,
) -> SurfacePriorFn:
    def _surface_prior(
        observation: ObservationBundle,
        atmo_prior: AtmosphericState | None,
        rt_model: Any,
        resolution: float,
    ) -> SurfacePrior:
        if atmo_prior is None:
            raise ValueError("Route-B monthly_database surface prior requires an atmospheric prior")

        logger.info("M3: Building monthly surface-prior database...")
        with observe_stage(
            "M3.monthly_database.prepare_runtime",
            details={"resolution_m": resolution},
        ):
            monthly_runtime = prepare_monthly_surface_prior_runtime(
                config,
                monthly_composite_provider,
                observation=observation,
                resolution=resolution,
            )

        observer = current_execution_observer()
        if observer is not None:
            observer.emit(
                "progress",
                stage="M3.monthly_database.prepare_runtime",
                message="Monthly surface-prior database ready.",
                query_band_count=len(monthly_runtime.query_bands),
                visible_band_count=len(monthly_runtime.visible_bands),
            )

        logger.info("M3: Querying monthly surface-prior database...")
        with observe_stage("M3.monthly_database.query"):
            prior = query_monthly_surface_prior(observation, atmo_prior, rt_model, monthly_runtime)
        if bool(np.asarray(prior.mask.values).any()):
            return prior

        logger.warning(
            "Route-B monthly database produced no valid surface-prior pixels; "
            "falling back to kernel-model BRDF priors."
        )
        if observer is not None:
            observer.emit(
                "warning",
                stage="M3.monthly_database.query",
                status="fallback",
                message="Monthly database produced no valid pixels; falling back to kernel-model priors.",
            )
        with observe_stage("M3.monthly_database.fallback"):
            if fallback_brdf_provider_factory is not None:
                brdf_prov = fallback_brdf_provider_factory()
            elif hasattr(monthly_composite_provider, "get_brdf_parameters"):
                brdf_prov = monthly_composite_provider
            else:
                raise ValueError(
                    "Route-B monthly_database fallback requires a BRDF provider factory"
                )
            brdf_weights = brdf_prov.get_brdf_parameters(
                **surface_prior_brdf_request(
                    observation,
                    brdf_provider=brdf_prov,
                    target_resolution=resolution,
                    temporal_window=config.providers.brdf.temporal_window,
                ),
            )
            fallback_deriver = KernelModelDeriver(
                psf_sigma_x=config.algorithms.surface_prior.psf_sigma_x,
                psf_sigma_y=config.algorithms.surface_prior.psf_sigma_y,
                apply_psf=config.algorithms.surface_prior.apply_psf,
            )
            return fallback_deriver.compute_surface_prior(
                brdf_weights,
                monthly_runtime.geometry,
                source_bands=brdf_prov.source_bands,
                target_bands=monthly_runtime.visible_bands,
                spectral_library=monthly_runtime.spectral_library,
                spectral_k_neighbors=monthly_runtime.spectral_k_neighbors,
            )

    return mark_surface_prior_metadata(_surface_prior, requires_atmo_prior=True)


def resolve_surface_prior_provider(
    config: Any, auth: CredentialManager | None = None
) -> SurfacePriorFn:
    method = config.algorithms.surface_prior.method
    if method == "monthly_database":
        monthly_provider = resolve_monthly_composite_provider(config, auth=auth)
        return cast(
            "SurfacePriorFn",
            _build_monthly_surface_prior(
                config,
                monthly_provider,
                fallback_brdf_provider_factory=lambda: resolve_brdf_provider(config, auth=auth),
            ),
        )

    brdf_prov = resolve_brdf_provider(config, auth=auth)
    return _build_surface_prior_component(method, config, brdf_prov)
