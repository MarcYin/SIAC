"""Surface-prior assembly helpers shared by runtime builders."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from siac.domain.sensors import SensorConfig
    from siac.runtime import AtmosphericState, ObservationBundle, SurfacePrior
    from siac.workflows.pipeline import SurfacePriorFn


def _select_surface_prior_bands(sensor_config: SensorConfig | None) -> list[Any]:
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
    preferred = [band for name in preferred_names for band in sensor_config.bands if band.name == name]
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


def _mark_surface_prior_metadata(provider: SurfacePriorFn, *, requires_atmo_prior: bool) -> SurfacePriorFn:
    provider.requires_atmo_prior = requires_atmo_prior  # type: ignore[attr-defined]
    return provider


def _surface_prior_mapping_state(
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


def _surface_prior_brdf_request(
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
class _MonthlySurfacePriorRuntime:
    visible_bands: list[Any]
    query_bands: list[Any]
    geometry: Any
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
) -> float:
    requested = float(requested_resolution)
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


def _prepare_monthly_surface_prior_runtime(
    config: Any,
    monthly_composite_provider: Any,
    *,
    observation: ObservationBundle,
    resolution: float,
    build_database_fn: Any | None = None,
    resample_geometry_fn: Any | None = None,
) -> _MonthlySurfacePriorRuntime:
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
    monthly_filter = getattr(getattr(config, "surface_prior", None), "monthly_database_filter", None)
    filter_enabled = bool(getattr(monthly_filter, "enabled", True))
    max_prediction_uncertainty = (
        float(getattr(monthly_filter, "max_prediction_uncertainty", 0.05))
        if filter_enabled
        else None
    )
    max_composite_quality = (
        float(getattr(monthly_filter, "max_composite_quality", 0.05))
        if filter_enabled
        else None
    )
    max_source_fit_rmse = (
        float(getattr(monthly_filter, "max_source_fit_rmse", 0.05))
        if filter_enabled
        else None
    )
    max_knn_feature_distance = (
        float(getattr(monthly_filter, "max_knn_feature_distance", 0.05))
        if filter_enabled
        else None
    )
    database_resolution = _resolve_monthly_database_resolution(
        monthly_composite_provider,
        resolution,
    )
    geometry = resample_geometry_fn(observation, resolution=database_resolution)
    get_monthly_composites = getattr(monthly_composite_provider, "get_monthly_composites", None)
    if callable(get_monthly_composites):
        monthly_composites = get_monthly_composites(observation, database_resolution)
        spectral_library, spectral_k_neighbors = _surface_prior_mapping_state(
            config,
            source_bands=tuple(monthly_composites.source_bands),
            target_bands=[*visible_bands, *query_bands],
            context="Route-B monthly-database surface priors",
        )
        database = build_database_fn(
            monthly_composites=monthly_composites,
            geometry=geometry,
            visible_bands=visible_bands,
            query_bands=query_bands,
            spectral_library=spectral_library,
            spectral_k_neighbors=spectral_k_neighbors,
            max_source_fit_rmse=max_source_fit_rmse,
        )
    else:
        source_bands = tuple(getattr(monthly_composite_provider, "source_bands", [*visible_bands, *query_bands]))
        spectral_library, spectral_k_neighbors = _surface_prior_mapping_state(
            config,
            source_bands=source_bands,
            target_bands=[*visible_bands, *query_bands],
            context="Route-B monthly-database surface priors",
        )
        database = build_database_fn(
            observation=observation,
            brdf_provider=monthly_composite_provider,
            resolution=database_resolution,
            geometry=geometry,
            visible_bands=visible_bands,
            query_bands=query_bands,
            spectral_library=spectral_library,
            spectral_k_neighbors=spectral_k_neighbors,
            max_source_fit_rmse=max_source_fit_rmse,
        )
    return _MonthlySurfacePriorRuntime(
        visible_bands=visible_bands,
        query_bands=query_bands,
        geometry=geometry,
        spectral_library=spectral_library,
        spectral_k_neighbors=spectral_k_neighbors,
        max_prediction_uncertainty=max_prediction_uncertainty,
        max_composite_quality=max_composite_quality,
        max_source_fit_rmse=max_source_fit_rmse,
        max_knn_feature_distance=max_knn_feature_distance,
        database=database,
    )


def _query_monthly_surface_prior(
    observation: ObservationBundle,
    atmo_prior: AtmosphericState,
    rt_model: Any,
    runtime: _MonthlySurfacePriorRuntime,
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


__all__ = [
    "_MonthlySurfacePriorRuntime",
    "_mark_surface_prior_metadata",
    "_prepare_monthly_surface_prior_runtime",
    "_query_monthly_surface_prior",
    "_select_route_b_query_bands",
    "_select_surface_prior_bands",
    "_select_visible_surface_prior_bands",
    "_surface_prior_brdf_request",
    "_surface_prior_mapping_state",
    "_surface_spectral_mapping_runtime",
]
