"""Surface-prior assembly helpers shared by runtime builders."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from siac.domain.sensors import SensorConfig
    from siac.runtime import AtmosphericState, ObservationBundle, SurfacePrior
    from siac.workflows.pipeline import SurfacePriorFn


def _select_surface_prior_bands(sensor_config: SensorConfig | None) -> list[Any]:
    if sensor_config is None:
        return list(range(1, 8))
    selected: list[Any] = list(sensor_config.select_bands_in_range(400.0, 520.0))
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
    database: Any


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
    geometry = resample_geometry_fn(observation, resolution=resolution)
    get_monthly_composites = getattr(monthly_composite_provider, "get_monthly_composites", None)
    if callable(get_monthly_composites):
        monthly_composites = get_monthly_composites(observation, resolution)
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
            resolution=resolution,
            geometry=geometry,
            visible_bands=visible_bands,
            query_bands=query_bands,
            spectral_library=spectral_library,
            spectral_k_neighbors=spectral_k_neighbors,
        )
    return _MonthlySurfacePriorRuntime(
        visible_bands=visible_bands,
        query_bands=query_bands,
        geometry=geometry,
        spectral_library=spectral_library,
        spectral_k_neighbors=spectral_k_neighbors,
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
