"""Assembly layer: resolve configured components into runtime callables.

This module remains the public compatibility facade so existing callers and
monkeypatch-based tests can keep targeting ``siac.app.assembly``.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from siac.adapters.s2_backend import build_s2_backend
from siac.algorithms.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.algorithms.surface.kernel_model import KernelModelDeriver
from siac.algorithms.surface.swir_refine import (
    build_monthly_surface_prior_database,
    query_surface_prior_from_monthly_database,
    resample_geometry_for_surface_prior,
)
from siac.app._assembly_providers import (
    _build_registered_component,
    resolve_atmo_provider,
    resolve_brdf_provider,
    resolve_monthly_composite_provider,
)
from siac.app._assembly_runtime import (
    PreprocessorRuntime,
    build_preprocessor_runtime,
    resolve_corrector,
    resolve_grid_assembler,
    resolve_output_writer,
    resolve_preprocessor,
    resolve_rt_model_for_pipeline,
    resolve_solver,
)
from siac.app._assembly_surface import (
    _mark_surface_prior_metadata,
    _MonthlySurfacePriorRuntime,
    _select_surface_prior_bands,
    _surface_prior_brdf_request,
    _surface_prior_mapping_state,
)
from siac.app._assembly_surface import (
    _prepare_monthly_surface_prior_runtime as _shared_prepare_monthly_surface_prior_runtime,
)
from siac.app._assembly_surface import (
    _query_monthly_surface_prior as _shared_query_monthly_surface_prior,
)
from siac.app.registry import S2_BACKEND_REGISTRY, SURFACE_PRIOR_METHOD_REGISTRY
from siac.observability import current_execution_observer, observe_stage

if TYPE_CHECKING:
    from siac.adapters.auth import CredentialManager
    from siac.runtime import AtmosphericState, ObservationBundle, SurfacePrior
    from siac.workflows.pipeline import SurfacePriorFn

logger = logging.getLogger(__name__)


def _build_s2_backend_with_auth(
    config: Any,
    auth: CredentialManager | None,
    *,
    use_auth: bool,
) -> Any:
    """Centralize the only behavioral branch in S2 backend assembly: auth usage."""
    return build_s2_backend(config, auth=auth if use_auth else None)


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
        target_bands = _select_surface_prior_bands(observation.sensor_config)
        spectral_library, spectral_k_neighbors = _surface_prior_mapping_state(
            config,
            source_bands=tuple(brdf_prov.source_bands),
            target_bands=target_bands,
            context="surface priors",
        )
        brdf_weights = brdf_prov.get_brdf_parameters(
            **_surface_prior_brdf_request(
                observation,
                brdf_provider=brdf_prov,
                target_resolution=resolution,
                temporal_window=config.brdf.temporal_window,
            ),
        )
        deriver = KernelModelDeriver(
            psf_sigma_x=config.surface_prior.psf_sigma_x,
            psf_sigma_y=config.surface_prior.psf_sigma_y,
            apply_psf=config.surface_prior.apply_psf,
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
    def _surface_prior(
        observation: ObservationBundle,
        atmo_prior: AtmosphericState | None,
        rt_model: Any,
        resolution: float,
    ) -> SurfacePrior:
        _ = (atmo_prior, rt_model)
        target_bands = _select_surface_prior_bands(observation.sensor_config)
        spectral_library, spectral_k_neighbors = _surface_prior_mapping_state(
            config,
            source_bands=tuple(brdf_prov.source_bands),
            target_bands=target_bands,
            context="Whittaker surface priors",
        )
        brdf_weights = brdf_prov.get_temporal_brdf_parameters(
            **_surface_prior_brdf_request(
                observation,
                brdf_provider=brdf_prov,
                target_resolution=resolution,
                temporal_window=config.brdf.temporal_window,
            ),
        )
        deriver = BRDFWhittakerDeriver(
            temporal_lambda=config.surface_prior.whittaker_lambda,
            psf_sigma_x=config.surface_prior.psf_sigma_x,
            psf_sigma_y=config.surface_prior.psf_sigma_y,
            apply_psf=config.surface_prior.apply_psf,
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


def _prepare_monthly_surface_prior_runtime(
    config: Any,
    monthly_composite_provider: Any,
    *,
    observation: ObservationBundle,
    resolution: float,
) -> _MonthlySurfacePriorRuntime:
    return _shared_prepare_monthly_surface_prior_runtime(
        config,
        monthly_composite_provider,
        observation=observation,
        resolution=resolution,
        build_database_fn=build_monthly_surface_prior_database,
        resample_geometry_fn=resample_geometry_for_surface_prior,
    )


def _query_monthly_surface_prior(
    observation: ObservationBundle,
    atmo_prior: AtmosphericState,
    rt_model: Any,
    runtime: _MonthlySurfacePriorRuntime,
) -> SurfacePrior:
    return _shared_query_monthly_surface_prior(
        observation,
        atmo_prior,
        rt_model,
        runtime,
        query_database_fn=query_surface_prior_from_monthly_database,
    )


@SURFACE_PRIOR_METHOD_REGISTRY.register("monthly_database")
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
            monthly_runtime = _prepare_monthly_surface_prior_runtime(
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
            prior = _query_monthly_surface_prior(observation, atmo_prior, rt_model, monthly_runtime)
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
                **_surface_prior_brdf_request(
                    observation,
                    brdf_provider=brdf_prov,
                    target_resolution=resolution,
                    temporal_window=config.brdf.temporal_window,
                ),
            )
            fallback_deriver = KernelModelDeriver(
                psf_sigma_x=config.surface_prior.psf_sigma_x,
                psf_sigma_y=config.surface_prior.psf_sigma_y,
                apply_psf=config.surface_prior.apply_psf,
            )
            return fallback_deriver.compute_surface_prior(
                brdf_weights,
                monthly_runtime.geometry,
                source_bands=brdf_prov.source_bands,
                target_bands=monthly_runtime.visible_bands,
                spectral_library=monthly_runtime.spectral_library,
                spectral_k_neighbors=monthly_runtime.spectral_k_neighbors,
            )

    return _mark_surface_prior_metadata(_surface_prior, requires_atmo_prior=True)


def resolve_surface_prior_provider(
    config: Any, auth: CredentialManager | None = None
) -> SurfacePriorFn:
    method = getattr(config.surface_prior, "method", "kernel_model")
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


@S2_BACKEND_REGISTRY.register("cdse")
def _build_cdse_s2_backend(config: Any, auth: CredentialManager | None = None) -> Any:
    return _build_s2_backend_with_auth(config, auth, use_auth=True)


@S2_BACKEND_REGISTRY.register("gcs")
def _build_gcs_s2_backend(config: Any, auth: CredentialManager | None = None) -> Any:
    return _build_s2_backend_with_auth(config, auth, use_auth=True)


@S2_BACKEND_REGISTRY.register("local")
def _build_local_s2_backend(config: Any, auth: CredentialManager | None = None) -> Any:
    return _build_s2_backend_with_auth(config, auth, use_auth=False)


def resolve_s2_backend(config: Any, *, auth: CredentialManager | None = None) -> Any:
    return _build_registered_component(S2_BACKEND_REGISTRY, config.s2_data.backend, config, auth)


__all__ = [
    "PreprocessorRuntime",
    "build_preprocessor_runtime",
    "resolve_atmo_provider",
    "resolve_brdf_provider",
    "resolve_monthly_composite_provider",
    "resolve_corrector",
    "resolve_grid_assembler",
    "resolve_output_writer",
    "resolve_preprocessor",
    "resolve_rt_model_for_pipeline",
    "resolve_s2_backend",
    "resolve_solver",
    "resolve_surface_prior_provider",
]
