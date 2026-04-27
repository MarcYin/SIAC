"""
Resolve system config and run requests into a fully merged runtime config.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from siac.config.schema import (
    ResolvedAlgorithmsConfig,
    ResolvedAtmoProviderConfig,
    ResolvedAuthConfig,
    ResolvedAWSAuthConfig,
    ResolvedBRDFProviderConfig,
    ResolvedCachePathsConfig,
    ResolvedCDSAuthConfig,
    ResolvedCDSEAuthConfig,
    ResolvedConfig,
    ResolvedEarthdataAuthConfig,
    ResolvedGCSAuthConfig,
    ResolvedMonthlyCompositeProviderConfig,
    ResolvedPathsConfig,
    ResolvedProvidersConfig,
    ResolvedRTAlgorithmConfig,
    ResolvedRunConfig,
    ResolvedS2ProviderConfig,
    ResolvedSpectralMappingConfig,
    ResolvedSurfacePriorAlgorithmConfig,
    RunRequest,
    SystemConfig,
)

if TYPE_CHECKING:
    from pathlib import Path


def _default_cache_path(root: Path | None, name: str) -> Path | None:
    if root is None:
        return None
    return root / name


def _resolve_cache_dir(explicit: Path | None, default: Path | None) -> Path | None:
    return explicit if explicit is not None else default


def resolve_auth(system: SystemConfig) -> ResolvedAuthConfig:
    return ResolvedAuthConfig(
        cdse=ResolvedCDSEAuthConfig(
            username=system.auth.cdse.username,
            password=system.auth.cdse.password,
        ),
        cds=ResolvedCDSAuthConfig(api_key=system.auth.cds.api_key),
        aws=ResolvedAWSAuthConfig(
            access_key_id=system.auth.aws.access_key_id,
            secret_access_key=system.auth.aws.secret_access_key,
        ),
        earthdata=ResolvedEarthdataAuthConfig(
            username=system.auth.earthdata.username,
            password=system.auth.earthdata.password,
        ),
        gcs=ResolvedGCSAuthConfig(credentials_file=system.auth.gcs.credentials_file),
    )


def resolve_paths(system: SystemConfig) -> ResolvedPathsConfig:
    cache_root = system.paths.cache_root
    defaults = ResolvedCachePathsConfig(
        atmo=_default_cache_path(cache_root, "atmo"),
        brdf=_default_cache_path(cache_root, "brdf"),
        s2=_default_cache_path(cache_root, "s2"),
    )
    return ResolvedPathsConfig(
        dem=system.paths.dem,
        water_mask=system.paths.water_mask,
        emulator_dir=system.paths.emulator_dir,
        lut_path=system.paths.lut_path,
        rsrf_root=system.paths.rsrf_root,
        cache_root=cache_root,
        caches=ResolvedCachePathsConfig(
            atmo=_resolve_cache_dir(system.paths.caches.atmo, defaults.atmo),
            brdf=_resolve_cache_dir(system.paths.caches.brdf, defaults.brdf),
            s2=_resolve_cache_dir(system.paths.caches.s2, defaults.s2),
        ),
    )


def resolve_config(system: SystemConfig, request: RunRequest) -> ResolvedConfig:
    paths = resolve_paths(system)
    auth = resolve_auth(system)
    sensor = request.sensor or getattr(system, "sensor", None) or "auto"
    aoi = request.aoi if request.aoi is not None else getattr(system, "aoi", None)

    providers = ResolvedProvidersConfig(
        atmo=ResolvedAtmoProviderConfig(
            **system.providers.atmo.model_dump(mode="python", exclude={"cache_dir"}),
            cache_dir=_resolve_cache_dir(system.providers.atmo.cache_dir, paths.caches.atmo),
        ),
        brdf=ResolvedBRDFProviderConfig(
            **system.providers.brdf.model_dump(mode="python", exclude={"cache_dir"}),
            cache_dir=_resolve_cache_dir(system.providers.brdf.cache_dir, paths.caches.brdf),
        ),
        s2=ResolvedS2ProviderConfig(
            **system.providers.s2.model_dump(mode="python", exclude={"cache_dir"}),
            cache_dir=_resolve_cache_dir(system.providers.s2.cache_dir, paths.caches.s2),
        ),
        monthly_composites=ResolvedMonthlyCompositeProviderConfig(
            **system.providers.monthly_composites.model_dump(mode="python"),
        ),
    )

    algorithms = ResolvedAlgorithmsConfig(
        surface_prior=ResolvedSurfacePriorAlgorithmConfig(
            **system.algorithms.surface_prior.model_dump(
                exclude={"spectral_mapping"}, mode="python"
            ),
            spectral_mapping=ResolvedSpectralMappingConfig(
                **system.algorithms.surface_prior.spectral_mapping.model_dump(mode="python"),
            ),
        ),
        rt=ResolvedRTAlgorithmConfig(
            **system.algorithms.rt.model_dump(mode="python"),
            emulator_dir=paths.emulator_dir,
            lut_path=paths.lut_path,
        ),
        solver=system.algorithms.solver.model_copy(deep=True),
        cloud_mask=system.algorithms.cloud_mask.model_copy(deep=True),
    )

    return ResolvedConfig(
        run=ResolvedRunConfig(
            input_path=request.input_path,
            output_path=request.output_path or system.output.defaults.output_dir,
            sensor=sensor,
            aoi=aoi,
            s2_query=request.s2_query,
        ),
        paths=paths,
        auth=auth,
        providers=providers,
        algorithms=algorithms,
        runtime=system.runtime.model_copy(deep=True),
        output=system.output.model_copy(deep=True),
    )
