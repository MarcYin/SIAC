"""Registry-backed provider assembly helpers."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from siac.adapters.atmo import CAMSProvider
from siac.adapters.earthdata import earthaccess_source_from_auth
from siac.app.registry import (
    ATMO_PROVIDER_REGISTRY,
    BRDF_PROVIDER_REGISTRY,
    MONTHLY_COMPOSITE_PROVIDER_REGISTRY,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from siac.adapters.auth import CredentialManager
    from siac.algorithms.surface.brdf_monthly_composite import MonthlyCompositeCollection
    from siac.domain.protocols import AtmosphericPriorProvider, MonthlyCompositeProvider
    from siac.runtime import ObservationBundle
    from siac.workflows.pipeline import AtmoPriorFn


def _build_registered_component(
    registry: Any,
    name: str,
    config: Any,
    auth: CredentialManager | None = None,
) -> Any:
    """Resolve and invoke a registry-backed factory with the common assembly signature."""
    return registry.get(name)(config, auth)


@ATMO_PROVIDER_REGISTRY.register("cams")
def _build_cams_provider(
    config: Any,
    auth: CredentialManager | None = None,
) -> AtmosphericPriorProvider:
    atmo_config = config.providers.atmo
    return cast(
        "AtmosphericPriorProvider",
        CAMSProvider(
            atmo_config.data_path,
            temporal_interp=atmo_config.temporal_interpolation == "linear",
            download_missing=atmo_config.download_missing,
            auth=auth,
            cache_dir=atmo_config.cache_dir,
        ),
    )


@ATMO_PROVIDER_REGISTRY.register("merra2")
def _build_merra2_provider(
    config: Any,
    auth: CredentialManager | None = None,
) -> AtmosphericPriorProvider:
    from siac.adapters.atmo.merra2 import MERRA2Provider

    return cast(
        "AtmosphericPriorProvider",
        MERRA2Provider(
            cache_dir=config.providers.atmo.cache_dir,
            source=earthaccess_source_from_auth(auth),
        ),
    )


@ATMO_PROVIDER_REGISTRY.register("mcd19")
def _build_mcd19_provider(
    config: Any,
    auth: CredentialManager | None = None,
) -> AtmosphericPriorProvider:
    from siac.adapters.atmo.mcd19_earthaccess import MCD19AODProvider

    return cast(
        "AtmosphericPriorProvider",
        MCD19AODProvider(
            cache_dir=config.providers.atmo.cache_dir,
            source=earthaccess_source_from_auth(auth),
            best_quality_qa=getattr(config.providers.atmo, "maiac_best_quality_qa", False),
            allow_default_prior=getattr(config.providers.atmo, "allow_default_prior", False),
        ),
    )


@ATMO_PROVIDER_REGISTRY.register("vnp19")
def _build_vnp19_provider(
    config: Any,
    auth: CredentialManager | None = None,
) -> AtmosphericPriorProvider:
    from siac.adapters.atmo.mcd19_earthaccess import VNP19AODProvider

    return cast(
        "AtmosphericPriorProvider",
        VNP19AODProvider(
            cache_dir=config.providers.atmo.cache_dir,
            source=earthaccess_source_from_auth(auth),
            best_quality_qa=getattr(config.providers.atmo, "maiac_best_quality_qa", False),
            allow_default_prior=getattr(config.providers.atmo, "allow_default_prior", False),
        ),
    )


def _resolve_reproject_cache_dir(config: Any) -> Any:
    """Pull ``paths.caches.reproject`` through to the BRDF providers (wave 19b).

    Returns None when the cache is unconfigured; the providers default to
    no cache (legacy behaviour).
    """
    paths = getattr(config, "paths", None)
    caches = getattr(paths, "caches", None) if paths is not None else None
    return getattr(caches, "reproject", None) if caches is not None else None


@BRDF_PROVIDER_REGISTRY.register("mcd43")
def _build_mcd43_provider(config: Any, auth: CredentialManager | None = None) -> Any:
    from siac.adapters.brdf.mcd43_earthaccess import MCD43EarthAccessProvider

    return MCD43EarthAccessProvider(
        cache_dir=config.providers.brdf.cache_dir,
        source=earthaccess_source_from_auth(auth),
        reproject_cache_dir=_resolve_reproject_cache_dir(config),
    )


@BRDF_PROVIDER_REGISTRY.register("mcd43_gee")
def _build_mcd43_gee_provider(config: Any, auth: CredentialManager | None = None) -> Any:
    from siac.adapters.brdf.mcd43_gee import MCD43GEEProvider

    del auth  # GEE auth is handled by edown via ~/.config/earthengine, not SIAC.
    brdf = config.providers.brdf
    cache_dir = brdf.cache_dir or (Path.home() / ".cache" / "siac" / "mcd43_gee")
    kwargs: dict[str, Any] = {"cache_dir": cache_dir}
    # Optional override of the edown executable via providers.brdf.data_path.
    if brdf.data_path is not None:
        kwargs["edown_executable"] = str(brdf.data_path)
    return MCD43GEEProvider(**kwargs)


@BRDF_PROVIDER_REGISTRY.register("vnp43")
def _build_vnp43_provider(config: Any, auth: CredentialManager | None = None) -> Any:
    from siac.adapters.brdf.mcd43_earthaccess import VNP43EarthAccessProvider

    return VNP43EarthAccessProvider(
        cache_dir=config.providers.brdf.cache_dir,
        source=earthaccess_source_from_auth(auth),
        reproject_cache_dir=_resolve_reproject_cache_dir(config),
    )


@BRDF_PROVIDER_REGISTRY.register("mcd19")
def _build_mcd19_brdf_provider(config: Any, auth: CredentialManager | None = None) -> Any:
    from siac.adapters.brdf.mcd43_earthaccess import MCD19EarthAccessProvider

    return MCD19EarthAccessProvider(
        cache_dir=config.providers.brdf.cache_dir,
        source=earthaccess_source_from_auth(auth),
        reproject_cache_dir=_resolve_reproject_cache_dir(config),
    )


def resolve_brdf_provider(config: Any, *, auth: CredentialManager | None = None) -> Any:
    provider_name = config.providers.brdf.kind
    return _build_registered_component(BRDF_PROVIDER_REGISTRY, provider_name, config, auth)


@MONTHLY_COMPOSITE_PROVIDER_REGISTRY.register("generated_brdf")
def _build_generated_brdf_monthly_composite_provider(
    config: Any,
    auth: CredentialManager | None = None,
) -> Any:
    from siac.algorithms.surface.swir_refine import generate_monthly_composites_from_brdf

    brdf_provider = resolve_brdf_provider(config, auth=auth)

    class _GeneratedBRDFMonthlyCompositeProvider:
        resolution_m: float | None = None

        @property
        def source_name(self) -> str:
            return f"{brdf_provider.source_name} monthly composites"

        @property
        def source_bands(self) -> Any:
            return brdf_provider.source_bands

        def get_monthly_composites(
            self,
            observation: ObservationBundle,
            resolution: float,
        ) -> MonthlyCompositeCollection:
            return generate_monthly_composites_from_brdf(
                observation=observation,
                brdf_provider=brdf_provider,
                resolution=resolution,
            )

    return _GeneratedBRDFMonthlyCompositeProvider()


@MONTHLY_COMPOSITE_PROVIDER_REGISTRY.register("user_callable")
def _build_user_callable_monthly_composite_provider(
    config: Any,
    _auth: CredentialManager | None = None,
) -> Any:
    provider = getattr(config.providers.monthly_composites, "user_callable", None)
    if provider is None:
        raise ValueError(
            "monthly_composites.user_callable must be provided when kind='user_callable'"
        )
    if hasattr(provider, "get_monthly_composites"):
        return provider
    if not callable(provider):
        raise TypeError(
            "monthly_composites.user_callable must be callable or define get_monthly_composites(...)"
        )
    provider_fn = cast("Callable[[ObservationBundle, float], MonthlyCompositeCollection]", provider)

    class _CallableMonthlyCompositeProvider:
        source_name = "user_callable"
        source_bands: tuple[Any, ...] = ()
        resolution_m = getattr(provider, "resolution_m", None)

        def get_monthly_composites(
            self,
            observation: ObservationBundle,
            resolution: float,
        ) -> MonthlyCompositeCollection:
            result = provider_fn(observation, resolution)
            source_bands = getattr(result, "source_bands", None)
            if source_bands is not None:
                self.source_bands = tuple(source_bands)
            return result

    return _CallableMonthlyCompositeProvider()


@MONTHLY_COMPOSITE_PROVIDER_REGISTRY.register("prepared_store")
def _build_prepared_store_monthly_composite_provider(
    config: Any,
    _auth: CredentialManager | None = None,
) -> Any:
    from siac.algorithms.surface.monthly_composite_store import (
        MonthlyCompositeStoreManifest,
        read_monthly_composite_collection,
        read_monthly_composite_store_manifest,
    )

    store_path = getattr(config.providers.monthly_composites, "store_path", None)
    if store_path is None:
        raise ValueError(
            "monthly_composites.store_path must be provided when kind='prepared_store'"
        )
    resolved_store_path = cast("str | Path", store_path)
    strict_coverage = bool(getattr(config.providers.monthly_composites, "strict_coverage", True))

    def _validate_store_grid(
        manifest: MonthlyCompositeStoreManifest,
        observation: ObservationBundle,
        resolution: float,
    ) -> None:
        grid = manifest.grid
        if grid is None:
            raise ValueError(
                f"Prepared monthly composite store {store_path} is missing grid metadata; "
                "set monthly_composites.strict_coverage=false to bypass grid checks."
            )
        _assert_matching_store_grid(
            grid,
            observation=observation,
            resolution=resolution,
            store_path=resolved_store_path,
        )

    class _PreparedStoreMonthlyCompositeProvider:
        def __init__(self) -> None:
            self._manifest: MonthlyCompositeStoreManifest | None = None
            self._collection: MonthlyCompositeCollection | None = None

        def _load_manifest(self) -> MonthlyCompositeStoreManifest:
            if self._manifest is None:
                self._manifest = read_monthly_composite_store_manifest(resolved_store_path)
            return self._manifest

        def _load_collection(self) -> MonthlyCompositeCollection:
            if self._collection is None:
                self._collection = read_monthly_composite_collection(resolved_store_path)
            return self._collection

        @property
        def source_name(self) -> str:
            manifest = self._load_manifest()
            return manifest.source_name or f"prepared monthly composites: {resolved_store_path}"

        @property
        def source_bands(self) -> Any:
            return self._load_manifest().source_bands

        @property
        def resolution_m(self) -> float | None:
            grid = self._load_manifest().grid
            if grid is None:
                return None
            return float(grid.resolution)

        def get_monthly_composites(
            self,
            observation: ObservationBundle,
            resolution: float,
        ) -> MonthlyCompositeCollection:
            manifest = self._load_manifest()
            if strict_coverage:
                _validate_store_grid(manifest, observation, float(resolution))
            return self._load_collection()

    return _PreparedStoreMonthlyCompositeProvider()


@MONTHLY_COMPOSITE_PROVIDER_REGISTRY.register("bestpixel")
def _build_bestpixel_monthly_composite_provider(
    config: Any,
    _auth: CredentialManager | None = None,
) -> Any:
    from siac.adapters.bestpixel import BestPixelMonthlyCompositeProvider

    mc = config.providers.monthly_composites
    disk_cache = getattr(mc, "bestpixel_disk_cache", None)
    return BestPixelMonthlyCompositeProvider(
        endpoint=mc.bestpixel_endpoint,
        bands=mc.bestpixel_bands,
        lookback_years=mc.bestpixel_lookback_years,
        months=mc.bestpixel_months,
        output_crs=mc.bestpixel_output_crs,
        top_k=mc.bestpixel_top_k,
        max_cloud_cover=mc.bestpixel_max_cloud_cover,
        resolution_m=mc.bestpixel_resolution_m,
        fetch_resolution_m=mc.bestpixel_fetch_resolution_m,
        disk_cache=str(disk_cache) if disk_cache is not None else None,
    )


def _assert_matching_store_grid(
    grid: Any,
    *,
    observation: ObservationBundle,
    resolution: float,
    store_path: Any,
) -> None:
    from siac.algorithms.surface.monthly_composite_store import MonthlyCompositeStoreGridSpec

    expected = MonthlyCompositeStoreGridSpec.from_bounds(
        observation.bounds,
        crs=str(observation.crs),
        resolution=float(resolution),
    )
    if str(grid.crs) != expected.crs:
        raise ValueError(
            f"Prepared monthly composite store {store_path} uses CRS {grid.crs!r}, "
            f"but the observation requires {expected.crs!r}."
        )
    tolerance = max(1e-6, expected.resolution * 1e-6)
    grid_resolution = float(grid.resolution)
    expected_resolution = float(expected.resolution)
    if grid_resolution > (expected_resolution + tolerance):
        raise ValueError(
            f"Prepared monthly composite store {store_path} uses resolution {grid.resolution}, "
            f"but the observation requires {expected.resolution} or finer."
        )
    resolution_ratio = expected_resolution / grid_resolution
    if resolution_ratio < (1.0 - 1e-6) or not np.isclose(
        resolution_ratio, round(resolution_ratio), rtol=0.0, atol=1e-6
    ):
        raise ValueError(
            f"Prepared monthly composite store {store_path} uses resolution {grid.resolution}, "
            f"which is not an integer finer subdivision of the observation resolution {expected.resolution}."
        )
    if not np.allclose(
        np.asarray(grid.bounds, dtype=np.float64),
        np.asarray(expected.bounds, dtype=np.float64),
        rtol=0.0,
        atol=tolerance,
    ):
        raise ValueError(
            f"Prepared monthly composite store {store_path} covers bounds {tuple(grid.bounds)!r}, "
            f"but the observation requires {expected.bounds!r}."
        )


def resolve_monthly_composite_provider(
    config: Any,
    *,
    auth: CredentialManager | None = None,
) -> MonthlyCompositeProvider:
    provider_name = config.providers.monthly_composites.kind
    return cast(
        "MonthlyCompositeProvider",
        _build_registered_component(
            MONTHLY_COMPOSITE_PROVIDER_REGISTRY, provider_name, config, auth
        ),
    )


def resolve_atmo_provider(config: Any, auth: CredentialManager | None = None) -> AtmoPriorFn:
    atmo_config = config.providers.atmo
    provider = cast(
        "AtmosphericPriorProvider",
        _build_registered_component(ATMO_PROVIDER_REGISTRY, atmo_config.kind, config, auth),
    )

    fuse_with = tuple(getattr(atmo_config, "fuse_aod_with", ()) or ())
    if fuse_with:
        from siac.adapters.atmo.fusion import FusedAODProvider

        sources = [
            cast(
                "AtmosphericPriorProvider",
                _build_registered_component(ATMO_PROVIDER_REGISTRY, kind, config, auth),
            )
            for kind in fuse_with
        ]
        provider = FusedAODProvider(
            provider, sources, op=getattr(atmo_config, "fuse_aod_op", "max")
        )

    prepared_path = getattr(atmo_config, "prepared_scalar_path", None)
    if prepared_path:
        from siac.adapters.atmo.prepared import PreparedScalarAODProvider

        # Applied last: a prepared AOD is authoritative over anything derived
        # live, including a fusion, since it was validated when it was prepared.
        provider = PreparedScalarAODProvider(
            provider,
            prepared_path,
            scene_key=getattr(atmo_config, "prepared_scalar_scene_key", None),
            required=bool(getattr(atmo_config, "prepared_scalar_required", True)),
            tcwv_cm=getattr(atmo_config, "prepared_scalar_tcwv_cm", None),
        )

    return provider.get_prior
