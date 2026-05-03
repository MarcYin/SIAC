"""Sentinel-2 query coercion and access helpers."""

from __future__ import annotations

from datetime import date, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from siac.adapters.auth import CredentialManager
from siac.app.planning import resolve_run_config
from siac.config import SIACConfig
from siac.errors import DataNotFoundError

if TYPE_CHECKING:
    from siac.adapters.data.s2_data_source import S2Product, S2Query
    from siac.app.requests import Sentinel2ResolveRequest, Sentinel2SearchRequest


def coerce_date(value: date | datetime | str | None) -> date | None:
    if value is None:
        return None
    if isinstance(value, datetime):
        return value.date()
    if isinstance(value, date):
        return value
    return datetime.strptime(str(value), "%Y-%m-%d").date()


def apply_s2_query_defaults(query: S2Query, *, config: Any) -> S2Query:
    s2_config = config.providers.s2
    if query.max_cloud_cover == 100.0:
        query.max_cloud_cover = s2_config.max_cloud_cover
    if query.processing_level == "L1C" and s2_config.processing_level != "L1C":
        query.processing_level = s2_config.processing_level
    return query


def coerce_s2_query(query: Any, *, config: Any) -> S2Query:
    from siac.adapters.data.s2_data_source import S2Query

    if isinstance(query, S2Query):
        resolved = S2Query(**query.__dict__)
    else:
        raw = str(query).strip()
        try:
            resolved = S2Query.from_tile_date(raw)
        except ValueError:
            resolved = S2Query.from_product_id(raw)
            if "MSIL2A" in raw:
                resolved.processing_level = "L2A"
            elif "MSIL1C" in raw:
                resolved.processing_level = "L1C"
    return apply_s2_query_defaults(resolved, config=config)


def _as_local_candidate(query: Any) -> Path:
    return Path(query).expanduser() if isinstance(query, Path) else Path(str(query)).expanduser()


def resolve_s2_input(
    request: Sentinel2ResolveRequest,
    *,
    resolve_s2_backend_fn: Any | None = None,
) -> Path:
    local_candidate = _as_local_candidate(request.query)
    if local_candidate.exists():
        return local_candidate

    if resolve_s2_backend_fn is None:
        from siac.app.s2_backend import resolve_s2_backend as resolve_s2_backend_fn
    resolve_s2_backend = cast("Any", resolve_s2_backend_fn)

    resolved_config = resolve_run_config(
        request.config,
        sensor=request.config.sensor if request.config.sensor != "auto" else "s2",
        s2_query=cast("Any", request.query),
    )
    auth_obj = request.auth or CredentialManager.from_config(resolved_config)
    backend = resolve_s2_backend(resolved_config, auth=auth_obj)
    if backend is None:
        raise DataNotFoundError(
            "S2 backend is 'local', but input path does not exist. "
            "Provide a local SAFE path or switch config.providers.s2.backend to 'cdse' or 'gcs'."
        )

    from siac.adapters.data.s2_data_source import S2DataAccess

    cache_dir = resolved_config.providers.s2.cache_dir
    accessor = S2DataAccess(backend=backend, cache_dir=cache_dir)
    query = coerce_s2_query(request.query, config=resolved_config)
    return cast("Path", accessor.get(query, dest_dir=cache_dir))


def _resolve_search_config(request: Sentinel2SearchRequest) -> Any:
    config = request.config or SIACConfig(sensor="s2")
    config = config.with_overrides(
        sensor="s2",
        providers={"s2": {"backend": request.backend, "max_cloud_cover": request.max_cloud_cover}},
    )
    return resolve_run_config(config, sensor="s2")


def search_sentinel2(
    request: Sentinel2SearchRequest,
    *,
    resolve_s2_backend_fn: Any | None = None,
) -> list[S2Product]:
    from siac.adapters.data.s2_data_source import S2Query, search_s2

    resolved_config = _resolve_search_config(request)
    auth_obj = request.auth or CredentialManager.from_config(resolved_config)

    if resolve_s2_backend_fn is None:
        from siac.app.s2_backend import resolve_s2_backend as resolve_s2_backend_fn
    resolve_s2_backend = cast("Any", resolve_s2_backend_fn)

    backend_obj = resolve_s2_backend(resolved_config, auth=auth_obj)
    if backend_obj is None:
        raise ValueError("search_sentinel2 does not support backend='local'.")

    query = S2Query(
        mgrs_tile=request.tile,
        date=coerce_date(request.date if request.date is not None else request.date_value),
        start_date=coerce_date(request.start_date),
        end_date=coerce_date(request.end_date),
        bbox=request.bbox,
        max_cloud_cover=request.max_cloud_cover,
        processing_level=resolved_config.providers.s2.processing_level,
    )
    return cast("list[S2Product]", search_s2(backend_obj, query))
