"""Sentinel-2 specific workflows."""

from __future__ import annotations

from datetime import date, datetime
from pathlib import Path

from siac.adapters.auth import CredentialManager
from siac.app.planning import resolve_run_config
from siac.errors import DataNotFoundError
from siac.workflows.scene import process_scene


def coerce_date(value: date | str | None) -> date | None:
    if value is None:
        return None
    if isinstance(value, datetime):
        return value.date()
    if isinstance(value, date):
        return value
    return datetime.strptime(str(value), "%Y-%m-%d").date()


def apply_s2_query_defaults(query, *, config):
    if query.max_cloud_cover == 100.0:
        query.max_cloud_cover = config.s2_data.max_cloud_cover
    if query.processing_level == "L1C" and config.s2_data.processing_level != "L1C":
        query.processing_level = config.s2_data.processing_level
    return query


def coerce_s2_query(query, *, config):
    from siac.adapters.data.s2_data_source import S2Query

    if isinstance(query, S2Query):
        q = S2Query(**query.__dict__)
    else:
        raw = str(query).strip()
        try:
            q = S2Query.from_tile_date(raw)
        except ValueError:
            q = S2Query.from_product_id(raw)
            if "MSIL2A" in raw:
                q.processing_level = "L2A"
            elif "MSIL1C" in raw:
                q.processing_level = "L1C"
    return apply_s2_query_defaults(q, config=config)


def resolve_s2_input(
    query,
    config,
    *,
    auth: CredentialManager | None = None,
    resolve_s2_backend_fn=None,
) -> Path:
    local_candidate = Path(query).expanduser() if isinstance(query, Path) else Path(str(query)).expanduser()
    if local_candidate.exists():
        return local_candidate

    if resolve_s2_backend_fn is None:
        from siac.app.assembly import resolve_s2_backend as resolve_s2_backend_fn

    resolved_config = resolve_run_config(
        config,
        sensor=config.sensor if config.sensor != "auto" else "s2",
        s2_query=query,
    )
    auth_obj = auth or CredentialManager.from_config(resolved_config)
    backend = resolve_s2_backend_fn(resolved_config, auth=auth_obj)
    if backend is None:
        raise DataNotFoundError(
            "S2 backend is 'local', but input path does not exist. "
            "Provide a local SAFE path or switch config.s2_data.backend to 'cdse' or 'gcs'."
        )

    from siac.adapters.data.s2_data_source import S2DataAccess

    cache_dir = resolved_config.s2_data.cache_dir
    accessor = S2DataAccess(backend=backend, cache_dir=cache_dir)
    q = coerce_s2_query(query, config=resolved_config)
    return accessor.get(q, dest_dir=cache_dir)


def search_sentinel2(
    *,
    tile: str | None = None,
    date: date | str | None = None,
    date_value: date | str | None = None,
    start_date: date | str | None = None,
    end_date: date | str | None = None,
    bbox: tuple[float, float, float, float] | None = None,
    max_cloud_cover: float = 80.0,
    backend: str = "cdse",
    config=None,
    auth: CredentialManager | None = None,
    resolve_s2_backend_fn=None,
):
    from siac.adapters.data.s2_data_source import S2Query, search_s2

    cfg = config or __import__("siac.config", fromlist=["SIACConfig"]).SIACConfig(sensor="s2")
    cfg = cfg.with_overrides(
        sensor="s2",
        providers={"s2": {"backend": backend, "max_cloud_cover": max_cloud_cover}},
    )
    resolved_config = resolve_run_config(cfg, sensor="s2")
    auth_obj = auth or CredentialManager.from_config(resolved_config)

    if resolve_s2_backend_fn is None:
        from siac.app.assembly import resolve_s2_backend as resolve_s2_backend_fn

    backend_obj = resolve_s2_backend_fn(resolved_config, auth=auth_obj)
    if backend_obj is None:
        raise ValueError("search_sentinel2 does not support backend='local'.")

    query = S2Query(
        mgrs_tile=tile,
        date=coerce_date(date if date is not None else date_value),
        start_date=coerce_date(start_date),
        end_date=coerce_date(end_date),
        bbox=bbox,
        max_cloud_cover=max_cloud_cover,
        processing_level=resolved_config.s2_data.processing_level,
    )
    return search_s2(backend_obj, query)


def process_s2(
    config,
    query,
    *,
    output_path: str | Path | None = None,
    aoi=None,
    auth: CredentialManager | None = None,
    resolve_s2_input_fn=None,
):
    resolved_config = resolve_run_config(
        config,
        sensor=config.sensor if config.sensor != "auto" else "s2",
        aoi=aoi if aoi is not None else config.aoi,
        s2_query=query,
    )
    auth_obj = auth or CredentialManager.from_config(resolved_config)
    if resolve_s2_input_fn is None:
        resolve_s2_input_fn = resolve_s2_input
    input_path = resolve_s2_input_fn(query, config, auth=auth_obj)
    return process_scene(
        config,
        input_path,
        output_path=output_path,
        aoi=aoi,
        auth=auth_obj,
    )


__all__ = [
    "apply_s2_query_defaults",
    "coerce_date",
    "coerce_s2_query",
    "process_s2",
    "resolve_s2_input",
    "search_sentinel2",
]
