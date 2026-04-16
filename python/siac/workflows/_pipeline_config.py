"""Pipeline configuration helpers.

Pure config-inspection functions extracted from :mod:`siac.workflows.pipeline`
to keep the orchestrator file focused on execution flow.  Everything here
is stateless: no observability, no futures, no I/O.  Safe to import from any
layer (including tests) without side-effects.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from siac.domain.sensors import SensorConfig

# Module-level defaults used by config accessors.
_MAX_SCATTER_POINTS_PER_BAND = 4096
_DEFAULT_AUX_RESOLUTION_M = 500.0

_EXECUTION_KEYS = (
    "backend",
    "max_workers",
    "retries",
    "stage_timeout_s",
    "dashboard",
    "dashboard_address",
    "performance_report_path",
    "show_progress",
    "profiling_sample_interval_s",
    "progress_heartbeat_s",
)


def _config_scatter_limit(config: Any) -> int:
    """Read max_scatter_points_per_band from config, falling back to module default."""
    execution = getattr(config, "execution", None)
    return int(getattr(execution, "max_scatter_points_per_band", _MAX_SCATTER_POINTS_PER_BAND))


def _config_aux_resolution(config: Any) -> float:
    """Read default_aux_resolution_m from config, falling back to module default."""
    execution = getattr(config, "execution", None)
    return float(getattr(execution, "default_aux_resolution_m", _DEFAULT_AUX_RESOLUTION_M))


def _execution_values(source: Any) -> dict[str, Any]:
    """Extract execution keys from dict/model-like objects."""
    if source is None:
        return {}
    if isinstance(source, dict):
        return {k: source[k] for k in _EXECUTION_KEYS if k in source}

    out: dict[str, Any] = {}
    for key in _EXECUTION_KEYS:
        if hasattr(source, key):
            out[key] = getattr(source, key)
    return out


def _resolve_execution_settings(
    config: Any,
    *,
    execution: Any | None,
    max_workers: int | None,
) -> dict[str, Any]:
    """Resolve execution settings from config + call overrides."""
    settings: dict[str, Any] = {
        "backend": "thread",
        "max_workers": 4,
        "retries": 2,
        "stage_timeout_s": None,
        "dashboard": False,
        "dashboard_address": None,
        "performance_report_path": None,
        "show_progress": False,
        "profiling_sample_interval_s": 5.0,
        "progress_heartbeat_s": 30.0,
    }

    settings.update(_execution_values(getattr(config, "execution", None)))
    settings.update(_execution_values(execution))
    if max_workers is not None:
        settings["max_workers"] = max_workers

    backend = str(settings.get("backend", "thread")).lower()
    if backend not in {"thread", "dask"}:
        raise ValueError(f"Unsupported execution backend: {backend!r}")
    settings["backend"] = backend

    workers = int(settings.get("max_workers", 4))
    if workers < 1:
        raise ValueError("max_workers must be >= 1")
    settings["max_workers"] = workers

    retries = int(settings.get("retries", 2))
    if retries < 0:
        raise ValueError("retries must be >= 0")
    settings["retries"] = retries

    timeout = settings.get("stage_timeout_s")
    if timeout is not None:
        timeout = float(timeout)
        if timeout <= 0:
            raise ValueError("stage_timeout_s must be > 0 when provided")
    settings["stage_timeout_s"] = timeout

    report_path = settings.get("performance_report_path")
    if report_path is not None:
        report_path = Path(report_path)
    settings["performance_report_path"] = report_path

    settings["dashboard"] = bool(settings.get("dashboard", False))
    settings["show_progress"] = bool(settings.get("show_progress", False))

    sample_interval = float(settings.get("profiling_sample_interval_s", 5.0))
    if sample_interval <= 0:
        raise ValueError("profiling_sample_interval_s must be > 0")
    settings["profiling_sample_interval_s"] = sample_interval

    heartbeat = float(settings.get("progress_heartbeat_s", 30.0))
    if heartbeat <= 0:
        raise ValueError("progress_heartbeat_s must be > 0")
    settings["progress_heartbeat_s"] = heartbeat
    return settings


def _aerosol_resolution(config: Any) -> float:
    solver_config = getattr(config, "solver", None)
    if solver_config is not None:
        return float(getattr(solver_config, "aerosol_resolution", 120.0))
    return 120.0


def _requested_solver_band_names(config: Any) -> tuple[str, ...] | None:
    """Return the union of stage-requested solver bands, if any."""
    solver_config = getattr(config, "solver", None)
    if solver_config is None:
        return None

    stages = tuple(getattr(solver_config, "stages", ()) or ())
    if not stages:
        return None

    requested: list[str] = []
    for stage in stages:
        raw_bands = getattr(stage, "bands", None)
        if raw_bands is None:
            continue
        bands = (raw_bands,) if isinstance(raw_bands, str) else tuple(raw_bands)
        for name in bands:
            normalized = str(name).strip()
            if normalized and normalized not in requested:
                requested.append(normalized)

    return tuple(requested) or None


def _select_solver_bands_for_preload(
    sensor_config: SensorConfig,
    requested_band_names: tuple[str, ...] | None = None,
) -> list[Any]:
    """Mirror M4 band-selection logic for LUT preloading hints."""
    default_selector = getattr(sensor_config, "default_aerosol_solver_bands", None)
    if callable(default_selector):
        default_bands = list(default_selector())
    else:
        range_selector = getattr(sensor_config, "select_bands_in_range", None)
        if callable(range_selector):
            bands = list(range_selector(400.0, 520.0))
            default_bands = bands or list(getattr(sensor_config, "bands", ())[:2])
        else:
            default_bands = list(getattr(sensor_config, "bands", ())[:2])

    if not requested_band_names:
        return default_bands

    available_by_name = {band.name: band for band in sensor_config.bands}
    missing = [name for name in requested_band_names if name not in available_by_name]
    if missing:
        raise ValueError(
            "Requested solver band(s) are not available for sensor "
            f"{getattr(sensor_config, 'sensor_id', 'unknown')}: {', '.join(missing)}"
        )

    selected_names = {band.name for band in default_bands}
    selected_names.update(requested_band_names)
    return [band for band in sensor_config.bands if band.name in selected_names]


def _should_skip_correction(config: Any) -> bool:
    output = getattr(config, "output", None)
    defaults = getattr(output, "defaults", None)
    if defaults is None:
        return False
    return bool(getattr(defaults, "skip_correction", False))


def _should_capture_aot_scatter(config: Any) -> bool:
    output = getattr(config, "output", None)
    defaults = getattr(output, "defaults", None)
    if defaults is None:
        return True
    return bool(getattr(defaults, "include_auxiliary", True))
