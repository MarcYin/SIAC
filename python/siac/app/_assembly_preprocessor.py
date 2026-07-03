"""Preprocessor assembly for the M1 observation stage."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from siac.adapters.satellite import detect_sensor, get_preprocessor
from siac.adapters.satellite.observation import preprocess_to_observation_bundle
from siac.algorithms.cloud.mask import preferred_cloud_mask_band_names
from siac.errors import SensorNotSupportedError

if TYPE_CHECKING:
    from siac.domain.aoi import AOI
    from siac.domain.sensors import SensorConfig
    from siac.runtime import ObservationBundle
    from siac.workflows.pipeline import PreprocessorFn


@dataclass(frozen=True)
class PreprocessorRuntime:
    """Preprocessor callable plus its sensor configuration."""

    preprocessor: PreprocessorFn
    sensor_config: SensorConfig


def _ordered_unique_band_names(names: list[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(name for name in names if name))


def _bestpixel_predict_visible_anchor_band_names(config: Any, sensor_config: SensorConfig) -> list[str]:
    surface_prior = getattr(getattr(config, "algorithms", None), "surface_prior", None)
    if surface_prior is None:
        return []
    if str(getattr(surface_prior, "method", "")) != "bestpixel":
        return []
    if not bool(getattr(surface_prior, "bestpixel_predict_visible", False)):
        return []
    available = {band.name for band in getattr(sensor_config, "bands", ())}
    preferred = ("B8A", "B11", "B12")
    return [name for name in preferred if name in available]


def _toa_preload_band_names(config: Any, sensor_config: SensorConfig) -> tuple[str, ...] | None:
    """Select bands to eagerly load at native resolution during preprocessing."""
    if not hasattr(sensor_config, "select_nearest_band") or not hasattr(sensor_config, "bands"):
        return None
    cloud_mask = getattr(getattr(config, "algorithms", None), "cloud_mask", None)
    cloud_mode = str(getattr(cloud_mask, "mode", "auto")).lower()
    if cloud_mode == "user_callable":
        return None

    names: list[str] = []
    ref_band = sensor_config.select_nearest_band(665.0, tolerance_nm=80.0)
    if ref_band is None:
        ref_band = sensor_config.bands[0]
    names.append(ref_band.name)

    if cloud_mode == "auto":
        for color_name, (wl_min, wl_max) in {
            "green": (530.0, 590.0),
            "red": (630.0, 690.0),
            "nir": (760.0, 900.0),
        }.items():
            preferred = preferred_cloud_mask_band_names(sensor_config, color_name)
            if preferred:
                names.extend(preferred)
                continue
            names.extend(band.name for band in sensor_config.select_bands_in_range(wl_min, wl_max))

    names.extend(_bestpixel_predict_visible_anchor_band_names(config, sensor_config))
    resolved = _ordered_unique_band_names(names)
    return resolved or None


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

    cloud_mask_config = config.algorithms.cloud_mask.model_dump(exclude={"user_callable"})
    paths = getattr(config, "paths", None)
    rsrf_root = getattr(paths, "rsrf_root", None)
    # Wave 18: thread the content-addressed cloud-mask cache directory
    # from ``paths.caches.cloud`` into the preprocessor so the OmniCloudMask
    # PyTorch inference can be short-circuited on repeated runs over the
    # same TOA inputs (~20-25 s saved per cache hit on a Sentinel-2 scene).
    caches = getattr(paths, "caches", None)
    cloud_cache_dir = getattr(caches, "cloud", None) if caches is not None else None
    if cloud_cache_dir is not None:
        cloud_mask_config["cache_dir"] = cloud_cache_dir
    preprocessor_config: dict[str, Any] = {"cloud_mask": cloud_mask_config}
    if rsrf_root is not None:
        preprocessor_config["rsrf_root"] = rsrf_root
    try:
        preprocessor_obj = get_preprocessor_fn(sensor_name, config=preprocessor_config)
    except KeyError as exc:
        raise SensorNotSupportedError(
            f"Sensor {sensor_name!r} has no registered scene preprocessor. "
            "Sentinel-2 ('s2') is the only end-to-end processing sensor currently available."
        ) from exc

    if input_path is not None and hasattr(preprocessor_obj, "get_metadata"):
        preprocessor_obj.get_metadata(Path(input_path))

    sensor_config = preprocessor_obj.sensor_config
    preload_toa_bands = _toa_preload_band_names(config, sensor_config)
    if hasattr(preprocessor_obj, "config") and isinstance(preprocessor_obj.config, dict):
        if rsrf_root is not None:
            preprocessor_obj.config.setdefault("rsrf_root", rsrf_root)
        preprocessor_obj.config.setdefault("cloud_mask", cloud_mask_config)
        if preload_toa_bands is not None:
            preprocessor_obj.config.setdefault("preload_toa_bands", preload_toa_bands)

    def _preprocessor(path: Path, aoi: AOI | None = None) -> ObservationBundle:
        return preprocess_to_observation_bundle(
            preprocessor_obj,
            path,
            aoi=aoi,
            default_aoi_resolver=default_aoi_resolver,
            fallback_sensor_config=sensor_config,
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
