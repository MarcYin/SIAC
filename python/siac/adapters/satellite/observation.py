"""Observation-bundle conversion helpers for satellite preprocessors."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Protocol, cast

from siac.domain.aoi import AOI
from siac.runtime import ObservationBundle

if TYPE_CHECKING:
    from pathlib import Path

    import xarray as xr

    from siac.domain.sensors import SensorConfig


class SatellitePreprocessorLike(Protocol):
    """Minimum protocol required to adapt a preprocessor to pipeline input."""

    sensor_config: SensorConfig

    def preprocess(self, input_path: str | Path) -> dict[str, Any]: ...


def _resolve_default_aoi(toa: xr.Dataset, default_aoi_resolver: Any | None) -> AOI:
    if callable(default_aoi_resolver):
        return cast("AOI", default_aoi_resolver(toa))
    return AOI.from_raster(toa[list(toa.data_vars)[0]])


def _clip_dataarray_to_aoi(
    field: object, *, bounds: tuple[float, float, float, float], crs: str
) -> object:
    rio = getattr(field, "rio", None)
    if rio is None:
        return field
    try:
        return rio.clip_box(*bounds, crs=crs)
    except Exception:
        return field


def clip_raw_preprocessor_output(raw: dict[str, Any], *, aoi: AOI) -> dict[str, Any]:
    """Clip raw preprocessor output to an AOI without changing its contract."""
    import xarray as xr

    bounds = aoi.get_bounds()
    crs = str(aoi.crs)
    clipped = dict(raw)

    toa = raw.get("toa")
    if isinstance(toa, xr.Dataset):
        clipped_vars = {
            name: cast("xr.DataArray", _clip_dataarray_to_aoi(data, bounds=bounds, crs=crs))
            for name, data in toa.data_vars.items()
        }
        clipped["toa"] = xr.Dataset(clipped_vars, attrs=toa.attrs)

    geometry = raw.get("geometry")
    if geometry is not None:
        field_names = ("sza", "saa", "vza", "vaa")
        if all(hasattr(geometry, name) for name in field_names):
            clipped["geometry"] = geometry.__class__(
                **{
                    name: _clip_dataarray_to_aoi(getattr(geometry, name), bounds=bounds, crs=crs)
                    for name in field_names
                }
            )

    cloud_mask = raw.get("cloud_mask")
    clipped["cloud_mask"] = _clip_dataarray_to_aoi(cloud_mask, bounds=bounds, crs=crs)
    cloud_classes = raw.get("cloud_classes")
    if cloud_classes is not None:
        clipped["cloud_classes"] = _clip_dataarray_to_aoi(cloud_classes, bounds=bounds, crs=crs)
    return clipped


def observation_bundle_from_raw(
    raw: dict[str, Any],
    *,
    sensor_config: SensorConfig,
    crs: str,
    bounds: tuple[float, float, float, float],
) -> ObservationBundle:
    """Convert a raw satellite preprocessor dictionary into the pipeline contract."""
    return ObservationBundle(
        toa=raw["toa"],
        geometry=raw["geometry"],
        cloud_mask=raw["cloud_mask"],
        sensor_config=sensor_config,
        metadata=raw["metadata"],
        crs=crs,
        bounds=bounds,
    )


def raw_output_to_observation_bundle(
    raw: dict[str, Any],
    *,
    sensor_config: SensorConfig,
    aoi: AOI | None,
    default_aoi_resolver: Any | None = None,
) -> ObservationBundle:
    """Resolve AOI metadata and adapt raw preprocessor output to an ObservationBundle."""
    toa = raw["toa"]
    resolved_aoi = aoi or _resolve_default_aoi(toa, default_aoi_resolver)
    if aoi is not None:
        raw = clip_raw_preprocessor_output(raw, aoi=resolved_aoi)
    return observation_bundle_from_raw(
        raw,
        sensor_config=sensor_config,
        crs=str(resolved_aoi.crs),
        bounds=resolved_aoi.get_bounds(),
    )


def preprocess_to_observation_bundle(
    preprocessor: SatellitePreprocessorLike,
    input_path: str | Path,
    *,
    aoi: AOI | None = None,
    default_aoi_resolver: Any | None = None,
    fallback_sensor_config: SensorConfig | None = None,
) -> ObservationBundle:
    """Run a satellite preprocessor and return the explicit M1 pipeline contract."""
    if aoi is not None:
        set_subset = getattr(preprocessor, "set_spatial_subset", None)
        if callable(set_subset):
            set_subset(aoi.get_bounds(), crs=str(aoi.crs))
    try:
        raw = preprocessor.preprocess(input_path)
    finally:
        clear_subset = getattr(preprocessor, "clear_spatial_subset", None)
        if callable(clear_subset):
            clear_subset()

    active_sensor_config = getattr(preprocessor, "sensor_config", fallback_sensor_config)
    if active_sensor_config is None:
        raise ValueError("Satellite preprocessor did not expose a sensor_config.")
    return raw_output_to_observation_bundle(
        raw,
        sensor_config=active_sensor_config,
        aoi=aoi,
        default_aoi_resolver=default_aoi_resolver,
    )
