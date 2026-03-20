"""Cloud/cloud-shadow mask orchestration and standardized class generation."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr

from siac.algorithms.cloud.mapping import apply_class_mapping
from siac.algorithms.cloud.providers.omnicloudmask import OmniCloudMaskProvider
from siac.geo import reproject_match
from siac.geo.reprojection import resample
from siac.storage.readers import read_raster

if TYPE_CHECKING:
    from collections.abc import Callable

    from siac.domain import SensorConfig

_EXPECTED_CLASS_VALUES = (0, 1, 2, 3)
_COLOR_WINDOWS = {
    "green": (530.0, 590.0),
    "red": (630.0, 690.0),
    "nir": (760.0, 900.0),
}


def classes_to_bool_mask(cloud_classes: xr.DataArray) -> xr.DataArray:
    """Convert standardized cloud classes to boolean cloud mask compatibility form."""
    classes = apply_class_mapping(cloud_classes, None, unmapped_to_missing=False)
    bool_mask = xr.where(classes == 1, False, True).astype(bool)
    bool_mask.name = "cloud_mask"
    return bool_mask


def _extract_band(da: xr.DataArray) -> xr.DataArray:
    if "band" in da.dims:
        if da.sizes["band"] == 1:
            return da.squeeze("band", drop=True)
        return da.isel(band=0, drop=True)
    return da


def _resample_continuous(
    da: xr.DataArray,
    *,
    target_resolution_m: float,
    resolution_policy: str,
    allow_upsample_to_target: bool,
) -> xr.DataArray:
    if target_resolution_m <= 0:
        raise ValueError("target_resolution_m must be > 0")
    if not hasattr(da, "rio") or da.rio.crs is None:
        return da

    current = abs(float(da.rio.resolution()[0]))

    should_resample = False
    if resolution_policy == "force":
        should_resample = abs(current - target_resolution_m) > 1e-6
    elif resolution_policy == "auto":
        if current < target_resolution_m - 1e-6:
            # downsample finer data to target (preferred default)
            should_resample = True
        elif current > target_resolution_m + 1e-6 and allow_upsample_to_target:
            should_resample = True
    else:
        raise ValueError("resolution_policy must be 'auto' or 'force'")

    if not should_resample:
        return da

    return resample(da, target_resolution=target_resolution_m, resampling="bilinear")


def _resample_classes(
    classes: xr.DataArray,
    *,
    target_resolution_m: float,
    resolution_policy: str,
    allow_upsample_to_target: bool,
) -> xr.DataArray:
    if target_resolution_m <= 0:
        raise ValueError("target_resolution_m must be > 0")
    if not hasattr(classes, "rio") or classes.rio.crs is None:
        return classes

    current = abs(float(classes.rio.resolution()[0]))
    should_resample = False
    if resolution_policy == "force":
        should_resample = abs(current - target_resolution_m) > 1e-6
    elif resolution_policy == "auto":
        if current < target_resolution_m - 1e-6 or current > target_resolution_m + 1e-6 and allow_upsample_to_target:
            should_resample = True
    else:
        raise ValueError("resolution_policy must be 'auto' or 'force'")

    if not should_resample:
        return classes

    return resample(classes, target_resolution=target_resolution_m, resampling="nearest")


def _group_band_names(
    toa: xr.Dataset,
    sensor_config: SensorConfig,
    color_name: str,
) -> list[str]:
    wl_min, wl_max = _COLOR_WINDOWS[color_name]
    names: list[str] = []

    for band in sensor_config.bands:
        if band.name not in toa.data_vars:
            continue
        if wl_min <= band.center_wavelength <= wl_max:
            names.append(band.name)

    if not names:
        raise ValueError(
            f"Could not find any {color_name} band in TOA data within "
            f"{wl_min:.0f}-{wl_max:.0f} nm"
        )

    return names


def _mean_group(
    toa: xr.Dataset,
    names: list[str],
    *,
    target_resolution_m: float,
    resolution_policy: str,
    allow_upsample_to_target: bool,
) -> xr.DataArray:
    arrays: list[xr.DataArray] = []
    for name in names:
        da = _extract_band(toa[name])
        da = _resample_continuous(
            da,
            target_resolution_m=target_resolution_m,
            resolution_policy=resolution_policy,
            allow_upsample_to_target=allow_upsample_to_target,
        )
        arrays.append(da)

    ref = arrays[0]
    aligned = [ref]
    for arr in arrays[1:]:
        if (
            hasattr(arr, "rio")
            and hasattr(ref, "rio")
            and arr.rio.crs is not None
            and ref.rio.crs is not None
        ):
            aligned.append(reproject_match(arr, ref, resampling="bilinear"))
        elif arr.dims == ref.dims and all(arr.sizes[d] == ref.sizes[d] for d in ref.dims):
            aligned.append(arr)
        else:
            aligned.append(arr.interp_like(ref))

    if len(aligned) == 1:
        out = aligned[0]
    else:
        out = xr.concat(aligned, dim="source_band").mean(dim="source_band")

    return out.astype(np.float32)


def _prepare_rgbnir(
    toa: xr.Dataset,
    sensor_config: SensorConfig,
    *,
    target_resolution_m: float,
    resolution_policy: str,
    allow_upsample_to_target: bool,
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    red_names = _group_band_names(toa, sensor_config, "red")
    green_names = _group_band_names(toa, sensor_config, "green")
    nir_names = _group_band_names(toa, sensor_config, "nir")

    red = _mean_group(
        toa,
        red_names,
        target_resolution_m=target_resolution_m,
        resolution_policy=resolution_policy,
        allow_upsample_to_target=allow_upsample_to_target,
    )
    green = _mean_group(
        toa,
        green_names,
        target_resolution_m=target_resolution_m,
        resolution_policy=resolution_policy,
        allow_upsample_to_target=allow_upsample_to_target,
    )
    nir = _mean_group(
        toa,
        nir_names,
        target_resolution_m=target_resolution_m,
        resolution_policy=resolution_policy,
        allow_upsample_to_target=allow_upsample_to_target,
    )

    if (
        hasattr(green, "rio")
        and hasattr(red, "rio")
        and green.rio.crs is not None
        and red.rio.crs is not None
    ):
        green = reproject_match(green, red, resampling="bilinear")
    elif not (green.dims == red.dims and all(green.sizes[d] == red.sizes[d] for d in red.dims)):
        green = green.interp_like(red)

    if (
        hasattr(nir, "rio")
        and hasattr(red, "rio")
        and nir.rio.crs is not None
        and red.rio.crs is not None
    ):
        nir = reproject_match(nir, red, resampling="bilinear")
    elif not (nir.dims == red.dims and all(nir.sizes[d] == red.sizes[d] for d in red.dims)):
        nir = nir.interp_like(red)

    return red, green, nir


def _external_classes(
    external_mask_path: Path,
    *,
    reference: xr.DataArray,
    class_mapping: dict[int, list[int]] | None,
    target_resolution_m: float,
    resolution_policy: str,
    allow_upsample_to_target: bool,
    unmapped_to_missing: bool,
) -> xr.DataArray:
    raw = read_raster(external_mask_path)
    raw = _extract_band(raw)

    # Align to reference grid first; nearest preserves class labels.
    if hasattr(raw, "rio") and raw.rio.crs is not None:
        raw = reproject_match(raw, reference, resampling="nearest")

    raw = _resample_classes(
        raw,
        target_resolution_m=target_resolution_m,
        resolution_policy=resolution_policy,
        allow_upsample_to_target=allow_upsample_to_target,
    )

    return apply_class_mapping(
        raw,
        class_mapping,
        unmapped_to_missing=unmapped_to_missing,
    )


def _call_user_callable(
    fn: Callable[..., Any],
    *,
    toa: xr.Dataset,
    sensor_config: SensorConfig,
    input_path: Path | None,
) -> Any:
    try:
        return fn(toa=toa, sensor_config=sensor_config, input_path=input_path)
    except TypeError:
        # Fallback for simpler callables.
        try:
            return fn(toa, sensor_config, input_path)
        except TypeError:
            return fn(toa)


def build_cloud_classes(
    toa: xr.Dataset,
    sensor_config: SensorConfig,
    *,
    mode: str = "auto",
    provider: str = "omnicloudmask",
    class_mapping: dict[int, list[int]] | None = None,
    external_mask_path: str | Path | None = None,
    user_callable: Callable[..., Any] | None = None,
    target_resolution_m: float = 10.0,
    resolution_policy: str = "auto",
    allow_upsample_to_target: bool = False,
    unmapped_to_missing: bool = True,
) -> xr.DataArray:
    """Build standardized cloud classes (0/1/2/3) from model/file/user inputs."""
    if len(toa.data_vars) == 0:
        raise ValueError("TOA dataset must contain at least one band")

    ref = _extract_band(next(iter(toa.data_vars.values())))

    if mode == "none":
        classes = xr.full_like(ref, 1, dtype=np.uint8)
        classes = classes.where(np.isfinite(ref), 0).astype(np.uint8)
        classes.name = "cloud_classes"
        return apply_class_mapping(classes, None, unmapped_to_missing=False)

    if mode == "external_file":
        if external_mask_path is None:
            raise ValueError("external_mask_path must be provided when mode='external_file'")
        return _external_classes(
            Path(external_mask_path),
            reference=ref,
            class_mapping=class_mapping,
            target_resolution_m=target_resolution_m,
            resolution_policy=resolution_policy,
            allow_upsample_to_target=allow_upsample_to_target,
            unmapped_to_missing=unmapped_to_missing,
        )

    if mode == "user_callable":
        if user_callable is None:
            raise ValueError("user_callable must be provided when mode='user_callable'")
        raw = _call_user_callable(
            user_callable,
            toa=toa,
            sensor_config=sensor_config,
            input_path=Path(external_mask_path) if external_mask_path else None,
        )
        if not isinstance(raw, xr.DataArray):
            raw = xr.DataArray(np.asarray(raw), dims=ref.dims, coords=ref.coords)
        raw = _extract_band(raw)
        if raw.shape != ref.shape:
            raw = reproject_match(raw, ref, resampling="nearest")
        raw = _resample_classes(
            raw,
            target_resolution_m=target_resolution_m,
            resolution_policy=resolution_policy,
            allow_upsample_to_target=allow_upsample_to_target,
        )
        return apply_class_mapping(raw, class_mapping, unmapped_to_missing=unmapped_to_missing)

    if mode != "auto":
        raise ValueError(f"Unsupported cloud mode: {mode!r}")

    red, green, nir = _prepare_rgbnir(
        toa,
        sensor_config,
        target_resolution_m=target_resolution_m,
        resolution_policy=resolution_policy,
        allow_upsample_to_target=allow_upsample_to_target,
    )

    if provider == "omnicloudmask":
        out = OmniCloudMaskProvider().predict(
            red,
            green,
            nir,
            class_mapping=class_mapping,
            unmapped_to_missing=unmapped_to_missing,
        )
    else:
        raise ValueError(f"Unsupported cloud provider: {provider!r}")

    out = apply_class_mapping(out, None, unmapped_to_missing=False)

    # Ensure valid standardized classes and expected naming.
    uniques = np.unique(out.values)
    if not {int(v) for v in uniques}.issubset(set(_EXPECTED_CLASS_VALUES)):
        raise ValueError(
            f"Cloud classes must be in {_EXPECTED_CLASS_VALUES}; got {uniques.tolist()}"
        )
    out.name = "cloud_classes"
    return out
