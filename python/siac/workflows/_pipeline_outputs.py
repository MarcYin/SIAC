"""Result-shaping helpers for SIAC pipeline outputs."""

from __future__ import annotations

from contextlib import suppress
from typing import Any

import numpy as np
import xarray as xr

from siac.runtime import MonthlyCompositeOutput
from siac.runtime.models import copy_spatial_metadata_like


def surface_template(data: xr.DataArray) -> xr.DataArray:
    if "band" in data.dims:
        band_coord = data.coords["band"].values[0] if "band" in data.coords else 0
        return data.sel(band=band_coord, drop=True)
    return data


def _band_name(value: object, index: int) -> str:
    if hasattr(value, "item"):
        with suppress(Exception):
            value = value.item()
    name = getattr(value, "name", None)
    if isinstance(name, str) and name:
        return name
    text = str(value)
    return text if text else f"band_{index + 1:02d}"


def banded_dataarray_to_dataset(
    data: xr.DataArray,
    *,
    default_name: str,
    template: xr.DataArray,
) -> xr.Dataset:
    if "band" not in data.dims:
        return xr.Dataset({default_name: copy_spatial_metadata_like(data, template)})

    band_values = (
        data.coords["band"].values if "band" in data.coords else np.arange(data.sizes["band"])
    )
    return xr.Dataset(
        {
            _band_name(band, index): copy_spatial_metadata_like(
                data.sel(band=band, drop=True),
                template,
            )
            for index, band in enumerate(band_values)
        }
    )


def monthly_composite_outputs(
    composites: tuple[Any, ...],
    *,
    template: xr.DataArray,
) -> dict[str, MonthlyCompositeOutput] | None:
    if not composites:
        return None

    outputs: dict[str, MonthlyCompositeOutput] = {}
    for composite in composites:
        label = f"{int(composite.year):04d}_{int(composite.month):02d}"
        outputs[label] = MonthlyCompositeOutput(
            reflectance=banded_dataarray_to_dataset(
                composite.reflectance,
                default_name="reflectance",
                template=template,
            ),
            quality=copy_spatial_metadata_like(composite.quality.astype(np.float32), template),
            sample_index=copy_spatial_metadata_like(
                composite.sample_index.astype(np.int16), template
            ),
        )
    return outputs
