"""Resampling/alignment helpers extracted from :mod:`swir_refine`.

Pure DataArray grid-alignment utilities that don't depend on the Route-B
orchestration state.  Kept private (underscore) because they are implementation
details; swir_refine re-exports the names it uses so existing call sites keep
working.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from siac.runtime import BRDFKernelWeights
from siac.runtime.models import copy_spatial_metadata_like

if TYPE_CHECKING:
    import collections.abc

    from siac.algorithms.surface.brdf_monthly_composite import (
        MonthlyBestPixelComposite,
        MonthlyKernelWeightComposite,
    )
    from siac.domain import SensorBand


def _resample_da_indirect(
    data: xr.DataArray,
    target_shape: tuple[int, int],
    method: str,
    *,
    template: xr.DataArray | None = None,
) -> xr.DataArray:
    """Call into swir_refine._resample_da so monkey-patches applied to the
    swir_refine module surface remain effective after this refactor.

    The indirection lives in its own helper rather than at module import time
    to avoid a circular import (swir_refine -> this module -> swir_refine).
    """
    from siac.algorithms.surface import swir_refine
    return swir_refine._resample_da(data, target_shape, method, template=template)


def _shares_spatial_grid(data: xr.DataArray, template: xr.DataArray) -> bool:
    """Return True when *data* already lives on the same x/y grid as *template*."""
    if data.sizes.get("y") != template.sizes.get("y") or data.sizes.get("x") != template.sizes.get("x"):
        return False
    for axis in ("y", "x"):
        if axis in data.coords and axis in template.coords and not np.array_equal(
            np.asarray(data.coords[axis].values),
            np.asarray(template.coords[axis].values),
        ):
            return False
    return True


def _is_coarser_target_grid(
    source: xr.DataArray,
    template: xr.DataArray,
) -> bool:
    """True when *template* has fewer rows or cols than *source* (i.e. coarser)."""
    source_height = int(source.sizes.get("y", 0))
    source_width = int(source.sizes.get("x", 0))
    target_height = int(template.sizes.get("y", 0))
    target_width = int(template.sizes.get("x", 0))
    return source_height > target_height or source_width > target_width


def _monthly_composite_downsample_method(
    composite: MonthlyBestPixelComposite | MonthlyKernelWeightComposite,
    template: xr.DataArray,
) -> str:
    """Pick 'area' for downsampling and 'bilinear' for up/equal-sampling."""
    # Local import to avoid a circular dependency at module-import time.
    from siac.algorithms.surface.brdf_monthly_composite import MonthlyKernelWeightComposite
    source = composite.kernels.f0 if isinstance(composite, MonthlyKernelWeightComposite) else composite.reflectance
    if _is_coarser_target_grid(source, template):
        return "area"
    return "bilinear"


def _resample_band_cube_to_template(
    data: xr.DataArray,
    template: xr.DataArray,
    method: str,
) -> xr.DataArray:
    """Resample a (band, y, x) cube onto the *template* grid."""
    if _shares_spatial_grid(data, template):
        return copy_spatial_metadata_like(data.astype(np.float32), template)

    target_shape = (int(template.sizes["y"]), int(template.sizes["x"]))
    band_coords = data.coords["band"].values if "band" in data.coords else np.arange(data.sizes["band"])
    resampled = xr.concat(
        [
            _resample_da_indirect(data.sel(band=band, drop=True), target_shape, method, template=template)
            for band in band_coords
        ],
        dim=xr.IndexVariable("band", band_coords),
    )
    coords: dict[str, object] = {"band": band_coords}
    if "y" in template.coords:
        coords["y"] = template.coords["y"]
    if "x" in template.coords:
        coords["x"] = template.coords["x"]
    return copy_spatial_metadata_like(resampled.assign_coords(**coords).astype(np.float32), template)


def _resample_spatial_field_to_template(
    data: xr.DataArray,
    template: xr.DataArray,
    method: str,
) -> xr.DataArray:
    """Resample a 2-D (y, x) field onto the *template* grid."""
    if _shares_spatial_grid(data, template):
        return copy_spatial_metadata_like(data.astype(np.float32), template)
    target_shape = (int(template.sizes["y"]), int(template.sizes["x"]))
    resampled = _resample_da_indirect(data, target_shape, method, template=template)
    coords: dict[str, object] = {}
    if "y" in template.coords:
        coords["y"] = template.coords["y"]
    if "x" in template.coords:
        coords["x"] = template.coords["x"]
    return copy_spatial_metadata_like(resampled.assign_coords(**coords).astype(np.float32), template)


def _resample_brdf_weights_to_template(
    weights: BRDFKernelWeights,
    template: xr.DataArray,
    *,
    method: str,
) -> BRDFKernelWeights:
    """Align all kernel weight cubes + their uncertainties to *template*."""
    return BRDFKernelWeights(
        f0=_resample_band_cube_to_template(weights.f0, template, method),
        f1=_resample_band_cube_to_template(weights.f1, template, method),
        f2=_resample_band_cube_to_template(weights.f2, template, method),
        f0_unc=_resample_band_cube_to_template(weights.f0_unc, template, method),
        f1_unc=_resample_band_cube_to_template(weights.f1_unc, template, method),
        f2_unc=_resample_band_cube_to_template(weights.f2_unc, template, method),
        reflectance_unc=(
            _resample_band_cube_to_template(weights.reflectance_unc, template, method)
            if weights.reflectance_unc is not None
            else None
        ),
    )


def _deduplicate_bands(bands: collections.abc.Sequence[SensorBand]) -> list[SensorBand]:
    """Preserve order while removing bands that share a name."""
    seen: set[str] = set()
    ordered: list[SensorBand] = []
    for band in bands:
        if band.name in seen:
            continue
        seen.add(band.name)
        ordered.append(band)
    return ordered
