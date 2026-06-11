"""Public helpers for driving pipeline stages standalone.

The comparison and diagnostic scripts under ``tools/`` assemble solver inputs
outside the full pipeline (e.g. to evaluate one RT backend or surface-prior
variant on a prepared scene). They previously reached into private
``siac.workflows.pipeline`` helpers, which broke whenever pipeline internals
moved. This module is the supported seam: the pipeline imports these same
functions, so tools and production cannot drift apart.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

import numpy as np

from siac.geo.resample import resample_field_to_template as resample_field_to_template

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    import xarray as xr

    from siac.runtime import (
        AtmosphericState,
        ObservationBundle,
        SolverInputBundle,
        SurfacePrior,
    )

    GridAssemblerFn = Callable[..., SolverInputBundle]

__all__ = [
    "aerosol_resolution",
    "call_grid_assembler",
    "resample_field_to_template",
    "select_band_slice",
]


def aerosol_resolution(config: Any) -> float:
    """Solver aerosol-grid resolution (m) from a resolved SIAC config."""
    solver_config = getattr(getattr(config, "algorithms", None), "solver", None)
    if solver_config is not None:
        return float(getattr(solver_config, "aerosol_resolution", 120.0))
    return 120.0


def select_band_slice(
    data: xr.DataArray,
    *,
    band_name: str,
    band_index: int,
) -> xr.DataArray | None:
    """Select one band from a banded array by name, falling back to index.

    Returns ``None`` when the array carries string band labels that do not
    include ``band_name`` (a positional fallback would silently pick the
    wrong band); band-less arrays pass through unchanged.
    """
    if "band" not in data.dims:
        return data
    band_coord = data.coords.get("band")
    if band_coord is not None:
        band_values = [str(value) for value in np.asarray(band_coord.values).tolist()]
        if band_name in band_values:
            return cast("xr.DataArray", data.sel(band=band_name, drop=True))
        if np.asarray(band_coord.values).dtype.kind in {"U", "S", "O"}:
            return None
    return cast("xr.DataArray", data.isel(band=band_index, drop=True))


def call_grid_assembler(
    grid_assembler: GridAssemblerFn,
    obs: ObservationBundle,
    atmo: AtmosphericState,
    surface: SurfacePrior,
    rt_model: Any,
    *,
    aerosol_resolution_m: float,
    sharp_transition_filter: Any | None = None,
    water_mask_path: str | Path | None = None,
    water_mask_cache_dir: str | Path | None = None,
    water_mask_buffer_pixels: int = 0,
    solver_band_names: tuple[str, ...] | None = None,
    reproject_cache_dir: str | Path | None = None,
) -> SolverInputBundle:
    """Call the grid assembler with the current standardized interface."""
    return grid_assembler(
        obs,
        atmo,
        surface,
        rt_model,
        aerosol_resolution_m=aerosol_resolution_m,
        sharp_transition_filter=sharp_transition_filter,
        water_mask_path=water_mask_path,
        water_mask_cache_dir=water_mask_cache_dir,
        water_mask_buffer_pixels=water_mask_buffer_pixels,
        solver_band_names=solver_band_names,
        reproject_cache_dir=reproject_cache_dir,
    )
