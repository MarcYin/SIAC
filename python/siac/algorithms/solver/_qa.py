"""QA-dataset assembly for the multigrid aerosol solver.

Builds the per-pixel QA layers attached to ``SolverResult``: exclusion
masks (invalid retrievals, missing observation support, water, sharp
transitions), parameter-boundary hits, the combined low-quality mask, and
the optional per-pixel fitting cost.

Extracted from ``multigrid.py`` to keep the solver module reviewable;
``MultiGridSolver._build_solver_qa_dataset`` delegates here.
"""

from __future__ import annotations

import numpy as np
import xarray as xr

from siac.runtime.models import copy_spatial_metadata_like


def _bound_tolerance(bounds: tuple[float, float]) -> float:
    span = max(float(bounds[1]) - float(bounds[0]), 0.0)
    return max(1.0e-6, span * 1.0e-5)


def _boundary_hit_masks(
    values: np.ndarray,
    bounds: tuple[float, float],
    valid_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    atol = _bound_tolerance(bounds)
    finite = np.isfinite(values)
    lower = valid_mask & finite & np.isclose(values, float(bounds[0]), rtol=0.0, atol=atol)
    upper = valid_mask & finite & np.isclose(values, float(bounds[1]), rtol=0.0, atol=atol)
    return lower, upper


def _mask_to_data_array(
    mask: np.ndarray,
    template: xr.DataArray,
) -> xr.DataArray:
    coords = {dim: template.coords[dim] for dim in template.dims if dim in template.coords}
    out = xr.DataArray(np.asarray(mask, dtype=bool), dims=template.dims, coords=coords)
    return copy_spatial_metadata_like(out, template)


def build_solver_qa_dataset(
    *,
    template: xr.DataArray,
    valid_mask: np.ndarray,
    aot: np.ndarray,
    tcwv: np.ndarray,
    invalid_mask: np.ndarray | None,
    zero_obs_mask: np.ndarray | None,
    insufficient_support_mask: np.ndarray | None,
    no_observation_mask: np.ndarray | None,
    sharp_transition_mask: np.ndarray | None,
    water_mask: np.ndarray | None,
    aot_bounds: tuple[float, float],
    tcwv_bounds: tuple[float, float],
    fitting_cost: np.ndarray | None = None,
) -> xr.Dataset:
    """Assemble the per-pixel solver QA dataset on the *template* grid."""
    valid = np.asarray(valid_mask, dtype=bool)
    invalid = (
        np.zeros_like(valid, dtype=bool)
        if invalid_mask is None
        else np.asarray(invalid_mask, dtype=bool) & valid
    )
    zero_obs = (
        np.zeros_like(valid, dtype=bool)
        if zero_obs_mask is None
        else np.asarray(zero_obs_mask, dtype=bool)
    )
    insufficient_support = (
        np.zeros_like(valid, dtype=bool)
        if insufficient_support_mask is None
        else np.asarray(insufficient_support_mask, dtype=bool)
    )
    no_observation = (
        np.zeros_like(valid, dtype=bool)
        if no_observation_mask is None
        else np.asarray(no_observation_mask, dtype=bool)
    )
    sharp_transition = (
        np.zeros_like(valid, dtype=bool)
        if sharp_transition_mask is None
        else np.asarray(sharp_transition_mask, dtype=bool)
    )
    water_excluded = (
        np.zeros_like(valid, dtype=bool)
        if water_mask is None
        else np.asarray(water_mask, dtype=bool)
    )
    aot_lower, aot_upper = _boundary_hit_masks(
        np.asarray(aot, dtype=np.float32),
        aot_bounds,
        valid,
    )
    tcwv_lower, tcwv_upper = _boundary_hit_masks(
        np.asarray(tcwv, dtype=np.float32),
        tcwv_bounds,
        valid,
    )
    parameter_boundary = aot_lower | aot_upper | tcwv_lower | tcwv_upper
    low_quality = (
        invalid
        | zero_obs
        | insufficient_support
        | no_observation
        | parameter_boundary
        | sharp_transition
        | water_excluded
    )

    qa_vars: dict[str, xr.DataArray] = {
        "invalid_retrieval": _mask_to_data_array(invalid, template),
        "zero_obs_support": _mask_to_data_array(zero_obs, template),
        "insufficient_observation_support": _mask_to_data_array(insufficient_support, template),
        "no_observation": _mask_to_data_array(no_observation, template),
        "sharp_transition_excluded": _mask_to_data_array(sharp_transition, template),
        "water_mask_excluded": _mask_to_data_array(water_excluded, template),
        "aot_lower_boundary": _mask_to_data_array(aot_lower, template),
        "aot_upper_boundary": _mask_to_data_array(aot_upper, template),
        "tcwv_lower_boundary": _mask_to_data_array(tcwv_lower, template),
        "tcwv_upper_boundary": _mask_to_data_array(tcwv_upper, template),
        "parameter_boundary": _mask_to_data_array(parameter_boundary, template),
        "low_quality": _mask_to_data_array(low_quality, template),
    }
    if fitting_cost is not None:
        cost_arr = np.asarray(fitting_cost, dtype=np.float32)
        if cost_arr.shape != template.shape:
            from scipy.ndimage import zoom

            cost_arr = zoom(
                cost_arr,
                (
                    template.shape[0] / cost_arr.shape[0],
                    template.shape[1] / cost_arr.shape[1],
                ),
                order=1,
            ).astype(np.float32)
            cost_arr = cost_arr[: template.shape[0], : template.shape[1]]
        solver_supported = valid & ~insufficient_support & ~no_observation
        cost_arr = np.where(
            solver_supported & np.isfinite(cost_arr),
            cost_arr,
            np.nan,
        ).astype(np.float32, copy=False)
        qa_vars["fitting_cost"] = xr.DataArray(
            cost_arr,
            dims=template.dims,
            coords=template.coords,
        )
    return xr.Dataset(qa_vars)
