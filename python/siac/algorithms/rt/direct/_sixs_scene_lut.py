"""Scene-LUT planning and interpolation helpers for the native 6S runner."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from siac.algorithms.rt.direct._sixs_session import _NativeBatchResult

_CASE_ARRAY_NAMES: tuple[str, ...] = (
    "sza_deg",
    "saa_deg",
    "vza_deg",
    "vaa_deg",
    "aot550",
    "tcwv_cm",
    "tco3_atmcm",
    "elevation_km",
)

#: Per-axis spread below which a scene-LUT axis collapses to a single (mean) node.
#: The 6S coefficients vary negligibly across these spans, but per-pixel *float*
#: geometry has thousands of distinct values, so without this the geometric axes
#: never collapse (``unique.size`` is large) and the LUT cross-products them — e.g.
#: a small AOI where the view/solar angles span <0.1 deg would still build
#: ``4**6`` geometric nodes (~500K cases) of full 6S evaluations. Collapsing each
#: near-constant axis to one node makes the per-build cost ~``N_aot * N_tcwv``.
#: Units: degrees for angles, atm-cm for ozone, km for elevation, unitless AOT, cm
#: TCWV. The grid-search aot550/tcwv_cm axes are supplied explicitly by the joint
#: LUT builder and are never routed through the tolerance collapse there.
_AXIS_COLLAPSE_TOLERANCE: dict[str, float] = {
    "sza_deg": 0.5,
    "saa_deg": 0.5,
    "vza_deg": 0.5,
    "vaa_deg": 0.5,
    "aot550": 0.005,
    "tcwv_cm": 0.05,
    "tco3_atmcm": 0.005,
    "elevation_km": 0.05,
}


@dataclass(frozen=True)
class _SceneLUTPlan:
    axes: dict[str, np.ndarray]
    grid_case_arrays: dict[str, np.ndarray]
    direct_case_count: int
    lut_case_count: int


def _build_scene_lut_axis(
    values: np.ndarray, target_size: int, *, tolerance: float = 0.0
) -> np.ndarray:
    finite = np.asarray(values[np.isfinite(values)], dtype=np.float64)
    if finite.size == 0:
        return np.zeros(1, dtype=np.float64)
    # Collapse a physically-near-constant axis to one representative node. Without
    # this, per-pixel float geometry (thousands of unique values, spread well below
    # the tolerance) keeps the full ``target_size`` nodes and the LUT cross-products
    # them needlessly. The single node uses the mean so per-pixel lookups land on it.
    if tolerance > 0.0 and float(finite.max() - finite.min()) <= tolerance:
        return np.array([float(finite.mean())], dtype=np.float64)
    unique = np.unique(finite)
    if unique.size <= max(1, target_size):
        return np.ascontiguousarray(unique, dtype=np.float64)
    if target_size <= 1:
        return np.ascontiguousarray(unique[:1], dtype=np.float64)
    quantiles = np.linspace(0.0, 1.0, target_size, dtype=np.float64)
    axis = np.quantile(finite, quantiles, method="linear")
    axis[0] = float(finite.min())
    axis[-1] = float(finite.max())
    axis = np.unique(np.asarray(axis, dtype=np.float64))
    if axis.size == 1 and unique.size > 1:
        axis = np.array([unique[0], unique[-1]], dtype=np.float64)
    return np.ascontiguousarray(axis, dtype=np.float64)


def _scene_lut_case_count(axes: dict[str, np.ndarray]) -> int:
    count = 1
    for axis in axes.values():
        count *= max(1, int(axis.size))
    return count


def _build_scene_lut_plan(
    case_arrays: dict[str, np.ndarray],
    *,
    max_nodes_per_axis: int,
    max_cases: int,
) -> _SceneLUTPlan:
    axes = {
        name: _build_scene_lut_axis(
            np.asarray(case_arrays[name], dtype=np.float64),
            max_nodes_per_axis,
            tolerance=_AXIS_COLLAPSE_TOLERANCE.get(name, 0.0),
        )
        for name in _CASE_ARRAY_NAMES
    }
    while _scene_lut_case_count(axes) > max_cases:
        reducible = [name for name, axis in axes.items() if axis.size > 1]
        if not reducible:
            break
        name = max(reducible, key=lambda item: axes[item].size)
        axes[name] = _build_scene_lut_axis(
            np.asarray(case_arrays[name], dtype=np.float64), axes[name].size - 1
        )

    mesh = np.meshgrid(*(axes[name] for name in _CASE_ARRAY_NAMES), indexing="ij")
    grid_case_arrays = {
        name: np.ascontiguousarray(mesh[idx].reshape(-1), dtype=np.float64)
        for idx, name in enumerate(_CASE_ARRAY_NAMES)
    }
    return _SceneLUTPlan(
        axes=axes,
        grid_case_arrays=grid_case_arrays,
        direct_case_count=int(np.asarray(case_arrays[_CASE_ARRAY_NAMES[0]]).size),
        lut_case_count=int(grid_case_arrays[_CASE_ARRAY_NAMES[0]].size),
    )


def _interpolate_scene_lut_outputs(
    plan: _SceneLUTPlan,
    native_outputs: _NativeBatchResult,
    case_arrays: dict[str, np.ndarray],
    selected_names: tuple[str, ...],
) -> dict[str, np.ndarray]:
    from scipy.interpolate import RegularGridInterpolator

    axes_order = _CASE_ARRAY_NAMES
    varying = [name for name in axes_order if plan.axes[name].size > 1]
    n_cases = int(np.asarray(case_arrays[axes_order[0]]).size)
    result: dict[str, np.ndarray] = {}
    full_shape = tuple(int(plan.axes[name].size) for name in axes_order)
    if not varying:
        for name in selected_names:
            value = float(np.asarray(native_outputs.outputs[name], dtype=np.float64).reshape(-1)[0])
            result[name] = np.full(n_cases, value, dtype=np.float64)
        return result

    sample_points = np.column_stack(
        [np.asarray(case_arrays[name], dtype=np.float64) for name in varying]
    )
    for name in selected_names:
        values = np.asarray(native_outputs.outputs[name], dtype=np.float64).reshape(full_shape)
        reduced = values
        for axis_index in reversed(range(len(axes_order))):
            if plan.axes[axes_order[axis_index]].size == 1:
                reduced = np.take(reduced, 0, axis=axis_index)
        interpolator = RegularGridInterpolator(
            tuple(plan.axes[axis_name] for axis_name in varying),
            reduced,
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        result[name] = np.ascontiguousarray(interpolator(sample_points), dtype=np.float64)
    return result


def _should_use_scene_lut(
    *,
    mode: str,
    direct_case_count: int,
    lut_case_count: int,
    min_pixels: int,
    required_speedup: float,
) -> bool:
    if mode == "direct":
        return False
    if mode == "scene_lut":
        return True
    if direct_case_count < min_pixels or lut_case_count <= 0 or lut_case_count >= direct_case_count:
        return False
    return (float(direct_case_count) / float(lut_case_count)) >= required_speedup


def _build_joint_grid_search_lut_plan(
    case_arrays: dict[str, np.ndarray],
    *,
    aot_axis: np.ndarray,
    tcwv_axis: np.ndarray,
    max_nodes_per_axis: int,
    max_cases: int,
) -> _SceneLUTPlan:
    """Build a scene-LUT plan with explicit aot/tcwv axes for joint grid-search reuse.

    Unlike :func:`_build_scene_lut_plan`, this builder takes the aot550 and
    tcwv_cm axes as inputs (rather than deriving them from per-pixel
    candidate values). The remaining six geometric/atmospheric axes are
    derived from the per-pixel ``case_arrays`` as usual. The trimming step
    that reduces total case count to fit ``max_cases`` only shrinks the
    geometric axes — the explicit aot/tcwv axes are preserved because their
    nodes must coincide with the grid-search candidate values for the
    block-grid-search reuse to be numerically exact at the grid points.
    """
    aot_axis_arr = np.ascontiguousarray(np.unique(np.asarray(aot_axis, dtype=np.float64)))
    tcwv_axis_arr = np.ascontiguousarray(np.unique(np.asarray(tcwv_axis, dtype=np.float64)))
    if aot_axis_arr.size == 0:
        aot_axis_arr = np.zeros(1, dtype=np.float64)
    if tcwv_axis_arr.size == 0:
        tcwv_axis_arr = np.zeros(1, dtype=np.float64)
    fixed_axes = {"aot550", "tcwv_cm"}
    axes: dict[str, np.ndarray] = {}
    for name in _CASE_ARRAY_NAMES:
        if name == "aot550":
            axes[name] = aot_axis_arr
        elif name == "tcwv_cm":
            axes[name] = tcwv_axis_arr
        else:
            axes[name] = _build_scene_lut_axis(
                np.asarray(case_arrays[name], dtype=np.float64),
                max_nodes_per_axis,
                tolerance=_AXIS_COLLAPSE_TOLERANCE.get(name, 0.0),
            )

    # Shrink only the geometric axes to fit the case budget.
    while _scene_lut_case_count(axes) > max_cases:
        reducible = [
            name for name in _CASE_ARRAY_NAMES if name not in fixed_axes and axes[name].size > 1
        ]
        if not reducible:
            break
        name = max(reducible, key=lambda item: axes[item].size)
        axes[name] = _build_scene_lut_axis(
            np.asarray(case_arrays[name], dtype=np.float64), axes[name].size - 1
        )

    mesh = np.meshgrid(*(axes[name] for name in _CASE_ARRAY_NAMES), indexing="ij")
    grid_case_arrays = {
        name: np.ascontiguousarray(mesh[idx].reshape(-1), dtype=np.float64)
        for idx, name in enumerate(_CASE_ARRAY_NAMES)
    }
    return _SceneLUTPlan(
        axes=axes,
        grid_case_arrays=grid_case_arrays,
        direct_case_count=int(np.asarray(case_arrays[_CASE_ARRAY_NAMES[0]]).size),
        lut_case_count=int(grid_case_arrays[_CASE_ARRAY_NAMES[0]].size),
    )
