"""Offline cost-field analysis for surface-driven AOD retrievals.

These helpers replay a solver cost cube around a point to explain WHY a scene
landed where it did — which band drove the minimum, how flat the curve was,
what a local window would have chosen. They are diagnostic tooling, not part of
the retrieval: nothing in :mod:`siac.algorithms.solver.surface_driven` calls
them, and they read dumped ``.npz`` cost cubes rather than live solver state.

They live here, beside their only callers, so the production solver module
carries only production code.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from siac._rust_compat import surface_driven_pool_argmin
from siac.algorithms.solver.surface_driven import (
    _cost_curve_diagnostics,
    _median_per_node,
)


def _prefix_surface_diagnostics(diagnostics: dict[str, Any], prefix: str) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in diagnostics.items():
        new_key = str(key)
        if new_key.startswith("surface_"):
            new_key = f"{prefix}_{new_key[len('surface_') :]}"
        else:
            new_key = f"{prefix}_{new_key}"
        out[new_key] = value
    return out


def _cost_field_window_indices(
    cube: np.ndarray,
    *,
    x: np.ndarray,
    y: np.ndarray,
    center_x: float,
    center_y: float,
    radius_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    cube_arr = np.asarray(cube, dtype=np.float64)
    if cube_arr.ndim != 3:
        raise ValueError(f"cost cube must have shape (aot, y, x); got {cube_arr.shape!r}")
    _n_aot, ny, nx = cube_arr.shape
    x_arr = np.asarray(x, dtype=np.float64).reshape(-1)
    y_arr = np.asarray(y, dtype=np.float64).reshape(-1)
    if x_arr.size != nx or y_arr.size != ny:
        raise ValueError(
            "cost-field x/y coordinates do not match cube shape: "
            f"x={x_arr.size}, y={y_arr.size}, cube_yx={(ny, nx)!r}"
        )
    if radius_m <= 0.0:
        raise ValueError("cost-field window radius must be positive.")

    x_idx = np.flatnonzero(np.abs(x_arr - float(center_x)) <= float(radius_m))
    y_idx = np.flatnonzero(np.abs(y_arr - float(center_y)) <= float(radius_m))
    if x_idx.size == 0 or y_idx.size == 0:
        raise ValueError(
            "cost-field window does not intersect the solver grid "
            f"(center=({center_x:.3f}, {center_y:.3f}), radius_m={radius_m:.1f})."
        )
    return cube_arr, x_arr, y_arr, x_idx, y_idx


def _local_cost_diagnostics(
    cube: np.ndarray,
    *,
    aot_axis: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    center_x: float,
    center_y: float,
    radius_m: float,
) -> dict[str, Any]:
    cube_arr, _x_arr, _y_arr, x_idx, y_idx = _cost_field_window_indices(
        cube,
        x=x,
        y=y,
        center_x=center_x,
        center_y=center_y,
        radius_m=radius_m,
    )
    cube_win = np.ascontiguousarray(cube_arr[:, y_idx][:, :, x_idx], dtype=np.float64)
    diagnostics = _prefix_surface_diagnostics(
        _cost_curve_diagnostics(cube_win, aot_axis),
        "local",
    )
    diagnostics.update(
        {
            "local_window_radius_m": float(radius_m),
            "local_window_shape_y": int(y_idx.size),
            "local_window_shape_x": int(x_idx.size),
        }
    )
    return diagnostics


def _local_band_diagnostics(
    *,
    band_names: list[str],
    aot_axis: np.ndarray,
    band_cost_cube: np.ndarray | None,
    band_residual_cube: np.ndarray | None,
    x: np.ndarray,
    y: np.ndarray,
    center_x: float,
    center_y: float,
    radius_m: float,
    final_aot: float,
) -> dict[str, Any]:
    if band_cost_cube is None:
        return {}
    band_cost = np.asarray(band_cost_cube, dtype=np.float64)
    if band_cost.ndim != 4 or band_cost.shape[0] != len(band_names):
        return {}
    first_band = band_cost[0]
    _cube_arr, _x_arr, _y_arr, x_idx, y_idx = _cost_field_window_indices(
        first_band,
        x=x,
        y=y,
        center_x=center_x,
        center_y=center_y,
        radius_m=radius_m,
    )
    band_residual = (
        np.asarray(band_residual_cube, dtype=np.float64)
        if band_residual_cube is not None
        and np.asarray(band_residual_cube).shape == band_cost.shape
        else None
    )
    final_node = int(np.nanargmin(np.abs(aot_axis - final_aot))) if np.isfinite(final_aot) else -1
    diagnostics: dict[str, Any] = {}
    if final_node >= 0:
        diagnostics["local_final_aot_node"] = float(aot_axis[final_node])

    argmins: list[float] = []
    for band_index, band_name in enumerate(band_names):
        cost_win = np.ascontiguousarray(
            band_cost[band_index][:, y_idx][:, :, x_idx],
            dtype=np.float64,
        )
        node_cost = _median_per_node(cost_win)
        finite = np.flatnonzero(np.isfinite(node_cost))
        if finite.size:
            best = int(finite[np.argmin(node_cost[finite])])
            argmin_aot = float(aot_axis[best])
            argmins.append(argmin_aot)
            diagnostics[f"local_band_{band_name}_argmin_aot"] = argmin_aot
            diagnostics[f"local_band_{band_name}_argmin_cost"] = float(node_cost[best])
        if final_node >= 0:
            diagnostics[f"local_band_{band_name}_cost_final_node"] = float(node_cost[final_node])
            if band_residual is not None:
                residual_win = np.ascontiguousarray(
                    band_residual[band_index][:, y_idx][:, :, x_idx],
                    dtype=np.float64,
                )
                node_residual = _median_per_node(residual_win)
                diagnostics[f"local_band_{band_name}_residual_final_node"] = float(
                    node_residual[final_node]
                )
    if len(argmins) >= 2:
        diagnostics["local_band_argmin_spread"] = float(max(argmins) - min(argmins))
    return diagnostics


def _scalar_int(value: np.ndarray, default: int) -> int:
    arr = np.asarray(value).reshape(-1)
    if arr.size == 0:
        return default
    try:
        return int(arr[0])
    except (TypeError, ValueError):
        return default


def _resolve_cost_field_window(
    cube: np.ndarray,
    *,
    aot_axis: np.ndarray,
    aot_prior: np.ndarray,
    aot_prior_unc: np.ndarray,
    solve_valid: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    center_x: float,
    center_y: float,
    radius_m: float,
    pool_window: int,
    min_count: int,
) -> dict[str, object]:
    cube_arr = np.asarray(cube, dtype=np.float64)
    if cube_arr.ndim != 3:
        raise ValueError(f"cost cube must have shape (aot, y, x); got {cube_arr.shape!r}")
    _n_aot, ny, nx = cube_arr.shape
    x_arr = np.asarray(x, dtype=np.float64).reshape(-1)
    y_arr = np.asarray(y, dtype=np.float64).reshape(-1)
    if x_arr.size != nx or y_arr.size != ny:
        raise ValueError(
            "cost-field x/y coordinates do not match cube shape: "
            f"x={x_arr.size}, y={y_arr.size}, cube_yx={(ny, nx)!r}"
        )
    if radius_m <= 0.0:
        raise ValueError("cost-field window radius must be positive.")

    x_idx = np.flatnonzero(np.abs(x_arr - float(center_x)) <= float(radius_m))
    y_idx = np.flatnonzero(np.abs(y_arr - float(center_y)) <= float(radius_m))
    if x_idx.size == 0 or y_idx.size == 0:
        raise ValueError(
            "cost-field window does not intersect the solver grid "
            f"(center=({center_x:.3f}, {center_y:.3f}), radius_m={radius_m:.1f})."
        )

    cube_win = np.ascontiguousarray(cube_arr[:, y_idx][:, :, x_idx], dtype=np.float64)
    prior_win = np.ascontiguousarray(np.asarray(aot_prior, dtype=np.float64)[y_idx][:, x_idx])
    unc_win = np.ascontiguousarray(np.asarray(aot_prior_unc, dtype=np.float64)[y_idx][:, x_idx])
    valid_win = np.ascontiguousarray(np.asarray(solve_valid, dtype=bool)[y_idx][:, x_idx])

    aod, unc, cost = surface_driven_pool_argmin(
        cube_win,
        np.asarray(aot_axis, dtype=np.float64),
        prior_win,
        unc_win,
        valid_win,
        int(pool_window),
        int(min_count),
    )
    aod_arr = np.asarray(aod, dtype=np.float64)
    unc_arr = np.asarray(unc, dtype=np.float64)
    cost_arr = np.asarray(cost, dtype=np.float64)
    finite = np.isfinite(aod_arr)
    if not np.any(finite):
        raise ValueError("integrated cost-field solve produced no finite pixels in the window.")

    cost_finite = np.isfinite(cost_arr) & finite
    unc_finite = np.isfinite(unc_arr) & finite
    return {
        "aod": float(np.nanmedian(aod_arr[finite])),
        "aod_unc": float(np.nanmedian(unc_arr[unc_finite])) if np.any(unc_finite) else None,
        "cost": float(np.nanmedian(cost_arr[cost_finite])) if np.any(cost_finite) else None,
        "n_finite": int(np.count_nonzero(finite)),
        "n_window": int(finite.size),
        "n_valid_input": int(np.count_nonzero(valid_win)),
        "window_shape": [int(y_idx.size), int(x_idx.size)],
        "x_range": [float(np.nanmin(x_arr[x_idx])), float(np.nanmax(x_arr[x_idx]))],
        "y_range": [float(np.nanmin(y_arr[y_idx])), float(np.nanmax(y_arr[y_idx]))],
    }


def integrated_cost_field_aod(
    *,
    cube: np.ndarray,
    aot_axis: np.ndarray,
    aot_prior: np.ndarray,
    aot_prior_unc: np.ndarray,
    solve_valid: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    center_x: float,
    center_y: float,
    radius_m: float = 1500.0,
    pool_window: int = 1,
    min_count: int = 1,
    cube_abs: np.ndarray | None = None,
    band_cost_cube: np.ndarray | None = None,
    band_residual_cube: np.ndarray | None = None,
    band_names: list[str] | None = None,
    clean_threshold: float = 0.15,
    high_threshold: float = 0.6,
) -> dict[str, object]:
    """Resolve a point retrieval by solving a cropped spatial cost field.

    The production raster solve returns a per-pixel AOT map and validation often
    samples the nearest pixel. The research harness instead crops the raw cost
    cube to the target window, runs the spatial median-pool/argmin there, and
    returns the median AOT over finite solved pixels. This helper implements
    that integrated site-window retrieval with the same Rust kernel used by the
    raster solver.
    """
    main = _resolve_cost_field_window(
        cube,
        aot_axis=aot_axis,
        aot_prior=aot_prior,
        aot_prior_unc=aot_prior_unc,
        solve_valid=solve_valid,
        x=x,
        y=y,
        center_x=center_x,
        center_y=center_y,
        radius_m=radius_m,
        pool_window=pool_window,
        min_count=min_count,
    )
    abs_result: dict[str, object] | None = None
    selected = main
    selected_pass = "main"
    mode = "single"

    if cube_abs is not None:
        cube_abs_arr = np.asarray(cube_abs)
        if cube_abs_arr.shape == np.asarray(cube).shape:
            abs_result = _resolve_cost_field_window(
                cube_abs_arr,
                aot_axis=aot_axis,
                aot_prior=aot_prior,
                aot_prior_unc=aot_prior_unc,
                solve_valid=solve_valid,
                x=x,
                y=y,
                center_x=center_x,
                center_y=center_y,
                radius_m=radius_m,
                pool_window=pool_window,
                min_count=min_count,
            )
            mode = "auto2"
            selected_pass = "shape"
            abs_aod = float(abs_result["aod"])
            if not (float(clean_threshold) <= abs_aod <= float(high_threshold)):
                selected = abs_result
                selected_pass = "abs"

    selected_aod = float(selected["aod"])
    local_diagnostics = _local_cost_diagnostics(
        cube,
        aot_axis=aot_axis,
        x=x,
        y=y,
        center_x=center_x,
        center_y=center_y,
        radius_m=radius_m,
    )
    if band_names is not None:
        local_diagnostics.update(
            _local_band_diagnostics(
                band_names=[str(name) for name in band_names],
                aot_axis=aot_axis,
                band_cost_cube=band_cost_cube,
                band_residual_cube=band_residual_cube,
                x=x,
                y=y,
                center_x=center_x,
                center_y=center_y,
                radius_m=radius_m,
                final_aot=selected_aod,
            )
        )

    return {
        "mode": mode,
        "selected_pass": selected_pass,
        "aod": selected_aod,
        "aod_unc": selected.get("aod_unc"),
        "cost": selected.get("cost"),
        "radius_m": float(radius_m),
        "pool_window": int(pool_window),
        "min_count": int(min_count),
        "clean_threshold": float(clean_threshold),
        "high_threshold": float(high_threshold),
        "main": main,
        "abs": abs_result,
        "diagnostics": local_diagnostics,
    }


def integrated_cost_field_aod_from_npz(
    path: str,
    *,
    center_x: float,
    center_y: float,
    radius_m: float = 1500.0,
    clean_threshold: float = 0.15,
    high_threshold: float = 0.6,
) -> dict[str, object]:
    """Load a ``SIAC_DUMP_COST_CUBE`` archive and resolve its target window."""
    with np.load(path, allow_pickle=False) as z:
        cube_abs = z["cube_abs"] if "cube_abs" in z.files and z["cube_abs"].size else None
        band_cost_cube = (
            z["band_cost_cube"]
            if "band_cost_cube" in z.files and z["band_cost_cube"].size
            else None
        )
        band_residual_cube = (
            z["band_residual_cube"]
            if "band_residual_cube" in z.files and z["band_residual_cube"].size
            else None
        )
        band_names = (
            [str(name) for name in z["band_names"].astype(str).tolist()]
            if "band_names" in z.files
            else None
        )
        return integrated_cost_field_aod(
            cube=z["cube"],
            cube_abs=cube_abs,
            band_cost_cube=band_cost_cube,
            band_residual_cube=band_residual_cube,
            band_names=band_names,
            aot_axis=z["aot_axis"],
            aot_prior=z["aot_prior"],
            aot_prior_unc=z["aot_prior_unc"],
            solve_valid=z["solve_valid"],
            x=z["x"],
            y=z["y"],
            center_x=center_x,
            center_y=center_y,
            radius_m=radius_m,
            pool_window=_scalar_int(
                z["pool_window"] if "pool_window" in z.files else np.array([]), 1
            ),
            min_count=_scalar_int(z["min_count"] if "min_count" in z.files else np.array([]), 1),
            clean_threshold=clean_threshold,
            high_threshold=high_threshold,
        )
