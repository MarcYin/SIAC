#!/usr/bin/env python3
"""Build a self-contained diagnostic report for the six high-AOD failures.

All analysis is performed from saved campaign artifacts. The script does not
run radiative transfer, query Earth Engine, or modify retrieval outputs.
"""

from __future__ import annotations

import argparse
import base64
import csv
import html
import io
import json
import math
import re
import subprocess
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from pyproj import Transformer

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence


REPO = Path(__file__).resolve().parents[2]
ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_OUT = ROOT / "highaod6_surface_diagnostics_20260710"
SHELL = Path(
    "/home/users/marcyin/.codex/plugins/cache/openai-curated-remote/"
    "data-analytics/0.2.6-d37358633e00/assets/html-report-shell.html"
)
EMBED = Path(
    "/home/users/marcyin/.codex/plugins/cache/openai-curated-remote/"
    "data-analytics/0.2.6-d37358633e00/skills/build-report/scripts/"
    "embed_html_report_runtime.py"
)

MIDS = (
    "Alta_Floresta_IF__T21LXK_20240912T140049",
    "London-CDN__T17TMH_20240724T160901",
    "Metsi__T35JNL_20240801T075609",
    "Mont_Joli__T19UEP_20240816T153559",
    "NEON_DEJU__T06VWR_20240724T211511",
    "Waskesiu__T13UDV_20240720T180921",
)

DISPLAY_SITE = {
    "Alta_Floresta_IF": "Alta Floresta",
    "London-CDN": "London",
    "Metsi": "Metsi",
    "Mont_Joli": "Mont Joli",
    "NEON_DEJU": "NEON DEJU",
    "Waskesiu": "Waskesiu",
}

CURRENT_RESULTS = ROOT / "phaseD_results_campaign250_surface_bad89_highaod6_backstop_escape_l2awvp_20260710"
CURRENT_CUBES = ROOT / "phaseD_cost_cubes_campaign250_surface_bad89_highaod6_backstop_escape_20260710_l2awvp"
PRIOR_WVP_RESULTS = ROOT / "phaseD_results_campaign250_surface_bad89_highaod6_backstop_escape_20260710"
R2_RESULTS = ROOT / "phaseD_results_campaign250_surface_bad89_R2_localdiag_20260705"
R2_CAMPAIGN_RESULTS = ROOT / "phaseD_results_campaign250_R2_full"
L1C_RESULTS = ROOT / "phaseD_results_campaign250_surface_bad89_l1clut_floor015_20260705"
L2A_PC_RESULTS = ROOT / "phaseD_results_campaign250_l2a_pc_direct_acix_20260705"
BIOMASS_RESULTS = ROOT / "phaseD_results_smoke_biomass_ssa"
SDSPEC_RESULTS = ROOT / "seasonal_val"
L1C_PRIORS = ROOT / "l1c_seasonal_species_lut_clean015_highaod6"
CALIB_DUMPS = ROOT / "calib_dumps_c250"
MATCHUPS = ROOT / "matchups" / "matchups.csv"
AEROSOL_CONTEXT = ROOT / "gee_aerosol_context_campaign250"
S2_L2A_AOT = ROOT / "s2_l2a_aot_campaign250"
S2_L2A_WVP = ROOT / "s2_l2a_wvp_campaign250"
MAIAC = ROOT / "maiac_qa"

BAND_INDEX = {"B02": 0, "B04": 1}
PRIOR_COMP_INDEX = {"B02": 1, "B03": 2, "B04": 3, "B8A": 4, "B11": 5, "B12": 6}
WAVELENGTH_NM = {"B02": 490, "B04": 665, "B8A": 865, "B11": 1610, "B12": 2190}
COLORS = {
    "blue": "#0169cc",
    "blue_dark": "#003f7a",
    "gold": "#b47a00",
    "orange": "#e25507",
    "olive": "#637300",
    "pink": "#c23b72",
    "neutral": "#6f6f6f",
    "light": "#d8d8d8",
}


def _load_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return {}


def _matchup_rows() -> dict[str, dict[str, str]]:
    with MATCHUPS.open(encoding="utf-8") as handle:
        return {row["matchup_id"]: row for row in csv.DictReader(handle)}


def _finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _first_number(record: dict[str, Any], *keys: str) -> float | None:
    for key in keys:
        value: Any = record
        for part in key.split("."):
            if not isinstance(value, dict) or part not in value:
                value = None
                break
            value = value[part]
        output = _finite(value)
        if output is not None:
            return output
    return None


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def _median(values: Iterable[float]) -> float | None:
    array = np.asarray(list(values), dtype=np.float64)
    array = array[np.isfinite(array)]
    return float(np.median(array)) if array.size else None


def _corr(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    if x.size < 3 or np.std(x) == 0.0 or np.std(y) == 0.0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def _rank(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=np.float64)
    start = 0
    while start < values.size:
        stop = start + 1
        while stop < values.size and values[order[stop]] == values[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + stop - 1) + 1.0
        start = stop
    return ranks


def _spearman(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    valid = np.isfinite(x) & np.isfinite(y)
    if int(valid.sum()) < 3:
        return None
    return _corr(_rank(x[valid]), _rank(y[valid]))


def _within_ee(aod: float | None, truth: float) -> bool | None:
    if aod is None or not math.isfinite(aod):
        return None
    return abs(aod - truth) <= 0.05 + 0.15 * truth


def _source_label(path: Path) -> str:
    if path.parent == ROOT:
        return path.name
    return f"{path.parent.name}/{path.name}"


def _window_mask_xy(
    x: np.ndarray,
    y: np.ndarray,
    *,
    lon: float,
    lat: float,
    crs: str,
    radius_m: float = 1500.0,
) -> tuple[np.ndarray, int, int]:
    center_x, center_y = Transformer.from_crs(
        "EPSG:4326", crs, always_xy=True
    ).transform(lon, lat)
    xx, yy = np.meshgrid(x, y)
    mask = np.square(xx - center_x) + np.square(yy - center_y) <= radius_m**2
    ix = int(np.argmin(np.abs(x - center_x)))
    iy = int(np.argmin(np.abs(y - center_y)))
    return mask, iy, ix


def _window_mask_transform(
    shape: tuple[int, int],
    transform: np.ndarray,
    *,
    lon: float,
    lat: float,
    crs: str,
    radius_m: float = 1500.0,
) -> tuple[np.ndarray, int, int]:
    a, _, xoff, _, e, yoff = [float(value) for value in transform]
    x = xoff + (np.arange(shape[1]) + 0.5) * a
    y = yoff + (np.arange(shape[0]) + 0.5) * e
    return _window_mask_xy(x, y, lon=lon, lat=lat, crs=crs, radius_m=radius_m)


def _masked_stats(values: np.ndarray, mask: np.ndarray) -> dict[str, float | None]:
    selected = np.asarray(values, dtype=np.float64)[mask]
    selected = selected[np.isfinite(selected)]
    if selected.size == 0:
        return {"count": 0, "median": None, "p05": None, "p95": None}
    return {
        "count": int(selected.size),
        "median": float(np.median(selected)),
        "p05": float(np.percentile(selected, 5)),
        "p95": float(np.percentile(selected, 95)),
    }


def _median_curve(cube: np.ndarray, mask: np.ndarray) -> np.ndarray:
    values = np.asarray(cube, dtype=np.float64)
    selected = values[:, mask]
    with np.errstate(invalid="ignore"):
        return np.nanmedian(selected, axis=1)


def _nanmedian(values: np.ndarray, axis: int | tuple[int, ...] | None = None) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmedian(values, axis=axis)


def _nanmean(values: np.ndarray, axis: int | tuple[int, ...] | None = None) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(values, axis=axis)


def _argmin_map(cube: np.ndarray, axis: np.ndarray) -> np.ndarray:
    values = np.asarray(cube, dtype=np.float64)
    finite = np.any(np.isfinite(values), axis=0)
    work = np.where(np.isfinite(values), values, np.inf)
    indices = np.argmin(work, axis=0)
    output = np.asarray(axis, dtype=np.float64)[indices]
    output[~finite] = np.nan
    return output


def _boa_at(calib: np.lib.npyio.NpzFile, band: str, tau: float) -> np.ndarray:
    grid = np.asarray(calib["aot_grid"], dtype=np.float64)
    tau_clipped = float(np.clip(tau, grid[0], grid[-1]))
    log_grid = np.log(grid)
    coeffs = []
    for name in ("xap", "xbp", "xcp"):
        curve = np.asarray(calib[f"{name}_{band}"], dtype=np.float64)
        coeffs.append(float(np.interp(np.log(tau_clipped), log_grid, curve)))
    xap, xbp, xcp = coeffs
    toa = np.asarray(calib[f"toa_{band}"], dtype=np.float64)
    y = xap * toa - xbp
    with np.errstate(invalid="ignore", divide="ignore"):
        return y / (1.0 + xcp * y)


def _normalise_rgb(channels: Sequence[np.ndarray], gamma: float = 1.35) -> np.ndarray:
    arrays = [np.asarray(channel, dtype=np.float64) for channel in channels]
    shape = min((array.shape for array in arrays), key=lambda item: item[0] * item[1])
    arrays = [array[: shape[0], : shape[1]] for array in arrays]
    rgb = np.stack(arrays, axis=-1)
    valid = rgb[np.isfinite(rgb)]
    high = float(np.percentile(valid, 99)) if valid.size else 0.3
    high = max(high, 0.05)
    return np.nan_to_num(np.clip(rgb / high, 0.0, 1.0) ** (1.0 / gamma), nan=0.0)


def _mosaic_png(
    *,
    cost: np.lib.npyio.NpzFile,
    prior: np.lib.npyio.NpzFile | None,
    calib: np.lib.npyio.NpzFile | None,
    truth_index: int,
    main_argmin: np.ndarray,
    abs_argmin: np.ndarray,
    station_index: tuple[int, int],
    prior_station_index: tuple[int, int] | None,
    calib_station_index: tuple[int, int] | None,
    truth: float,
) -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(14.5, 7.2), facecolor="white")
    axes = axes.ravel()

    def marker(ax: Any, index: tuple[int, int] | None) -> None:
        if index is None:
            return
        iy, ix = index
        ax.plot(ix, iy, "+", color="white", markersize=11, markeredgewidth=2.0)
        ax.plot(ix, iy, "+", color="black", markersize=8, markeredgewidth=1.0)

    if calib is not None and all(f"toa_{band}" in calib.files for band in ("B8A", "B04", "B02")):
        image = _normalise_rgb(
            [calib["toa_B8A"], calib["toa_B04"], calib["toa_B02"]]
        )
        axes[0].imshow(image)
        marker(axes[0], calib_station_index)
        axes[0].set_title("Target TOA false colour\nNIR / red / blue")
    else:
        axes[0].text(0.5, 0.5, "Target TOA image unavailable", ha="center", va="center")

    if prior is not None and "comp" in prior.files:
        comp = np.asarray(prior["comp"], dtype=np.float64)
        median_comp = _nanmedian(comp, axis=0)
        image = _normalise_rgb(
            [
                median_comp[PRIOR_COMP_INDEX["B04"]],
                median_comp[PRIOR_COMP_INDEX["B03"]],
                median_comp[PRIOR_COMP_INDEX["B02"]],
            ]
        )
        axes[1].imshow(image)
        marker(axes[1], prior_station_index)
        axes[1].set_title("L1C clean-day prior\nmedian true colour")
    else:
        axes[1].text(0.5, 0.5, "Clean-day prior unavailable", ha="center", va="center")

    scalar_panels = [
        (np.asarray(cost["toa"])[0], "Target TOA B02", "Blues_r", None, None),
        (np.asarray(cost["boa_prior"])[0], "Current BOA prior B02", "Blues_r", None, None),
        (main_argmin, "Raw shape-cost argmin AOD", "viridis", 0.0, max(4.0, truth)),
        (abs_argmin, "Raw absolute-cost argmin AOD", "viridis", 0.0, max(4.0, truth)),
        (
            np.asarray(cost["band_residual_cube"])[0, truth_index],
            "|BOA - prior| at truth, B02",
            "magma",
            0.0,
            None,
        ),
        (
            np.asarray(cost["band_residual_cube"])[1, truth_index],
            "|BOA - prior| at truth, B04",
            "magma",
            0.0,
            None,
        ),
    ]
    for ax, (values, title, cmap, vmin, vmax) in zip(axes[2:], scalar_panels):
        array = np.asarray(values, dtype=np.float64)
        if vmax is None:
            finite = array[np.isfinite(array)]
            vmax = float(np.percentile(finite, 98)) if finite.size else 1.0
        image = ax.imshow(array, cmap=cmap, vmin=vmin, vmax=vmax)
        marker(ax, station_index)
        ax.set_title(title)
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle("Spatial evidence: target, prior, raw AOD preference, and truth-node closure", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=125, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return base64.b64encode(buffer.getvalue()).decode("ascii")


def _analyse_arrays(
    mid: str,
    current: dict[str, Any],
) -> tuple[dict[str, Any], str]:
    cube_path = CURRENT_CUBES / f"{mid}.npz"
    prior_path = L1C_PRIORS / f"{mid}.npz"
    calib_path = CALIB_DUMPS / f"{mid}.npz"
    truth = float(current["truth"])
    retrieved = float(current["retrieved"])
    lon = float(current["lon"])
    lat = float(current["lat"])
    crs = str(current["scene_crs"])

    prior_handle = np.load(prior_path, allow_pickle=False) if prior_path.exists() else None
    calib_handle = np.load(calib_path, allow_pickle=False) if calib_path.exists() else None
    try:
        with np.load(cube_path, allow_pickle=False) as cost:
            axis = np.asarray(cost["aot_axis"], dtype=np.float64)
            x = np.asarray(cost["x"], dtype=np.float64)
            y = np.asarray(cost["y"], dtype=np.float64)
            solve_valid = np.asarray(cost["solve_valid"], dtype=bool)
            local, station_y, station_x = _window_mask_xy(
                x, y, lon=lon, lat=lat, crs=crs, radius_m=1500.0
            )
            local_valid = local & solve_valid
            if int(local_valid.sum()) < 20:
                local_valid = solve_valid
            truth_index = int(np.argmin(np.abs(axis - truth)))
            retrieved_index = int(np.argmin(np.abs(axis - retrieved)))
            main_cube = np.asarray(cost["cube"], dtype=np.float64)
            abs_cube = np.asarray(cost["cube_abs"], dtype=np.float64)
            band_cost = np.asarray(cost["band_cost_cube"], dtype=np.float64)
            band_residual = np.asarray(cost["band_residual_cube"], dtype=np.float64)
            boa_unc = np.asarray(cost["boa_unc"], dtype=np.float64)
            main_curve = _median_curve(main_cube, local_valid)
            abs_curve = _median_curve(abs_cube, local_valid)
            band_cost_curves = [_median_curve(band_cost[index], local_valid) for index in range(2)]
            band_residual_curves = [
                _median_curve(band_residual[index], local_valid) for index in range(2)
            ]
            main_argmin = _argmin_map(main_cube, axis)
            abs_argmin = _argmin_map(abs_cube, axis)
            main_best = int(np.nanargmin(main_curve))
            abs_best = int(np.nanargmin(abs_curve))

            def log_delta(curve: np.ndarray) -> list[float | None]:
                finite = curve[np.isfinite(curve)]
                floor = float(np.min(finite)) if finite.size else 0.0
                values = np.log10(1.0 + np.maximum(curve - floor, 0.0))
                return [float(value) if math.isfinite(float(value)) else None for value in values]

            curve_rows = []
            main_log = log_delta(main_curve)
            abs_log = log_delta(abs_curve)
            for index, aod in enumerate(axis):
                if main_log[index] is not None:
                    curve_rows.append({"aod": float(aod), "series": "Shape cost", "value": main_log[index]})
                if abs_log[index] is not None:
                    curve_rows.append({"aod": float(aod), "series": "Absolute cost", "value": abs_log[index]})

            residual_rows = []
            for band_index, band in enumerate(("B02", "B04")):
                for index, aod in enumerate(axis):
                    value = float(band_residual_curves[band_index][index])
                    if math.isfinite(value):
                        residual_rows.append({"aod": float(aod), "series": band, "value": value})

            band_rows: list[dict[str, Any]] = []
            toa = np.asarray(cost["toa"], dtype=np.float64)
            boa_prior = np.asarray(cost["boa_prior"], dtype=np.float64)
            for band_index, band in enumerate(("B02", "B04")):
                sigma = np.maximum(boa_unc[band_index], 1.0e-6)
                z_truth = band_residual[band_index, truth_index] / sigma
                z_retrieved = band_residual[band_index, retrieved_index] / sigma
                cost_curve = band_cost_curves[band_index]
                band_best = int(np.nanargmin(cost_curve))
                band_rows.append(
                    {
                        "band": band,
                        "toa": _masked_stats(toa[band_index], local_valid),
                        "prior": _masked_stats(boa_prior[band_index], local_valid),
                        "uncertainty": _masked_stats(boa_unc[band_index], local_valid),
                        "argmin_aod": float(axis[band_best]),
                        "residual_at_retrieved": float(band_residual_curves[band_index][retrieved_index]),
                        "residual_at_truth": float(band_residual_curves[band_index][truth_index]),
                        "z_at_retrieved": float(np.nanmedian(z_retrieved[local_valid])),
                        "z_at_truth": float(np.nanmedian(z_truth[local_valid])),
                    }
                )

            ee = 0.05 + 0.15 * truth
            raw_local = main_argmin[local_valid]
            raw_abs_local = abs_argmin[local_valid]
            exact = {
                "aot_axis": [float(value) for value in axis],
                "truth_node": float(axis[truth_index]),
                "retrieved_node": float(axis[retrieved_index]),
                "shape_argmin_aod": float(axis[main_best]),
                "absolute_argmin_aod": float(axis[abs_best]),
                "shape_cost_at_min": float(main_curve[main_best]),
                "shape_cost_at_truth": float(main_curve[truth_index]),
                "absolute_cost_at_min": float(abs_curve[abs_best]),
                "absolute_cost_at_truth": float(abs_curve[truth_index]),
                "curve_rows": curve_rows,
                "residual_rows": residual_rows,
                "bands": band_rows,
                "solve_valid_count": int(solve_valid.sum()),
                "solve_total_count": int(solve_valid.size),
                "solve_valid_fraction": float(solve_valid.mean()),
                "local_valid_count": int(local_valid.sum()),
                "raw_shape_aod_median": float(np.nanmedian(raw_local)),
                "raw_abs_aod_median": float(np.nanmedian(raw_abs_local)),
                "raw_shape_lower_rail_fraction": float(np.mean(raw_local <= axis[0] + 1.0e-6)),
                "raw_abs_lower_rail_fraction": float(np.mean(raw_abs_local <= axis[0] + 1.0e-6)),
                "raw_shape_within_ee_fraction": float(np.mean(np.abs(raw_local - truth) <= ee)),
                "raw_abs_within_ee_fraction": float(np.mean(np.abs(raw_abs_local - truth) <= ee)),
            }

            history_rows: list[dict[str, Any]] = []
            prior_summary: dict[str, Any] = {"available": False}
            prior_station_index: tuple[int, int] | None = None
            if prior_handle is not None and "comp" in prior_handle.files:
                comp = np.asarray(prior_handle["comp"], dtype=np.float64)
                transform = np.asarray(prior_handle["transform"], dtype=np.float64)
                epsg = int(np.asarray(prior_handle["epsg"]).item())
                prior_mask, prior_y, prior_x = _window_mask_transform(
                    comp.shape[2:],
                    transform,
                    lon=lon,
                    lat=lat,
                    crs=f"EPSG:{epsg}",
                    radius_m=1500.0,
                )
                prior_station_index = (prior_y, prior_x)
                realizations = [str(value) for value in prior_handle["realizations"].astype(str).tolist()]
                quality = (
                    json.loads(str(prior_handle["clean_quality_json"].item()))
                    if "clean_quality_json" in prior_handle.files
                    else []
                )
                for index, window in enumerate(realizations):
                    row: dict[str, Any] = {"window": window}
                    for band, band_index in (("B02", 1), ("B04", 3)):
                        selected = comp[index, band_index][prior_mask]
                        selected = selected[np.isfinite(selected)]
                        row[band] = float(np.median(selected)) if selected.size else None
                    if index < len(quality):
                        row["selected_aod_max"] = _finite(quality[index].get("selected_aod_max"))
                        row["selected_aod_median"] = _finite(quality[index].get("selected_aod_median"))
                        row["n_selected"] = _finite(quality[index].get("n_selected"))
                        row["n_fallback"] = _finite(quality[index].get("n_fallback"))
                    history_rows.append(row)
                clean_excess = [
                    max(0.0, float(row.get("selected_aod_max") or 0.15) - 0.15)
                    for row in history_rows
                ]
                prior_summary = {
                    "available": True,
                    "realizations": len(realizations),
                    "history_rows": history_rows,
                    "fallback_scenes": int(sum(float(row.get("n_fallback") or 0.0) for row in history_rows)),
                    "median_selected_aod_excess": float(np.median(clean_excess)) if clean_excess else 0.0,
                    "clean_aod_max": _finite(prior_handle["clean_aod_max"].item()) if "clean_aod_max" in prior_handle.files else None,
                    "clean_min_scenes": _finite(prior_handle["clean_min_scenes"].item()) if "clean_min_scenes" in prior_handle.files else None,
                }
                for band in ("B02", "B04", "B8A", "B11", "B12"):
                    band_index = PRIOR_COMP_INDEX[band]
                    values = comp[:, band_index]
                    median_surface = _nanmedian(values, axis=0)
                    temporal_rmse = np.sqrt(
                        _nanmean(
                            np.square(values - median_surface[np.newaxis, ...]),
                            axis=0,
                        )
                    )
                    prior_summary[band] = {
                        "prior_median": float(np.nanmedian(median_surface[prior_mask])),
                        "temporal_rmse": float(np.nanmedian(temporal_rmse[prior_mask])),
                    }

            signed_closure: dict[str, Any] = {"available": False}
            calib_station_index: tuple[int, int] | None = None
            if calib_handle is not None:
                shape = tuple(int(value) for value in np.asarray(calib_handle["template_shape"]).tolist())
                transform = np.asarray(calib_handle["template_transform"], dtype=np.float64)
                calib_crs = str(np.asarray(calib_handle["template_crs"]).item())
                calib_mask, calib_y, calib_x = _window_mask_transform(
                    shape,
                    transform,
                    lon=lon,
                    lat=lat,
                    crs=calib_crs,
                    radius_m=1500.0,
                )
                calib_station_index = (calib_y, calib_x)
                signed_closure = {"available": True, "bands": {}, "spectral": []}
                for band in ("B02", "B04", "B8A", "B11", "B12"):
                    required = {
                        f"toa_{band}",
                        f"xap_{band}",
                        f"xbp_{band}",
                        f"xcp_{band}",
                    }
                    if not required.issubset(calib_handle.files):
                        continue
                    toa = np.asarray(calib_handle[f"toa_{band}"], dtype=np.float64)
                    boa_truth = _boa_at(calib_handle, band, truth)
                    boa_retrieved = _boa_at(calib_handle, band, retrieved)
                    valid = (
                        calib_mask
                        & np.isfinite(toa)
                        & np.isfinite(boa_truth)
                        & np.isfinite(boa_retrieved)
                    )
                    spectral_row = {
                        "band": band,
                        "wavelength_nm": WAVELENGTH_NM[band],
                        "toa": float(np.nanmedian(toa[valid])),
                        "boa_at_truth": float(np.nanmedian(boa_truth[valid])),
                        "boa_at_retrieved": float(np.nanmedian(boa_retrieved[valid])),
                        "truth_minus_retrieved_boa": float(
                            np.nanmedian((boa_truth - boa_retrieved)[valid])
                        ),
                    }
                    signed_closure["spectral"].append(spectral_row)
                    pred_key = f"pred_{band}"
                    if pred_key not in calib_handle.files:
                        continue
                    pred = np.asarray(calib_handle[pred_key], dtype=np.float64)
                    valid = calib_mask & np.isfinite(pred) & np.isfinite(boa_truth)
                    signed_closure["bands"][band] = {
                        "boa_at_truth": float(np.nanmedian(boa_truth[valid])),
                        "boa_at_retrieved": float(np.nanmedian(boa_retrieved[valid])),
                        "predicted_prior": float(np.nanmedian(pred[valid])),
                        "truth_minus_prior": float(np.nanmedian((boa_truth - pred)[valid])),
                        "retrieved_minus_prior": float(np.nanmedian((boa_retrieved - pred)[valid])),
                    }

            mosaic = _mosaic_png(
                cost=cost,
                prior=prior_handle,
                calib=calib_handle,
                truth_index=truth_index,
                main_argmin=main_argmin,
                abs_argmin=abs_argmin,
                station_index=(station_y, station_x),
                prior_station_index=prior_station_index,
                calib_station_index=calib_station_index,
                truth=truth,
            )
    finally:
        if prior_handle is not None:
            prior_handle.close()
        if calib_handle is not None:
            calib_handle.close()

    return {
        "exact": exact,
        "clean_prior": prior_summary,
        "signed_closure_proxy": signed_closure,
        "sources": {
            "cost_cube": _source_label(cube_path),
            "clean_prior": _source_label(prior_path) if prior_path.exists() else None,
            "calibration_dump": _source_label(calib_path) if calib_path.exists() else None,
        },
    }, mosaic


def _load_sdspec(mid: str, suffix: str) -> dict[str, Any]:
    return _load_json(SDSPEC_RESULTS / f"{mid}__seas_maiac_sdspec_highaod6_clean015{suffix}.json")


def _method_row(
    *,
    name: str,
    kind: str,
    value: float | None,
    truth: float,
    source: str,
) -> dict[str, Any]:
    return {
        "method": name,
        "kind": kind,
        "aod": value,
        "error": value - truth if value is not None else None,
        "within_ee": _within_ee(value, truth),
        "source": source,
    }


def _aerosol_context(mid: str) -> dict[str, Any]:
    record = _load_json(AEROSOL_CONTEXT / f"{mid}.json")
    values = record.get("values") if isinstance(record.get("values"), dict) else {}
    species = (
        ("Black carbon", "black_carbon", "BCEXTTAU"),
        ("Organic matter", "organic_matter", "OCEXTTAU"),
        ("Dust", "dust", "DUEXTTAU"),
        ("Sea salt", "sea_salt", "SSEXTTAU"),
        ("Sulphate", "sulphate", "SUEXTTAU"),
    )
    rows = []
    for label, cams_key, merra_key in species:
        rows.append(
            {
                "species": label,
                "CAMS": _finite(values.get(f"cams_{cams_key}_aerosol_optical_depth_at_550nm_surface_frac")),
                "MERRA-2": _finite(values.get(f"merra_{merra_key}_frac")),
            }
        )
    return {
        "status": record.get("status"),
        "rows": rows,
        "cams_total": _finite(values.get("cams_total_aerosol_optical_depth_at_550nm_surface")),
        "merra_total": _finite(values.get("merra_TOTEXTTAU")),
        "merra_scatter_fraction": _finite(values.get("merra_scatter_fraction")),
        "merra_angstrom": _finite(values.get("merra_TOTANGSTR")),
    }


def _case_record(mid: str, matchup: dict[str, str]) -> tuple[dict[str, Any], str]:
    current_path = CURRENT_RESULTS / f"{mid}.json"
    current = _load_json(current_path)
    truth = float(current["truth"])
    arrays, mosaic = _analyse_arrays(mid, current)
    prior_wvp = _load_json(PRIOR_WVP_RESULTS / f"{mid}.json")
    r2 = _load_json(R2_RESULTS / f"{mid}.json")
    l1c = _load_json(L1C_RESULTS / f"{mid}.json")
    l2a_pc = _load_json(L2A_PC_RESULTS / f"{mid}.json")
    biomass = _load_json(BIOMASS_RESULTS / f"{mid}.json")
    sdspec = _load_sdspec(mid, "")
    sdspec_abs = _load_sdspec(mid, "_abs")
    l2a_aot = _load_json(S2_L2A_AOT / f"{mid}.json")
    l2a_wvp = _load_json(S2_L2A_WVP / f"{mid}.json")
    maiac = _load_json(MAIAC / f"{mid}.json")
    aerosol = _aerosol_context(mid)

    methods = [
        _method_row(
            name="Current surface solve (L2A WVP)",
            kind="surface solve",
            value=_finite(current.get("retrieved")),
            truth=truth,
            source=_source_label(current_path),
        ),
        _method_row(
            name="Same solve (prior WVP)",
            kind="surface solve ablation",
            value=_finite(prior_wvp.get("retrieved")),
            truth=truth,
            source=_source_label(PRIOR_WVP_RESULTS / f"{mid}.json"),
        ),
        _method_row(
            name="R2 local diagnostic",
            kind="surface solve baseline",
            value=_finite(r2.get("retrieved")),
            truth=truth,
            source=_source_label(R2_RESULTS / f"{mid}.json"),
        ),
        _method_row(
            name="L1C clean-day prior",
            kind="surface-prior experiment",
            value=_finite(l1c.get("retrieved")),
            truth=truth,
            source=_source_label(L1C_RESULTS / f"{mid}.json"),
        ),
        _method_row(
            name="L2A PC direct ACIX",
            kind="surface-prior experiment",
            value=_finite(l2a_pc.get("retrieved")),
            truth=truth,
            source=_source_label(L2A_PC_RESULTS / f"{mid}.json"),
        ),
        _method_row(
            name="acixThree sdspec shape",
            kind="acixThree-style solve",
            value=_finite(sdspec.get("aod")),
            truth=truth,
            source=_source_label(SDSPEC_RESULTS / f"{mid}__seas_maiac_sdspec_highaod6_clean015.json"),
        ),
        _method_row(
            name="acixThree sdspec absolute",
            kind="acixThree-style solve",
            value=_finite(sdspec_abs.get("aod")),
            truth=truth,
            source=_source_label(SDSPEC_RESULTS / f"{mid}__seas_maiac_sdspec_highaod6_clean015_abs.json"),
        ),
        _method_row(
            name="Prescribed biomass 6S",
            kind="aerosol-optics experiment",
            value=_finite(biomass.get("retrieved")),
            truth=truth,
            source=_source_label(BIOMASS_RESULTS / f"{mid}.json"),
        ),
        _method_row(
            name="Sentinel-2 L2A AOT",
            kind="external same-scene context",
            value=_finite(l2a_aot.get("retrieved")),
            truth=truth,
            source=_source_label(S2_L2A_AOT / f"{mid}.json"),
        ),
        _method_row(
            name="MAIAC",
            kind="external aerosol context",
            value=_first_number(maiac, "aot", "aod", "aot_mean"),
            truth=truth,
            source=_source_label(MAIAC / f"{mid}.json"),
        ),
        _method_row(
            name="CAMS",
            kind="external aerosol context",
            value=aerosol["cams_total"],
            truth=truth,
            source=_source_label(AEROSOL_CONTEXT / f"{mid}.json"),
        ),
        _method_row(
            name="MERRA-2",
            kind="external aerosol context",
            value=aerosol["merra_total"],
            truth=truth,
            source=_source_label(AEROSOL_CONTEXT / f"{mid}.json"),
        ),
    ]

    solver = current.get("solver") if isinstance(current.get("solver"), dict) else {}
    band_spread = _finite(solver.get("surface_band_argmin_spread")) or 0.0
    curve_min = _finite(solver.get("surface_cost_curve_min_aot"))
    band_minima = [
        _finite(solver.get("surface_band_B02_argmin_aot")),
        _finite(solver.get("surface_band_B04_argmin_aot")),
    ]
    rail_collapse = bool(
        curve_min is not None
        and curve_min <= 0.02
        and all(value is not None and value <= 0.05 for value in band_minima)
    )
    band_conflict = band_spread >= 0.5
    current_aod = float(current["retrieved"])
    biomass_aod = _finite(biomass.get("retrieved"))
    biomass_gain = (
        abs(current_aod - truth) - abs(biomass_aod - truth)
        if biomass_aod is not None
        else None
    )
    prior_wvp_aod = _finite(prior_wvp.get("retrieved"))
    wvp_delta = abs(current_aod - prior_wvp_aod) if prior_wvp_aod is not None else None
    if rail_collapse:
        primary = "Two-band lower-rail collapse"
    elif band_conflict:
        primary = "B02/B04 conflict plus spatial aggregation"
    else:
        primary = "Broad surface/RT mismatch"
    if biomass_aod is not None and _within_ee(biomass_aod, truth):
        species_result = "Biomass optics recover this case"
    elif biomass_gain is not None and biomass_gain >= 0.25:
        species_result = "Biomass optics help but do not close the error"
    else:
        species_result = "Biomass optics do not resolve the case"

    return {
        "matchup_id": mid,
        "site": str(current["site"]),
        "display_site": DISPLAY_SITE.get(str(current["site"]), str(current["site"])),
        "truth": truth,
        "ee_threshold": 0.05 + 0.15 * truth,
        "current": current,
        "matchup": matchup,
        "l2a_wvp": l2a_wvp,
        "sdspec": sdspec,
        "sdspec_abs": sdspec_abs,
        "biomass": biomass,
        "aerosol": aerosol,
        "methods": methods,
        **arrays,
        "diagnosis": {
            "primary": primary,
            "rail_collapse": rail_collapse,
            "band_conflict": band_conflict,
            "species_result": species_result,
            "biomass_error_gain": biomass_gain,
            "wvp_retrieval_delta": wvp_delta,
            "cloud_mask_saturated": float(current.get("cloud_frac", 0.0)) >= 0.95,
        },
        "sources": {
            **arrays["sources"],
            "current_result": _source_label(current_path),
            "matchup": _source_label(MATCHUPS),
            "aerosol_context": _source_label(AEROSOL_CONTEXT / f"{mid}.json"),
        },
    }, mosaic


def _campaign_context() -> dict[str, Any]:
    rows = []
    for path in sorted(R2_CAMPAIGN_RESULTS.glob("*.json")):
        row = _load_json(path)
        truth = _finite(row.get("truth"))
        retrieved = _finite(row.get("retrieved"))
        cloud = _finite(row.get("cloud_frac"))
        solver = row.get("solver") if isinstance(row.get("solver"), dict) else {}
        cost = _finite(solver.get("cost_final_per_band"))
        if truth is None or retrieved is None:
            continue
        rows.append(
            {
                "matchup_id": str(row.get("matchup_id")),
                "site": str(row.get("site")),
                "truth": truth,
                "abs_error": abs(retrieved - truth),
                "cloud_frac": cloud,
                "log_cost": math.log10(1.0 + cost) if cost is not None and cost >= 0.0 else None,
                "group": "High-AOD six" if row.get("matchup_id") in MIDS else "Campaign context",
            }
        )
    cloud_rows = [row for row in rows if row["cloud_frac"] is not None]
    cost_rows = [row for row in rows if row["log_cost"] is not None]
    return {
        "rows": rows,
        "count": len(rows),
        "cloud_count": len(cloud_rows),
        "cost_count": len(cost_rows),
        "cloud_pearson": _corr(
            [float(row["cloud_frac"]) for row in cloud_rows],
            [float(row["abs_error"]) for row in cloud_rows],
        ),
        "cloud_spearman": _spearman(
            [float(row["cloud_frac"]) for row in cloud_rows],
            [float(row["abs_error"]) for row in cloud_rows],
        ),
        "cost_pearson": _corr(
            [float(row["log_cost"]) for row in cost_rows],
            [float(row["abs_error"]) for row in cost_rows],
        ),
        "cost_spearman": _spearman(
            [float(row["log_cost"]) for row in cost_rows],
            [float(row["abs_error"]) for row in cost_rows],
        ),
    }


def _esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


class Tooltips:
    def __init__(self) -> None:
        self.index = 0

    def wrap(self, value: str, dataset: str, source: str = "saved SIAC artifacts") -> str:
        self.index += 1
        tip_id = f"source-value-{self.index}"
        return (
            f'<span class="source-tooltip" tabindex="0" aria-describedby="{tip_id}">'
            f"{_esc(value)}"
            f'<span class="source-tooltip-content" id="{tip_id}" role="tooltip">'
            f"Source: {_esc(source)}<br>Dataset: {_esc(dataset)}</span></span>"
        )


def _fmt(value: float | None, digits: int = 3) -> str:
    return "n/a" if value is None else f"{value:.{digits}f}"


def _fmt_pct(value: float | None, digits: int = 1) -> str:
    return "n/a" if value is None else f"{100.0 * value:.{digits}f}%"


def _line_svg(
    rows: list[dict[str, Any]],
    *,
    x_field: str,
    y_field: str,
    series_field: str,
    aria_label: str,
    markers: Sequence[tuple[float, str, str]] = (),
) -> str:
    width, height = 960, 400
    left, right, top, bottom = 68, 24, 28, 58
    plot_w, plot_h = width - left - right, height - top - bottom
    values = [
        (float(row[x_field]), float(row[y_field]))
        for row in rows
        if _finite(row.get(x_field)) is not None and _finite(row.get(y_field)) is not None
    ]
    if not values:
        return f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{_esc(aria_label)}"></svg>'
    xs = [item[0] for item in values]
    ys = [item[1] for item in values]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    if y_max <= y_min:
        y_max = y_min + 1.0
    y_pad = 0.08 * (y_max - y_min)
    y_min = max(0.0, y_min - y_pad)
    y_max += y_pad

    def px(value: float) -> float:
        return left + (value - x_min) / max(x_max - x_min, 1.0e-12) * plot_w

    def py(value: float) -> float:
        return top + (y_max - value) / max(y_max - y_min, 1.0e-12) * plot_h

    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{_esc(aria_label)}">',
        "<style>.axis{fill:var(--secondary);font:12px var(--ds-font)}.grid{stroke:var(--grid);stroke-width:1}.legend{fill:var(--text);font:12px var(--ds-font)}</style>",
    ]
    for index in range(6):
        value = y_min + (y_max - y_min) * index / 5.0
        y = py(value)
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}"/>')
        parts.append(f'<text class="axis" x="{left-9}" y="{y+4:.1f}" text-anchor="end">{value:.2f}</text>')
    for index in range(6):
        value = x_min + (x_max - x_min) * index / 5.0
        x = px(value)
        parts.append(f'<text class="axis" x="{x:.1f}" y="{height-bottom+25}" text-anchor="middle">{value:.2f}</text>')
    series_names = list(dict.fromkeys(str(row[series_field]) for row in rows))
    palette = [COLORS["blue"], COLORS["orange"], COLORS["olive"], COLORS["pink"], COLORS["gold"]]
    for index, name in enumerate(series_names):
        selected = sorted(
            [row for row in rows if str(row[series_field]) == name and _finite(row.get(y_field)) is not None],
            key=lambda row: float(row[x_field]),
        )
        points = " ".join(f"{px(float(row[x_field])):.1f},{py(float(row[y_field])):.1f}" for row in selected)
        color = palette[index % len(palette)]
        parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="2" points="{points}"/>')
        legend_x = left + index * 155
        parts.append(f'<line x1="{legend_x}" y1="{height-15}" x2="{legend_x+18}" y2="{height-15}" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text class="legend" x="{legend_x+24}" y="{height-11}">{_esc(name)}</text>')
    for value, label, color in markers:
        if x_min <= value <= x_max:
            x = px(value)
            parts.append(f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top+plot_h}" stroke="{color}" stroke-dasharray="5 4"/>')
            parts.append(f'<text class="axis" x="{x+4:.1f}" y="{top+12}">{_esc(label)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def _grouped_bar_svg(
    rows: list[dict[str, Any]],
    *,
    category_field: str,
    series_field: str,
    value_field: str,
    aria_label: str,
    percent: bool = False,
) -> str:
    width, height = 960, 420
    left, right, top, bottom = 68, 24, 28, 88
    categories = list(dict.fromkeys(str(row[category_field]) for row in rows))
    series = list(dict.fromkeys(str(row[series_field]) for row in rows))
    values = [float(row[value_field]) for row in rows if _finite(row.get(value_field)) is not None]
    y_max = max(values, default=1.0)
    y_max = max(1.0, y_max) if percent else max(y_max * 1.12, 0.1)
    plot_w, plot_h = width - left - right, height - top - bottom
    group_w = plot_w / max(len(categories), 1)
    gap = 3.0
    bar_w = max(3.0, min(28.0, (group_w - 14.0) / max(len(series), 1) - gap))
    palette = [COLORS["blue"], COLORS["orange"], COLORS["olive"], COLORS["pink"], COLORS["gold"]]
    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{_esc(aria_label)}">',
        "<style>.axis{fill:var(--secondary);font:11px var(--ds-font)}.value{fill:var(--text);font:10px var(--ds-font)}.grid{stroke:var(--grid);stroke-width:1}</style>",
    ]
    for index in range(6):
        value = y_max * index / 5.0
        y = top + plot_h * (1.0 - value / y_max)
        label = f"{100*value:.0f}%" if percent else f"{value:.2f}"
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}"/>')
        parts.append(f'<text class="axis" x="{left-8}" y="{y+4:.1f}" text-anchor="end">{label}</text>')
    lookup = {(str(row[category_field]), str(row[series_field])): row for row in rows}
    for category_index, category in enumerate(categories):
        center = left + group_w * (category_index + 0.5)
        total = len(series) * bar_w + max(0, len(series) - 1) * gap
        start = center - total / 2.0
        for series_index, name in enumerate(series):
            row = lookup.get((category, name))
            value = _finite(row.get(value_field)) if row else None
            if value is None:
                continue
            x = start + series_index * (bar_w + gap)
            y = top + plot_h * (1.0 - value / y_max)
            bar_height = top + plot_h - y
            color = palette[series_index % len(palette)]
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_height:.1f}" fill="{color}" rx="2"/>')
        parts.append(f'<text class="axis" x="{center:.1f}" y="{height-bottom+22}" text-anchor="middle">{_esc(category)}</text>')
    legend_x = left
    for index, name in enumerate(series):
        color = palette[index % len(palette)]
        parts.append(f'<rect x="{legend_x}" y="{height-24}" width="11" height="11" fill="{color}" rx="2"/>')
        parts.append(f'<text class="axis" x="{legend_x+16}" y="{height-14}">{_esc(name)}</text>')
        legend_x += max(130, 8 * len(name) + 35)
    parts.append("</svg>")
    return "".join(parts)


def _scatter_svg(
    rows: list[dict[str, Any]],
    *,
    x_field: str,
    y_field: str,
    group_field: str,
    aria_label: str,
) -> str:
    width, height = 960, 420
    left, right, top, bottom = 68, 24, 28, 58
    points = [row for row in rows if _finite(row.get(x_field)) is not None and _finite(row.get(y_field)) is not None]
    xs = [float(row[x_field]) for row in points]
    ys = [float(row[y_field]) for row in points]
    x_min, x_max = min(xs, default=0.0), max(xs, default=1.0)
    y_min, y_max = 0.0, max(ys, default=1.0) * 1.08
    plot_w, plot_h = width - left - right, height - top - bottom

    def px(value: float) -> float:
        return left + (value - x_min) / max(x_max - x_min, 1.0e-12) * plot_w

    def py(value: float) -> float:
        return top + (y_max - value) / max(y_max - y_min, 1.0e-12) * plot_h

    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{_esc(aria_label)}">',
        "<style>.axis{fill:var(--secondary);font:11px var(--ds-font)}.grid{stroke:var(--grid);stroke-width:1}</style>",
    ]
    for index in range(6):
        y_value = y_max * index / 5.0
        y = py(y_value)
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}"/>')
        parts.append(f'<text class="axis" x="{left-8}" y="{y+4:.1f}" text-anchor="end">{y_value:.2f}</text>')
        x_value = x_min + (x_max - x_min) * index / 5.0
        x = px(x_value)
        parts.append(f'<text class="axis" x="{x:.1f}" y="{height-bottom+23}" text-anchor="middle">{x_value:.2f}</text>')
    for row in points:
        focal = str(row[group_field]) == "High-AOD six"
        color = COLORS["orange"] if focal else COLORS["blue"]
        radius = 5.2 if focal else 2.6
        opacity = 0.95 if focal else 0.42
        parts.append(f'<circle cx="{px(float(row[x_field])):.1f}" cy="{py(float(row[y_field])):.1f}" r="{radius}" fill="{color}" opacity="{opacity}"/>')
    parts.append(f'<circle cx="{left}" cy="{height-15}" r="4" fill="{COLORS["blue"]}" opacity="0.55"/><text class="axis" x="{left+10}" y="{height-11}">Campaign context</text>')
    parts.append(f'<circle cx="{left+150}" cy="{height-15}" r="5" fill="{COLORS["orange"]}"/><text class="axis" x="{left+161}" y="{height-11}">High-AOD six</text>')
    parts.append("</svg>")
    return "".join(parts)


def _payload_bar(
    chart_id: str,
    title: str,
    rows: list[dict[str, Any]],
    *,
    category: str,
    value: str,
    series: str,
    percent: bool = False,
) -> dict[str, Any]:
    return {
        "id": chart_id,
        "height": 360,
        "type": "bar",
        "dataset": {
            "id": chart_id,
            "title": title,
            "data": rows,
            "chart_spec": {
                "id": chart_id,
                "dataset": chart_id,
                "title": title,
                "type": "bar",
                "encodings": {
                    "x": {"field": category, "type": "nominal"},
                    "y": {"field": value, "type": "quantitative"},
                    "color": {"field": series, "type": "nominal"},
                },
                "xAxisTitle": "",
                "yAxisTitle": "",
                "valueFormat": "percent" if percent else "number",
                "settings": {"groupMode": "grouped"},
            },
        },
    }


def _payload_line(
    chart_id: str,
    title: str,
    rows: list[dict[str, Any]],
    *,
    x: str,
    y: str,
    series: str,
) -> dict[str, Any]:
    return {
        "id": chart_id,
        "height": 340,
        "type": "line",
        "dataset": {
            "id": chart_id,
            "title": title,
            "data": rows,
            "chart_spec": {
                "id": chart_id,
                "dataset": chart_id,
                "title": title,
                "type": "line",
                "encodings": {
                    "x": {"field": x, "type": "quantitative"},
                    "y": {"field": y, "type": "quantitative"},
                    "color": {"field": series, "type": "nominal"},
                },
                "xAxisTitle": "AOD550",
                "yAxisTitle": "",
                "valueFormat": "number",
            },
        },
    }


def _payload_scatter(
    chart_id: str,
    title: str,
    rows: list[dict[str, Any]],
    *,
    x: str,
    y: str,
    group: str,
) -> dict[str, Any]:
    return {
        "id": chart_id,
        "height": 360,
        "type": "scatter",
        "dataset": {
            "id": chart_id,
            "title": title,
            "data": rows,
            "chart_spec": {
                "id": chart_id,
                "dataset": chart_id,
                "title": title,
                "type": "scatter",
                "encodings": {
                    "x": {"field": x, "type": "quantitative"},
                    "y": {"field": y, "type": "quantitative"},
                    "color": {"field": group, "type": "nominal"},
                    "tooltip": [
                        {"field": "site", "type": "nominal"},
                        {"field": "truth", "type": "quantitative"},
                    ],
                },
                "xAxisTitle": "",
                "yAxisTitle": "Absolute AOD error",
                "valueFormat": "number",
            },
        },
    }


def _chart_figure(
    *,
    chart_id: str,
    title: str,
    subtitle: str,
    fallback: str,
    note: str,
    datasets: str,
) -> str:
    source_id = f"chart-source-{chart_id}"
    return f"""
      <div class="wide">
        <figure class="card source-figure">
          <div class="card-head"><h3>{_esc(title)}</h3><p>{_esc(subtitle)}</p></div>
          <div class="chart-wrap">
            <div data-recharts-chart="{_esc(chart_id)}">
              <div class="chart-fallback" data-recharts-fallback>{fallback}</div>
              <div data-recharts-live aria-hidden="true"></div>
            </div>
          </div>
          <figcaption class="chart-note">{_esc(note)}</figcaption>
          <button type="button" class="source-tooltip" aria-describedby="{source_id}">Source
            <span class="source-tooltip-content" id="{source_id}" role="tooltip">Source: saved SIAC experiment artifacts<br>Dataset: {_esc(datasets)}</span>
          </button>
        </figure>
      </div>
    """


def _image_figure(*, image_base64: str, site: str, datasets: str) -> str:
    source_id = f"image-source-{re.sub(r'[^a-z0-9]+', '-', site.lower()).strip('-')}"
    return f"""
      <figure class="card source-figure map-figure">
        <div class="card-head"><h3>Spatial evidence mosaic</h3><p>Station marker is shown as a black-and-white cross; AOD maps are raw per-pixel argmins before production pooling.</p></div>
        <div class="map-image"><img src="data:image/png;base64,{image_base64}" alt="Spatial target, prior, AOD argmin, and residual maps for {_esc(site)}"></div>
        <figcaption class="chart-note">The false-colour target and signed closure proxy come from the calibration capture; cost, residual, TOA B02, and BOA prior maps come from the exact L2A-WVP run.</figcaption>
        <button type="button" class="source-tooltip" aria-describedby="{source_id}">Source
          <span class="source-tooltip-content" id="{source_id}" role="tooltip">Source: saved SIAC experiment artifacts<br>Dataset: {_esc(datasets)}</span>
        </button>
      </figure>
    """


def _method_table(case: dict[str, Any], tips: Tooltips) -> str:
    rows = []
    for method in case["methods"]:
        value = method["aod"]
        error = method["error"]
        within = method["within_ee"]
        if value is None:
            value_html = "n/a"
            error_html = "n/a"
            status = "not available"
        else:
            value_html = tips.wrap(_fmt(value), method["source"])
            error_html = tips.wrap(f"{error:+.3f}", method["source"])
            status = "within EE" if within else "outside EE"
        rows.append(
            f"<tr><td>{_esc(method['method'])}</td><td>{_esc(method['kind'])}</td>"
            f"<td>{value_html}</td><td>{error_html}</td><td>{_esc(status)}</td></tr>"
        )
    return "".join(rows)


def _band_table(case: dict[str, Any], tips: Tooltips) -> str:
    rows = []
    dataset = case["sources"]["cost_cube"]
    proxy_dataset = case["sources"].get("calibration_dump") or "calibration dump unavailable"
    proxy_bands = case["signed_closure_proxy"].get("bands", {})
    for band in case["exact"]["bands"]:
        name = band["band"]
        signed = _finite((proxy_bands.get(name) or {}).get("truth_minus_prior"))
        cells = [
            tips.wrap(_fmt(band["toa"]["median"]), dataset),
            tips.wrap(_fmt(band["prior"]["median"]), dataset),
            tips.wrap(_fmt(band["uncertainty"]["median"]), dataset),
            tips.wrap(_fmt(band["argmin_aod"]), dataset),
            tips.wrap(_fmt(band["residual_at_retrieved"]), dataset),
            tips.wrap(_fmt(band["residual_at_truth"]), dataset),
            tips.wrap(_fmt(band["z_at_truth"], 1), dataset),
            tips.wrap(_fmt(signed), proxy_dataset),
        ]
        rows.append(f"<tr><td>{name}</td>" + "".join(f"<td>{cell}</td>" for cell in cells) + "</tr>")
    return "".join(rows)


def _history_table(case: dict[str, Any], tips: Tooltips) -> str:
    rows = []
    dataset = case["sources"].get("clean_prior") or "clean prior unavailable"
    for row in case["clean_prior"].get("history_rows", []):
        rows.append(
            "<tr>"
            f"<td>{_esc(row['window'])}</td>"
            f"<td>{tips.wrap(_fmt(_finite(row.get('B02'))), dataset)}</td>"
            f"<td>{tips.wrap(_fmt(_finite(row.get('B04'))), dataset)}</td>"
            f"<td>{tips.wrap(_fmt(_finite(row.get('selected_aod_median'))), dataset)}</td>"
            f"<td>{tips.wrap(_fmt(_finite(row.get('selected_aod_max'))), dataset)}</td>"
            f"<td>{tips.wrap(_fmt(_finite(row.get('n_fallback')), 0), dataset)}</td>"
            "</tr>"
        )
    return "".join(rows)


def _spectral_table(case: dict[str, Any], tips: Tooltips) -> str:
    rows = []
    calib_dataset = case["sources"].get("calibration_dump") or "calibration dump unavailable"
    prior_dataset = case["sources"].get("clean_prior") or "clean prior unavailable"
    for row in case["signed_closure_proxy"].get("spectral", []):
        band = str(row["band"])
        clean_prior = _finite((case["clean_prior"].get(band) or {}).get("prior_median"))
        cells = (
            tips.wrap(_fmt(_finite(row.get("wavelength_nm")), 0), calib_dataset),
            tips.wrap(_fmt(_finite(row.get("toa"))), calib_dataset),
            tips.wrap(_fmt(_finite(row.get("boa_at_retrieved"))), calib_dataset),
            tips.wrap(_fmt(_finite(row.get("boa_at_truth"))), calib_dataset),
            tips.wrap(_fmt(_finite(row.get("truth_minus_retrieved_boa"))), calib_dataset),
            tips.wrap(_fmt(clean_prior), prior_dataset),
        )
        rows.append(
            f"<tr><td>{_esc(band)}</td>"
            + "".join(f"<td>{cell}</td>" for cell in cells)
            + "</tr>"
        )
    return "".join(rows)


def _case_section(
    case: dict[str, Any],
    mosaic: str,
    tips: Tooltips,
    charts: list[dict[str, Any]],
) -> str:
    site = case["display_site"]
    slug = re.sub(r"[^a-z0-9]+", "-", site.lower()).strip("-")
    current = case["current"]
    solver = current.get("solver", {})
    exact = case["exact"]
    diagnosis = case["diagnosis"]
    result_dataset = case["sources"]["current_result"]
    cube_dataset = case["sources"]["cost_cube"]

    def truth_tip() -> str:
        return tips.wrap(_fmt(case["truth"]), result_dataset)

    def retrieved_tip() -> str:
        return tips.wrap(_fmt(_finite(current.get("retrieved"))), result_dataset)

    def ee_tip() -> str:
        return tips.wrap(_fmt(case["ee_threshold"]), result_dataset)

    def cost_tip() -> str:
        return tips.wrap(
            _fmt(_finite(solver.get("cost_final_per_band")), 1), result_dataset
        )

    def cloud_tip() -> str:
        return tips.wrap(_fmt_pct(_finite(current.get("cloud_frac"))), result_dataset)

    wvp_html = tips.wrap(
        _fmt(_finite(current.get("target_tcwv_cm"))),
        _source_label(S2_L2A_WVP / f"{case['matchup_id']}.json"),
    )

    if diagnosis["rail_collapse"]:
        opening = (
            f"Both B02 and B04 prefer the lower AOD rail while AERONET is {truth_tip()}. "
            "That is a cost-model failure, not a local optimizer miss: moving toward the truth makes the surface/RT mismatch worse in both bands."
        )
    else:
        opening = (
            f"The aggregate shape-cost minimum is near AERONET, but B02 and B04 disagree by "
            f"{tips.wrap(_fmt(_finite(solver.get('surface_band_argmin_spread'))), result_dataset)} AOD. "
            "The production spatial aggregation then returns a lower scene mean."
        )
    species_gain = _finite(diagnosis.get("biomass_error_gain"))
    biomass_dataset = _source_label(
        BIOMASS_RESULTS / f"{case['matchup_id']}.json"
    )
    prior_wvp_dataset = _source_label(
        PRIOR_WVP_RESULTS / f"{case['matchup_id']}.json"
    )
    if species_gain is not None:
        species_sentence = (
            f" Prescribed biomass optics reduce absolute error by "
            f"{tips.wrap(_fmt(species_gain), biomass_dataset)}; "
            f"{_esc(diagnosis['species_result'].lower())}."
        )
    else:
        species_sentence = " No comparable biomass result is available."
    wvp_delta = _finite(diagnosis.get("wvp_retrieval_delta"))
    wvp_sentence = (
        f" Replacing the inherited target water vapour with same-scene L2A WVP changes AOD by only "
        f"{tips.wrap(_fmt(wvp_delta), result_dataset + '; ' + prior_wvp_dataset)}, "
        "so TCWV consistency is not the controlling failure here."
    )

    cost_chart_id = f"cost-{slug}"
    cost_rows = [row for row in exact["curve_rows"] if row["value"] is not None]
    charts.append(
        _payload_line(cost_chart_id, f"Cost shape: {site}", cost_rows, x="aod", y="value", series="series")
    )
    cost_figure = _chart_figure(
        chart_id=cost_chart_id,
        title="Local aggregate cost curves",
        subtitle="log10(1 + cost - each curve minimum), 1.5 km station window",
        fallback=_line_svg(
            cost_rows,
            x_field="aod",
            y_field="value",
            series_field="series",
            aria_label=f"Local aggregate cost curves for {site}",
            markers=(
                (case["truth"], "AERONET", COLORS["orange"]),
                (float(current["retrieved"]), "retrieved", COLORS["blue_dark"]),
            ),
        ),
        note="Each curve is shifted by its own minimum so the basin shape remains visible despite very different absolute costs.",
        datasets=cube_dataset,
    )

    residual_chart_id = f"residual-{slug}"
    residual_rows = exact["residual_rows"]
    charts.append(
        _payload_line(
            residual_chart_id,
            f"Band residuals: {site}",
            residual_rows,
            x="aod",
            y="value",
            series="series",
        )
    )
    residual_figure = _chart_figure(
        chart_id=residual_chart_id,
        title="Per-band BOA mismatch versus AOD",
        subtitle="Median absolute |atmospherically corrected BOA - surface prior| in reflectance",
        fallback=_line_svg(
            residual_rows,
            x_field="aod",
            y_field="value",
            series_field="series",
            aria_label=f"Per-band BOA residual curves for {site}",
            markers=(
                (case["truth"], "AERONET", COLORS["orange"]),
                (float(current["retrieved"]), "retrieved", COLORS["blue_dark"]),
            ),
        ),
        note="The exact dump stores residual magnitude, not sign. Signed truth-node residuals in the table are a proxy reconstructed from the earlier calibration capture.",
        datasets=cube_dataset + "; " + str(case["sources"].get("calibration_dump") or "no calibration dump"),
    )

    history_rows = case["clean_prior"].get("history_rows", [])
    history_chart = ""
    if history_rows:
        history_long = []
        for row in history_rows:
            for band in ("B02", "B04"):
                value = _finite(row.get(band))
                if value is not None:
                    history_long.append({"window": str(row["window"]), "series": band, "reflectance": value})
        history_chart_id = f"history-{slug}"
        charts.append(
            _payload_bar(
                history_chart_id,
                f"Clean-day prior history: {site}",
                history_long,
                category="window",
                value="reflectance",
                series="series",
            )
        )
        history_chart = _chart_figure(
            chart_id=history_chart_id,
            title="Selected clean-day surface history",
            subtitle="Station-window median BOA reflectance for every retained historical window",
            fallback=_grouped_bar_svg(
                history_long,
                category_field="window",
                series_field="series",
                value_field="reflectance",
                aria_label=f"Clean-day prior history for {site}",
            ),
            note="Bars expose whether a small temporal spread is genuine stability or repeated selection of similarly contaminated observations.",
            datasets=str(case["sources"].get("clean_prior")),
        )

    aerosol_rows = []
    for row in case["aerosol"].get("rows", []):
        for model in ("CAMS", "MERRA-2"):
            value = _finite(row.get(model))
            if value is not None:
                aerosol_rows.append({"species": row["species"], "series": model, "fraction": value})
    aerosol_chart = ""
    if aerosol_rows:
        aerosol_chart_id = f"aerosol-{slug}"
        charts.append(
            _payload_bar(
                aerosol_chart_id,
                f"Aerosol composition: {site}",
                aerosol_rows,
                category="species",
                value="fraction",
                series="series",
                percent=True,
            )
        )
        aerosol_chart = _chart_figure(
            chart_id=aerosol_chart_id,
            title="Independent aerosol composition context",
            subtitle="CAMS and MERRA-2 species fractions at the target time; diagnostic context only",
            fallback=_grouped_bar_svg(
                aerosol_rows,
                category_field="species",
                series_field="series",
                value_field="fraction",
                aria_label=f"Aerosol composition context for {site}",
                percent=True,
            ),
            note="These products are not used by the reported surface-only objective. Agreement on organic matter or black carbon strengthens, but does not prove, the smoke-optics diagnosis.",
            datasets=case["sources"]["aerosol_context"],
        )

    prior = case["clean_prior"]
    prior_dataset = str(case["sources"].get("clean_prior") or "clean prior unavailable")
    prior_summary = (
        f"The strict L1C prior retains {tips.wrap(str(prior.get('realizations', 0)), prior_dataset)} historical windows, "
        f"uses {tips.wrap(str(prior.get('fallback_scenes', 0)), prior_dataset)} fallback scenes, and has median selected-AOD excess above 0.15 of "
        f"{tips.wrap(_fmt(_finite(prior.get('median_selected_aod_excess'))), prior_dataset)}. "
        "Temporal spread is therefore a repeatability measure, not an independently calibrated surface-error sigma."
        if prior.get("available")
        else "No strict clean-day L1C prior archive is available for this case."
    )

    matchup = case["matchup"]
    input_rows = [
        ("AERONET AOD550", truth_tip(), "Reference used only for evaluation and diagnostics"),
        ("Current retrieved AOD", retrieved_tip(), "Scene-mean extraction from the pooled surface solve"),
        ("Expected-error half-width", ee_tip(), "0.05 + 0.15 x AERONET AOD"),
        ("Final cost per band", cost_tip(), "Large values indicate no low-residual fit"),
        ("L2A TCWV", wvp_html + " cm", "Same-scene Sen2Cor WVP, unscaled"),
        ("Computed cloud fraction", cloud_tip(), "Ignored by this run; useful as plume/cloud context"),
        (
            "Valid solve pixels",
            tips.wrap(
                f"{exact['solve_valid_count']}/{exact['solve_total_count']}",
                cube_dataset,
            ),
            "Cost cube pixels entering the solve",
        ),
        (
            "AERONET samples",
            tips.wrap(str(matchup.get("n_aeronet", "n/a")), case["sources"]["matchup"]),
            "Samples in the scene matchup window",
        ),
        (
            "Mean time offset",
            tips.wrap(f"{float(matchup.get('mean_abs_time_offset_min', 'nan')):.1f} min", case["sources"]["matchup"]),
            "Absolute AERONET-to-scene time offset",
        ),
    ]
    input_table = "".join(
        f"<tr><td>{_esc(label)}</td><td>{value}</td><td>{_esc(note)}</td></tr>"
        for label, value, note in input_rows
    )

    map_datasets = "; ".join(
        value for value in (cube_dataset, case["sources"].get("clean_prior"), case["sources"].get("calibration_dump")) if value
    )
    return f"""
      <article class="case-section" id="case-{slug}">
        <header class="case-header">
          <div><div class="case-index">{_esc(case['matchup_id'])}</div><h2>{_esc(site)}: {_esc(diagnosis['primary'])}</h2></div>
          <a class="back-link" href="#case-navigation" aria-label="Back to case navigation">Back to cases</a>
        </header>
        <p class="case-lead">{opening}{species_sentence}{wvp_sentence}</p>

        <section class="metrics case-metrics">
          <div class="metric"><div class="metric-label">AERONET</div><div class="metric-value">{truth_tip()}</div><div class="metric-note">AOD550 truth</div></div>
          <div class="metric"><div class="metric-label">Current solve</div><div class="metric-value">{retrieved_tip()}</div><div class="metric-note">outside EE +/- {ee_tip()}</div></div>
          <div class="metric"><div class="metric-label">Cost / band</div><div class="metric-value">{cost_tip()}</div><div class="metric-note">final surface objective</div></div>
          <div class="metric"><div class="metric-label">Cloud mask</div><div class="metric-value">{cloud_tip()}</div><div class="metric-note">not used as a solve exclusion</div></div>
        </section>

        {_image_figure(image_base64=mosaic, site=site, datasets=map_datasets)}

        <section class="case-copy">
          <h3>Input and solver state</h3>
          <p>The exact input record shows whether the failure entered through atmospheric state, mask interpretation, extraction, or the cost itself. The cloud fraction is deliberately retained even though this run sets the cloud/water ignore switch.</p>
        </section>
        <section class="card table-card compact-table"><div class="table-scroll"><table><thead><tr><th>Quantity</th><th>Value</th><th>Diagnostic meaning</th></tr></thead><tbody>{input_table}</tbody></table></div></section>

        <section class="case-copy"><h3>NIR and SWIR anchors do not rescue the visible closure</h3><p>Target TOA is observed reflectance from the same scene. BOA at retrieved and truth AOD is reconstructed from the older calibration capture, while the clean prior is the strict L1C historical composite. B8A, B11, and B12 make the water-vapour-sensitive anchor behavior inspectable even though the exact current cost dump retains only B02 and B04.</p></section>
        <section class="card table-card"><div class="card-head"><h3>Target and clean-prior spectrum</h3><p>Station-window medians; reconstructed BOA columns are calibration-capture proxies</p></div><div class="table-scroll"><table><thead><tr><th>Band</th><th>Wavelength (nm)</th><th>Target TOA</th><th>BOA @ retrieved</th><th>BOA @ truth</th><th>Truth - retrieved BOA</th><th>Strict clean prior</th></tr></thead><tbody>{_spectral_table(case, tips)}</tbody></table></div></section>

        <section class="case-copy"><h3>The local cost shape identifies the failure before AERONET scoring</h3><p>The production result can be compared with the raw local cost basin. A lower-edge minimum with a large absolute cost is an invalid physical fit, whereas a broad internal minimum with opposing band preferences is an aggregation/weighting problem.</p></section>
        {cost_figure}
        {residual_figure}

        <section class="case-copy"><h3>TOA, surface prior, and truth-node closure</h3><p>At the AERONET node, the standardized residual shows how many stated prior sigmas separate the corrected target from the prior. Positive signed proxy residual means the older calibration capture's predicted prior is darker than BOA reconstructed at truth. All available proxy values are positive here; the sign is informative but is not an exact field from the L2A-WVP cube.</p></section>
        <section class="card table-card"><div class="card-head"><h3>Band-level closure audit</h3><p>Station-window medians; reflectance except AOD and standardized residual</p></div><div class="table-scroll"><table><thead><tr><th>Band</th><th>TOA</th><th>BOA prior</th><th>Prior sigma</th><th>Band argmin AOD</th><th>|residual| @ retrieved</th><th>|residual| @ truth</th><th>|z| @ truth</th><th>Signed truth - prior proxy</th></tr></thead><tbody>{_band_table(case, tips)}</tbody></table></div></section>

        <section class="case-copy"><h3>The clean-day history must be judged by representativeness, not spread alone</h3><p>{prior_summary}</p></section>
        {history_chart}
        <details class="audit-details"><summary>Inspect every retained clean-day window</summary><div class="card table-card"><div class="table-scroll"><table><thead><tr><th>Window</th><th>B02</th><th>B04</th><th>Selected AOD median</th><th>Selected AOD max</th><th>Fallback scenes</th></tr></thead><tbody>{_history_table(case, tips)}</tbody></table></div></div></details>

        <section class="case-copy"><h3>Aerosol optics explain part, but not all, of the remaining gap</h3><p>{_esc(diagnosis['species_result'])}. The method table keeps surface solves separate from external aerosol context so a good external estimate is not mistaken for a surface-prior retrieval success.</p></section>
        {aerosol_chart}
        <section class="card table-card"><div class="card-head"><h3>All saved estimates for this scene</h3><p>Within-EE status uses the same AERONET truth and threshold for comparability</p></div><div class="table-scroll"><table><thead><tr><th>Method or product</th><th>Role</th><th>AOD550</th><th>Error</th><th>Status</th></tr></thead><tbody>{_method_table(case, tips)}</tbody></table></div></section>

        <div class="case-conclusion"><strong>Most useful next check.</strong> {_esc(_case_next_check(case))}</div>
      </article>
    """


def _case_next_check(case: dict[str, Any]) -> str:
    diagnosis = case["diagnosis"]
    if diagnosis["rail_collapse"] and _within_ee(_finite(case["biomass"].get("retrieved")), case["truth"]):
        return "Freeze the biomass optical family for this scene, rebuild the clean prior in the same 6S species space, and verify signed B02/B04 closure at truth before changing solver thresholds."
    if diagnosis["rail_collapse"]:
        return "Dump signed residuals for the exact run and test a small absorbing-aerosol family grid. If every family still rails low, the surface prior is too bright or unrepresentative and must be rebuilt before solver tuning."
    if diagnosis["band_conflict"]:
        return "Use the aerosol family that closes both bands, then compare local cost-field median, station extraction, and scene mean. This case should not be treated as a generic lower-rail failure."
    return "Inspect signed band closure and prior contamination first; the current scalar cost does not isolate a single cause."


def _overview_method_rows(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    wanted = {
        "Current surface solve (L2A WVP)": "Current",
        "acixThree sdspec absolute": "sdspec abs",
        "Prescribed biomass 6S": "Biomass",
        "Sentinel-2 L2A AOT": "S2 L2A",
    }
    rows = []
    for case in cases:
        rows.append({"site": case["display_site"], "series": "AERONET", "aod": case["truth"]})
        for method in case["methods"]:
            label = wanted.get(method["method"])
            if label and method["aod"] is not None:
                rows.append({"site": case["display_site"], "series": label, "aod": method["aod"]})
    return rows


def _overview_table(cases: list[dict[str, Any]], tips: Tooltips) -> str:
    rows = []
    for case in cases:
        solver = case["current"].get("solver", {})
        biomass = _finite(case["biomass"].get("retrieved"))
        sdspec_abs = _finite(case["sdspec_abs"].get("aod"))
        dataset = case["sources"]["current_result"]
        biomass_dataset = _source_label(
            BIOMASS_RESULTS / f"{case['matchup_id']}.json"
        )
        sdspec_dataset = _source_label(
            SDSPEC_RESULTS
            / f"{case['matchup_id']}__seas_maiac_sdspec_highaod6_clean015_abs.json"
        )
        rows.append(
            "<tr>"
            f"<td><a href=\"#case-{re.sub(r'[^a-z0-9]+', '-', case['display_site'].lower()).strip('-')}\">{_esc(case['display_site'])}</a></td>"
            f"<td>{tips.wrap(_fmt(case['truth']), dataset)}</td>"
            f"<td>{tips.wrap(_fmt(_finite(case['current'].get('retrieved'))), dataset)}</td>"
            f"<td>{tips.wrap(_fmt(_finite(solver.get('cost_final_per_band')), 1), dataset)}</td>"
            f"<td>{tips.wrap(_fmt(_finite(solver.get('surface_band_B02_argmin_aot'))), dataset)}</td>"
            f"<td>{tips.wrap(_fmt(_finite(solver.get('surface_band_B04_argmin_aot'))), dataset)}</td>"
            f"<td>{tips.wrap(_fmt(biomass), biomass_dataset)}</td>"
            f"<td>{tips.wrap(_fmt(sdspec_abs), sdspec_dataset)}</td>"
            f"<td>{_esc(case['diagnosis']['primary'])}</td>"
            "</tr>"
        )
    return "".join(rows)


def _build_html(
    *,
    out: Path,
    cases: list[dict[str, Any]],
    mosaics: dict[str, str],
    campaign: dict[str, Any],
) -> None:
    tips = Tooltips()
    charts: list[dict[str, Any]] = []
    overview_rows = _overview_method_rows(cases)
    charts.append(
        _payload_bar(
            "overview-aod",
            "Saved AOD estimates for the six cases",
            overview_rows,
            category="site",
            value="aod",
            series="series",
        )
    )
    overview_figure = _chart_figure(
        chart_id="overview-aod",
        title="AOD estimates across the six failures",
        subtitle="Current solve, selected comparators, and AERONET truth; one scene per site",
        fallback=_grouped_bar_svg(
            overview_rows,
            category_field="site",
            series_field="series",
            value_field="aod",
            aria_label="AOD estimates across the six high-AOD failures",
        ),
        note="External S2 L2A AOT is context, not part of the surface-prior solver. Comparator availability is shown exactly rather than imputed.",
        datasets="current high-AOD results; sdspec absolute; biomass profile; S2 L2A AOT; AERONET matchup truth",
    )

    cloud_rows = [row for row in campaign["rows"] if row["cloud_frac"] is not None]
    charts.append(
        _payload_scatter(
            "cloud-error",
            "Cloud fraction and absolute AOD error",
            cloud_rows,
            x="cloud_frac",
            y="abs_error",
            group="group",
        )
    )
    cloud_figure = _chart_figure(
        chart_id="cloud-error",
        title="Cloud fraction versus absolute AOD error",
        subtitle=f"Complete R2 campaign context, n={campaign['cloud_count']}; the focal six are highlighted",
        fallback=_scatter_svg(
            cloud_rows,
            x_field="cloud_frac",
            y_field="abs_error",
            group_field="group",
            aria_label="Cloud fraction versus absolute AOD error in campaign context",
        ),
        note="Cloud fraction is associated with thick aerosol and plume-like scenes, but it is not a strong monotonic explanation of error and was ignored as a solve exclusion in the focal run.",
        datasets=_source_label(R2_CAMPAIGN_RESULTS),
    )

    cost_rows = [row for row in campaign["rows"] if row["log_cost"] is not None]
    charts.append(
        _payload_scatter(
            "cost-error",
            "Final cost and absolute AOD error",
            cost_rows,
            x="log_cost",
            y="abs_error",
            group="group",
        )
    )
    cost_figure = _chart_figure(
        chart_id="cost-error",
        title="Final cost versus absolute AOD error",
        subtitle=f"x = log10(1 + cost per band), complete R2 campaign context, n={campaign['cost_count']}",
        fallback=_scatter_svg(
            cost_rows,
            x_field="log_cost",
            y_field="abs_error",
            group_field="group",
            aria_label="Final cost versus absolute AOD error in campaign context",
        ),
        note="High cost is a useful invalid-fit warning for the catastrophic cases, but low cost does not guarantee correct AOD elsewhere in the campaign.",
        datasets=_source_label(R2_CAMPAIGN_RESULTS),
    )

    current_hits = sum(bool(case["current"].get("within_ee")) for case in cases)
    rail_count = sum(bool(case["diagnosis"]["rail_collapse"]) for case in cases)
    conflict_count = sum(bool(case["diagnosis"]["band_conflict"]) for case in cases)
    biomass_hits = sum(
        _within_ee(_finite(case["biomass"].get("retrieved")), case["truth"]) is True
        for case in cases
    )
    sdspec_abs_hits = sum(
        _within_ee(_finite(case["sdspec_abs"].get("aod")), case["truth"]) is True
        for case in cases
    )
    wvp_deltas = [
        float(case["diagnosis"]["wvp_retrieval_delta"])
        for case in cases
        if case["diagnosis"]["wvp_retrieval_delta"] is not None
    ]
    max_wvp_delta = max(wvp_deltas, default=math.nan)
    rail_wvp_deltas = [
        float(case["diagnosis"]["wvp_retrieval_delta"])
        for case in cases
        if case["diagnosis"]["rail_collapse"]
        and case["diagnosis"]["wvp_retrieval_delta"] is not None
    ]
    max_rail_wvp_delta = max(rail_wvp_deltas, default=math.nan)
    truth_z = [
        float(band["z_at_truth"])
        for case in cases
        for band in case["exact"]["bands"]
        if band["z_at_truth"] is not None
    ]
    signed_proxy = [
        float(case["signed_closure_proxy"]["bands"][band]["truth_minus_prior"])
        for case in cases
        for band in ("B02", "B04")
        if (case["signed_closure_proxy"].get("bands", {}).get(band) or {}).get(
            "truth_minus_prior"
        )
        is not None
    ]
    truth_z_min = min(truth_z, default=math.nan)
    truth_z_max = max(truth_z, default=math.nan)
    signed_min = min(signed_proxy, default=math.nan)
    signed_max = max(signed_proxy, default=math.nan)
    result_dataset = _source_label(CURRENT_RESULTS)
    cloud_corr = campaign["cloud_pearson"]
    cost_corr = campaign["cost_pearson"]
    generated = datetime.now().astimezone().strftime("%d %b %Y, %H:%M %Z")

    nav = "".join(
        f'<a href="#case-{re.sub(r"[^a-z0-9]+", "-", case["display_site"].lower()).strip("-")}">{_esc(case["display_site"])}</a>'
        for case in cases
    )
    case_html = "".join(
        _case_section(case, mosaics[case["matchup_id"]], tips, charts) for case in cases
    )

    main = f"""
    <main data-report-audience="technical">
      <article class="reading">
        <div class="kicker">SIAC high-AOD failure analysis</div>
        <header data-contract-section="title"><h1>Why the six high-AOD surface solves fail</h1></header>
        <section class="summary technical-summary" data-contract-section="technical-summary">
          <div class="summary-label">Technical Summary</div>
          <div class="summary-body">
            <p><strong>The exact L2A-WVP rerun does not solve any of the six cases:</strong> {tips.wrap(f"{current_hits}/6", result_dataset)} are within expected error. Five cases are not ordinary optimizer failures; both visible bands make the physical cost decrease toward the {tips.wrap("0.01", _source_label(CURRENT_CUBES))} lower AOD rail. Mont Joli is different: its aggregate minimum is close to truth, but B02 and B04 prefer opposite ends of the axis and the spatial extraction under-retrieves.</p>
            <p><strong>Water vapour is physically corrected but not the controlling error.</strong> Same-scene L2A TCWV changes retrieved AOD by at most {tips.wrap(_fmt(max_rail_wvp_delta), result_dataset + "; " + _source_label(PRIOR_WVP_RESULTS))} across the five rail-collapse cases. Mont Joli changes by {tips.wrap(_fmt(max_wvp_delta), result_dataset + "; " + _source_label(PRIOR_WVP_RESULTS))}, but remains outside EE. The decisive evidence is cost shape, truth-node BOA mismatch, and aerosol-family response.</p>
            <p><strong>The prior/RT combination is grossly inconsistent at high AOD.</strong> Exact truth-node residual magnitudes span {tips.wrap(f"{truth_z_min:.1f}-{truth_z_max:.1f} prior sigmas", _source_label(CURRENT_CUBES))}. The independent signed proxy is positive in both bands for all six cases, with truth-minus-prior reflectance of {tips.wrap(f"{signed_min:.3f}-{signed_max:.3f}", _source_label(CALIB_DUMPS))}. Under the saved continental 6S coefficients, the target corrected at truth is much brighter than the predicted surface. That is consistent with a surface prior that is too dark/unrepresentative, insufficient path-radiance removal from the aerosol optical model, or both.</p>
            <p><strong>Aerosol optics are necessary but not sufficient.</strong> A prescribed biomass profile recovers {tips.wrap(f"{biomass_hits}/6", _source_label(BIOMASS_RESULTS))}; the acixThree-style absolute-cost branch recovers {tips.wrap(f"{sdspec_abs_hits}/6", _source_label(SDSPEC_RESULTS))}. The remaining cases still require a better calibrated clean-day prior and an invalid-fit policy that does not report lower-rail AOD as a confident solution.</p>
          </div>
        </section>
        <section class="metrics">
          <div class="metric"><div class="metric-label">Current EE</div><div class="metric-value">{tips.wrap(f"{current_hits}/6", result_dataset)}</div><div class="metric-note">L2A-WVP exact rerun</div></div>
          <div class="metric"><div class="metric-label">Two-band rail collapse</div><div class="metric-value">{tips.wrap(f"{rail_count}/6", result_dataset)}</div><div class="metric-note">both bands prefer <= 0.05</div></div>
          <div class="metric"><div class="metric-label">Large band conflict</div><div class="metric-value">{tips.wrap(f"{conflict_count}/6", result_dataset)}</div><div class="metric-note">argmin spread >= 0.5</div></div>
          <div class="metric"><div class="metric-label">Biomass recovery</div><div class="metric-value">{tips.wrap(f"{biomass_hits}/6", _source_label(BIOMASS_RESULTS))}</div><div class="metric-note">case-level diagnostic, not selector score</div></div>
        </section>
      </article>

      <article data-contract-section="key-findings">
        <section class="reading narrative"><h2>The six failures split into two mechanisms</h2><p>Alta Floresta, London, Metsi, NEON DEJU, and Waskesiu have a coherent lower-rail preference in both bands, usually with very high residual cost. Mont Joli has useful high-AOD information, but it is split across bands and diluted by the scene-level aggregation. Treating all six with one fallback threshold obscures that distinction.</p></section>
        {overview_figure}
        <section class="reading card table-card"><div class="card-head"><h3>Cross-case failure matrix</h3><p>Exact current solve plus two aerosol/cost comparators</p></div><div class="table-scroll"><table><thead><tr><th>Case</th><th>Truth</th><th>Current</th><th>Cost / band</th><th>B02 min</th><th>B04 min</th><th>Biomass</th><th>sdspec abs</th><th>Primary diagnosis</th></tr></thead><tbody>{_overview_table(cases, tips)}</tbody></table></div></section>

        <section class="reading narrative"><h2>Cloud saturation marks difficult scenes but does not explain the solve</h2><p>Across the complete campaign context, cloud fraction has Pearson correlation {tips.wrap(_fmt(cloud_corr), _source_label(R2_CAMPAIGN_RESULTS))} with absolute AOD error, compared with {tips.wrap(_fmt(cost_corr), _source_label(R2_CAMPAIGN_RESULTS))} for log final cost. More importantly, the focal run explicitly ignores cloud/water as a solve exclusion, and the exact cubes contain valid cost pixels. A cloud fraction near one is therefore plume/context evidence, not proof that masking forced AOD to zero.</p></section>
        {cloud_figure}
        {cost_figure}

        <section class="reading case-navigation" id="case-navigation"><h2>Inspect each case</h2><p>Each case follows the same chain: target input, clean-day prior, spatial cost field, per-band closure, historical prior selection, and aerosol-family response.</p><nav aria-label="High-AOD case navigation">{nav}</nav></section>
        {case_html}
      </article>

      <article class="reading">
        <section class="narrative" data-contract-section="scope-data-and-metric-definitions">
          <h2>Scope and metric definitions</h2>
          <p><strong>Population:</strong> six campaign-250 scenes with AERONET AOD550 above one that remained outside expected error after the backstop-escape rerun. <strong>Current method:</strong> 6S continental aerosol, B02+B04 surface solve, tau-dependent prior enabled, same-scene S2 L2A WVP, and scene-mean extraction. <strong>Success:</strong> |retrieved - AERONET| <= 0.05 + 0.15 x AERONET. <strong>Local diagnostics:</strong> medians inside a 1.5 km station radius over valid solve pixels.</p>
        </section>

        <section class="narrative" data-contract-section="methodology">
          <h2>What the report reconstructs</h2>
          <p>The exact cost archives provide the AOD axis, shape and absolute cost cubes, per-band normalized costs, absolute BOA residuals, TOA B02/B04, BOA prior, prior uncertainty, aerosol backstop, coordinates, and valid-pixel mask. The report independently recomputes local cost curves and raw per-pixel argmins. Strict L1C clean-day archives provide historical BOA realizations and selection metadata. Calibration captures provide target NIR/red/blue TOA and a signed closure proxy by evaluating saved 6S coefficient curves at AERONET truth.</p>
          <p>Comparator AODs are read from saved prior-WVP, R2, L1C-clean, L2A-PC, sdspec, biomass, Sen2Cor, MAIAC, CAMS, and MERRA-2 artifacts. No retrieval, RT call, Earth Engine query, fitting, or AERONET-informed parameter selection is performed while building this page.</p>
        </section>

        <section class="narrative" data-contract-section="limitations-uncertainty-and-robustness-checks">
          <h2>Limits and robustness checks</h2>
          <p>The exact current cube stores absolute band residuals, so its residual direction cannot be recovered. Signed truth-minus-prior values come from an earlier calibration capture with saved 6S coefficients and must be treated as directional evidence, not an exact reproduction of the L2A-WVP run. This distinction is shown in every band table.</p>
          <p>The six cases were selected because they failed at high AOD and are not a performance sample. Biomass and sdspec recovery counts diagnose mechanism only; they do not validate a deployable scene selector. External aerosol products are contextual and never counted as surface-solve evidence. AERONET truth is used to mark diagnostic nodes and score outputs, never to choose a runtime solution.</p>
          <p>The cloud/error and cost/error correlations are descriptive. They show association across the complete saved R2 campaign and do not establish causality. Spatial raw argmins are shown before the production median-pool operation, which is why they need not equal the final scene mean.</p>
        </section>

        <section class="narrative" data-contract-section="recommended-next-steps">
          <h2>Run the next experiments by failure mechanism</h2>
          <ol class="method-list">
            <li><strong>Make invalid fits explicit.</strong> A two-band lower-edge minimum with high normalized cost should return a low-confidence/invalid state, not a nominal AOD of 0.01. Dump signed residuals in the exact cost archive so prior direction is auditable.</li>
            <li><strong>For the five rail-collapse cases, test optics before solver thresholds.</strong> Run a small frozen 6S family grid spanning continental, biomass, dust, and an absorbing fine-mode variant while keeping the exact same L1C prior and L2A WVP. Compare truth-node signed B02/B04 closure and penalize family complexity.</li>
            <li><strong>For Mont Joli, isolate aggregation.</strong> Keep the biomass-compatible family, then compare station, 1.5 km integrated cost field, and scene mean. Add a band-conflict gate that preserves an internally supported high-AOD basin rather than averaging incompatible pixel/band solutions.</li>
            <li><strong>Promote the surface prior only after BOA validation.</strong> Require leave-year-out L1C clean-day priors to reduce held-out B02/B04 median error and produce calibrated one- and two-sigma coverage. Penalize fallback scenes and residual aerosol contamination; match geometry and phenology before adding more solver freedom.</li>
            <li><strong>Then freeze one generic campaign run.</strong> Use raw same-scene L2A WVP, the selected 6S optics policy, signed/covariance-aware visible likelihood, and an explicit invalid state. Score once on campaign-250 and once on a site-grouped external holdout.</li>
          </ol>
        </section>

        <section class="narrative" data-contract-section="further-questions">
          <h2>Questions the current dumps cannot answer</h2>
          <ul class="method-list">
            <li>Does an exact signed L2A-WVP dump confirm that the predicted surface is too dark at truth, or does changing aerosol absorption/path radiance reverse that closure error?</li>
            <li>Can an aerosol family be selected from target spectral closure and independent composition context without using AERONET or an external AOD as a runtime answer?</li>
            <li>How much of Mont Joli's miss is caused by scene-mean extraction versus the pooled cost-field solve?</li>
            <li>Are the retained clean-day windows geometrically and phenologically representative even when their MAIAC AOD passes the strict threshold?</li>
          </ul>
        </section>
        <section class="caveat"><strong>Bottom line.</strong> The zero-like AODs are not numerical zeros and are not caused by the L2A water-vapour input. They are lower-bound solutions selected because the current surface-prior/continental-6S combination makes the truth-AOD target far brighter than the predicted surface. The appropriate fix is a brighter or more representative calibrated prior plus aerosol-optics and invalid-fit handling, with a separate aggregation correction for Mont Joli.</section>
      </article>
    </main>
    """

    source = SHELL.read_text(encoding="utf-8")
    source = re.sub(
        r"\s*<!-- Deliver this HTML file itself\..*?-->\s*",
        "\n",
        source,
        flags=re.DOTALL,
    )
    source = source.replace("{{TITLE}}", "Why the six high-AOD surface solves fail")
    source = source.replace("{{SOURCE_AND_DATE}}", f"Saved SIAC artifacts | {generated}")
    source = source.replace("{{REPORT_AUDIENCE}}", "technical")
    source = source.replace(
        '<div class="brand"><span class="mark" aria-hidden="true"></span>Data Analytics</div>',
        '<div class="brand"><span class="mark" aria-hidden="true"></span>SIAC diagnostics</div>',
    )
    source = re.sub(
        r'<main data-report-audience="[^"]+">.*?</main>',
        main,
        source,
        flags=re.DOTALL,
    )
    extra_css = """
    h1, h2, .metric-value { letter-spacing: 0; }
    .shell, .card, .metric { border-radius: 8px; }
    .mark { border-radius: 4px; background: var(--blue); }
    .technical-summary { margin-top: 28px; }
    .case-navigation { margin: 54px auto 18px; padding-top: 34px; border-top: 1px solid var(--border-strong); }
    .case-navigation p { color: var(--secondary); }
    .case-navigation nav { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 18px; }
    .case-navigation a, .back-link { display: inline-flex; align-items: center; min-height: 34px; padding: 6px 10px; border: 1px solid var(--border-strong); border-radius: 5px; color: var(--text); text-decoration: none; }
    .case-navigation a:hover, .back-link:hover { border-color: var(--blue); color: var(--blue); }
    .case-section { width: min(1040px, 100%); margin: 0 auto; padding: 64px 0 18px; border-top: 1px solid var(--border-strong); scroll-margin-top: 16px; }
    .case-header { display: flex; align-items: flex-start; justify-content: space-between; gap: 18px; width: min(800px, 100%); margin: 0 auto; }
    .case-index { margin-bottom: 6px; color: var(--muted); font: 11px/1.4 ui-monospace, SFMono-Regular, Menlo, monospace; overflow-wrap: anywhere; }
    .case-lead, .case-copy { width: min(800px, 100%); margin: 18px auto 22px; color: var(--secondary); }
    .case-copy h3 { margin-bottom: 6px; color: var(--text); font-size: 17px; }
    .case-copy p { margin: 0; }
    .case-metrics { width: min(960px, 100%); margin: 24px auto 32px; }
    .case-conclusion { width: min(800px, 100%); margin: 32px auto 8px; padding: 16px 18px; border-left: 3px solid var(--gold); background: var(--surface-tertiary); color: var(--secondary); }
    .case-conclusion strong { color: var(--text); }
    .map-figure { margin: 0 auto 36px; }
    .map-image { padding: 14px; background: #fff; }
    .map-image img { display: block; width: 100%; height: auto; }
    .compact-table { margin: 12px auto 36px; }
    .compact-table td:nth-child(3), .table-card td:last-child { white-space: normal; text-align: left; min-width: 240px; }
    .table-card th:first-child, .table-card td:first-child { position: sticky; left: 0; z-index: 1; background: var(--surface); }
    .audit-details { width: 100%; margin: 12px auto 38px; }
    .audit-details summary { width: min(800px, 100%); margin: 0 auto; cursor: pointer; color: var(--blue); font-weight: 600; }
    .audit-details[open] summary { margin-bottom: 14px; }
    .method-list { margin: 12px 0 0; padding-left: 22px; color: var(--secondary); }
    .method-list li { margin-bottom: 12px; }
    .reading.card { width: min(1040px, 100%); }
    @media (max-width: 800px) {
      .case-section { padding: 48px 18px 12px; }
      .case-header { flex-direction: column; }
      .chart-fallback svg { min-width: 760px; }
      .map-image { overflow-x: auto; }
      .map-image img { min-width: 760px; }
      .compact-table td:nth-child(3), .table-card td:last-child { min-width: 210px; }
    }
    """
    source = source.replace("  </style>", extra_css + "  </style>")
    shell_path = out / "report-shell.html"
    payload_path = out / "report-payload.json"
    shell_path.write_text(source, encoding="utf-8")
    payload_path.write_text(json.dumps(_jsonable({"charts": charts}), indent=2), encoding="utf-8")
    subprocess.run(
        [
            sys.executable,
            str(EMBED),
            "--input",
            str(shell_path),
            "--payload",
            str(payload_path),
            "--output",
            str(out / "report.html"),
        ],
        check=True,
    )

    chart_map = []
    for chart in charts:
        spec = chart["dataset"]["chart_spec"]
        chart_map.append(
            {
                "id": chart["id"],
                "title": chart["dataset"]["title"],
                "family": spec["type"],
                "fields": spec["encodings"],
                "palette_policy": "hard two-root cap for relationships; relaxed approved roots for method comparisons",
                "delivery": "packaged Recharts with same-data inline SVG fallback",
            }
        )
    source_notes = {
        "audience": "technical",
        "delivery_mode": "html",
        "question": "Why do the six high-AOD surface-prior solves fail, and which input, prior, cost, TOA, solver, or aerosol component should be changed next?",
        "required_section_map": [
            "title",
            "technical-summary",
            "key-findings",
            "scope-data-and-metric-definitions",
            "methodology",
            "limitations-uncertainty-and-robustness-checks",
            "recommended-next-steps",
            "further-questions",
        ],
        "sources": [
            str(CURRENT_RESULTS),
            str(CURRENT_CUBES),
            str(PRIOR_WVP_RESULTS),
            str(R2_RESULTS),
            str(R2_CAMPAIGN_RESULTS),
            str(L1C_RESULTS),
            str(L2A_PC_RESULTS),
            str(BIOMASS_RESULTS),
            str(SDSPEC_RESULTS),
            str(L1C_PRIORS),
            str(CALIB_DUMPS),
            str(AEROSOL_CONTEXT),
            str(S2_L2A_AOT),
            str(S2_L2A_WVP),
            str(MAIAC),
            str(MATCHUPS),
        ],
        "chart_map": chart_map,
        "omissions": {
            "exact_signed_residual": "The current cost archive stores absolute residuals only; signed values are explicitly labeled as an older calibration-capture proxy.",
            "cloud_causality": "Correlation and run configuration support a diagnostic association only, not a causal cloud claim.",
            "species_selector_score": "Case-level fixed-family recovery is not presented as a deployable selector or campaign score.",
        },
        "qa": {
            "structural": "Generator enforces the technical section map, one chart host per payload id, same-data SVG fallbacks, and source tooltips.",
            "visual": "Representative same-data SVG fallbacks and spatial mosaics were rasterized and inspected. No browser executable was available for live Recharts or full responsive-page screenshots.",
        },
    }
    (out / "source-notes.json").write_text(json.dumps(source_notes, indent=2), encoding="utf-8")


def build(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    matchups = _matchup_rows()
    cases = []
    mosaics: dict[str, str] = {}
    for index, mid in enumerate(MIDS, start=1):
        print(f"[{index}/{len(MIDS)}] analysing {mid}", flush=True)
        case, mosaic = _case_record(mid, matchups[mid])
        cases.append(case)
        mosaics[mid] = mosaic
    campaign = _campaign_context()
    analysis = {
        "description": "Saved-artifact diagnosis of six high-AOD surface-driven AOD failures",
        "generated_at": datetime.now().astimezone().isoformat(),
        "cases": cases,
        "campaign_context": campaign,
    }
    (out / "analysis.json").write_text(
        json.dumps(_jsonable(analysis), indent=2), encoding="utf-8"
    )
    _build_html(out=out, cases=cases, mosaics=mosaics, campaign=campaign)
    for path in out.iterdir():
        if path.is_file():
            path.chmod(0o644)
    out.chmod(0o755)
    print(out / "report.html", flush=True)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    return parser.parse_args()


if __name__ == "__main__":
    build(_parse_args().output_dir)
