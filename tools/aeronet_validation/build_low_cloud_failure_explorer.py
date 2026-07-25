"""Build an interactive spatial/spectral/cost explorer for low-cloud AOD failures."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from pyproj import Transformer

if TYPE_CHECKING:
    from collections.abc import Sequence

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.aeronet_validation.build_low_cloud_aod_report import (  # noqa: E402
    SINGLE_PRIOR_LABELS,
)
from tools.aeronet_validation.summarize_single_prior_adaptive import OPTIONS  # noqa: E402

from siac._rust_compat import surface_driven_pool_argmin  # noqa: E402

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
CURRENT_RESULTS = ROOT / "phaseD_results_lowcloud20_singleprior_b03_chi2_20260711"
DIAGNOSTIC_RESULTS = ROOT / "phaseD_results_lowcloud20_failurediag_b03_chi2_20260712"
COST_CUBES = ROOT / "phaseD_cost_cubes_lowcloud20_failurediag_b03_chi2_20260712"
BASELINE_RESULTS = ROOT / "phaseD_results_campaign250_R2_full_localdiag_20260705"
PRIOR_RESULTS = ROOT / "prior_quality"
DEFAULT_OUTPUT = ROOT / "reports/aod-low-cloud-failure-explorer-20260712"
WEB_ASSETS = Path(__file__).with_name("low_cloud_failure_explorer")

BANDS = ("B02", "B03", "B04")
WAVELENGTHS = {"B02": 490.0, "B03": 560.0, "B04": 665.0}
BAND_COLORS = {"B02": "#1565c0", "B03": "#2e7d32", "B04": "#c62828"}
GROUP_COLORS = {
    "A": "#1769aa",
    "B": "#6d4c41",
    "C": "#7b1fa2",
    "D": "#00838f",
    "E": "#ad6800",
    "F": "#546e7a",
    "G": "#9e3d64",
}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_records(directory: Path) -> dict[str, dict[str, Any]]:
    records = {}
    for path in sorted(directory.glob("*.json")):
        record = _load_json(path)
        records[str(record.get("matchup_id") or path.stem)] = record
    return records


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
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _within_ee(value: float, truth: float) -> bool:
    return abs(value - truth) <= 0.05 + 0.15 * truth + 1.0e-12


def _clarify_cams_label(value: object) -> str:
    return str(value).replace("CAMS", "surface-anchor CAMS")


def _window_mask(
    x: np.ndarray,
    y: np.ndarray,
    *,
    lon: float,
    lat: float,
    crs: str,
    radius_m: float = 1500.0,
) -> tuple[np.ndarray, int, int]:
    center_x, center_y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon, lat)
    xx, yy = np.meshgrid(x, y)
    mask = np.square(xx - center_x) + np.square(yy - center_y) <= radius_m**2
    ix = int(np.argmin(np.abs(x - center_x)))
    iy = int(np.argmin(np.abs(y - center_y)))
    return mask, iy, ix


def _argmin_map(cube: np.ndarray, axis: np.ndarray) -> np.ndarray:
    finite = np.any(np.isfinite(cube), axis=0)
    indices = np.argmin(np.where(np.isfinite(cube), cube, np.inf), axis=0)
    values = np.asarray(axis, dtype=np.float64)[indices]
    values[~finite] = np.nan
    return values


def _median_curve(cube: np.ndarray, mask: np.ndarray) -> np.ndarray:
    selected = np.asarray(cube, dtype=np.float64)[:, mask]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(selected, axis=1)


def _normalise_curve(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    finite = array[np.isfinite(array)]
    floor = float(np.min(finite)) if finite.size else 0.0
    return np.log10(1.0 + np.maximum(array - floor, 0.0))


def _normalise_rgb(channels: Sequence[np.ndarray]) -> np.ndarray:
    rgb = np.stack([np.asarray(channel, dtype=np.float64) for channel in channels], axis=-1)
    finite = rgb[np.isfinite(rgb)]
    low = float(np.percentile(finite, 1)) if finite.size else 0.0
    high = float(np.percentile(finite, 99)) if finite.size else 0.3
    scale = max(high - low, 0.03)
    output = np.clip((rgb - low) / scale, 0.0, 1.0) ** (1.0 / 1.25)
    return np.nan_to_num(output, nan=0.0, posinf=1.0, neginf=0.0)


def _mark_station(ax: Any, iy: int, ix: int, radius_px: float) -> None:
    from matplotlib.patches import Circle

    ax.plot(ix, iy, "+", color="white", markersize=12, markeredgewidth=2.4)
    ax.plot(ix, iy, "+", color="black", markersize=8, markeredgewidth=1.1)
    ax.add_patch(
        Circle(
            (ix, iy),
            radius=radius_px,
            edgecolor="white",
            facecolor="none",
            linewidth=1.1,
            linestyle="--",
            alpha=0.9,
        )
    )


def _panel(
    fig: Any,
    ax: Any,
    values: np.ndarray,
    title: str,
    *,
    iy: int,
    ix: int,
    radius_px: float,
    cmap: str,
    vmin: float | None = None,
    vmax: float | None = None,
    colorbar: bool = True,
) -> None:
    image = ax.imshow(values, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    _mark_station(ax, iy, ix, radius_px)
    ax.set_title(title, fontsize=9.5)
    ax.set_xticks([])
    ax.set_yticks([])
    if colorbar:
        fig.colorbar(image, ax=ax, fraction=0.042, pad=0.02)


def _pooled_maps(
    cost: np.lib.npyio.NpzFile,
    *,
    prior_unc: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    uncertainty = cost["aot_prior_unc"] if prior_unc is None else prior_unc
    return tuple(
        np.asarray(value, dtype=np.float64)
        for value in surface_driven_pool_argmin(
            np.ascontiguousarray(cost["cube"], dtype=np.float64),
            np.ascontiguousarray(cost["aot_axis"], dtype=np.float64),
            np.ascontiguousarray(cost["aot_prior"], dtype=np.float64),
            np.ascontiguousarray(uncertainty, dtype=np.float64),
            np.ascontiguousarray(cost["solve_valid"], dtype=bool),
            int(np.asarray(cost["pool_window"]).item()),
            int(np.asarray(cost["min_count"]).item()),
        )
    )


def _save_spatial_figure(
    path: Path,
    *,
    cost: np.lib.npyio.NpzFile,
    record: dict[str, Any],
    pooled_aod: np.ndarray,
    pooled_unc: np.ndarray,
    no_backstop_aod: np.ndarray,
    iy: int,
    ix: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    axis = np.asarray(cost["aot_axis"], dtype=np.float64)
    cube = np.asarray(cost["cube"], dtype=np.float64)
    band_cost = np.asarray(cost["band_cost_cube"], dtype=np.float64)
    band_residual = np.asarray(cost["band_residual_cube"], dtype=np.float64)
    truth = float(record["truth"])
    truth_index = int(np.argmin(np.abs(axis - truth)))
    total_argmin = _argmin_map(cube, axis)
    band_argmins = [_argmin_map(band_cost[index], axis) for index in range(3)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        local_min = np.nanmin(cube, axis=0)
    truth_penalty = cube[truth_index] - local_min
    truth_penalty = np.log10(1.0 + np.maximum(truth_penalty, 0.0))
    valid = np.asarray(cost["solve_valid"], dtype=bool)
    toa = np.asarray(cost["toa"], dtype=np.float64)
    prior = np.asarray(cost["boa_prior"], dtype=np.float64)
    surface_unc = np.asarray(cost["boa_unc"], dtype=np.float64)
    atmo_prior = np.asarray(cost["aot_prior"], dtype=np.float64)
    atmo_prior_unc = np.asarray(cost["aot_prior_unc"], dtype=np.float64)
    backstop_delta = pooled_aod - no_backstop_aod
    finite_delta = np.abs(backstop_delta[np.isfinite(backstop_delta)])
    delta_vmax = max(float(np.percentile(finite_delta, 98)), 0.05) if finite_delta.size else 0.05
    x = np.asarray(cost["x"], dtype=np.float64)
    pixel_size = float(np.median(np.abs(np.diff(x)))) if x.size > 1 else 60.0
    radius_px = 1500.0 / max(pixel_size, 1.0)
    aod_max = min(4.0, max(1.0, math.ceil(truth * 1.4 * 10.0) / 10.0))

    fig, axes = plt.subplots(3, 6, figsize=(22.2, 11.5), facecolor="white")
    axes = axes.ravel()
    axes[0].imshow(_normalise_rgb([toa[2], toa[1], toa[0]]), interpolation="nearest")
    _mark_station(axes[0], iy, ix, radius_px)
    axes[0].set_title("Target TOA true colour\nB04 / B03 / B02", fontsize=9.5)
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[1].imshow(_normalise_rgb([prior[2], prior[1], prior[0]]), interpolation="nearest")
    _mark_station(axes[1], iy, ix, radius_px)
    axes[1].set_title("Predicted surface prior\nB04 / B03 / B02", fontsize=9.5)
    axes[1].set_xticks([])
    axes[1].set_yticks([])
    _panel(
        fig,
        axes[2],
        surface_unc[1],
        "Predicted B03 surface uncertainty",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="cividis",
        vmin=0.0,
        vmax=float(np.nanpercentile(surface_unc[1], 98)),
    )
    _panel(
        fig,
        axes[3],
        valid.astype(float),
        f"Solver support mask\n{100 * valid.mean():.1f}% valid",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="Greys",
        vmin=0.0,
        vmax=1.0,
    )
    _panel(
        fig,
        axes[4],
        atmo_prior,
        "Solver MAIAC AOD backstop",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="viridis",
        vmin=0.0,
        vmax=aod_max,
    )
    _panel(
        fig,
        axes[5],
        atmo_prior_unc,
        "Solver MAIAC backstop sigma",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="cividis",
        vmin=0.0,
        vmax=float(np.nanpercentile(atmo_prior_unc, 98)),
    )
    _panel(
        fig,
        axes[6],
        pooled_aod,
        "Operational pooled AOD",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="viridis",
        vmin=0.0,
        vmax=aod_max,
    )
    _panel(
        fig,
        axes[7],
        no_backstop_aod,
        "No-backstop pooled AOD replay",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="viridis",
        vmin=0.0,
        vmax=aod_max,
    )
    _panel(
        fig,
        axes[8],
        backstop_delta,
        "Backstop effect\noperational - no-backstop AOD",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="coolwarm",
        vmin=-delta_vmax,
        vmax=delta_vmax,
    )
    _panel(
        fig,
        axes[9],
        pooled_unc,
        "Operational pooled AOD uncertainty",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="cividis",
        vmin=0.0,
        vmax=float(np.nanpercentile(pooled_unc, 98)),
    )
    _panel(
        fig,
        axes[10],
        total_argmin,
        "Raw total-cost argmin AOD",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="viridis",
        vmin=0.0,
        vmax=aod_max,
    )
    _panel(
        fig,
        axes[11],
        truth_penalty,
        "Truth-node cost penalty\nlog10(1 + cost - local min)",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="magma",
        vmin=0.0,
        vmax=float(np.nanpercentile(truth_penalty, 98)),
    )
    for index, band in enumerate(BANDS):
        _panel(
            fig,
            axes[12 + index],
            band_argmins[index],
            f"{band} cost argmin AOD",
            iy=iy,
            ix=ix,
            radius_px=radius_px,
            cmap="viridis",
            vmin=0.0,
            vmax=aod_max,
        )
    residual_vmax = float(np.nanpercentile(band_residual[:, truth_index], 98))
    for index, band in enumerate(BANDS):
        _panel(
            fig,
            axes[15 + index],
            band_residual[index, truth_index],
            f"{band} |BOA - prior| at truth",
            iy=iy,
            ix=ix,
            radius_px=radius_px,
            cmap="magma",
            vmin=0.0,
            vmax=residual_vmax,
        )
    fig.suptitle(
        f"{record['matchup_id']}  |  spatial evidence on the 60 m solver grid",
        fontsize=13,
        x=0.02,
        ha="left",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=125, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _save_diagnostic_figure(
    path: Path,
    *,
    cost: np.lib.npyio.NpzFile,
    record: dict[str, Any],
    prior_record: dict[str, Any],
    local_mask: np.ndarray,
    pooled_aod: np.ndarray,
    pooled_unc: np.ndarray,
    candidates: list[dict[str, Any]],
) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    axis = np.asarray(cost["aot_axis"], dtype=np.float64)
    cube = np.asarray(cost["cube"], dtype=np.float64)
    band_cost = np.asarray(cost["band_cost_cube"], dtype=np.float64)
    residual = np.asarray(cost["band_residual_cube"], dtype=np.float64)
    valid = np.asarray(cost["solve_valid"], dtype=bool)
    window_valid = local_mask & valid
    curve_valid = window_valid
    curve_scope = "1.5 km station window"
    if int(window_valid.sum()) < 20:
        curve_valid = valid
        curve_scope = "scene support; station window has fewer than 20 valid pixels"
    if not np.any(curve_valid):
        raise ValueError(f"No valid cost-cube support for {record['matchup_id']}")
    truth = float(record["truth"])
    retrieved = float(record["retrieved"])
    surface_anchor_cams = float(prior_record["cams_aot"])
    solver_cams = float(np.nanmedian(np.asarray(cost["aot_prior"])[curve_valid]))
    threshold = 0.05 + 0.15 * truth
    truth_index = int(np.argmin(np.abs(axis - truth)))
    retrieved_index = int(np.argmin(np.abs(axis - retrieved)))
    local_total = _median_curve(cube, curve_valid)
    scene_total = _median_curve(cube, valid)
    local_band_cost = [_median_curve(band_cost[index], curve_valid) for index in range(3)]
    local_residual = [_median_curve(residual[index], curve_valid) for index in range(3)]
    toa = np.asarray(cost["toa"], dtype=np.float64)
    boa_prior = np.asarray(cost["boa_prior"], dtype=np.float64)
    boa_unc = np.asarray(cost["boa_unc"], dtype=np.float64)

    fig, axes = plt.subplots(2, 3, figsize=(15.8, 9.2), facecolor="white")
    axes = axes.ravel()

    def aod_context(ax: Any) -> None:
        ax.axvspan(truth - threshold, truth + threshold, color="#eceff1", alpha=0.9, label="EE")
        ax.axvline(truth, color="#111111", linewidth=1.6, label="AERONET")
        ax.axvline(retrieved, color="#0b67a3", linewidth=1.4, linestyle="--", label="Retrieved")
        ax.axvline(
            solver_cams,
            color="#737373",
            linewidth=1.2,
            linestyle=":",
            label="Solver MAIAC backstop",
        )
        ax.axvline(
            surface_anchor_cams,
            color="#ad6800",
            linewidth=1.0,
            linestyle="-.",
            label="Surface-anchor CAMS",
        )

    axes[0].plot(axis, _normalise_curve(scene_total), color="#111111", label="Scene median")
    axes[0].plot(axis, _normalise_curve(local_total), color="#0b67a3", label=curve_scope)
    aod_context(axes[0])
    axes[0].set_title("Total chi-square cost profile")
    axes[0].set_xlabel("AOD")
    axes[0].set_ylabel("log10(1 + cost above minimum)")
    axes[0].set_xlim(axis[0], min(axis[-1], max(1.5, truth * 1.45)))
    axes[0].legend(fontsize=7.5, ncol=2)
    axes[0].grid(alpha=0.18)

    for index, band in enumerate(BANDS):
        axes[1].plot(
            axis,
            _normalise_curve(local_band_cost[index]),
            color=BAND_COLORS[band],
            label=band,
        )
    aod_context(axes[1])
    axes[1].set_title("Per-band chi-square cost profiles")
    axes[1].set_xlabel("AOD")
    axes[1].set_ylabel("log10(1 + cost above minimum)")
    axes[1].set_xlim(axis[0], min(axis[-1], max(1.5, truth * 1.45)))
    axes[1].legend(fontsize=7.5, ncol=3)
    axes[1].grid(alpha=0.18)

    for index, band in enumerate(BANDS):
        axes[2].plot(axis, local_residual[index], color=BAND_COLORS[band], label=band)
    aod_context(axes[2])
    axes[2].set_title("Per-band absolute BOA-prior residual")
    axes[2].set_xlabel("AOD")
    axes[2].set_ylabel("Absolute reflectance residual")
    axes[2].set_xlim(axis[0], min(axis[-1], max(1.5, truth * 1.45)))
    axes[2].legend(fontsize=7.5, ncol=3)
    axes[2].grid(alpha=0.18)

    wavelengths = np.asarray([WAVELENGTHS[band] for band in BANDS])
    toa_median = np.asarray([np.nanmedian(toa[index][curve_valid]) for index in range(3)])
    prior_median = np.asarray([np.nanmedian(boa_prior[index][curve_valid]) for index in range(3)])
    unc_median = np.asarray([np.nanmedian(boa_unc[index][curve_valid]) for index in range(3)])
    axes[3].plot(wavelengths, toa_median, "o-", color="#303030", label="TOA")
    axes[3].errorbar(
        wavelengths,
        prior_median,
        yerr=unc_median,
        fmt="o-",
        color="#0b67a3",
        capsize=4,
        label="Predicted BOA prior",
    )
    axes[3].set_xticks(wavelengths, BANDS)
    axes[3].set_title(f"Median reflectance spectrum\n{curve_scope}")
    axes[3].set_ylabel("Reflectance")
    axes[3].legend(fontsize=8)
    axes[3].grid(alpha=0.18)

    x_positions = np.arange(3)
    width = 0.34
    current_residual = np.asarray([local_residual[index][retrieved_index] for index in range(3)])
    truth_residual = np.asarray([local_residual[index][truth_index] for index in range(3)])
    axes[4].bar(
        x_positions - width / 2, current_residual, width, color="#0b67a3", label="Retrieved node"
    )
    axes[4].bar(x_positions + width / 2, truth_residual, width, color="#737373", label="Truth node")
    axes[4].plot(x_positions, unc_median, "D", color="#111111", label="Prior uncertainty")
    axes[4].set_xticks(x_positions, BANDS)
    axes[4].set_title("Spectral closure at retrieved and truth AOD")
    axes[4].set_ylabel("Absolute reflectance residual")
    axes[4].legend(fontsize=8)
    axes[4].grid(axis="y", alpha=0.18)

    plotted = [item for item in candidates if item.get("value") is not None]
    plotted = plotted[:24]
    y_positions = np.arange(len(plotted))
    values = np.asarray([float(item["value"]) for item in plotted])
    colors = ["#0b67a3" if item.get("within_ee") else "#777777" for item in plotted]
    axes[5].axvspan(truth - threshold, truth + threshold, color="#eceff1", alpha=0.9)
    axes[5].axvline(truth, color="#111111", linewidth=1.5)
    axes[5].scatter(values, y_positions, c=colors, s=30, zorder=3)
    for ypos, value in zip(y_positions, values):
        axes[5].plot([truth, value], [ypos, ypos], color="#cccccc", linewidth=0.8, zorder=1)
    axes[5].set_yticks(y_positions, [item["label"] for item in plotted], fontsize=7)
    axes[5].invert_yaxis()
    axes[5].set_title("Saved scalar outputs and diagnostic anchors")
    axes[5].set_xlabel("AOD")
    axes[5].grid(axis="x", alpha=0.18)

    for ax in axes:
        ax.tick_params(axis="both", labelsize=8)
        ax.title.set_fontsize(10)
        ax.xaxis.label.set_fontsize(9)
        ax.yaxis.label.set_fontsize(9)

    fig.suptitle(
        f"{record['matchup_id']}  |  spectral, cost, and scalar-output evidence",
        fontsize=13,
        x=0.02,
        ha="left",
    )
    fig.tight_layout(rect=(0.015, 0.01, 0.995, 0.955), pad=1.2)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=125, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    pooled_local = pooled_aod[window_valid]
    pooled_unc_local = pooled_unc[window_valid]

    def finite_median(values: np.ndarray) -> float | None:
        finite = np.asarray(values, dtype=np.float64)
        return float(np.nanmedian(finite)) if np.isfinite(finite).any() else None

    def finite_mean(values: np.ndarray) -> float | None:
        finite = np.asarray(values, dtype=np.float64)
        return float(np.nanmean(finite)) if np.isfinite(finite).any() else None

    return {
        "aod_axis": axis,
        "total_scene": _normalise_curve(scene_total),
        "total_local": _normalise_curve(local_total),
        "band_cost_local": {
            band: _normalise_curve(local_band_cost[index]) for index, band in enumerate(BANDS)
        },
        "band_residual_local": {band: local_residual[index] for index, band in enumerate(BANDS)},
        "curve_support_scope": curve_scope,
        "curve_support_count": int(curve_valid.sum()),
        "station_window_valid_count": int(window_valid.sum()),
        "pooled_local_median": finite_median(pooled_local),
        "pooled_local_mean": finite_mean(pooled_local),
        "pooled_local_unc_median": finite_median(pooled_unc_local),
        "toa_local_median": dict(zip(BANDS, toa_median)),
        "prior_local_median": dict(zip(BANDS, prior_median)),
        "prior_unc_local_median": dict(zip(BANDS, unc_median)),
        "residual_retrieved": dict(zip(BANDS, current_residual)),
        "residual_truth": dict(zip(BANDS, truth_residual)),
    }


def _candidate_outputs(
    matchup_id: str,
    *,
    failure_row: dict[str, Any],
    current: dict[str, Any],
    baseline: dict[str, Any],
    prior: dict[str, Any],
    variants: dict[str, dict[str, dict[str, Any]]],
    solver_cams: float | None = None,
    no_backstop: float | None = None,
) -> list[dict[str, Any]]:
    truth = float(current["truth"])
    extraction = current.get("retrieval_extraction") or {}
    solver = current.get("solver") or {}
    values: list[tuple[str, str, float | None]] = [
        ("current", "Current scene mean", current.get("retrieved")),
        ("extraction", "Station pixel", extraction.get("station")),
        ("extraction", "1.5 km median", extraction.get("winmed")),
        ("reference", "Historical R2", baseline.get("retrieved")),
        ("reference", "Surface-prior anchor CAMS", prior.get("cams_aot")),
        ("backstop", "Solver CAMS backstop", solver_cams),
        ("sensitivity", "No-backstop replay, scene mean", no_backstop),
        ("cost", "Global total-cost min", solver.get("surface_cost_curve_min_aot")),
        ("cost", "B02 cost min", solver.get("surface_band_B02_argmin_aot")),
        ("cost", "B03 cost min", solver.get("surface_band_B03_argmin_aot")),
        ("cost", "B04 cost min", solver.get("surface_band_B04_argmin_aot")),
        ("audit", "ExtraTrees OOF audit", failure_row.get("oof_aod")),
    ]
    bands = np.asarray(
        [solver.get(f"surface_band_{band}_argmin_aot") for band in BANDS],
        dtype=np.float64,
    )
    cams = float(prior["cams_aot"])
    values.append(
        (
            "consensus",
            "Fixed median posterior 0.275",
            0.275 * float(np.median(bands)) + 0.725 * cams,
        )
    )
    for option, records in variants.items():
        record = records.get(matchup_id)
        if record is None or str(record.get("status", "")).upper() != "OK":
            continue
        values.append(("variant", SINGLE_PRIOR_LABELS.get(option, option), record.get("retrieved")))
    output = []
    for kind, label, raw_value in values:
        try:
            value = float(raw_value) if raw_value is not None else None
        except (TypeError, ValueError):
            value = None
        if value is not None and not math.isfinite(value):
            value = None
        output.append(
            {
                "kind": kind,
                "label": label,
                "value": value,
                "within_ee": _within_ee(value, truth) if value is not None else None,
            }
        )
    return output


def _save_global_figures(
    output: Path,
    *,
    all_current: dict[str, dict[str, Any]],
    failures: list[dict[str, Any]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    all_rows = [
        record
        for record in all_current.values()
        if str(record.get("status", "")).upper() == "OK"
        and record.get("truth") is not None
        and record.get("retrieved") is not None
    ]
    truth = np.asarray([row["truth"] for row in all_rows], dtype=np.float64)
    retrieved = np.asarray([row["retrieved"] for row in all_rows], dtype=np.float64)
    ee = 0.05 + 0.15 * truth
    hit = np.abs(retrieved - truth) <= ee
    max_aod = max(1.8, float(max(np.max(truth), np.max(retrieved))) * 1.04)

    fig, axes = plt.subplots(2, 2, figsize=(14.8, 10.5), facecolor="white")
    axes = axes.ravel()
    grid = np.linspace(0.0, max_aod, 300)
    axes[0].fill_between(
        grid,
        grid - (0.05 + 0.15 * grid),
        grid + (0.05 + 0.15 * grid),
        color="#e8edf1",
        label="Expected error",
    )
    axes[0].scatter(
        truth[hit], retrieved[hit], s=18, color="#9aa5ad", alpha=0.65, label="Within EE"
    )
    for code, color in GROUP_COLORS.items():
        rows = [row for row in failures if row["mechanism_code"] == code]
        if rows:
            axes[0].scatter(
                [row["truth_aod"] for row in rows],
                [row["retrieved_aod"] for row in rows],
                s=30,
                color=color,
                label=f"Group {code}",
            )
    axes[0].plot(grid, grid, color="#111111", linewidth=1.0)
    axes[0].set(xlabel="AERONET AOD", ylabel="Retrieved AOD", xlim=(0, max_aod), ylim=(0, max_aod))
    axes[0].set_title("All 152 current retrievals")
    axes[0].legend(fontsize=7, ncol=2)
    axes[0].grid(alpha=0.16)

    axes[1].scatter(
        [row["truth_aod"] for row in failures],
        [row["error"] / row["ee_threshold"] for row in failures],
        c=[GROUP_COLORS[row["mechanism_code"]] for row in failures],
        s=34,
    )
    axes[1].axhline(-1.0, color="#777777", linestyle="--")
    axes[1].axhline(1.0, color="#777777", linestyle="--")
    axes[1].set(xlabel="AERONET AOD", ylabel="Signed error / EE")
    axes[1].set_title("Failure direction and normalized magnitude")
    axes[1].grid(alpha=0.16)

    axes[2].scatter(
        [row["band_spread"] for row in failures],
        [row["error"] / row["ee_threshold"] for row in failures],
        c=[GROUP_COLORS[row["mechanism_code"]] for row in failures],
        s=34,
    )
    axes[2].axvline(0.2, color="#777777", linestyle="--")
    axes[2].axhline(0.0, color="#111111", linewidth=0.8)
    axes[2].set(xlabel="B02/B03/B04 argmin spread", ylabel="Signed error / EE")
    axes[2].set_title("Band preference spread")
    axes[2].grid(alpha=0.16)

    cams_offset = [
        (row["solver_cams_aod"] - row["truth_aod"]) / row["ee_threshold"] for row in failures
    ]
    axes[3].scatter(
        cams_offset,
        [row["error"] / row["ee_threshold"] for row in failures],
        c=[GROUP_COLORS[row["mechanism_code"]] for row in failures],
        s=34,
    )
    axes[3].axvline(0.0, color="#777777", linestyle="--")
    axes[3].axhline(0.0, color="#111111", linewidth=0.8)
    axes[3].set(xlabel="(Solver CAMS - AERONET) / EE", ylabel="Signed retrieval error / EE")
    axes[3].set_title("Operational atmospheric-backstop offset")
    axes[3].grid(alpha=0.16)

    fig.suptitle(
        "Cross-case evidence; colours are rule-based filter groups, not diagnoses", fontsize=13
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(output / "assets" / "cross_case_overview.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    signed_fields = (
        ("Current", lambda row: row["error"] / row["ee_threshold"]),
        (
            "No backstop",
            lambda row: (row["no_backstop_aod"] - row["truth_aod"]) / row["ee_threshold"],
        ),
        (
            "Solver CAMS",
            lambda row: (row["solver_cams_aod"] - row["truth_aod"]) / row["ee_threshold"],
        ),
        (
            "Surface-anchor CAMS",
            lambda row: (row["surface_anchor_cams_aod"] - row["truth_aod"]) / row["ee_threshold"],
        ),
        ("Curve min", lambda row: (row["curve_min_aod"] - row["truth_aod"]) / row["ee_threshold"]),
        ("B02 min", lambda row: (row["band_b02_min_aod"] - row["truth_aod"]) / row["ee_threshold"]),
        ("B03 min", lambda row: (row["band_b03_min_aod"] - row["truth_aod"]) / row["ee_threshold"]),
        ("B04 min", lambda row: (row["band_b04_min_aod"] - row["truth_aod"]) / row["ee_threshold"]),
    )
    context_fields = (
        ("Cloud", lambda row: row["cloud_fraction"]),
        ("Valid support", lambda row: row["valid_support_fraction"]),
        ("Cost / band", lambda row: row["cost_per_band"]),
        ("Band spread", lambda row: row["band_spread"]),
        ("|Backstop shift|", lambda row: abs(row["backstop_shift"])),
    )
    ordered = sorted(failures, key=lambda row: row["severity_rank"])
    signed = np.asarray([[function(row) for _, function in signed_fields] for row in ordered])
    context = np.asarray([[function(row) for _, function in context_fields] for row in ordered])
    context_rank = np.zeros_like(context)
    for column in range(context.shape[1]):
        order = np.argsort(np.argsort(context[:, column]))
        context_rank[:, column] = order / max(context.shape[0] - 1, 1)

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(14.5, 14.5),
        gridspec_kw={"width_ratios": [1.45, 1.0]},
        facecolor="white",
    )
    image1 = ax1.imshow(signed, aspect="auto", cmap="RdBu_r", vmin=-3.0, vmax=3.0)
    ax1.set_xticks(
        range(len(signed_fields)), [label for label, _ in signed_fields], rotation=35, ha="right"
    )
    row_labels = [f"{row['site']} | {row['matchup_id'].rsplit('_', 1)[-1][:8]}" for row in ordered]
    ax1.set_yticks(range(len(ordered)), row_labels, fontsize=7)
    ax1.set_title("AOD evidence offset from AERONET, in EE units")
    fig.colorbar(image1, ax=ax1, fraction=0.035, pad=0.02, label="Signed offset / EE")
    image2 = ax2.imshow(context_rank, aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    ax2.set_xticks(
        range(len(context_fields)), [label for label, _ in context_fields], rotation=35, ha="right"
    )
    ax2.set_yticks([])
    ax2.set_title("Within-column percentile for context fields")
    fig.colorbar(image2, ax=ax2, fraction=0.05, pad=0.02, label="Column percentile")
    fig.suptitle("Failure evidence matrix; rows are ordered by absolute error / EE", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(output / "assets" / "evidence_matrix.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def build(root: Path, output: Path) -> dict[str, Any]:
    analysis = _load_json(root / "reports/aod-low-cloud-20260711/low-cloud-failure-analysis.json")
    failures = list(analysis["failure_cases"])
    current = _load_records(root / CURRENT_RESULTS.name)
    diagnostic = _load_records(root / DIAGNOSTIC_RESULTS.name)
    baseline = _load_records(root / BASELINE_RESULTS.name)
    priors = _load_records(root / PRIOR_RESULTS.name)
    variants = {
        option: _load_records(root / f"phaseD_results_lowcloud20_singleprior_{option}_20260711")
        for option in OPTIONS
        if option != "b03_chi2"
    }
    expected = {row["matchup_id"] for row in failures}
    missing = {
        "diagnostic_results": sorted(expected - set(diagnostic)),
        "cost_cubes": sorted(
            matchup_id
            for matchup_id in expected
            if not (root / COST_CUBES.name / f"{matchup_id}.npz").exists()
        ),
    }
    if any(missing.values()):
        raise ValueError(f"Diagnostic artifacts are incomplete: {missing}")

    output.mkdir(parents=True, exist_ok=True)
    (output / "assets" / "spatial").mkdir(parents=True, exist_ok=True)
    (output / "assets" / "diagnostic").mkdir(parents=True, exist_ok=True)
    (output / "data").mkdir(parents=True, exist_ok=True)
    shutil.copy2(WEB_ASSETS / "app.css", output / "app.css")
    shutil.copy2(WEB_ASSETS / "app.js", output / "app.js")

    case_rows = []
    for index, failure in enumerate(failures, start=1):
        matchup_id = failure["matchup_id"]
        record = current[matchup_id]
        prior_record = priors[matchup_id]
        cube_path = root / COST_CUBES.name / f"{matchup_id}.npz"
        with np.load(cube_path, allow_pickle=False) as cost:
            x = np.asarray(cost["x"], dtype=np.float64)
            y = np.asarray(cost["y"], dtype=np.float64)
            local_mask, iy, ix = _window_mask(
                x,
                y,
                lon=float(record["lon"]),
                lat=float(record["lat"]),
                crs=str(record["scene_crs"]),
            )
            solve_valid = np.asarray(cost["solve_valid"], dtype=bool)
            pooled_aod, pooled_unc, pooled_cost = _pooled_maps(cost)
            no_backstop_aod, _no_backstop_unc, _no_backstop_cost = _pooled_maps(
                cost,
                prior_unc=np.full(solve_valid.shape, np.inf, dtype=np.float64),
            )
            no_backstop_values = no_backstop_aod[np.isfinite(no_backstop_aod)]
            no_backstop_scene_mean = (
                float(np.mean(no_backstop_values)) if no_backstop_values.size else None
            )
            solver_cams_values = np.asarray(cost["aot_prior"], dtype=np.float64)[solve_valid]
            solver_cams_unc_values = np.asarray(cost["aot_prior_unc"], dtype=np.float64)[
                solve_valid
            ]
            solver_cams_aod = float(np.nanmedian(solver_cams_values))
            solver_cams_unc = float(np.nanmedian(solver_cams_unc_values))
            candidates = _candidate_outputs(
                matchup_id,
                failure_row=failure,
                current=record,
                baseline=baseline[matchup_id],
                prior=prior_record,
                variants=variants,
                solver_cams=solver_cams_aod,
                no_backstop=no_backstop_scene_mean,
            )
            spatial_path = output / "assets" / "spatial" / f"{matchup_id}.png"
            diagnostic_path = output / "assets" / "diagnostic" / f"{matchup_id}.png"
            _save_spatial_figure(
                spatial_path,
                cost=cost,
                record=record,
                pooled_aod=pooled_aod,
                pooled_unc=pooled_unc,
                no_backstop_aod=no_backstop_aod,
                iy=iy,
                ix=ix,
            )
            curve_data = _save_diagnostic_figure(
                diagnostic_path,
                cost=cost,
                record=record,
                prior_record=prior_record,
                local_mask=local_mask,
                pooled_aod=pooled_aod,
                pooled_unc=pooled_unc,
                candidates=candidates,
            )
            local_valid = local_mask & solve_valid
            local_cost = pooled_cost[local_valid]
            cube_meta = {
                "shape": list(np.asarray(cost["cube"]).shape),
                "bands": [str(value) for value in np.asarray(cost["band_names"])],
                "aod_nodes": int(np.asarray(cost["aot_axis"]).size),
                "aod_min": float(np.min(cost["aot_axis"])),
                "aod_max": float(np.max(cost["aot_axis"])),
                "pool_window": int(np.asarray(cost["pool_window"]).item()),
                "pool_min_count": int(np.asarray(cost["min_count"]).item()),
                "local_pooled_cost_median": (
                    float(np.nanmedian(local_cost)) if np.isfinite(local_cost).any() else None
                ),
            }

        case_rows.append(
            {
                **failure,
                "mechanism": _clarify_cams_label(failure["mechanism"]),
                "diagnostic_evidence": _clarify_cams_label(failure.get("diagnostic_evidence", "")),
                "case_index": index,
                "spatial_image": f"assets/spatial/{matchup_id}.png",
                "diagnostic_image": f"assets/diagnostic/{matchup_id}.png",
                "canonical_result_url": (f"../../{CURRENT_RESULTS.name}/{matchup_id}.json"),
                "diagnostic_result_url": (f"../../{DIAGNOSTIC_RESULTS.name}/{matchup_id}.json"),
                "cost_cube_url": f"../../{COST_CUBES.name}/{matchup_id}.npz",
                "candidates": candidates,
                "curves": curve_data,
                "cube": cube_meta,
                "canonical_solver": record.get("solver") or {},
                "retrieval_extraction": record.get("retrieval_extraction") or {},
                "prior_boa": record.get("prior_boa") or {},
                "prior_unc": record.get("prior_unc") or {},
                "diagnostic_retrieved": diagnostic[matchup_id].get("retrieved"),
                "surface_anchor_cams_aod": float(prior_record["cams_aot"]),
                "solver_cams_aod": solver_cams_aod,
                "solver_cams_unc": solver_cams_unc,
                "no_backstop_aod": no_backstop_scene_mean,
                "backstop_shift": (
                    float(record["retrieved"]) - no_backstop_scene_mean
                    if no_backstop_scene_mean is not None
                    else None
                ),
            }
        )

    _save_global_figures(output, all_current=current, failures=case_rows)
    data = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "title": "Low-cloud AOD failure evidence explorer",
        "method": "B02/B03/B04 chi-square with one S2 monthly SWIR/NIR-anchored ExtraTree surface prior",
        "cohort_count": 152,
        "failure_count": len(case_rows),
        "stats": {
            "status_ok": len(
                [
                    record
                    for record in current.values()
                    if str(record.get("status", "")).upper() == "OK"
                ]
            ),
            "within_ee": len(current) - len(case_rows),
            "within_ee_rate": (len(current) - len(case_rows)) / len(current),
            "outside_ee": len(case_rows),
            "underreads": sum(row["direction"] == "Underread" for row in case_rows),
            "overreads": sum(row["direction"] == "Overread" for row in case_rows),
        },
        "cases": _jsonable(case_rows),
        "group_definitions": {
            row["mechanism_code"]: {
                "label": _clarify_cams_label(row["mechanism"]),
                "evidence": _clarify_cams_label(row["diagnostic_evidence"]),
                "count": row["count"],
                "share_of_misses": row["share_of_misses"],
                "underreads": row["underreads"],
                "overreads": row["overreads"],
            }
            for row in analysis["mechanisms"]
        },
        "sources": {
            "failure_csv": "data/failure-cases.csv",
            "failure_analysis": "../aod-low-cloud-20260711/low-cloud-failure-analysis.json",
            "current_results": f"../../{CURRENT_RESULTS.name}/",
            "diagnostic_results": f"../../{DIAGNOSTIC_RESULTS.name}/",
            "cost_cubes": f"../../{COST_CUBES.name}/",
        },
        "notes": [
            "Every case shown is status OK and outside expected error.",
            "Rule-based groups are optional filters, not diagnoses or decisions.",
            "Spatial maps and curves come from the exact current-method diagnostic cost cubes.",
            "Surface-prior anchor CAMS and solver CAMS backstop are recorded separately.",
            "The no-backstop replay keeps the same surface cost, spatial pooling, and scene-mean extraction.",
            "AERONET is used only as the evaluation reference shown in the explorer.",
        ],
    }
    (output / "data" / "cases.json").write_text(
        json.dumps(data, separators=(",", ":")), encoding="utf-8"
    )
    shutil.copy2(
        root / "reports/aod-low-cloud-20260711/low-cloud-failure-cases.csv",
        output / "data" / "failure-cases.csv",
    )
    index = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="light dark">
  <title>Low-cloud AOD failure evidence explorer</title>
  <link rel="icon" href="data:,">
  <link rel="stylesheet" href="app.css">
</head>
<body>
  <div id="app" class="app-shell" aria-live="polite">
    <div class="loading-state">Loading 42-case diagnostic snapshot...</div>
  </div>
  <noscript>This investigation workspace requires JavaScript. The source CSV and PNG assets remain available in this directory.</noscript>
  <script src="app.js"></script>
</body>
</html>
"""
    (output / "index.html").write_text(index, encoding="utf-8")
    receipt = {
        "output": str(output),
        "cases": len(case_rows),
        "spatial_images": len(list((output / "assets" / "spatial").glob("*.png"))),
        "diagnostic_images": len(list((output / "assets" / "diagnostic").glob("*.png"))),
        "missing": missing,
    }
    (output / "build-receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    receipt = build(args.root, args.output)
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
