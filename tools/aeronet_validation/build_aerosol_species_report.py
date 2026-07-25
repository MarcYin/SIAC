"""Build the detailed low-cloud aerosol-species experiment web report."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

from tools.aeronet_validation.analyze_aerosol_species_experiment import (
    CCI_COSTS,
    COHORT,
    METHODS,
    PRIMARY_METHODS,
    ROOT,
    _load_records,
    _month,
    _retrieved,
    analyze,
)
from tools.aeronet_validation.build_low_cloud_failure_explorer import (
    BAND_COLORS,
    BANDS,
    WAVELENGTHS,
    _jsonable,
    _mark_station,
    _median_curve,
    _normalise_rgb,
    _window_mask,
)

from siac._rust_compat import surface_driven_pool_argmin
from siac.algorithms.rt.aerosol_species import (
    candidate_fraction_sets,
    climatology_fraction_percentages,
)

DEFAULT_OUTPUT = ROOT / "reports/aod-aerosol-species-20260712"
FIXED_COSTS = (
    ROOT / "phaseD_cost_cubes_lowcloud20_aerosol_species_sixs_continental_full_t4_20260712"
)
EXACT_CCI_COSTS = (
    ROOT / "phaseD_cost_cubes_lowcloud20_aerosol_species_sixs_cci_exact_exact_smoke_t4_20260712"
)
DIRECT_COSTS = ROOT / "phaseD_cost_cubes_lowcloud20_aerosol_species_libradtran_continental_20260712"
FAILURE_ANALYSIS = ROOT / "reports/aod-low-cloud-20260711/low-cloud-failure-analysis.json"
WEB_ASSETS = Path(__file__).with_name("aerosol_species_report")
SMOKE_COHORT = Path(__file__).with_name("aerosol_species_smoke_mids.txt")

METHOD_LABELS = {
    "lut_continental": "libRadTran continental Zarr LUT",
    "sixs_continental": "Native 6S fixed continental",
    "sixs_cci3": "Native 6S CCI nearest-3 selection",
    "sixs_cci_exact_smoke": "Native 6S exact monthly CCI mixture (smoke)",
    "libradtran_continental_smoke": "Direct libRadTran fixed continental (smoke)",
}
SPECIES_REPLAY_DELTA_TOLERANCE = 1e-4
METHOD_SHORT = {
    "lut_continental": "LUT continental",
    "sixs_continental": "6S continental",
    "sixs_cci3": "6S CCI-3",
    "sixs_cci_exact_smoke": "6S exact CCI",
    "libradtran_continental_smoke": "Direct libRadTran",
}
METHOD_COLORS = {
    "lut_continental": "#59636a",
    "sixs_continental": "#1565a8",
    "sixs_cci3": "#2d7a4b",
    "sixs_cci_exact_smoke": "#7b4f9d",
    "libradtran_continental_smoke": "#a65d00",
}
COMPONENTS = ("dust", "sea_salt", "fine_strong", "fine_weak")
COMPONENT_LABELS = {
    "dust": "Dust",
    "sea_salt": "Sea salt",
    "fine_strong": "Fine, strong absorption",
    "fine_weak": "Fine, weak absorption",
}
COMPONENT_COLORS = {
    "dust": "#bd7b2c",
    "sea_salt": "#3b83bd",
    "fine_strong": "#9b3f52",
    "fine_weak": "#4b8c61",
}


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _within_ee(value: float | None, truth: float) -> bool:
    return value is not None and abs(value - truth) <= 0.05 + 0.15 * truth + 1e-12


def _wilson_interval(hits: int, count: int, z: float = 1.959963984540054) -> list[float] | None:
    if count <= 0:
        return None
    p = hits / count
    denominator = 1.0 + z * z / count
    centre = (p + z * z / (2.0 * count)) / denominator
    half = z * math.sqrt(p * (1.0 - p) / count + z * z / (4.0 * count * count)) / denominator
    return [max(0.0, centre - half), min(1.0, centre + half)]


def _metric_summary(
    records: dict[str, dict[str, Any]],
    cohort: Sequence[str],
    truths: dict[str, float],
) -> dict[str, Any]:
    values: list[float] = []
    reference: list[float] = []
    runtimes: list[float] = []
    hit_ids: list[str] = []
    for matchup_id in cohort:
        value = _retrieved(records.get(matchup_id))
        if value is None:
            continue
        truth = truths[matchup_id]
        values.append(value)
        reference.append(truth)
        if _within_ee(value, truth):
            hit_ids.append(matchup_id)
        runtime = _finite(records[matchup_id].get("runtime_s"))
        if runtime is not None:
            runtimes.append(runtime)
    retrieved = np.asarray(values, dtype=np.float64)
    truth_values = np.asarray(reference, dtype=np.float64)
    errors = retrieved - truth_values
    if retrieved.size >= 2 and np.std(retrieved) > 0.0 and np.std(truth_values) > 0.0:
        correlation = float(np.corrcoef(truth_values, retrieved)[0, 1])
    else:
        correlation = None
    denominator = float(np.sum(np.square(truth_values - np.mean(truth_values))))
    r2 = (
        float(1.0 - np.sum(np.square(errors)) / denominator)
        if retrieved.size >= 2 and denominator > 0.0
        else None
    )
    hits = len(hit_ids)
    target_hits = math.floor(0.87 * len(cohort)) + 1 if cohort else 0
    return {
        "expected": len(cohort),
        "valid": int(retrieved.size),
        "hits": hits,
        "target_hits_above_87": target_hits,
        "additional_hits_to_above_87": max(0, target_hits - hits),
        "strict_rate": hits / len(cohort) if cohort else 0.0,
        "strict_rate_ci95": _wilson_interval(hits, len(cohort)),
        "valid_rate": hits / retrieved.size if retrieved.size else 0.0,
        "rmse": float(np.sqrt(np.mean(np.square(errors)))) if errors.size else None,
        "mae": float(np.mean(np.abs(errors))) if errors.size else None,
        "bias": float(np.mean(errors)) if errors.size else None,
        "median_abs_error": float(np.median(np.abs(errors))) if errors.size else None,
        "correlation": correlation,
        "r2": r2,
        "median_runtime_s": float(np.median(runtimes)) if runtimes else None,
        "hit_matchup_ids": hit_ids,
    }


def _bootstrap_abs_error_delta(
    baseline: dict[str, dict[str, Any]],
    candidate: dict[str, dict[str, Any]],
    cohort: Sequence[str],
    truths: dict[str, float],
    *,
    draws: int = 20_000,
) -> dict[str, Any]:
    deltas = []
    for matchup_id in cohort:
        base = _retrieved(baseline.get(matchup_id))
        test = _retrieved(candidate.get(matchup_id))
        if base is None or test is None:
            continue
        truth = truths[matchup_id]
        deltas.append(abs(test - truth) - abs(base - truth))
    values = np.asarray(deltas, dtype=np.float64)
    if values.size == 0:
        return {"count": 0, "mean": None, "ci95": None, "p_mean_below_zero": None}
    rng = np.random.default_rng(20260712)
    chunk = 2_000
    means = []
    for start in range(0, draws, chunk):
        size = min(chunk, draws - start)
        indices = rng.integers(0, values.size, size=(size, values.size))
        means.append(np.mean(values[indices], axis=1))
    bootstrap = np.concatenate(means)
    return {
        "count": int(values.size),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "ci95": [float(value) for value in np.quantile(bootstrap, [0.025, 0.975])],
        "p_mean_below_zero": float(np.mean(bootstrap < 0.0)),
    }


def _pooled(cost: np.lib.npyio.NpzFile) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    result = surface_driven_pool_argmin(
        np.ascontiguousarray(cost["cube"], dtype=np.float64),
        np.ascontiguousarray(cost["aot_axis"], dtype=np.float64),
        np.ascontiguousarray(cost["aot_prior"], dtype=np.float64),
        np.ascontiguousarray(cost["aot_prior_unc"], dtype=np.float64),
        np.ascontiguousarray(cost["solve_valid"], dtype=bool),
        int(np.asarray(cost["pool_window"]).item()),
        int(np.asarray(cost["min_count"]).item()),
    )
    return tuple(np.asarray(value, dtype=np.float64) for value in result)  # type: ignore[return-value]


def _normalise_curve(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    finite = array[np.isfinite(array)]
    floor = float(np.min(finite)) if finite.size else 0.0
    return np.log10(1.0 + np.maximum(array - floor, 0.0))


def _curve(cube: np.ndarray, mask: np.ndarray) -> np.ndarray:
    if not np.any(mask):
        return np.full(cube.shape[0], np.nan, dtype=np.float64)
    return _median_curve(np.asarray(cube, dtype=np.float64), mask)


def _component_text(mixture: dict[str, float]) -> str:
    return ", ".join(
        f"{COMPONENT_LABELS[name]} {100.0 * float(mixture[name]):.0f}%" for name in COMPONENTS
    )


def _candidate_label(index: int, mixture: dict[str, float]) -> str:
    pieces = sorted(mixture.items(), key=lambda item: item[1], reverse=True)
    return f"CCI-{index + 1}: " + " / ".join(
        f"{COMPONENT_LABELS[name]} {100.0 * value:.0f}%" for name, value in pieces if value > 0.0
    )


def _compact_mixture_label(mixture: dict[str, float]) -> str:
    return " ".join(
        [
            f"D{100.0 * mixture['dust']:.0f}",
            f"SS{100.0 * mixture['sea_salt']:.0f}",
            f"FS{100.0 * mixture['fine_strong']:.0f}",
            f"FW{100.0 * mixture['fine_weak']:.0f}",
        ]
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
    ticks: Sequence[float] | None = None,
    ticklabels: Sequence[str] | None = None,
) -> None:
    image = ax.imshow(values, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    _mark_station(ax, iy, ix, radius_px)
    ax.set_title(title, fontsize=9.3)
    ax.set_xticks([])
    ax.set_yticks([])
    if colorbar:
        bar = fig.colorbar(image, ax=ax, fraction=0.042, pad=0.02, ticks=ticks)
        if ticklabels is not None:
            bar.ax.set_yticklabels(ticklabels, fontsize=7)


def _candidate_residual_at_solution(
    residual_cube: np.ndarray,
    aot_axis: np.ndarray,
    aod_map: np.ndarray,
) -> np.ndarray:
    indices = np.abs(aot_axis[:, None, None] - aod_map[None, :, :]).argmin(axis=0)
    rows, columns = np.indices(aod_map.shape)
    output = residual_cube[:, indices, rows, columns]
    output[:, ~np.isfinite(aod_map)] = np.nan
    return np.asarray(output, dtype=np.float64)


def _save_spatial_figure(
    path: Path,
    *,
    matchup_id: str,
    truth: float,
    fixed: dict[str, Any],
    candidates: list[dict[str, Any]],
    selected_index: np.ndarray,
    selected_aod: np.ndarray,
    selected_cost: np.ndarray,
    confidence: np.ndarray,
    iy: int,
    ix: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    x = fixed["x"]
    pixel_size = float(np.median(np.abs(np.diff(x)))) if x.size > 1 else 60.0
    radius_px = 1500.0 / max(pixel_size, 1.0)
    all_aod = np.concatenate(
        [fixed["aod"][np.isfinite(fixed["aod"])]]
        + [item["aod"][np.isfinite(item["aod"])] for item in candidates]
    )
    upper = float(np.nanpercentile(all_aod, 99)) if all_aod.size else 1.0
    aod_max = min(4.0, max(0.5, truth * 1.35, upper))
    delta = selected_aod - fixed["aod"]
    finite_delta = np.abs(delta[np.isfinite(delta)])
    delta_max = max(0.05, float(np.nanpercentile(finite_delta, 98))) if finite_delta.size else 0.05
    cost_advantage = selected_cost - fixed["cost"]
    finite_advantage = np.abs(cost_advantage[np.isfinite(cost_advantage)])
    advantage_max = (
        max(0.1, float(np.nanpercentile(finite_advantage, 98))) if finite_advantage.size else 0.1
    )
    confidence_max = (
        max(0.1, float(np.nanpercentile(confidence[np.isfinite(confidence)], 98)))
        if np.any(np.isfinite(confidence))
        else 0.1
    )

    fig, axes = plt.subplots(3, 5, figsize=(19.5, 11.6), facecolor="white")
    axes = axes.ravel()
    axes[0].imshow(_normalise_rgb([fixed["toa"][2], fixed["toa"][1], fixed["toa"][0]]))
    _mark_station(axes[0], iy, ix, radius_px)
    axes[0].set_title("Target TOA true colour\nB04 / B03 / B02", fontsize=9.3)
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[1].imshow(
        _normalise_rgb([fixed["boa_prior"][2], fixed["boa_prior"][1], fixed["boa_prior"][0]])
    )
    _mark_station(axes[1], iy, ix, radius_px)
    axes[1].set_title("Fixed-6S final surface prior\nB04 / B03 / B02", fontsize=9.3)
    axes[1].set_xticks([])
    axes[1].set_yticks([])
    _panel(
        fig,
        axes[2],
        fixed["valid"].astype(float),
        f"Solver support mask\n{100.0 * fixed['valid'].mean():.1f}% valid",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="Greys",
        vmin=0.0,
        vmax=1.0,
    )
    _panel(
        fig,
        axes[3],
        fixed["aot_prior"],
        "Solver CAMS/MAIAC backstop AOD",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="viridis",
        vmin=0.0,
        vmax=aod_max,
    )
    _panel(
        fig,
        axes[4],
        fixed["aod"],
        "Fixed continental 6S AOD",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="viridis",
        vmin=0.0,
        vmax=aod_max,
    )
    for index, candidate in enumerate(candidates):
        _panel(
            fig,
            axes[5 + index],
            candidate["aod"],
            f"CCI candidate {index + 1} AOD\n{candidate['short_label']}",
            iy=iy,
            ix=ix,
            radius_px=radius_px,
            cmap="viridis",
            vmin=0.0,
            vmax=aod_max,
        )
    species_display = np.where(selected_index >= 0, selected_index, np.nan)
    _panel(
        fig,
        axes[8],
        species_display,
        "Selected CCI candidate rank\nminimum pooled surface cost",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap=ListedColormap(["#a65d00", "#2879b7", "#2d7a4b"]),
        vmin=-0.5,
        vmax=2.5,
        ticks=[0, 1, 2],
        ticklabels=["CCI-1", "CCI-2", "CCI-3"],
    )
    _panel(
        fig,
        axes[9],
        selected_aod,
        "CCI per-pixel selected AOD",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="viridis",
        vmin=0.0,
        vmax=aod_max,
    )
    _panel(
        fig,
        axes[10],
        delta,
        "CCI selected - fixed 6S AOD",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="RdBu_r",
        vmin=-delta_max,
        vmax=delta_max,
    )
    _panel(
        fig,
        axes[11],
        selected_cost,
        "Selected minimum pooled\nsurface-only cost",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="magma",
        vmin=0.0,
        vmax=float(np.nanpercentile(selected_cost, 98)),
    )
    _panel(
        fig,
        axes[12],
        cost_advantage,
        "CCI minimum - fixed 6S\nsurface-only cost",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="RdBu_r",
        vmin=-advantage_max,
        vmax=advantage_max,
    )
    _panel(
        fig,
        axes[13],
        confidence,
        "Candidate separation\nsecond-best - best cost",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="cividis",
        vmin=0.0,
        vmax=confidence_max,
    )
    selected_strong = np.full(selected_index.shape, np.nan, dtype=np.float64)
    for index, candidate in enumerate(candidates):
        selected_strong[selected_index == index] = candidate["mixture"]["fine_strong"]
    _panel(
        fig,
        axes[14],
        selected_strong,
        "Selected strong-absorbing\nfine-mode fraction",
        iy=iy,
        ix=ix,
        radius_px=radius_px,
        cmap="inferno",
        vmin=0.0,
        vmax=max(0.5, float(np.nanmax(selected_strong))),
    )
    fig.suptitle(
        f"{matchup_id} | exact final-pass spatial species evidence",
        fontsize=13,
        x=0.02,
        ha="left",
    )
    fig.tight_layout(rect=(0.005, 0.005, 1, 0.965))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=125, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _save_diagnostic_figure(
    path: Path,
    *,
    matchup_id: str,
    truth: float,
    records: dict[str, dict[str, Any]],
    fixed: dict[str, Any],
    candidates: list[dict[str, Any]],
    climatology_target: dict[str, float],
    selected_residual: np.ndarray,
    local_mask: np.ndarray,
) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = fixed["valid"]
    local_valid = local_mask & valid
    curve_mask = local_valid if int(np.count_nonzero(local_valid)) >= 20 else valid
    curve_scope = "1.5 km station window" if curve_mask is local_valid else "scene support"
    axis = fixed["axis"]
    threshold = 0.05 + 0.15 * truth

    fig, axes = plt.subplots(3, 3, figsize=(16.4, 13.3), facecolor="white")
    axes = axes.ravel()
    axes[0].plot(
        axis,
        _normalise_curve(fixed["scene_curve"]),
        color=METHOD_COLORS["sixs_continental"],
        label="Fixed 6S",
    )
    axes[1].plot(
        axis,
        _normalise_curve(fixed["local_curve"]),
        color=METHOD_COLORS["sixs_continental"],
        label="Fixed 6S",
    )
    for index, candidate in enumerate(candidates):
        color = ["#a65d00", "#2879b7", "#2d7a4b"][index]
        axes[0].plot(
            axis, _normalise_curve(candidate["scene_curve"]), color=color, label=f"CCI-{index + 1}"
        )
        axes[1].plot(
            axis, _normalise_curve(candidate["local_curve"]), color=color, label=f"CCI-{index + 1}"
        )
    for ax, title in zip(
        axes[:2],
        ["Scene median surface cost", f"Local median surface cost: {curve_scope}"],
        strict=True,
    ):
        ax.axvline(truth, color="#111111", linewidth=1.3, label="AERONET")
        ax.set_title(title)
        ax.set_xlabel("AOD")
        ax.set_ylabel("log10(1 + cost above minimum)")
        ax.set_xlim(axis[0], min(axis[-1], max(1.5, truth * 1.45)))
        ax.grid(alpha=0.18)
        ax.legend(fontsize=7.5, ncol=2)

    mixture_rows = [climatology_target, *[item["mixture"] for item in candidates]]
    x = np.arange(len(mixture_rows))
    bottoms = np.zeros(len(mixture_rows))
    for component in COMPONENTS:
        values = np.asarray([item[component] for item in mixture_rows])
        axes[2].bar(
            x,
            values,
            bottom=bottoms,
            color=COMPONENT_COLORS[component],
            label=COMPONENT_LABELS[component],
        )
        bottoms += values
    for index, candidate in enumerate(candidates, start=1):
        axes[2].text(
            index,
            1.025,
            f"selected {100.0 * candidate['selected_fraction']:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    axes[2].set_xticks(x, ["CCI target", "CCI-1", "CCI-2", "CCI-3"])
    axes[2].set_ylim(0.0, 1.14)
    axes[2].set_ylabel("Aerosol component fraction")
    axes[2].set_title("Monthly CCI target, candidate mixtures, and selected share")
    axes[2].legend(fontsize=7, ncol=2, loc="lower center")
    axes[2].grid(axis="y", alpha=0.18)

    method_residuals = [fixed["solution_residual"]]
    residual_labels = ["Fixed 6S"]
    residual_colors = [METHOD_COLORS["sixs_continental"]]
    selected_local = selected_residual[:, curve_mask]
    method_residuals.append(np.nanmedian(selected_local, axis=1))
    residual_labels.append("CCI selected")
    residual_colors.append(METHOD_COLORS["sixs_cci3"])
    positions = np.arange(len(BANDS))
    width = 0.34
    for index, (values, label, color) in enumerate(
        zip(method_residuals, residual_labels, residual_colors, strict=True)
    ):
        axes[3].bar(positions + (index - 0.5) * width, values, width, color=color, label=label)
    axes[3].set_xticks(positions, BANDS)
    axes[3].set_ylabel("Absolute BOA - prior residual")
    axes[3].set_title(f"Spectral closure at each pixel's solution\n{curve_scope}")
    axes[3].legend(fontsize=8)
    axes[3].grid(axis="y", alpha=0.18)

    method_order = list(PRIMARY_METHODS)
    output_rows = []
    for method in method_order:
        record = records[method]
        pass1 = _finite((record.get("anchor_iterate") or {}).get("pass1_scene_mean"))
        final = _retrieved(record)
        if pass1 is not None:
            output_rows.append(
                (f"{METHOD_SHORT[method]} pass 1", pass1, METHOD_COLORS[method], "o")
            )
        if final is not None:
            output_rows.append((f"{METHOD_SHORT[method]} final", final, METHOD_COLORS[method], "s"))
    for y, (_label, value, color, marker) in enumerate(output_rows):
        axes[4].plot([truth, value], [y, y], color="#c8ced2", linewidth=0.9)
        axes[4].scatter([value], [y], color=color, marker=marker, s=36, zorder=3)
    axes[4].axvspan(truth - threshold, truth + threshold, color="#eceff1", alpha=0.9)
    axes[4].axvline(truth, color="#111111", linewidth=1.4)
    axes[4].set_yticks(range(len(output_rows)), [row[0] for row in output_rows], fontsize=7.5)
    axes[4].invert_yaxis()
    axes[4].set_xlabel("AOD")
    axes[4].set_title("Pass-1 anchor and final retrieval outputs")
    axes[4].grid(axis="x", alpha=0.18)

    wavelengths = np.asarray([WAVELENGTHS[band] for band in BANDS], dtype=np.float64)
    toa = np.asarray([np.nanmedian(fixed["toa"][index][curve_mask]) for index in range(3)])
    axes[5].plot(wavelengths, toa, "o-", color="#222222", label="Target TOA")
    for method, source in [
        ("sixs_continental", fixed),
        ("sixs_cci3", candidates[0]),
    ]:
        prior = np.asarray(
            [np.nanmedian(source["boa_prior"][index][curve_mask]) for index in range(3)]
        )
        uncertainty = np.asarray(
            [np.nanmedian(source["boa_unc"][index][curve_mask]) for index in range(3)]
        )
        axes[5].errorbar(
            wavelengths,
            prior,
            yerr=uncertainty,
            fmt="o-",
            capsize=3,
            color=METHOD_COLORS[method],
            label=f"{METHOD_SHORT[method]} BOA prior",
        )
    axes[5].set_xticks(wavelengths, BANDS)
    axes[5].set_ylabel("Reflectance")
    axes[5].set_title("Final-pass surface-prior spectrum")
    axes[5].legend(fontsize=7.5)
    axes[5].grid(alpha=0.18)

    fixed_final = _retrieved(records["sixs_continental"])
    for index, band in enumerate(BANDS):
        axes[6].plot(
            axis,
            _normalise_curve(fixed["local_band_cost"][index]),
            color=BAND_COLORS[band],
            label=band,
        )
    axes[6].axvline(truth, color="#111111", linewidth=1.3, label="AERONET")
    if fixed_final is not None:
        axes[6].axvline(
            fixed_final,
            color=METHOD_COLORS["sixs_continental"],
            linestyle="--",
            linewidth=1.2,
            label="Fixed final",
        )
    axes[6].set_xlim(axis[0], min(axis[-1], max(1.5, truth * 1.45)))
    axes[6].set_xlabel("AOD")
    axes[6].set_ylabel("log10(1 + cost above minimum)")
    axes[6].set_title(f"Fixed-6S per-band cost: {curve_scope}")
    axes[6].legend(fontsize=7.5, ncol=2)
    axes[6].grid(alpha=0.18)

    candidate_colors = ["#a65d00", "#2879b7", "#2d7a4b"]
    band_styles = ["-", "--", ":"]
    for candidate_index, candidate in enumerate(candidates):
        for band_index, _band in enumerate(BANDS):
            axes[7].plot(
                axis,
                _normalise_curve(candidate["local_band_cost"][band_index]),
                color=candidate_colors[candidate_index],
                linestyle=band_styles[band_index],
                linewidth=1.25,
                label=f"CCI-{candidate_index + 1}" if band_index == 0 else None,
            )
    for band_index, band in enumerate(BANDS):
        axes[7].plot([], [], color="#333333", linestyle=band_styles[band_index], label=band)
    axes[7].axvline(truth, color="#111111", linewidth=1.3, label="AERONET")
    axes[7].set_xlim(axis[0], min(axis[-1], max(1.5, truth * 1.45)))
    axes[7].set_xlabel("AOD")
    axes[7].set_ylabel("log10(1 + cost above minimum)")
    axes[7].set_title(f"CCI-candidate per-band costs: {curve_scope}")
    axes[7].legend(fontsize=7, ncol=3)
    axes[7].grid(alpha=0.18)

    truth_residuals = [fixed["truth_residual"]] + [
        candidate["truth_residual"] for candidate in candidates
    ]
    truth_labels = ["Fixed 6S", "CCI-1", "CCI-2", "CCI-3"]
    truth_colors = [METHOD_COLORS["sixs_continental"], *candidate_colors]
    positions = np.arange(len(BANDS))
    width = 0.2
    for index, (values, label, color) in enumerate(
        zip(truth_residuals, truth_labels, truth_colors, strict=True)
    ):
        axes[8].bar(
            positions + (index - 1.5) * width,
            values,
            width,
            color=color,
            label=label,
        )
    axes[8].set_xticks(positions, BANDS)
    axes[8].set_ylabel("Absolute BOA - prior residual")
    axes[8].set_title(f"Spectral closure at AERONET truth node\n{curve_scope}")
    axes[8].legend(fontsize=7, ncol=2)
    axes[8].grid(axis="y", alpha=0.18)

    for ax in axes:
        ax.tick_params(axis="both", labelsize=8)
        ax.title.set_fontsize(10)
        ax.xaxis.label.set_fontsize(9)
        ax.yaxis.label.set_fontsize(9)
    fig.suptitle(
        f"{matchup_id} | spectral, cost, iteration, and mixture evidence",
        fontsize=13,
        x=0.02,
        ha="left",
    )
    fig.tight_layout(rect=(0.01, 0.01, 0.995, 0.97), pad=1.2)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=125, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return {
        "curve_scope": curve_scope,
        "curve_support_count": int(np.count_nonzero(curve_mask)),
        "station_window_valid_count": int(np.count_nonzero(local_valid)),
    }


def _load_case_costs(
    *,
    matchup_id: str,
    record: dict[str, Any],
    fixed_path: Path,
    cci_paths: Sequence[Path],
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    mixtures = candidate_fraction_sets(
        float(record["lon"]),
        float(record["lat"]),
        _month(matchup_id),
        n=len(cci_paths),
    )
    target_percentages = climatology_fraction_percentages(
        float(record["lon"]),
        float(record["lat"]),
        _month(matchup_id),
    )
    climatology_target = {
        component: float(target_percentages[index] / 100.0)
        for index, component in enumerate(COMPONENTS)
    }
    with np.load(fixed_path, allow_pickle=False) as cost:
        fixed_aod, fixed_unc, fixed_cost = _pooled(cost)
        axis = np.asarray(cost["aot_axis"], dtype=np.float64)
        valid = np.asarray(cost["solve_valid"], dtype=bool)
        fixed_curve_scene = _curve(cost["cube"], valid)
        x = np.asarray(cost["x"], dtype=np.float64)
        y = np.asarray(cost["y"], dtype=np.float64)
        local_mask, iy, ix = _window_mask(
            x,
            y,
            lon=float(record["lon"]),
            lat=float(record["lat"]),
            crs=str(record["scene_crs"]),
        )
        local_valid = local_mask & valid
        curve_mask = local_valid if int(np.count_nonzero(local_valid)) >= 20 else valid
        fixed_curve_local = _curve(cost["cube"], curve_mask)
        fixed_band_cost = np.asarray(cost["band_cost_cube"], dtype=np.float64)
        fixed_band_residual = np.asarray(cost["band_residual_cube"], dtype=np.float64)
        truth_index = int(np.argmin(np.abs(axis - float(record["truth"]))))
        fixed_residual = _candidate_residual_at_solution(
            fixed_band_residual,
            axis,
            fixed_aod,
        )
        fixed = {
            "aod": fixed_aod,
            "unc": fixed_unc,
            "cost": fixed_cost,
            "axis": axis,
            "valid": valid,
            "aot_prior": np.asarray(cost["aot_prior"], dtype=np.float64),
            "aot_prior_unc": np.asarray(cost["aot_prior_unc"], dtype=np.float64),
            "boa_prior": np.asarray(cost["boa_prior"], dtype=np.float64),
            "boa_unc": np.asarray(cost["boa_unc"], dtype=np.float64),
            "toa": np.asarray(cost["toa"], dtype=np.float64),
            "x": x,
            "y": y,
            "scene_curve": fixed_curve_scene,
            "local_curve": fixed_curve_local,
            "local_band_cost": [_curve(fixed_band_cost[index], curve_mask) for index in range(3)],
            "solution_residual": np.nanmedian(fixed_residual[:, curve_mask], axis=1),
            "truth_residual": np.nanmedian(
                fixed_band_residual[:, truth_index][:, curve_mask], axis=1
            ),
            "pool_window": int(np.asarray(cost["pool_window"]).item()),
            "min_count": int(np.asarray(cost["min_count"]).item()),
        }

    candidates = []
    for index, (path, mixture) in enumerate(zip(cci_paths, mixtures, strict=True)):
        with np.load(path, allow_pickle=False) as cost:
            aod, uncertainty, obs_cost = _pooled(cost)
            candidate_band_cost = np.asarray(cost["band_cost_cube"], dtype=np.float64)
            candidate_band_residual = np.asarray(cost["band_residual_cube"], dtype=np.float64)
            candidate = {
                "index": index,
                "mixture": mixture,
                "label": _candidate_label(index, mixture),
                "short_label": _compact_mixture_label(mixture),
                "climatology_l1_distance": float(
                    sum(abs(mixture[name] - climatology_target[name]) for name in COMPONENTS)
                ),
                "aod": aod,
                "unc": uncertainty,
                "cost": obs_cost,
                "scene_curve": _curve(cost["cube"], valid),
                "local_curve": _curve(cost["cube"], curve_mask),
                "local_band_cost": [
                    _curve(candidate_band_cost[band_index], curve_mask) for band_index in range(3)
                ],
                "truth_residual": np.nanmedian(
                    candidate_band_residual[:, truth_index][:, curve_mask], axis=1
                ),
                "boa_prior": np.asarray(cost["boa_prior"], dtype=np.float64),
                "boa_unc": np.asarray(cost["boa_unc"], dtype=np.float64),
            }
        candidates.append(candidate)

    costs = np.stack([item["cost"] for item in candidates], axis=0)
    finite_cost = np.where(np.isfinite(costs), costs, np.inf)
    selected_index = np.argmin(finite_cost, axis=0).astype(np.int16)
    valid_selection = np.any(np.isfinite(costs), axis=0)
    selected_index[~valid_selection] = -1
    rows, columns = np.indices(selected_index.shape)
    safe_index = np.maximum(selected_index, 0)
    aod_stack = np.stack([item["aod"] for item in candidates], axis=0)
    selected_aod = aod_stack[safe_index, rows, columns]
    selected_aod[~valid_selection] = np.nan
    selected_cost = costs[safe_index, rows, columns]
    selected_cost[~valid_selection] = np.nan
    sorted_cost = np.sort(finite_cost, axis=0)
    with np.errstate(invalid="ignore"):
        confidence = sorted_cost[1] - sorted_cost[0]
    confidence[~np.isfinite(confidence)] = np.nan

    selected_residual = np.full((3, *selected_index.shape), np.nan, dtype=np.float64)
    for index, path in enumerate(cci_paths):
        with np.load(path, allow_pickle=False) as cost:
            residual = _candidate_residual_at_solution(
                np.asarray(cost["band_residual_cube"], dtype=np.float64),
                np.asarray(cost["aot_axis"], dtype=np.float64),
                candidates[index]["aod"],
            )
        use = selected_index == index
        selected_residual[:, use] = residual[:, use]

    selected_count = int(np.count_nonzero(valid_selection))
    for index, candidate in enumerate(candidates):
        count = int(np.count_nonzero(selected_index == index))
        candidate["selected_count"] = count
        candidate["selected_fraction"] = count / selected_count if selected_count else 0.0
        candidate["scene_mean_aod"] = float(np.nanmean(candidate["aod"]))
        candidate["median_surface_cost"] = float(np.nanmedian(candidate["cost"]))
        candidate["local_selected_count"] = int(
            np.count_nonzero((selected_index == index) & local_mask)
        )

    replayed = float(np.nanmean(selected_aod)) if np.any(np.isfinite(selected_aod)) else None
    diagnostics = {
        "selected_index": selected_index,
        "selected_aod": selected_aod,
        "selected_cost": selected_cost,
        "selected_residual": selected_residual,
        "confidence": confidence,
        "local_mask": local_mask,
        "iy": iy,
        "ix": ix,
        "selected_count": selected_count,
        "replayed_scene_mean_aod": replayed,
        "replay_delta": (
            replayed - float(record["retrieved"])
            if replayed is not None and record.get("retrieved") is not None
            else None
        ),
        "solver_backstop_median": float(np.nanmedian(fixed["aot_prior"][valid])),
        "solver_backstop_unc_median": float(np.nanmedian(fixed["aot_prior_unc"][valid])),
        "selection_confidence_median": float(np.nanmedian(confidence)),
        "selection_confidence_p10": float(np.nanpercentile(confidence, 10)),
        "climatology_target": climatology_target,
        "cci_minus_fixed_map_mean": float(np.nanmean(selected_aod - fixed["aod"])),
        "cci_minus_fixed_map_median": float(np.nanmedian(selected_aod - fixed["aod"])),
    }
    return fixed, candidates, diagnostics


def _method_case(record: dict[str, Any] | None, truth: float) -> dict[str, Any]:
    value = _retrieved(record)
    anchor = (record or {}).get("anchor_iterate") or {}
    error = value - truth if value is not None else None
    threshold = 0.05 + 0.15 * truth
    return {
        "status": str((record or {}).get("status", "MISSING")),
        "aod": value,
        "error": error,
        "abs_error": abs(error) if error is not None else None,
        "error_ee": error / threshold if error is not None else None,
        "within_ee": _within_ee(value, truth),
        "pass1_aod": _finite(anchor.get("pass1_scene_mean")),
        "pass2_aod": _finite(anchor.get("pass2_scene_mean")),
        "runtime_s": _finite((record or {}).get("runtime_s")),
        "station_aod": _finite(((record or {}).get("retrieval_extraction") or {}).get("station")),
        "window_median_aod": _finite(
            ((record or {}).get("retrieval_extraction") or {}).get("winmed")
        ),
        "quality_score": _finite((record or {}).get("aod_quality_score")),
        "solver": (record or {}).get("solver") or {},
        "prior_boa": (record or {}).get("prior_boa") or {},
        "prior_unc": (record or {}).get("prior_unc") or {},
    }


def _transition(base_hit: bool, candidate_hit: bool) -> str:
    if not base_hit and candidate_hit:
        return "gain"
    if base_hit and not candidate_hit:
        return "loss"
    return "both_hit" if base_hit else "both_miss"


def _case_cache_path(output: Path, matchup_id: str) -> Path:
    return output / "data" / "case-cache" / f"{matchup_id}.json"


def _build_case(
    *,
    output: Path,
    matchup_id: str,
    method_records: dict[str, dict[str, dict[str, Any]]],
    failure_groups: dict[str, dict[str, Any]],
    force: bool,
) -> dict[str, Any]:
    cache_path = _case_cache_path(output, matchup_id)
    spatial_path = output / "assets" / "spatial" / f"{matchup_id}.png"
    diagnostic_path = output / "assets" / "diagnostic" / f"{matchup_id}.png"
    if not force and cache_path.exists() and spatial_path.exists() and diagnostic_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))

    records = {method: rows.get(matchup_id) for method, rows in method_records.items()}
    baseline = records["lut_continental"]
    fixed_record = records["sixs_continental"]
    cci_record = records["sixs_cci3"]
    if baseline is None or fixed_record is None or cci_record is None:
        raise ValueError(f"Primary record missing for {matchup_id}")
    truth = float(baseline["truth"])
    fixed_path = FIXED_COSTS / f"{matchup_id}.npz"
    cci_paths = sorted(CCI_COSTS.glob(f"{matchup_id}.species*.npz"))
    if not fixed_path.exists() or len(cci_paths) != 3:
        raise ValueError(
            f"Cost artifacts incomplete for {matchup_id}: fixed={fixed_path.exists()} cci={len(cci_paths)}"
        )
    fixed, candidates, diagnostics = _load_case_costs(
        matchup_id=matchup_id,
        record=cci_record,
        fixed_path=fixed_path,
        cci_paths=cci_paths,
    )
    replay_delta = float(diagnostics["replay_delta"])
    if abs(replay_delta) > 1e-2:
        raise ValueError(
            f"Species replay mismatch for {matchup_id}: {replay_delta:.9g}"
        )
    if abs(replay_delta) > SPECIES_REPLAY_DELTA_TOLERANCE:
        print(
            f"Species replay tolerance warning for {matchup_id}: {replay_delta:.9g} > "
            f"{SPECIES_REPLAY_DELTA_TOLERANCE:.1e}"
        )
    _save_spatial_figure(
        spatial_path,
        matchup_id=matchup_id,
        truth=truth,
        fixed=fixed,
        candidates=candidates,
        selected_index=diagnostics["selected_index"],
        selected_aod=diagnostics["selected_aod"],
        selected_cost=diagnostics["selected_cost"],
        confidence=diagnostics["confidence"],
        iy=diagnostics["iy"],
        ix=diagnostics["ix"],
    )
    curve_metadata = _save_diagnostic_figure(
        diagnostic_path,
        matchup_id=matchup_id,
        truth=truth,
        records={method: record for method, record in records.items() if record is not None},
        fixed=fixed,
        candidates=candidates,
        climatology_target=diagnostics["climatology_target"],
        selected_residual=diagnostics["selected_residual"],
        local_mask=diagnostics["local_mask"],
    )
    method_data = {method: _method_case(records.get(method), truth) for method in METHODS}
    fixed_hit = method_data["sixs_continental"]["within_ee"]
    cci_hit = method_data["sixs_cci3"]["within_ee"]
    lut_hit = method_data["lut_continental"]["within_ee"]
    group = failure_groups.get(matchup_id)
    cost_urls = {
        "fixed_6s": f"../../{FIXED_COSTS.name}/{matchup_id}.npz",
        "cci_1": f"../../{CCI_COSTS.name}/{matchup_id}.species0.npz",
        "cci_2": f"../../{CCI_COSTS.name}/{matchup_id}.species1.npz",
        "cci_3": f"../../{CCI_COSTS.name}/{matchup_id}.species2.npz",
    }
    if (EXACT_CCI_COSTS / f"{matchup_id}.species0.npz").exists():
        cost_urls["exact_cci"] = f"../../{EXACT_CCI_COSTS.name}/{matchup_id}.species0.npz"
    if (DIRECT_COSTS / f"{matchup_id}.npz").exists():
        cost_urls["direct_libradtran"] = f"../../{DIRECT_COSTS.name}/{matchup_id}.npz"
    case = {
        "matchup_id": matchup_id,
        "site": str(baseline.get("site", matchup_id.split("__", 1)[0])),
        "lon": float(baseline["lon"]),
        "lat": float(baseline["lat"]),
        "truth_aod": truth,
        "ee_threshold": 0.05 + 0.15 * truth,
        "cloud_fraction": float(baseline.get("cloud_frac", 0.0)),
        "invalid_fraction": float(baseline.get("invalid_frac", 0.0)),
        "regime": str(baseline.get("regime", "")),
        "baseline_failure_group": str(group.get("mechanism_code")) if group else "H",
        "baseline_failure_label": str(group.get("mechanism")) if group else "Baseline within EE",
        "mask_fallback": bool(
            method_data["sixs_continental"]["solver"].get("surface_cloud_mask_bypassed")
            or method_data["sixs_continental"]["solver"].get("surface_water_mask_bypassed")
            or method_data["sixs_cci3"]["solver"].get("surface_cloud_mask_bypassed")
            or method_data["sixs_cci3"]["solver"].get("surface_water_mask_bypassed")
        ),
        "methods": method_data,
        "transition_cci_vs_fixed": _transition(fixed_hit, cci_hit),
        "transition_fixed_vs_lut": _transition(lut_hit, fixed_hit),
        "transition_cci_vs_lut": _transition(lut_hit, cci_hit),
        "cci_minus_fixed_abs_error": (
            method_data["sixs_cci3"]["abs_error"] - method_data["sixs_continental"]["abs_error"]
        ),
        "fixed_minus_lut_abs_error": (
            method_data["sixs_continental"]["abs_error"]
            - method_data["lut_continental"]["abs_error"]
        ),
        "species": {
            "candidate_count": len(candidates),
            "selected_pixel_count": diagnostics["selected_count"],
            "replayed_scene_mean_aod": diagnostics["replayed_scene_mean_aod"],
            "reported_scene_mean_aod": cci_record["retrieved"],
            "replay_delta": diagnostics["replay_delta"],
            "selection_confidence_median": diagnostics["selection_confidence_median"],
            "selection_confidence_p10": diagnostics["selection_confidence_p10"],
            "climatology_target": diagnostics["climatology_target"],
            "climatology_target_text": _component_text(diagnostics["climatology_target"]),
            "cci_minus_fixed_map_mean": diagnostics["cci_minus_fixed_map_mean"],
            "cci_minus_fixed_map_median": diagnostics["cci_minus_fixed_map_median"],
            "candidates": [
                {
                    "index": item["index"],
                    "label": item["label"],
                    "mixture": item["mixture"],
                    "mixture_text": _component_text(item["mixture"]),
                    "climatology_l1_distance": item["climatology_l1_distance"],
                    "selected_count": item["selected_count"],
                    "selected_fraction": item["selected_fraction"],
                    "local_selected_count": item["local_selected_count"],
                    "scene_mean_aod": item["scene_mean_aod"],
                    "median_surface_cost": item["median_surface_cost"],
                }
                for item in candidates
            ],
        },
        "cost_grid": {
            "shape": list(fixed["aod"].shape),
            "aod_nodes": int(fixed["axis"].size),
            "aod_min": float(fixed["axis"].min()),
            "aod_max": float(fixed["axis"].max()),
            "pool_window": fixed["pool_window"],
            "pool_min_count": fixed["min_count"],
            "solver_backstop_median": diagnostics["solver_backstop_median"],
            "solver_backstop_unc_median": diagnostics["solver_backstop_unc_median"],
            **curve_metadata,
        },
        "spatial_image": f"assets/spatial/{matchup_id}.png",
        "diagnostic_image": f"assets/diagnostic/{matchup_id}.png",
        "source_urls": {
            method: f"../../{Path(METHODS[method]).name}/{matchup_id}.json"
            for method in METHODS
            if records.get(method) is not None
        },
        "cost_urls": cost_urls,
    }
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(_jsonable(case), separators=(",", ":")), encoding="utf-8")
    return _jsonable(case)


def _stratum_rows(
    cases: Sequence[dict[str, Any]],
    field: str,
    labels: Sequence[tuple[str, Any]],
) -> list[dict[str, Any]]:
    output = []
    for label, predicate in labels:
        rows = [row for row in cases if predicate(row[field])]
        item = {"label": label, "count": len(rows), "methods": {}}
        for method in PRIMARY_METHODS:
            hits = sum(bool(row["methods"][method]["within_ee"]) for row in rows)
            item["methods"][method] = {
                "hits": hits,
                "rate": hits / len(rows) if rows else None,
            }
        output.append(item)
    return output


def _stratification(cases: Sequence[dict[str, Any]]) -> dict[str, Any]:
    truth_bins = [
        ("<0.15", lambda value: value < 0.15),
        ("0.15-0.30", lambda value: 0.15 <= value < 0.30),
        ("0.30-0.60", lambda value: 0.30 <= value < 0.60),
        ("0.60-1.00", lambda value: 0.60 <= value < 1.00),
        (">=1.00", lambda value: value >= 1.00),
    ]
    cloud_bins = [
        ("<2%", lambda value: value < 0.02),
        ("2-5%", lambda value: 0.02 <= value < 0.05),
        ("5-10%", lambda value: 0.05 <= value < 0.10),
        ("10-20%", lambda value: 0.10 <= value < 0.20),
    ]
    groups = sorted({row["baseline_failure_group"] for row in cases})
    group_rows = []
    for group in groups:
        rows = [row for row in cases if row["baseline_failure_group"] == group]
        group_rows.append(
            {
                "label": group,
                "description": rows[0]["baseline_failure_label"],
                "count": len(rows),
                "methods": {
                    method: {
                        "hits": sum(row["methods"][method]["within_ee"] for row in rows),
                        "rate": (
                            sum(row["methods"][method]["within_ee"] for row in rows) / len(rows)
                            if rows
                            else None
                        ),
                    }
                    for method in PRIMARY_METHODS
                },
            }
        )
    return {
        "truth": _stratum_rows(cases, "truth_aod", truth_bins),
        "cloud": _stratum_rows(cases, "cloud_fraction", cloud_bins),
        "baseline_groups": group_rows,
    }


def _value_summary(values: Sequence[float], truths: Sequence[float]) -> dict[str, Any]:
    retrieved = np.asarray(values, dtype=np.float64)
    reference = np.asarray(truths, dtype=np.float64)
    errors = retrieved - reference
    thresholds = 0.05 + 0.15 * reference
    hits = int(np.count_nonzero(np.abs(errors) <= thresholds + 1e-12))
    return {
        "count": int(retrieved.size),
        "hits": hits,
        "rate": hits / retrieved.size if retrieved.size else None,
        "rmse": float(np.sqrt(np.mean(np.square(errors)))) if errors.size else None,
        "mae": float(np.mean(np.abs(errors))) if errors.size else None,
        "bias": float(np.mean(errors)) if errors.size else None,
    }


def _final_pass_policy_replays(cases: Sequence[dict[str, Any]]) -> dict[str, Any]:
    truths = [float(row["truth_aod"]) for row in cases]
    values: dict[str, list[float]] = {
        "nearest_cci": [],
        "scene_min_surface_cost": [],
        "pixel_min_surface_cost": [],
        "cci_candidate_oracle": [],
        "fixed_plus_cci_oracle": [],
    }
    for row in cases:
        truth = float(row["truth_aod"])
        candidate_rows = row["species"]["candidates"]
        candidate_aod = [float(candidate["scene_mean_aod"]) for candidate in candidate_rows]
        scene_min = min(candidate_rows, key=lambda candidate: candidate["median_surface_cost"])
        fixed = float(row["methods"]["sixs_continental"]["aod"])
        values["nearest_cci"].append(candidate_aod[0])
        values["scene_min_surface_cost"].append(float(scene_min["scene_mean_aod"]))
        values["pixel_min_surface_cost"].append(float(row["methods"]["sixs_cci3"]["aod"]))
        values["cci_candidate_oracle"].append(
            min(candidate_aod, key=lambda value: abs(value - truth))
        )
        values["fixed_plus_cci_oracle"].append(
            min([fixed, *candidate_aod], key=lambda value: abs(value - truth))
        )
    labels = {
        "nearest_cci": "Nearest climatological CCI mixture only",
        "scene_min_surface_cost": "One scene-wide CCI mixture by minimum median surface cost",
        "pixel_min_surface_cost": "Current per-pixel CCI minimum surface cost",
        "cci_candidate_oracle": "Diagnostic best of three CCI candidates by AERONET",
        "fixed_plus_cci_oracle": "Diagnostic best of fixed 6S plus three CCI candidates by AERONET",
    }
    operational = {
        "nearest_cci": True,
        "scene_min_surface_cost": True,
        "pixel_min_surface_cost": True,
        "cci_candidate_oracle": False,
        "fixed_plus_cci_oracle": False,
    }
    return {
        key: {
            "label": labels[key],
            "operational_without_aeronet": operational[key],
            **_value_summary(policy_values, truths),
        }
        for key, policy_values in values.items()
    }


def _save_global_figures(output: Path, report: dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cases = report["cases"]
    figure_dir = output / "assets" / "global"
    figure_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.6), facecolor="white")
    axes = axes.ravel()
    max_aod = (
        max(
            1.5,
            max(row["truth_aod"] for row in cases),
            max(
                row["methods"][method]["aod"] or 0.0 for row in cases for method in PRIMARY_METHODS
            ),
        )
        * 1.04
    )
    grid = np.linspace(0.0, max_aod, 300)
    for ax, method in zip(axes[:3], PRIMARY_METHODS, strict=True):
        truth = np.asarray([row["truth_aod"] for row in cases])
        retrieved = np.asarray([row["methods"][method]["aod"] for row in cases])
        hit = np.asarray([row["methods"][method]["within_ee"] for row in cases])
        ax.fill_between(
            grid,
            grid - (0.05 + 0.15 * grid),
            grid + (0.05 + 0.15 * grid),
            color="#e8ecef",
        )
        ax.scatter(truth[hit], retrieved[hit], s=18, color=METHOD_COLORS[method], alpha=0.65)
        ax.scatter(
            truth[~hit],
            retrieved[~hit],
            s=24,
            facecolor="white",
            edgecolor="#a33a3a",
            linewidth=0.9,
        )
        ax.plot(grid, grid, color="#111111", linewidth=1.0)
        summary = report["metrics"][method]
        ax.set_title(
            f"{METHOD_SHORT[method]}\n{summary['hits']}/{summary['expected']} within EE; "
            f"RMSE {summary['rmse']:.3f}"
        )
        ax.set(xlabel="AERONET AOD", ylabel="Retrieved AOD", xlim=(0, max_aod), ylim=(0, max_aod))
        ax.grid(alpha=0.16)

    comparisons = [
        ("fixed_minus_lut_abs_error", "6S fixed - LUT |error|", "transition_fixed_vs_lut"),
        ("cci_minus_fixed_abs_error", "CCI-3 - fixed 6S |error|", "transition_cci_vs_fixed"),
    ]
    transition_colors = {
        "gain": "#2d7a4b",
        "loss": "#a33a3a",
        "both_hit": "#2879b7",
        "both_miss": "#7b858b",
    }
    for ax, (field, title, transition_field) in zip(axes[3:5], comparisons, strict=True):
        for transition, color in transition_colors.items():
            rows = [row for row in cases if row[transition_field] == transition]
            ax.scatter(
                [row["truth_aod"] for row in rows],
                [row[field] for row in rows],
                s=25,
                color=color,
                alpha=0.75,
                label=transition.replace("_", " "),
            )
        ax.axhline(0.0, color="#111111", linewidth=1.0)
        ax.set(xlabel="AERONET AOD", ylabel="Change in absolute error")
        ax.set_title(f"{title}\nnegative values improve absolute error")
        ax.grid(alpha=0.16)
        ax.legend(fontsize=7, ncol=2)

    order = ["gain", "loss", "both_hit", "both_miss"]
    fixed_counts = Counter(row["transition_fixed_vs_lut"] for row in cases)
    cci_counts = Counter(row["transition_cci_vs_fixed"] for row in cases)
    x = np.arange(len(order))
    width = 0.36
    axes[5].bar(
        x - width / 2,
        [fixed_counts[key] for key in order],
        width,
        color="#1565a8",
        label="6S fixed vs LUT",
    )
    axes[5].bar(
        x + width / 2,
        [cci_counts[key] for key in order],
        width,
        color="#2d7a4b",
        label="CCI-3 vs fixed 6S",
    )
    axes[5].set_xticks(x, [value.replace("_", "\n") for value in order])
    axes[5].set_ylabel("Matchups")
    axes[5].set_title("Expected-error transitions")
    axes[5].legend(fontsize=8)
    axes[5].grid(axis="y", alpha=0.16)
    fig.suptitle("Full-cohort retrieval performance and paired effects", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(figure_dir / "performance_overview.png", dpi=135, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(15.5, 10.8), facecolor="white")
    strata = report["stratification"]
    for ax, key, title in [
        (axes[0, 0], "truth", "Within-EE rate by AERONET AOD"),
        (axes[0, 1], "cloud", "Within-EE rate by image cloud fraction"),
    ]:
        rows = strata[key]
        x = np.arange(len(rows))
        width = 0.25
        for index, method in enumerate(PRIMARY_METHODS):
            ax.bar(
                x + (index - 1) * width,
                [row["methods"][method]["rate"] for row in rows],
                width,
                color=METHOD_COLORS[method],
                label=METHOD_SHORT[method],
            )
        ax.set_xticks(x, [f"{row['label']}\n(n={row['count']})" for row in rows])
        ax.set_ylim(0.0, 1.05)
        ax.set_ylabel("Within-EE fraction")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.16)
        ax.legend(fontsize=7)

    group_rows = strata["baseline_groups"]
    x = np.arange(len(group_rows))
    width = 0.25
    for index, method in enumerate(PRIMARY_METHODS):
        axes[1, 0].bar(
            x + (index - 1) * width,
            [row["methods"][method]["rate"] for row in group_rows],
            width,
            color=METHOD_COLORS[method],
            label=METHOD_SHORT[method],
        )
    axes[1, 0].set_xticks(x, [f"{row['label']}\n(n={row['count']})" for row in group_rows])
    axes[1, 0].set_ylim(0.0, 1.05)
    axes[1, 0].set_ylabel("Within-EE fraction")
    axes[1, 0].set_title("Performance by prior baseline failure-evidence group")
    axes[1, 0].grid(axis="y", alpha=0.16)
    axes[1, 0].legend(fontsize=7)

    colors = [transition_colors[row["transition_cci_vs_fixed"]] for row in cases]
    sizes = [25 + 130 * min(abs(row["cci_minus_fixed_abs_error"]), 0.4) / 0.4 for row in cases]
    axes[1, 1].scatter(
        [row["lon"] for row in cases],
        [row["lat"] for row in cases],
        c=colors,
        s=sizes,
        alpha=0.78,
        edgecolor="white",
        linewidth=0.4,
    )
    axes[1, 1].set(xlim=(-180, 180), ylim=(-65, 80), xlabel="Longitude", ylabel="Latitude")
    axes[1, 1].set_xticks(np.arange(-180, 181, 60))
    axes[1, 1].set_yticks(np.arange(-60, 81, 30))
    axes[1, 1].set_title(
        "Spatial distribution of CCI-3 vs fixed-6S transitions\nmarker size follows |absolute-error change|"
    )
    axes[1, 1].grid(alpha=0.2)
    for transition, color in transition_colors.items():
        axes[1, 1].scatter([], [], color=color, label=transition.replace("_", " "), s=30)
    axes[1, 1].legend(fontsize=7, ncol=2)
    fig.suptitle(
        "Stratified and spatial evidence; rates retain the fixed 152-case denominator", fontsize=13
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(figure_dir / "stratified_spatial.png", dpi=135, bbox_inches="tight")
    plt.close(fig)

    selected_fractions = np.asarray(
        [
            [candidate["selected_fraction"] for candidate in row["species"]["candidates"]]
            for row in cases
        ]
    )
    mixtures = np.asarray(
        [
            [
                [candidate["mixture"][component] for component in COMPONENTS]
                for candidate in row["species"]["candidates"]
            ]
            for row in cases
        ]
    )
    climatology_targets = np.asarray(
        [
            [row["species"]["climatology_target"][component] for component in COMPONENTS]
            for row in cases
        ]
    )
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.2), facecolor="white")
    axes[0, 0].boxplot(
        [selected_fractions[:, index] for index in range(3)],
        tick_labels=["CCI-1", "CCI-2", "CCI-3"],
        patch_artist=True,
        boxprops={"facecolor": "#dce6eb"},
        medianprops={"color": "#111111"},
    )
    axes[0, 0].set_ylabel("Selected-pixel fraction per scene")
    axes[0, 0].set_title("Distribution of candidate-rank selection shares")
    axes[0, 0].grid(axis="y", alpha=0.18)
    mean_mixtures = np.vstack([climatology_targets.mean(axis=0), mixtures.mean(axis=0)])
    x = np.arange(4)
    bottom = np.zeros(4)
    for component_index, component in enumerate(COMPONENTS):
        values = mean_mixtures[:, component_index]
        axes[0, 1].bar(
            x,
            values,
            bottom=bottom,
            color=COMPONENT_COLORS[component],
            label=COMPONENT_LABELS[component],
        )
        bottom += values
    axes[0, 1].set_xticks(x, ["CCI target", "CCI-1", "CCI-2", "CCI-3"])
    axes[0, 1].set_ylabel("Mean candidate component fraction")
    axes[0, 1].set_title("Cohort-average CCI target and candidate composition")
    axes[0, 1].legend(fontsize=7, ncol=2)
    axes[0, 1].grid(axis="y", alpha=0.18)

    cci_delta = np.asarray([row["cci_minus_fixed_abs_error"] for row in cases])
    confidence = np.asarray([row["species"]["selection_confidence_median"] for row in cases])
    dominant = np.max(selected_fractions, axis=1)
    transition = [row["transition_cci_vs_fixed"] for row in cases]
    axes[1, 0].scatter(
        confidence,
        cci_delta,
        c=[transition_colors[value] for value in transition],
        s=28,
        alpha=0.75,
    )
    axes[1, 0].axhline(0.0, color="#111111", linewidth=1.0)
    axes[1, 0].set(
        xlabel="Median second-best - best surface cost", ylabel="CCI-3 - fixed 6S |error|"
    )
    axes[1, 0].set_title("Selection separation versus retrieval effect")
    axes[1, 0].grid(alpha=0.18)
    axes[1, 1].scatter(
        dominant,
        cci_delta,
        c=[transition_colors[value] for value in transition],
        s=28,
        alpha=0.75,
    )
    axes[1, 1].axhline(0.0, color="#111111", linewidth=1.0)
    axes[1, 1].set(
        xlabel="Largest candidate selected-pixel fraction", ylabel="CCI-3 - fixed 6S |error|"
    )
    axes[1, 1].set_title("Spatial dominance versus retrieval effect")
    axes[1, 1].grid(alpha=0.18)
    fig.suptitle("CCI candidate composition and selection behavior", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(figure_dir / "species_selection.png", dpi=135, bbox_inches="tight")
    plt.close(fig)

    smoke_ids = set(report["direct_libradtran_smoke"]["matchup_ids"])
    smoke_cases = [row for row in cases if row["matchup_id"] in smoke_ids]
    smoke_methods = (
        "lut_continental",
        "sixs_continental",
        "sixs_cci3",
        "sixs_cci_exact_smoke",
        "libradtran_continental_smoke",
    )
    if smoke_cases:
        fig, axes = plt.subplots(2, 2, figsize=(15.5, 10.5), facecolor="white")
        x = np.arange(len(smoke_cases))
        truth = np.asarray([row["truth_aod"] for row in smoke_cases])
        thresholds = 0.05 + 0.15 * truth
        axes[0, 0].errorbar(
            x,
            truth,
            yerr=thresholds,
            fmt="*",
            color="#111111",
            capsize=4,
            markersize=8,
            label="AERONET and EE",
        )
        offsets = np.linspace(-0.28, 0.28, len(smoke_methods))
        for offset, method in zip(offsets, smoke_methods, strict=True):
            axes[0, 0].scatter(
                x + offset,
                [row["methods"][method]["aod"] for row in smoke_cases],
                color=METHOD_COLORS[method],
                s=30,
                label=METHOD_SHORT[method],
            )
        axes[0, 0].set_ylabel("AOD")
        axes[0, 0].set_title("Retrieval output by predefined smoke case")
        axes[0, 0].grid(axis="y", alpha=0.18)
        axes[0, 0].legend(fontsize=7, ncol=2)

        width = 0.15
        for method_index, method in enumerate(smoke_methods):
            normalized = [row["methods"][method]["error_ee"] for row in smoke_cases]
            axes[0, 1].bar(
                x + (method_index - 2) * width,
                normalized,
                width,
                color=METHOD_COLORS[method],
                label=METHOD_SHORT[method],
            )
        axes[0, 1].axhline(-1.0, color="#777777", linestyle="--", linewidth=0.9)
        axes[0, 1].axhline(1.0, color="#777777", linestyle="--", linewidth=0.9)
        axes[0, 1].axhline(0.0, color="#111111", linewidth=0.8)
        axes[0, 1].set_ylabel("Signed error / EE")
        axes[0, 1].set_title("Normalized retrieval error")
        axes[0, 1].grid(axis="y", alpha=0.18)
        axes[0, 1].legend(fontsize=7, ncol=2)

        iterative_methods = (
            "sixs_continental",
            "sixs_cci3",
            "sixs_cci_exact_smoke",
            "libradtran_continental_smoke",
        )
        maxima = [1.0]
        for method in iterative_methods:
            pass1 = np.asarray(
                [row["methods"][method]["pass1_aod"] for row in smoke_cases], dtype=np.float64
            )
            final = np.asarray(
                [row["methods"][method]["aod"] for row in smoke_cases], dtype=np.float64
            )
            axes[1, 0].scatter(
                pass1,
                final,
                color=METHOD_COLORS[method],
                s=34,
                alpha=0.8,
                label=METHOD_SHORT[method],
            )
            maxima.extend(pass1[np.isfinite(pass1)].tolist())
            maxima.extend(final[np.isfinite(final)].tolist())
        maximum = max(maxima) * 1.08
        axes[1, 0].plot([0.0, maximum], [0.0, maximum], color="#111111", linewidth=0.9)
        axes[1, 0].set(
            xlabel="Pass-1 scene AOD",
            ylabel="Final scene AOD",
            xlim=(0.0, maximum),
            ylim=(0.0, maximum),
        )
        axes[1, 0].set_title("Anchor-iteration movement")
        axes[1, 0].grid(alpha=0.18)
        axes[1, 0].legend(fontsize=7)

        runtime_methods = (
            "lut_continental",
            "sixs_continental",
            "sixs_cci3",
            "sixs_cci_exact_smoke",
            "libradtran_continental_smoke",
        )
        for method_index, method in enumerate(runtime_methods):
            axes[1, 1].bar(
                x + (method_index - 2) * width,
                [row["methods"][method]["runtime_s"] for row in smoke_cases],
                width,
                color=METHOD_COLORS[method],
                label=METHOD_SHORT[method],
            )
        axes[1, 1].set_yscale("log")
        axes[1, 1].set_ylabel("Runtime (seconds, log scale)")
        axes[1, 1].set_title("End-to-end runtime")
        axes[1, 1].grid(axis="y", alpha=0.18)
        axes[1, 1].legend(fontsize=7, ncol=2)

        labels = [row["site"] for row in smoke_cases]
        for ax in (axes[0, 0], axes[0, 1], axes[1, 1]):
            ax.set_xticks(x, labels, rotation=22, ha="right")
        fig.suptitle("Predefined six-case RT and species smoke controls", fontsize=13)
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        fig.savefig(figure_dir / "smoke_controls.png", dpi=135, bbox_inches="tight")
        plt.close(fig)


def _write_index(output: Path) -> None:
    index = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="light dark">
  <title>SIAC aerosol-species experiment</title>
  <link rel="icon" href="data:,">
  <link rel="stylesheet" href="app.css">
</head>
<body>
  <div id="app" class="app-shell" aria-live="polite">
    <div class="loading-state">Loading the aerosol-species experiment snapshot...</div>
  </div>
  <noscript>The report requires JavaScript; the CSV, JSON, and PNG evidence remains available in this directory.</noscript>
  <script src="app.js"></script>
</body>
</html>
"""
    (output / "index.html").write_text(index, encoding="utf-8")


def _write_detailed_case_csv(output: Path, cases: Sequence[dict[str, Any]]) -> Path:
    rows = []
    for case in cases:
        row = {
            "matchup_id": case["matchup_id"],
            "site": case["site"],
            "lon": case["lon"],
            "lat": case["lat"],
            "truth_aod": case["truth_aod"],
            "ee_threshold": case["ee_threshold"],
            "cloud_fraction": case["cloud_fraction"],
            "baseline_failure_group": case["baseline_failure_group"],
            "mask_fallback": case["mask_fallback"],
            "transition_fixed_vs_lut": case["transition_fixed_vs_lut"],
            "transition_cci_vs_fixed": case["transition_cci_vs_fixed"],
            "transition_cci_vs_lut": case["transition_cci_vs_lut"],
            "fixed_minus_lut_abs_error": case["fixed_minus_lut_abs_error"],
            "cci_minus_fixed_abs_error": case["cci_minus_fixed_abs_error"],
            "solver_backstop_median": case["cost_grid"]["solver_backstop_median"],
            "solver_backstop_unc_median": case["cost_grid"]["solver_backstop_unc_median"],
            "species_selection_confidence_median": case["species"]["selection_confidence_median"],
            "species_replay_delta": case["species"]["replay_delta"],
        }
        for component in COMPONENTS:
            row[f"cci_climatology_target_{component}_fraction"] = case["species"][
                "climatology_target"
            ][component]
        for method in METHODS:
            values = case["methods"][method]
            row[f"{method}_aod"] = values["aod"]
            row[f"{method}_error"] = values["error"]
            row[f"{method}_within_ee"] = values["within_ee"]
            row[f"{method}_pass1_aod"] = values["pass1_aod"]
            row[f"{method}_runtime_s"] = values["runtime_s"]
        for candidate in case["species"]["candidates"]:
            prefix = f"cci_candidate_{candidate['index'] + 1}"
            row[f"{prefix}_scene_mean_aod"] = candidate["scene_mean_aod"]
            row[f"{prefix}_median_surface_cost"] = candidate["median_surface_cost"]
            row[f"{prefix}_selected_fraction"] = candidate["selected_fraction"]
            row[f"{prefix}_climatology_l1_distance"] = candidate["climatology_l1_distance"]
            for component in COMPONENTS:
                row[f"{prefix}_{component}_fraction"] = candidate["mixture"][component]
        rows.append(row)
    path = output / "data" / "detailed-case-results.csv"
    if rows:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    return path


def build(*, root: Path, output: Path, force: bool) -> dict[str, Any]:
    cohort_path = root / COHORT.name
    cohort = [
        line.strip()
        for line in cohort_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    method_dirs = {method: root / path.name for method, path in METHODS.items()}
    method_records = {method: _load_records(path) for method, path in method_dirs.items()}
    incomplete = {
        method: sorted(set(cohort) - set(method_records[method]))
        for method in PRIMARY_METHODS
        if set(cohort) - set(method_records[method])
    }
    if incomplete:
        raise ValueError(
            "Primary experiment results are incomplete: "
            + ", ".join(
                f"{method}={len(missing)} missing" for method, missing in incomplete.items()
            )
        )
    non_ok = {
        method: [
            matchup_id
            for matchup_id in cohort
            if str(method_records[method][matchup_id].get("status", "")).upper() != "OK"
        ]
        for method in PRIMARY_METHODS
    }
    non_ok = {method: matchup_ids for method, matchup_ids in non_ok.items() if matchup_ids}
    if non_ok:
        raise ValueError(
            "Primary experiment has non-OK results: "
            + ", ".join(f"{method}={len(matchup_ids)}" for method, matchup_ids in non_ok.items())
        )
    smoke_cohort = [
        line.strip()
        for line in SMOKE_COHORT.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    for method in ("libradtran_continental_smoke", "sixs_cci_exact_smoke"):
        smoke_records = method_records[method]
        smoke_missing = [
            matchup_id for matchup_id in smoke_cohort if matchup_id not in smoke_records
        ]
        smoke_non_ok = [
            matchup_id
            for matchup_id in smoke_cohort
            if matchup_id in smoke_records
            and str(smoke_records[matchup_id].get("status", "")).upper() != "OK"
        ]
        if smoke_missing or smoke_non_ok:
            raise ValueError(
                f"{method} outputs are incomplete: missing={smoke_missing}, non_ok={smoke_non_ok}"
            )
    missing_costs = [
        matchup_id
        for matchup_id in cohort
        if not (root / FIXED_COSTS.name / f"{matchup_id}.npz").exists()
        or len(list((root / CCI_COSTS.name).glob(f"{matchup_id}.species*.npz"))) != 3
    ]
    if missing_costs:
        raise ValueError(f"Cost artifacts are incomplete for {len(missing_costs)} matchups")
    smoke_cost_missing = {
        "sixs_cci_exact_smoke": [
            matchup_id
            for matchup_id in smoke_cohort
            if not (root / EXACT_CCI_COSTS.name / f"{matchup_id}.species0.npz").exists()
        ],
        "libradtran_continental_smoke": [
            matchup_id
            for matchup_id in smoke_cohort
            if not (root / DIRECT_COSTS.name / f"{matchup_id}.npz").exists()
        ],
    }
    smoke_cost_missing = {
        method: matchup_ids for method, matchup_ids in smoke_cost_missing.items() if matchup_ids
    }
    if smoke_cost_missing:
        raise ValueError(f"Smoke cost artifacts are incomplete: {smoke_cost_missing}")

    output.mkdir(parents=True, exist_ok=True)
    for directory in [
        output / "assets" / "spatial",
        output / "assets" / "diagnostic",
        output / "assets" / "global",
        output / "data" / "case-cache",
    ]:
        directory.mkdir(parents=True, exist_ok=True)
    shutil.copy2(WEB_ASSETS / "app.css", output / "app.css")
    shutil.copy2(WEB_ASSETS / "app.js", output / "app.js")

    base_analysis = analyze(
        root=root,
        cohort_path=cohort_path,
        output=output,
        require_complete=True,
        replay_species=False,
    )
    baseline = method_records["lut_continental"]
    truths = {matchup_id: float(baseline[matchup_id]["truth"]) for matchup_id in cohort}
    metrics = {
        method: _metric_summary(records, cohort, truths)
        for method, records in method_records.items()
    }
    smoke_truths = {matchup_id: truths[matchup_id] for matchup_id in smoke_cohort}
    for method in ("libradtran_continental_smoke", "sixs_cci_exact_smoke"):
        metrics[method] = _metric_summary(
            method_records[method],
            smoke_cohort,
            smoke_truths,
        )
    failure_analysis = json.loads(
        (root / FAILURE_ANALYSIS.relative_to(ROOT)).read_text(encoding="utf-8")
    )
    failure_groups = {row["matchup_id"]: row for row in failure_analysis["failure_cases"]}

    cases = []
    for index, matchup_id in enumerate(cohort, start=1):
        cases.append(
            _build_case(
                output=output,
                matchup_id=matchup_id,
                method_records=method_records,
                failure_groups=failure_groups,
                force=force,
            )
        )
        if index % 5 == 0 or index == len(cohort):
            print(f"case evidence {index}/{len(cohort)}")

    comparisons = dict(base_analysis["comparisons"])
    comparisons["sixs_continental_vs_lut"]["bootstrap_abs_error_delta"] = (
        _bootstrap_abs_error_delta(
            method_records["lut_continental"],
            method_records["sixs_continental"],
            cohort,
            truths,
        )
    )
    comparisons["sixs_cci3_vs_lut"]["bootstrap_abs_error_delta"] = _bootstrap_abs_error_delta(
        method_records["lut_continental"],
        method_records["sixs_cci3"],
        cohort,
        truths,
    )
    comparisons["sixs_cci3_vs_sixs_continental"]["bootstrap_abs_error_delta"] = (
        _bootstrap_abs_error_delta(
            method_records["sixs_continental"],
            method_records["sixs_cci3"],
            cohort,
            truths,
        )
    )
    smoke_ids = smoke_cohort
    direct_comparison = None
    exact_comparisons = None
    if smoke_ids:
        from tools.aeronet_validation.analyze_aerosol_species_experiment import paired_comparison

        direct_comparison = paired_comparison(
            method_records["lut_continental"],
            method_records["libradtran_continental_smoke"],
            smoke_ids,
            smoke_truths,
        )
        exact_comparisons = {
            "vs_sixs_continental": paired_comparison(
                method_records["sixs_continental"],
                method_records["sixs_cci_exact_smoke"],
                smoke_ids,
                smoke_truths,
            ),
            "vs_sixs_cci3": paired_comparison(
                method_records["sixs_cci3"],
                method_records["sixs_cci_exact_smoke"],
                smoke_ids,
                smoke_truths,
            ),
        }

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "title": "Low-cloud AOD aerosol-species experiment",
        "cohort_count": len(cohort),
        "cohort_definition": "152 Sentinel-2/AERONET matchups with image cloud fraction below 20%",
        "expected_error": "abs(retrieved - AERONET) <= 0.05 + 0.15 x AERONET",
        "method_labels": METHOD_LABELS,
        "method_short_labels": METHOD_SHORT,
        "method_colors": METHOD_COLORS,
        "metrics": metrics,
        "comparisons": comparisons,
        "direct_libradtran_smoke": {
            "case_count": len(smoke_ids),
            "matchup_ids": smoke_ids,
            "comparison_vs_lut": direct_comparison,
        },
        "exact_cci_smoke": {
            "case_count": len(smoke_ids),
            "matchup_ids": smoke_ids,
            "comparisons": exact_comparisons,
        },
        "stratification": _stratification(cases),
        "final_pass_policy_replays": _final_pass_policy_replays(cases),
        "cases": cases,
        "method_contract": [
            {
                "method": "lut_continental",
                "rt": "Precomputed libRadTran continental_average Zarr LUT",
                "species": "One fixed continental-average optical model",
                "scope": "Full 152-case control",
            },
            {
                "method": "sixs_continental",
                "rt": "Native 6S direct coefficient generation",
                "species": "One fixed 6S continental profile",
                "scope": "Full 152-case backend control",
            },
            {
                "method": "sixs_cci3",
                "rt": "Native 6S direct coefficient generation",
                "species": "Three nearest monthly 1-degree Aerosol_cci mixtures; per-pixel minimum pooled surface-only cost",
                "scope": "Full 152-case species arm",
            },
            {
                "method": "sixs_cci_exact_smoke",
                "rt": "Native 6S direct coefficient generation",
                "species": "One continuous monthly 1-degree Aerosol_cci mixture for the scene centre; no candidate selector",
                "scope": f"{len(smoke_ids)}-case simple-species smoke arm",
            },
            {
                "method": "libradtran_continental_smoke",
                "rt": "Direct libRadTran uvspec; 16-stream DISORT with reptran medium windows and fine deep-water segments",
                "species": "One fixed continental-average profile",
                "scope": f"{len(smoke_ids)}-case implementation smoke control",
            },
        ],
        "fixed_recipe": {
            "surface_prior": "One S2 monthly NIR/SWIR-anchored ExtraTrees surface prior",
            "surface_predictors": "B8A, B11, B12",
            "predicted_and_solved_bands": "B02, B03, B04",
            "cost": "Three-band chi-square; visible uncertainty floor 0.006",
            "surface_model": "ExtraTrees, 20 estimators",
            "iteration": "One anchor iteration; pass 1 changes the final surface prediction",
            "atmospheric_backstop": "CAMS/MAIAC Gaussian AOD backstop, tau gate 8, calibrated backstop",
            "masking": "Cloud/water masks by default; the canonical recorded broad-mask fallback is used only when masked support is zero",
            "spatial_output": "Scene mean of valid 60 m solver pixels",
            "evaluation": "AERONET is used only after retrieval for scoring and visualization",
        },
        "source_urls": {
            "case_csv": "case-results.csv",
            "detailed_case_csv": "data/detailed-case-results.csv",
            "analysis_json": "analysis.json",
            "policy_replay_analysis": "../aod-aerosol-species-policy-replay-20260712/analysis.json",
            "baseline_failure_explorer": "../aod-low-cloud-failure-explorer-20260712/",
            **{method: f"../../{path.name}/" for method, path in method_dirs.items()},
            "fixed_costs": f"../../{FIXED_COSTS.name}/",
            "cci_costs": f"../../{CCI_COSTS.name}/",
            "exact_cci_costs": f"../../{EXACT_CCI_COSTS.name}/",
            "direct_libradtran_costs": f"../../{DIRECT_COSTS.name}/",
        },
        "audit_notes": [
            "The fixed-6S versus CCI-3 pair holds the RT backend, bands, cost, masks, pooling, backstop settings, surface-prior type, and extraction fixed.",
            "Tai Ping, MCO Hanimaadhoo, and Churchill have zero support under the normal mask and use the same recorded cloud/water-mask fallback as the LUT baseline in every primary arm.",
            "Because anchor iteration is enabled, aerosol optics affect pass 1 and therefore the final predicted visible surface reflectance; the species comparison is an end-to-end operational effect, not a frozen-surface sensitivity.",
            "The CCI arm replaces the fixed continental candidate with the three nearest CCI mixtures; fixed continental is not a fourth candidate in this arm.",
            "The exact-CCI smoke arm uses one continuous location/month climatological mixture and has no per-pixel or scene-level species competition.",
            "Candidate selection minimizes pooled observation-only surface cost per pixel. The CAMS/MAIAC penalty affects each candidate's AOD node but not the candidate-rank comparison returned by the kernel.",
            "Candidate ranks are nearest-mixture ranks and their physical component fractions vary by scene and month.",
            "Nearest-mixture and scene-wide policy replays use the CCI arm's saved final-pass cubes. They are diagnostics of the selection policy, not complete reruns with that policy applied in pass 1.",
            "The direct-libRadTran smoke arm uses the runtime adaptive medium/fine molecular-absorption setup; it is an implementation control and not a claim of bit-identical coefficients to the precomputed Zarr LUT.",
            "AERONET does not participate in species selection, iteration, or retrieval; it is the held-out scoring reference shown here.",
            "The page reports gains and losses symmetrically and does not make an operational method decision.",
        ],
    }
    _save_global_figures(output, report)
    _write_detailed_case_csv(output, cases)
    (output / "data" / "report.json").write_text(
        json.dumps(_jsonable(report), separators=(",", ":")), encoding="utf-8"
    )
    _write_index(output)
    receipt = {
        "output": str(output),
        "cases": len(cases),
        "spatial_images": len(list((output / "assets" / "spatial").glob("*.png"))),
        "diagnostic_images": len(list((output / "assets" / "diagnostic").glob("*.png"))),
        "global_images": len(list((output / "assets" / "global").glob("*.png"))),
        "direct_libradtran_smoke_cases": len(smoke_ids),
        "exact_cci_smoke_cases": len(smoke_ids),
    }
    (output / "build-receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    receipt = build(root=args.root, output=args.output, force=args.force)
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
