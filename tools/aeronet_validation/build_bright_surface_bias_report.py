"""Build a detailed bright-surface attribution report for the medium-AOD cohort.

This is an investigation artifact, not a retrieval recipe.  It combines saved
surface-driven cost cubes with exact same-day L2A/L1C-current-RT pairs to show
where bright-target error enters the current chain.  AERONET is used only to
score and visualise the already-saved retrievals; no AERONET value is used to
fit the pair model or select a correction.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from pyproj import Transformer

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
CUBES = ROOT / "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
FRESH = ROOT / "phaseD_results_lowcloud20_mediumphysics_baseline_costcube_mediumdev_20260713"
CURRENT = ROOT / "analysis/medium_aod_current_end_to_end_prior_conflict_z2576_development_20260713"
PAIR_DIR = ROOT / "analysis/l2a_l1c_exact_pairs_mediumdev_20260713"
DIRECT_HARMONIZED = ROOT / "phaseD_results_lowcloud20_l2aharm_full_a100_all_cap0p030_direct_static_mediumdev_20260713"
DEFAULT_OUTPUT = ROOT / "reports/bright-surface-bias-20260714"
WEB_ASSETS = Path(__file__).with_name("bright_surface_bias_report")

BANDS = ("B02", "B03", "B04")
PAIR_BANDS = ("blue", "green", "red")
BAND_COLORS = {"B02": "#2d6db4", "B03": "#2a8c67", "B04": "#c64e42"}
LUMINANCE = np.asarray([0.0722, 0.7152, 0.2126], dtype=np.float64)
PAIR_TRAINING_CUTOFF = "2023-12-31"
PAIR_CORRECTION_CAP = 0.03


def _finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _ee_threshold(truth: float) -> float:
    return 0.05 + 0.15 * truth


def _within_ee(value: float | None, truth: float) -> bool | None:
    return None if value is None else abs(value - truth) <= _ee_threshold(truth)


def _outcome(value: float | None, truth: float) -> str:
    if value is None:
        return "missing"
    if _within_ee(value, truth):
        return "hit"
    return "under" if value < truth else "over"


def _window_mask(
    x: np.ndarray,
    y: np.ndarray,
    *,
    lon: float,
    lat: float,
    crs: str,
    radius_m: float = 1500.0,
) -> np.ndarray:
    center_x, center_y = Transformer.from_crs(
        "EPSG:4326", crs, always_xy=True
    ).transform(lon, lat)
    xx, yy = np.meshgrid(x, y)
    return np.square(xx - center_x) + np.square(yy - center_y) <= radius_m**2


def _median_curve(cube: np.ndarray, mask: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore"):
        return np.nanmedian(np.asarray(cube, dtype=np.float64)[:, mask], axis=1)


def _normalise_curve(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.zeros(values.shape, dtype=np.float64)
    return np.log10(1.0 + np.maximum(values - np.min(finite), 0.0))


def _luminance(values: np.ndarray) -> np.ndarray:
    return np.einsum("b,b...->...", LUMINANCE, np.asarray(values, dtype=np.float64))


def _correlation(left: np.ndarray, right: np.ndarray) -> float | None:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    good = np.isfinite(left) & np.isfinite(right)
    if int(np.count_nonzero(good)) < 20:
        return None
    x = left[good]
    y = right[good]
    if np.std(x) <= 1.0e-12 or np.std(y) <= 1.0e-12:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def _linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float | None, float | None]:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    good = np.isfinite(x) & np.isfinite(y)
    if int(np.count_nonzero(good)) < 20 or np.std(x[good]) <= 1.0e-12:
        return None, None
    slope, intercept = np.polyfit(x[good], y[good], 1)
    return float(slope), float(intercept)


def _safe_median(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.median(finite)) if finite.size else None


def _safe_mean(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.mean(finite)) if finite.size else None


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def _metrics(cases: list[dict[str, Any]], key: str) -> dict[str, Any]:
    pairs = [
        (float(case[key]), float(case["truth"]))
        for case in cases
        if _finite(case.get(key)) is not None
    ]
    if not pairs:
        return {"n": 0, "within_ee": 0, "within_ee_percent": None, "bias": None, "mae": None}
    errors = np.asarray([value - truth for value, truth in pairs], dtype=np.float64)
    return {
        "n": len(pairs),
        "within_ee": int(sum(abs(value - truth) <= _ee_threshold(truth) for value, truth in pairs)),
        "within_ee_percent": float(
            100.0 * sum(abs(value - truth) <= _ee_threshold(truth) for value, truth in pairs) / len(pairs)
        ),
        "bias": float(np.mean(errors)),
        "mae": float(np.mean(np.abs(errors))),
        "rmse": float(np.sqrt(np.mean(np.square(errors)))),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({field for row in rows for field in row})
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(_jsonable(row) for row in rows)


def _pair_training_statistics() -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Return source-L2A -> current-RT affine sufficient statistics by site.

    Each row is ``count, sum(x), sum(x^2), sum(y), sum(x*y), sum(y^2)`` for
    ``y = current_RT_surface - L2A_surface``.  The cutoff makes this a
    production-valid historical training screen rather than future leakage.
    """

    by_site: dict[str, np.ndarray] = {}
    for path in sorted(PAIR_DIR.glob("*.npz")):
        site = path.stem.split("__", 1)[0]
        with np.load(path, allow_pickle=False) as pair:
            scenes = json.loads(str(pair["scenes_json"].item()))
            scene_keep = np.asarray(
                [str(scene.get("day", "")) <= PAIR_TRAINING_CUTOFF for scene in scenes],
                dtype=bool,
            )
            keep = scene_keep[np.asarray(pair["scene_index"], dtype=np.int64)]
            l2a = np.asarray(pair["l2a"], dtype=np.float64)[keep]
            target = np.asarray(pair["siac"], dtype=np.float64)[keep]
            names = [str(value) for value in pair["band_names"].tolist()]
        if not l2a.size:
            continue
        stats = by_site.setdefault(site, np.zeros((3, 6), dtype=np.float64))
        for index, band in enumerate(PAIR_BANDS):
            column = names.index(band)
            x = l2a[:, column]
            y = target[:, column] - x
            good = (
                np.isfinite(x)
                & np.isfinite(y)
                & (x > 0.001)
                & (x < 0.8)
                & (np.abs(y) < 0.2)
            )
            x = x[good]
            y = y[good]
            stats[index] += np.asarray(
                [x.size, np.sum(x), np.dot(x, x), np.sum(y), np.dot(x, y), np.dot(y, y)],
                dtype=np.float64,
            )
    if not by_site:
        raise RuntimeError(f"no usable exact pairs found under {PAIR_DIR}")
    return by_site, np.sum(np.stack(list(by_site.values()), axis=0), axis=0)


def _fit_held_site_pair_model(
    *,
    all_stats: np.ndarray,
    by_site: dict[str, np.ndarray],
    site: str,
) -> tuple[np.ndarray, np.ndarray]:
    stats = np.asarray(all_stats, dtype=np.float64) - by_site.get(site, 0.0)
    coefficients = np.zeros((3, 2), dtype=np.float64)
    residual_sigma = np.zeros(3, dtype=np.float64)
    for index in range(3):
        count, sx, sx2, sy, sxy, sy2 = stats[index]
        if count < 100.0:
            raise RuntimeError(f"too few held-site pair samples for {site!r}")
        matrix = np.asarray([[count, sx], [sx, sx2]], dtype=np.float64)
        rhs = np.asarray([sy, sxy], dtype=np.float64)
        intercept, slope = np.linalg.solve(matrix, rhs)
        coefficients[index] = [intercept, slope]
        sse = (
            sy2
            - 2.0 * intercept * sy
            - 2.0 * slope * sxy
            + intercept * intercept * count
            + 2.0 * intercept * slope * sx
            + slope * slope * sx2
        )
        residual_sigma[index] = math.sqrt(max(sse / count, 1.0e-10))
    return coefficients, residual_sigma


def _pair_case_metrics(path: Path) -> tuple[dict[str, Any], tuple[np.ndarray, np.ndarray] | None]:
    if not path.exists():
        return {"available": False}, None
    with np.load(path, allow_pickle=False) as pair:
        names = [str(value) for value in pair["band_names"].tolist()]
        l2a = np.asarray(pair["l2a"], dtype=np.float64)
        target = np.asarray(pair["siac"], dtype=np.float64)
    columns = [names.index(band) for band in PAIR_BANDS]
    source = l2a[:, columns].T
    current_rt = target[:, columns].T
    residual = source - current_rt
    source_luminance = _luminance(source)
    current_rt_luminance = _luminance(current_rt)
    residual_luminance = _luminance(residual)
    midpoint_luminance = 0.5 * (source_luminance + current_rt_luminance)
    good = (
        np.isfinite(source_luminance)
        & np.isfinite(current_rt_luminance)
        & np.isfinite(residual_luminance)
    )
    source_luminance = source_luminance[good]
    midpoint_luminance = midpoint_luminance[good]
    residual_luminance = residual_luminance[good]
    source_slope, source_intercept = _linear_fit(source_luminance, residual_luminance)
    midpoint_slope, midpoint_intercept = _linear_fit(midpoint_luminance, residual_luminance)
    if midpoint_luminance.size:
        midpoint_low = _safe_median(
            residual_luminance[midpoint_luminance <= np.quantile(midpoint_luminance, 0.25)]
        )
        midpoint_high = _safe_median(
            residual_luminance[midpoint_luminance >= np.quantile(midpoint_luminance, 0.75)]
        )
    else:
        midpoint_low = midpoint_high = None
    sample = None
    if source_luminance.size:
        count = min(6000, source_luminance.size)
        indices = np.linspace(0, source_luminance.size - 1, count, dtype=np.int64)
        sample = (midpoint_luminance[indices], residual_luminance[indices])
    return (
        {
            "available": True,
            "samples": int(source_luminance.size),
            "l2a_minus_current_rt_slope": source_slope,
            "l2a_minus_current_rt_intercept": source_intercept,
            "l2a_minus_current_rt_correlation": _correlation(source_luminance, residual_luminance),
            "l2a_minus_current_rt_vs_mean_slope": midpoint_slope,
            "l2a_minus_current_rt_vs_mean_intercept": midpoint_intercept,
            "l2a_minus_current_rt_vs_mean_correlation": _correlation(
                midpoint_luminance, residual_luminance
            ),
            "l2a_minus_current_rt_median": _safe_median(residual_luminance),
            "l2a_minus_current_rt_q1": midpoint_low,
            "l2a_minus_current_rt_q4": midpoint_high,
            "l2a_minus_current_rt_q4_minus_q1": None
            if midpoint_low is None or midpoint_high is None
            else midpoint_high - midpoint_low,
        },
        sample,
    )


def _screen_solution(
    cost: np.ndarray,
    axis: np.ndarray,
    aot_prior: np.ndarray,
    aot_unc: np.ndarray,
    valid: np.ndarray,
) -> float | None:
    """Cheap independent-pixel screen, intentionally separate from production pooling."""

    with np.errstate(invalid="ignore", divide="ignore"):
        total = cost + np.square(
            (axis[:, np.newaxis, np.newaxis] - aot_prior[np.newaxis])
            / np.maximum(aot_unc[np.newaxis], 1.0e-6)
        )
    index = np.argmin(np.where(valid[np.newaxis] & np.isfinite(total), total, np.inf), axis=0)
    output = np.asarray(axis, dtype=np.float64)[index]
    output[~valid] = np.nan
    return _safe_mean(output)


def _rgb(values: np.ndarray) -> np.ndarray:
    rgb = np.stack([values[2], values[1], values[0]], axis=-1).astype(np.float64)
    return np.clip(rgb / 0.35, 0.0, 1.0) ** (1.0 / 1.25)


def _case_figure(path: Path, case: dict[str, Any], data: dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    axis = data["axis"]
    truth = float(case["truth"])
    fig, axes = plt.subplots(2, 4, figsize=(19, 10), constrained_layout=True, facecolor="white")
    flat = axes.ravel()
    for index, (item, title) in enumerate(zip(
        (data["toa"], data["prior"], data["truth_boa"]),
        ("Target TOA", "Operating BOA prior", "BOA at AERONET node"),
        strict=True,
    )):
        ax = flat[index]
        ax.imshow(_rgb(item), interpolation="nearest")
        ax.set(title=f"{title}\nB04 / B03 / B02")
        ax.set_xticks([])
        ax.set_yticks([])
    residual_map = data["signed_luminance"].copy()
    residual_map[~data["local_mask"]] = np.nan
    finite = residual_map[np.isfinite(residual_map)]
    limit = max(float(np.percentile(np.abs(finite), 98)) if finite.size else 0.01, 0.01)
    image = flat[3].imshow(residual_map, cmap="RdBu_r", vmin=-limit, vmax=limit, interpolation="nearest")
    flat[3].set(title="Diagnostic BOA minus prior\nwithin 1.5 km station window")
    flat[3].set_xticks([])
    flat[3].set_yticks([])
    fig.colorbar(image, ax=flat[3], fraction=0.045, pad=0.02, label="Reflectance")

    total_curve = data["curve"]
    flat[4].plot(axis, _normalise_curve(total_curve), color="#111827", linewidth=2.0, label="Joint surface")
    for name, curve in zip(BANDS, data["band_curves"], strict=True):
        flat[4].plot(axis, _normalise_curve(curve), color=BAND_COLORS[name], linewidth=1.25, label=name)
    for value, color, label, style in (
        (truth, "#111111", "AERONET", "-"),
        (case["retrieved"], "#2364aa", "Current", "--"),
        (case["maiac_aod"], "#b34d2e", "MAIAC", ":"),
        (case["local_surface_min"], "#6855a5", "Surface min", "-."),
    ):
        if value is not None:
            flat[4].axvline(value, color=color, linestyle=style, linewidth=1.1, label=label)
    flat[4].set(title="Station-window surface cost", xlabel="AOD550", ylabel="log10(1 + delta cost)")
    flat[4].legend(fontsize=7, ncol=2)

    for name, residual in zip(BANDS, data["signed_curves"], strict=True):
        flat[5].plot(axis, residual, color=BAND_COLORS[name], linewidth=1.4, label=name)
    flat[5].axhline(0.0, color="#6b7280", linewidth=0.8)
    flat[5].axvline(truth, color="#111827", linewidth=1.1)
    flat[5].set(title="Median signed BOA residual", xlabel="AOD550", ylabel="BOA - prior")
    flat[5].legend(fontsize=8)

    pair_sample = data.get("pair_sample")
    if pair_sample is not None:
        pair_x, pair_y = pair_sample
        flat[6].scatter(pair_x, pair_y, s=3, alpha=0.18, color="#4f7cac", rasterized=True)
        slope = case["pair"].get("l2a_minus_current_rt_vs_mean_slope")
        intercept = case["pair"].get("l2a_minus_current_rt_vs_mean_intercept")
        if slope is not None and intercept is not None:
            xline = np.linspace(np.nanpercentile(pair_x, 1), np.nanpercentile(pair_x, 99), 100)
            flat[6].plot(xline, intercept + slope * xline, color="#b34d2e", linewidth=1.8)
        flat[6].axhline(0.0, color="#6b7280", linewidth=0.8)
        flat[6].set(
            title="Exact same-day pair diagnostic",
            xlabel="Mean L2A /\ncurrent-RT visible",
            ylabel="L2A - current RT",
        )
    else:
        flat[6].text(0.5, 0.5, "No exact-pair archive", ha="center", va="center")
        flat[6].set_axis_off()

    target_x = data["target_midpoint_luminance"]
    target_y = data["signed_luminance"][data["local_mask"]]
    finite = np.isfinite(target_x) & np.isfinite(target_y)
    target_x = target_x[finite]
    target_y = target_y[finite]
    if target_x.size:
        count = min(5000, target_x.size)
        take = np.linspace(0, target_x.size - 1, count, dtype=np.int64)
        flat[7].scatter(target_x[take], target_y[take], s=3, alpha=0.2, color="#3f7f5f", rasterized=True)
        slope = case["target_residual_vs_midpoint_slope"]
        intercept = case["target_residual_vs_midpoint_intercept"]
        if slope is not None and intercept is not None:
            xline = np.linspace(np.nanpercentile(target_x, 1), np.nanpercentile(target_x, 99), 100)
            flat[7].plot(xline, intercept + slope * xline, color="#b34d2e", linewidth=1.8)
    flat[7].axhline(0.0, color="#6b7280", linewidth=0.8)
    flat[7].set(
        title="Target residual diagnostic",
        xlabel="Mean operating-prior /\nRT BOA",
        ylabel="BOA - prior",
    )
    fig.suptitle(
        f"{case['site']} | current {case['retrieved']:.3f}, AERONET {truth:.3f}, {case['outcome']}",
        fontsize=14,
        fontweight="bold",
    )
    fig.savefig(path, dpi=155)
    plt.close(fig)


def _summary_figures(output: Path, cases: list[dict[str, Any]]) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output.mkdir(parents=True, exist_ok=True)
    brightness = np.asarray([case["brightness"] for case in cases], dtype=np.float64)
    information = np.asarray([case["aod_information"] for case in cases], dtype=np.float64)
    spread = np.asarray([case["band_minimum_spread"] for case in cases], dtype=np.float64)
    colors = {"hit": "#2f7d5e", "under": "#b44b3f", "over": "#365f9d"}
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.1), constrained_layout=True, facecolor="white")
    for case in cases:
        axes[0].scatter(case["brightness"], case["error"], color=colors[case["outcome"]], s=48, alpha=0.85)
        if case["outcome"] != "hit":
            axes[0].annotate(case["site"], (case["brightness"], case["error"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
    axes[0].axhline(0.0, color="#6b7280", linewidth=0.8)
    axes[0].set(title="Current retrieval error versus operating-prior brightness", xlabel="Visible prior luminance", ylabel="Retrieved - AERONET AOD")
    axes[1].scatter(brightness, information, color="#3f7f5f", s=48, alpha=0.8, label="AOD information")
    twin = axes[1].twinx()
    twin.scatter(brightness, spread, color="#b34d2e", s=40, marker="s", alpha=0.72, label="Band-minimum spread")
    axes[1].set(title="Bright targets have weaker and less coherent surface evidence", xlabel="Visible prior luminance", ylabel="AOD information (1/AOD)")
    twin.set_ylabel("Local per-band AOD-minimum spread")
    fig.savefig(output / "brightness-error-and-information.png", dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.1), constrained_layout=True, facecolor="white")
    pair_slope = np.asarray(
        [case["pair"].get("l2a_minus_current_rt_vs_mean_slope") for case in cases],
        dtype=np.float64,
    )
    pair_delta = np.asarray([case["pair"].get("l2a_minus_current_rt_q4_minus_q1") for case in cases], dtype=np.float64)
    axes[0].scatter(brightness, pair_slope, color="#6f5aa6", s=48, alpha=0.82)
    axes[0].axhline(0.0, color="#6b7280", linewidth=0.8)
    axes[0].set(
        title="Exact-pair difference-versus-midpoint slope",
        xlabel="Operating prior luminance",
        ylabel="L2A - current RT slope",
    )
    axes[1].scatter(brightness, pair_delta, color="#b34d2e", s=48, alpha=0.82)
    axes[1].axhline(0.0, color="#6b7280", linewidth=0.8)
    axes[1].set(
        title="Exact-pair bright-minus-dark residual",
        xlabel="Operating prior luminance",
        ylabel="Q4 - Q1 L2A - current RT",
    )
    fig.savefig(output / "exact-pair-brightness-bias.png", dpi=170)
    plt.close(fig)

    groups = ["low", "mid", "high"]
    metrics = [_metrics([case for case in cases if case["brightness_group"] == group], "retrieved") for group in groups]
    screen = [_metrics([case for case in cases if case["brightness_group"] == group], "screen_mapped_pair_calibrated") for group in groups]
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.6), constrained_layout=True, facecolor="white")
    x = np.arange(len(groups))
    axes[0].bar(x - 0.18, [metric["within_ee_percent"] for metric in metrics], width=0.36, color="#2f7d5e", label="Current pooled retrieval")
    axes[0].bar(x + 0.18, [metric["within_ee_percent"] for metric in screen], width=0.36, color="#aa6a45", label="Raw-pair post-prior screen")
    axes[0].set_xticks(x, ["Low", "Middle", "High"])
    axes[0].set_ylim(0, 105)
    axes[0].set(title="Raw L2A-pair correction screen", ylabel="Within EE (%)")
    axes[0].legend(fontsize=8)
    axes[1].bar(x - 0.18, [metric["bias"] for metric in metrics], width=0.36, color="#2f7d5e", label="Current pooled retrieval")
    axes[1].bar(x + 0.18, [metric["bias"] for metric in screen], width=0.36, color="#aa6a45", label="Raw-pair post-prior screen")
    axes[1].axhline(0.0, color="#6b7280", linewidth=0.8)
    axes[1].set_xticks(x, ["Low", "Middle", "High"])
    axes[1].set(title="The direct correction shifts high-brightness cases lower", ylabel="AOD bias")
    axes[1].legend(fontsize=8)
    fig.savefig(output / "post-prior-screen.png", dpi=170)
    plt.close(fig)
    return {
        "brightness_error": "assets/brightness-error-and-information.png",
        "pair_bias": "assets/exact-pair-brightness-bias.png",
        "screen": "assets/post-prior-screen.png",
    }


def _build_case(
    cube_path: Path,
    *,
    fresh: dict[str, Any],
    current: dict[str, Any],
    direct: dict[str, Any] | None,
    pair_coefficients: np.ndarray,
    pair_sigma: np.ndarray,
    output: Path,
    reuse_figures: bool,
) -> dict[str, Any]:
    matchup_id = cube_path.stem
    truth = float(current["truth"])
    pair, pair_sample = _pair_case_metrics(PAIR_DIR / cube_path.name)
    with np.load(cube_path, allow_pickle=False) as archive:
        axis = np.asarray(archive["aot_axis"], dtype=np.float64)
        valid = np.asarray(archive["solve_valid"], dtype=bool)
        prior = np.asarray(archive["boa_prior"], dtype=np.float64)
        uncertainty = np.asarray(archive["boa_unc"], dtype=np.float64)
        toa = np.asarray(archive["toa"], dtype=np.float64)
        signed_cube = np.asarray(archive["band_signed_residual_cube"], dtype=np.float64)
        joint_cube = np.asarray(archive["cube"], dtype=np.float64)
        band_cost_cube = np.asarray(archive["band_cost_cube"], dtype=np.float64)
        aot_prior = np.asarray(archive["aot_prior"], dtype=np.float64)
        aot_unc = np.asarray(archive["aot_prior_unc"], dtype=np.float64)
        x = np.asarray(archive["x"], dtype=np.float64)
        y = np.asarray(archive["y"], dtype=np.float64)
    station_mask = _window_mask(
        x,
        y,
        lon=float(current["lon"]),
        lat=float(current["lat"]),
        crs=str(current["scene_crs"]),
    )
    local_mask = station_mask & valid
    if int(np.count_nonzero(local_mask)) < 20:
        local_mask = valid
    truth_index = int(np.argmin(np.abs(axis - truth)))
    left_index = max(0, truth_index - 1)
    right_index = min(axis.size - 1, truth_index + 1)
    signed = signed_cube[:, truth_index]
    signed_luminance = _luminance(signed)
    local_prior = prior[:, local_mask]
    local_uncertainty = uncertainty[:, local_mask]
    local_signed = signed[:, local_mask]
    local_prior_luminance = _luminance(local_prior)
    local_signed_luminance = _luminance(local_signed)
    local_truth_luminance = _luminance(local_prior + local_signed)
    local_midpoint_luminance = 0.5 * (local_prior_luminance + local_truth_luminance)
    common = local_signed_luminance
    spectral = np.sqrt(np.nanmean(np.square(local_signed - common[np.newaxis]), axis=0))
    derivative = (
        signed_cube[:, right_index, local_mask] - signed_cube[:, left_index, local_mask]
    ) / max(float(axis[right_index] - axis[left_index]), 1.0e-6)
    information = np.sqrt(np.nansum(np.square(derivative / np.maximum(local_uncertainty, 1.0e-6)), axis=0))
    curve = _median_curve(joint_cube, local_mask)
    band_curves = [_median_curve(band_cost_cube[index], local_mask) for index in range(3)]
    band_minima = [float(axis[int(np.nanargmin(values))]) for values in band_curves]
    local_surface_min = float(axis[int(np.nanargmin(curve))])
    curve_truth_penalty = float(curve[truth_index] - np.nanmin(curve))
    target_prior_slope, target_prior_intercept = _linear_fit(
        local_prior_luminance, local_signed_luminance
    )
    target_midpoint_slope, target_midpoint_intercept = _linear_fit(
        local_midpoint_luminance, local_signed_luminance
    )
    # The pair correction is intentionally applied only in an independent-pixel
    # screen.  This checks directionality without claiming operational pooling.
    correction = np.stack(
        [
            np.clip(pair_coefficients[index, 0] + pair_coefficients[index, 1] * prior[index], -PAIR_CORRECTION_CAP, PAIR_CORRECTION_CAP)
            for index in range(3)
        ],
        axis=0,
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        mapped_cost = np.sum(
            np.square(signed_cube - correction[:, np.newaxis])
            / np.square(np.maximum(uncertainty[:, np.newaxis], 1.0e-6)),
            axis=0,
        )
        mapped_pair_cost = np.sum(
            np.square(signed_cube - correction[:, np.newaxis])
            / np.square(
                np.maximum(
                    np.sqrt(
                        np.square(uncertainty[:, np.newaxis])
                        + np.square(pair_sigma[:, np.newaxis, np.newaxis, np.newaxis])
                    ),
                    1.0e-6,
                )
            ),
            axis=0,
        )
    retrieved = _finite(current.get("retrieved"))
    fresh_retrieved = _finite(fresh.get("retrieved"))
    direct_retrieved = _finite((direct or {}).get("retrieved"))
    current_error = None if retrieved is None else retrieved - truth
    conflict = current.get("prior_conflict_replay") or {}
    solver = current.get("solver") or {}
    atmospheric = current.get("atmo_prior") or {}
    case = {
        "matchup_id": matchup_id,
        "site": str(current.get("site", matchup_id.split("__", 1)[0])),
        "truth": truth,
        "retrieved": retrieved,
        "retrieved_station": _finite(current.get("retrieved_station")),
        "retrieved_winmed": _finite(current.get("retrieved_winmed")),
        "fresh_retrieved": fresh_retrieved,
        "direct_harmonized_retrieved": direct_retrieved,
        "error": current_error,
        "ee": _ee_threshold(truth),
        "outcome": _outcome(retrieved, truth),
        "brightness": _safe_median(local_prior_luminance),
        "brightness_p90": float(np.nanpercentile(local_prior_luminance, 90)),
        "target_brightness_at_truth": _safe_median(_luminance(local_prior + local_signed)),
        "prior_uncertainty_absolute": _safe_median(np.nanmean(local_uncertainty, axis=0)),
        "prior_uncertainty_relative": _safe_median(np.nanmean(local_uncertainty / np.maximum(local_prior, 1.0e-4), axis=0)),
        "truth_common_residual": _safe_median(common),
        "truth_common_residual_z": _safe_median(
            _luminance(local_signed / np.maximum(local_uncertainty, 1.0e-6))
        ),
        "truth_spectral_residual": _safe_median(spectral),
        "truth_spatial_residual_std": float(np.nanstd(common)),
        "aod_information": _safe_median(information),
        "aod_resolution_proxy": _safe_median(1.0 / np.maximum(information, 1.0e-6)),
        "band_minima": dict(zip(BANDS, band_minima, strict=True)),
        "band_minimum_spread": float(np.ptp(band_minima)),
        "local_surface_min": local_surface_min,
        "surface_truth_penalty": curve_truth_penalty,
        "maiac_aod": _finite(atmospheric.get("aot_median")),
        "maiac_uncertainty": _finite(atmospheric.get("aot_unc_median")),
        "cloud_fraction": _finite(current.get("cloud_frac")),
        "invalid_fraction": _finite(current.get("invalid_frac")),
        "valid_aot_fraction": _finite(current.get("valid_aot_fraction")),
        "valid_support_fraction": _finite(solver.get("surface_valid_observation_fraction")),
        "aod_quality_score": _finite(current.get("aod_quality_score")),
        "atmospheric_context": {
            "tcwv_median": _finite(atmospheric.get("tcwv_median")),
            "elevation_median": _finite(atmospheric.get("elevation_median")),
            "tco3_median": _finite(atmospheric.get("tco3_median")),
        },
        "saved_solver_band_minima": {
            band: _finite(solver.get(f"surface_band_{band}_argmin_aot")) for band in BANDS
        },
        "saved_solver_band_final_residual": {
            band: _finite(solver.get(f"surface_band_{band}_residual_final_node"))
            for band in BANDS
        },
        "saved_prior_boa": {
            band: _finite((current.get("prior_boa") or {}).get(band)) for band in BANDS
        },
        "saved_prior_uncertainty": {
            band: _finite((current.get("prior_unc") or {}).get(band)) for band in BANDS
        },
        "surface_cost_per_band": _finite(solver.get("surface_static_cost_per_band")),
        "surface_curve_curvature": _finite(solver.get("surface_cost_curve_curvature")),
        "target_residual_vs_brightness_slope": target_prior_slope,
        "target_residual_vs_brightness_intercept": target_prior_intercept,
        "target_residual_vs_brightness_correlation": _correlation(
            local_prior_luminance, local_signed_luminance
        ),
        "target_residual_vs_midpoint_slope": target_midpoint_slope,
        "target_residual_vs_midpoint_intercept": target_midpoint_intercept,
        "target_residual_vs_midpoint_correlation": _correlation(
            local_midpoint_luminance, local_signed_luminance
        ),
        "pair": pair,
        "pair_model": {
            "training_cutoff": PAIR_TRAINING_CUTOFF,
            "coefficients_target_minus_l2a": {
                band: {"intercept": float(pair_coefficients[index, 0]), "slope": float(pair_coefficients[index, 1]), "rmse": float(pair_sigma[index])}
                for index, band in enumerate(BANDS)
            },
            "local_median_correction": {
                band: float(np.nanmedian(correction[index][local_mask]))
                for index, band in enumerate(BANDS)
            },
            "correction_cap": PAIR_CORRECTION_CAP,
        },
        "screen_baseline_calibrated": _screen_solution(joint_cube, axis, aot_prior, aot_unc, valid),
        "screen_mapped_calibrated": _screen_solution(mapped_cost, axis, aot_prior, aot_unc, valid),
        "screen_mapped_pair_calibrated": _screen_solution(mapped_pair_cost, axis, aot_prior, aot_unc, valid),
        "screen_mapped_pair_product": _screen_solution(
            mapped_pair_cost, axis, aot_prior, np.maximum(0.5 * aot_prior, 0.02), valid
        ),
        "prior_conflict_decision": conflict.get("decision"),
        "prior_conflict_z": _finite(conflict.get("standardized_positive_conflict")),
        "urls": {
            "current_result": f"../../analysis/{CURRENT.name}/{matchup_id}.json",
            "fresh_result": f"../../{FRESH.name}/{matchup_id}.json",
            "cost_cube": f"../../{CUBES.name}/{matchup_id}.npz",
            "exact_pair": f"../../analysis/{PAIR_DIR.name}/{matchup_id}.npz",
            "direct_harmonized_result": f"../../{DIRECT_HARMONIZED.name}/{matchup_id}.json" if direct is not None else None,
        },
    }
    figure_path = output / "assets/cases" / f"{matchup_id}.png"
    if not (reuse_figures and figure_path.exists()):
        _case_figure(
            figure_path,
            case,
            {
                "axis": axis,
                "toa": toa,
                "prior": prior,
                "truth_boa": prior + signed,
                "signed_luminance": signed_luminance,
                "local_mask": local_mask,
                "curve": curve,
                "band_curves": band_curves,
                "signed_curves": [_median_curve(signed_cube[index], local_mask) for index in range(3)],
                "pair_sample": pair_sample,
                "target_midpoint_luminance": local_midpoint_luminance,
            },
        )
    case["figure"] = f"assets/cases/{matchup_id}.png"
    return case


def build(output: Path, *, reuse_figures: bool) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=True)
    (output / "assets").mkdir(exist_ok=True)
    (output / "data").mkdir(exist_ok=True)
    (output / "downloads").mkdir(exist_ok=True)
    for name in ("app.css", "app.js"):
        shutil.copy2(WEB_ASSETS / name, output / name)
    by_site, all_stats = _pair_training_statistics()
    cases: list[dict[str, Any]] = []
    for cube_path in sorted(CUBES.glob("*.npz")):
        matchup_id = cube_path.stem
        fresh_path = FRESH / f"{matchup_id}.json"
        current_path = CURRENT / f"{matchup_id}.json"
        if not fresh_path.exists() or not current_path.exists():
            continue
        fresh = json.loads(fresh_path.read_text(encoding="utf-8"))
        current = json.loads(current_path.read_text(encoding="utf-8"))
        direct_path = DIRECT_HARMONIZED / f"{matchup_id}.json"
        direct = json.loads(direct_path.read_text(encoding="utf-8")) if direct_path.exists() else None
        coefficients, sigma = _fit_held_site_pair_model(
            all_stats=all_stats,
            by_site=by_site,
            site=matchup_id.split("__", 1)[0],
        )
        cases.append(
            _build_case(
                cube_path,
                fresh=fresh,
                current=current,
                direct=direct,
                pair_coefficients=coefficients,
                pair_sigma=sigma,
                output=output,
                reuse_figures=reuse_figures,
            )
        )
    if len(cases) != 36:
        raise RuntimeError(f"expected 36 complete bright-surface cases, found {len(cases)}")
    cuts = np.quantile(np.asarray([case["brightness"] for case in cases], dtype=np.float64), [1.0 / 3.0, 2.0 / 3.0])
    for case in cases:
        brightness = float(case["brightness"])
        case["brightness_group"] = "low" if brightness <= cuts[0] else "mid" if brightness <= cuts[1] else "high"
    cases.sort(key=lambda case: (case["brightness_group"], -float(case["brightness"])))
    summary_figures = _summary_figures(output / "assets", cases)
    groups = {
        group: _metrics([case for case in cases if case["brightness_group"] == group], "retrieved")
        for group in ("low", "mid", "high")
    }
    extraction_metrics = {
        "scene_mean": _metrics(cases, "retrieved"),
        "station": _metrics(cases, "retrieved_station"),
        "station_window_median": _metrics(cases, "retrieved_winmed"),
    }
    correlations = {
        key: _correlation(
            np.asarray([case[key] for case in cases], dtype=np.float64),
            np.asarray([case["error"] for case in cases], dtype=np.float64),
        )
        for key in (
            "brightness",
            "truth_common_residual",
            "truth_common_residual_z",
            "aod_information",
            "band_minimum_spread",
            "surface_truth_penalty",
        )
    }
    screen_metrics = {
        "baseline_independent_pixel": _metrics(cases, "screen_baseline_calibrated"),
        "raw_pair_post_prior": _metrics(cases, "screen_mapped_pair_calibrated"),
        "raw_pair_post_prior_product": _metrics(cases, "screen_mapped_pair_product"),
    }
    report = {
        "schema_version": 1,
        "title": "Bright-Surface AOD Attribution Explorer",
        "subtitle": "Saved 36-case medium-AOD investigation: spatial, spectral, cost and exact L2A/L1C-current-RT evidence.",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "scope": {
            "cases": len(cases),
            "brightness_definition": "Station-window median of the operating B02/B03/B04 BOA-prior luminance, weights 0.0722/0.7152/0.2126.",
            "window": "1.5 km radius around the AERONET station",
            "target": "AERONET is used only for retrospective scoring and the truth-node diagnostics.",
            "pair_target": "Exact same-day operational L2A versus L1C corrected with MAIAC and the current continental libRadtran LUT.",
            "pair_training_cutoff": PAIR_TRAINING_CUTOFF,
        },
        "summary": {
            "current": extraction_metrics["scene_mean"],
            "extraction_modes": extraction_metrics,
            "brightness_cutpoints": {"low_to_mid": float(cuts[0]), "mid_to_high": float(cuts[1])},
            "brightness_groups": groups,
            "correlations_with_current_error": correlations,
            "findings": [
                "Brightness is associated with current retrieval error: the highest operating-prior brightness tercile has materially lower within-EE coverage than the lowest tercile.",
                "The exact, AERONET-free same-day pairs retain a positive L2A-minus-current-RT trend when assessed against their mean brightness. This is an upstream frame mismatch, not merely a shared-axis correlation.",
                "The saved current retrieval does not reduce to one common bright-prior offset: bright failures include low surface-only minima, high surface-only minima, MAIAC-dominated cases, and per-band AOD disagreement.",
                "A direct post-predictor application of the raw pair affine mapping is a negative independent-pixel screen. It moves the already under-retrieved bright cases lower, so it is not promoted to a pooled operational test.",
                "The raw-history harmonizer also cannot be treated as proof against the pair effect: it alters NIR/SWIR anchors before the seasonal predictor, changing the predictor domain as well as the visible surface level.",
                f"Replacing the canonical scene mean with saved station or station-window extraction does not improve this cohort ({extraction_metrics['station']['within_ee']}/{extraction_metrics['station']['n']} and {extraction_metrics['station_window_median']['within_ee']}/{extraction_metrics['station_window_median']['n']} within EE, respectively).",
            ],
            "limitations": [
                "The truth-node BOA diagnostics are retrospective and must not be used as an operational decision variable.",
                "Difference-versus-one-input correlations are partly algebraically coupled. The figures and headline pair metrics therefore use the mean of the two surfaces on the x-axis; raw L2A/prior-axis slopes remain available only as secondary diagnostics.",
                "The post-prior mapping screen intentionally omits the production spatial pooling. It is direction-of-effect evidence, not an end-to-end retrieval result.",
                "The exact-pair bias identifies an L2A/current-RT frame mismatch, but it does not by itself establish that this is the sole cause of a given AOD miss.",
            ],
        },
        "screen": {
            "description": "Held-site affine L2A-to-current-RT pair mapping applied after the existing predicted BOA prior in an independent-pixel likelihood screen. Pair coefficients use only scenes on or before 2023-12-31; AERONET is not used by the mapping.",
            "result": screen_metrics,
            "not_promoted": True,
        },
        "figures": summary_figures,
        "cases": cases,
        "sources": {
            "current_results": f"../../analysis/{CURRENT.name}/",
            "fresh_results": f"../../{FRESH.name}/",
            "cost_cubes": f"../../{CUBES.name}/",
            "exact_pairs": f"../../analysis/{PAIR_DIR.name}/",
            "previous_harmonization_report": "../l2a-harmonization-20260714/",
            "cases_csv": "downloads/cases.csv",
            "report_json": "data/report.json",
        },
    }
    (output / "data/report.json").write_text(
        json.dumps(_jsonable(report), separators=(",", ":")) + "\n", encoding="utf-8"
    )
    _write_csv(
        output / "downloads/cases.csv",
        [
            {
                "matchup_id": case["matchup_id"],
                "site": case["site"],
                "brightness_group": case["brightness_group"],
                "brightness": case["brightness"],
                "truth": case["truth"],
                "current": case["retrieved"],
                "current_error": case["error"],
                "outcome": case["outcome"],
                "maiac_aod": case["maiac_aod"],
                "local_surface_min": case["local_surface_min"],
                "band_minimum_spread": case["band_minimum_spread"],
                "aod_information": case["aod_information"],
                "relative_surface_uncertainty": case["prior_uncertainty_relative"],
                "truth_common_residual": case["truth_common_residual"],
                "truth_spectral_residual": case["truth_spectral_residual"],
                "pair_brightness_slope_raw_l2a_axis": case["pair"].get("l2a_minus_current_rt_slope"),
                "pair_brightness_slope_mean_axis": case["pair"].get(
                    "l2a_minus_current_rt_vs_mean_slope"
                ),
                "pair_q4_minus_q1": case["pair"].get("l2a_minus_current_rt_q4_minus_q1"),
                "screen_mapped_pair_calibrated": case["screen_mapped_pair_calibrated"],
                "direct_harmonized_retrieved": case["direct_harmonized_retrieved"],
            }
            for case in cases
        ],
    )
    (output / "index.html").write_text(
        """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="light">
  <title>Bright-Surface AOD Attribution Explorer</title>
  <link rel="icon" href="data:,">
  <link rel="stylesheet" href="app.css">
</head>
<body>
  <div id="app"><div class="loading">Loading bright-surface evidence...</div></div>
  <noscript>This report requires JavaScript. The JSON and CSV downloads remain available.</noscript>
  <script src="app.js"></script>
</body>
</html>
""",
        encoding="utf-8",
    )
    receipt = {
        "output": str(output),
        "cases": len(cases),
        "case_figures": len(list((output / "assets/cases").glob("*.png"))),
        "current_hits": report["summary"]["current"]["within_ee"],
        "high_brightness_hits": groups["high"]["within_ee"],
        "low_brightness_hits": groups["low"]["within_ee"],
        "raw_pair_screen_promoted": False,
    }
    (output / "build-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--reuse-figures", action="store_true")
    args = parser.parse_args()
    print(json.dumps(build(args.output, reuse_figures=args.reuse_figures), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
