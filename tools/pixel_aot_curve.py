#!/usr/bin/env python3
"""Inspect corrected BOA as a function of AOT for selected scene pixels.

This tool rebuilds the same solver-grid inputs used by the retrieval, then
evaluates RT coefficients on a small 1x1 grid at user-selected coordinates for
multiple candidate AOT values. It is intended for diagnosing suspicious pixels
such as lower-bound AOT hits.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from tools.compare_rt_backends import (
    LegacySplitEmulator,
    _geometry_on_template,
    _load_dem_elevation,
    _reconstruct_run_atmosphere,
    _surface_prior_from_reference,
    _with_elevation,
)

from siac.adapters.rt import build_rt_model
from siac.app.planning import build_execution_plan
from siac.app.requests import SceneProcessRequest
from siac.config import SIACConfig
from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients
from siac.workflows.pipeline import _aerosol_resolution, _call_grid_assembler, _select_band_slice

LOGGER = logging.getLogger("pixel_aot_curve")
DEFAULT_AOT_VALUES = "0.001,0.05,0.1,0.2,0.4,0.8,1.2,1.6,2.0"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep AOT for selected pixels and report corrected BOA-vs-AOT "
            "curves for the LUT backend and, when available, the emulator backend."
        )
    )
    parser.add_argument("--config", type=Path, required=True, help="SIAC config TOML.")
    parser.add_argument("--input", type=Path, required=True, help="Scene SAFE path.")
    parser.add_argument(
        "--reference-output-dir",
        type=Path,
        required=True,
        help="Existing SIAC output directory containing auxiliary.nc and surface_prior.nc.",
    )
    parser.add_argument(
        "--pixel",
        nargs=2,
        type=float,
        action="append",
        metavar=("X", "Y"),
        required=True,
        help="Pixel coordinate in the solver-grid CRS. Repeat for multiple pixels.",
    )
    parser.add_argument(
        "--band",
        action="append",
        default=None,
        help="Band to inspect. Repeat for multiple bands. Defaults to all solver bands.",
    )
    parser.add_argument(
        "--aot-values",
        default=DEFAULT_AOT_VALUES,
        help=f"Comma-separated AOT sweep values. Defaults to {DEFAULT_AOT_VALUES!r}.",
    )
    parser.add_argument(
        "--emulator-dir",
        type=Path,
        default=None,
        help="Optional emulator weights directory. When omitted, emulator curves are skipped.",
    )
    parser.add_argument(
        "--surface-prior-comparison",
        type=Path,
        default=None,
        help=(
            "Optional surface_prior_comparison.nc to sample monthly/direct priors. "
            "If omitted, the tool auto-detects the common experiments/surface_prior_compare_sync location."
        ),
    )
    parser.add_argument(
        "--sensor",
        default="s2",
        help="Sensor override for config resolution. Defaults to s2.",
    )
    parser.add_argument(
        "--elevation-source",
        default="provider",
        choices=("provider", "dem"),
        help=(
            "Which elevation to feed into the RT backends. "
            "'dem' matches compare_rt_backends.py when a DEM is configured."
        ),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional JSON output path. When omitted, the JSON summary is printed to stdout.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="Optional CSV output path for flattened curve samples.",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=None,
        help="Optional directory where PNG figures will be written.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="Logging verbosity.",
    )
    return parser.parse_args()


def _setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def _parse_aot_values(raw: str) -> list[float]:
    values: list[float] = []
    for token in raw.split(","):
        stripped = token.strip()
        if not stripped:
            continue
        value = float(stripped)
        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"Invalid AOT value: {stripped!r}")
        values.append(value)
    if not values:
        raise ValueError("At least one AOT value is required.")
    return values


def _classify_curve(values: list[float], *, tol: float = 1.0e-6) -> str:
    finite = np.asarray([value for value in values if np.isfinite(value)], dtype=np.float64)
    if finite.size < 2:
        return "insufficient"
    diffs = np.diff(finite)
    if np.all(diffs >= -tol):
        return "increasing"
    if np.all(diffs <= tol):
        return "decreasing"
    return "mixed"


def _prior_position(prior: float | None, values: list[float], *, tol: float = 1.0e-6) -> str | None:
    if prior is None or not np.isfinite(prior):
        return None
    finite = np.asarray([value for value in values if np.isfinite(value)], dtype=np.float64)
    if finite.size == 0:
        return None
    lower = float(np.min(finite))
    upper = float(np.max(finite))
    if prior < lower - tol:
        return "below_curve"
    if prior > upper + tol:
        return "above_curve"
    return "within_curve"


def _nearest_curve_point(prior: float | None, points: list[dict[str, float]]) -> dict[str, float] | None:
    if prior is None or not np.isfinite(prior):
        return None
    finite_points = [point for point in points if np.isfinite(point["corrected_boa"])]
    if not finite_points:
        return None
    best = min(finite_points, key=lambda point: abs(point["corrected_boa"] - prior))
    return {
        "aot": float(best["aot"]),
        "curve_boa": float(best["corrected_boa"]),
        "abs_difference": float(abs(best["corrected_boa"] - prior)),
    }


def _scalar_da(value: float, *, x: float, y: float) -> xr.DataArray:
    return xr.DataArray(
        np.array([[np.float32(value)]], dtype=np.float32),
        dims=("y", "x"),
        coords={"y": [np.float32(y)], "x": [np.float32(x)]},
    )


def _pixelize(da: xr.DataArray, *, x: float, y: float) -> xr.DataArray:
    selected = da.sel(x=x, y=y, method="nearest")
    return _scalar_da(
        float(selected.values),
        x=float(selected.coords["x"].values),
        y=float(selected.coords["y"].values),
    )


def _sample_optional_prior(
    ds: xr.Dataset | None,
    *,
    variable: str,
    x: float,
    y: float,
) -> float | None:
    if ds is None or variable not in ds:
        return None
    value = ds[variable].sel(x=x, y=y, method="nearest").values
    scalar = float(value)
    return scalar if np.isfinite(scalar) else None


def _auto_surface_prior_comparison(reference_output_dir: Path) -> Path | None:
    candidate = reference_output_dir.parent / "experiments" / "surface_prior_compare_sync" / "surface_prior_comparison.nc"
    return candidate if candidate.exists() else None


def _build_optional_emulator_backend(
    *,
    config: SIACConfig,
    plan: Any,
    obs: Any,
    emulator_dir: Path | None,
) -> tuple[Any | None, str | None]:
    if emulator_dir is None:
        return None, "No emulator directory provided."
    if not emulator_dir.exists():
        return None, f"Emulator directory does not exist: {emulator_dir}"

    if any(emulator_dir.glob("*_xap.npz")) and any(emulator_dir.glob("*_xbp.npz")) and any(
        emulator_dir.glob("*_xcp.npz")
    ):
        return (
            LegacySplitEmulator(
                emulator_dir=emulator_dir,
                sensor_id=obs.sensor_config.sensor_id,
                satellite_id=obs.sensor_config.satellite_id,
            ),
            None,
        )

    emulator_config = config.model_copy(deep=True)
    emulator_config.rt_model.backend = "emulator"
    emulator_config.rt_model.fallback_to_lut = False
    emulator_config.paths.emulator_dir = emulator_dir
    try:
        return build_rt_model(emulator_config, auth=plan.auth, sensor_config=obs.sensor_config), None
    except Exception as exc:  # pragma: no cover - depends on local emulator installation
        return None, f"Failed to build emulator backend: {exc}"


def _collect_curve(
    *,
    rt_backend: Any,
    geometry: GeometryAngles,
    atmo_base: AtmosphericState,
    band: Any,
    obs_toa: xr.DataArray,
    aot_values: list[float],
    x: float,
    y: float,
) -> list[dict[str, float]]:
    points: list[dict[str, float]] = []
    for aot in aot_values:
        atmo_state = atmo_base.with_updated_aot_tcwv(
            aot=_scalar_da(aot, x=x, y=y),
            tcwv=atmo_base.tcwv,
        )
        coeffs: RTCoefficients = rt_backend.compute_coefficients(geometry, atmo_state, band, False)
        corrected = coeffs.apply_correction(obs_toa)
        points.append(
            {
                "aot": float(aot),
                "corrected_boa": float(corrected.values[0, 0]),
                "xap": float(coeffs.xap.values[0, 0]),
                "xbp": float(coeffs.xbp.values[0, 0]),
                "xcp": float(coeffs.xcp.values[0, 0]),
            }
        )
    return points


def _summarize_backend_curve(
    *,
    points: list[dict[str, float]],
    priors: dict[str, float | None],
) -> dict[str, Any]:
    curve_values = [point["corrected_boa"] for point in points]
    summary: dict[str, Any] = {
        "direction": _classify_curve(curve_values),
        "curve_min_boa": float(np.nanmin(curve_values)) if np.isfinite(curve_values).any() else float("nan"),
        "curve_max_boa": float(np.nanmax(curve_values)) if np.isfinite(curve_values).any() else float("nan"),
        "points": points,
        "priors": {},
    }
    for name, value in priors.items():
        summary["priors"][name] = {
            "value": value,
            "position": _prior_position(value, curve_values),
            "nearest_curve_point": _nearest_curve_point(value, points),
        }
    return summary


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _figure_basename(label: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in label.strip())
    return safe or "pixel"


def _curve_value_bounds(_pixel_result: dict[str, Any], band_result: dict[str, Any]) -> tuple[float, float]:
    values: list[float] = []
    values.extend(point["corrected_boa"] for point in band_result["lut"]["points"] if np.isfinite(point["corrected_boa"]))
    emulator = band_result.get("emulator", {})
    if emulator.get("available"):
        values.extend(point["corrected_boa"] for point in emulator["points"] if np.isfinite(point["corrected_boa"]))
    for key in ("surface_prior_input", "monthly_prior", "direct_prior"):
        value = band_result.get(key)
        if value is not None and np.isfinite(value):
            values.append(float(value))
    if not values:
        return 0.0, 1.0
    ymin = float(min(values))
    ymax = float(max(values))
    pad = 0.05 * max(1.0, abs(ymin)) if np.isclose(ymin, ymax) else 0.1 * (ymax - ymin)
    return ymin - pad, ymax + pad


def _draw_dash_line(
    draw: Any,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str,
    width: int = 2,
    dash: int = 8,
    gap: int = 5,
) -> None:
    x0, y0 = start
    x1, y1 = end
    dx = x1 - x0
    dy = y1 - y0
    length = float(np.hypot(dx, dy))
    if length <= 0.0:
        return
    step = dash + gap
    n = int(np.ceil(length / step))
    for idx in range(n):
        s0 = min(length, idx * step)
        s1 = min(length, s0 + dash)
        if s1 <= s0:
            continue
        t0 = s0 / length
        t1 = s1 / length
        draw.line(
            (
                x0 + dx * t0,
                y0 + dy * t0,
                x0 + dx * t1,
                y0 + dy * t1,
            ),
            fill=color,
            width=width,
        )


def _plot_band_panel(ax: Any, *, pixel_result: dict[str, Any], band_result: dict[str, Any]) -> None:
    lut = band_result["lut"]
    lut_x = [point["aot"] for point in lut["points"]]
    lut_y = [point["corrected_boa"] for point in lut["points"]]
    ax.plot(lut_x, lut_y, color="#1f77b4", marker="o", linewidth=2.0, label="LUT")

    emulator = band_result.get("emulator", {})
    if emulator.get("available"):
        emu_x = [point["aot"] for point in emulator["points"]]
        emu_y = [point["corrected_boa"] for point in emulator["points"]]
        ax.plot(emu_x, emu_y, color="#ff7f0e", marker="s", linewidth=1.8, linestyle="--", label="Emulator")

    prior_styles = [
        ("surface_prior_input", band_result.get("surface_prior_input"), "#2ca02c", "Surface prior"),
        ("monthly_prior", band_result.get("monthly_prior"), "#d62728", "Monthly prior"),
        ("direct_prior", band_result.get("direct_prior"), "#9467bd", "Direct prior"),
    ]
    for _, value, color, label in prior_styles:
        if value is not None and np.isfinite(value):
            ax.axhline(float(value), color=color, linewidth=1.4, linestyle=":", label=label)

    solved_aot = pixel_result["reference_run_state"]["aot"]
    if solved_aot is not None and np.isfinite(solved_aot):
        ax.axvline(float(solved_aot), color="#4d4d4d", linewidth=1.4, linestyle="-.", label="Solved AOT")

    ax.axvline(0.001, color="#9a9a9a", linewidth=1.0, linestyle="--", alpha=0.9, label="AOT floor")
    ax.set_xscale("log")
    ax.set_xlabel("AOT")
    ax.set_ylabel("Corrected BOA")
    ax.grid(True, alpha=0.25)
    ax.set_title(
        f"{band_result['band']}  LUT: {lut['direction']}"
        + (
            f"  Emulator: {emulator['direction']}"
            if emulator.get("available") and "direction" in emulator
            else ""
        )
    )

    handles, labels = ax.get_legend_handles_labels()
    dedup: dict[str, Any] = {}
    for handle, label in zip(handles, labels):
        dedup[label] = handle
    ax.legend(dedup.values(), dedup.keys(), fontsize=8, loc="best")


def _draw_band_panel_pil(
    draw: Any,
    *,
    box: tuple[int, int, int, int],
    pixel_result: dict[str, Any],
    band_result: dict[str, Any],
    font: Any,
    small_font: Any,
) -> None:
    left, top, right, bottom = box
    chart_left = left + 58
    chart_right = right - 12
    chart_top = top + 26
    chart_bottom = bottom - 38
    draw.rectangle((left, top, right, bottom), outline="#c9c9c9", width=1)
    draw.rectangle((chart_left, chart_top, chart_right, chart_bottom), outline="#4a4a4a", width=1)

    all_aot = [point["aot"] for point in band_result["lut"]["points"] if point["aot"] > 0]
    emulator = band_result.get("emulator", {})
    if emulator.get("available"):
        all_aot.extend(point["aot"] for point in emulator["points"] if point["aot"] > 0)
    min_log = float(np.log10(min(all_aot)))
    max_log = float(np.log10(max(all_aot)))
    ymin, ymax = _curve_value_bounds(pixel_result, band_result)

    def xmap(aot: float) -> float:
        if max_log == min_log:
            return float((chart_left + chart_right) / 2)
        return chart_left + (np.log10(aot) - min_log) / (max_log - min_log) * (chart_right - chart_left)

    def ymap(value: float) -> float:
        if ymax == ymin:
            return float((chart_top + chart_bottom) / 2)
        return chart_bottom - (value - ymin) / (ymax - ymin) * (chart_bottom - chart_top)

    grid_values = np.linspace(ymin, ymax, 5)
    for grid in grid_values:
        yy = ymap(float(grid))
        draw.line((chart_left, yy, chart_right, yy), fill="#ececec", width=1)
        draw.text((left + 4, yy - 6), f"{grid:.3f}", fill="#444444", font=small_font)

    for aot in sorted(set(all_aot)):
        xx = xmap(float(aot))
        draw.line((xx, chart_top, xx, chart_bottom), fill="#f3f3f3", width=1)
        draw.text((xx - 12, chart_bottom + 6), f"{aot:g}", fill="#444444", font=small_font)

    prior_styles = [
        (band_result.get("surface_prior_input"), "#2ca02c", "Surface prior"),
        (band_result.get("monthly_prior"), "#d62728", "Monthly prior"),
        (band_result.get("direct_prior"), "#9467bd", "Direct prior"),
    ]
    for value, color, label in prior_styles:
        if value is None or not np.isfinite(value):
            continue
        yy = ymap(float(value))
        _draw_dash_line(draw, (chart_left, yy), (chart_right, yy), color=color, width=2, dash=10, gap=6)
        draw.text((chart_right - 110, yy - 10), label, fill=color, font=small_font)

    solved_aot = pixel_result["reference_run_state"]["aot"]
    if solved_aot is not None and np.isfinite(solved_aot) and solved_aot > 0:
        xx = xmap(float(solved_aot))
        _draw_dash_line(draw, (xx, chart_top), (xx, chart_bottom), color="#4d4d4d", width=2, dash=12, gap=4)
        draw.text((min(xx + 4, chart_right - 70), chart_top + 4), "Solved AOT", fill="#4d4d4d", font=small_font)

    if min(all_aot) <= 0.001 <= max(all_aot):
        xx = xmap(0.001)
        _draw_dash_line(draw, (xx, chart_top), (xx, chart_bottom), color="#8c8c8c", width=1, dash=4, gap=4)

    lut_points = band_result["lut"]["points"]
    lut_xy = [(xmap(point["aot"]), ymap(point["corrected_boa"])) for point in lut_points if np.isfinite(point["corrected_boa"])]
    if len(lut_xy) >= 2:
        draw.line(lut_xy, fill="#1f77b4", width=3)
    for xx, yy in lut_xy:
        draw.ellipse((xx - 3, yy - 3, xx + 3, yy + 3), fill="#1f77b4", outline="#1f77b4")

    if emulator.get("available"):
        emu_points = emulator["points"]
        emu_xy = [(xmap(point["aot"]), ymap(point["corrected_boa"])) for point in emu_points if np.isfinite(point["corrected_boa"])]
        if len(emu_xy) >= 2:
            _draw_dash_line(draw, emu_xy[0], emu_xy[1], color="#ff7f0e", width=2, dash=5, gap=4)
            for start, end in zip(emu_xy[:-1], emu_xy[1:]):
                _draw_dash_line(draw, start, end, color="#ff7f0e", width=2, dash=8, gap=5)
        for xx, yy in emu_xy:
            draw.rectangle((xx - 3, yy - 3, xx + 3, yy + 3), fill="#ff7f0e", outline="#ff7f0e")

    title = f"{band_result['band']}  LUT: {band_result['lut']['direction']}"
    if emulator.get("available") and "direction" in emulator:
        title += f"  Emulator: {emulator['direction']}"
    draw.text((left + 8, top + 4), title, fill="#111111", font=font)
    draw.text((chart_left, chart_bottom + 22), "AOT (log scale)", fill="#222222", font=small_font)
    draw.text((left + 4, chart_top - 18), "Corrected BOA", fill="#222222", font=small_font)


def _plot_band_panel(ax: Any, *, pixel_result: dict[str, Any], band_result: dict[str, Any]) -> None:
    lut = band_result["lut"]
    lut_x = [point["aot"] for point in lut["points"]]
    lut_y = [point["corrected_boa"] for point in lut["points"]]
    ax.plot(lut_x, lut_y, color="#1f77b4", marker="o", linewidth=2.0, label="LUT")

    emulator = band_result.get("emulator", {})
    if emulator.get("available"):
        emu_x = [point["aot"] for point in emulator["points"]]
        emu_y = [point["corrected_boa"] for point in emulator["points"]]
        ax.plot(emu_x, emu_y, color="#ff7f0e", marker="s", linewidth=1.8, linestyle="--", label="Emulator")

    prior_styles = [
        ("surface_prior_input", band_result.get("surface_prior_input"), "#2ca02c", "Surface prior"),
        ("monthly_prior", band_result.get("monthly_prior"), "#d62728", "Monthly prior"),
        ("direct_prior", band_result.get("direct_prior"), "#9467bd", "Direct prior"),
    ]
    for _, value, color, label in prior_styles:
        if value is not None and np.isfinite(value):
            ax.axhline(float(value), color=color, linewidth=1.4, linestyle=":", label=label)

    solved_aot = pixel_result["reference_run_state"]["aot"]
    if solved_aot is not None and np.isfinite(solved_aot):
        ax.axvline(float(solved_aot), color="#4d4d4d", linewidth=1.4, linestyle="-.", label="Solved AOT")

    ax.axvline(0.001, color="#9a9a9a", linewidth=1.0, linestyle="--", alpha=0.9, label="AOT floor")
    ax.set_xscale("log")
    ax.set_xlabel("AOT")
    ax.set_ylabel("Corrected BOA")
    ax.grid(True, alpha=0.25)
    ax.set_title(
        f"{band_result['band']}  LUT: {lut['direction']}"
        + (
            f"  Emulator: {emulator['direction']}"
            if emulator.get("available") and "direction" in emulator
            else ""
        )
    )

    handles, labels = ax.get_legend_handles_labels()
    dedup: dict[str, Any] = {}
    for handle, label in zip(handles, labels):
        dedup[label] = handle
    ax.legend(dedup.values(), dedup.keys(), fontsize=8, loc="best")


def _plot_pixel_figure(path: Path, *, pixel_result: dict[str, Any]) -> None:
    bands = pixel_result["bands"]
    if not bands:
        return
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        from PIL import Image, ImageDraw, ImageFont

        width = 1080
        header_h = 50
        panel_h = 260
        image = Image.new("RGB", (width, header_h + panel_h * len(bands) + 20), "white")
        draw = ImageDraw.Draw(image)
        font = ImageFont.load_default()
        small_font = ImageFont.load_default()
        snapped = pixel_result["snapped"]
        req = pixel_result["requested"]
        ref = pixel_result["reference_run_state"]
        header = (
            f"{pixel_result['label']}  requested=({req['x']:.1f}, {req['y']:.1f})  "
            f"snapped=({snapped['x']:.1f}, {snapped['y']:.1f})  solved_aot={ref['aot']:.4f}"
        )
        draw.text((12, 12), header, fill="#111111", font=font)
        for index, band_result in enumerate(bands):
            top = header_h + index * panel_h
            _draw_band_panel_pil(
                draw,
                box=(10, top, width - 10, top + panel_h - 10),
                pixel_result=pixel_result,
                band_result=band_result,
                font=font,
                small_font=small_font,
            )
        image.save(path)
        return

    fig, axes = plt.subplots(
        nrows=len(bands),
        ncols=1,
        figsize=(9.0, max(3.5, 3.4 * len(bands))),
        constrained_layout=True,
    )
    if len(bands) == 1:
        axes = [axes]

    for ax, band_result in zip(axes, bands):
        _plot_band_panel(ax, pixel_result=pixel_result, band_result=band_result)

    snapped = pixel_result["snapped"]
    req = pixel_result["requested"]
    ref = pixel_result["reference_run_state"]
    fig.suptitle(
        (
            f"{pixel_result['label']}  "
            f"requested=({req['x']:.1f}, {req['y']:.1f})  "
            f"snapped=({snapped['x']:.1f}, {snapped['y']:.1f})  "
            f"solved_aot={ref['aot']:.4f}"
        ),
        fontsize=11,
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_overview_figure(path: Path, *, result: dict[str, Any]) -> None:
    pixels = result["pixels"]
    if not pixels:
        return

    max_bands = max(len(pixel["bands"]) for pixel in pixels)
    if max_bands == 0:
        return
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        from PIL import Image, ImageDraw, ImageFont

        panel_w = 520
        panel_h = 250
        header_h = 45
        width = panel_w * max_bands + 20
        height = header_h + panel_h * len(pixels) + 20
        image = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(image)
        font = ImageFont.load_default()
        small_font = ImageFont.load_default()
        draw.text((12, 12), "Pixel AOT Curve Overview", fill="#111111", font=font)
        for row, pixel_result in enumerate(pixels):
            for col in range(max_bands):
                if col >= len(pixel_result["bands"]):
                    continue
                left = 10 + col * panel_w
                top = header_h + row * panel_h
                _draw_band_panel_pil(
                    draw,
                    box=(left, top, left + panel_w - 10, top + panel_h - 10),
                    pixel_result=pixel_result,
                    band_result=pixel_result["bands"][col],
                    font=font,
                    small_font=small_font,
                )
                draw.text((left + 8, top + 6), pixel_result["label"], fill="#111111", font=small_font)
        image.save(path)
        return

    fig, axes = plt.subplots(
        nrows=len(pixels),
        ncols=max_bands,
        figsize=(5.2 * max_bands, 3.6 * len(pixels)),
        constrained_layout=True,
        squeeze=False,
    )

    for row, pixel_result in enumerate(pixels):
        for col in range(max_bands):
            ax = axes[row][col]
            if col < len(pixel_result["bands"]):
                _plot_band_panel(ax, pixel_result=pixel_result, band_result=pixel_result["bands"][col])
                ax.text(
                    0.02,
                    0.98,
                    pixel_result["label"],
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=9,
                    bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
                )
            else:
                ax.axis("off")

    fig.suptitle("Pixel AOT Curve Overview", fontsize=13)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = _parse_args()
    _setup_logging(args.log_level)

    aot_values = _parse_aot_values(args.aot_values)

    config = SIACConfig.from_file(args.config).with_overrides(sensor=args.sensor)
    request = SceneProcessRequest(config=config, input_path=args.input)
    plan = build_execution_plan(request)

    LOGGER.info("Preprocessing observation from %s", args.input)
    obs = plan.preprocessor(plan.input_path, plan.runtime_aoi)
    aerosol_resolution = _aerosol_resolution(plan.config)

    LOGGER.info("Fetching atmospheric prior")
    atmo_prior = plan.atmo_provider(obs.bounds, obs.crs, obs.metadata["observation_time"], aerosol_resolution)
    run_atmo = _reconstruct_run_atmosphere(args.reference_output_dir, atmo_prior)
    surface_prior = _surface_prior_from_reference(args.reference_output_dir / "surface_prior.nc")
    auxiliary = xr.open_dataset(args.reference_output_dir / "auxiliary.nc")

    LOGGER.info("Building LUT backend")
    lut_config = config.model_copy(deep=True)
    lut_config.rt_model.backend = "lut"
    lut_rt = build_rt_model(lut_config, auth=plan.auth, sensor_config=obs.sensor_config)

    solver_inputs = _call_grid_assembler(
        plan.grid_assembler,
        obs,
        run_atmo,
        surface_prior,
        lut_rt,
        aerosol_resolution_m=float(aerosol_resolution),
    )

    rt_atmo = solver_inputs.atmo_prior
    if args.elevation_source == "dem":
        dem_path = getattr(config.paths, "dem", None)
        if dem_path is None:
            raise ValueError("elevation_source='dem' requires config.paths.dem to be configured.")
        rt_atmo = _with_elevation(rt_atmo, _load_dem_elevation(rt_atmo.aot, dem_path))

    geometry = _geometry_on_template(solver_inputs.geometry, rt_atmo.aot)
    band_lookup = {band.name: (index, band) for index, band in enumerate(solver_inputs.bands)}

    selected_bands = args.band or list(band_lookup)
    unknown_bands = [band for band in selected_bands if band not in band_lookup]
    if unknown_bands:
        raise ValueError(f"Unknown solver band(s): {unknown_bands}")

    emulator_backend, emulator_error = _build_optional_emulator_backend(
        config=config,
        plan=plan,
        obs=obs,
        emulator_dir=args.emulator_dir,
    )

    prior_comparison_path = args.surface_prior_comparison
    if prior_comparison_path is None:
        prior_comparison_path = _auto_surface_prior_comparison(args.reference_output_dir)
    prior_comparison = xr.open_dataset(prior_comparison_path) if prior_comparison_path is not None else None

    result: dict[str, Any] = {
        "config_path": str(args.config.resolve()),
        "input_path": str(args.input.resolve()),
        "reference_output_dir": str(args.reference_output_dir.resolve()),
        "aot_values": [float(value) for value in aot_values],
        "elevation_source": args.elevation_source,
        "surface_prior_comparison": str(prior_comparison_path.resolve()) if prior_comparison_path is not None else None,
        "emulator": {
            "requested_dir": str(args.emulator_dir.resolve()) if args.emulator_dir is not None else None,
            "available": emulator_backend is not None,
            "error": emulator_error,
        },
        "pixels": [],
    }
    if args.figure_dir is not None:
        result["figures"] = {
            "directory": str(args.figure_dir.resolve()),
            "overview": None,
            "per_pixel": [],
        }
    csv_rows: list[dict[str, Any]] = []

    for pixel_index, (req_x, req_y) in enumerate(args.pixel, start=1):
        snapped = rt_atmo.aot.sel(x=req_x, y=req_y, method="nearest")
        snapped_x = float(snapped.coords["x"].values)
        snapped_y = float(snapped.coords["y"].values)

        pixel_geometry = GeometryAngles(
            sza=_pixelize(geometry.sza, x=snapped_x, y=snapped_y),
            saa=_pixelize(geometry.saa, x=snapped_x, y=snapped_y),
            vza=_pixelize(geometry.vza, x=snapped_x, y=snapped_y),
            vaa=_pixelize(geometry.vaa, x=snapped_x, y=snapped_y),
        )
        pixel_atmo = AtmosphericState(
            aot=_pixelize(rt_atmo.aot, x=snapped_x, y=snapped_y),
            tcwv=_pixelize(rt_atmo.tcwv, x=snapped_x, y=snapped_y),
            tco3=_pixelize(rt_atmo.tco3, x=snapped_x, y=snapped_y),
            aot_unc=_pixelize(rt_atmo.aot_unc, x=snapped_x, y=snapped_y),
            tcwv_unc=_pixelize(rt_atmo.tcwv_unc, x=snapped_x, y=snapped_y),
            tco3_unc=_pixelize(rt_atmo.tco3_unc, x=snapped_x, y=snapped_y),
            elevation=_pixelize(rt_atmo.elevation, x=snapped_x, y=snapped_y),
        )

        pixel_result: dict[str, Any] = {
            "label": f"pixel_{pixel_index}",
            "requested": {"x": float(req_x), "y": float(req_y)},
            "snapped": {
                "x": snapped_x,
                "y": snapped_y,
                "dx": float(snapped_x - req_x),
                "dy": float(snapped_y - req_y),
            },
            "reference_run_state": {
                "aot": float(auxiliary["aot"].sel(x=snapped_x, y=snapped_y, method="nearest").values),
                "tcwv": float(auxiliary["tcwv"].sel(x=snapped_x, y=snapped_y, method="nearest").values),
                "elevation_km": float(pixel_atmo.elevation.values[0, 0]),
            },
            "bands": [],
        }

        for band_name in selected_bands:
            band_index, band = band_lookup[band_name]
            toa_band = _select_band_slice(solver_inputs.toa, band_name=band_name, band_index=band_index)
            surface_band = _select_band_slice(
                solver_inputs.surface_prior.boa,
                band_name=band_name,
                band_index=band_index,
            )
            surface_unc_band = _select_band_slice(
                solver_inputs.surface_prior.boa_unc,
                band_name=band_name,
                band_index=band_index,
            )
            if toa_band is None or surface_band is None:
                continue

            obs_toa = _pixelize(toa_band, x=snapped_x, y=snapped_y)
            surface_prior_input = float(_pixelize(surface_band, x=snapped_x, y=snapped_y).values[0, 0])
            surface_prior_unc = (
                float(_pixelize(surface_unc_band, x=snapped_x, y=snapped_y).values[0, 0])
                if surface_unc_band is not None
                else None
            )
            monthly_prior = _sample_optional_prior(
                prior_comparison,
                variable=f"monthly_{band_name}",
                x=snapped_x,
                y=snapped_y,
            )
            direct_prior = _sample_optional_prior(
                prior_comparison,
                variable=f"direct_{band_name}",
                x=snapped_x,
                y=snapped_y,
            )

            priors = {
                "surface_prior_input": surface_prior_input,
                "monthly_prior": monthly_prior,
                "direct_prior": direct_prior,
            }

            band_result: dict[str, Any] = {
                "band": band_name,
                "observed_toa": float(obs_toa.values[0, 0]),
                "surface_prior_input": surface_prior_input,
                "surface_prior_unc": surface_prior_unc,
                "monthly_prior": monthly_prior,
                "direct_prior": direct_prior,
            }

            lut_points = _collect_curve(
                rt_backend=lut_rt,
                geometry=pixel_geometry,
                atmo_base=pixel_atmo,
                band=band,
                obs_toa=obs_toa,
                aot_values=aot_values,
                x=snapped_x,
                y=snapped_y,
            )
            band_result["lut"] = _summarize_backend_curve(points=lut_points, priors=priors)

            for point in lut_points:
                csv_rows.append(
                    {
                        "pixel_label": pixel_result["label"],
                        "requested_x": pixel_result["requested"]["x"],
                        "requested_y": pixel_result["requested"]["y"],
                        "snapped_x": snapped_x,
                        "snapped_y": snapped_y,
                        "band": band_name,
                        "backend": "lut",
                        "observed_toa": band_result["observed_toa"],
                        "surface_prior_input": surface_prior_input,
                        "monthly_prior": monthly_prior,
                        "direct_prior": direct_prior,
                        **point,
                    }
                )

            if emulator_backend is not None:
                try:
                    emulator_points = _collect_curve(
                        rt_backend=emulator_backend,
                        geometry=pixel_geometry,
                        atmo_base=pixel_atmo,
                        band=band,
                        obs_toa=obs_toa,
                        aot_values=aot_values,
                        x=snapped_x,
                        y=snapped_y,
                    )
                except Exception as exc:
                    band_result["emulator"] = {"available": False, "error": str(exc)}
                else:
                    band_result["emulator"] = {
                        "available": True,
                        **_summarize_backend_curve(points=emulator_points, priors=priors),
                    }
                    for point in emulator_points:
                        csv_rows.append(
                            {
                                "pixel_label": pixel_result["label"],
                                "requested_x": pixel_result["requested"]["x"],
                                "requested_y": pixel_result["requested"]["y"],
                                "snapped_x": snapped_x,
                                "snapped_y": snapped_y,
                                "band": band_name,
                                "backend": "emulator",
                                "observed_toa": band_result["observed_toa"],
                                "surface_prior_input": surface_prior_input,
                                "monthly_prior": monthly_prior,
                                "direct_prior": direct_prior,
                                **point,
                            }
                        )
            else:
                band_result["emulator"] = {
                    "available": False,
                    "error": emulator_error,
                }

            pixel_result["bands"].append(band_result)

        result["pixels"].append(pixel_result)

    if args.figure_dir is not None:
        args.figure_dir.mkdir(parents=True, exist_ok=True)
        overview_path = args.figure_dir / "overview.png"
        _plot_overview_figure(overview_path, result=result)
        result["figures"]["overview"] = str(overview_path.resolve())
        for pixel_result in result["pixels"]:
            figure_path = args.figure_dir / f"{_figure_basename(pixel_result['label'])}.png"
            _plot_pixel_figure(figure_path, pixel_result=pixel_result)
            pixel_result["figure_png"] = str(figure_path.resolve())
            result["figures"]["per_pixel"].append(str(figure_path.resolve()))

    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        LOGGER.info("Wrote JSON summary to %s", args.output_json)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))

    if args.output_csv is not None:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        _write_csv(args.output_csv, csv_rows)
        LOGGER.info("Wrote CSV curve samples to %s", args.output_csv)


if __name__ == "__main__":
    main()
