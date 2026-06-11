#!/usr/bin/env python3
"""Run a detailed A/B experiment for SIAC surface-prior derivation paths."""

from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr
from PIL import Image, ImageDraw, ImageFont

from siac.algorithms.solver import build_solver_valid_mask
from siac.app._assembly_surface import resolve_surface_prior_provider
from siac.app.planning import build_execution_plan
from siac.app.requests import SceneProcessRequest
from siac.config import SIACConfig
from siac.runtime import AtmosphericState, SolverInputBundle, SurfacePrior
from siac.workflows.scene_setup import aerosol_resolution

LOGGER = logging.getLogger("compare_surface_prior_experiment")
PRIOR_BANDS = ("B01", "B02", "B04")
SCATTER_MAX_POINTS = 4096


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare monthly_database/SWIR_refine and direct kernel_model "
            "surface priors for the same SIAC scene."
        )
    )
    parser.add_argument("--config", type=Path, required=True, help="SIAC TOML config.")
    parser.add_argument("--input-path", type=Path, required=True, help="Scene SAFE path.")
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Directory for experiment outputs."
    )
    parser.add_argument(
        "--stored-surface-prior",
        type=Path,
        default=None,
        help="Optional stored surface_prior.nc from an existing run for reproducibility checks.",
    )
    parser.add_argument(
        "--sensor",
        default="s2",
        help="Sensor override passed to SIACConfig.with_overrides(...).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="Logging level.",
    )
    return parser.parse_args()


def setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _sample_pairs(
    x: np.ndarray, y: np.ndarray, max_points: int = SCATTER_MAX_POINTS
) -> tuple[np.ndarray, np.ndarray]:
    if x.size <= max_points:
        return x, y
    order = np.argsort(x, kind="mergesort")
    x = x[order]
    y = y[order]
    indices = np.linspace(0, x.size - 1, max_points, dtype=np.int64)
    return x[indices], y[indices]


def _safe_float(value: Any) -> float | None:
    array = np.asarray(value)
    if array.size == 0:
        return None
    scalar = float(array.reshape(-1)[0])
    if not np.isfinite(scalar):
        return None
    return scalar


def _corrcoef(x: np.ndarray, y: np.ndarray) -> float | None:
    if x.size < 2:
        return None
    if float(np.nanstd(x)) == 0.0 or float(np.nanstd(y)) == 0.0:
        return None
    return _safe_float(np.corrcoef(x, y)[0, 1])


def _linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float | None, float | None]:
    if x.size < 2:
        return None, None
    if float(np.nanstd(x)) == 0.0:
        return None, None
    slope, intercept = np.polyfit(x, y, 1)
    return _safe_float(slope), _safe_float(intercept)


def _flatten_valid(
    a: xr.DataArray, b: xr.DataArray, mask: xr.DataArray | np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    valid = np.asarray(mask, dtype=bool) & np.isfinite(a.values) & np.isfinite(b.values)
    return a.values[valid].astype(np.float32, copy=False), b.values[valid].astype(
        np.float32, copy=False
    )


def _pair_metrics(
    a: xr.DataArray, b: xr.DataArray, mask: xr.DataArray | np.ndarray
) -> dict[str, Any]:
    av, bv = _flatten_valid(a, b, mask)
    if av.size == 0:
        return {"valid_count": 0}
    delta = av - bv
    slope, intercept = _linear_fit(bv, av)
    quantiles = np.quantile(delta, [0.05, 0.25, 0.5, 0.75, 0.95])
    return {
        "valid_count": int(av.size),
        "mean_a": _safe_float(np.mean(av)),
        "mean_b": _safe_float(np.mean(bv)),
        "std_a": _safe_float(np.std(av)),
        "std_b": _safe_float(np.std(bv)),
        "bias_a_minus_b": _safe_float(np.mean(delta)),
        "median_bias_a_minus_b": _safe_float(np.median(delta)),
        "rmse": _safe_float(np.sqrt(np.mean(delta**2))),
        "mae": _safe_float(np.mean(np.abs(delta))),
        "corr": _corrcoef(av, bv),
        "fit_slope": slope,
        "fit_intercept": intercept,
        "delta_q05": _safe_float(quantiles[0]),
        "delta_q25": _safe_float(quantiles[1]),
        "delta_q50": _safe_float(quantiles[2]),
        "delta_q75": _safe_float(quantiles[3]),
        "delta_q95": _safe_float(quantiles[4]),
    }


def _toa_metrics(
    observed: xr.DataArray,
    simulated: xr.DataArray,
    mask: xr.DataArray | np.ndarray,
) -> dict[str, Any]:
    obs, sim = _flatten_valid(observed, simulated, mask)
    if obs.size == 0:
        return {"valid_count": 0}
    resid = sim - obs
    return {
        "valid_count": int(obs.size),
        "observed_mean": _safe_float(np.mean(obs)),
        "simulated_mean": _safe_float(np.mean(sim)),
        "bias": _safe_float(np.mean(resid)),
        "median_bias": _safe_float(np.median(resid)),
        "rmse": _safe_float(np.sqrt(np.mean(resid**2))),
        "mae": _safe_float(np.mean(np.abs(resid))),
        "corr_obs_vs_sim": _corrcoef(obs, sim),
        "residual_q05": _safe_float(np.quantile(resid, 0.05)),
        "residual_q50": _safe_float(np.quantile(resid, 0.50)),
        "residual_q95": _safe_float(np.quantile(resid, 0.95)),
    }


def _plot_scatter(
    x: np.ndarray,
    y: np.ndarray,
    *,
    xlabel: str,
    ylabel: str,
    title: str,
    path: Path,
    annotate_lines: list[str],
    diagonal: bool = True,
) -> None:
    xs, ys = _sample_pairs(x, y)
    width = 900
    height = 900
    margin_left = 90
    margin_right = 30
    margin_top = 70
    margin_bottom = 90
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    font = ImageFont.load_default()

    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw.rectangle(
        [margin_left, margin_top, margin_left + plot_width, margin_top + plot_height],
        outline="black",
        width=1,
    )

    if xs.size:
        if diagonal:
            lo = float(np.nanmin([xs.min(), ys.min()]))
            hi = float(np.nanmax([xs.max(), ys.max()]))
            if hi <= lo:
                hi = lo + 1.0
            x_min = y_min = lo
            x_max = y_max = hi
            draw.line(
                [
                    (margin_left, margin_top + plot_height),
                    (margin_left + plot_width, margin_top),
                ],
                fill=(80, 80, 80),
                width=1,
            )
        else:
            x_min = float(xs.min())
            x_max = float(xs.max())
            y_min = float(ys.min())
            y_max = float(ys.max())
            if x_max <= x_min:
                x_max = x_min + 1.0
            if y_max <= y_min:
                y_max = y_min + 1.0
        for xv, yv in zip(xs.tolist(), ys.tolist(), strict=True):
            px = margin_left + int(round((xv - x_min) / (x_max - x_min) * plot_width))
            py = margin_top + plot_height - int(round((yv - y_min) / (y_max - y_min) * plot_height))
            draw.rectangle([px, py, px + 1, py + 1], fill=(31, 119, 180))

    draw.text((margin_left, 15), title, fill="black", font=font)
    draw.text((margin_left, height - 22), xlabel, fill="black", font=font)
    draw.text((10, margin_top), ylabel, fill="black", font=font)
    if annotate_lines:
        draw.multiline_text(
            (margin_left + 8, margin_top + 8),
            "\n".join(annotate_lines),
            fill="black",
            font=font,
            spacing=2,
        )
    image.save(path)


def _plot_map(data: xr.DataArray, *, title: str, path: Path, cmap: str = "RdBu_r") -> None:
    _ = cmap
    values = np.asarray(data.values, dtype=np.float32)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        vlim = 1.0
    else:
        vlim = float(np.nanquantile(np.abs(finite), 0.98))
        if vlim == 0.0:
            vlim = float(np.nanmax(np.abs(finite))) or 1.0
    norm = np.clip(values / vlim, -1.0, 1.0)
    rgb = np.full((*norm.shape, 3), 245, dtype=np.uint8)
    finite_mask = np.isfinite(norm)
    red = np.clip(norm, 0.0, 1.0)
    blue = np.clip(-norm, 0.0, 1.0)
    white = 1.0 - np.abs(norm)
    rgb[..., 0] = np.where(finite_mask, np.round(255.0 * (red + white)).clip(0, 255), 220).astype(
        np.uint8
    )
    rgb[..., 1] = np.where(finite_mask, np.round(255.0 * white).clip(0, 255), 220).astype(np.uint8)
    rgb[..., 2] = np.where(finite_mask, np.round(255.0 * (blue + white)).clip(0, 255), 220).astype(
        np.uint8
    )

    image = Image.fromarray(rgb, mode="RGB")
    if max(image.size) < 600:
        try:
            resample = Image.Resampling.NEAREST
        except AttributeError:
            resample = Image.NEAREST
        scale = max(1, 600 // max(image.size))
        image = image.resize((image.width * scale, image.height * scale), resample=resample)

    font = ImageFont.load_default()
    canvas = Image.new("RGB", (image.width, image.height + 28), "white")
    canvas.paste(image, (0, 28))
    draw = ImageDraw.Draw(canvas)
    draw.text((8, 8), title, fill="black", font=font)
    canvas.save(path)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _surface_prior_band(prior: SurfacePrior, band_name: str) -> xr.DataArray:
    return prior.boa.sel(band=band_name, drop=True).astype(np.float32)


def _surface_prior_unc_band(prior: SurfacePrior, band_name: str) -> xr.DataArray:
    return prior.boa_unc.sel(band=band_name, drop=True).astype(np.float32)


def _build_plan(config_path: Path, input_path: Path, sensor: str) -> tuple[Any, Any, Any]:
    config = SIACConfig.from_file(config_path).with_overrides(sensor=sensor)
    request = SceneProcessRequest(input_path=input_path, output_path=None, config=config)
    plan = build_execution_plan(request)
    observation = plan.preprocessor(request.input_path, None)
    return config, plan, observation


def _compute_direct_prior(
    config: Any, plan: Any, observation: Any, resolution: float
) -> SurfacePrior:
    kernel_config = config.model_copy(deep=True)
    kernel_config.surface_prior.method = "kernel_model"
    provider = resolve_surface_prior_provider(kernel_config, auth=plan.auth)
    return provider(observation, None, plan.rt_model, resolution)


def _simulate_toa(
    solver_inputs: SolverInputBundle,
    atmo_state: AtmosphericState,
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    valid_mask = build_solver_valid_mask(
        solver_inputs.cloud_mask,
        solver_inputs.toa,
        solver_inputs.surface_prior,
    ).astype(bool)

    simulated_bands: list[xr.DataArray] = []
    residual_bands: list[xr.DataArray] = []
    band_masks: list[xr.DataArray] = []
    for band in solver_inputs.bands:
        observed = solver_inputs.toa.sel(band=band.name, drop=True)
        surface = solver_inputs.surface_prior.boa.sel(band=band.name, drop=True)
        coeffs = solver_inputs.rt_model.compute_coefficients(
            solver_inputs.geometry,
            atmo_state,
            band,
            False,
        )
        simulated = coeffs.simulate_toa(surface).astype(np.float32)
        residual = (simulated - observed).astype(np.float32)
        band_valid = (
            valid_mask & np.isfinite(observed) & np.isfinite(surface) & np.isfinite(simulated)
        ).astype(bool)
        simulated_bands.append(simulated.expand_dims(band=[band.name]))
        residual_bands.append(residual.expand_dims(band=[band.name]))
        band_masks.append(band_valid.expand_dims(band=[band.name]))

    return (
        xr.concat(simulated_bands, dim="band").astype(np.float32),
        xr.concat(residual_bands, dim="band").astype(np.float32),
        xr.concat(band_masks, dim="band").astype(bool),
    )


def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Building execution plan and observation bundle")
    config, plan, observation = _build_plan(args.config, args.input_path, args.sensor)
    resolution = aerosol_resolution(config)

    LOGGER.info("Fetching atmospheric prior at %.1f m", resolution)
    atmo_prior = plan.atmo_provider(
        observation.bounds,
        observation.crs,
        observation.metadata["observation_time"],
        resolution,
    )

    LOGGER.info("Computing monthly_database surface prior")
    monthly_prior = plan.surface_prior_provider(observation, atmo_prior, plan.rt_model, resolution)
    LOGGER.info("Computing direct kernel_model surface prior")
    direct_prior = _compute_direct_prior(config, plan, observation, resolution)

    LOGGER.info("Assembling solver grids")
    monthly_inputs = plan.grid_assembler(
        observation,
        atmo_prior,
        monthly_prior,
        plan.rt_model,
        aerosol_resolution_m=resolution,
    )
    direct_inputs = plan.grid_assembler(
        observation,
        atmo_prior,
        direct_prior,
        plan.rt_model,
        aerosol_resolution_m=resolution,
    )

    LOGGER.info("Running baseline solver with monthly_database prior")
    monthly_solved = plan.solver(monthly_inputs, config)
    LOGGER.info("Running counterfactual solver with direct kernel_model prior")
    direct_solved = plan.solver(direct_inputs, config)

    LOGGER.info("Computing TOA diagnostics")
    toa_monthly_prior_atmo, resid_monthly_prior_atmo, valid_monthly_prior_atmo = _simulate_toa(
        monthly_inputs,
        monthly_inputs.atmo_prior,
    )
    toa_direct_prior_atmo, resid_direct_prior_atmo, valid_direct_prior_atmo = _simulate_toa(
        direct_inputs,
        direct_inputs.atmo_prior,
    )
    toa_monthly_fixed, resid_monthly_fixed, valid_monthly_fixed = _simulate_toa(
        monthly_inputs,
        monthly_solved.atmo_state,
    )
    toa_direct_fixed, resid_direct_fixed, valid_direct_fixed = _simulate_toa(
        direct_inputs,
        monthly_solved.atmo_state,
    )
    toa_direct_own, resid_direct_own, valid_direct_own = _simulate_toa(
        direct_inputs,
        direct_solved.atmo_state,
    )

    summary: dict[str, Any] = {
        "scene": {
            "config": str(args.config),
            "input_path": str(args.input_path),
            "output_dir": str(args.output_dir),
            "stored_surface_prior": str(args.stored_surface_prior)
            if args.stored_surface_prior
            else None,
            "observation_time": str(observation.metadata["observation_time"]),
            "bounds": [float(v) for v in observation.bounds],
            "crs": observation.crs,
            "aerosol_resolution_m": float(resolution),
            "solver_bands": [band.name for band in monthly_inputs.bands],
        },
        "prior_metrics": {},
        "stored_repro_metrics": {},
        "toa_metrics": {
            "atmo_prior": {},
            "monthly_solved_fixed_atmo": {},
            "end_to_end": {},
        },
        "attribution_metrics": {
            "atmo_prior": {},
            "monthly_solved_fixed_atmo": {},
            "end_to_end": {},
        },
    }

    prior_rows: list[dict[str, Any]] = []
    for band_name in PRIOR_BANDS:
        monthly_band = _surface_prior_band(monthly_prior, band_name)
        direct_band = _surface_prior_band(direct_prior, band_name)
        common_mask = (
            monthly_prior.mask.astype(bool)
            & direct_prior.mask.astype(bool)
            & np.isfinite(monthly_band)
            & np.isfinite(direct_band)
        )
        metrics = _pair_metrics(monthly_band, direct_band, common_mask)
        summary["prior_metrics"][band_name] = metrics
        prior_rows.append({"band": band_name, **metrics})

        x, y = _flatten_valid(direct_band, monthly_band, common_mask)
        annotate_lines = [
            f"n={metrics.get('valid_count', 0)}",
            f"bias={metrics.get('bias_a_minus_b'):.5f}"
            if metrics.get("bias_a_minus_b") is not None
            else "bias=n/a",
            f"rmse={metrics.get('rmse'):.5f}" if metrics.get("rmse") is not None else "rmse=n/a",
            f"corr={metrics.get('corr'):.3f}" if metrics.get("corr") is not None else "corr=n/a",
        ]
        _plot_scatter(
            x,
            y,
            xlabel=f"Direct MCD43 {band_name}",
            ylabel=f"Monthly database {band_name}",
            title=f"Surface Prior Scatter {band_name}",
            path=args.output_dir / f"prior_scatter_{band_name}.png",
            annotate_lines=annotate_lines,
        )
        _plot_map(
            (monthly_band - direct_band).astype(np.float32),
            title=f"Monthly - Direct Surface Prior {band_name}",
            path=args.output_dir / f"prior_delta_{band_name}.png",
        )

    if args.stored_surface_prior is not None and args.stored_surface_prior.exists():
        stored = xr.open_dataset(args.stored_surface_prior)
        for band_name in PRIOR_BANDS:
            if band_name not in stored.data_vars:
                continue
            stored_band = stored[band_name].astype(np.float32)
            monthly_band = _surface_prior_band(monthly_prior, band_name)
            mask = np.isfinite(stored_band) & np.isfinite(monthly_band)
            summary["stored_repro_metrics"][band_name] = _pair_metrics(
                monthly_band, stored_band, mask
            )
        stored.close()

    toa_rows: list[dict[str, Any]] = []
    simulated_scenarios = {
        "atmo_prior": (
            toa_monthly_prior_atmo,
            toa_direct_prior_atmo,
            resid_monthly_prior_atmo,
            resid_direct_prior_atmo,
            valid_monthly_prior_atmo,
            valid_direct_prior_atmo,
        ),
        "monthly_solved_fixed_atmo": (
            toa_monthly_fixed,
            toa_direct_fixed,
            resid_monthly_fixed,
            resid_direct_fixed,
            valid_monthly_fixed,
            valid_direct_fixed,
        ),
        "end_to_end": (
            toa_monthly_fixed,
            toa_direct_own,
            resid_monthly_fixed,
            resid_direct_own,
            valid_monthly_fixed,
            valid_direct_own,
        ),
    }

    for scenario_name, (
        sim_monthly,
        sim_direct,
        resid_monthly,
        resid_direct,
        valid_monthly,
        valid_direct,
    ) in simulated_scenarios.items():
        for band_name in sim_monthly.coords["band"].values.tolist():
            observed = monthly_inputs.toa.sel(band=band_name, drop=True).astype(np.float32)
            monthly_band = sim_monthly.sel(band=band_name, drop=True)
            direct_band = sim_direct.sel(band=band_name, drop=True)
            monthly_mask = valid_monthly.sel(band=band_name, drop=True)
            direct_mask = valid_direct.sel(band=band_name, drop=True)
            common_mask = monthly_mask.astype(bool) & direct_mask.astype(bool)

            monthly_metrics = _toa_metrics(observed, monthly_band, common_mask)
            direct_metrics = _toa_metrics(observed, direct_band, common_mask)
            monthly_rmse = monthly_metrics.get("rmse")
            direct_rmse = direct_metrics.get("rmse")
            rmse_improvement_pct = None
            if monthly_rmse not in (None, 0.0) and direct_rmse is not None:
                rmse_improvement_pct = (
                    100.0 * (float(monthly_rmse) - float(direct_rmse)) / float(monthly_rmse)
                )

            summary["toa_metrics"][scenario_name][band_name] = {
                "monthly_database": monthly_metrics,
                "direct_kernel_model": direct_metrics,
                "rmse_improvement_pct": rmse_improvement_pct,
            }
            toa_rows.append(
                {
                    "scenario": scenario_name,
                    "band": band_name,
                    "monthly_bias": monthly_metrics.get("bias"),
                    "direct_bias": direct_metrics.get("bias"),
                    "monthly_rmse": monthly_metrics.get("rmse"),
                    "direct_rmse": direct_metrics.get("rmse"),
                    "rmse_improvement_pct": rmse_improvement_pct,
                }
            )

            monthly_resid = resid_monthly.sel(band=band_name, drop=True)
            direct_resid = resid_direct.sel(band=band_name, drop=True)
            direct_surface = direct_inputs.surface_prior.boa.sel(band=band_name, drop=True)
            monthly_surface = monthly_inputs.surface_prior.boa.sel(band=band_name, drop=True)
            delta_prior = (monthly_surface - direct_surface).astype(np.float32)
            delta_resid = (monthly_resid - direct_resid).astype(np.float32)
            dx, dy = _flatten_valid(delta_prior, delta_resid, common_mask)
            slope, intercept = _linear_fit(dx, dy)
            attribution = {
                "valid_count": int(dx.size),
                "corr_delta_prior_vs_delta_residual": _corrcoef(dx, dy),
                "fit_slope": slope,
                "fit_intercept": intercept,
            }
            summary["attribution_metrics"][scenario_name][band_name] = attribution

            obs_vals, monthly_vals = _flatten_valid(observed, monthly_band, common_mask)
            _, direct_vals = _flatten_valid(observed, direct_band, common_mask)
            _plot_scatter(
                obs_vals,
                monthly_vals,
                xlabel=f"Observed TOA {band_name}",
                ylabel="Simulated TOA",
                title=f"{scenario_name} Monthly TOA Scatter {band_name}",
                path=args.output_dir / f"toa_scatter_{scenario_name}_monthly_{band_name}.png",
                annotate_lines=[
                    f"bias={monthly_metrics.get('bias'):.5f}"
                    if monthly_metrics.get("bias") is not None
                    else "bias=n/a",
                    f"rmse={monthly_metrics.get('rmse'):.5f}"
                    if monthly_metrics.get("rmse") is not None
                    else "rmse=n/a",
                ],
            )
            _plot_scatter(
                obs_vals,
                direct_vals,
                xlabel=f"Observed TOA {band_name}",
                ylabel="Simulated TOA",
                title=f"{scenario_name} Direct TOA Scatter {band_name}",
                path=args.output_dir / f"toa_scatter_{scenario_name}_direct_{band_name}.png",
                annotate_lines=[
                    f"bias={direct_metrics.get('bias'):.5f}"
                    if direct_metrics.get("bias") is not None
                    else "bias=n/a",
                    f"rmse={direct_metrics.get('rmse'):.5f}"
                    if direct_metrics.get("rmse") is not None
                    else "rmse=n/a",
                ],
            )
            _plot_map(
                (np.abs(monthly_resid) - np.abs(direct_resid)).astype(np.float32),
                title=f"|Monthly resid| - |Direct resid| {scenario_name} {band_name}",
                path=args.output_dir
                / f"toa_abs_residual_improvement_{scenario_name}_{band_name}.png",
            )
            _plot_scatter(
                dx,
                dy,
                xlabel="Delta surface prior (monthly - direct)",
                ylabel="Delta TOA residual (monthly - direct)",
                title=f"Attribution {scenario_name} {band_name}",
                path=args.output_dir / f"toa_attribution_{scenario_name}_{band_name}.png",
                annotate_lines=[
                    f"n={attribution['valid_count']}",
                    (
                        f"corr={attribution['corr_delta_prior_vs_delta_residual']:.3f}"
                        if attribution["corr_delta_prior_vs_delta_residual"] is not None
                        else "corr=n/a"
                    ),
                    f"slope={attribution['fit_slope']:.3f}"
                    if attribution["fit_slope"] is not None
                    else "slope=n/a",
                ],
                diagonal=False,
            )

    summary_path = args.output_dir / "summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    _write_csv(args.output_dir / "prior_metrics.csv", prior_rows)
    _write_csv(args.output_dir / "toa_metrics.csv", toa_rows)

    priors_ds = xr.Dataset(
        {
            "monthly_boa": monthly_prior.boa.astype(np.float32),
            "monthly_boa_unc": monthly_prior.boa_unc.astype(np.float32),
            "monthly_mask": monthly_prior.mask.astype(np.uint8),
            "direct_boa": direct_prior.boa.astype(np.float32),
            "direct_boa_unc": direct_prior.boa_unc.astype(np.float32),
            "direct_mask": direct_prior.mask.astype(np.uint8),
            "delta_boa": (monthly_prior.boa - direct_prior.boa).astype(np.float32),
        }
    )
    priors_ds.to_netcdf(args.output_dir / "surface_prior_comparison.nc")

    toa_ds = xr.Dataset(
        {
            "observed_toa": monthly_inputs.toa.astype(np.float32),
            "monthly_surface_solver": monthly_inputs.surface_prior.boa.astype(np.float32),
            "direct_surface_solver": direct_inputs.surface_prior.boa.astype(np.float32),
            "toa_monthly_atmo_prior": toa_monthly_prior_atmo.astype(np.float32),
            "toa_direct_atmo_prior": toa_direct_prior_atmo.astype(np.float32),
            "toa_monthly_fixed": toa_monthly_fixed.astype(np.float32),
            "toa_direct_fixed": toa_direct_fixed.astype(np.float32),
            "toa_direct_own": toa_direct_own.astype(np.float32),
            "resid_monthly_atmo_prior": resid_monthly_prior_atmo.astype(np.float32),
            "resid_direct_atmo_prior": resid_direct_prior_atmo.astype(np.float32),
            "resid_monthly_fixed": resid_monthly_fixed.astype(np.float32),
            "resid_direct_fixed": resid_direct_fixed.astype(np.float32),
            "resid_direct_own": resid_direct_own.astype(np.float32),
            "valid_monthly_fixed": valid_monthly_fixed.astype(np.uint8),
            "valid_direct_fixed": valid_direct_fixed.astype(np.uint8),
            "monthly_aot_solved": monthly_solved.atmo_state.aot.astype(np.float32),
            "monthly_tcwv_solved": monthly_solved.atmo_state.tcwv.astype(np.float32),
            "direct_aot_solved": direct_solved.atmo_state.aot.astype(np.float32),
            "direct_tcwv_solved": direct_solved.atmo_state.tcwv.astype(np.float32),
        }
    )
    toa_ds.to_netcdf(args.output_dir / "toa_comparison.nc")

    LOGGER.info("Experiment complete")
    LOGGER.info("Summary written to %s", summary_path)


if __name__ == "__main__":
    main()
