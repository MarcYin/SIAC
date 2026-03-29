#!/usr/bin/env python3
"""Run a pixelwise RT closure analysis from stored LUT/emulator comparison outputs."""

from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

LOGGER = logging.getLogger("pixelwise_rt_closure")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyse pixelwise RT closure using stored compare_rt_backends outputs, "
            "surface priors, and sampled per-pixel coefficient tables."
        )
    )
    parser.add_argument(
        "--reference-output-dir",
        type=Path,
        required=True,
        help="Scene output directory containing surface_prior.nc.",
    )
    parser.add_argument(
        "--rt-output-dir",
        type=Path,
        required=True,
        help="Directory produced by compare_rt_backends.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where closure summaries and CSV samples will be written.",
    )
    parser.add_argument(
        "--samples-per-band",
        type=int,
        default=21,
        help="Number of quantile-stratified pixel samples to export per band and scenario.",
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


def _corrcoef(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 2 or b.size < 2:
        return float("nan")
    if np.allclose(a, a[0]) or np.allclose(b, b[0]):
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _linear_fit(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    if x.size < 2 or y.size < 2:
        return {"slope": float("nan"), "intercept": float("nan")}
    slope, intercept = np.polyfit(x, y, 1)
    return {"slope": float(slope), "intercept": float(intercept)}


def _stats(values: np.ndarray) -> dict[str, float]:
    return {
        "min": float(np.min(values)),
        "p05": float(np.quantile(values, 0.05)),
        "p50": float(np.quantile(values, 0.50)),
        "p95": float(np.quantile(values, 0.95)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
    }


def _simulate_toa(
    xap: np.ndarray,
    xbp: np.ndarray,
    xcp: np.ndarray,
    boa: np.ndarray,
) -> np.ndarray:
    denom = 1.0 - xcp * boa
    stable = np.isfinite(denom) & (np.abs(denom) > 1.0e-6) & np.isfinite(xap) & (np.abs(xap) > 1.0e-12)
    y = boa / denom
    toa = (y + xbp) / xap
    return np.where(stable, toa, np.nan).astype(np.float32)


def _scenario_names(rt_output_dir: Path) -> list[str]:
    scenarios = []
    for path in sorted(rt_output_dir.glob("toa_*.nc")):
        scenarios.append(path.stem.removeprefix("toa_"))
    if not scenarios:
        raise FileNotFoundError(f"No toa_*.nc files found under {rt_output_dir}")
    return scenarios


def _band_names(toa_ds: xr.Dataset) -> list[str]:
    bands = []
    prefix = "observed_toa_"
    for name in toa_ds.data_vars:
        if name.startswith(prefix):
            bands.append(name.removeprefix(prefix))
    if not bands:
        raise ValueError("Could not infer bands from TOA dataset.")
    return sorted(bands)


def _positive_mask(obs: np.ndarray, valid: np.ndarray) -> np.ndarray:
    return valid & np.isfinite(obs) & (obs > 0.0)


def _summary_for_backend(
    obs: np.ndarray,
    sim: np.ndarray,
    valid: np.ndarray,
) -> dict[str, Any]:
    obs_vals = obs[valid]
    sim_vals = sim[valid]
    residual = sim_vals - obs_vals
    fit = _linear_fit(obs_vals, sim_vals)
    return {
        "stats": _stats(sim_vals),
        "bias_simulated_minus_observed": float(np.mean(residual)),
        "rmse_simulated_minus_observed": float(np.sqrt(np.mean(residual**2))),
        "mae_simulated_minus_observed": float(np.mean(np.abs(residual))),
        "corr": _corrcoef(obs_vals, sim_vals),
        "linear_fit_observed_to_simulated": fit,
    }


def _select_quantile_samples(obs: np.ndarray, valid: np.ndarray, count: int) -> list[tuple[int, float]]:
    candidate_mask = _positive_mask(obs, valid)
    candidate_indices = np.flatnonzero(candidate_mask.ravel())
    if candidate_indices.size == 0:
        candidate_indices = np.flatnonzero(valid.ravel())
    if candidate_indices.size == 0:
        return []

    flat_obs = obs.ravel()
    candidate_values = flat_obs[candidate_indices]
    sample_count = min(int(count), int(candidate_indices.size))
    if sample_count <= 0:
        return []

    quantiles = np.linspace(0.02, 0.98, sample_count)
    targets = np.quantile(candidate_values, quantiles)

    chosen: list[tuple[int, float]] = []
    remaining_indices = candidate_indices.copy()
    remaining_values = candidate_values.copy()
    for quantile, target in zip(quantiles, targets):
        pick = int(np.argmin(np.abs(remaining_values - target)))
        chosen.append((int(remaining_indices[pick]), float(quantile)))
        remaining_indices = np.delete(remaining_indices, pick)
        remaining_values = np.delete(remaining_values, pick)
        if remaining_indices.size == 0:
            break
    return chosen


def _write_samples_csv(
    path: Path,
    *,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
    elevation: np.ndarray | None,
    obs: np.ndarray,
    boa: np.ndarray,
    lut_xap: np.ndarray,
    lut_xbp: np.ndarray,
    lut_xcp: np.ndarray,
    lut_toa: np.ndarray,
    emu_xap: np.ndarray,
    emu_xbp: np.ndarray,
    emu_xcp: np.ndarray,
    emu_toa: np.ndarray,
    valid: np.ndarray,
    samples: list[tuple[int, float]],
) -> None:
    fieldnames = [
        "sample_index",
        "target_quantile",
        "y_index",
        "x_index",
        "y_coord",
        "x_coord",
        "elevation_km",
        "observed_toa",
        "surface_boa",
        "lut_xap",
        "lut_xbp",
        "lut_xcp",
        "lut_forward_toa",
        "lut_residual",
        "emu_xap",
        "emu_xbp",
        "emu_xcp",
        "emu_forward_toa",
        "emu_residual",
        "emu_minus_lut_toa",
        "valid",
    ]

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for sample_index, (flat_index, quantile) in enumerate(samples, start=1):
            y_index, x_index = np.unravel_index(flat_index, obs.shape)
            writer.writerow(
                {
                    "sample_index": sample_index,
                    "target_quantile": quantile,
                    "y_index": int(y_index),
                    "x_index": int(x_index),
                    "y_coord": float(y_coords[y_index]),
                    "x_coord": float(x_coords[x_index]),
                    "elevation_km": float(elevation[y_index, x_index]) if elevation is not None else float("nan"),
                    "observed_toa": float(obs[y_index, x_index]),
                    "surface_boa": float(boa[y_index, x_index]),
                    "lut_xap": float(lut_xap[y_index, x_index]),
                    "lut_xbp": float(lut_xbp[y_index, x_index]),
                    "lut_xcp": float(lut_xcp[y_index, x_index]),
                    "lut_forward_toa": float(lut_toa[y_index, x_index]),
                    "lut_residual": float(lut_toa[y_index, x_index] - obs[y_index, x_index]),
                    "emu_xap": float(emu_xap[y_index, x_index]),
                    "emu_xbp": float(emu_xbp[y_index, x_index]),
                    "emu_xcp": float(emu_xcp[y_index, x_index]),
                    "emu_forward_toa": float(emu_toa[y_index, x_index]),
                    "emu_residual": float(emu_toa[y_index, x_index] - obs[y_index, x_index]),
                    "emu_minus_lut_toa": float(emu_toa[y_index, x_index] - lut_toa[y_index, x_index]),
                    "valid": bool(valid[y_index, x_index]),
                }
            )


def main() -> None:
    args = _parse_args()
    _setup_logging(args.log_level)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    surface_prior = xr.open_dataset(args.reference_output_dir / "surface_prior.nc")
    scenarios = _scenario_names(args.rt_output_dir)

    summary: dict[str, Any] = {
        "reference_output_dir": str(args.reference_output_dir.resolve()),
        "rt_output_dir": str(args.rt_output_dir.resolve()),
        "samples_per_band": int(args.samples_per_band),
        "scenarios": {},
    }

    for scenario in scenarios:
        LOGGER.info("Analysing scenario %s", scenario)
        toa_ds = xr.open_dataset(args.rt_output_dir / f"toa_{scenario}.nc")
        coeff_ds = xr.open_dataset(args.rt_output_dir / f"coefficients_{scenario}.nc")
        bands = _band_names(toa_ds)

        y_coords = np.asarray(toa_ds.coords["y"].values)
        x_coords = np.asarray(toa_ds.coords["x"].values)
        elevation = None
        if "elevation_km" in coeff_ds:
            elevation = coeff_ds["elevation_km"].transpose("y", "x").values.astype(np.float32)

        scenario_summary: dict[str, Any] = {
            "elevation_km": _stats(elevation[np.isfinite(elevation)]) if elevation is not None else None
        }

        for band in bands:
            obs = toa_ds[f"observed_toa_{band}"].transpose("y", "x").values.astype(np.float32)
            lut_toa = toa_ds[f"simulated_toa_lut_{band}"].transpose("y", "x").values.astype(np.float32)
            emu_toa = toa_ds[f"simulated_toa_emulator_{band}"].transpose("y", "x").values.astype(np.float32)
            boa = surface_prior[band].transpose("y", "x").values.astype(np.float32)

            lut_xap = coeff_ds[f"{band}_lut_xap"].transpose("y", "x").values.astype(np.float32)
            lut_xbp = coeff_ds[f"{band}_lut_xbp"].transpose("y", "x").values.astype(np.float32)
            lut_xcp = coeff_ds[f"{band}_lut_xcp"].transpose("y", "x").values.astype(np.float32)
            emu_xap = coeff_ds[f"{band}_emulator_xap"].transpose("y", "x").values.astype(np.float32)
            emu_xbp = coeff_ds[f"{band}_emulator_xbp"].transpose("y", "x").values.astype(np.float32)
            emu_xcp = coeff_ds[f"{band}_emulator_xcp"].transpose("y", "x").values.astype(np.float32)

            lut_toa_rebuilt = _simulate_toa(lut_xap, lut_xbp, lut_xcp, boa)
            emu_toa_rebuilt = _simulate_toa(emu_xap, emu_xbp, emu_xcp, boa)

            valid = (
                np.isfinite(obs)
                & np.isfinite(boa)
                & np.isfinite(lut_toa)
                & np.isfinite(emu_toa)
                & np.isfinite(lut_toa_rebuilt)
                & np.isfinite(emu_toa_rebuilt)
            )
            positive_valid = _positive_mask(obs, valid)
            if not np.any(valid):
                raise ValueError(f"No valid pixels found for {scenario}/{band}")

            obs_vals = obs[valid]
            lut_vals = lut_toa[valid]
            emu_vals = emu_toa[valid]
            pos_obs_vals = obs[positive_valid] if np.any(positive_valid) else obs_vals
            pos_lut_vals = lut_toa[positive_valid] if np.any(positive_valid) else lut_vals
            pos_emu_vals = emu_toa[positive_valid] if np.any(positive_valid) else emu_vals

            scenario_summary[band] = {
                "valid_count": int(np.count_nonzero(valid)),
                "positive_observed_count": int(np.count_nonzero(positive_valid)),
                "observed_all": _stats(obs_vals),
                "observed_positive": _stats(pos_obs_vals),
                "lut": _summary_for_backend(obs, lut_toa, valid),
                "emulator": _summary_for_backend(obs, emu_toa, valid),
                "lut_positive": _summary_for_backend(obs, lut_toa, positive_valid if np.any(positive_valid) else valid),
                "emulator_positive": _summary_for_backend(
                    obs, emu_toa, positive_valid if np.any(positive_valid) else valid
                ),
                "range_ratios_positive": {
                    "lut_p95_over_observed_p95": float(
                        np.quantile(pos_lut_vals, 0.95) / np.quantile(pos_obs_vals, 0.95)
                    ),
                    "emulator_p95_over_observed_p95": float(
                        np.quantile(pos_emu_vals, 0.95) / np.quantile(pos_obs_vals, 0.95)
                    ),
                },
                "forward_formula_check": {
                    "lut_max_abs_diff": float(np.nanmax(np.abs(lut_toa_rebuilt[valid] - lut_toa[valid]))),
                    "lut_rmse_diff": float(np.sqrt(np.nanmean((lut_toa_rebuilt[valid] - lut_toa[valid]) ** 2))),
                    "emulator_max_abs_diff": float(np.nanmax(np.abs(emu_toa_rebuilt[valid] - emu_toa[valid]))),
                    "emulator_rmse_diff": float(np.sqrt(np.nanmean((emu_toa_rebuilt[valid] - emu_toa[valid]) ** 2))),
                },
            }

            samples = _select_quantile_samples(obs, valid, args.samples_per_band)
            csv_path = args.output_dir / f"samples_{scenario}_{band}.csv"
            _write_samples_csv(
                csv_path,
                x_coords=x_coords,
                y_coords=y_coords,
                elevation=elevation,
                obs=obs,
                boa=boa,
                lut_xap=lut_xap,
                lut_xbp=lut_xbp,
                lut_xcp=lut_xcp,
                lut_toa=lut_toa_rebuilt,
                emu_xap=emu_xap,
                emu_xbp=emu_xbp,
                emu_xcp=emu_xcp,
                emu_toa=emu_toa_rebuilt,
                valid=valid,
                samples=samples,
            )

        summary["scenarios"][scenario] = scenario_summary

    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    LOGGER.info("Pixelwise closure analysis written to %s", args.output_dir)


if __name__ == "__main__":
    main()
