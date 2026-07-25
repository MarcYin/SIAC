"""Evaluate seasonal L1C priors against saved truth-AOD BOA references."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from pyproj import Transformer

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
BAND_INDEX = {"B02": 1, "B04": 3}


def _parse_prior(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("prior must be LABEL=PATH")
    label, path = value.split("=", 1)
    return label, Path(path)


def _window_mask(
    *,
    shape: tuple[int, int],
    transform: np.ndarray,
    epsg: int,
    lon: float,
    lat: float,
    radius_m: float,
) -> np.ndarray:
    a, _, xoff, _, e, yoff = [float(value) for value in transform]
    y_count, x_count = shape
    x = xoff + (np.arange(x_count) + 0.5) * a
    y = yoff + (np.arange(y_count) + 0.5) * e
    center_x, center_y = Transformer.from_crs(
        "EPSG:4326", f"EPSG:{epsg}", always_xy=True
    ).transform(lon, lat)
    return (np.abs(x[np.newaxis, :] - center_x) <= radius_m) & (
        np.abs(y[:, np.newaxis] - center_y) <= radius_m
    )


def _prior_window(
    path: Path,
    *,
    lon: float,
    lat: float,
    radius_m: float,
) -> dict[str, dict[str, float]]:
    with np.load(path, allow_pickle=False) as archive:
        comp = np.asarray(archive["comp"], dtype=np.float64)
        epsg = int(archive["epsg"])
        transform = np.asarray(archive["transform"], dtype=np.float64)
        clean_quality = (
            json.loads(str(archive["clean_quality_json"].item()))
            if "clean_quality_json" in archive.files
            else []
        )
    mask = _window_mask(
        shape=comp.shape[2:],
        transform=transform,
        epsg=epsg,
        lon=lon,
        lat=lat,
        radius_m=radius_m,
    )
    output = {}
    for band, band_index in BAND_INDEX.items():
        values = comp[:, band_index]
        median_surface = np.nanmedian(values, axis=0)
        temporal_rmse = np.sqrt(
            np.nanmean(np.square(values - median_surface[np.newaxis, ...]), axis=0)
        )
        output[band] = {
            "prior": float(np.nanmean(median_surface[mask])),
            "temporal_rmse": float(np.nanmedian(temporal_rmse[mask])),
        }
    selected_excess = [
        max(0.0, float(row.get("selected_aod_max", 0.15)) - 0.15) for row in clean_quality
    ]
    output["quality"] = {
        "windows": float(len(clean_quality)),
        "fallback_scenes": float(sum(int(row.get("n_fallback", 0)) for row in clean_quality)),
        "median_aod_excess": float(np.median(selected_excess)) if selected_excess else 0.0,
    }
    return output


def _summarize(rows: list[dict[str, Any]], model_error: float) -> dict[str, Any]:
    bands = {}
    for band in BAND_INDEX:
        residual = np.asarray([row[band]["residual"] for row in rows], dtype=float)
        temporal = np.asarray([row[band]["temporal_rmse"] for row in rows], dtype=float)
        sigma = np.sqrt(np.square(temporal) + model_error**2)
        z = np.abs(residual) / np.maximum(sigma, 1e-6)
        valid = np.isfinite(residual) & np.isfinite(z)
        bands[band] = {
            "count": int(valid.sum()),
            "bias": float(np.mean(residual[valid])),
            "median_abs_error": float(np.median(np.abs(residual[valid]))),
            "rmse": float(np.sqrt(np.mean(np.square(residual[valid])))),
            "median_temporal_rmse": float(np.median(temporal[valid])),
            "median_z": float(np.median(z[valid])),
            "one_sigma_coverage_pct": float(100.0 * np.mean(z[valid] <= 1.0)),
            "two_sigma_coverage_pct": float(100.0 * np.mean(z[valid] <= 2.0)),
        }
    return {
        "count": len(rows),
        "model_error": model_error,
        "bands": bands,
        "median_fallback_scenes": float(
            np.median([row["quality"]["fallback_scenes"] for row in rows])
        ),
        "median_clean_aod_excess": float(
            np.median([row["quality"]["median_aod_excess"] for row in rows])
        ),
    }


def _paired_summary(
    baseline_label: str,
    baseline_rows: list[dict[str, Any]],
    candidate_label: str,
    candidate_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    baseline = {row["matchup_id"]: row for row in baseline_rows}
    candidate = {row["matchup_id"]: row for row in candidate_rows}
    common = sorted(set(baseline) & set(candidate))
    bands = {}
    for band in BAND_INDEX:
        baseline_abs = np.asarray(
            [abs(float(baseline[mid][band]["residual"])) for mid in common], dtype=float
        )
        candidate_abs = np.asarray(
            [abs(float(candidate[mid][band]["residual"])) for mid in common], dtype=float
        )
        delta = candidate_abs - baseline_abs
        bands[band] = {
            "count": len(common),
            "baseline_median_abs_error": float(np.median(baseline_abs)),
            "candidate_median_abs_error": float(np.median(candidate_abs)),
            "median_abs_error_delta": float(np.median(delta)),
            "baseline_rmse": float(np.sqrt(np.mean(np.square(baseline_abs)))),
            "candidate_rmse": float(np.sqrt(np.mean(np.square(candidate_abs)))),
            "improved": int(np.sum(delta < -1e-12)),
            "tied": int(np.sum(np.abs(delta) <= 1e-12)),
            "worse": int(np.sum(delta > 1e-12)),
        }
    return {
        "baseline": baseline_label,
        "candidate": candidate_label,
        "common": len(common),
        "bands": bands,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prior", type=_parse_prior, action="append", required=True)
    parser.add_argument("--mids", type=Path, default=ROOT / "surface_e1_matched_mids.txt")
    parser.add_argument("--radius-m", type=float, default=1500.0)
    parser.add_argument("--model-errors", default="0,0.010,0.013,0.015")
    parser.add_argument("--output", type=Path, default=ROOT / "l1c_prior_quality_e2_20260709.json")
    args = parser.parse_args()

    matchup_rows = {
        row["matchup_id"]: row
        for row in csv.DictReader((ROOT / "matchups" / "matchups.csv").open())
    }
    matchup_ids = [line.strip() for line in args.mids.read_text().splitlines() if line.strip()]
    references = {
        matchup_id: json.loads((ROOT / "prior_quality" / f"{matchup_id}.json").read_text())
        for matchup_id in matchup_ids
        if (ROOT / "prior_quality" / f"{matchup_id}.json").exists()
    }
    model_errors = [float(value) for value in args.model_errors.split(",")]
    result = {
        "mids": str(args.mids),
        "radius_m": args.radius_m,
        "priors": {},
        "paired": [],
    }
    rows_by_label: dict[str, list[dict[str, Any]]] = {}
    for label, directory in args.prior:
        rows = []
        for matchup_id in matchup_ids:
            path = directory / f"{matchup_id}.npz"
            if matchup_id not in references or not path.exists():
                continue
            match = matchup_rows[matchup_id]
            prior = _prior_window(
                path,
                lon=float(match["longitude"]),
                lat=float(match["latitude"]),
                radius_m=args.radius_m,
            )
            reference = references[matchup_id]["ref"]
            row: dict[str, Any] = {
                "matchup_id": matchup_id,
                "truth": float(references[matchup_id]["truth"]),
                "quality": prior["quality"],
            }
            for band in BAND_INDEX:
                row[band] = {
                    **prior[band],
                    "reference": float(reference[band]),
                    "residual": float(prior[band]["prior"] - reference[band]),
                }
            rows.append(row)
        summaries = [_summarize(rows, model_error) for model_error in model_errors]
        result["priors"][label] = {
            "directory": str(directory),
            "available": len(rows),
            "summaries": summaries,
            "rows": rows,
        }
        rows_by_label[label] = rows
        print(f"{label}: {len(rows)} priors", flush=True)
    baseline_label = args.prior[0][0]
    result["paired"] = [
        _paired_summary(
            baseline_label,
            rows_by_label[baseline_label],
            candidate_label,
            rows_by_label[candidate_label],
        )
        for candidate_label, _ in args.prior[1:]
    ]
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
