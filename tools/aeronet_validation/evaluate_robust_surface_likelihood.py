"""Offline screen of robust surface-only scene likelihoods.

The saved R2 diagnostic cubes contain absolute per-band BOA residuals at every
AOD node. This script removes the aerosol-product backstop and recomputes one
scene likelihood using diagonal chi-square, Huber, or Student-t losses. It is a
screen for E3, not an independent validation: AERONET is used only after each
curve is minimized to score the frozen candidates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
RESULTS = ROOT / "phaseD_results_campaign250_R2_full_localdiag_20260705"
CUBES = ROOT / "phaseD_cost_cubes_campaign250_R2_full_localdiag_20260705"
DEFAULT_OUTPUT = ROOT / "robust_surface_likelihood_e3_20260709.json"
AOD_BINS = (
    (-math.inf, 0.1, "<0.1"),
    (0.1, 0.2, "0.1-0.2"),
    (0.2, 0.4, "0.2-0.4"),
    (0.4, 0.6, "0.4-0.6"),
    (0.6, 1.0, "0.6-1.0"),
    (1.0, 1.5, "1.0-1.5"),
    (1.5, math.inf, ">=1.5"),
)


@dataclass(frozen=True)
class Candidate:
    loss: str
    aggregation: str
    model_error: float
    interpolate: bool

    @property
    def name(self) -> str:
        interpolation = "interp" if self.interpolate else "node"
        return f"{self.loss}_{self.aggregation}_model{self.model_error:.3f}_{interpolation}"


def robust_loss(z: np.ndarray, mode: str) -> np.ndarray:
    """Return a monotonic robust loss for standardized absolute residuals."""

    if mode == "chi2":
        return 0.5 * np.square(z)
    if mode == "huber1p5":
        delta = 1.5
        return np.where(z <= delta, 0.5 * np.square(z), delta * (z - 0.5 * delta))
    if mode == "student3":
        degrees = 3.0
        return 0.5 * (degrees + 1.0) * np.log1p(np.square(z) / degrees)
    raise ValueError(f"unknown robust loss {mode!r}")


def aggregate_curve(pixel_loss: np.ndarray, mode: str) -> np.ndarray:
    """Aggregate a ``(node, pixel)`` loss matrix into one scene curve."""

    if mode == "mean":
        return np.nanmean(pixel_loss, axis=1)
    if mode == "median":
        return np.nanmedian(pixel_loss, axis=1)
    raise ValueError(f"unknown aggregation {mode!r}")


def interpolate_curve_minimum(axis: np.ndarray, curve: np.ndarray, best: int) -> float:
    """Fit a local quadratic on the non-uniform AOD axis when it is convex."""

    node = float(axis[best])
    if best <= 0 or best >= axis.size - 1:
        return node
    x = np.asarray(axis[best - 1 : best + 2], dtype=np.float64)
    y = np.asarray(curve[best - 1 : best + 2], dtype=np.float64)
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        return node
    a, b, _ = np.polyfit(x, y, 2)
    if not math.isfinite(a) or not math.isfinite(b) or a <= 0.0:
        return node
    vertex = -b / (2.0 * a)
    return float(vertex) if x[0] <= vertex <= x[-1] else node


def _load_results(directory: Path) -> dict[str, dict[str, Any]]:
    records = {}
    for path in directory.glob("*.json"):
        try:
            row = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        matchup_id = str(row.get("matchup_id") or "")
        if matchup_id:
            records[matchup_id] = row
    return records


def _fold(site: str, folds: int) -> int:
    digest = hashlib.sha256(site.encode()).digest()
    return int.from_bytes(digest[:4], "big") % folds


def _aod_bin(value: float) -> str:
    return next(label for lower, upper, label in AOD_BINS if lower <= value < upper)


def _within(retrieved: float, truth: float) -> bool:
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


def _summarize(rows: list[dict[str, Any]], campaign_count: int, folds: int) -> dict[str, Any]:
    errors = np.asarray([row["retrieved"] - row["truth"] for row in rows], dtype=np.float64)
    by_bin = []
    for _, _, label in AOD_BINS:
        selected = [row for row in rows if row["aod_bin"] == label]
        by_bin.append(
            {
                "aod_bin": label,
                "hits": sum(row["within_ee"] for row in selected),
                "count": len(selected),
                "pct": 100.0 * sum(row["within_ee"] for row in selected) / len(selected)
                if selected
                else None,
            }
        )
    by_fold = []
    for fold in range(folds):
        selected = [row for row in rows if row["fold"] == fold]
        by_fold.append(
            {
                "fold": fold,
                "hits": sum(row["within_ee"] for row in selected),
                "count": len(selected),
                "pct": 100.0 * sum(row["within_ee"] for row in selected) / len(selected)
                if selected
                else None,
            }
        )
    hits = sum(row["within_ee"] for row in rows)
    return {
        "hits": hits,
        "available": len(rows),
        "raw_pct": 100.0 * hits / campaign_count,
        "available_pct": 100.0 * hits / len(rows),
        "rmse": float(np.sqrt(np.mean(np.square(errors)))),
        "bias": float(np.mean(errors)),
        "edge_minima": sum(row["minimum_at_edge"] for row in rows),
        "aod_bins": by_bin,
        "site_folds": by_fold,
    }


def evaluate(args: argparse.Namespace) -> dict[str, Any]:
    records = _load_results(args.results_dir)
    model_errors = [float(value) for value in args.model_errors.split(",")]
    base_candidates = [
        (loss, aggregation, model_error)
        for model_error in model_errors
        for loss in ("chi2", "huber1p5", "student3")
        for aggregation in ("mean", "median")
    ]
    candidates = [
        Candidate(loss, aggregation, model_error, interpolate)
        for loss, aggregation, model_error in base_candidates
        for interpolate in (False, True)
    ]
    rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    missing = []
    for index, (matchup_id, record) in enumerate(sorted(records.items()), start=1):
        if str(record.get("status")).upper() != "OK":
            continue
        cube_path = args.cubes_dir / f"{matchup_id}.npz"
        if not cube_path.exists():
            missing.append(matchup_id)
            continue
        with np.load(cube_path, allow_pickle=False) as archive:
            residual = np.asarray(archive["band_residual_cube"], dtype=np.float64)
            uncertainty = np.asarray(archive["boa_unc"], dtype=np.float64)
            solve_valid = np.asarray(archive["solve_valid"], dtype=bool).reshape(-1)
            axis = np.asarray(archive["aot_axis"], dtype=np.float64)
        if residual.ndim != 4 or uncertainty.shape != residual.shape[:1] + residual.shape[2:]:
            raise ValueError(f"unexpected residual/uncertainty shape in {cube_path}")
        residual = residual.reshape(residual.shape[0], residual.shape[1], -1)
        truth = float(record["truth"])
        site = str(record.get("site") or matchup_id.split("__")[0])
        fold = _fold(site, args.folds)
        for model_error in model_errors:
            sigma = np.sqrt(np.square(uncertainty) + model_error**2).reshape(
                uncertainty.shape[0], 1, -1
            )
            standardized = residual / np.maximum(sigma, 1e-6)
            finite = solve_valid & np.all(np.isfinite(standardized), axis=(0, 1))
            if int(np.count_nonzero(finite)) < args.min_pixels:
                finite = solve_valid
            for loss_name in ("chi2", "huber1p5", "student3"):
                loss = robust_loss(standardized, loss_name)
                pixel_loss = np.sum(loss, axis=0)
                pixel_loss[:, ~finite] = np.nan
                for aggregation in ("mean", "median"):
                    curve = aggregate_curve(pixel_loss, aggregation)
                    valid_nodes = np.flatnonzero(np.isfinite(curve))
                    if not valid_nodes.size:
                        continue
                    best = int(valid_nodes[np.argmin(curve[valid_nodes])])
                    for interpolate in (False, True):
                        candidate = Candidate(loss_name, aggregation, model_error, interpolate)
                        retrieved = (
                            interpolate_curve_minimum(axis, curve, best)
                            if interpolate
                            else float(axis[best])
                        )
                        rows[candidate.name].append(
                            {
                                "matchup_id": matchup_id,
                                "site": site,
                                "fold": fold,
                                "truth": truth,
                                "retrieved": retrieved,
                                "within_ee": _within(retrieved, truth),
                                "aod_bin": _aod_bin(truth),
                                "minimum_at_edge": best in {0, axis.size - 1},
                            }
                        )
        if index % 25 == 0:
            print(f"processed {index}/{len(records)} records", flush=True)

    summaries = []
    for candidate in candidates:
        candidate_rows = rows[candidate.name]
        summary = _summarize(candidate_rows, args.campaign_count, args.folds)
        summaries.append(
            {
                "candidate": candidate.name,
                "loss": candidate.loss,
                "aggregation": candidate.aggregation,
                "model_error": candidate.model_error,
                "interpolate": candidate.interpolate,
                **summary,
            }
        )
    summaries.sort(key=lambda row: (-row["hits"], row["rmse"], abs(row["bias"])))
    best_name = summaries[0]["candidate"]
    return {
        "description": "Surface-only robust likelihood screen; no aerosol-product backstop",
        "results_dir": str(args.results_dir),
        "cubes_dir": str(args.cubes_dir),
        "campaign_count": args.campaign_count,
        "available_records": len(rows[best_name]),
        "missing_cubes": missing,
        "signed_covariance_available": False,
        "candidates": summaries,
        "best_candidate": summaries[0],
        "best_rows": rows[best_name],
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=RESULTS)
    parser.add_argument("--cubes-dir", type=Path, default=CUBES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--model-errors", default="0,0.010,0.015")
    parser.add_argument("--min-pixels", type=int, default=100)
    parser.add_argument("--campaign-count", type=int, default=250)
    parser.add_argument("--folds", type=int, default=5)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    result = evaluate(args)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result["best_candidate"], indent=2), flush=True)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
