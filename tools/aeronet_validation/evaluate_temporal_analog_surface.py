"""Evaluate an AERONET-free same-pixel temporal-analogue surface prior.

For each candidate AOD, current-scene NIR/SWIR TOA is corrected to BOA using
the saved RT coefficient curves. Historical clear-surface realizations at the
same pixel are ranked by their NIR/SWIR distance to that corrected anchor. The
visible prior is the median of the closest realizations, and its uncertainty is
their robust spread. No AERONET value is used by the retrieval.

The script scores only the locked development partition unless
``--unlock-holdout`` is explicitly supplied.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from analyze_medium_aod_physics import (
    DEFAULT_BASELINE,
    DEFAULT_EXCLUDED,
    DEFAULT_HOLDOUT_FOLDS,
    DEFAULT_MIDS,
    DEFAULT_SPLIT_SEED,
    _load_results,
    _site_group_folds,
    within_ee,
)

DEFAULT_DUMPS = Path(
    "/gws/ssde/j25a/nceo_isp/public/siac_refactor/calib_dumps_c250"
)
DEFAULT_OUTPUT = Path(
    "/gws/ssde/j25a/nceo_isp/public/siac_refactor/"
    "analysis/medium_aod_temporal_analog_development_20260713"
)
ANCHOR_COLUMNS = (4, 5, 6)
VISIBLE_COLUMNS = {"B02": 1, "B04": 3}
MAD_TO_STD = 1.4826


@dataclass(frozen=True)
class AnalogConfig:
    neighbors: int
    uncertainty_floor: float
    spatial_stride: int


def _correct(toa: np.ndarray, xap: float, xbp: float, xcp: float) -> np.ndarray:
    y = float(xap) * toa - float(xbp)
    denominator = 1.0 + float(xcp) * y
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(np.abs(denominator) > 1.0e-12, y / denominator, np.nan)


def _calibrated_backstop_sigma(aot_prior: float) -> float:
    loose = max(0.5 * aot_prior, 0.02)
    mid = max(0.07, 0.5 * aot_prior / (1.0 + math.exp(-(aot_prior - 0.5) / 0.15)))
    return loose if aot_prior < 0.15 else mid


def _nearest_analogue_prior(
    composites: np.ndarray,
    target_anchor: np.ndarray,
    *,
    neighbors: int,
    uncertainty_floor: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    historical_anchor = np.moveaxis(composites[:, ANCHOR_COLUMNS], 1, -1)
    finite_anchor = np.all(np.isfinite(historical_anchor), axis=-1)
    with np.errstate(invalid="ignore"):
        distance2 = np.mean(
            np.square(historical_anchor - target_anchor[np.newaxis]),
            axis=-1,
        )
    distance2 = np.where(finite_anchor, distance2, np.inf)
    k = min(max(int(neighbors), 1), composites.shape[0])
    nearest = np.argpartition(distance2, kth=k - 1, axis=0)[:k]
    finite_count = np.count_nonzero(np.isfinite(distance2), axis=0)

    priors: list[np.ndarray] = []
    uncertainties: list[np.ndarray] = []
    for column in VISIBLE_COLUMNS.values():
        values = np.take_along_axis(composites[:, column], nearest, axis=0)
        with np.errstate(invalid="ignore"):
            centre = np.nanmedian(values, axis=0)
            spread = np.nanmedian(np.abs(values - centre[np.newaxis]), axis=0) * MAD_TO_STD
        priors.append(centre)
        uncertainties.append(np.maximum(spread, float(uncertainty_floor)))
    return np.stack(priors), np.stack(uncertainties), finite_count


def _evaluate_case(
    dump_path: Path,
    record: dict[str, Any],
    config: AnalogConfig,
) -> dict[str, Any]:
    with np.load(dump_path, allow_pickle=False) as source:
        aot_axis = np.asarray(source["aot_grid"], dtype=np.float64)
        composites = np.asarray(source["comp"], dtype=np.float64)[
            :, :, :: config.spatial_stride, :: config.spatial_stride
        ]
        toa = {
            name: np.asarray(source[f"toa_{name}"], dtype=np.float64)[
                :: config.spatial_stride, :: config.spatial_stride
            ]
            for name in (*VISIBLE_COLUMNS, "B8A", "B11", "B12")
        }
        coefficients = {
            name: tuple(
                np.asarray(source[f"{term}_{name}"], dtype=np.float64)
                for term in ("xap", "xbp", "xcp")
            )
            for name in (*VISIBLE_COLUMNS, "B8A", "B11", "B12")
        }

    node_cost = np.full(aot_axis.shape, np.nan, dtype=np.float64)
    node_support = np.zeros(aot_axis.shape, dtype=np.int64)
    for node in range(aot_axis.size):
        target_anchor = np.stack(
            [
                _correct(toa[name], *(values[node] for values in coefficients[name]))
                for name in ("B8A", "B11", "B12")
            ],
            axis=-1,
        )
        prior, uncertainty, analogue_count = _nearest_analogue_prior(
            composites,
            target_anchor,
            neighbors=config.neighbors,
            uncertainty_floor=config.uncertainty_floor,
        )
        corrected_visible = np.stack(
            [
                _correct(toa[name], *(values[node] for values in coefficients[name]))
                for name in VISIBLE_COLUMNS
            ]
        )
        with np.errstate(invalid="ignore", divide="ignore"):
            band_cost = np.square((corrected_visible - prior) / uncertainty)
        valid = (
            np.all(np.isfinite(band_cost), axis=0)
            & np.all(np.isfinite(target_anchor), axis=-1)
            & (analogue_count >= config.neighbors)
        )
        if np.count_nonzero(valid) >= 50:
            node_cost[node] = float(np.nanmedian(np.sum(band_cost[:, valid], axis=0)))
            node_support[node] = int(np.count_nonzero(valid))

    aot_prior = float(record["atmo_prior"]["aot_median"])
    prior_sigma = _calibrated_backstop_sigma(aot_prior)
    total_cost = node_cost + np.square((aot_axis - aot_prior) / prior_sigma)
    finite = np.isfinite(total_cost)
    if not np.any(finite):
        raise ValueError("no AOD node has sufficient temporal-analogue support")
    best = int(np.nanargmin(total_cost))
    surface_best = int(np.nanargmin(node_cost))
    return {
        "retrieved": float(aot_axis[best]),
        "surface_min_aot": float(aot_axis[surface_best]),
        "surface_min_cost": float(node_cost[surface_best]),
        "aot_prior": aot_prior,
        "aot_prior_sigma": prior_sigma,
        "best_support": int(node_support[best]),
        "aot_axis": aot_axis.tolist(),
        "surface_cost": [float(value) if math.isfinite(value) else None for value in node_cost],
        "total_cost": [float(value) if math.isfinite(value) else None for value in total_cost],
        "node_support": node_support.tolist(),
    }


def _summary(rows: list[dict[str, Any]], subset: str) -> dict[str, Any]:
    selected = [row for row in rows if row["subset"] == subset]
    errors = np.asarray([row["retrieved"] - row["truth"] for row in selected])
    hits = sum(bool(row["within_ee"]) for row in selected)
    return {
        "n": len(selected),
        "hits": hits,
        "within_ee_pct": 100.0 * hits / len(selected) if selected else None,
        "bias": float(np.mean(errors)) if errors.size else None,
        "mae": float(np.mean(np.abs(errors))) if errors.size else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--mids", type=Path, default=DEFAULT_MIDS)
    parser.add_argument("--dumps", type=Path, default=DEFAULT_DUMPS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--neighbors", type=int, default=3)
    parser.add_argument("--uncertainty-floor", type=float, default=0.006)
    parser.add_argument("--spatial-stride", type=int, default=2)
    parser.add_argument("--all-regimes", action="store_true")
    parser.add_argument("--unlock-holdout", action="store_true")
    args = parser.parse_args()

    if args.neighbors < 1:
        raise SystemExit("--neighbors must be positive")
    if args.uncertainty_floor <= 0.0:
        raise SystemExit("--uncertainty-floor must be positive")
    if args.spatial_stride < 1:
        raise SystemExit("--spatial-stride must be positive")

    mids = [line.strip() for line in args.mids.read_text().splitlines() if line.strip()]
    baseline, statuses = _load_results(args.baseline, mids)
    if len(baseline) != len(mids):
        raise SystemExit(f"baseline is incomplete: {len(baseline)}/{len(mids)} ({statuses})")
    folds = _site_group_folds(baseline, mids, seed=DEFAULT_SPLIT_SEED)
    selected = [
        matchup_id
        for matchup_id in mids
        if args.unlock_holdout or folds[matchup_id] not in set(DEFAULT_HOLDOUT_FOLDS)
    ]
    if not args.all_regimes:
        selected = [
            matchup_id
            for matchup_id in selected
            if 0.25 <= float(baseline[matchup_id]["truth"]) <= 0.85
        ]

    config = AnalogConfig(
        neighbors=args.neighbors,
        uncertainty_floor=args.uncertainty_floor,
        spatial_stride=args.spatial_stride,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for index, matchup_id in enumerate(selected, start=1):
        dump_path = args.dumps / f"{matchup_id}.npz"
        if not dump_path.exists():
            missing.append(matchup_id)
            continue
        record = baseline[matchup_id]
        result = _evaluate_case(dump_path, record, config)
        truth = float(record["truth"])
        retrieved = float(result["retrieved"])
        if truth < 0.25:
            subset = "low"
        elif truth <= 0.85:
            subset = "medium"
        elif matchup_id in DEFAULT_EXCLUDED:
            subset = "excluded_extreme"
        else:
            subset = "high"
        row = {
            "matchup_id": matchup_id,
            "site": record.get("site"),
            "fold": folds[matchup_id],
            "subset": subset,
            "truth": truth,
            "baseline": float(record["retrieved"]),
            "retrieved": retrieved,
            "error": retrieved - truth,
            "within_ee": within_ee(retrieved, truth),
            "baseline_within_ee": within_ee(float(record["retrieved"]), truth),
            **result,
        }
        rows.append(row)
        (args.output_dir / f"{matchup_id}.json").write_text(
            json.dumps(row, indent=2, allow_nan=False) + "\n"
        )
        print(
            f"[{index:02d}/{len(selected):02d}] {matchup_id}: "
            f"truth={truth:.3f} baseline={row['baseline']:.3f} analog={retrieved:.3f} "
            f"hit={row['within_ee']}"
        )

    fieldnames = [
        "matchup_id",
        "site",
        "fold",
        "subset",
        "truth",
        "baseline",
        "retrieved",
        "error",
        "within_ee",
        "baseline_within_ee",
        "surface_min_aot",
        "surface_min_cost",
        "aot_prior",
        "aot_prior_sigma",
        "best_support",
    ]
    with (args.output_dir / "cases.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "scope": "all_unlocked" if args.unlock_holdout else "development_locked",
        "config": vars(config),
        "requested": len(selected),
        "completed": len(rows),
        "missing_dumps": missing,
        "subsets": {name: _summary(rows, name) for name in ("low", "medium", "high")},
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    medium = summary["subsets"]["medium"]
    print(
        f"medium: {medium['hits']}/{medium['n']} "
        f"({medium['within_ee_pct']:.2f}%), bias={medium['bias']:+.4f}; "
        f"missing={len(missing)}"
    )


if __name__ == "__main__":
    main()
