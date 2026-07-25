"""Analyze surface-solver cost and cloud relationships with AOD error."""

from __future__ import annotations

import argparse
import json
import math
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
FULL_R2_DIR = ROOT / "phaseD_results_campaign250_R2_full"
BAD_DIAG_DIR = ROOT / "phaseD_results_campaign250_surface_bad89_R2_diag_20260705"

COST_BINS = (
    ("<1", 0.0, 1.0),
    ("1-2", 1.0, 2.0),
    ("2-5", 2.0, 5.0),
    ("5-20", 5.0, 20.0),
    ("20-100", 20.0, 100.0),
    (">=100", 100.0, math.inf),
)

CLOUD_BINS = (
    ("0-5%", 0.0, 0.05),
    ("5-25%", 0.05, 0.25),
    ("25-75%", 0.25, 0.75),
    ("75-100%", 0.75, math.inf),
)


def _safe_float(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def _load_records(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for json_path in sorted(path.glob("*.json")):
        with json_path.open("r", encoding="utf-8") as handle:
            record = json.load(handle)
        record["_path"] = str(json_path)
        records.append(record)
    return records


def _solver(record: dict[str, Any]) -> dict[str, Any]:
    solver = record.get("solver")
    return solver if isinstance(solver, dict) else {}


def _ok_rows(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        if str(record.get("status")).upper() != "OK":
            continue
        truth = _safe_float(record.get("truth"))
        retrieved = _safe_float(record.get("retrieved"))
        if not (math.isfinite(truth) and math.isfinite(retrieved)):
            continue
        solver = _solver(record)
        err = retrieved - truth
        cloud = _safe_float(record.get("cloud_frac"))
        rows.append(
            {
                "mid": str(record.get("matchup_id") or Path(str(record.get("_path"))).stem),
                "truth": truth,
                "retrieved": retrieved,
                "err": err,
                "abs_err": abs(err),
                "within": bool(record.get("within_ee")),
                "cloud_frac": cloud,
                "cloud_at_station": _safe_float(record.get("cloud_at_station")),
                "cost": _safe_float(solver.get("cost_final_per_band")),
                "aot_std": _safe_float(solver.get("aot_std")),
                "curve_delta": _safe_float(
                    solver.get("surface_cost_curve_relative_second_delta")
                ),
                "band_spread": _safe_float(solver.get("surface_band_argmin_spread")),
                "curve_min": _safe_float(solver.get("surface_cost_curve_min_aot")),
                "final_node": _safe_float(solver.get("surface_final_aot_median")),
                "tau_gate": bool(solver.get("surface_tau_gate_fired")),
            }
        )
    return rows


def _values(rows: list[dict[str, Any]], field: str) -> np.ndarray:
    return np.asarray([_safe_float(row.get(field)) for row in rows], dtype=np.float64)


def _finite_pair(rows: list[dict[str, Any]], x_field: str, y_field: str) -> tuple[np.ndarray, np.ndarray]:
    x = _values(rows, x_field)
    y = _values(rows, y_field)
    mask = np.isfinite(x) & np.isfinite(y)
    return x[mask], y[mask]


def _rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=np.float64)
    i = 0
    while i < values.size:
        j = i + 1
        while j < values.size and values[order[j]] == values[order[i]]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + j - 1) + 1.0
        i = j
    return ranks


def _corr(x: np.ndarray, y: np.ndarray, *, spearman: bool = False) -> float:
    if x.size < 3 or y.size < 3:
        return math.nan
    if spearman:
        x = _rankdata(x)
        y = _rankdata(y)
    if np.nanstd(x) == 0.0 or np.nanstd(y) == 0.0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def _partial_corr(
    rows: list[dict[str, Any]],
    x_field: str,
    y_field: str,
    control_fields: tuple[str, ...],
) -> tuple[int, float]:
    fields = (x_field, y_field, *control_fields)
    matrix = np.vstack([_values(rows, field) for field in fields]).T
    mask = np.all(np.isfinite(matrix), axis=1)
    matrix = matrix[mask]
    if matrix.shape[0] < len(control_fields) + 4:
        return int(matrix.shape[0]), math.nan
    x = matrix[:, 0]
    y = matrix[:, 1]
    controls = matrix[:, 2:]
    design = np.column_stack([np.ones(matrix.shape[0]), controls])
    x_resid = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    y_resid = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return int(matrix.shape[0]), _corr(x_resid, y_resid)


def _median(values: Iterable[float]) -> float:
    arr = np.asarray([value for value in values if math.isfinite(value)], dtype=np.float64)
    if arr.size == 0:
        return math.nan
    return float(np.nanmedian(arr))


def _pct(numerator: int, denominator: int) -> float:
    return numerator / denominator * 100.0 if denominator else 0.0


def _print_correlations(title: str, rows: list[dict[str, Any]], fields: tuple[str, ...]) -> None:
    print(f"\n{title}")
    failure = np.asarray([not bool(row["within"]) for row in rows], dtype=np.float64)
    for field in fields:
        x = _values(rows, field)
        mask = np.isfinite(x)
        if not np.any(mask):
            continue
        abs_err = _values(rows, "abs_err")
        err = _values(rows, "err")
        abs_mask = mask & np.isfinite(abs_err)
        err_mask = mask & np.isfinite(err)
        fail_mask = mask & np.isfinite(failure)
        print(
            f"  {field:16s} n={int(mask.sum()):3d} "
            f"pearson(abs_err)={_corr(x[abs_mask], abs_err[abs_mask]):+0.3f} "
            f"spearman(abs_err)={_corr(x[abs_mask], abs_err[abs_mask], spearman=True):+0.3f} "
            f"pearson(err)={_corr(x[err_mask], err[err_mask]):+0.3f} "
            f"pearson(fail)={_corr(x[fail_mask], failure[fail_mask]):+0.3f}"
        )


def _bin_name(value: float, bins: tuple[tuple[str, float, float], ...]) -> str:
    if not math.isfinite(value):
        return "nan"
    for name, lo, hi in bins:
        if lo <= value < hi:
            return name
    return "nan"


def _print_bin_table(
    title: str,
    rows: list[dict[str, Any]],
    *,
    field: str,
    bins: tuple[tuple[str, float, float], ...],
) -> None:
    print(f"\n{title}")
    print("  bin       n  within%  median_err  median_abs_err  under  over")
    order = [name for name, _lo, _hi in bins] + ["nan"]
    for name in order:
        subset = [row for row in rows if _bin_name(_safe_float(row.get(field)), bins) == name]
        if not subset:
            continue
        within = sum(bool(row["within"]) for row in subset)
        under = sum(float(row["err"]) < 0.0 for row in subset)
        over = sum(float(row["err"]) > 0.0 for row in subset)
        print(
            f"  {name:8s} {len(subset):3d} {_pct(within, len(subset)):7.1f} "
            f"{_median(float(row['err']) for row in subset):+11.3f} "
            f"{_median(float(row['abs_err']) for row in subset):14.3f} "
            f"{under:6d} {over:5d}"
        )


def _print_threshold_table(
    title: str,
    rows: list[dict[str, Any]],
    *,
    field: str,
    thresholds: tuple[float, ...],
) -> None:
    total_fail = sum(not bool(row["within"]) for row in rows)
    print(f"\n{title}")
    print("  threshold      flagged  failures  precision  recall")
    for threshold in thresholds:
        flagged = [
            row
            for row in rows
            if math.isfinite(_safe_float(row.get(field))) and _safe_float(row.get(field)) >= threshold
        ]
        failures = sum(not bool(row["within"]) for row in flagged)
        print(
            f"  {field}>={threshold:<6g} {len(flagged):8d} {failures:9d} "
            f"{_pct(failures, len(flagged)):9.1f} {_pct(failures, total_fail):7.1f}"
        )


def _print_bad_curve_summary(rows: list[dict[str, Any]]) -> None:
    rows = [row for row in rows if math.isfinite(float(row["curve_min"]))]
    if not rows:
        return
    curve_better = sum(
        abs(float(row["curve_min"]) - float(row["truth"])) < float(row["abs_err"])
        for row in rows
    )
    curve_within = sum(
        abs(float(row["curve_min"]) - float(row["truth"]))
        <= 0.05 + 0.15 * float(row["truth"])
        for row in rows
    )
    final_within = sum(
        math.isfinite(float(row["final_node"]))
        and abs(float(row["final_node"]) - float(row["truth"]))
        <= 0.05 + 0.15 * float(row["truth"])
        for row in rows
    )
    flat = sum(math.isfinite(float(row["curve_delta"])) and float(row["curve_delta"]) <= 0.05 for row in rows)
    spread = sum(math.isfinite(float(row["band_spread"])) and float(row["band_spread"]) >= 0.2 for row in rows)
    print("\nBad-suite cost-curve checks")
    print(f"  curve min closer than retrieved: {curve_better}/{len(rows)}")
    print(f"  curve min within EE: {curve_within}/{len(rows)}")
    print(f"  final node within EE: {final_within}/{len(rows)}")
    print(f"  relative second delta <= 5%: {flat}/{len(rows)}")
    print(f"  B02/B04 argmin spread >= 0.2: {spread}/{len(rows)}")


def _print_cloud_control_checks(rows: list[dict[str, Any]]) -> None:
    cloud, cost = _finite_pair(rows, "cloud_frac", "cost")
    n_abs, partial_abs = _partial_corr(rows, "cloud_frac", "abs_err", ("cost",))
    n_signed, partial_signed = _partial_corr(rows, "cloud_frac", "err", ("cost",))
    print("\nCloud/cost control checks")
    print(
        f"  cloud_frac vs cost: n={cloud.size} "
        f"pearson={_corr(cloud, cost):+0.3f} "
        f"spearman={_corr(cloud, cost, spearman=True):+0.3f}"
    )
    print(
        f"  cloud_frac vs abs_err controlling cost: "
        f"n={n_abs} partial_pearson={partial_abs:+0.3f}"
    )
    print(
        f"  cloud_frac vs signed err controlling cost: "
        f"n={n_signed} partial_pearson={partial_signed:+0.3f}"
    )


def analyze_dataset(label: str, rows: list[dict[str, Any]], *, include_curve: bool) -> None:
    total = len(rows)
    within = sum(bool(row["within"]) for row in rows)
    print(f"\n=== {label} ===")
    print(f"OK rows: {total}; within EE: {within}/{total} = {_pct(within, total):.1f}%")
    print(
        "Error summary: "
        f"median_err={_median(float(row['err']) for row in rows):+.3f}; "
        f"median_abs_err={_median(float(row['abs_err']) for row in rows):.3f}; "
        f"under={sum(float(row['err']) < 0.0 for row in rows)}; "
        f"over={sum(float(row['err']) > 0.0 for row in rows)}"
    )
    _print_correlations(
        "Correlations",
        rows,
        (
            "cost",
            "aot_std",
            "cloud_frac",
            "cloud_at_station",
            "curve_delta",
            "band_spread",
        )
        if include_curve
        else ("cost", "aot_std", "cloud_frac", "cloud_at_station"),
    )
    _print_bin_table("Cost bins", rows, field="cost", bins=COST_BINS)
    _print_threshold_table(
        "Cost thresholds as failure detector",
        rows,
        field="cost",
        thresholds=(1.0, 2.0, 5.0, 20.0, 100.0),
    )
    _print_bin_table("Cloud-fraction bins", rows, field="cloud_frac", bins=CLOUD_BINS)
    _print_threshold_table(
        "Cloud thresholds as failure detector",
        rows,
        field="cloud_frac",
        thresholds=(0.25, 0.5, 0.75, 0.95),
    )
    _print_cloud_control_checks(rows)
    if include_curve:
        _print_bad_curve_summary(rows)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze cost/cloud relationships with AERONET AOD error."
    )
    parser.add_argument("--full-r2-dir", type=Path, default=FULL_R2_DIR)
    parser.add_argument("--bad-diag-dir", type=Path, default=BAD_DIAG_DIR)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    full_rows = _ok_rows(_load_records(args.full_r2_dir))
    bad_rows = _ok_rows(_load_records(args.bad_diag_dir))
    analyze_dataset("Full R2 campaign", full_rows, include_curve=False)
    analyze_dataset("Bad89 diagnostic suite", bad_rows, include_curve=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
