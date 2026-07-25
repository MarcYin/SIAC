"""Analyze local station-window cost-shape diagnostics from dumped cost cubes."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from pyproj import Transformer

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.aeronet_validation.analyze_surface_driven_failures import (  # noqa: E402
    DEFAULT_ROOT,
    SURFACE_RESULT_DIRS,
    _action_bucket,
    _load_dir,
    _row_for_mid,
)

from siac.algorithms.solver.surface_driven import integrated_cost_field_aod_from_npz  # noqa: E402

DEFAULT_RESULTS_DIR = (
    DEFAULT_ROOT / "phaseD_results_campaign250_surface_bad89_R2_localdiag_20260705"
)
DEFAULT_COST_DIR = (
    DEFAULT_ROOT / "phaseD_cost_cubes_campaign250_surface_bad89_R2_localdiag_20260705"
)
DEFAULT_MIDS_FILE = DEFAULT_ROOT / "campaign250_surface_bad89_mids.txt"


def _safe_float(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _read_mids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _result_records(path: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for json_path in sorted(path.glob("*.json")):
        record = _load_json(json_path)
        records[str(record.get("matchup_id") or json_path.stem)] = record
    return records


def _cost_cube_paths(cost_dir: Path, mid: str) -> list[Path]:
    base = cost_dir / f"{mid}.npz"
    if base.exists():
        return [base]
    return sorted(cost_dir.glob(f"{mid}.species*.npz"))


def _candidate_cost(candidate: dict[str, Any]) -> float:
    cost = _safe_float(candidate.get("cost"))
    return cost if math.isfinite(cost) else math.inf


def _resolve_local_candidate(
    *,
    record: dict[str, Any],
    cost_paths: list[Path],
    radius_m: float,
) -> dict[str, Any]:
    crs = record.get("scene_crs")
    lon = _safe_float(record.get("lon"))
    lat = _safe_float(record.get("lat"))
    if not crs or not (math.isfinite(lon) and math.isfinite(lat)):
        raise ValueError("record is missing scene_crs/lon/lat")
    transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    sx, sy = transformer.transform(lon, lat)

    candidates: list[dict[str, Any]] = []
    for path in cost_paths:
        local = integrated_cost_field_aod_from_npz(
            str(path),
            center_x=float(sx),
            center_y=float(sy),
            radius_m=radius_m,
        )
        candidates.append({"path": str(path), **local})
    if not candidates:
        raise ValueError("no local candidates")
    selected = min(candidates, key=_candidate_cost)
    return {
        "selected_index": int(candidates.index(selected)),
        "candidate_count": int(len(candidates)),
        "selected_path": selected["path"],
        **selected,
    }


def _within(aod: float, truth: float) -> bool:
    return math.isfinite(aod) and math.isfinite(truth) and abs(aod - truth) <= 0.05 + 0.15 * truth


def _median(values: list[float]) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    mid = len(finite) // 2
    if len(finite) % 2:
        return finite[mid]
    return 0.5 * (finite[mid - 1] + finite[mid])


def aggregation_consensus_aod(
    retrieved: float,
    band_argmins: list[float],
    *,
    spread_max: float = 0.10,
    deviation_min: float = 0.10,
) -> tuple[float, bool]:
    """Use band consensus when agreeing band curves contradict aggregation."""

    finite = [float(value) for value in band_argmins if math.isfinite(float(value))]
    if len(finite) < 2 or not math.isfinite(retrieved):
        return retrieved, False
    consensus = _median(finite)
    spread = max(finite) - min(finite)
    use = spread <= spread_max and abs(retrieved - consensus) >= deviation_min
    return (consensus if use else retrieved), use


def _corr(x_values: list[float], y_values: list[float]) -> float:
    x = np.asarray(x_values, dtype=np.float64)
    y = np.asarray(y_values, dtype=np.float64)
    mask = np.isfinite(x) & np.isfinite(y)
    if int(np.count_nonzero(mask)) < 3:
        return math.nan
    x = x[mask]
    y = y[mask]
    if float(np.nanstd(x)) == 0.0 or float(np.nanstd(y)) == 0.0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def _case_group(row: dict[str, Any]) -> str:
    if row["within_r2"]:
        return "R R2 within EE"
    if not row["surface_oracle"]:
        return _action_bucket(row)
    return "G R2 miss recoverable by another saved surface-only arm"


def analyze(args: argparse.Namespace) -> int:
    root = Path(args.root)
    mids = _read_mids(Path(args.mids_file))
    records = _result_records(Path(args.results_dir))
    surface_records = {
        label: _load_dir(root, dirname) for label, dirname in SURFACE_RESULT_DIRS.items()
    }

    rows: list[dict[str, Any]] = []
    missing_results: list[str] = []
    missing_cubes: list[str] = []
    failures: list[str] = []
    for mid in mids:
        record = records.get(mid)
        if record is None:
            missing_results.append(mid)
            continue
        if str(record.get("status")).upper() != "OK":
            failures.append(mid)
            continue
        paths = _cost_cube_paths(Path(args.cost_dir), mid)
        if not paths:
            missing_cubes.append(mid)
            continue
        truth = _safe_float(record.get("truth"))
        retrieved = _safe_float(record.get("retrieved"))
        local = _resolve_local_candidate(
            record=record,
            cost_paths=paths,
            radius_m=float(args.radius_m),
        )
        diagnostics = local.get("diagnostics") if isinstance(local.get("diagnostics"), dict) else {}
        local_aod = _safe_float(local.get("aod"))
        local_curve_min = _safe_float(diagnostics.get("local_cost_curve_min_aot"))
        old_row = _row_for_mid(mid, surface_records)
        row = {
            "mid": mid,
            "group": _case_group(old_row),
            "truth": truth,
            "retrieved": retrieved,
            "within": bool(record.get("within_ee")),
            "err": retrieved - truth,
            "abs_err": abs(retrieved - truth),
            "local_aod": local_aod,
            "local_err": local_aod - truth,
            "local_within": _within(local_aod, truth),
            "local_curve_min": local_curve_min,
            "local_curve_min_within": _within(local_curve_min, truth),
            "local_curve_delta": _safe_float(
                diagnostics.get("local_cost_curve_relative_second_delta")
            ),
            "local_curve_min_at_edge": bool(diagnostics.get("local_cost_curve_min_at_edge")),
            "local_band_spread": _safe_float(diagnostics.get("local_band_argmin_spread")),
            "local_b02_argmin": _safe_float(diagnostics.get("local_band_B02_argmin_aot")),
            "local_b04_argmin": _safe_float(diagnostics.get("local_band_B04_argmin_aot")),
            "local_cost": _safe_float(local.get("cost")),
            "selected_pass": local.get("selected_pass"),
            "candidate_count": int(local.get("candidate_count", 1)),
        }
        consensus_aod, consensus_used = aggregation_consensus_aod(
            retrieved,
            [float(row["local_b02_argmin"]), float(row["local_b04_argmin"])],
        )
        row.update(
            aggregation_consensus_aod=consensus_aod,
            aggregation_consensus_used=consensus_used,
            aggregation_consensus_within=_within(consensus_aod, truth),
        )
        rows.append(row)

    print("Local cost diagnostic coverage")
    print(f"  expected mids: {len(mids)}")
    print(f"  result records: {len(records)}")
    print(f"  analyzed OK rows with cost cubes: {len(rows)}")
    print(f"  missing results: {len(missing_results)}")
    print(f"  non-OK result records: {len(failures)}")
    print(f"  missing cost cubes: {len(missing_cubes)}")

    if not rows:
        return 2

    local_within = sum(bool(row["local_within"]) for row in rows)
    curve_min_within = sum(bool(row["local_curve_min_within"]) for row in rows)
    flat = sum(
        math.isfinite(float(row["local_curve_delta"]))
        and float(row["local_curve_delta"]) <= float(args.flat_delta)
        for row in rows
    )
    spread = sum(
        math.isfinite(float(row["local_band_spread"]))
        and float(row["local_band_spread"]) >= float(args.spread_threshold)
        for row in rows
    )
    edge = sum(bool(row["local_curve_min_at_edge"]) for row in rows)
    improved = sum(abs(float(row["local_err"])) < float(row["abs_err"]) for row in rows)
    print("\nLocal cost-shape summary")
    print(f"  local integrated AOD within EE: {local_within}/{len(rows)}")
    print(f"  local curve-min within EE: {curve_min_within}/{len(rows)}")
    print(f"  local AOD closer than reported R2 retrieved: {improved}/{len(rows)}")
    print(f"  local curve min at AOD-axis edge: {edge}/{len(rows)}")
    print(f"  local relative second delta <= {float(args.flat_delta):.3f}: {flat}/{len(rows)}")
    print(
        f"  local B02/B04 argmin spread >= {float(args.spread_threshold):.3f}: {spread}/{len(rows)}"
    )
    consensus_rows = [row for row in rows if bool(row["aggregation_consensus_used"])]
    consensus_gains = [
        row["mid"]
        for row in consensus_rows
        if not bool(row["within"]) and bool(row["aggregation_consensus_within"])
    ]
    consensus_losses = [
        row["mid"]
        for row in consensus_rows
        if bool(row["within"]) and not bool(row["aggregation_consensus_within"])
    ]
    consensus_hits = sum(bool(row["aggregation_consensus_within"]) for row in rows)
    print(
        "  aggregation-consensus rule: "
        f"{consensus_hits}/{len(rows)} within EE; used={len(consensus_rows)}; "
        f"gains={len(consensus_gains)}; losses={len(consensus_losses)}"
    )
    print(
        "  medians: "
        f"reported_err={_median([float(row['err']) for row in rows]):+.3f}; "
        f"local_err={_median([float(row['local_err']) for row in rows]):+.3f}; "
        f"local_abs_err={_median([abs(float(row['local_err'])) for row in rows]):.3f}; "
        f"local_cost={_median([float(row['local_cost']) for row in rows]):.3f}; "
        f"local_delta={_median([float(row['local_curve_delta']) for row in rows]):.3f}; "
        f"local_spread={_median([float(row['local_band_spread']) for row in rows]):.3f}"
    )
    print("\nLocal diagnostic correlations with reported absolute error")
    print(
        f"  local_cost: {_corr([float(row['local_cost']) for row in rows], [float(row['abs_err']) for row in rows]):+.3f}"
    )
    print(
        f"  local_curve_delta: {_corr([float(row['local_curve_delta']) for row in rows], [float(row['abs_err']) for row in rows]):+.3f}"
    )
    print(
        f"  local_band_spread: {_corr([float(row['local_band_spread']) for row in rows], [float(row['abs_err']) for row in rows]):+.3f}"
    )

    total_failures = sum(not bool(row["within"]) for row in rows)
    print("\nGate diagnostic table against reported R2 failures")
    print(
        "  gate                                      flagged  failures  hits  precision  recall  hit_risk"
    )

    def _print_gate(name: str, flagged: list[dict[str, Any]]) -> None:
        failures_count = sum(not bool(row["within"]) for row in flagged)
        hits = len(flagged) - failures_count
        precision = failures_count / len(flagged) * 100.0 if flagged else 0.0
        recall = failures_count / total_failures * 100.0 if total_failures else 0.0
        hit_risk = hits / len(flagged) * 100.0 if flagged else 0.0
        print(
            f"  {name:40s} {len(flagged):7d} {failures_count:9d} {hits:5d} "
            f"{precision:9.1f} {recall:7.1f} {hit_risk:8.1f}"
        )

    for threshold in (0.1, 0.2, 0.3, 0.5, 1.0):
        _print_gate(
            f"local_spread >= {threshold:g}",
            [
                row
                for row in rows
                if math.isfinite(float(row["local_band_spread"]))
                and float(row["local_band_spread"]) >= threshold
            ],
        )
    for threshold in (0.001, 0.005, 0.01, 0.02, 0.05):
        _print_gate(
            f"local_delta <= {threshold:g}",
            [
                row
                for row in rows
                if math.isfinite(float(row["local_curve_delta"]))
                and float(row["local_curve_delta"]) <= threshold
            ],
        )
    for delta, spread_threshold in ((0.05, 0.2), (0.02, 0.2), (0.01, 0.2), (0.05, 0.5)):
        _print_gate(
            f"delta <= {delta:g} and spread >= {spread_threshold:g}",
            [
                row
                for row in rows
                if math.isfinite(float(row["local_curve_delta"]))
                and math.isfinite(float(row["local_band_spread"]))
                and float(row["local_curve_delta"]) <= delta
                and float(row["local_band_spread"]) >= spread_threshold
            ],
        )
    for cost_threshold, spread_threshold in ((5.0, 0.2), (20.0, 0.2), (5.0, 0.5)):
        _print_gate(
            f"cost >= {cost_threshold:g} and spread >= {spread_threshold:g}",
            [
                row
                for row in rows
                if math.isfinite(float(row["local_cost"]))
                and math.isfinite(float(row["local_band_spread"]))
                and float(row["local_cost"]) >= cost_threshold
                and float(row["local_band_spread"]) >= spread_threshold
            ],
        )

    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["group"])].append(row)
    print("\nBy bad-case group")
    for group, items in sorted(grouped.items()):
        print(
            f"  {group}: n={len(items)} "
            f"local_within={sum(bool(item['local_within']) for item in items)} "
            f"curve_min_within={sum(bool(item['local_curve_min_within']) for item in items)} "
            f"flat={sum(math.isfinite(float(item['local_curve_delta'])) and float(item['local_curve_delta']) <= float(args.flat_delta) for item in items)} "
            f"spread={sum(math.isfinite(float(item['local_band_spread'])) and float(item['local_band_spread']) >= float(args.spread_threshold) for item in items)} "
            f"median_local_err={_median([float(item['local_err']) for item in items]):+.3f}"
        )

    print(f"\nTop {int(args.top)} local opportunities by reported-to-local |error| improvement")
    for row in sorted(
        rows,
        key=lambda item: float(item["abs_err"]) - abs(float(item["local_err"])),
        reverse=True,
    )[: int(args.top)]:
        gain = float(row["abs_err"]) - abs(float(row["local_err"]))
        print(
            f"  {str(row['mid'])[:54]:54s} "
            f"gain={gain:+.3f} reported_err={float(row['err']):+.3f} "
            f"local_err={float(row['local_err']):+.3f} "
            f"local_aod={float(row['local_aod']):.3f} "
            f"curve_min={float(row['local_curve_min']):.3f} "
            f"delta={float(row['local_curve_delta']):.3f} "
            f"spread={float(row['local_band_spread']):.3f} "
            f"group={row['group']}"
        )

    if args.output is not None:
        Path(args.output).write_text(
            json.dumps(
                {
                    "results_dir": str(args.results_dir),
                    "cost_dir": str(args.cost_dir),
                    "radius_m": float(args.radius_m),
                    "available": len(rows),
                    "aggregation_consensus": {
                        "spread_max": 0.10,
                        "deviation_min": 0.10,
                        "hits": consensus_hits,
                        "used": len(consensus_rows),
                        "gains": consensus_gains,
                        "losses": consensus_losses,
                    },
                    "rows": rows,
                    "missing_results": missing_results,
                    "missing_cubes": missing_cubes,
                    "failures": failures,
                },
                indent=2,
            ),
            encoding="utf-8",
        )

    return 2 if missing_results or missing_cubes else 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze local cost-curve and band-disagreement diagnostics."
    )
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--mids-file", type=Path, default=DEFAULT_MIDS_FILE)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--cost-dir", type=Path, default=DEFAULT_COST_DIR)
    parser.add_argument("--radius-m", type=float, default=1500.0)
    parser.add_argument("--flat-delta", type=float, default=0.05)
    parser.add_argument("--spread-threshold", type=float, default=0.2)
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    return analyze(_parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
