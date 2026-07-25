"""Check diagnostic reruns for campaign-250 surface-driven bad cases."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.aeronet_validation.analyze_surface_driven_failures import (
    DEFAULT_ROOT,
    R2_LABEL,
    SURFACE_RESULT_DIRS,
    _action_bucket,
    _load_dir,
    _row_for_mid,
    _safe_float,
    _within_ee,
)


DEFAULT_MIDS_FILE = DEFAULT_ROOT / "campaign250_surface_bad89_mids.txt"
DEFAULT_RESULTS_DIR = DEFAULT_ROOT / "phaseD_results_campaign250_surface_bad89_R2_diag_20260705"

REQUIRED_DIAGNOSTIC_KEYS = (
    "surface_cost_curve_min_aot",
    "surface_cost_curve_min_cost",
    "surface_cost_curve_second_aot",
    "surface_cost_curve_second_delta",
    "surface_cost_curve_relative_second_delta",
    "surface_cost_curve_curvature",
    "surface_final_aot_median",
    "surface_final_aot_node",
    "surface_band_B02_argmin_aot",
    "surface_band_B02_argmin_cost",
    "surface_band_B02_cost_final_node",
    "surface_band_B02_residual_final_node",
    "surface_band_B04_argmin_aot",
    "surface_band_B04_argmin_cost",
    "surface_band_B04_cost_final_node",
    "surface_band_B04_residual_final_node",
    "surface_band_argmin_spread",
    "surface_tau_gate_configured",
    "surface_tau_available",
    "surface_tau_gate_fired",
)


def _read_mid_list(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _load_records(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    records: dict[str, dict[str, Any]] = {}
    for json_path in sorted(path.glob("*.json")):
        with json_path.open("r", encoding="utf-8") as handle:
            record = json.load(handle)
        records[str(record.get("matchup_id") or json_path.stem)] = record
    return records


def _result_dir(root: Path, value: Path) -> Path:
    if value.is_absolute():
        return value
    return root / value


def _median(values: list[float]) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    mid = len(finite) // 2
    if len(finite) % 2:
        return finite[mid]
    return 0.5 * (finite[mid - 1] + finite[mid])


def _quantile(values: list[float], q: float) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    if len(finite) == 1:
        return finite[0]
    pos = (len(finite) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return finite[lo]
    return finite[lo] * (hi - pos) + finite[hi] * (pos - lo)


def _fmt(value: float, digits: int = 3) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.{digits}f}"


def _ee(record: dict[str, Any], truth: float) -> float:
    value = _safe_float(record.get("ee_threshold"))
    if math.isfinite(value):
        return value
    return 0.05 + 0.15 * truth


def _solver(record: dict[str, Any]) -> dict[str, Any]:
    solver = record.get("solver")
    return solver if isinstance(solver, dict) else {}


def _curve_min_is_edge(solver: dict[str, Any]) -> bool:
    if solver.get("surface_cost_curve_min_at_edge") is True:
        return True
    min_aot = _safe_float(solver.get("surface_cost_curve_min_aot"))
    axis_min = _safe_float(solver.get("surface_aot_axis_min"))
    axis_max = _safe_float(solver.get("surface_aot_axis_max"))
    return math.isfinite(min_aot) and (
        math.isclose(min_aot, axis_min, rel_tol=0.0, abs_tol=1e-9)
        or math.isclose(min_aot, axis_max, rel_tol=0.0, abs_tol=1e-9)
    )


def _missing_required_keys(solver: dict[str, Any]) -> list[str]:
    missing_keys = [key for key in REQUIRED_DIAGNOSTIC_KEYS if key not in solver]
    if "surface_cost_curve_curvature" in missing_keys and _curve_min_is_edge(solver):
        missing_keys.remove("surface_cost_curve_curvature")
    return missing_keys


def _case_group(row: dict[str, Any]) -> str:
    if not row["surface_oracle"]:
        return _action_bucket(row)
    return "G R2 miss recoverable by another saved surface-only arm"


def _print_counter(title: str, counts: Counter[str]) -> None:
    print(f"\n{title}")
    for key, count in counts.most_common():
        print(f"  {key}: {count}")


def _print_group_table(rows: list[dict[str, Any]]) -> None:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["group"])].append(row)

    print("\nBy bad-case group")
    for group, items in sorted(grouped.items()):
        total = len(items)
        within = sum(bool(item["within"]) for item in items)
        miss = total - within
        curve_min_ok = sum(bool(item["curve_min_within"]) for item in items)
        final_node_ok = sum(bool(item["final_node_within"]) for item in items)
        print(
            f"  {group}: {within}/{total} within EE; "
            f"misses={miss}; curve-min-within={curve_min_ok}; "
            f"final-node-within={final_node_ok}; "
            f"median_err={_fmt(_median([float(item['err']) for item in items]))}"
        )


def check(args: argparse.Namespace) -> int:
    root = Path(args.root)
    mids_file = Path(args.mids_file)
    result_dir = _result_dir(root, Path(args.result_dir))
    baseline_dir = _result_dir(root, Path(args.baseline_dir))

    expected_mids = _read_mid_list(mids_file)
    expected = set(expected_mids)
    records = _load_records(result_dir)
    baseline = _load_records(baseline_dir)
    surface_records = {
        label: _load_dir(root, dirname) for label, dirname in SURFACE_RESULT_DIRS.items()
    }

    missing = sorted(expected - set(records))
    extra = sorted(set(records) - expected)
    present = [records[mid] for mid in expected_mids if mid in records]
    ok_records = [record for record in present if str(record.get("status")).upper() == "OK"]
    failed_records = [record for record in present if str(record.get("status")).upper() != "OK"]
    known_baseline_failed = [
        record
        for record in failed_records
        if str(
            baseline.get(str(record.get("matchup_id") or ""), {}).get("status")
        ).upper()
        != "OK"
    ]
    new_failed_records = [
        record
        for record in failed_records
        if str(
            baseline.get(str(record.get("matchup_id") or ""), {}).get("status")
        ).upper()
        == "OK"
    ]
    within = sum(_within_ee(record) for record in ok_records)

    print("Diagnostic result coverage")
    print(f"  results_dir: {result_dir}")
    print(f"  expected: {len(expected_mids)}")
    print(f"  present: {len(present)}")
    print(f"  missing: {len(missing)}")
    print(f"  extra: {len(extra)}")
    print(f"  OK: {len(ok_records)}")
    print(f"  FAILED/non-OK: {len(failed_records)}")
    print(f"  known baseline non-OK: {len(known_baseline_failed)}")
    print(f"  new FAILED/non-OK: {len(new_failed_records)}")
    if ok_records:
        print(f"  within EE: {within}/{len(ok_records)} = {within / len(ok_records) * 100.0:.1f}%")

    missing_key_counts: Counter[str] = Counter()
    complete_diagnostics = 0
    for record in ok_records:
        solver = _solver(record)
        missing_keys = _missing_required_keys(solver)
        if missing_keys:
            missing_key_counts.update(missing_keys)
        else:
            complete_diagnostics += 1
    print(
        f"  complete required diagnostics: "
        f"{complete_diagnostics}/{len(ok_records)}"
    )
    if missing_key_counts:
        _print_counter("Missing diagnostic keys", missing_key_counts)

    rows: list[dict[str, Any]] = []
    baseline_deltas: list[float] = []
    baseline_status_changes: Counter[str] = Counter()
    for mid in expected_mids:
        record = records.get(mid)
        if record is None or str(record.get("status")).upper() != "OK":
            continue
        old_row = _row_for_mid(mid, surface_records)
        solver = _solver(record)
        truth = _safe_float(record.get("truth"))
        retrieved = _safe_float(record.get("retrieved"))
        err = retrieved - truth
        ee = _ee(record, truth)
        curve_min = _safe_float(solver.get("surface_cost_curve_min_aot"))
        final_node = _safe_float(solver.get("surface_final_aot_median"))
        old = baseline.get(mid, {})
        old_retrieved = _safe_float(old.get("retrieved"))
        delta = retrieved - old_retrieved
        if math.isfinite(delta):
            baseline_deltas.append(abs(delta))
        old_within = _within_ee(old) if old else False
        new_within = _within_ee(record)
        baseline_status_changes[f"{old_within}->{new_within}"] += 1
        rows.append(
            {
                "mid": mid,
                "group": _case_group(old_row),
                "truth": truth,
                "retrieved": retrieved,
                "err": err,
                "within": new_within,
                "curve_min": curve_min,
                "curve_min_err": curve_min - truth,
                "curve_min_within": math.isfinite(curve_min)
                and abs(curve_min - truth) <= ee,
                "final_node": final_node,
                "final_node_err": final_node - truth,
                "final_node_within": math.isfinite(final_node)
                and abs(final_node - truth) <= ee,
                "cost": _safe_float(solver.get("cost_final_per_band")),
                "curve_delta": _safe_float(
                    solver.get("surface_cost_curve_relative_second_delta")
                ),
                "curvature": _safe_float(solver.get("surface_cost_curve_curvature")),
                "band_spread": _safe_float(solver.get("surface_band_argmin_spread")),
                "tau_gate": bool(solver.get("surface_tau_gate_fired")),
                "tau_available": bool(solver.get("surface_tau_available")),
            }
        )

    if rows:
        print("\nBaseline R2 comparison")
        print(f"  status transitions: {dict(baseline_status_changes)}")
        print(f"  median |retrieved_new - retrieved_R2|: {_fmt(_median(baseline_deltas), 6)}")
        print(f"  max |retrieved_new - retrieved_R2|: {_fmt(max(baseline_deltas), 6)}")

        _print_group_table(rows)

        misses = [row for row in rows if not row["within"]]
        curve_min_ok_misses = sum(bool(row["curve_min_within"]) for row in misses)
        final_node_ok_misses = sum(bool(row["final_node_within"]) for row in misses)
        tau_available = sum(bool(row["tau_available"]) for row in rows)
        tau_fired = sum(bool(row["tau_gate"]) for row in rows)
        curve_edge = sum(_curve_min_is_edge(_solver(records[str(row["mid"])])) for row in rows)
        high_spread = sum(
            math.isfinite(float(row["band_spread"])) and float(row["band_spread"]) >= 0.2
            for row in misses
        )
        flat_curve = sum(
            math.isfinite(float(row["curve_delta"])) and float(row["curve_delta"]) <= 0.05
            for row in misses
        )
        print("\nDiagnostic clues")
        print(f"  tau prior available: {tau_available}/{len(rows)}")
        print(f"  tau gate fired: {tau_fired}/{len(rows)}")
        print(f"  cost-curve minimum on AOD-axis edge: {curve_edge}/{len(rows)}")
        print(f"  misses where cost-curve min would be within EE: {curve_min_ok_misses}/{len(misses)}")
        print(f"  misses where final solver node would be within EE: {final_node_ok_misses}/{len(misses)}")
        print(f"  misses with B02/B04 argmin spread >= 0.2 AOD: {high_spread}/{len(misses)}")
        print(f"  misses with flat/ambiguous cost curve delta <= 5%: {flat_curve}/{len(misses)}")
        print(
            "  median miss diagnostics: "
            f"err={_fmt(_median([float(row['err']) for row in misses]))}, "
            f"curve_min_err={_fmt(_median([float(row['curve_min_err']) for row in misses]))}, "
            f"band_spread={_fmt(_median([float(row['band_spread']) for row in misses]))}, "
            f"curve_delta={_fmt(_median([float(row['curve_delta']) for row in misses]))}, "
            f"cost={_fmt(_median([float(row['cost']) for row in misses]))}"
        )

        print(f"\nTop {int(args.top)} current misses by |error|")
        for row in sorted(misses, key=lambda item: abs(float(item["err"])), reverse=True)[
            : int(args.top)
        ]:
            print(
                f"  {str(row['mid'])[:54]:54s} "
                f"err={float(row['err']):+6.3f} truth={float(row['truth']):5.3f} "
                f"ret={float(row['retrieved']):5.3f} curve_min={_fmt(float(row['curve_min']))} "
                f"final_node={_fmt(float(row['final_node']))} cost={_fmt(float(row['cost']), 2)} "
                f"spread={_fmt(float(row['band_spread']))} "
                f"delta={_fmt(float(row['curve_delta']))} tau_gate={row['tau_gate']} "
                f"group={row['group']}"
            )

    if missing or new_failed_records or missing_key_counts:
        return 2
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check a surface-driven bad-case diagnostic rerun."
    )
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--mids-file", type=Path, default=DEFAULT_MIDS_FILE)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument(
        "--baseline-dir",
        type=Path,
        default=Path("phaseD_results_campaign250_R2_full"),
    )
    parser.add_argument("--top", type=int, default=30)
    return parser.parse_args()


def main() -> int:
    return check(_parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
