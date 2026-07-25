"""Analyze campaign-250 failures for surface-driven AOD retrievals.

This is a read-only diagnostic over saved Phase-D JSON outputs. It deliberately
does not use external AOD products or write selector outputs; the goal is to
understand which failures a genuine surface-prior solve still needs to fix.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")

SURFACE_RESULT_DIRS: dict[str, str] = {
    "K2": "phaseD_results_campaign250_K2",
    "N1_f006": "phaseD_results_campaign250_N1_f006",
    "N2_loose": "phaseD_results_campaign250_N2_f006_loose",
    "N3_iter": "phaseD_results_campaign250_N3_f006_iter",
    "O3_et20": "phaseD_results_campaign250_O3_et20_iter",
    "O4_debias": "phaseD_results_campaign250_O4_debias_iter",
    "O5_shape": "phaseD_results_campaign250_O5_shape",
    "O6_b11b12": "phaseD_results_campaign250_O6_b11b12",
    "P1_tau": "phaseD_results_campaign250_P1_tau",
    "Q1_static_ext": "phaseD_results_campaign250_Q1_static_ext",
    "Q2_tau_ext": "phaseD_results_campaign250_Q2_tau_ext",
    "R1_gated": "phaseD_results_campaign250_R1_gated",
    "R2_full": "phaseD_results_campaign250_R2_full",
    "monthly_median3": "phaseD_results_campaign250_l2a_monthly_median3_scene_mean",
    "pc_mean3": "phaseD_results_campaign250_l2a_pc_production_mean3_scene_mean",
}

R2_LABEL = "R2_full"


def _safe_float(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def _within_ee(record: dict[str, Any]) -> bool:
    return bool(record.get("within_ee"))


def _load_dir(root: Path, dirname: str) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for path in sorted((root / dirname).glob("*.json")):
        with path.open("r", encoding="utf-8") as handle:
            record = json.load(handle)
        records[str(record.get("matchup_id") or path.stem)] = record
    return records


def _aod_bin(truth: float) -> str:
    if truth < 0.1:
        return "<0.1"
    if truth < 0.2:
        return "0.1-0.2"
    if truth < 0.4:
        return "0.2-0.4"
    if truth < 0.6:
        return "0.4-0.6"
    if truth < 1.0:
        return "0.6-1.0"
    if truth < 1.5:
        return "1.0-1.5"
    return ">1.5"


def _cost_bucket(cost: float) -> str:
    if not math.isfinite(cost):
        return "nan"
    if cost < 1.0:
        return "<1"
    if cost < 2.0:
        return "1-2"
    if cost < 5.0:
        return "2-5"
    if cost < 20.0:
        return "5-20"
    if cost < 100.0:
        return "20-100"
    return ">=100"


def _action_bucket(row: dict[str, Any]) -> str:
    truth = float(row["truth"])
    err = float(row["err"])
    cost = float(row["cost"])
    if cost >= 20.0 and err < 0.0:
        return "A high-cost underread: RT/smoke-cloud optical mismatch"
    if truth >= 0.6 and cost < 2.0 and err < 0.0:
        return "B low-cost high-AOD underread: wrong aerosol optics / weak spectral constraint"
    if 0.4 <= truth < 0.6 and err < 0.0:
        return "C moderate underread: surface prior/uncertainty or backstop too low"
    if err > 0.0:
        return "D overread: dark/contaminated prior or cloud/adjacency"
    if truth < 0.4 and err < 0.0:
        return "E clean/moderate underread: overbright prior or tight backstop"
    return "F other"


def _row_for_mid(
    mid: str,
    records: dict[str, dict[str, dict[str, Any]]],
) -> dict[str, Any]:
    r2 = records[R2_LABEL].get(mid, {})
    fallback = next((arm[mid] for arm in records.values() if mid in arm), {})
    source = r2 or fallback
    solver = r2.get("solver") if isinstance(r2.get("solver"), dict) else {}
    truth = _safe_float(source.get("truth"))
    retrieved = _safe_float(r2.get("retrieved"))
    winners = tuple(
        label for label, arm in records.items() if mid in arm and _within_ee(arm[mid])
    )
    return {
        "mid": mid,
        "site": source.get("site") or mid.split("__")[0],
        "truth": truth,
        "retrieved": retrieved,
        "err": retrieved - truth,
        "bin": _aod_bin(truth) if math.isfinite(truth) else "nan",
        "within_r2": _within_ee(r2),
        "surface_oracle": bool(winners),
        "winners": winners,
        "cost": _safe_float(solver.get("cost_final_per_band")),
        "std": _safe_float(solver.get("aot_std")),
        "cloud": _safe_float(r2.get("cloud_frac")),
    }


def _print_counter(title: str, counter: Counter[str], order: tuple[str, ...] = ()) -> None:
    print(f"\n{title}")
    seen = set()
    for key in order:
        print(f"  {key}: {counter.get(key, 0)}")
        seen.add(key)
    for key, value in counter.most_common():
        if key not in seen:
            print(f"  {key}: {value}")


def _write_mid_list(path: Path, mids: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{mid}\n" for mid in mids), encoding="utf-8")
    print(f"Wrote {len(mids)} mids -> {path}")


def analyze(
    root: Path,
    *,
    top: int,
    write_unsolved_mids: Path | None = None,
    write_r2_fixed_mids: Path | None = None,
    write_r2_fail_mids: Path | None = None,
    write_bad_suite_mids: Path | None = None,
) -> int:
    records = {label: _load_dir(root, dirname) for label, dirname in SURFACE_RESULT_DIRS.items()}
    mids = sorted(set().union(*(set(arm) for arm in records.values())))
    if not mids:
        raise ValueError(f"No saved campaign records found under {root}")

    print("Surface-only campaign scores")
    for label, arm in records.items():
        total = len(arm)
        within = sum(_within_ee(record) for record in arm.values())
        pct = within / total * 100.0 if total else 0.0
        print(f"  {label:16s} {within:3d}/{total:<3d} {pct:5.1f}%")

    rows = [_row_for_mid(mid, records) for mid in mids]
    oracle = sum(bool(row["surface_oracle"]) for row in rows)
    print(f"\nSurface-only oracle: {oracle}/{len(rows)} = {oracle / len(rows) * 100.0:.1f}%")

    r2_rows = [row for row in rows if math.isfinite(float(row["err"]))]
    r2_fail = [row for row in r2_rows if not row["within_r2"]]
    print(f"R2 misses: {len(r2_fail)}")
    print(
        "  under/over:",
        Counter("under" if float(row["err"]) < 0.0 else "over" for row in r2_fail),
    )
    _print_counter(
        "R2 misses by truth-AOD bin",
        Counter(str(row["bin"]) for row in r2_fail),
        ("<0.1", "0.1-0.2", "0.2-0.4", "0.4-0.6", "0.6-1.0", "1.0-1.5", ">1.5"),
    )
    _print_counter(
        "R2 misses by cost_final_per_band",
        Counter(_cost_bucket(float(row["cost"])) for row in r2_fail),
        ("<1", "1-2", "2-5", "5-20", "20-100", ">=100", "nan"),
    )

    r2_fixed = [row for row in r2_fail if row["surface_oracle"]]
    print(f"\nR2 misses fixed by another saved surface-only arm: {len(r2_fixed)}")
    fixed_by: Counter[str] = Counter()
    for row in r2_fixed:
        fixed_by.update(str(label) for label in row["winners"])
    for label, count in fixed_by.most_common():
        print(f"  {label:16s} {count}")

    unsolved = [row for row in rows if not row["surface_oracle"]]
    print(f"\nUnsolved by every saved surface-only arm: {len(unsolved)}")
    _print_counter(
        "Unsolved cases by action bucket",
        Counter(_action_bucket(row) for row in unsolved),
    )

    print(f"\nTop {top} unsolved cases by |R2 error|")
    for row in sorted(
        (row for row in unsolved if math.isfinite(float(row["err"]))),
        key=lambda item: abs(float(item["err"])),
        reverse=True,
    )[:top]:
        print(
            f"  {str(row['mid'])[:54]:54s} "
            f"bin={str(row['bin']):7s} truth={float(row['truth']):5.2f} "
            f"R2={float(row['retrieved']):5.2f} err={float(row['err']):+6.2f} "
            f"cost={float(row['cost']):8.2f} std={float(row['std']):5.2f} "
            f"cloud={float(row['cloud']):4.2f}"
        )

    if write_unsolved_mids is not None:
        _write_mid_list(write_unsolved_mids, sorted(str(row["mid"]) for row in unsolved))
    if write_r2_fixed_mids is not None:
        _write_mid_list(write_r2_fixed_mids, sorted(str(row["mid"]) for row in r2_fixed))
    if write_r2_fail_mids is not None:
        _write_mid_list(write_r2_fail_mids, sorted(str(row["mid"]) for row in r2_fail))
    if write_bad_suite_mids is not None:
        bad_mids = sorted({str(row["mid"]) for row in unsolved + r2_fixed})
        _write_mid_list(write_bad_suite_mids, bad_mids)
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze saved campaign-250 surface-driven AOD failures."
    )
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument(
        "--write-unsolved-mids",
        type=Path,
        default=None,
        help="Optional path for the 61 mids unsolved by every saved surface-only arm.",
    )
    parser.add_argument(
        "--write-r2-fixed-mids",
        type=Path,
        default=None,
        help="Optional path for R2 misses fixed by another saved surface-only arm.",
    )
    parser.add_argument(
        "--write-r2-fail-mids",
        type=Path,
        default=None,
        help="Optional path for all R2 misses.",
    )
    parser.add_argument(
        "--write-bad-suite-mids",
        type=Path,
        default=None,
        help="Optional path for unsolved mids plus R2 misses fixed by another surface-only arm.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    return analyze(
        args.root,
        top=int(args.top),
        write_unsolved_mids=args.write_unsolved_mids,
        write_r2_fixed_mids=args.write_r2_fixed_mids,
        write_r2_fail_mids=args.write_r2_fail_mids,
        write_bad_suite_mids=args.write_bad_suite_mids,
    )


if __name__ == "__main__":
    raise SystemExit(main())
