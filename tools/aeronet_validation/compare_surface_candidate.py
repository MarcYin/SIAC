"""Compare a surface-driven AOD candidate against a baseline result set.

The campaign notes use this for quick, auditable checks of a new surface-prior
candidate against the current R2 baseline: total score, fixes, breaks, closer
retrievals, and simple truth-AOD/error-regime splits.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any


DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")


@dataclass(frozen=True)
class ScoredRecord:
    matchup_id: str
    status_ok: bool
    within_ee: bool
    aod: float
    truth: float
    err: float
    ee: float
    raw: dict[str, Any]


def _finite_float(record: dict[str, Any], *keys: str) -> float:
    for key in keys:
        value = record.get(key)
        if value is None:
            continue
        try:
            out = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(out):
            return out
    return math.nan


def _within_ee(record: dict[str, Any], err: float, ee: float) -> bool:
    if "within_ee" in record:
        return bool(record["within_ee"])
    if "flag" in record:
        return str(record["flag"]).upper() == "OK"
    return math.isfinite(err) and math.isfinite(ee) and abs(err) <= ee


def _load_record(directory: Path, matchup_id: str) -> ScoredRecord | None:
    path = directory / f"{matchup_id}.json"
    if not path.exists():
        return None
    record = json.loads(path.read_text())
    status_ok = str(record.get("status", "OK")).upper() == "OK"
    aod = _finite_float(record, "aod", "retrieved_aod", "retrieved")
    truth = _finite_float(record, "truth", "aeronet_aod", "reference_aod")
    err = _finite_float(record, "err", "error")
    if not math.isfinite(err) and math.isfinite(aod) and math.isfinite(truth):
        err = aod - truth
    ee = _finite_float(record, "ee", "expected_error")
    within = status_ok and _within_ee(record, err, ee)
    return ScoredRecord(
        matchup_id=matchup_id,
        status_ok=status_ok,
        within_ee=within,
        aod=aod,
        truth=truth,
        err=err,
        ee=ee,
        raw=record,
    )


def _truth_bin(truth: float) -> str:
    if not math.isfinite(truth):
        return "unknown"
    bins = (
        (0.1, "<0.1"),
        (0.2, "0.1-0.2"),
        (0.4, "0.2-0.4"),
        (0.6, "0.4-0.6"),
        (1.0, "0.6-1.0"),
    )
    for limit, label in bins:
        if truth < limit:
            return label
    return ">=1.0"


def _direction(err: float) -> str:
    if not math.isfinite(err):
        return "unknown"
    return "under" if err < 0 else "over"


def _fmt(value: float, digits: int = 3) -> str:
    return "nan" if not math.isfinite(value) else f"{value:.{digits}f}"


def _print_score(label: str, records: list[ScoredRecord | None]) -> None:
    present = [record for record in records if record is not None]
    ok = [record for record in present if record.status_ok]
    hits = sum(record.within_ee for record in present)
    ok_hits = sum(record.within_ee for record in ok)
    pct = hits / len(records) * 100.0 if records else 0.0
    ok_pct = ok_hits / len(ok) * 100.0 if ok else 0.0
    missing = len(records) - len(present)
    print(
        f"{label}: {hits}/{len(records)} raw = {pct:.1f}% within EE; "
        f"{ok_hits}/{len(ok)} OK = {ok_pct:.1f}%"
        + (f"; missing {missing}" if missing else "")
    )


def compare(
    *,
    mids: list[str],
    baseline_dir: Path,
    candidate_dir: Path,
    label: str,
    top_n: int,
) -> int:
    baseline = [_load_record(baseline_dir, mid) for mid in mids]
    candidate = [_load_record(candidate_dir, mid) for mid in mids]

    _print_score("baseline", baseline)
    _print_score(label, candidate)

    paired = [
        (base, cand)
        for base, cand in zip(baseline, candidate, strict=True)
        if base is not None
        and cand is not None
        and base.status_ok
        and cand.status_ok
        and math.isfinite(base.err)
        and math.isfinite(cand.err)
    ]
    fixes = [(base, cand) for base, cand in paired if not base.within_ee and cand.within_ee]
    breaks = [(base, cand) for base, cand in paired if base.within_ee and not cand.within_ee]
    both_miss = [(base, cand) for base, cand in paired if not base.within_ee and not cand.within_ee]
    closer = [(base, cand) for base, cand in paired if abs(cand.err) < abs(base.err)]
    worse = [(base, cand) for base, cand in paired if abs(cand.err) > abs(base.err)]

    print("\nCandidate versus baseline")
    print(f"  paired OK records: {len(paired)}")
    print(f"  fixes: {len(fixes)}")
    print(f"  breaks: {len(breaks)}")
    print(f"  both miss: {len(both_miss)}")
    print(f"  closer abs error: {len(closer)}")
    print(f"  worse abs error: {len(worse)}")

    by_truth: Counter[str] = Counter()
    by_direction: Counter[str] = Counter()
    for base, cand in fixes:
        by_truth[_truth_bin(base.truth)] += 1
        by_direction[_direction(base.err)] += 1
    if fixes:
        print("\nFixes by baseline truth-AOD bin")
        for key in ("<0.1", "0.1-0.2", "0.2-0.4", "0.4-0.6", "0.6-1.0", ">=1.0", "unknown"):
            if by_truth[key]:
                print(f"  {key}: {by_truth[key]}")
        print("Fixes by baseline error direction")
        for key in ("under", "over", "unknown"):
            if by_direction[key]:
                print(f"  {key}: {by_direction[key]}")

    print(f"\nTop {top_n} fixes by candidate error")
    for base, cand in sorted(fixes, key=lambda pair: abs(pair[1].err))[:top_n]:
        print(
            f"  {base.matchup_id:58s} truth={_fmt(base.truth)} "
            f"base={_fmt(base.aod)} err={_fmt(base.err)} "
            f"cand={_fmt(cand.aod)} err={_fmt(cand.err)}"
        )

    print(f"\nTop {top_n} breaks by baseline margin")
    for base, cand in sorted(breaks, key=lambda pair: abs(pair[0].err - pair[0].ee))[:top_n]:
        print(
            f"  {base.matchup_id:58s} truth={_fmt(base.truth)} "
            f"base={_fmt(base.aod)} err={_fmt(base.err)} ee={_fmt(base.ee)} "
            f"cand={_fmt(cand.aod)} err={_fmt(cand.err)}"
        )

    print(f"\nTop {top_n} remaining candidate misses by |error|")
    for base, cand in sorted(both_miss, key=lambda pair: abs(pair[1].err), reverse=True)[:top_n]:
        print(
            f"  {base.matchup_id:58s} truth={_fmt(base.truth)} "
            f"base={_fmt(base.aod)} err={_fmt(base.err)} "
            f"cand={_fmt(cand.aod)} err={_fmt(cand.err)}"
        )
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--mids-file", type=Path, required=True)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_ROOT / "phaseD_results_campaign250_R2_full")
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--label", default="candidate")
    parser.add_argument("--top-n", type=int, default=12)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    mids_path = args.mids_file
    if not mids_path.is_absolute():
        mids_path = args.root / mids_path
    candidate_dir = args.candidate_dir
    if not candidate_dir.is_absolute():
        candidate_dir = args.root / candidate_dir
    baseline_dir = args.baseline_dir
    if not baseline_dir.is_absolute():
        baseline_dir = args.root / baseline_dir
    mids = [line.strip() for line in mids_path.read_text().splitlines() if line.strip()]
    return compare(
        mids=mids,
        baseline_dir=baseline_dir,
        candidate_dir=candidate_dir,
        label=args.label,
        top_n=max(0, args.top_n),
    )


if __name__ == "__main__":
    raise SystemExit(main())
