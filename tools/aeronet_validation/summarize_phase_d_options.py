"""Rank Phase-D result directories without hiding retrieval-coverage losses."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")


@dataclass(frozen=True)
class Outcome:
    matchup_id: str
    present: bool
    status: str
    valid: bool
    within_ee: bool
    aod: float
    truth: float
    error: float
    ee: float


def _finite(record: dict[str, Any], *keys: str) -> float:
    for key in keys:
        value = record.get(key)
        if value is None:
            continue
        try:
            result = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(result):
            return result
    return math.nan


def load_outcome(
    directory: Path,
    matchup_id: str,
    *,
    retrieved_key: str | None = None,
) -> Outcome:
    path = directory / f"{matchup_id}.json"
    if not path.exists():
        return Outcome(matchup_id, False, "MISSING", False, False, *(math.nan,) * 4)
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return Outcome(matchup_id, True, "MALFORMED", False, False, *(math.nan,) * 4)

    status = str(record.get("status", "OK")).upper()
    aod = (
        _finite(record, retrieved_key)
        if retrieved_key is not None
        else _finite(record, "retrieved", "aod", "retrieved_aod")
    )
    truth = _finite(record, "truth", "aeronet_aod", "reference_aod")
    error = _finite(record, "err", "error")
    if not math.isfinite(error) and math.isfinite(aod) and math.isfinite(truth):
        error = aod - truth
    ee = _finite(record, "ee", "ee_threshold", "expected_error")
    if not math.isfinite(ee) and math.isfinite(truth):
        ee = 0.05 + 0.15 * truth
    valid = status == "OK" and math.isfinite(aod) and math.isfinite(truth)
    explicit_hit = record.get("within_ee") if retrieved_key is None else None
    if isinstance(explicit_hit, bool):
        within_ee = valid and explicit_hit
    else:
        within_ee = valid and math.isfinite(error) and math.isfinite(ee) and abs(error) <= ee
    return Outcome(matchup_id, True, status, valid, within_ee, aod, truth, error, ee)


def summarize(label: str, outcomes: list[Outcome]) -> dict[str, Any]:
    expected = len(outcomes)
    present = sum(row.present for row in outcomes)
    valid = sum(row.valid for row in outcomes)
    hits = sum(row.within_ee for row in outcomes)
    errors = [row.error for row in outcomes if row.valid and math.isfinite(row.error)]
    return {
        "label": label,
        "hits": hits,
        "expected": expected,
        "strict_within_ee_pct": 100.0 * hits / expected if expected else None,
        "present": present,
        "valid": valid,
        "coverage_pct": 100.0 * valid / expected if expected else None,
        "valid_within_ee_pct": 100.0 * hits / valid if valid else None,
        "no_valid_observation": sum(row.status == "NO_VALID_OBSERVATION" for row in outcomes),
        "failed": sum(row.status == "FAILED" for row in outcomes),
        "malformed": sum(row.status == "MALFORMED" for row in outcomes),
        "missing": expected - present,
        "rmse": math.sqrt(sum(value * value for value in errors) / len(errors))
        if errors
        else None,
        "bias": sum(errors) / len(errors) if errors else None,
    }


def paired_summary(
    baseline_label: str,
    baseline: list[Outcome],
    candidate_label: str,
    candidate: list[Outcome],
) -> dict[str, Any]:
    paired_present = [
        (base, cand)
        for base, cand in zip(baseline, candidate, strict=True)
        if base.present and cand.present
    ]
    gains = [
        base.matchup_id
        for base, cand in paired_present
        if not base.within_ee and cand.within_ee
    ]
    losses = [
        base.matchup_id
        for base, cand in paired_present
        if base.within_ee and not cand.within_ee
    ]
    shared_valid = [
        (base, cand)
        for base, cand in paired_present
        if base.valid and cand.valid
    ]
    return {
        "baseline": baseline_label,
        "candidate": candidate_label,
        "gains": len(gains),
        "losses": len(losses),
        "net": len(gains) - len(losses),
        "paired_present": len(paired_present),
        "shared_valid": len(shared_valid),
        "gain_ids": gains,
        "loss_ids": losses,
    }


def evaluate(
    *,
    mids: list[str],
    baseline_label: str,
    baseline_dir: Path,
    candidates: list[tuple[str, Path]],
    retrieved_key: str | None = None,
) -> dict[str, Any]:
    baseline = [
        load_outcome(baseline_dir, matchup_id, retrieved_key=retrieved_key)
        for matchup_id in mids
    ]
    summaries = [summarize(baseline_label, baseline)]
    paired = []
    for label, directory in candidates:
        outcomes = [
            load_outcome(directory, matchup_id, retrieved_key=retrieved_key)
            for matchup_id in mids
        ]
        summaries.append(summarize(label, outcomes))
        paired.append(paired_summary(baseline_label, baseline, label, outcomes))
    return {"retrieved_key": retrieved_key or "reported", "summaries": summaries, "paired": paired}


def _resolve(root: Path, value: Path) -> Path:
    return value if value.is_absolute() else root / value


def _labelled_path(value: str) -> tuple[str, Path]:
    label, separator, path = value.partition("=")
    if not separator or not label or not path:
        raise argparse.ArgumentTypeError("expected LABEL=RESULT_DIRECTORY")
    return label, Path(path)


def _pct(value: float | None) -> str:
    return "-" if value is None else f"{value:.1f}%"


def _print_table(result: dict[str, Any]) -> None:
    paired = {row["candidate"]: row for row in result["paired"]}
    print("| option | strict EE | valid coverage | EE of valid | no valid | failed | missing | net |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    for row in result["summaries"]:
        comparison = paired.get(row["label"])
        net = "-" if comparison is None else f'{comparison["net"]:+d}'
        print(
            f'| {row["label"]} | {row["hits"]}/{row["expected"]} '
            f'({_pct(row["strict_within_ee_pct"])}) | {row["valid"]}/{row["expected"]} '
            f'({_pct(row["coverage_pct"])}) | {_pct(row["valid_within_ee_pct"])} | '
            f'{row["no_valid_observation"]} | {row["failed"]} | {row["missing"]} | {net} |'
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--mids-file", type=Path, required=True)
    parser.add_argument("--baseline", type=_labelled_path, required=True)
    parser.add_argument("--candidate", type=_labelled_path, action="append", default=[])
    parser.add_argument(
        "--retrieved-key",
        help="Score an alternative saved AOD field, e.g. retrieved_winmed or retrieved_station.",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    mids_path = _resolve(args.root, args.mids_file)
    mids = [line.strip() for line in mids_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    baseline_label, baseline_value = args.baseline
    candidates = [(label, _resolve(args.root, path)) for label, path in args.candidate]
    result = evaluate(
        mids=mids,
        baseline_label=baseline_label,
        baseline_dir=_resolve(args.root, baseline_value),
        candidates=candidates,
        retrieved_key=args.retrieved_key,
    )
    _print_table(result)
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
