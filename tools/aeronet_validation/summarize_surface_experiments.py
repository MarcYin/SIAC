"""Summarize and pair tagged surface-retrieval experiments."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
RESULTS = ROOT / "seasonal_val"


def _load_tag(tag: str) -> dict[str, dict[str, Any]]:
    rows = {}
    for path in RESULTS.glob(f"*__seas_maiac_{tag}.json"):
        try:
            row = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        rows[str(row["matchup_id"])] = row
    return rows


def _hit(row: dict[str, Any]) -> bool:
    if row.get("within_ee") is not None:
        return bool(row["within_ee"])
    if row.get("flag") is not None:
        return str(row["flag"]).upper() == "OK"
    return abs(float(row["err"])) <= float(row["ee"])


def _aod(row: dict[str, Any]) -> float:
    return float(row.get("aod", row.get("retrieved")))


def _truth(row: dict[str, Any]) -> float:
    return float(row["truth"])


def _pearson(xs: list[float], ys: list[float]) -> float | None:
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if int(valid.sum()) < 3 or np.std(x[valid]) == 0.0 or np.std(y[valid]) == 0.0:
        return None
    return float(np.corrcoef(x[valid], y[valid])[0, 1])


def _summary(
    tag: str,
    records: dict[str, dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    expected = list(metadata)
    rows = [records[matchup_id] for matchup_id in expected if matchup_id in records]
    hits = sum(_hit(row) for row in rows)
    errors = np.asarray([_aod(row) - _truth(row) for row in rows], dtype=float)
    by_role = {}
    for role in sorted({str(row["role"]) for row in metadata.values()}):
        selected = [
            records[matchup_id]
            for matchup_id, meta in metadata.items()
            if meta["role"] == role and matchup_id in records
        ]
        by_role[role] = {
            "hits": sum(_hit(row) for row in selected),
            "available": len(selected),
            "expected": sum(meta["role"] == role for meta in metadata.values()),
        }
    by_bin = {}
    for label in sorted({str(row["aod_bin"]) for row in metadata.values()}):
        selected = [
            records[matchup_id]
            for matchup_id, meta in metadata.items()
            if meta["aod_bin"] == label and matchup_id in records
        ]
        by_bin[label] = {
            "hits": sum(_hit(row) for row in selected),
            "available": len(selected),
            "expected": sum(meta["aod_bin"] == label for meta in metadata.values()),
        }
    return {
        "tag": tag,
        "hits": hits,
        "available": len(rows),
        "expected": len(expected),
        "expected_pct": 100.0 * hits / len(expected),
        "available_pct": 100.0 * hits / len(rows) if rows else None,
        "rmse": float(np.sqrt(np.mean(np.square(errors)))) if rows else None,
        "bias": float(np.mean(errors)) if rows else None,
        "by_role": by_role,
        "by_aod_bin": by_bin,
    }


def _paired(
    baseline_tag: str,
    baseline: dict[str, dict[str, Any]],
    tag: str,
    records: dict[str, dict[str, Any]],
    metadata: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    common = [
        matchup_id for matchup_id in metadata if matchup_id in baseline and matchup_id in records
    ]
    gains = [
        matchup_id
        for matchup_id in common
        if not _hit(baseline[matchup_id]) and _hit(records[matchup_id])
    ]
    losses = [
        matchup_id
        for matchup_id in common
        if _hit(baseline[matchup_id]) and not _hit(records[matchup_id])
    ]
    deltas = [_aod(records[mid]) - _aod(baseline[mid]) for mid in common]

    def _tcwv(row: dict[str, Any]) -> float:
        value = row.get("target_tcwv_cm")
        try:
            return float(value) if value is not None else math.nan
        except (TypeError, ValueError):
            return math.nan

    tcwv = [_tcwv(records[mid]) - _tcwv(baseline[mid]) for mid in common]
    return {
        "baseline": baseline_tag,
        "candidate": tag,
        "common": len(common),
        "gains": len(gains),
        "losses": len(losses),
        "net": len(gains) - len(losses),
        "median_aod_delta": float(np.median(deltas)) if deltas else None,
        "tcwv_delta_vs_aod_delta_pearson": _pearson(tcwv, deltas),
        "gain_ids": gains,
        "loss_ids": losses,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("tags", nargs="+")
    parser.add_argument("--manifest", type=Path, default=ROOT / "surface_e1_matched_manifest.json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    metadata = {row["matchup_id"]: row for row in manifest["rows"]}
    records = {tag: _load_tag(tag) for tag in args.tags}
    summaries = [_summary(tag, records[tag], metadata) for tag in args.tags]
    baseline_tag = args.tags[0]
    paired = [
        _paired(baseline_tag, records[baseline_tag], tag, records[tag], metadata)
        for tag in args.tags[1:]
    ]
    result = {
        "manifest": str(args.manifest),
        "baseline": baseline_tag,
        "summaries": summaries,
        "paired": paired,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
