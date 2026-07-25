#!/usr/bin/env python3
"""Summarize acixThree-compatible S2 AOD replicate outputs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import median

import numpy as np


def _load(path: Path) -> dict[str, object] | None:
    try:
        with path.open("r", encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        return None
    except json.JSONDecodeError as exc:
        return {"matchup_id": path.parent.name, "status": "bad_json", "error": str(exc)}


def _records(root: Path, mids: list[str] | None) -> list[dict[str, object]]:
    ids = mids if mids is not None else sorted(p.name for p in root.iterdir() if p.is_dir())
    out: list[dict[str, object]] = []
    for mid in ids:
        rec = _load(root / mid / "summary.json")
        if rec is None:
            rec = {"matchup_id": mid, "status": "missing_summary"}
        out.append(rec)
    return out


def _as_float(value: object) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")


def summarize(records: list[dict[str, object]]) -> dict[str, object]:
    usable = []
    for rec in records:
        retrieved = _as_float(rec.get("acixthree_s2_aod550_window_mean"))
        truth = _as_float(rec.get("truth_aod550"))
        err = retrieved - truth
        if np.isfinite(retrieved) and np.isfinite(truth):
            usable.append((rec, err, abs(err)))

    within = [
        rec
        for rec, _, _ in usable
        if bool(rec.get("within_expected_error"))
    ]
    abs_err = [ae for _, _, ae in usable]
    bias = [err for _, err, _ in usable]
    return {
        "total_requested": len(records),
        "usable": len(usable),
        "within_expected_error": len(within),
        "within_expected_error_pct": 100.0 * len(within) / len(usable) if usable else None,
        "mae": float(np.mean(abs_err)) if abs_err else None,
        "median_abs_error": float(median(abs_err)) if abs_err else None,
        "bias": float(np.mean(bias)) if bias else None,
        "missing_or_failed": len(records) - len(usable),
    }


def write_csv(records: list[dict[str, object]], path: Path) -> None:
    fields = [
        "matchup_id",
        "site",
        "truth_aod550",
        "acixthree_s2_aod550_window_mean",
        "error",
        "expected_error",
        "within_expected_error",
        "hyper_status",
        "status",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=Path)
    p.add_argument("--mids-file", type=Path)
    p.add_argument("--csv", type=Path)
    p.add_argument("--json", type=Path)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    mids = None
    if args.mids_file is not None:
        mids = [line.strip() for line in args.mids_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    records = _records(args.root, mids)
    summary = summarize(records)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    if args.csv is not None:
        write_csv(records, args.csv)
    if args.json is not None:
        args.json.write_text(json.dumps({"summary": summary, "records": records}, indent=2, default=float), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
