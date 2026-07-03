"""Summarize saved AERONET validation JSON records.

This is intentionally lightweight: it scores an existing result set without
rerunning retrievals. It supports both the Phase-D JSON shape
(``within_ee``/``status``) and the older ``seasonal_val`` shape
(``flag == "OK"`` for within-EE).
"""

from __future__ import annotations

import argparse
import glob
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class SavedResultSummary:
    total: int
    within_ee: int
    pct: float
    out_sites: tuple[str, ...]


def _record_id(record: dict[str, Any], path: Path) -> str:
    return str(record.get("matchup_id") or record.get("site") or path.stem)


def _is_within_ee(record: dict[str, Any]) -> bool:
    if "within_ee" in record:
        return bool(record["within_ee"])
    if "flag" in record:
        return str(record["flag"]).upper() == "OK"
    raise ValueError("Result JSON must contain either 'within_ee' or 'flag'.")


def expand_result_patterns(patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        expanded = sorted(Path(match) for match in glob.glob(pattern))  # noqa: PTH207
        if expanded:
            paths.extend(expanded)
            continue
        path = Path(pattern)
        if path.exists():
            paths.append(path)
    return paths


def summarize_result_files(paths: list[Path]) -> SavedResultSummary:
    if not paths:
        raise ValueError("No result JSON files matched.")

    within = 0
    out_sites: list[str] = []
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            record = json.load(handle)
        if _is_within_ee(record):
            within += 1
        else:
            out_sites.append(_record_id(record, path))

    total = len(paths)
    return SavedResultSummary(
        total=total,
        within_ee=within,
        pct=within / total * 100.0,
        out_sites=tuple(out_sites),
    )


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize saved SIAC/AERONET validation result JSON files."
    )
    parser.add_argument("patterns", nargs="+", help="JSON file path or glob pattern.")
    parser.add_argument(
        "--target-pct",
        type=float,
        default=85.0,
        help="Minimum within-EE percentage required for exit 0.",
    )
    parser.add_argument(
        "--expect-count",
        type=int,
        default=None,
        help="Optional exact number of JSON files expected.",
    )
    parser.add_argument("--label", default="saved-results", help="Label shown in output.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    paths = expand_result_patterns(list(args.patterns))
    summary = summarize_result_files(paths)

    print(
        f"{args.label}: {summary.within_ee}/{summary.total} = {summary.pct:.1f}% "
        f"within EE"
    )
    if summary.out_sites:
        print("OUT:", ", ".join(summary.out_sites))

    if args.expect_count is not None and summary.total != args.expect_count:
        print(f"FAILED: expected {args.expect_count} records, found {summary.total}")
        return 2
    if summary.pct < float(args.target_pct):
        print(f"FAILED: {summary.pct:.1f}% is below target {float(args.target_pct):.1f}%")
        return 1
    print(f"PASSED: {summary.pct:.1f}% >= target {float(args.target_pct):.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
