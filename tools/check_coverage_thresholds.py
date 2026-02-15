#!/usr/bin/env python3
"""Validate coverage thresholds from coverage.json.

This gate enforces:
- overall coverage >= min_total
- per-file coverage > min_file
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check total and per-file coverage thresholds from coverage.json."
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=Path("coverage.json"),
        help="Coverage JSON path (default: coverage.json).",
    )
    parser.add_argument(
        "--min-total",
        type=float,
        default=95.0,
        help="Minimum overall coverage percentage (default: 95.0).",
    )
    parser.add_argument(
        "--min-file",
        type=float,
        default=90.0,
        help="Per-file coverage must be strictly greater than this value (default: 90.0).",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="python/siac/",
        help="Only validate files under this path prefix (default: python/siac/).",
    )
    return parser.parse_args()


def _read_file_rates(json_path: Path, prefix: str) -> tuple[float, dict[str, float]]:
    with json_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    total = float(data.get("totals", {}).get("percent_covered", 0.0))

    files: dict[str, dict[str, object]] = data.get("files", {})
    per_file: dict[str, float] = {}
    for filename, payload in files.items():
        if not filename.endswith(".py"):
            continue
        if prefix and not filename.startswith(prefix):
            continue
        summary = payload.get("summary", {})
        per_file[filename] = float(summary.get("percent_covered", 0.0))

    return total, per_file


def main() -> int:
    args = _parse_args()
    if not args.json.exists():
        print(
            f"[coverage-gate] Missing coverage report: {args.json}. "
            "Run pytest with JSON report first.",
            file=sys.stderr,
        )
        return 2

    total, per_file = _read_file_rates(args.json, args.prefix)
    failures: list[str] = []

    if total < args.min_total:
        failures.append(
            f"overall coverage {total:.2f}% is below required {args.min_total:.2f}%"
        )

    low_files = sorted(
        [(name, rate) for name, rate in per_file.items() if rate <= args.min_file],
        key=lambda item: item[1],
    )
    if low_files:
        failures.append(
            f"{len(low_files)} file(s) are not above {args.min_file:.2f}% per-file coverage"
        )

    if failures:
        print("[coverage-gate] FAILED")
        for msg in failures:
            print(f"- {msg}")
        if low_files:
            print("- files below threshold:")
            for name, rate in low_files:
                print(f"  - {name}: {rate:.2f}%")
        return 1

    print("[coverage-gate] PASSED")
    print(f"- overall coverage: {total:.2f}%")
    print(f"- per-file coverage: all {len(per_file)} file(s) are above {args.min_file:.2f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
