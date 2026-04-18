#!/usr/bin/env python3
"""Run original-vs-native 6SV2.1 parity validation and write a report."""

from __future__ import annotations

import argparse
from pathlib import Path

from siac.sixs_upstream_parity import run_upstream_parity_suite, write_upstream_parity_report


def _default_source_dir() -> Path:
    return Path("tmp/6s_upstream")


def _default_output_dir() -> Path:
    return Path("tmp/6s_upstream_parity")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=_default_source_dir())
    parser.add_argument("--upstream-build-dir", type=Path, default=_default_output_dir() / "upstream_exec")
    parser.add_argument("--native-build-dir", type=Path, default=_default_output_dir() / "native_build")
    parser.add_argument("--output-dir", type=Path, default=_default_output_dir() / "cases")
    parser.add_argument("--report-path", type=Path, default=_default_output_dir() / "report.json")
    parser.add_argument("--compiler", default="gfortran")
    args = parser.parse_args()

    report = run_upstream_parity_suite(
        source_dir=args.source_dir,
        upstream_build_dir=args.upstream_build_dir,
        native_build_dir=args.native_build_dir,
        output_dir=args.output_dir,
        compiler=args.compiler,
    )
    write_upstream_parity_report(report, args.report_path)

    case_reports = report["cases"]
    print(f"Report: {args.report_path.resolve()}")
    print(f"Original executable: {report['upstream_executable']}")
    print(f"Native module: {report['native_module']}")
    for case in case_reports:
        status = "MATCH" if case["matched"] else "MISMATCH"
        print(
            f"{status:8s} {case['name']}: "
            f"{case['matched_variable_count']}/{case['compared_variable_count']} variables matched"
        )
        if case["mismatches"]:
            first = case["mismatches"][0]
            print(
                "         "
                f"first mismatch: {first['name']} "
                f"original={first['original']} native={first['native']} tol={first['tolerance']}"
            )
    return 0 if report["all_cases_matched"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
