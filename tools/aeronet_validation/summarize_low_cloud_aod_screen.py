"""Summarize low-cloud AOD methods against a fixed campaign cohort."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_COHORT = "campaign250_lowcloud20_mids.txt"
DEFAULT_OUTPUT = "reports/aod-low-cloud-20260711/screen-analysis.json"
DEFAULT_BENCHMARKS = (
    ("historical_r2", "phaseD_results_campaign250_R2_full_localdiag_20260705"),
    ("masked_r2", "phaseD_results_campaign250_masked_r2_l2awvp_6s_20260710"),
    ("multisource_tree_v1", "phaseD_results_campaign250_multisource_tree_v1"),
    ("multisource_tree_v2", "phaseD_results_campaign250_multisource_tree_v2"),
)
DEFAULT_CANDIDATES = (
    (
        "modis_single_day_followup",
        "phaseD_results_masked_modis_prior_modis_single_day_20260711",
    ),
    (
        "modis_timeseries_smooth_followup",
        "phaseD_results_masked_modis_prior_modis_timeseries_smooth_20260711",
    ),
    (
        "s2_swir_nir_anchor_followup",
        "phaseD_results_masked_modis_prior_swir_nir_anchored_20260711",
    ),
    (
        "modis_monthly_anchor_calibrated",
        "phaseD_results_masked_modis_prior_modis_monthly_swir_nir_anchored_lowcloud20_20260711",
    ),
    (
        "modis_monthly_anchor_loose_backstop",
        "phaseD_results_masked_modis_prior_modis_monthly_swir_nir_anchored_"
        "loose_backstop_lowcloud20_20260711",
    ),
    (
        "modis_monthly_anchor_spread_only",
        "phaseD_results_masked_modis_prior_modis_monthly_swir_nir_anchored_"
        "spread_only_lowcloud20_20260711",
    ),
    (
        "modis_monthly_multigrid_current",
        "experiments/lowcloud20_modis_monthly_database_20260711/runs/monthly_database",
    ),
    (
        "modis_monthly_multigrid_june",
        "runs/monthly_database",
    ),
    (
        "modis_monthly_multigrid_spread_only",
        "experiments/lowcloud20_modis_monthly_database_20260711/runs/monthly_database_spread_only",
    ),
    (
        "modis_monthly_multigrid_spread_only_aot25",
        "experiments/lowcloud20_modis_monthly_database_20260711/"
        "runs/monthly_database_spread_only_aot25",
    ),
    (
        "modis_monthly_multigrid_legacy_resample_control",
        "experiments/lowcloud20_modis_monthly_database_20260711/"
        "runs/monthly_database_spread_only_aot25_legacy_resample",
    ),
)
EXCLUDED_SURFACE_DIRECTORY_TOKENS = (
    "highaod",
    "multisource",
    "surface_bad",
)


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _valid(record: dict[str, Any] | None) -> bool:
    return bool(
        record
        and str(record.get("status", "")).upper() == "OK"
        and _finite(record.get("truth")) is not None
        and _finite(record.get("retrieved")) is not None
    )


def _hit(record: dict[str, Any] | None) -> bool:
    if not _valid(record):
        return False
    assert record is not None
    truth = float(record["truth"])
    retrieved = float(record["retrieved"])
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth + 1.0e-12


def _normalize_record(
    record: dict[str, Any],
    *,
    truth_by_id: dict[str, float] | None = None,
) -> dict[str, Any]:
    normalized = dict(record)
    matchup_id = str(normalized.get("matchup_id", ""))
    if _finite(normalized.get("retrieved")) is None:
        normalized["retrieved"] = _finite(normalized.get("aot_window_mean"))
    if _finite(normalized.get("truth")) is None and truth_by_id is not None:
        normalized["truth"] = truth_by_id.get(matchup_id)
    if str(normalized.get("status", "")).lower() == "ok":
        normalized["status"] = "OK"
    return normalized


def _load_records(
    directory: Path,
    *,
    truth_by_id: dict[str, float] | None = None,
) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    paths = sorted(directory.glob("*.json"))
    if not paths:
        paths = sorted(directory.glob("**/result.json"))
    for path in paths:
        record = _normalize_record(
            json.loads(path.read_text(encoding="utf-8")),
            truth_by_id=truth_by_id,
        )
        matchup_id = str(record.get("matchup_id") or path.stem)
        records[matchup_id] = record
    return records


def _summary(
    records: dict[str, dict[str, Any]],
    cohort: list[str],
) -> dict[str, Any]:
    present = [matchup_id for matchup_id in cohort if matchup_id in records]
    valid = [matchup_id for matchup_id in cohort if _valid(records.get(matchup_id))]
    hits = [matchup_id for matchup_id in cohort if _hit(records.get(matchup_id))]
    errors = np.asarray(
        [records[mid]["retrieved"] - records[mid]["truth"] for mid in valid],
        dtype=np.float64,
    )
    statuses: dict[str, int] = {}
    for matchup_id in present:
        status = str(records[matchup_id].get("status", "MISSING_STATUS")).upper()
        statuses[status] = statuses.get(status, 0) + 1
    return {
        "expected": len(cohort),
        "present": len(present),
        "valid": len(valid),
        "hits": len(hits),
        "strict_rate": len(hits) / len(cohort) if cohort else None,
        "valid_hit_rate": len(hits) / len(valid) if valid else None,
        "rmse": float(np.sqrt(np.mean(np.square(errors)))) if errors.size else None,
        "bias": float(np.mean(errors)) if errors.size else None,
        "statuses": statuses,
        "hit_matchup_ids": hits,
    }


def _surface_campaign_directories(root: Path) -> list[Path]:
    directories = []
    for directory in sorted(root.glob("phaseD_results_campaign250*")):
        if not directory.is_dir():
            continue
        if any(token in directory.name for token in EXCLUDED_SURFACE_DIRECTORY_TOKENS):
            continue
        if sum(1 for _ in directory.glob("*.json")) < 244:
            continue
        directories.append(directory)
    return directories


def _parse_candidate(value: str) -> tuple[str, str]:
    try:
        name, directory = value.split("=", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("candidate must be NAME=DIRECTORY") from exc
    if not name.strip() or not directory.strip():
        raise argparse.ArgumentTypeError("candidate must be NAME=DIRECTORY")
    return name.strip(), directory.strip()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--cohort", default=DEFAULT_COHORT)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--expect-count", type=int, default=152)
    parser.add_argument("--target-rate", type=float, default=0.87)
    parser.add_argument(
        "--candidate",
        action="append",
        type=_parse_candidate,
        dest="candidates",
        help="Additional or replacement candidate as NAME=DIRECTORY.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    root = args.root
    cohort = [
        line.strip()
        for line in (root / args.cohort).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(cohort) != args.expect_count or len(cohort) != len(set(cohort)):
        raise SystemExit(f"expected {args.expect_count} unique cohort IDs, found {len(cohort)}")

    surface_directories = _surface_campaign_directories(root)
    surface_records = {
        directory.name: _load_records(directory) for directory in surface_directories
    }
    surface_method_summaries = {
        name: _summary(records, cohort) for name, records in surface_records.items()
    }
    truth_by_id = {
        matchup_id: truth
        for records in surface_records.values()
        for matchup_id, record in records.items()
        if (truth := _finite(record.get("truth"))) is not None
    }
    surface_hits = {
        matchup_id
        for matchup_id in cohort
        if any(_hit(records.get(matchup_id)) for records in surface_records.values())
    }
    screen_ids = [matchup_id for matchup_id in cohort if matchup_id not in surface_hits]

    candidate_specs = args.candidates or list(DEFAULT_CANDIDATES)
    candidates: dict[str, Any] = {}
    candidate_hit_union: set[str] = set()
    for name, directory_name in candidate_specs:
        records = _load_records(root / directory_name, truth_by_id=truth_by_id)
        summary = _summary(records, cohort)
        screen_summary = _summary(records, screen_ids)
        unique_hits = [matchup_id for matchup_id in screen_ids if _hit(records.get(matchup_id))]
        candidate_hit_union.update(summary["hit_matchup_ids"])
        candidates[name] = {
            "directory": str(root / directory_name),
            "cohort": summary,
            "screen": screen_summary,
            "unique_vs_existing_surface_oracle": len(unique_hits),
            "unique_hit_matchup_ids": unique_hits,
        }

    required_hits = math.ceil(float(args.target_rate) * len(cohort) - 1.0e-12)
    expanded_hits = surface_hits | candidate_hit_union
    benchmarks = {
        name: {
            "directory": str(root / directory_name),
            "cohort": _summary(
                _load_records(root / directory_name, truth_by_id=truth_by_id),
                cohort,
            ),
        }
        for name, directory_name in DEFAULT_BENCHMARKS
    }
    result = {
        "cohort_file": str(root / args.cohort),
        "cohort_count": len(cohort),
        "target_rate": float(args.target_rate),
        "required_hits": required_hits,
        "surface_campaign_directory_count": len(surface_directories),
        "surface_campaign_directories": [path.name for path in surface_directories],
        "surface_method_summaries": surface_method_summaries,
        "existing_surface_oracle_hits": len(surface_hits),
        "existing_surface_oracle_rate": len(surface_hits) / len(cohort),
        "screen_count": len(screen_ids),
        "screen_matchup_ids": screen_ids,
        "benchmarks": benchmarks,
        "candidates": candidates,
        "expanded_surface_oracle_hits": len(expanded_hits),
        "expanded_surface_oracle_rate": len(expanded_hits) / len(cohort),
        "expanded_oracle_meets_target": len(expanded_hits) >= required_hits,
    }
    output_path = root / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        f"existing surface oracle: {len(surface_hits)}/{len(cohort)}; "
        f"expanded: {len(expanded_hits)}/{len(cohort)}; target: {required_hits}/{len(cohort)}"
    )
    for name, candidate in candidates.items():
        summary = candidate["cohort"]
        print(
            f"{name}: present={summary['present']} valid={summary['valid']} "
            f"hits={summary['hits']} unique={candidate['unique_vs_existing_surface_oracle']}"
        )
    for name, benchmark in benchmarks.items():
        summary = benchmark["cohort"]
        print(
            f"{name}: present={summary['present']} valid={summary['valid']} "
            f"hits={summary['hits']} strict={100.0 * summary['strict_rate']:.1f}%"
        )
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
