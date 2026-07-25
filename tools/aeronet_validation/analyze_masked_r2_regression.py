"""Decompose the apparent score regression from historical to masked R2."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")


def _finite(value: object) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan


def _solver(record: dict[str, Any]) -> dict[str, Any]:
    value = record.get("solver")
    return value if isinstance(value, dict) else {}


def _prior(record: dict[str, Any], band: str) -> float:
    value = record.get("prior_boa")
    return _finite(value.get(band)) if isinstance(value, dict) else math.nan


def _valid(record: dict[str, Any] | None) -> bool:
    return bool(
        record is not None
        and str(record.get("status", "OK")).upper() == "OK"
        and math.isfinite(_finite(record.get("truth")))
        and math.isfinite(_finite(record.get("retrieved")))
    )


def _hit(record: dict[str, Any] | None) -> bool:
    if not _valid(record):
        return False
    explicit = record.get("within_ee")
    if isinstance(explicit, bool):
        return explicit
    truth = _finite(record.get("truth"))
    error = _finite(record.get("retrieved")) - truth
    return abs(error) <= 0.05 + 0.15 * truth


def _load(directory: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for path in directory.glob("*.json"):
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        records[str(record.get("matchup_id") or path.stem)] = record
    return records


def _corr(rows: list[dict[str, Any]], x_key: str, y_key: str) -> float | None:
    x = np.asarray([_finite(row.get(x_key)) for row in rows], dtype=np.float64)
    y = np.asarray([_finite(row.get(y_key)) for row in rows], dtype=np.float64)
    keep = np.isfinite(x) & np.isfinite(y)
    if int(keep.sum()) < 3 or np.std(x[keep]) == 0.0 or np.std(y[keep]) == 0.0:
        return None
    return float(np.corrcoef(x[keep], y[keep])[0, 1])


def _median(rows: list[dict[str, Any]], key: str) -> float | None:
    values = np.asarray([_finite(row.get(key)) for row in rows], dtype=np.float64)
    values = values[np.isfinite(values)]
    return float(np.median(values)) if values.size else None


def _score(records: dict[str, dict[str, Any]], mids: list[str]) -> dict[str, Any]:
    valid = sum(_valid(records.get(mid)) for mid in mids)
    hits = sum(_hit(records.get(mid)) for mid in mids)
    statuses: dict[str, int] = {}
    for mid in mids:
        record = records.get(mid)
        status = "MISSING" if record is None else str(record.get("status", "OK")).upper()
        statuses[status] = statuses.get(status, 0) + 1
    return {
        "hits": hits,
        "expected": len(mids),
        "strict_pct": 100.0 * hits / len(mids) if mids else None,
        "valid": valid,
        "coverage_pct": 100.0 * valid / len(mids) if mids else None,
        "valid_pct": 100.0 * hits / valid if valid else None,
        "statuses": statuses,
    }


def _abstention_class(record: dict[str, Any] | None) -> str:
    if record is None or str(record.get("status", "")).upper() != "NO_VALID_OBSERVATION":
        return ""
    solver = _solver(record)
    support_count = _finite(solver.get("surface_valid_observation_count"))
    solved_count = _finite(solver.get("surface_solved_pixel_count"))
    cloud_fraction = _finite(record.get("cloud_frac"))
    if math.isfinite(support_count) and support_count <= 0.0:
        if math.isfinite(cloud_fraction) and cloud_fraction >= 0.999:
            return "fully_cloud_masked"
        return "zero_admissible_support_other"
    if (
        math.isfinite(support_count)
        and support_count > 0.0
        and math.isfinite(solved_count)
        and solved_count <= 0.0
    ):
        return "valid_support_but_pool_unsolved"
    return "nonfinite_retrieval_other"


def _comparison_rows(
    old: dict[str, dict[str, Any]],
    new: dict[str, dict[str, Any]],
    mids: list[str],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for mid in mids:
        old_record = old.get(mid)
        new_record = new.get(mid)
        old_valid = _valid(old_record)
        new_valid = _valid(new_record)
        old_hit = _hit(old_record)
        new_hit = _hit(new_record)
        if old_valid and new_valid:
            if not old_hit and new_hit:
                outcome = "fix"
            elif old_hit and not new_hit:
                outcome = "break"
            elif old_hit:
                outcome = "retained_hit"
            else:
                outcome = "retained_miss"
        elif old_valid:
            outcome = "new_abstention"
        elif new_valid:
            outcome = "new_only"
        else:
            outcome = "neither_valid"

        old_aod = _finite(old_record.get("retrieved")) if old_record else math.nan
        new_aod = _finite(new_record.get("retrieved")) if new_record else math.nan
        old_b02 = _prior(old_record, "B02") if old_record else math.nan
        old_b04 = _prior(old_record, "B04") if old_record else math.nan
        new_b02 = _prior(new_record, "B02") if new_record else math.nan
        new_b04 = _prior(new_record, "B04") if new_record else math.nan
        solver = _solver(new_record) if new_record else {}
        rows.append(
            {
                "matchup_id": mid,
                "site": str((new_record or old_record or {}).get("site", "")),
                "truth": _finite((new_record or old_record or {}).get("truth")),
                "outcome": outcome,
                "old_status": str((old_record or {}).get("status", "MISSING")).upper(),
                "new_status": str((new_record or {}).get("status", "MISSING")).upper(),
                "old_hit": old_hit,
                "new_hit": new_hit,
                "old_aod": old_aod,
                "new_aod": new_aod,
                "aod_delta": new_aod - old_aod,
                "old_prior_b02": old_b02,
                "new_prior_b02": new_b02,
                "prior_b02_delta": new_b02 - old_b02,
                "old_prior_b04": old_b04,
                "new_prior_b04": new_b04,
                "prior_b04_delta": new_b04 - old_b04,
                "cloud_frac": _finite((new_record or {}).get("cloud_frac")),
                "invalid_frac": _finite((new_record or {}).get("invalid_frac")),
                "valid_support_frac": _finite(solver.get("surface_valid_observation_fraction")),
                "valid_support_count": _finite(solver.get("surface_valid_observation_count")),
                "solved_pixel_count": _finite(solver.get("surface_solved_pixel_count")),
                "cost_curve_valid_nodes": _finite(solver.get("surface_cost_curve_valid_nodes")),
                "abstention_class": _abstention_class(new_record),
                "target_tcwv_cm": _finite((new_record or {}).get("target_tcwv_cm")),
            }
        )
    return rows


def _mask_control(
    masked: dict[str, dict[str, Any]],
    bypass: dict[str, dict[str, Any]],
    mids: list[str],
) -> dict[str, Any]:
    common = [mid for mid in mids if _valid(masked.get(mid)) and _valid(bypass.get(mid))]
    fixes = [mid for mid in common if not _hit(masked[mid]) and _hit(bypass[mid])]
    breaks = [mid for mid in common if _hit(masked[mid]) and not _hit(bypass[mid])]
    bypass_only_hits = [
        mid
        for mid in mids
        if not _valid(masked.get(mid)) and _valid(bypass.get(mid)) and _hit(bypass[mid])
    ]
    return {
        "masked": _score(masked, mids),
        "cloud_bypass": _score(bypass, mids),
        "common_valid": len(common),
        "bypass_fixes_on_common": len(fixes),
        "bypass_breaks_on_common": len(breaks),
        "bypass_common_net": len(fixes) - len(breaks),
        "bypass_only_hits": len(bypass_only_hits),
        "bypass_fix_ids": fixes,
        "bypass_break_ids": breaks,
        "bypass_only_hit_ids": bypass_only_hits,
    }


def analyze(
    *,
    mids: list[str],
    old: dict[str, dict[str, Any]],
    new: dict[str, dict[str, Any]],
    control_mids: list[str],
    cloud_bypass: dict[str, dict[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = _comparison_rows(old, new, mids)
    common = [row for row in rows if row["old_status"] == "OK" and row["new_status"] == "OK"]
    fixes = [row for row in common if row["outcome"] == "fix"]
    breaks = [row for row in common if row["outcome"] == "break"]
    old_only_hits = [
        row for row in rows if row["old_hit"] and not _valid(new.get(row["matchup_id"]))
    ]
    new_only_hits = [
        row for row in rows if row["new_hit"] and not _valid(old.get(row["matchup_id"]))
    ]
    old_score = _score(old, mids)
    new_score = _score(new, mids)
    common_net = len(fixes) - len(breaks)
    coverage_net = len(new_only_hits) - len(old_only_hits)
    abstentions = [row for row in rows if row["new_status"] == "NO_VALID_OBSERVATION"]
    abstention_classes = {
        name: sum(row["abstention_class"] == name for row in abstentions)
        for name in (
            "fully_cloud_masked",
            "zero_admissible_support_other",
            "valid_support_but_pool_unsolved",
            "nonfinite_retrieval_other",
        )
    }
    common_old_hits = sum(bool(row["old_hit"]) for row in common)
    common_new_hits = sum(bool(row["new_hit"]) for row in common)
    result = {
        "old": old_score,
        "new": new_score,
        "strict_hit_delta": new_score["hits"] - old_score["hits"],
        "decomposition": {
            "common_valid": len(common),
            "common_old_hits": common_old_hits,
            "common_new_hits": common_new_hits,
            "common_old_pct": 100.0 * common_old_hits / len(common) if common else None,
            "common_new_pct": 100.0 * common_new_hits / len(common) if common else None,
            "fixes": len(fixes),
            "breaks": len(breaks),
            "common_net": common_net,
            "old_hits_lost_to_abstention": len(old_only_hits),
            "new_hits_from_old_invalid": len(new_only_hits),
            "coverage_net": coverage_net,
            "reconciled_delta": common_net + coverage_net,
        },
        "common_diagnostics": {
            "median_aod_delta": _median(common, "aod_delta"),
            "median_prior_b02_delta": _median(common, "prior_b02_delta"),
            "median_prior_b04_delta": _median(common, "prior_b04_delta"),
            "corr_aod_delta_prior_b02_delta": _corr(common, "aod_delta", "prior_b02_delta"),
            "corr_aod_delta_prior_b04_delta": _corr(common, "aod_delta", "prior_b04_delta"),
            "corr_aod_delta_cloud_frac": _corr(common, "aod_delta", "cloud_frac"),
            "corr_aod_delta_valid_support_frac": _corr(common, "aod_delta", "valid_support_frac"),
        },
        "break_diagnostics": {
            "median_aod_delta": _median(breaks, "aod_delta"),
            "median_prior_b02_delta": _median(breaks, "prior_b02_delta"),
            "median_prior_b04_delta": _median(breaks, "prior_b04_delta"),
            "median_cloud_frac": _median(breaks, "cloud_frac"),
            "median_valid_support_frac": _median(breaks, "valid_support_frac"),
        },
        "mask_only_control": _mask_control(new, cloud_bypass, control_mids),
        "new_abstention_diagnostics": {
            "total": len(abstentions),
            "classes": abstention_classes,
            "rows": abstentions,
            "pool_window_cells": 20,
            "pool_min_count": 80,
            "pool_basis": (
                "600 m configured half-width at 60 m resolution gives a 20-cell window; "
                "the solver requires max(20, window_area / 5) = 80 finite samples"
            ),
        },
        "breaks": breaks,
        "old_hits_lost_to_abstention": old_only_hits,
        "configuration_changes": [
            "historical: libRadtran LUT backend; new: native 6S backend",
            "historical: atmospheric-provider TCWV; new: fixed same-scene S2 L2A WVP",
            "historical: cloud/water masking bypassed and invalid output filled; new: cloud/water/full-support masking with NaN abstention",
            "historical and new runs use different solver-grid extents in some scenes after intervening grid/QA changes",
        ],
    }
    return result, rows


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _truth_bin(value: float) -> str:
    if value < 0.2:
        return "clean"
    if value < 0.6:
        return "moderate"
    if value < 1.0:
        return "high"
    return "very_high"


def _diagnostic_subset(rows: list[dict[str, Any]]) -> list[str]:
    selected = [
        row["matchup_id"] for row in rows if row["outcome"] in {"fix", "break", "new_abstention"}
    ]
    selected_set = set(selected)
    for outcome in ("retained_hit", "retained_miss"):
        for truth_bin in ("clean", "moderate", "high", "very_high"):
            controls = [
                row["matchup_id"]
                for row in rows
                if row["outcome"] == outcome
                and _truth_bin(_finite(row["truth"])) == truth_bin
                and row["matchup_id"] not in selected_set
            ][:3]
            selected.extend(controls)
            selected_set.update(controls)
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--mids-file", default="campaign250_mids.txt")
    parser.add_argument("--control-mids-file", default="surface_e1_matched_mids.txt")
    parser.add_argument(
        "--old-dir", default="phaseD_results_campaign250_R2_full_localdiag_20260705"
    )
    parser.add_argument(
        "--new-dir", default="phaseD_results_campaign250_masked_r2_l2awvp_6s_20260710"
    )
    parser.add_argument(
        "--cloud-bypass-dir", default="phaseD_results_surface_e1_masked_cloud_bypass_20260710"
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--rows-output", type=Path)
    parser.add_argument("--subset-output", type=Path)
    args = parser.parse_args()

    mids = [
        line.strip()
        for line in (args.root / args.mids_file).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    control_mids = [
        line.strip()
        for line in (args.root / args.control_mids_file).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    result, rows = analyze(
        mids=mids,
        old=_load(args.root / args.old_dir),
        new=_load(args.root / args.new_dir),
        control_mids=control_mids,
        cloud_bypass=_load(args.root / args.cloud_bypass_dir),
    )
    text = json.dumps(result, indent=2)
    print(text)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    if args.rows_output is not None:
        _write_rows(args.rows_output, rows)
    if args.subset_output is not None:
        args.subset_output.parent.mkdir(parents=True, exist_ok=True)
        subset = _diagnostic_subset(rows)
        args.subset_output.write_text("\n".join(subset) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
