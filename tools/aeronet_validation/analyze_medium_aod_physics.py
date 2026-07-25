"""Score physical AOD experiments on a fixed low-cloud AERONET cohort."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterable

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_BASELINE = DEFAULT_ROOT / "phaseD_results_lowcloud20_native_maiac_adaptive_b03_chi2_20260713"
DEFAULT_MIDS = DEFAULT_ROOT / "campaign250_lowcloud20_mids.txt"
DEFAULT_MATCHUPS = (
    DEFAULT_ROOT
    / "reports/aod-production-reproduction-spec-20260713/downloads/target-matchups.csv"
)
DEFAULT_EXCLUDED = (
    "Dakar_Belair__T28PBB_20240630T113321",
    "IER_Cinzana__T29PRQ_20240302T103839",
    "HESS__T33KXQ_20240803T083559",
)
DEFAULT_SPLIT_SEED = 20260713
DEFAULT_HOLDOUT_FOLDS = (2, 3)


@dataclass(frozen=True)
class Metrics:
    n: int
    within_ee: int
    within_ee_rate: float | None
    bias: float | None
    mae: float | None
    rmse: float | None
    median_error: float | None
    under_ee: int
    over_ee: int

    def as_dict(self) -> dict[str, Any]:
        return self.__dict__.copy()


def _finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def ee_threshold(truth: float) -> float:
    return 0.05 + 0.15 * truth


def within_ee(retrieved: float, truth: float) -> bool:
    return abs(retrieved - truth) <= ee_threshold(truth)


def _percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("percentile requires at least one value")
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def compute_metrics(records: Iterable[dict[str, Any]]) -> Metrics:
    pairs: list[tuple[float, float]] = []
    for record in records:
        truth = _finite(record.get("truth"))
        retrieved = _finite(record.get("retrieved"))
        if truth is not None and retrieved is not None:
            pairs.append((retrieved, truth))
    if not pairs:
        return Metrics(0, 0, None, None, None, None, None, 0, 0)

    errors = [retrieved - truth for retrieved, truth in pairs]
    hits = [within_ee(retrieved, truth) for retrieved, truth in pairs]
    under = sum(error < -ee_threshold(truth) for error, (_, truth) in zip(errors, pairs))
    over = sum(error > ee_threshold(truth) for error, (_, truth) in zip(errors, pairs))
    return Metrics(
        n=len(pairs),
        within_ee=sum(hits),
        within_ee_rate=sum(hits) / len(pairs),
        bias=sum(errors) / len(errors),
        mae=sum(abs(error) for error in errors) / len(errors),
        rmse=math.sqrt(sum(error * error for error in errors) / len(errors)),
        median_error=_percentile(errors, 0.5),
        under_ee=under,
        over_ee=over,
    )


def _load_results(directory: Path, mids: list[str]) -> tuple[dict[str, dict[str, Any]], dict[str, int]]:
    records: dict[str, dict[str, Any]] = {}
    statuses: dict[str, int] = {}
    for matchup_id in mids:
        path = directory / f"{matchup_id}.json"
        if not path.exists():
            statuses["MISSING"] = statuses.get("MISSING", 0) + 1
            continue
        try:
            record = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            statuses["MALFORMED"] = statuses.get("MALFORMED", 0) + 1
            continue
        status = str(record.get("status", "UNKNOWN")).upper()
        statuses[status] = statuses.get(status, 0) + 1
        if status == "OK" and _finite(record.get("retrieved")) is not None:
            records[matchup_id] = record
    return records, statuses


def _load_matchups(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            matchup_id = row.get("matchup_id", "")
            if matchup_id and matchup_id != "matchup_id":
                rows[matchup_id] = row
    return rows


def _subset_ids(truth: dict[str, float], excluded: set[str]) -> dict[str, list[str]]:
    return {
        "low": [mid for mid, value in truth.items() if value < 0.25],
        "medium": [mid for mid, value in truth.items() if 0.25 <= value <= 0.85],
        "medium_025_040": [mid for mid, value in truth.items() if 0.25 <= value < 0.40],
        "medium_040_055": [mid for mid, value in truth.items() if 0.40 <= value < 0.55],
        "medium_055_070": [mid for mid, value in truth.items() if 0.55 <= value < 0.70],
        "medium_070_085": [mid for mid, value in truth.items() if 0.70 <= value <= 0.85],
        "high": [mid for mid, value in truth.items() if value > 0.85],
        "non_extreme": [mid for mid in truth if mid not in excluded],
        "all": list(truth),
    }


def _site_group_folds(
    baseline: dict[str, dict[str, Any]],
    mids: list[str],
    *,
    seed: int,
) -> dict[str, int]:
    import numpy as np
    from sklearn.model_selection import StratifiedGroupKFold

    labels: list[str] = []
    groups: list[str] = []
    for matchup_id in mids:
        record = baseline[matchup_id]
        truth = float(record["truth"])
        regime = "low" if truth < 0.25 else "medium" if truth <= 0.85 else "high"
        outcome = "hit" if within_ee(float(record["retrieved"]), truth) else "miss"
        labels.append(regime if regime == "high" else f"{regime}_{outcome}")
        groups.append(str(record["site"]))

    assignments: dict[str, int] = {}
    splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
    placeholder = np.zeros(len(mids), dtype=np.uint8)
    for fold, (_, test_indices) in enumerate(splitter.split(placeholder, labels, groups)):
        for index in test_indices:
            assignments[mids[int(index)]] = fold
    if len(assignments) != len(mids):
        raise RuntimeError("site-group split did not assign every matchup")
    return assignments


def _paired_transitions(
    baseline: dict[str, dict[str, Any]],
    candidate: dict[str, dict[str, Any]],
    ids: Iterable[str],
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    gains = losses = stable_hits = stable_misses = 0
    for matchup_id in ids:
        if matchup_id not in baseline or matchup_id not in candidate:
            continue
        base = baseline[matchup_id]
        cand = candidate[matchup_id]
        truth = float(base["truth"])
        base_value = float(base["retrieved"])
        cand_value = float(cand["retrieved"])
        base_hit = within_ee(base_value, truth)
        cand_hit = within_ee(cand_value, truth)
        if not base_hit and cand_hit:
            transition = "gain"
            gains += 1
        elif base_hit and not cand_hit:
            transition = "loss"
            losses += 1
        elif base_hit:
            transition = "stable_hit"
            stable_hits += 1
        else:
            transition = "stable_miss"
            stable_misses += 1
        rows.append(
            {
                "matchup_id": matchup_id,
                "truth": truth,
                "baseline": base_value,
                "candidate": cand_value,
                "baseline_error": base_value - truth,
                "candidate_error": cand_value - truth,
                "delta_retrieved": cand_value - base_value,
                "transition": transition,
            }
        )
    return {
        "n_paired": len(rows),
        "gains": gains,
        "losses": losses,
        "net_gains": gains - losses,
        "stable_hits": stable_hits,
        "stable_misses": stable_misses,
        "rows": rows,
    }


def _case_row(
    matchup_id: str,
    baseline: dict[str, Any],
    candidate: dict[str, Any] | None,
    matchup: dict[str, str] | None,
    variant: str,
    excluded: set[str],
) -> dict[str, Any]:
    truth = float(baseline["truth"])
    base_value = float(baseline["retrieved"])
    candidate_value = _finite((candidate or {}).get("retrieved"))
    solver = (candidate or baseline).get("solver") or {}
    atmo = (candidate or baseline).get("atmo_prior") or {}
    source = matchup or {}
    return {
        "variant": variant,
        "matchup_id": matchup_id,
        "site": baseline.get("site"),
        "truth": truth,
        "baseline": base_value,
        "candidate": candidate_value,
        "baseline_error": base_value - truth,
        "candidate_error": None if candidate_value is None else candidate_value - truth,
        "baseline_within_ee": within_ee(base_value, truth),
        "candidate_within_ee": (
            None if candidate_value is None else within_ee(candidate_value, truth)
        ),
        "excluded_extreme": matchup_id in excluded,
        "elevation_m": _finite(source.get("elevation_m")),
        "angstrom": _finite(source.get("aeronet_angstrom_mean")),
        "scene_cloud_cover": _finite(source.get("scene_cloud_cover")),
        "atmo_aot": _finite(atmo.get("aot_median")),
        "atmo_tcwv_cm": _finite(atmo.get("tcwv_median")),
        "surface_curve_min_aot": _finite(solver.get("surface_cost_curve_min_aot")),
        "surface_band_argmin_spread": _finite(solver.get("surface_band_argmin_spread")),
        "surface_cost_per_band": _finite(solver.get("surface_static_cost_per_band")),
        "surface_b02_argmin_aot": _finite(solver.get("surface_band_B02_argmin_aot")),
        "surface_b03_argmin_aot": _finite(solver.get("surface_band_B03_argmin_aot")),
        "surface_b04_argmin_aot": _finite(solver.get("surface_band_B04_argmin_aot")),
    }


def _variant_arg(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("variant must be NAME=RESULT_DIRECTORY")
    name, path = value.split("=", 1)
    if not name.strip() or not path.strip():
        raise argparse.ArgumentTypeError("variant must be NAME=RESULT_DIRECTORY")
    return name.strip(), Path(path).expanduser()


def _fmt(value: float | None, *, percent: bool = False) -> str:
    if value is None:
        return "-"
    return f"{100.0 * value:.2f}%" if percent else f"{value:+.4f}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--mids", type=Path, default=DEFAULT_MIDS)
    parser.add_argument("--matchups", type=Path, default=DEFAULT_MATCHUPS)
    parser.add_argument("--variant", action="append", type=_variant_arg, default=[])
    parser.add_argument("--exclude", action="append", default=list(DEFAULT_EXCLUDED))
    parser.add_argument("--split-seed", type=int, default=DEFAULT_SPLIT_SEED)
    parser.add_argument(
        "--holdout-fold",
        action="append",
        type=int,
        dest="holdout_folds",
        help="Held-out site-group fold (repeatable); defaults to 0 and 1.",
    )
    parser.add_argument(
        "--unlock-holdout",
        action="store_true",
        help="Score held-out cases after the physical recipe has been frozen.",
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    mids = [line.strip() for line in args.mids.read_text().splitlines() if line.strip()]
    baseline, baseline_statuses = _load_results(args.baseline, mids)
    if len(baseline) != len(mids):
        raise SystemExit(
            f"baseline must have {len(mids)} OK finite records; got {len(baseline)} "
            f"({baseline_statuses})"
        )
    truth = {mid: float(baseline[mid]["truth"]) for mid in mids}
    excluded = set(args.exclude)
    folds = _site_group_folds(baseline, mids, seed=args.split_seed)
    holdout_folds = set(
        DEFAULT_HOLDOUT_FOLDS if args.holdout_folds is None else args.holdout_folds
    )
    invalid_folds = holdout_folds - set(range(5))
    if invalid_folds:
        raise SystemExit(f"holdout folds must be in 0..4; got {sorted(invalid_folds)}")
    development_ids = [mid for mid in mids if folds[mid] not in holdout_folds]
    holdout_ids = [mid for mid in mids if folds[mid] in holdout_folds]
    scoring_ids = mids if args.unlock_holdout else development_ids
    scoring_truth = {mid: truth[mid] for mid in scoring_ids}
    subsets = _subset_ids(scoring_truth, excluded)
    if args.unlock_holdout:
        development_subsets = _subset_ids({mid: truth[mid] for mid in development_ids}, excluded)
        holdout_subsets = _subset_ids({mid: truth[mid] for mid in holdout_ids}, excluded)
        subsets.update({f"development_{name}": ids for name, ids in development_subsets.items()})
        subsets.update({f"holdout_{name}": ids for name, ids in holdout_subsets.items()})
    matchups = _load_matchups(args.matchups)

    all_results: dict[str, dict[str, dict[str, Any]]] = {"baseline": baseline}
    status_counts: dict[str, dict[str, int]] = {"baseline": baseline_statuses}
    for name, directory in args.variant:
        records, statuses = _load_results(directory, mids)
        all_results[name] = records
        status_counts[name] = statuses

    summary: dict[str, Any] = {
        "cohort": {
            "mids": str(args.mids),
            "n": len(mids),
            "ee": "abs(retrieved-truth) <= 0.05 + 0.15*truth",
            "medium_range": [0.25, 0.85],
            "excluded_extreme_ids": sorted(excluded),
            "score_scope": "all_unlocked" if args.unlock_holdout else "development_locked",
            "split_seed": args.split_seed,
            "holdout_folds": sorted(holdout_folds),
            "development_n": len(development_ids),
            "holdout_n": len(holdout_ids),
            "fold_by_matchup_id": folds,
        },
        "variants": {},
    }
    case_rows: list[dict[str, Any]] = []
    for name, records in all_results.items():
        subset_metrics = {
            subset_name: compute_metrics(records[mid] for mid in ids if mid in records).as_dict()
            for subset_name, ids in subsets.items()
        }
        transitions = {
            subset_name: _paired_transitions(baseline, records, ids)
            for subset_name, ids in subsets.items()
            if name != "baseline"
        }
        summary["variants"][name] = {
            "status_counts": status_counts[name],
            "metrics": subset_metrics,
            "transitions": transitions,
        }
        for matchup_id in scoring_ids:
            case_rows.append(
                _case_row(
                    matchup_id,
                    baseline[matchup_id],
                    records.get(matchup_id),
                    matchups.get(matchup_id),
                    name,
                    excluded,
                )
            )

    scope = "all" if args.unlock_holdout else "development"
    print(f"score scope: {scope}; holdout folds={sorted(holdout_folds)}")
    print("variant                 OK/152  medium within    medium bias   low within      non-extreme")
    for name, payload in summary["variants"].items():
        metrics = payload["metrics"]
        medium = metrics["medium"]
        low = metrics["low"]
        non_extreme = metrics["non_extreme"]
        ok = payload["status_counts"].get("OK", 0)
        print(
            f"{name:23s} {ok:3d}/152  "
            f"{medium['within_ee']:2d}/{medium['n']:2d} "
            f"{_fmt(medium['within_ee_rate'], percent=True):>8s}  "
            f"{_fmt(medium['bias']):>11s}  "
            f"{low['within_ee']:2d}/{low['n']:2d}  "
            f"{non_extreme['within_ee']:3d}/{non_extreme['n']:3d}"
        )

    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        if case_rows:
            with (args.output_dir / "cases.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(case_rows[0]))
                writer.writeheader()
                writer.writerows(case_rows)


if __name__ == "__main__":
    main()
