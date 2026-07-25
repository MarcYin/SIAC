"""Analyze controlled low-cloud aerosol-species retrieval experiments."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
COHORT = ROOT / "campaign250_lowcloud20_mids.txt"
DEFAULT_OUTPUT = ROOT / "reports/aod-aerosol-species-20260712"
METHODS = {
    "lut_continental": ROOT / "phaseD_results_lowcloud20_singleprior_b03_chi2_20260711",
    "sixs_continental": (
        ROOT / "phaseD_results_lowcloud20_aerosol_species_sixs_continental_full_t4_20260712"
    ),
    "sixs_cci3": ROOT / "phaseD_results_lowcloud20_aerosol_species_sixs_cci3_full_t4_20260712",
    "sixs_cci_exact_smoke": (
        ROOT / "phaseD_results_lowcloud20_aerosol_species_sixs_cci_exact_exact_smoke_t4_20260712"
    ),
    "libradtran_continental_smoke": (
        ROOT / "phaseD_results_lowcloud20_aerosol_species_libradtran_continental_20260712"
    ),
}
CCI_COSTS = ROOT / "phaseD_cost_cubes_lowcloud20_aerosol_species_sixs_cci3_full_t4_20260712"
PRIMARY_METHODS = ("lut_continental", "sixs_continental", "sixs_cci3")
SMOKE_METHODS = ("sixs_cci_exact_smoke", "libradtran_continental_smoke")
SMOKE_COHORT = Path(__file__).with_name("aerosol_species_smoke_mids.txt")
_DATE_RE = re.compile(r"_(\d{8})T\d{6}$")


def _load_records(directory: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.glob("*.json")):
        try:
            row = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        records[str(row.get("matchup_id") or path.stem)] = row
    return records


def _retrieved(record: dict[str, Any] | None) -> float | None:
    if record is None or str(record.get("status", "")).upper() != "OK":
        return None
    try:
        value = float(record["retrieved"])
    except (KeyError, TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _within_ee(retrieved: float | None, truth: float) -> bool:
    return retrieved is not None and abs(retrieved - truth) <= 0.05 + 0.15 * truth + 1e-12


def summarize(
    records: dict[str, dict[str, Any]],
    cohort: list[str],
    truths: dict[str, float],
) -> dict[str, Any]:
    status_counts = Counter(
        str(records[mid].get("status", "MISSING")).upper() if mid in records else "MISSING"
        for mid in cohort
    )
    valid_rows: list[tuple[float, float]] = []
    hit_ids: list[str] = []
    runtimes: list[float] = []
    for mid in cohort:
        result = _retrieved(records.get(mid))
        truth = truths[mid]
        if result is not None:
            valid_rows.append((result, truth))
        if _within_ee(result, truth):
            hit_ids.append(mid)
        runtime = records.get(mid, {}).get("runtime_s")
        try:
            runtime_value = float(runtime)
        except (TypeError, ValueError):
            continue
        if math.isfinite(runtime_value):
            runtimes.append(runtime_value)

    errors = np.asarray([value - truth for value, truth in valid_rows], dtype=np.float64)
    valid = len(valid_rows)
    return {
        "expected": len(cohort),
        "present": sum(mid in records for mid in cohort),
        "valid": valid,
        "hits": len(hit_ids),
        "strict_rate": len(hit_ids) / len(cohort) if cohort else 0.0,
        "valid_hit_rate": len(hit_ids) / valid if valid else 0.0,
        "rmse": float(np.sqrt(np.mean(np.square(errors)))) if errors.size else None,
        "mae": float(np.mean(np.abs(errors))) if errors.size else None,
        "bias": float(np.mean(errors)) if errors.size else None,
        "median_runtime_s": float(np.median(runtimes)) if runtimes else None,
        "statuses": dict(sorted(status_counts.items())),
        "hit_matchup_ids": hit_ids,
    }


def paired_comparison(
    baseline: dict[str, dict[str, Any]],
    candidate: dict[str, dict[str, Any]],
    cohort: list[str],
    truths: dict[str, float],
) -> dict[str, Any]:
    baseline_hits: set[str] = set()
    candidate_hits: set[str] = set()
    common_valid: list[str] = []
    common_baseline_hits: set[str] = set()
    common_candidate_hits: set[str] = set()
    abs_error_delta: list[float] = []
    for mid in cohort:
        truth = truths[mid]
        base_value = _retrieved(baseline.get(mid))
        candidate_value = _retrieved(candidate.get(mid))
        if _within_ee(base_value, truth):
            baseline_hits.add(mid)
        if _within_ee(candidate_value, truth):
            candidate_hits.add(mid)
        if base_value is not None and candidate_value is not None:
            common_valid.append(mid)
            abs_error_delta.append(abs(candidate_value - truth) - abs(base_value - truth))
            if _within_ee(base_value, truth):
                common_baseline_hits.add(mid)
            if _within_ee(candidate_value, truth):
                common_candidate_hits.add(mid)

    gains = sorted(candidate_hits - baseline_hits)
    losses = sorted(baseline_hits - candidate_hits)
    discordant = len(gains) + len(losses)
    if discordant:
        from scipy.stats import binomtest

        mcnemar_p = float(binomtest(min(len(gains), len(losses)), discordant, 0.5).pvalue)
    else:
        mcnemar_p = 1.0
    deltas = np.asarray(abs_error_delta, dtype=np.float64)
    common_gains = sorted(common_candidate_hits - common_baseline_hits)
    common_losses = sorted(common_baseline_hits - common_candidate_hits)
    return {
        "expected": len(cohort),
        "baseline_hits": len(baseline_hits),
        "candidate_hits": len(candidate_hits),
        "gains": len(gains),
        "losses": len(losses),
        "net_hits": len(gains) - len(losses),
        "mcnemar_exact_p": mcnemar_p,
        "common_valid": len(common_valid),
        "common_valid_gains": len(common_gains),
        "common_valid_losses": len(common_losses),
        "common_valid_net_hits": len(common_gains) - len(common_losses),
        "common_valid_gain_matchup_ids": common_gains,
        "common_valid_loss_matchup_ids": common_losses,
        "mean_abs_error_delta": float(np.mean(deltas)) if deltas.size else None,
        "median_abs_error_delta": float(np.median(deltas)) if deltas.size else None,
        "gain_matchup_ids": gains,
        "loss_matchup_ids": losses,
    }


def _month(matchup_id: str) -> int:
    match = _DATE_RE.search(matchup_id)
    if match is None:
        raise ValueError(f"Cannot derive month from matchup_id {matchup_id!r}")
    return int(match.group(1)[4:6])


def replay_species_selection(
    matchup_id: str,
    record: dict[str, Any],
    cost_dir: Path,
) -> dict[str, Any] | None:
    from siac._rust_compat import surface_driven_pool_argmin
    from siac.algorithms.rt.aerosol_species import candidate_fraction_sets

    paths = sorted(cost_dir.glob(f"{matchup_id}.species*.npz"))
    if not paths:
        return None
    maps: list[tuple[np.ndarray, np.ndarray]] = []
    for path in paths:
        with np.load(path, allow_pickle=False) as cost:
            aod, _unc, obs_cost = surface_driven_pool_argmin(
                np.ascontiguousarray(cost["cube"], dtype=np.float64),
                np.ascontiguousarray(cost["aot_axis"], dtype=np.float64),
                np.ascontiguousarray(cost["aot_prior"], dtype=np.float64),
                np.ascontiguousarray(cost["aot_prior_unc"], dtype=np.float64),
                np.ascontiguousarray(cost["solve_valid"], dtype=bool),
                int(np.asarray(cost["pool_window"]).item()),
                int(np.asarray(cost["min_count"]).item()),
            )
        maps.append((np.asarray(aod, dtype=np.float64), np.asarray(obs_cost, dtype=np.float64)))

    shape = maps[0][0].shape
    selected = np.full(shape, -1, dtype=np.int16)
    selected_cost = np.full(shape, np.inf, dtype=np.float64)
    selected_aod = np.full(shape, np.nan, dtype=np.float64)
    for index, (aod, obs_cost) in enumerate(maps):
        better = np.isfinite(obs_cost) & (obs_cost < selected_cost)
        selected[better] = index
        selected_cost[better] = obs_cost[better]
        selected_aod[better] = aod[better]

    valid = selected >= 0
    mixtures = candidate_fraction_sets(
        float(record["lon"]),
        float(record["lat"]),
        _month(matchup_id),
        n=len(paths),
    )
    candidates = []
    for index, ((aod, obs_cost), mixture) in enumerate(zip(maps, mixtures, strict=True)):
        count = int(np.count_nonzero(selected == index))
        candidates.append(
            {
                "candidate_index": index,
                "selected_pixel_count": count,
                "selected_pixel_fraction": count / int(np.count_nonzero(valid))
                if np.any(valid)
                else 0.0,
                "candidate_scene_mean_aod": float(np.nanmean(aod)),
                "candidate_median_surface_cost": float(np.nanmedian(obs_cost)),
                "mixture": mixture,
            }
        )
    reported = _retrieved(record)
    replayed = float(np.nanmean(selected_aod)) if np.any(np.isfinite(selected_aod)) else None
    return {
        "candidate_count": len(paths),
        "selected_pixel_count": int(np.count_nonzero(valid)),
        "replayed_scene_mean_aod": replayed,
        "reported_scene_mean_aod": reported,
        "replay_delta": replayed - reported
        if replayed is not None and reported is not None
        else None,
        "candidates": candidates,
    }


def _case_rows(
    cohort: list[str],
    truths: dict[str, float],
    methods: dict[str, dict[str, dict[str, Any]]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for mid in cohort:
        row: dict[str, Any] = {"matchup_id": mid, "truth_aod": truths[mid]}
        for method, records in methods.items():
            record = records.get(mid)
            value = _retrieved(record)
            row[f"{method}_status"] = str(record.get("status", "MISSING")) if record else "MISSING"
            row[f"{method}_aod"] = value
            row[f"{method}_within_ee"] = _within_ee(value, truths[mid])
        rows.append(row)
    return rows


def _value_summary(values: list[float], truths: list[float]) -> dict[str, Any]:
    retrieved = np.asarray(values, dtype=np.float64)
    reference = np.asarray(truths, dtype=np.float64)
    errors = retrieved - reference
    thresholds = 0.05 + 0.15 * reference
    hits = int(np.count_nonzero(np.abs(errors) <= thresholds + 1e-12))
    return {
        "count": int(retrieved.size),
        "hits": hits,
        "rate": hits / retrieved.size if retrieved.size else None,
        "rmse": float(np.sqrt(np.mean(np.square(errors)))) if errors.size else None,
        "mae": float(np.mean(np.abs(errors))) if errors.size else None,
        "bias": float(np.mean(errors)) if errors.size else None,
    }


def selection_policy_replays(
    species_replay: dict[str, dict[str, Any]],
    records: dict[str, dict[str, dict[str, Any]]],
    cohort: list[str],
    truths: dict[str, float],
) -> dict[str, Any]:
    values: dict[str, list[float]] = {
        "nearest_cci": [],
        "scene_min_surface_cost": [],
        "pixel_min_surface_cost": [],
        "cci_candidate_oracle": [],
        "fixed_plus_cci_oracle": [],
    }
    references: list[float] = []
    for matchup_id in cohort:
        replay = species_replay.get(matchup_id)
        fixed = _retrieved(records["sixs_continental"].get(matchup_id))
        pixel = _retrieved(records["sixs_cci3"].get(matchup_id))
        if replay is None or fixed is None or pixel is None:
            continue
        truth = truths[matchup_id]
        candidates = replay["candidates"]
        candidate_aod = [float(row["candidate_scene_mean_aod"]) for row in candidates]
        scene_min = min(candidates, key=lambda row: row["candidate_median_surface_cost"])
        values["nearest_cci"].append(candidate_aod[0])
        values["scene_min_surface_cost"].append(float(scene_min["candidate_scene_mean_aod"]))
        values["pixel_min_surface_cost"].append(pixel)
        values["cci_candidate_oracle"].append(
            min(candidate_aod, key=lambda value: abs(value - truth))
        )
        values["fixed_plus_cci_oracle"].append(
            min([fixed, *candidate_aod], key=lambda value: abs(value - truth))
        )
        references.append(truth)
    operational = {
        "nearest_cci": True,
        "scene_min_surface_cost": True,
        "pixel_min_surface_cost": True,
        "cci_candidate_oracle": False,
        "fixed_plus_cci_oracle": False,
    }
    return {
        key: {
            "operational_without_aeronet": operational[key],
            **_value_summary(policy_values, references),
        }
        for key, policy_values in values.items()
    }


def analyze(
    *,
    root: Path,
    cohort_path: Path,
    output: Path,
    require_complete: bool,
    replay_species: bool,
) -> dict[str, Any]:
    cohort = [
        line.strip()
        for line in cohort_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    smoke_cohort = [
        line.strip()
        for line in SMOKE_COHORT.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    method_dirs = {key: root / path.name for key, path in METHODS.items()}
    records = {key: _load_records(path) for key, path in method_dirs.items()}
    baseline = records["lut_continental"]
    missing_truth = [
        mid for mid in cohort if mid not in baseline or baseline[mid].get("truth") is None
    ]
    if missing_truth:
        raise ValueError(f"Baseline truth is missing for {len(missing_truth)} matchups")
    truths = {mid: float(baseline[mid]["truth"]) for mid in cohort}

    summaries = {
        key: summarize(rows, smoke_cohort if key in SMOKE_METHODS else cohort, truths)
        for key, rows in records.items()
    }
    if require_complete:
        incomplete = {
            key: summaries[key]["present"]
            for key in PRIMARY_METHODS
            if summaries[key]["present"] != len(cohort)
        }
        if incomplete:
            raise ValueError(f"Primary experiment outputs are incomplete: {incomplete}")

    comparisons = {
        "sixs_continental_vs_lut": paired_comparison(
            baseline, records["sixs_continental"], cohort, truths
        ),
        "sixs_cci3_vs_lut": paired_comparison(baseline, records["sixs_cci3"], cohort, truths),
        "sixs_cci3_vs_sixs_continental": paired_comparison(
            records["sixs_continental"], records["sixs_cci3"], cohort, truths
        ),
    }
    species_replay: dict[str, Any] = {}
    if replay_species:
        cci_costs = root / CCI_COSTS.name
        for index, mid in enumerate(cohort, start=1):
            record = records["sixs_cci3"].get(mid)
            if _retrieved(record) is not None and record is not None:
                replay = replay_species_selection(mid, record, cci_costs)
                if replay is not None:
                    species_replay[mid] = replay
            if index % 10 == 0:
                print(f"species replay {index}/{len(cohort)}")
    policy_replays = (
        selection_policy_replays(species_replay, records, cohort, truths) if species_replay else {}
    )

    output.mkdir(parents=True, exist_ok=True)
    cases = _case_rows(cohort, truths, records)
    if cases:
        with (output / "case-results.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(cases[0]))
            writer.writeheader()
            writer.writerows(cases)
    result = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cohort_count": len(cohort),
        "smoke_cohort_count": len(smoke_cohort),
        "expected_error": "abs(retrieved - truth) <= 0.05 + 0.15 * truth",
        "method_directories": {key: str(path) for key, path in method_dirs.items()},
        "summaries": summaries,
        "comparisons": comparisons,
        "species_replay": species_replay,
        "selection_policy_replays": policy_replays,
        "case_results_csv": str(output / "case-results.csv"),
    }
    (output / "analysis.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--cohort", type=Path, default=COHORT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument("--replay-species", action="store_true")
    args = parser.parse_args()
    result = analyze(
        root=args.root,
        cohort_path=args.cohort,
        output=args.output,
        require_complete=args.require_complete,
        replay_species=args.replay_species,
    )
    print(
        json.dumps(
            {"summaries": result["summaries"], "comparisons": result["comparisons"]}, indent=2
        )
    )


if __name__ == "__main__":
    main()
