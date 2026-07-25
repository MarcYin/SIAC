"""Select and audit one global shifted-log AOD offset.

The base retrieval calibrator is fixed before this analysis. Target folds 0-3
select one scalar offset, while fold 4 remains untouched confirmation. The
same offset is then replayed across every seed/tree variant and the external
holdout. No scene, site, geography, or AERONET value is used at inference.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from tools.aeronet_validation.aod_residual_calibration import expected_error, metrics

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

SHIFT = 1.0 / 3.0
TARGET_CONFIRMATION_FOLD = 4
TARGET_HITS = 133
OFFSETS = tuple(step / 400.0 for step in range(-20, 21))


def apply_shifted_log_offset(prediction: np.ndarray, offset: float) -> np.ndarray:
    """Apply a scalar correction in the shifted-log AOD domain."""
    values = np.asarray(prediction, dtype=float)
    return np.clip((values + SHIFT) * np.exp(offset) - SHIFT, 0.0, 4.0)


def evaluate_prediction(
    rows: Sequence[dict[str, Any]],
    prediction: np.ndarray,
    indices: np.ndarray,
    *,
    include_predictions: bool = False,
) -> dict[str, Any]:
    """Evaluate one prediction vector on explicitly supplied row indices."""
    truth_all = np.asarray([row["truth"] for row in rows], dtype=float)
    baseline_all = np.asarray([row["baseline"] for row in rows], dtype=float)
    truth = truth_all[indices]
    baseline = baseline_all[indices]
    candidate = np.asarray(prediction, dtype=float)[indices]
    baseline_hit = np.abs(baseline - truth) <= expected_error(truth)
    candidate_hit = np.abs(candidate - truth) <= expected_error(truth)
    result: dict[str, Any] = {
        "count": int(len(indices)),
        "baseline": metrics(truth, baseline),
        "candidate": metrics(truth, candidate),
        "gains": int((candidate_hit & ~baseline_hit).sum()),
        "losses": int((~candidate_hit & baseline_hit).sum()),
    }
    if include_predictions:
        result["predictions"] = [
            {
                "matchup_id": rows[index]["matchup_id"],
                "site": rows[index]["site"],
                "site_fold": int(rows[index]["site_fold"]),
                "truth": float(truth_all[index]),
                "baseline": float(baseline_all[index]),
                "candidate": float(prediction[index]),
                "baseline_within_ee": bool(baseline_hit[position]),
                "candidate_within_ee": bool(candidate_hit[position]),
            }
            for position, index in enumerate(indices)
        ]
    return result


def rank_offsets(
    rows: Sequence[dict[str, Any]],
    prediction: np.ndarray,
    indices: np.ndarray,
    offsets: Iterable[float] = OFFSETS,
) -> list[dict[str, Any]]:
    """Rank scalar offsets using only the requested development indices."""
    ranking: list[dict[str, Any]] = []
    for offset in offsets:
        adjusted = apply_shifted_log_offset(prediction, offset)
        summary = evaluate_prediction(rows, adjusted, indices)["candidate"]
        ranking.append(
            {
                "offset": float(offset),
                "metrics": summary,
                "score": [
                    int(summary["hits"]),
                    -float(summary["mae"]),
                    -float(summary["rmse"]),
                    -abs(float(offset)),
                    -float(offset),
                ],
            }
        )
    ranking.sort(key=lambda row: tuple(row["score"]), reverse=True)
    return ranking


def _site_from_matchup_id(matchup_id: str) -> str:
    return matchup_id.split("__", 1)[0]


def _external_sites(path: Path) -> set[str]:
    return {
        _site_from_matchup_id(line.strip())
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


def _assert_ids(label: str, expected: Sequence[str], actual: Sequence[str]) -> None:
    if list(expected) != list(actual):
        raise ValueError(f"{label} matchup IDs or order do not agree")


def _case_receipts(
    rows: Sequence[dict[str, Any]],
    unadjusted: np.ndarray,
    adjusted: np.ndarray,
    seed_predictions: Sequence[np.ndarray],
    seen_sites: set[str],
) -> list[dict[str, Any]]:
    truth = np.asarray([row["truth"] for row in rows], dtype=float)
    baseline = np.asarray([row["baseline"] for row in rows], dtype=float)
    ee = expected_error(truth)
    baseline_hit = np.abs(baseline - truth) <= ee
    unadjusted_hit = np.abs(unadjusted - truth) <= ee
    adjusted_hit = np.abs(adjusted - truth) <= ee
    seed_stack = np.stack(seed_predictions)
    receipts: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        if adjusted_hit[index] and not baseline_hit[index]:
            transition = "gain"
        elif not adjusted_hit[index] and baseline_hit[index]:
            transition = "loss"
        elif adjusted_hit[index]:
            transition = "retained_hit"
        else:
            transition = "remaining_miss"
        receipts.append(
            {
                "matchup_id": row["matchup_id"],
                "site": row["site"],
                "site_fold": int(row["site_fold"]),
                "external_site_seen": bool(row["site"] in seen_sites),
                "truth": float(truth[index]),
                "expected_error": float(ee[index]),
                "baseline": float(baseline[index]),
                "unadjusted": float(unadjusted[index]),
                "adjusted": float(adjusted[index]),
                "baseline_error": float(baseline[index] - truth[index]),
                "unadjusted_error": float(unadjusted[index] - truth[index]),
                "adjusted_error": float(adjusted[index] - truth[index]),
                "baseline_error_over_ee": float(
                    (baseline[index] - truth[index]) / ee[index]
                ),
                "adjusted_error_over_ee": float(
                    (adjusted[index] - truth[index]) / ee[index]
                ),
                "baseline_within_ee": bool(baseline_hit[index]),
                "unadjusted_within_ee": bool(unadjusted_hit[index]),
                "adjusted_within_ee": bool(adjusted_hit[index]),
                "transition": transition,
                "seed300_adjusted_min": float(seed_stack[:, index].min()),
                "seed300_adjusted_max": float(seed_stack[:, index].max()),
                "seed300_adjusted_std": float(seed_stack[:, index].std()),
            }
        )
    return receipts


def _compact(summary: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in summary.items() if key != "predictions"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-analysis", type=Path, required=True)
    parser.add_argument("--seed-analysis", type=Path, required=True)
    parser.add_argument("--external-train-mids", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    selection = json.loads(args.selection_analysis.read_text(encoding="utf-8"))
    seed_analysis = json.loads(args.seed_analysis.read_text(encoding="utf-8"))
    target_rows = selection["target_all_descriptive"]["predictions"]
    external_rows = selection["external_holdout"]["predictions"]
    _assert_ids(
        "Target",
        [row["matchup_id"] for row in target_rows],
        seed_analysis["target_matchup_ids"],
    )
    _assert_ids(
        "External holdout",
        [row["matchup_id"] for row in external_rows],
        seed_analysis["external_holdout_matchup_ids"],
    )

    primary_name = seed_analysis["primary_stability_variant"]
    variants = seed_analysis["variants"]
    primary = np.asarray(variants[primary_name]["target_prediction"], dtype=float)
    folds = np.asarray([row["site_fold"] for row in target_rows], dtype=int)
    development_indices = np.flatnonzero(folds != TARGET_CONFIRMATION_FOLD)
    confirmation_indices = np.flatnonzero(folds == TARGET_CONFIRMATION_FOLD)
    all_indices = np.arange(len(target_rows))
    ranking = rank_offsets(target_rows, primary, development_indices)
    selected_offset = float(ranking[0]["offset"])
    adjusted_primary = apply_shifted_log_offset(primary, selected_offset)

    nested_prediction = np.full(len(target_rows), np.nan, dtype=float)
    nested_folds: list[dict[str, Any]] = []
    for fold in sorted(set(folds.tolist())):
        fold_development = np.flatnonzero(folds != fold)
        fold_confirmation = np.flatnonzero(folds == fold)
        fold_ranking = rank_offsets(target_rows, primary, fold_development)
        fold_offset = float(fold_ranking[0]["offset"])
        fold_prediction = apply_shifted_log_offset(primary, fold_offset)
        nested_prediction[fold_confirmation] = fold_prediction[fold_confirmation]
        nested_folds.append(
            {
                "held_out_fold": int(fold),
                "selected_offset": fold_offset,
                "development": evaluate_prediction(
                    target_rows, fold_prediction, fold_development
                ),
                "confirmation": evaluate_prediction(
                    target_rows, fold_prediction, fold_confirmation
                ),
            }
        )
    if not np.all(np.isfinite(nested_prediction)):
        raise ValueError("Nested offset predictions are incomplete")

    seen_sites = _external_sites(args.external_train_mids)
    seen = np.asarray([row["site"] in seen_sites for row in target_rows], dtype=bool)
    seed_names = sorted(
        name for name in variants if name.startswith("seed") and name.endswith("_n300")
    )
    adjusted_seed_predictions = [
        apply_shifted_log_offset(
            np.asarray(variants[name]["target_prediction"], dtype=float),
            selected_offset,
        )
        for name in seed_names
    ]

    replay: dict[str, Any] = {}
    for name, variant in variants.items():
        target_prediction = apply_shifted_log_offset(
            np.asarray(variant["target_prediction"], dtype=float), selected_offset
        )
        external_prediction = apply_shifted_log_offset(
            np.asarray(variant["external_holdout_prediction"], dtype=float),
            selected_offset,
        )
        replay[name] = {
            "target_all": evaluate_prediction(
                target_rows, target_prediction, all_indices
            ),
            "target_development": evaluate_prediction(
                target_rows, target_prediction, development_indices
            ),
            "target_confirmation": evaluate_prediction(
                target_rows, target_prediction, confirmation_indices
            ),
            "target_seen_sites": evaluate_prediction(
                target_rows, target_prediction, np.flatnonzero(seen)
            ),
            "target_unseen_sites": evaluate_prediction(
                target_rows, target_prediction, np.flatnonzero(~seen)
            ),
            "external_holdout": evaluate_prediction(
                external_rows, external_prediction, np.arange(len(external_rows))
            ),
        }

    primary_external = apply_shifted_log_offset(
        np.asarray(variants[primary_name]["external_holdout_prediction"], dtype=float),
        selected_offset,
    )
    target_final = evaluate_prediction(
        target_rows, adjusted_primary, all_indices, include_predictions=True
    )
    analysis = {
        "schema_version": "siac-global-aod-offset-v1",
        "selection_policy": {
            "base_recipe_fixed_before_offset_selection": seed_analysis["recipe"],
            "primary_prediction_variant": primary_name,
            "target_development_folds": [0, 1, 2, 3],
            "target_confirmation_fold": TARGET_CONFIRMATION_FOLD,
            "confirmation_labels_used_for_selection": False,
            "offset_grid": list(OFFSETS),
            "ranking": ["within_ee_hits", "mae", "rmse", "smallest_abs_offset"],
            "target_or_site_label_used_at_inference": False,
            "geography_or_case_routing_used_at_inference": False,
            "formula": "clip((aod + 1/3) * exp(offset) - 1/3, 0, 4)",
            "selection_scope_note": (
                "The nested replay validates scalar-offset selection for the already-fixed "
                "base recipe; recipe-family selection was audited separately."
            ),
        },
        "selected_offset": selected_offset,
        "target_hits_required": TARGET_HITS,
        "target_threshold_met": bool(target_final["candidate"]["hits"] >= TARGET_HITS),
        "development_ranking": ranking,
        "target_unadjusted": evaluate_prediction(target_rows, primary, all_indices),
        "target_development": evaluate_prediction(
            target_rows, adjusted_primary, development_indices
        ),
        "target_confirmation": evaluate_prediction(
            target_rows, adjusted_primary, confirmation_indices
        ),
        "target_all": target_final,
        "target_seen_sites": evaluate_prediction(
            target_rows, adjusted_primary, np.flatnonzero(seen)
        ),
        "target_unseen_sites": evaluate_prediction(
            target_rows, adjusted_primary, np.flatnonzero(~seen)
        ),
        "target_nested_offset_selection": {
            "aggregate": evaluate_prediction(
                target_rows, nested_prediction, all_indices, include_predictions=True
            ),
            "folds": nested_folds,
        },
        "external_holdout": evaluate_prediction(
            external_rows,
            primary_external,
            np.arange(len(external_rows)),
            include_predictions=True,
        ),
        "seed_tree_replay": replay,
        "seed300_target_hits": {
            name: replay[name]["target_all"]["candidate"]["hits"]
            for name in seed_names
        },
        "all_replayed_variants_meet_target": all(
            row["target_all"]["candidate"]["hits"] >= TARGET_HITS
            for row in replay.values()
        ),
        "case_receipts": _case_receipts(
            target_rows,
            primary,
            adjusted_primary,
            adjusted_seed_predictions,
            seen_sites,
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "selected_offset": selected_offset,
                "target_unadjusted": analysis["target_unadjusted"],
                "target_development": analysis["target_development"],
                "target_confirmation": analysis["target_confirmation"],
                "target_all": _compact(analysis["target_all"]),
                "target_nested_offset_selection": _compact(
                    analysis["target_nested_offset_selection"]["aggregate"]
                ),
                "target_seen_sites": analysis["target_seen_sites"],
                "target_unseen_sites": analysis["target_unseen_sites"],
                "external_holdout": _compact(analysis["external_holdout"]),
                "seed300_target_hits": analysis["seed300_target_hits"],
                "all_replayed_variants_meet_target": analysis[
                    "all_replayed_variants_meet_target"
                ],
            },
            indent=2,
        )
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
