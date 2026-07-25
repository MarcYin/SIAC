from __future__ import annotations

from tools.aeronet_validation.analyze_aerosol_species_experiment import (
    paired_comparison,
    selection_policy_replays,
    summarize,
)


def _record(value: float, truth: float, *, status: str = "OK") -> dict[str, object]:
    return {"status": status, "retrieved": value, "truth": truth, "runtime_s": 12.0}


def test_summary_uses_fixed_denominator_and_non_ok_as_miss() -> None:
    cohort = ["a", "b", "c"]
    truths = {"a": 0.2, "b": 0.4, "c": 0.6}
    records = {
        "a": _record(0.2, 0.2),
        "b": _record(0.9, 0.4),
        "c": _record(0.6, 0.6, status="FAILED"),
    }
    result = summarize(records, cohort, truths)
    assert result["valid"] == 2
    assert result["hits"] == 1
    assert result["strict_rate"] == 1 / 3
    assert result["statuses"] == {"FAILED": 1, "OK": 2}


def test_paired_comparison_counts_gains_and_losses() -> None:
    cohort = ["a", "b", "c"]
    truths = {"a": 0.2, "b": 0.4, "c": 0.6}
    baseline = {
        "a": _record(0.2, 0.2),
        "b": _record(0.9, 0.4),
        "c": _record(0.6, 0.6),
    }
    candidate = {
        "a": _record(0.8, 0.2),
        "b": _record(0.4, 0.4),
        "c": _record(0.6, 0.6),
    }
    result = paired_comparison(baseline, candidate, cohort, truths)
    assert result["gains"] == 1
    assert result["losses"] == 1
    assert result["net_hits"] == 0
    assert result["gain_matchup_ids"] == ["b"]
    assert result["loss_matchup_ids"] == ["a"]


def test_selection_policy_replays_separates_operational_and_oracle() -> None:
    cohort = ["a"]
    truths = {"a": 0.4}
    records = {
        "sixs_continental": {"a": _record(0.8, 0.4)},
        "sixs_cci3": {"a": _record(0.5, 0.4)},
    }
    replay = {
        "a": {
            "candidates": [
                {"candidate_scene_mean_aod": 0.7, "candidate_median_surface_cost": 2.0},
                {"candidate_scene_mean_aod": 0.42, "candidate_median_surface_cost": 1.0},
                {"candidate_scene_mean_aod": 0.9, "candidate_median_surface_cost": 3.0},
            ]
        }
    }

    result = selection_policy_replays(replay, records, cohort, truths)

    assert result["scene_min_surface_cost"]["hits"] == 1
    assert result["cci_candidate_oracle"]["hits"] == 1
    assert not result["cci_candidate_oracle"]["operational_without_aeronet"]


def test_paired_comparison_separates_strict_missing_from_common_valid() -> None:
    result = paired_comparison(
        {"a": _record(0.2, 0.2), "b": _record(0.4, 0.4)},
        {"a": _record(0.2, 0.2)},
        ["a", "b"],
        {"a": 0.2, "b": 0.4},
    )

    assert result["losses"] == 1
    assert result["common_valid"] == 1
    assert result["common_valid_losses"] == 0
