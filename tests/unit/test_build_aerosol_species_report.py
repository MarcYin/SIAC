from __future__ import annotations

import pytest
from tools.aeronet_validation.build_aerosol_species_report import (
    _compact_mixture_label,
    _final_pass_policy_replays,
    _transition,
    _wilson_interval,
)


@pytest.mark.parametrize(
    ("baseline", "candidate", "expected"),
    [
        (False, True, "gain"),
        (True, False, "loss"),
        (True, True, "both_hit"),
        (False, False, "both_miss"),
    ],
)
def test_transition_labels_all_paired_outcomes(
    baseline: bool,
    candidate: bool,
    expected: str,
) -> None:
    assert _transition(baseline, candidate) == expected


def test_wilson_interval_contains_observed_rate() -> None:
    low, high = _wilson_interval(110, 152) or (0.0, 0.0)

    assert low < 110 / 152 < high
    assert 0.0 <= low < high <= 1.0


def test_compact_mixture_label_is_stable() -> None:
    label = _compact_mixture_label(
        {"dust": 0.25, "sea_salt": 0.25, "fine_strong": 0.0, "fine_weak": 0.5}
    )

    assert label == "D25 SS25 FS0 FW50"


def test_final_pass_policy_replays_keep_oracle_non_operational() -> None:
    cases = [
        {
            "truth_aod": 0.4,
            "methods": {"sixs_continental": {"aod": 0.8}, "sixs_cci3": {"aod": 0.5}},
            "species": {
                "candidates": [
                    {"scene_mean_aod": 0.7, "median_surface_cost": 2.0},
                    {"scene_mean_aod": 0.42, "median_surface_cost": 1.0},
                    {"scene_mean_aod": 0.9, "median_surface_cost": 3.0},
                ]
            },
        }
    ]

    result = _final_pass_policy_replays(cases)

    assert result["scene_min_surface_cost"]["hits"] == 1
    assert result["cci_candidate_oracle"]["hits"] == 1
    assert not result["cci_candidate_oracle"]["operational_without_aeronet"]
