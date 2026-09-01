"""Tests for the local Cloud Score+ quality-mosaic composer."""

from __future__ import annotations

import numpy as np
import pytest
from tools.aeronet_validation.cloudscore_local_mosaic import (
    NO_WINNER,
    Candidate,
    agreement,
    compose_winners,
    mosaic_weight,
    winner_day_plane,
)


def _candidate(day: str, cs: float | np.ndarray, *, day_weight=0.0, coverage=0.0, shape=(4, 4)):
    plane = np.full(shape, cs, dtype=np.float64) if np.isscalar(cs) else np.asarray(cs, float)
    return Candidate(day=day, score=plane, day_weight=day_weight, coverage_ratio=coverage)


def test_weight_matches_the_committed_formula() -> None:
    # (cs_cdf - 0.6)/0.4 + day_weight + coverage_ratio
    got = mosaic_weight(np.array([0.6, 1.0]), day_weight=0.25, coverage_ratio=0.5)
    np.testing.assert_allclose(got, [0.0 + 0.75, 1.0 + 0.75])


def test_clear_span_must_be_positive() -> None:
    with pytest.raises(ValueError, match="clear_span"):
        mosaic_weight(np.zeros(3), day_weight=0.0, coverage_ratio=0.0, clear_span=0.0)


def test_highest_quality_acquisition_wins() -> None:
    winners, days = compose_winners(
        [
            _candidate("2024-01-01", 0.7),
            _candidate("2024-02-01", 0.95),
            _candidate("2024-03-01", 0.8),
        ]
    )
    assert days == ("2024-01-01", "2024-02-01", "2024-03-01")
    assert np.all(winners == 1)


def test_day_weight_can_outrank_a_better_pixel_score() -> None:
    # A cleaner pixel on a day with a poor AOD/coverage weight must lose to a
    # slightly less clear pixel on a strongly weighted day.
    winners, _ = compose_winners(
        [
            _candidate("2024-01-01", 1.0, day_weight=0.0),
            _candidate("2024-02-01", 0.9, day_weight=1.0),
        ]
    )
    assert np.all(winners == 1)


def test_ties_follow_earth_engine_mosaic_semantics() -> None:
    # cs_cdf saturates at 1.0 over confidently clear pixels, so ties are common.
    # qualityMosaic keeps the LAST maximum; numpy argmax would keep the first.
    tied = [_candidate("2024-01-01", 1.0), _candidate("2024-02-01", 1.0)]
    last, _ = compose_winners(tied)
    first, _ = compose_winners(tied, tie_breaking="first")
    assert np.all(last == 1)
    assert np.all(first == 0)


def test_rejects_unknown_tie_breaking() -> None:
    with pytest.raises(ValueError, match="tie_breaking"):
        compose_winners([_candidate("2024-01-01", 1.0)], tie_breaking="random")


def test_pixels_with_no_usable_acquisition_get_no_winner() -> None:
    a = _candidate("2024-01-01", np.array([[np.nan, 0.9], [0.9, np.nan]]))
    b = _candidate("2024-02-01", np.array([[np.nan, 0.7], [0.5, np.nan]]))
    winners, _ = compose_winners([a, b])
    assert winners[0, 0] == NO_WINNER
    assert winners[1, 1] == NO_WINNER
    assert winners[0, 1] == 0
    assert winners[1, 0] == 0


def test_mismatched_grids_are_rejected() -> None:
    with pytest.raises(ValueError, match="one grid"):
        compose_winners(
            [_candidate("2024-01-01", 0.9), _candidate("2024-02-01", 0.9, shape=(3, 3))]
        )


def test_empty_candidate_list_is_rejected() -> None:
    with pytest.raises(ValueError, match="no candidate"):
        compose_winners([])


def test_winner_day_plane_maps_indices_to_day_codes() -> None:
    winners = np.array([[0, 1], [NO_WINNER, 1]], dtype=np.int32)
    plane = winner_day_plane(winners, ["2024-01-05", "2024-02-17"])
    assert plane[0, 0] == 20240105
    assert plane[0, 1] == 20240217
    assert plane[1, 0] == -1


def test_agreement_reports_perfect_and_partial_matches() -> None:
    ref = np.array([[0, 1], [2, NO_WINNER]])
    assert agreement(ref, ref.copy())["identical_fraction"] == pytest.approx(1.0)
    other = np.array([[0, 1], [1, NO_WINNER]])
    stats = agreement(ref, other)
    assert stats["identical_fraction"] == pytest.approx(0.75)
    assert stats["coverage_disagreement"] == pytest.approx(0.0)


def test_erosion_removes_cloud_adjacent_pixels() -> None:
    from tools.aeronet_validation.cloudscore_local_mosaic import erode_valid

    valid = np.ones((7, 7), dtype=bool)
    valid[3, 3] = False
    eroded = erode_valid(valid, radius_px=1)
    # The invalid pixel and its 4-neighbours go, and the border is dropped.
    assert not eroded[3, 3]
    assert not eroded[2, 3] and not eroded[4, 3]
    assert not eroded[0, :].any() and not eroded[-1, :].any()
    assert eroded[1, 1]


def test_erosion_is_a_noop_at_zero_radius() -> None:
    from tools.aeronet_validation.cloudscore_local_mosaic import erode_valid

    valid = np.array([[True, False], [True, True]])
    np.testing.assert_array_equal(erode_valid(valid, radius_px=0), valid)


def test_monthly_composition_indexes_within_each_month() -> None:
    from tools.aeronet_validation.cloudscore_local_mosaic import compose_monthly_winners

    months, winners, ordering = compose_monthly_winners(
        {
            "2019-02": [_candidate("2019-02-03", 0.7), _candidate("2019-02-18", 0.95)],
            "2019-01": [_candidate("2019-01-09", 0.9)],
        }
    )
    # Months are ordered, and indices are local to each month's candidate list.
    assert months == ("2019-01", "2019-02")
    assert winners.shape == (2, 4, 4)
    assert winners.dtype == np.int16
    assert np.all(winners[0] == 0)
    assert np.all(winners[1] == 1)
    assert ordering["2019-02"] == ("2019-02-03", "2019-02-18")


def test_monthly_composition_rejects_mismatched_grids() -> None:
    from tools.aeronet_validation.cloudscore_local_mosaic import compose_monthly_winners

    with pytest.raises(ValueError, match="one grid"):
        compose_monthly_winners(
            {
                "2019-01": [_candidate("2019-01-09", 0.9)],
                "2019-02": [_candidate("2019-02-09", 0.9, shape=(5, 5))],
            }
        )


def test_monthly_composition_rejects_empty_input() -> None:
    from tools.aeronet_validation.cloudscore_local_mosaic import compose_monthly_winners

    with pytest.raises(ValueError, match="no months"):
        compose_monthly_winners({})
