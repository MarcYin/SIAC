"""Tests for global current-date scene selection."""

from __future__ import annotations

import pytest
from tools.aeronet_validation.global_scene_selection import (
    SceneOption,
    assign_target_months,
    choose_scene,
    matchup_id,
    season_balance,
)


def _option(date: str, cloud: float = 5.0, tile: str = "T31UDQ", pid: str | None = None):
    return SceneOption(
        product_id=pid or f"S2A_{tile}_{date}",
        datetime=f"{date}T10:30:00Z",
        mgrs_tile=tile,
        cloud_cover=cloud,
    )


def test_target_months_are_evenly_spread() -> None:
    assigned = assign_target_months([f"s{i:04d}" for i in range(120)])
    counts = {m: sum(1 for v in assigned.values() if v == m) for m in range(1, 13)}
    assert set(counts.values()) == {10}


def test_target_months_are_deterministic_and_seed_sensitive() -> None:
    ids = [f"s{i:03d}" for i in range(48)]
    assert assign_target_months(ids, seed=3) == assign_target_months(ids, seed=3)
    assert assign_target_months(ids, seed=3) != assign_target_months(ids, seed=4)


def test_target_month_is_preferred_over_a_clearer_other_month() -> None:
    # Season balance must not be traded away for a marginally clearer scene.
    chosen = choose_scene(
        [_option("2020-06-15", cloud=0.0), _option("2020-03-10", cloud=15.0)],
        target_month=3,
    )
    assert chosen is not None and chosen.month == 3


def test_clearest_wins_within_the_target_month() -> None:
    chosen = choose_scene(
        [_option("2020-03-05", cloud=18.0), _option("2020-03-20", cloud=2.0)], target_month=3
    )
    assert chosen is not None and chosen.cloud_cover == 2.0


def test_nearest_month_is_used_when_the_target_is_unavailable() -> None:
    # A polar AOI with no usable winter acquisitions must stay in the catalogue.
    chosen = choose_scene([_option("2020-07-15"), _option("2020-11-15")], target_month=12)
    assert chosen is not None and chosen.month == 11


def test_month_distance_wraps_across_the_year_boundary() -> None:
    chosen = choose_scene([_option("2020-06-15"), _option("2020-12-15")], target_month=1)
    assert chosen is not None and chosen.month == 12


def test_cloudy_and_out_of_range_scenes_are_rejected() -> None:
    assert choose_scene([_option("2020-03-10", cloud=55.0)], target_month=3) is None
    assert choose_scene([_option("2017-03-10", cloud=1.0)], target_month=3) is None


def test_no_usable_option_returns_none_rather_than_raising() -> None:
    assert choose_scene([], target_month=6) is None


def test_invalid_target_month_is_rejected() -> None:
    with pytest.raises(ValueError, match="target_month"):
        choose_scene([_option("2020-03-10")], target_month=13)


def test_matchup_id_matches_the_existing_convention() -> None:
    scene = _option("2020-12-24", tile="T37PDL")
    assert matchup_id("t0042_bare_sparse_00120_00311", scene).endswith("__T37PDL_20201224T103000")


def test_season_balance_reports_realised_shares() -> None:
    selected = {"a": _option("2020-01-05"), "b": _option("2020-01-09"), "c": _option("2020-07-09")}
    shares = season_balance(selected)
    assert shares[1] == pytest.approx(2 / 3)
    assert shares[7] == pytest.approx(1 / 3)
    assert sum(shares.values()) == pytest.approx(1.0)
