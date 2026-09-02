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
    # A moderately clouded scene is now admitted on purpose: the model is asked
    # to correct whatever scene it is given, so a near-clear-only corpus would
    # not cover inference. The ceiling still bites above it.
    assert choose_scene([_option("2020-03-10", cloud=55.0)], target_month=3) is not None
    assert choose_scene([_option("2020-03-10", cloud=95.0)], target_month=3) is None
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


def test_raising_the_ceiling_alone_would_not_spread_cloud_cover() -> None:
    """The reason a target level exists at all.

    Selection minimises cloud cover, so a higher ceiling only admits scenes it
    then declines to pick. Without a target the clearest still wins, however
    high the ceiling goes.
    """

    options = [
        SceneOption("clear", "2021-06-10T00:00:00", "T31UDQ", 2.0),
        SceneOption("mid", "2021-06-15T00:00:00", "T31UDQ", 45.0),
        SceneOption("cloudy", "2021-06-20T00:00:00", "T31UDQ", 80.0),
    ]
    assert choose_scene(options, target_month=6, max_cloud=90.0).product_id == "clear"
    assert choose_scene(options, target_month=6, target_cloud=45.0).product_id == "mid"
    assert choose_scene(options, target_month=6, target_cloud=75.0).product_id == "cloudy"


def test_a_target_level_is_approached_not_enforced() -> None:
    # An AOI that is never cloudy has no 75% scene to give; refusing it would
    # drop the driest regions, so it contributes its closest scene instead.
    options = [
        SceneOption("a", "2021-06-10T00:00:00", "T31UDQ", 1.0),
        SceneOption("b", "2021-06-15T00:00:00", "T31UDQ", 6.0),
    ]
    assert choose_scene(options, target_month=6, target_cloud=75.0).product_id == "b"


def test_cloud_quota_is_exact_and_sums_to_the_catalogue() -> None:
    from tools.aeronet_validation.global_scene_selection import CLOUD_TARGETS, cloud_quota

    quota = cloud_quota(5000)
    assert sum(quota.values()) == 5000
    for level, share in CLOUD_TARGETS:
        assert quota[level] == pytest.approx(share * 5000, abs=1)
    # An awkward total must still sum exactly, not lose a scene to rounding.
    assert sum(cloud_quota(37).values()) == 37


def test_targeting_follows_the_deficit_not_a_fixed_assignment() -> None:
    """Stratification has to hold on what was *selected*, not what was wished for.

    A pre-assigned target is only a wish: an AOI with no cloudy acquisition
    returns a clear one instead, and under a fixed assignment that shortfall is
    never made up, so the cloudy strata finish short.
    """

    from tools.aeronet_validation.global_scene_selection import cloud_quota, next_cloud_target

    quota = cloud_quota(100)
    realised = dict.fromkeys(quota, 0)
    # Nothing selected yet: the largest quota leads.
    assert next_cloud_target(realised, quota) == 5.0
    # Clear stratum already satisfied, so the deficit moves on.
    realised[5.0] = quota[5.0]
    assert next_cloud_target(realised, quota) != 5.0
    # Everything satisfied except the cloudiest: that is what gets targeted.
    for level in quota:
        realised[level] = quota[level]
    realised[75.0] = 0
    assert next_cloud_target(realised, quota) == 75.0


def test_deficit_targeting_fills_every_stratum_when_scenes_allow() -> None:
    from tools.aeronet_validation.global_scene_selection import (
        cloud_quota,
        cloud_stratum,
        next_cloud_target,
    )

    quota = cloud_quota(200)
    realised = dict.fromkeys(quota, 0)
    # Every AOI can offer any cloud level, so the quota should be met exactly.
    for _ in range(200):
        target = next_cloud_target(realised, quota)
        realised[cloud_stratum(target)] += 1
    assert realised == quota


def test_cloud_stratum_assigns_each_scene_to_one_target() -> None:
    from tools.aeronet_validation.global_scene_selection import cloud_stratum

    assert cloud_stratum(0.0) == 5.0
    assert cloud_stratum(11.0) == 5.0
    assert cloud_stratum(30.0) == 20.0
    assert cloud_stratum(60.0) == 45.0
    assert cloud_stratum(88.0) == 75.0


def test_the_ceiling_still_rejects_beyond_it() -> None:
    options = [SceneOption("a", "2021-06-10T00:00:00", "T31UDQ", 95.0)]
    assert choose_scene(options, target_month=6, target_cloud=75.0) is None
