"""Tests for the local Cloud Score+ index policy."""

from __future__ import annotations

import numpy as np
import pytest
from tools.aeronet_validation.cloudscore_index_policy import (
    clean_coverage,
    daily_weights,
    index_policy,
    seasonal_windows,
    select_candidate_days,
)


def test_seasonal_windows_centre_on_the_calendar_month_midpoint() -> None:
    # December's midpoint is the 15th, so the window crosses the year boundary.
    windows = seasonal_windows("2020-12-24", library_years=(2019, 2020))
    assert windows == (("2019-10-31", "2020-01-29"), ("2020-10-31", "2021-01-29"))


def test_seasonal_windows_span_the_declared_width() -> None:
    ((start, end),) = seasonal_windows("2021-06-02", library_years=(2021,))
    import datetime as dt

    assert (dt.date.fromisoformat(end) - dt.date.fromisoformat(start)).days == 90


def test_clean_coverage_counts_pixels_at_or_above_threshold() -> None:
    plane = np.array([0.2, 0.6, 0.9, np.nan])
    # 0.6 and 0.9 qualify out of three finite pixels.
    assert clean_coverage(plane) == pytest.approx(2 / 3)


def test_clean_coverage_of_an_empty_plane_is_zero() -> None:
    assert clean_coverage(np.full(4, np.nan)) == 0.0


def test_candidate_rule_compares_against_the_calendar_month_mean() -> None:
    # Two Junes and two Decembers; December is cloudier and must be judged
    # against its own norm rather than against June.
    coverage = {
        "2019-06-10": 0.9,
        "2020-06-10": 0.5,  # below the June mean of 0.7
        "2019-12-10": 0.3,
        "2020-12-10": 0.1,  # below the December mean of 0.2
    }
    assert select_candidate_days(coverage) == ("2019-06-10", "2019-12-10")


def test_candidate_rule_drops_zero_coverage_days() -> None:
    assert select_candidate_days({"2019-06-10": 0.0, "2019-06-20": 0.0}) == ()


def test_daily_weights_substitute_the_locked_value_for_maiac_gaps() -> None:
    weights = daily_weights(
        {"2019-06-10": 0.2, "2019-06-20": None},
        {"2019-06-10": 0.9, "2019-06-20": 0.8},
    )
    assert weights["2019-06-10"]["aod_source"] == "maiac"
    assert weights["2019-06-20"]["aod_source"] == "locked_constant_0p1"
    assert weights["2019-06-20"]["aod"] == pytest.approx(0.1)


def test_daily_weight_rises_with_clean_coverage() -> None:
    weights = daily_weights(
        {"2019-06-10": 0.2, "2019-06-20": 0.2},
        {"2019-06-10": 0.2, "2019-06-20": 0.95},
    )
    assert weights["2019-06-20"]["weight"] > weights["2019-06-10"]["weight"]


def test_policy_does_not_claim_the_earth_engine_winner_source() -> None:
    policy = index_policy()
    # The ordering matches the committed index but the argmax is local; an
    # archive must not misattribute where its winners came from.
    assert policy["winner_source"] == "local_mosaic_from_edown_cloud_score_plus"
    assert policy["schema"] == "siac_l1c_cloudscore_winner_index_v1"
    assert policy["quality_order"] == (
        "normalized_cs + locked_daily_aod_weight + aoi_overlap_ratio"
    )
    assert policy["cloud_score_clear_threshold"] == 0.6


def test_aod_weight_uses_raw_maiac_integer_units() -> None:
    # Physical-scale AOD collapses the sigmoid to a near-constant ~0.5, making
    # the aerosol term inert. Raw units must produce a real spread.
    from tools.aeronet_validation.cloudscore_index_policy import MAIAC_RAW_SCALE

    assert MAIAC_RAW_SCALE == 1000.0
    coverage = {f"2019-06-{d:02d}": 0.9 for d in range(1, 11)}
    aod = {day: 0.05 + 0.08 * i for i, day in enumerate(sorted(coverage))}
    weights = daily_weights(aod, coverage)
    spread = max(w["weight"] for w in weights.values()) - min(w["weight"] for w in weights.values())
    assert spread > 0.3


def test_default_band_is_cs_not_the_cdf_form() -> None:
    from tools.aeronet_validation.cloudscore_index_policy import CLOUD_SCORE_BAND

    # cs is the direct clearness estimate and the more aggressive masker; the
    # committed index used cs_cdf, this pipeline deliberately does not.
    assert CLOUD_SCORE_BAND == "cs"
    policy = index_policy()
    assert policy["cloud_score_band"] == "cs"
    assert policy["quality_order"].startswith("normalized_cs +")
