"""Tests for the free-signal prefilter that shortlists Cloud Score+ requests."""

from __future__ import annotations

import pytest
from tools.aeronet_validation.scout_prefilter import (
    ScoutRecord,
    request_saving,
    shortlist,
    shortlist_month,
)


def _record(
    day: str,
    *,
    clear: float = 1.0,
    coverage: float = 1.0,
    aod: float | None = 0.1,
    image_id: str | None = None,
) -> ScoutRecord:
    return ScoutRecord(
        image_id=image_id or f"img_{day}",
        day=day,
        clear_fraction=clear,
        coverage_fraction=coverage,
        aod=aod,
    )


def test_partial_coverage_is_rejected_outright() -> None:
    # A swath-edge acquisition cannot contribute over the missing part, so it is
    # rejected rather than merely ranked down.
    records = [_record("2020-06-01", coverage=0.5), _record("2020-06-02", coverage=1.0)]
    assert [r.day for r in shortlist_month(records, top_k=5)] == ["2020-06-02"]


def test_lowest_aerosol_wins_among_clear_acquisitions() -> None:
    records = [
        _record("2020-06-01", clear=1.00, aod=0.40),
        _record("2020-06-02", clear=0.99, aod=0.05),
        _record("2020-06-03", clear=1.00, aod=0.20),
    ]
    assert [r.day for r in shortlist_month(records, top_k=2)] == ["2020-06-02", "2020-06-03"]


def test_shortlist_is_capped_at_top_k() -> None:
    records = [_record(f"2020-06-{d:02d}", aod=0.01 * d) for d in range(1, 21)]
    assert len(shortlist_month(records, top_k=5)) == 5


def test_cloudy_month_falls_back_rather_than_returning_nothing() -> None:
    # A persistently cloudy AOI must still contribute, or the corpus silently
    # loses the regions that are hardest to sample.
    records = [_record(f"2020-06-{d:02d}", clear=0.30 + 0.05 * d) for d in range(1, 5)]
    chosen = shortlist_month(records, top_k=3)
    assert len(chosen) == 3
    # Fallback ranks by clarity first, since none are clear enough to compare on aerosol.
    assert chosen[0].clear_fraction == pytest.approx(0.50)


def test_missing_aod_is_treated_as_the_locked_constant() -> None:
    records = [_record("2020-06-01", aod=None), _record("2020-06-02", aod=0.05)]
    assert [r.day for r in shortlist_month(records, top_k=2)] == ["2020-06-02", "2020-06-01"]
    assert records[0].effective_aod == pytest.approx(0.1)


def test_all_partial_coverage_yields_nothing_rather_than_failing() -> None:
    assert shortlist_month([_record("2020-06-01", coverage=0.2)], top_k=5) == ()


def test_top_k_must_be_positive() -> None:
    with pytest.raises(ValueError, match="top_k"):
        shortlist_month([_record("2020-06-01")], top_k=0)


def test_months_are_shortlisted_independently() -> None:
    records = [_record(f"2020-0{m}-{d:02d}", aod=0.01 * d) for m in (5, 6) for d in range(1, 9)]
    got = shortlist(records, top_k=3)
    assert set(got) == {"2020-05", "2020-06"}
    assert all(len(v) == 3 for v in got.values())


def test_saving_report_quantifies_avoided_requests() -> None:
    records = [_record(f"2020-06-{d:02d}") for d in range(1, 31)]
    got = shortlist(records, top_k=5)
    saving = request_saving(records, got)
    assert saving["acquisitions_scouted"] == 30
    assert saving["cloud_score_requests"] == 5
    assert saving["saving_fraction"] == pytest.approx(25 / 30)


def test_fully_clouded_acquisitions_never_reach_the_shortlist() -> None:
    # The clarity fallback exists so a cloudy AOI still gets a composite, but an
    # acquisition with no clear pixel cannot win one in a per-pixel mosaic, so
    # scoring it is a wasted Cloud Score+ request however scarce candidates are.
    records = [
        _record("2020-11-06", clear=0.0),
        _record("2020-11-11", clear=0.0),
        _record("2020-11-16", clear=0.03),
    ]
    chosen = shortlist_month(records, top_k=5)
    assert [r.day for r in chosen] == ["2020-11-16"]


def test_a_month_with_nothing_usable_shortlists_nothing() -> None:
    assert shortlist_month([_record("2020-11-06", clear=0.0)], top_k=5) == ()
