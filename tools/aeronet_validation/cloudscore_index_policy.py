#!/usr/bin/env python3
"""Candidate selection and daily weighting for a locally composed Cloud Score+ index.

These are the parts of the committed ``siac_l1c_cloudscore_winner_index_v1``
policy that decide *which* acquisitions compete and how strongly, as opposed to
``cloudscore_local_mosaic``, which decides which one wins per pixel. Everything
here is pure NumPy so the policy can be inspected and retuned without an Earth
Engine round-trip -- the whole point of moving the mosaic off the server.
"""

from __future__ import annotations

import datetime as dt
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
from tools.aeronet_validation.acix3_surface_prior import _aod_weight, _coverage_weight

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

#: Committed policy constants, mirrored from ``index_policy_json``.
#: Cloud Score+ band. ``cs`` is the direct per-pixel clearness estimate and is
#: the more aggressive masker; ``cs_cdf`` is its CDF-normalised form. The
#: committed index used ``cs_cdf`` throughout. This pipeline uses ``cs``.
CLOUD_SCORE_BAND = "cs"

#: Clear threshold, and the offset in the ``(score - threshold) / span``
#: normalisation. 0.6 is retained for ``cs`` by decision, not inherited from the
#: committed ``cs_cdf`` index: it is Google's conventional clear cut for both
#: bands. Note ``cs`` is distributed differently, so the effective spread of the
#: ordering term is not identical to the committed one even at the same cut.
CLEAR_THRESHOLD = 0.6
SEASONAL_HALF_WIDTH_DAYS = 45
#: S2A+S2B give full 5-day revisit from 2018, so the earliest year is as dense
#: as the latest and every added year is genuine diversity rather than thin
#: coverage. Widening the range is affordable because Cloud Score+ is now
#: fetched only for a scouted shortlist, not for every acquisition in a window.
LIBRARY_YEARS = tuple(range(2018, 2026))
LOCKED_AOD_GAP_VALUE = 0.1

#: MCD19 AOD is stored as scaled integers and the committed weighting is a
#: ``locked_raw_sigmoid`` over those raw units. Feeding physical reflectance-
#: scale AOD into ``_aod_weight`` collapses its range to 0.495-0.501 instead of
#: 0.00-0.96, which would make the aerosol term nearly inert. Verified against
#: the committed ``day_scalars`` of a reference index.
MAIAC_RAW_SCALE = 1000.0


def seasonal_windows(
    reference_day: str,
    *,
    library_years: Sequence[int] = LIBRARY_YEARS,
    half_width_days: int = SEASONAL_HALF_WIDTH_DAYS,
) -> tuple[tuple[str, str], ...]:
    """Calendar-month midpoint +/- ``half_width_days`` in each library year.

    Mirrors ``seasonal_window`` in the committed policy. Windows are returned
    as inclusive ``(start, end)`` ISO dates suitable for one edown call each.
    """

    anchor = dt.date.fromisoformat(str(reference_day))
    windows = []
    for year in library_years:
        try:
            midpoint = dt.date(year, anchor.month, 15)
        except ValueError:  # pragma: no cover - month 15th always exists
            continue
        delta = dt.timedelta(days=int(half_width_days))
        windows.append(((midpoint - delta).isoformat(), (midpoint + delta).isoformat()))
    return tuple(windows)


def clean_coverage(plane: np.ndarray, *, clear_threshold: float = CLEAR_THRESHOLD) -> float:
    """Fraction of finite pixels at or above the Cloud Score+ clear threshold.

    Note the candidate rule that consumes this is *relative* -- coverage against
    the calendar-month mean -- so a more aggressive band lowers every day
    proportionally and selects the same days. Measured: cs mean coverage 0.751
    versus cs_cdf 0.798, identical day set. The band choice bites on the mosaic
    ordering, which is an absolute argmax, not here.

    This is the statistic the candidate rule thresholds. Note it is *not* used
    to mask the mosaic input: the score stays continuous there.
    """

    values = np.asarray(plane, dtype=np.float64)
    finite = np.isfinite(values)
    if not finite.any():
        return 0.0
    return float(np.mean(values[finite] >= float(clear_threshold)))


def select_candidate_days(coverage_by_day: Mapping[str, float]) -> tuple[str, ...]:
    """Apply ``daily_clean_coverage_ge_calendar_month_mean_and_gt_zero``.

    The comparison is against the mean over that *calendar month* (across
    library years), not the whole archive, so a persistently cloudy season is
    judged against its own seasonal norm rather than against clear months.
    """

    by_month: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for day, coverage in coverage_by_day.items():
        by_month[dt.date.fromisoformat(str(day)).strftime("%m")].append((str(day), float(coverage)))
    selected: list[str] = []
    for entries in by_month.values():
        values = np.asarray([coverage for _, coverage in entries], dtype=np.float64)
        threshold = float(np.mean(values))
        selected.extend(
            day for day, coverage in entries if coverage >= threshold and coverage > 0.0
        )
    return tuple(sorted(selected))


def daily_weights(
    aod_by_day: Mapping[str, float | None],
    coverage_by_day: Mapping[str, float],
    *,
    locked_gap_value: float = LOCKED_AOD_GAP_VALUE,
) -> dict[str, dict[str, object]]:
    """Per-day ``aod_weight + coverage_weight``, matching the committed order.

    Days without a MAIAC retrieval take the locked constant rather than being
    dropped (``maiac_gap_policy: locked_constant_0p1``); the source is recorded
    per day so a later audit can separate measured from substituted values.
    """

    days = sorted(coverage_by_day)
    sources = [
        "maiac" if aod_by_day.get(day) is not None else "locked_constant_0p1" for day in days
    ]
    aod = np.asarray(
        [
            float(aod_by_day[day]) if aod_by_day.get(day) is not None else float(locked_gap_value)
            for day in days
        ],
        dtype=np.float64,
    )
    coverage = np.asarray([float(coverage_by_day[day]) for day in days], dtype=np.float64)
    aod_component = np.asarray(_aod_weight(aod * MAIAC_RAW_SCALE), dtype=np.float64)
    coverage_component = np.asarray(_coverage_weight(coverage), dtype=np.float64)
    total = np.where(np.isfinite(aod_component), aod_component, 0.0) + coverage_component
    return {
        day: {
            "day": day,
            "aod": float(aod[index]),
            "aod_source": sources[index],
            "weight": float(total[index]),
        }
        for index, day in enumerate(days)
    }


def index_policy(*, library_years: Sequence[int] = LIBRARY_YEARS) -> dict[str, object]:
    """Declare how this index was produced, mirroring the committed schema.

    ``winner_source`` deliberately differs from the committed
    ``earth_engine_cloud_score_plus_index_only``: the ordering is the same but
    the argmax happened locally, and an archive must not claim otherwise.
    """

    return {
        "schema": "siac_l1c_cloudscore_winner_index_v1",
        "cloud_score_collection": "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED",
        "cloud_score_band": CLOUD_SCORE_BAND,
        "cloud_score_clear_threshold": CLEAR_THRESHOLD,
        "quality_order": (
            f"normalized_{CLOUD_SCORE_BAND} + locked_daily_aod_weight + aoi_overlap_ratio"
        ),
        "candidate_rule": "daily_clean_coverage_ge_calendar_month_mean_and_gt_zero",
        "seasonal_window": "calendar month midpoint +/-45 days in each library year",
        "library_years": list(library_years),
        "day_aod_source": "maiac",
        "maiac_gap_policy": "locked_constant_0p1",
        "aod_quality_mode": "locked_raw_sigmoid",
        "winner_source": "local_mosaic_from_edown_cloud_score_plus",
        "tie_breaking": "last",
            "score_resampling": "nearest",
    }
