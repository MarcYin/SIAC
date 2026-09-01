#!/usr/bin/env python3
"""Rank library acquisitions from free signals before spending Earth Engine quota.

Cloud Score+ is the constrained resource: every acquisition costs a request, and
compositing a whole +/-45 day window fetches roughly 35 per month-bucket when
only a handful can plausibly win a pixel. L2A SCL COGs and MAIAC day AOD are
public HTTP reads under no quota at all, so they can scout the whole window for
free and Cloud Score+ then only has to score the shortlist.

That inverts the cost structure: the expensive band is fetched for candidates
already known to be clear, fully covering the AOI, and lightly loaded with
aerosol, rather than for every acquisition including the obviously cloudy ones.
It is the bestpixel selection idea used purely as a prefilter -- the actual
per-pixel weighting still comes from Cloud Score+.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

#: A partially covered AOI cannot contribute a usable composite pixel over the
#: missing part, so swath-edge acquisitions are rejected outright rather than
#: ranked down.
MIN_COVERAGE = 0.999

#: Clear-fraction floor for the preferred pool. Deliberately high: with 8
#: library years there are enough acquisitions that near-perfect ones can be
#: insisted upon, which is the quality argument for widening the year range.
PREFERRED_CLEAR = 0.98

#: Shortlist size per month-bucket.
TOP_K = 5

#: An acquisition with no clear pixel over the AOI cannot win one in the
#: per-pixel mosaic, so scoring it buys nothing at any candidate scarcity. This
#: is a floor on the clarity *fallback*, not a second gate: partial clarity is
#: genuinely useful, because the mosaic composes per pixel rather than per
#: image. Measured at Aosta in November, where nothing reaches the preferred
#: gate and the fallback was shortlisting 0.000-clear acquisitions.
MIN_USABLE_CLEAR = 0.0

#: Substituted when MAIAC has no retrieval, matching ``maiac_gap_policy``.
LOCKED_AOD_GAP_VALUE = 0.1


@dataclass(frozen=True)
class ScoutRecord:
    """One acquisition summarised from free signals only."""

    image_id: str
    day: str
    clear_fraction: float
    coverage_fraction: float
    aod: float | None = None

    @property
    def month_bucket(self) -> str:
        return str(self.day)[:7]

    @property
    def effective_aod(self) -> float:
        return LOCKED_AOD_GAP_VALUE if self.aod is None else float(self.aod)


def rank_key(record: ScoutRecord) -> tuple[float, float, str]:
    """Lowest aerosol first, then clearest, then a stable id.

    Aerosol leads because among acquisitions that are already essentially cloud
    free, the residual quality difference is atmospheric: the teacher corrects
    at the retrieved AOD, so a lightly loaded day needs less correction and
    carries less of that correction's error into the surface label.
    """

    return (record.effective_aod, -float(record.clear_fraction), record.image_id)


def shortlist_month(
    records: Sequence[ScoutRecord],
    *,
    top_k: int = TOP_K,
    min_coverage: float = MIN_COVERAGE,
    preferred_clear: float = PREFERRED_CLEAR,
) -> tuple[ScoutRecord, ...]:
    """Shortlist one month-bucket, relaxing clarity only if too few qualify.

    The coverage gate is absolute. The clear-fraction gate is not: a
    persistently cloudy AOI would otherwise contribute no composite at all,
    which silently removes exactly the regions that are hardest to sample. When
    the preferred pool is short it is topped up with the next-clearest
    acquisitions rather than failing.

    Relaxing stops at ``MIN_USABLE_CLEAR``. A month can therefore return fewer
    than ``top_k``, which is the honest answer when fewer than ``top_k``
    acquisitions can contribute a pixel.
    """

    if top_k <= 0:
        raise ValueError("top_k must be positive")
    covered = [
        r
        for r in records
        if float(r.coverage_fraction) >= float(min_coverage)
        and float(r.clear_fraction) > MIN_USABLE_CLEAR
    ]
    if not covered:
        return ()
    preferred = sorted(
        (r for r in covered if float(r.clear_fraction) >= float(preferred_clear)), key=rank_key
    )
    if len(preferred) >= top_k:
        return tuple(preferred[:top_k])
    chosen = list(preferred)
    fallback = sorted(
        (r for r in covered if r not in set(preferred)),
        key=lambda r: (-float(r.clear_fraction), r.effective_aod, r.image_id),
    )
    chosen.extend(fallback[: top_k - len(chosen)])
    return tuple(chosen)


def shortlist(
    records: Sequence[ScoutRecord],
    *,
    top_k: int = TOP_K,
    min_coverage: float = MIN_COVERAGE,
    preferred_clear: float = PREFERRED_CLEAR,
) -> dict[str, tuple[ScoutRecord, ...]]:
    """Shortlist every month-bucket independently."""

    by_month: dict[str, list[ScoutRecord]] = defaultdict(list)
    for record in records:
        by_month[record.month_bucket].append(record)
    return {
        month: shortlist_month(
            entries, top_k=top_k, min_coverage=min_coverage, preferred_clear=preferred_clear
        )
        for month, entries in sorted(by_month.items())
    }


def request_saving(records: Sequence[ScoutRecord], shortlisted: dict[str, tuple[ScoutRecord, ...]]):
    """Cloud Score+ requests avoided by scouting first."""

    scouted = len(records)
    selected = sum(len(v) for v in shortlisted.values())
    return {
        "acquisitions_scouted": scouted,
        "cloud_score_requests": selected,
        "requests_avoided": scouted - selected,
        "saving_fraction": (scouted - selected) / scouted if scouted else 0.0,
    }
