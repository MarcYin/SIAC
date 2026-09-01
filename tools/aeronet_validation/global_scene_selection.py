#!/usr/bin/env python3
"""Choose one current-date Sentinel-2 acquisition per catalogue AOI.

The AERONET corpus selected scenes by matchup: a photometer reading within
+/-30 minutes. A global catalogue has no such anchor, so selection needs its own
policy. Two properties drive it.

*Season must be balanced deliberately.* The seasonal library window follows the
scene's calendar month, so an unbalanced season draw silently unbalances the
library too. Worse, cloud filtering favours dry seasons, which correlates season
with surface brightness -- the exact axis the catalogue stratifies on -- and
would quietly undo that stratification.

*One scene per AOI.* A second scene at the same AOI in a different month needs a
different library window, so it saves no fetches while producing a strongly
correlated sample. Reuse only pays for same-month repeats, which are the least
informative samples available.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

#: Scene cloud cover ceiling, matching the committed cloud<20 cohort.
MAX_SCENE_CLOUD_COVER = 20.0

#: Eligible acquisition years. The library spans 2018-2023, and a current-date
#: scene may sit inside or after it: the library is a seasonal climatology, not
#: a strictly prior history.
ELIGIBLE_YEARS = tuple(range(2018, 2026))


@dataclass(frozen=True)
class SceneOption:
    """One candidate acquisition returned by a STAC search."""

    product_id: str
    datetime: str
    mgrs_tile: str
    cloud_cover: float

    @property
    def month(self) -> int:
        return int(self.datetime[5:7])

    @property
    def year(self) -> int:
        return int(self.datetime[:4])


def assign_target_months(sample_ids: Sequence[str], *, seed: int = 20260901) -> dict[str, int]:
    """Give each AOI a target calendar month, evenly spread and shuffled.

    Round-robin rather than random so the draw is exactly balanced at any
    catalogue size, then shuffled so month does not correlate with the
    catalogue's land-cover ordering.
    """

    if not sample_ids:
        raise ValueError("no sample ids supplied")
    months = [1 + index % 12 for index in range(len(sample_ids))]
    rng = np.random.default_rng(seed)
    rng.shuffle(months)
    return dict(zip(sorted(sample_ids), months, strict=True))


def choose_scene(
    options: Sequence[SceneOption],
    *,
    target_month: int,
    max_cloud: float = MAX_SCENE_CLOUD_COVER,
    eligible_years: Sequence[int] = ELIGIBLE_YEARS,
) -> SceneOption | None:
    """Pick the clearest acquisition in the target month, else the nearest month.

    Falling back to a neighbouring month keeps an AOI in the catalogue when its
    target month is unusable -- persistently cloudy tropics, or polar winter
    with no acquisitions at all -- rather than silently thinning the catalogue
    in exactly the regions hardest to sample.
    """

    if not 1 <= int(target_month) <= 12:
        raise ValueError(f"target_month must be 1-12, got {target_month}")
    usable = [
        option
        for option in options
        if option.cloud_cover <= float(max_cloud) and option.year in set(eligible_years)
    ]
    if not usable:
        return None

    def month_distance(option: SceneOption) -> int:
        gap = abs(option.month - int(target_month))
        return min(gap, 12 - gap)

    # Nearest month first, then clearest, then a stable id for determinism.
    return min(usable, key=lambda o: (month_distance(o), o.cloud_cover, o.product_id))


def matchup_id(sample_id: str, scene: SceneOption) -> str:
    """Identifier in the existing ``SITE__TILE_STAMP`` form used by every stage."""

    stamp = scene.datetime.replace("-", "").replace(":", "").replace("Z", "")
    stamp = stamp.replace(" ", "T").split(".")[0]
    return f"{sample_id}__{scene.mgrs_tile}_{stamp}"


def season_balance(selected: Mapping[str, SceneOption]) -> dict[int, float]:
    """Realised share per calendar month, for auditing the draw."""

    if not selected:
        return {}
    counts: dict[int, int] = {}
    for scene in selected.values():
        counts[scene.month] = counts.get(scene.month, 0) + 1
    total = len(selected)
    return {month: counts.get(month, 0) / total for month in range(1, 13)}
