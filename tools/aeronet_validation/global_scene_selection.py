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

#: Scene cloud cover ceiling. The committed cohort used 20, which is a poor
#: match for inference: the model is asked to correct whatever scene it is given,
#: and a corpus of near-clear scenes under-represents broken cloud, cloud shadow
#: and adjacency effects entirely.
MAX_SCENE_CLOUD_COVER = 90.0

#: Cloud-cover targets and their shares of the catalogue. Raising the ceiling
#: alone would change nothing, because selection *minimises* cloud cover and so
#: keeps returning the clearest scene available; the corpus only spans the range
#: if scenes are drawn towards deliberately different cloud levels.
#:
#: Weighted towards the clear end on purpose. The teacher label comes from the
#: M5 solver on this scene, and the solver constrains AOD from clear pixels, so
#: a heavily clouded scene yields both less label area and a weaker retrieval.
#: The spread buys realism; the weighting keeps label quality.
CLOUD_TARGETS: tuple[tuple[float, float], ...] = (
    (5.0, 0.40),
    (20.0, 0.25),
    (45.0, 0.20),
    (75.0, 0.15),
)

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


def cloud_stratum(
    cloud_cover: float, *, targets: Sequence[tuple[float, float]] = CLOUD_TARGETS
) -> float:
    """The target level a realised cloud cover counts towards.

    Strata are the midpoints between adjacent targets, so every scene lands in
    exactly one and the realised distribution is comparable with the quota.
    """

    levels = [float(level) for level, _ in targets]
    return min(levels, key=lambda level: (abs(float(cloud_cover) - level), level))


def cloud_quota(
    total: int, *, targets: Sequence[tuple[float, float]] = CLOUD_TARGETS
) -> dict[float, int]:
    """How many scenes each cloud stratum should receive."""

    if total <= 0:
        raise ValueError("total must be positive")
    quota = {float(level): int(round(float(share) * total)) for level, share in targets}
    # Rounding must not lose or invent scenes; settle the remainder on the
    # largest stratum so the quota sums to the catalogue exactly.
    largest = max(quota, key=lambda level: (quota[level], -level))
    quota[largest] += total - sum(quota.values())
    return quota


def next_cloud_target(
    realised: Mapping[float, int],
    quota: Mapping[float, int],
) -> float:
    """The stratum currently furthest below its quota.

    Targets are chosen against what has actually been *selected* so far rather
    than assigned up front. A pre-assigned target is only a wish: an AOI that
    has no cloudy acquisition to offer silently returns a clear one, and under a
    fixed assignment that shortfall is never made up, so the cloudy strata
    finish under-filled. Re-deriving the target each time redirects the deficit
    onto AOIs that can still satisfy it.
    """

    return max(quota, key=lambda level: (quota[level] - realised.get(level, 0), -level))


def choose_scene(
    options: Sequence[SceneOption],
    *,
    target_month: int,
    target_cloud: float | None = None,
    max_cloud: float = MAX_SCENE_CLOUD_COVER,
    eligible_years: Sequence[int] = ELIGIBLE_YEARS,
) -> SceneOption | None:
    """Pick the acquisition nearest the target month and cloud level.

    Falling back to a neighbouring month keeps an AOI in the catalogue when its
    target month is unusable -- persistently cloudy tropics, or polar winter
    with no acquisitions at all -- rather than silently thinning the catalogue
    in exactly the regions hardest to sample.

    ``target_cloud`` is a level to approach, not a band to enforce. An AOI that
    is simply never cloudy has no 75%-cloud acquisition to offer, and refusing
    it would drop the driest regions from the corpus; it contributes its
    closest available scene instead. Left at ``None`` the choice reverts to the
    clearest scene.
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

    def cloud_distance(option: SceneOption) -> float:
        if target_cloud is None:
            return float(option.cloud_cover)
        return abs(float(option.cloud_cover) - float(target_cloud))

    # Nearest month first, then nearest the target cloud level, then a stable id.
    return min(usable, key=lambda o: (month_distance(o), cloud_distance(o), o.product_id))


def matchup_id(sample_id: str, scene: SceneOption) -> str:
    """Identifier in the existing ``SITE__TILE_STAMP`` form used by every stage."""

    stamp = scene.datetime.replace("-", "").replace(":", "").replace("Z", "")
    stamp = stamp.replace(" ", "T").split(".")[0]
    return f"{sample_id}__{scene.mgrs_tile}_{stamp}"


def cloud_balance(
    selected: Mapping[str, SceneOption], *, targets: Sequence[tuple[float, float]] = CLOUD_TARGETS
) -> dict[str, float]:
    """Realised share per cloud stratum, for auditing the draw."""

    if not selected:
        return {}
    counts = dict.fromkeys((float(level) for level, _ in targets), 0)
    for scene in selected.values():
        counts[cloud_stratum(scene.cloud_cover, targets=targets)] += 1
    total = len(selected)
    return {f"{level:g}%_target": counts[level] / total for level in counts}


def season_balance(selected: Mapping[str, SceneOption]) -> dict[int, float]:
    """Realised share per calendar month, for auditing the draw."""

    if not selected:
        return {}
    counts: dict[int, int] = {}
    for scene in selected.values():
        counts[scene.month] = counts.get(scene.month, 0) + 1
    total = len(selected)
    return {month: counts.get(month, 0) / total for month in range(1, 13)}
