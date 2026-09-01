#!/usr/bin/env python3
"""Stratified global sampling of AOI centres for the expanded training corpus.

The existing corpus is anchored on AERONET sites, which exists for validation
convenience rather than method necessity: ``build_optimal_m5_teacher.py`` never
reads AERONET. That anchoring is also what starves the corpus of the surfaces
that actually fail -- ``bare_sparse`` is 8.4% of scenes while the bright octile
carries about a third of the error, and only 189 of 523 sites contribute any
bright scene, because sun photometers are not sited in deserts.

Sampling is therefore over land globally, stratified toward the failure mode.
Selection is deterministic given a seed so a catalogue can be regenerated and
audited.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

import numpy as np

#: Share of the catalogue reserved for each land-cover class. Everything not
#: named here competes for the remainder in proportion to its availability.
#: ``bare_sparse`` is lifted from its natural ~8% because that is where the
#: residual error concentrates; it is not a claim about global land area.
DEFAULT_TARGETS: dict[str, float] = {
    "bare_sparse": 0.25,
    "shrubland": 0.10,
    "grassland": 0.10,
}

#: Latitude bands used to keep illumination geometry diverse.
LATITUDE_BANDS: tuple[tuple[str, float, float], ...] = (
    ("tropical", 0.0, 23.5),
    ("subtropical", 23.5, 45.0),
    ("temperate", 45.0, 60.0),
    ("polar", 60.0, 90.0),
)


@dataclass(frozen=True)
class Candidate:
    """One land location eligible for the catalogue."""

    sample_id: str
    longitude: float
    latitude: float
    land_cover: str
    continent: str


def latitude_band(latitude: float) -> str:
    magnitude = abs(float(latitude))
    for name, lower, upper in LATITUDE_BANDS:
        if lower <= magnitude < upper:
            return name
    return LATITUDE_BANDS[-1][0]


def _quota(
    targets: Mapping[str, float], available: Mapping[str, int], total: int
) -> dict[str, int]:
    """Allocate per-class counts, capped by what is actually available.

    A target that cannot be met from the candidate pool is filled as far as it
    goes and the shortfall returned to the open remainder, rather than silently
    producing a short catalogue.
    """

    quota: dict[str, int] = {}
    for name, share in targets.items():
        quota[name] = min(int(round(share * total)), int(available.get(name, 0)))
    remaining = total - sum(quota.values())
    others = {k: v for k, v in available.items() if k not in targets and v > 0}
    pool = sum(others.values())
    if remaining > 0 and pool > 0:
        for name, count in others.items():
            quota[name] = min(int(round(remaining * count / pool)), count)
    return quota


def sample_catalog(
    candidates: Sequence[Candidate],
    *,
    total: int,
    targets: Mapping[str, float] | None = None,
    seed: int = 20260901,
) -> tuple[Candidate, ...]:
    """Draw a stratified catalogue, balancing latitude and continent within class.

    Within each land-cover quota the draw is round-robin over
    ``(continent, latitude band)`` cells, so a class dominated by one region
    cannot monopolise its allocation.
    """

    if total <= 0:
        raise ValueError("total must be positive")
    if not candidates:
        raise ValueError("no candidates supplied")
    targets = dict(DEFAULT_TARGETS if targets is None else targets)
    if any(share < 0 for share in targets.values()):
        raise ValueError("targets must be non-negative")
    if sum(targets.values()) > 1.0 + 1e-9:
        raise ValueError("targets must not exceed 1.0 in total")

    by_class: dict[str, list[Candidate]] = {}
    for candidate in candidates:
        by_class.setdefault(candidate.land_cover, []).append(candidate)
    quota = _quota(targets, {k: len(v) for k, v in by_class.items()}, total)

    rng = np.random.default_rng(seed)
    selected: list[Candidate] = []
    for land_cover, count in sorted(quota.items()):
        if count <= 0:
            continue
        cells: dict[tuple[str, str], list[Candidate]] = {}
        for candidate in by_class[land_cover]:
            cells.setdefault((candidate.continent, latitude_band(candidate.latitude)), []).append(
                candidate
            )
        for members in cells.values():
            rng.shuffle(members)
        order = sorted(cells)
        taken = 0
        while taken < count and any(cells[key] for key in order):
            for key in order:
                if taken >= count:
                    break
                if cells[key]:
                    selected.append(cells[key].pop())
                    taken += 1
    selected.sort(key=lambda c: c.sample_id)
    return tuple(selected)


def composition(selected: Sequence[Candidate]) -> dict[str, dict[str, float]]:
    """Realised shares, for auditing a catalogue against its targets."""

    if not selected:
        return {}
    total = len(selected)
    out: dict[str, dict[str, float]] = {}
    for field, key in (
        ("land_cover", lambda c: c.land_cover),
        ("continent", lambda c: c.continent),
        ("latitude_band", lambda c: latitude_band(c.latitude)),
    ):
        counts: dict[str, int] = {}
        for candidate in selected:
            counts[key(candidate)] = counts.get(key(candidate), 0) + 1
        out[field] = {name: count / total for name, count in sorted(counts.items())}
    return out
