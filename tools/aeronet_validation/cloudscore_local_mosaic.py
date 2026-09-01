#!/usr/bin/env python3
"""Compose Cloud Score+ winner indices locally, without Earth Engine compute.

Two index builders exist already and neither has the properties needed for a
global corpus:

* ``build_l1c_cloudscore_index20.py`` keeps the Cloud Score+ ordering but runs
  ``qualityMosaic`` server-side, so the weighting cannot be changed without
  re-running Earth Engine over the whole corpus;
* ``build_l2a_scl_index.py`` is Earth-Engine-free but replaces Cloud Score+
  with a binary SCL mask, which cannot reproduce a continuous ordering -- every
  clear pixel ties -- and it depends on L2A coverage.

This module takes the middle path: Earth Engine (via ``edown``) supplies only
per-date ``cs_cdf`` rasters, and the quality mosaic is composed here. The
weighting then becomes an offline, editable choice and no L2A is required.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

#: Cloud Score+ clear threshold and span used by the committed weight formula
#: ``(cs_cdf - CLEAR_THRESHOLD) / CLEAR_SPAN``.
CLEAR_THRESHOLD = 0.6
CLEAR_SPAN = 0.4

#: Value written where no candidate acquisition is usable.
NO_WINNER = -1


@dataclass(frozen=True)
class Candidate:
    """One candidate acquisition contributing to a monthly mosaic."""

    day: str
    score: np.ndarray
    day_weight: float
    coverage_ratio: float


def mosaic_weight(
    score: np.ndarray,
    *,
    day_weight: float,
    coverage_ratio: float,
    clear_threshold: float = CLEAR_THRESHOLD,
    clear_span: float = CLEAR_SPAN,
) -> np.ndarray:
    """Per-pixel quality score.

    Band-agnostic: ``score`` is whichever Cloud Score+ band the index is built
    on. ``cs`` is the direct clearness estimate and is the more aggressive
    masker; ``cs_cdf`` is its CDF-normalised form. The committed index used
    ``cs_cdf``; this pipeline uses ``cs``.

    ``(score - threshold) / span + day_weight + coverage_ratio``. Only the
    score term is per-pixel; the other two are scene-day scalars, which is
    why they can be recomputed locally without touching Earth Engine.
    """

    if clear_span <= 0.0:
        raise ValueError("clear_span must be positive")
    scores = (np.asarray(score, dtype=np.float64) - float(clear_threshold)) / float(clear_span)
    return scores + float(day_weight) + float(coverage_ratio)


def compose_winners(
    candidates: Sequence[Candidate],
    *,
    clear_threshold: float = CLEAR_THRESHOLD,
    clear_span: float = CLEAR_SPAN,
    tie_breaking: str = "last",
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Argmax the per-pixel quality score across candidate acquisitions.

    ``tie_breaking`` defaults to ``"last"`` because Earth Engine's
    ``qualityMosaic`` uses mosaic semantics: among pixels sharing the maximum
    quality, the last image in collection order wins. NumPy's ``argmax``
    returns the *first* maximum, so reproducing the server-side index requires
    the reversed convention. This is not academic -- ``cs_cdf`` saturates at
    1.0 over confidently clear pixels, so ties are common rather than rare.

    Returns the winner plane (``NO_WINNER`` where nothing is usable) and the
    ordered day labels the indices refer to.
    """

    if tie_breaking not in {"first", "last"}:
        raise ValueError(f"unsupported tie_breaking {tie_breaking!r}")
    if not candidates:
        raise ValueError("no candidate acquisitions supplied")
    shapes = {candidate.score.shape for candidate in candidates}
    if len(shapes) != 1:
        raise ValueError(f"candidate rasters must share one grid; got {sorted(shapes)}")

    stack = np.stack(
        [
            mosaic_weight(
                candidate.score,
                day_weight=candidate.day_weight,
                coverage_ratio=candidate.coverage_ratio,
                clear_threshold=clear_threshold,
                clear_span=clear_span,
            )
            for candidate in candidates
        ]
    )
    usable = np.isfinite(stack)
    if tie_breaking == "last":
        # Reverse, take the first maximum, then map the index back.
        reversed_stack = stack[::-1]
        filled = np.where(np.isfinite(reversed_stack), reversed_stack, -np.inf)
        winners = (len(candidates) - 1) - np.argmax(filled, axis=0)
    else:
        filled = np.where(usable, stack, -np.inf)
        winners = np.argmax(filled, axis=0)
    winners = winners.astype(np.int32)
    winners[~usable.any(axis=0)] = NO_WINNER
    return winners, tuple(candidate.day for candidate in candidates)


def winner_day_plane(winners: np.ndarray, days: Sequence[str]) -> np.ndarray:
    """Map a winner index plane to integer day codes (``YYYYMMDD``), -1 if none."""

    codes = np.asarray([int(str(day).replace("-", "")) for day in days], dtype=np.int64)
    plane = np.full(winners.shape, -1, dtype=np.int64)
    valid = winners >= 0
    plane[valid] = codes[winners[valid]]
    return plane


def agreement(reference: np.ndarray, replica: np.ndarray) -> dict[str, float]:
    """Compare a local winner plane against a server-side reference plane."""

    reference = np.asarray(reference)
    replica = np.asarray(replica)
    if reference.shape != replica.shape:
        raise ValueError("winner planes must share one grid")
    both = (reference >= 0) & (replica >= 0)
    either = (reference >= 0) | (replica >= 0)
    return {
        "identical_fraction": float(np.mean(reference == replica)),
        "identical_fraction_where_both_valid": (
            float(np.mean(reference[both] == replica[both])) if both.any() else float("nan")
        ),
        "coverage_reference": float(np.mean(reference >= 0)),
        "coverage_replica": float(np.mean(replica >= 0)),
        "coverage_disagreement": float(np.mean(either & ~both)),
    }


def erode_valid(valid: np.ndarray, *, radius_px: int) -> np.ndarray:
    """Shrink a validity mask by ``radius_px`` using a square structuring element.

    The committed index applies a 500 m erosion to the Cloud Score+ clear mask
    before compositing (``cloud_score_erosion_radius_m``), which removes
    cloud-adjacent pixels whose score is high but whose neighbours are not.
    Reproducing the server-side winners requires reproducing this step.
    """

    mask = np.asarray(valid, dtype=bool)
    if radius_px <= 0:
        return mask.copy()
    eroded = mask.copy()
    for shift in range(1, int(radius_px) + 1):
        for axis in (0, 1):
            eroded &= np.roll(mask, shift, axis=axis)
            eroded &= np.roll(mask, -shift, axis=axis)
        # Rolling wraps at the border; treat out-of-array neighbours as invalid.
    eroded[:radius_px, :] = False
    eroded[-radius_px:, :] = False
    eroded[:, :radius_px] = False
    eroded[:, -radius_px:] = False
    return eroded


def compose_monthly_winners(
    candidates_by_month: Mapping[str, Sequence[Candidate]],
    **kwargs: object,
) -> tuple[tuple[str, ...], np.ndarray, dict[str, tuple[str, ...]]]:
    """Compose one winner plane per calendar month.

    The committed index runs ``qualityMosaic`` per calendar month of candidate
    days, and its ``winners`` array is indexed within that month's candidate
    list rather than into the global image table. Returns the ordered months,
    a stacked ``(month, y, x)`` winner array, and the per-month day ordering
    the indices refer to.
    """

    if not candidates_by_month:
        raise ValueError("no months supplied")
    months = tuple(sorted(candidates_by_month))
    planes = []
    ordering: dict[str, tuple[str, ...]] = {}
    for month in months:
        winners, days = compose_winners(list(candidates_by_month[month]), **kwargs)  # type: ignore[arg-type]
        planes.append(winners)
        ordering[month] = days
    shapes = {plane.shape for plane in planes}
    if len(shapes) != 1:
        raise ValueError(f"monthly winner planes must share one grid; got {sorted(shapes)}")
    return months, np.stack(planes).astype(np.int16), ordering
