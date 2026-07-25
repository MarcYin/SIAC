"""Adaptive clean-day selection for L1C surface-prior sidecars."""

from __future__ import annotations

import json
import math
from pathlib import Path
from statistics import median
from typing import Any


def select_clean_scenes(
    scenes: dict[str, dict[str, Any]],
    *,
    max_aod: float,
    min_scenes: int,
) -> tuple[dict[str, dict[str, Any]], dict[str, float | int]]:
    """Prefer an absolute clean threshold, then backfill by ascending AOD."""

    if not math.isfinite(max_aod) or max_aod <= 0.0:
        raise ValueError("max_aod must be positive and finite")
    if min_scenes < 1:
        raise ValueError("min_scenes must be at least one")
    ranked = sorted(
        (
            (scene_id, scene, float(scene["maiac_aod"]))
            for scene_id, scene in scenes.items()
            if math.isfinite(float(scene["maiac_aod"]))
        ),
        key=lambda item: (item[2], item[0]),
    )
    if not ranked:
        raise ValueError("sidecar contains no finite MAIAC AOD values")
    preferred_count = sum(aod <= max_aod for _, _, aod in ranked)
    selected_count = min(len(ranked), max(preferred_count, min_scenes))
    selected_rows = ranked[:selected_count]
    selected = {scene_id: scene for scene_id, scene, _ in selected_rows}
    aods = [aod for _, _, aod in selected_rows]
    quality: dict[str, float | int] = {
        "n_total": len(ranked),
        "n_preferred": preferred_count,
        "n_selected": selected_count,
        "n_fallback": max(0, selected_count - preferred_count),
        "max_aod_threshold": float(max_aod),
        "selected_aod_min": min(aods),
        "selected_aod_median": median(aods),
        "selected_aod_max": max(aods),
    }
    return selected, quality


def select_clean_windows(
    windows: list[dict[str, Any]],
    *,
    max_median_aod: float,
    min_windows: int,
) -> tuple[list[int], dict[str, float | int]]:
    """Prefer clean monthly composites, then backfill globally by median AOD.

    Per-month scene backfilling is needed to form a spatial composite, but it
    can leave every month represented by a polluted realization. This second
    stage ranks the completed monthly composites together so those fallback
    months do not receive equal weight in the seasonal prior.
    """

    if not math.isfinite(max_median_aod) or max_median_aod <= 0.0:
        raise ValueError("max_median_aod must be positive and finite")
    if min_windows < 1:
        raise ValueError("min_windows must be at least one")
    ranked = sorted(
        (
            (index, float(window["selected_aod_median"]))
            for index, window in enumerate(windows)
            if math.isfinite(float(window["selected_aod_median"]))
        ),
        key=lambda item: (item[1], item[0]),
    )
    if not ranked:
        raise ValueError("window metadata contains no finite median AOD values")
    preferred_count = sum(aod <= max_median_aod for _, aod in ranked)
    selected_count = min(len(ranked), max(preferred_count, min_windows))
    selected_rows = ranked[:selected_count]
    selected_indices = sorted(index for index, _ in selected_rows)
    aods = [aod for _, aod in selected_rows]
    quality: dict[str, float | int] = {
        "n_total_windows": len(ranked),
        "n_preferred_windows": preferred_count,
        "n_selected_windows": selected_count,
        "n_fallback_windows": max(0, selected_count - preferred_count),
        "max_median_aod_threshold": float(max_median_aod),
        "selected_window_aod_min": min(aods),
        "selected_window_aod_median": median(aods),
        "selected_window_aod_max": max(aods),
    }
    return selected_indices, quality


def write_filtered_sidecar(
    source: str | Path,
    destination: str | Path,
    *,
    max_aod: float,
    min_scenes: int,
) -> dict[str, float | int]:
    """Write an AtmoSidecar containing only adaptively selected scenes."""

    source_path = Path(source)
    destination_path = Path(destination)
    payload = json.loads(source_path.read_text(encoding="utf-8"))
    selected, quality = select_clean_scenes(
        payload["scenes"], max_aod=max_aod, min_scenes=min_scenes
    )
    output = dict(payload)
    output["scenes"] = selected
    output["clean_day_selection"] = quality
    destination_path.write_text(json.dumps(output), encoding="utf-8")
    return quality
