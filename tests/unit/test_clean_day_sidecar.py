from __future__ import annotations

import json

from tools.aeronet_validation.clean_day_sidecar import (
    select_clean_scenes,
    select_clean_windows,
    write_filtered_sidecar,
)


def _scenes() -> dict[str, dict[str, float]]:
    return {
        "d": {"maiac_aod": 0.40},
        "b": {"maiac_aod": 0.12},
        "a": {"maiac_aod": 0.05},
        "c": {"maiac_aod": 0.20},
    }


def test_select_clean_scenes_prefers_absolute_threshold():
    selected, quality = select_clean_scenes(_scenes(), max_aod=0.15, min_scenes=2)
    assert list(selected) == ["a", "b"]
    assert quality["n_preferred"] == 2
    assert quality["n_fallback"] == 0
    assert quality["selected_aod_max"] == 0.12


def test_select_clean_scenes_backfills_to_minimum():
    selected, quality = select_clean_scenes(_scenes(), max_aod=0.08, min_scenes=3)
    assert list(selected) == ["a", "b", "c"]
    assert quality["n_preferred"] == 1
    assert quality["n_fallback"] == 2
    assert quality["selected_aod_max"] == 0.20


def test_write_filtered_sidecar_preserves_metadata(tmp_path):
    source = tmp_path / "source.json"
    destination = tmp_path / "filtered.json"
    source.write_text(json.dumps({"bands": ["blue"], "scenes": _scenes()}), encoding="utf-8")
    quality = write_filtered_sidecar(source, destination, max_aod=0.15, min_scenes=2)
    result = json.loads(destination.read_text(encoding="utf-8"))
    assert result["bands"] == ["blue"]
    assert list(result["scenes"]) == ["a", "b"]
    assert result["clean_day_selection"] == quality


def test_select_clean_windows_keeps_preferred_then_global_lowest():
    windows = [
        {"window": "2020-06", "selected_aod_median": 0.32},
        {"window": "2021-06", "selected_aod_median": 0.11},
        {"window": "2022-06", "selected_aod_median": 0.19},
        {"window": "2023-06", "selected_aod_median": 0.14},
        {"window": "2024-06", "selected_aod_median": 0.25},
    ]

    selected, quality = select_clean_windows(
        windows,
        max_median_aod=0.15,
        min_windows=3,
    )

    assert selected == [1, 2, 3]
    assert quality["n_preferred_windows"] == 2
    assert quality["n_fallback_windows"] == 1
    assert quality["selected_window_aod_max"] == 0.19


def test_select_clean_windows_keeps_all_clean_windows():
    windows = [
        {"selected_aod_median": 0.05},
        {"selected_aod_median": 0.08},
        {"selected_aod_median": 0.12},
        {"selected_aod_median": 0.30},
    ]

    selected, quality = select_clean_windows(
        windows,
        max_median_aod=0.15,
        min_windows=2,
    )

    assert selected == [0, 1, 2]
    assert quality["n_selected_windows"] == 3
    assert quality["n_fallback_windows"] == 0
