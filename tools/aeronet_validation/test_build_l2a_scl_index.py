from __future__ import annotations

import json
from typing import TYPE_CHECKING

import numpy as np
import pytest
from tools.aeronet_validation.build_l2a_scl_index import (
    _aod_selection_weights,
    _cams_day_aod,
    _candidate_mask,
    _clear_distance_quality,
    _constant_aod_for_unresolved_days,
    _library_years,
    _load_dump,
    _locked_aod_selection_weights,
    _ocm_mask_and_quality,
    _reference_candidates,
    _resolve_day_aod,
    _session,
    _valid_existing_index,
)


def test_library_year_modes_keep_frozen_and_production_contracts_distinct() -> None:
    import datetime as dt

    observation = dt.date(2024, 7, 1)
    assert _library_years(observation, "fixed_2018_2023") == tuple(range(2018, 2024))
    assert _library_years(observation, "previous_5") == tuple(range(2019, 2024))
    with pytest.raises(ValueError, match="Unsupported library-year mode"):
        _library_years(observation, "future")


def test_stac_session_retries_transient_parallel_forbidden_responses() -> None:
    import requests

    retry = _session(requests).get_adapter("https://").max_retries
    assert retry.total == 7
    assert 403 in retry.status_forcelist
    assert 429 in retry.status_forcelist
    assert retry.backoff_jitter > 0


if TYPE_CHECKING:
    from pathlib import Path


def test_load_dump_accepts_mixed_resolution_local_grid(tmp_path: Path) -> None:
    from affine import Affine
    from pyproj import Transformer

    archive = tmp_path / "scene.npz"
    np.savez_compressed(
        archive,
        local60_shape=np.asarray([128, 96]),
        local60_transform=np.asarray([60.0, 0.0, 500000.0, 0.0, -60.0, 6000000.0, 0, 0, 1]),
        crs=np.asarray("EPSG:32631"),
    )
    loaded = _load_dump(archive, {"Affine": Affine, "Transformer": Transformer})
    assert (loaded["height"], loaded["width"]) == (128, 96)
    assert tuple(loaded["transform"])[:6] == (60.0, 0.0, 500000.0, 0.0, -60.0, 6000000.0)
    assert loaded["crs"] == "EPSG:32631"


def test_aod_selection_weights_are_month_local_and_prefer_clean_days() -> None:
    days = ["2020-01-01", "2020-01-06", "2020-01-11", "2020-02-01", "2020-02-06"]
    weights = _aod_selection_weights(
        days,
        {
            "2020-01-01": 0.5,
            "2020-01-06": 0.1,
            "2020-01-11": 0.3,
            "2020-02-01": 0.2,
        },
    )
    assert np.allclose(weights[:3], [0.0, 1.0, 0.5])
    assert weights[3] == 1.0
    assert np.isnan(weights[4])


def test_locked_aod_weights_reproduce_global_raw_maiac_sigmoid() -> None:
    days = ["2020-01-01", "2020-01-06", "2020-02-01", "2020-02-06"]
    physical = np.asarray([0.02, 0.04, 0.08, np.nan])
    weights = _locked_aod_selection_weights(
        days,
        dict(zip(days[:3], physical[:3], strict=True)),
    )

    raw = physical * 1000.0
    aod_min = np.nanmin(raw)
    aod_max = np.nanpercentile(raw, 60.0)
    expected = 1.0 - 1.0 / (
        1.0 + np.exp(-0.2 * (np.minimum(raw, aod_max) - (aod_max - aod_min) / 2.0))
    )
    assert np.allclose(weights[:3], expected[:3])
    assert weights[0] > weights[1] > weights[2]
    assert np.isnan(weights[3])


def test_candidates_are_ranked_within_measured_day_pool() -> None:
    coverage = np.asarray([0.9, 0.2, 0.1, 0.0], dtype=np.float32)
    months = np.asarray(["2020-01", "2020-01", "2020-01", "2020-02"])
    measured = np.asarray([False, True, True, True])
    selected = _candidate_mask(coverage, months, eligible=measured)
    assert selected.tolist() == [False, True, False, False]


def test_clear_distance_quality_rewards_cloud_interior_without_masking_edges() -> None:
    clear = np.ones((7, 7), dtype=bool)
    clear[[0, -1], :] = False
    clear[:, [0, -1]] = False
    quality = _clear_distance_quality(clear, pixel_size_m=20.0, radius_m=60.0)

    assert quality[0, 0] == 0.0
    assert quality[1, 1] == pytest.approx(1.0 / 3.0)
    assert quality[3, 3] == 1.0
    assert clear[1, 1]  # near-edge support is ranked, not hard-eroded


def test_ocm_quality_rejects_thick_and_shadow_but_retains_thin() -> None:
    confidence = np.zeros((4, 2, 3), dtype=np.float32)
    classes = np.asarray([[0, 1, 2], [3, 0, 2]])
    for class_id in range(4):
        confidence[class_id][classes == class_id] = 1.0

    valid, quality = _ocm_mask_and_quality(
        confidence,
        np.ones((2, 3), dtype=bool),
        thin_quality_weight=0.25,
        erosion=None,
    )

    assert valid.tolist() == [[True, False, True], [False, True, True]]
    assert quality.tolist() == [[1.0, 0.0, 0.25], [0.0, 1.0, 0.25]]


def test_cams_daily_fallback_reads_scaled_multitime_tiff(tmp_path: Path) -> None:
    import rasterio
    from rasterio.transform import from_origin

    token = "2020_01_02"
    folder = tmp_path / token
    folder.mkdir()
    path = folder / f"{token}_aod550.tif"
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=3,
        height=3,
        count=2,
        dtype="int16",
        crs="EPSG:4326",
        transform=from_origin(-1.5, 1.5, 1.0, 1.0),
        nodata=-32767,
    ) as target:
        target.write(np.full((2, 3, 3), 10, dtype=np.int16))
        target.scales = (0.01, 0.01)
        target.offsets = (0.0, 0.0)

    resolved = _cams_day_aod(
        ["2020-01-02", "2020-01-03"],
        scene_bounds=(-0.1, -0.1, 0.1, 0.1),
        crs="EPSG:4326",
        cams_dir=tmp_path,
    )
    assert resolved == {"2020-01-02": pytest.approx(0.1)}


def test_unresolved_cams_days_receive_locked_rt_fallback() -> None:
    fallback = _constant_aod_for_unresolved_days(
        ["2020-01-01", "2020-01-02", "2020-01-03"],
        {"2020-01-01": 0.2, "2020-01-03": 0.4},
    )
    assert fallback == {"2020-01-02": 0.1}


def test_missing_reference_never_invents_aod(tmp_path: Path) -> None:
    assert _reference_candidates(tmp_path / "missing.npz") == ([], [], {})
    resolved = _resolve_day_aod(
        source="reference",
        days=["2020-01-01"],
        reference_scalars={},
        scene_bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        cache_dir=tmp_path,
    )
    assert resolved == {}


def test_existing_index_is_reused_only_when_complete_and_readable(tmp_path: Path) -> None:
    valid = tmp_path / "valid.npz"
    np.savez_compressed(
        valid,
        months=np.asarray(["2020-01"]),
        winners=np.zeros((1, 2, 3), dtype=np.int16),
        image_table=np.asarray(["{}"]),
        day_scalars=np.asarray(["{}"]),
    )
    assert _valid_existing_index(valid)

    policy = {"ocm_confidence_quality": True, "scl_mode": "standard"}
    policy_index = tmp_path / "policy.npz"
    np.savez_compressed(
        policy_index,
        months=np.asarray(["2020-01"]),
        winners=np.zeros((1, 2, 3), dtype=np.int16),
        image_table=np.asarray(["{}"]),
        day_scalars=np.asarray(["{}"]),
        index_policy_json=np.asarray(json.dumps(policy, sort_keys=True)),
    )
    assert _valid_existing_index(policy_index, expected_policy=policy)
    assert not _valid_existing_index(
        policy_index,
        expected_policy={"ocm_confidence_quality": False, "scl_mode": "standard"},
    )
    assert not _valid_existing_index(valid, expected_policy=policy)

    incomplete = tmp_path / "incomplete.npz"
    np.savez_compressed(incomplete, months=np.asarray(["2020-01"]))
    assert not _valid_existing_index(incomplete)

    truncated = tmp_path / "truncated.npz"
    truncated.write_bytes(valid.read_bytes()[:64])
    assert not _valid_existing_index(truncated)


def test_coverage_ratio_measures_real_footprint_overlap() -> None:
    """Regression: CoordinateTransformation is on osr, not ogr.

    Calling it on ``ogr`` raises AttributeError, which the function's broad
    guard turned into a 0.0 overlap for every acquisition -- making the
    ``aoi_overlap_ratio`` term of the quality ordering identically zero rather
    than down-weighting partial tiles.
    """

    from pyproj import Transformer
    from tools.aeronet_validation.build_l2a_scl_index import _coverage_ratio, _imports

    imports = _imports()
    crs = "EPSG:32631"
    bounds = (400000.0, 4988000.0, 412000.0, 5000000.0)
    west, south, east, north = bounds
    to_wgs84 = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform

    def footprint(x_share: float) -> dict[str, object]:
        limit = west + (east - west) * x_share
        ring = [
            list(to_wgs84(*point))
            for point in [(west, south), (limit, south), (limit, north), (west, north)]
        ]
        ring.append(ring[0])
        return {"type": "Polygon", "coordinates": [ring]}

    full = _coverage_ratio({"geometry": footprint(1.0)}, bounds, crs, imports)
    half = _coverage_ratio({"geometry": footprint(0.5)}, bounds, crs, imports)
    assert full == pytest.approx(1.0, abs=0.01)
    assert half == pytest.approx(0.5, abs=0.01)
