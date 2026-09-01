"""Tests for the local Cloud Score+ index builder."""

from __future__ import annotations

import json

import numpy as np
import pytest
from tools.aeronet_validation.build_cloudscore_index_local import (
    build_index,
    image_day,
    write_index,
)


def _planes(spec: dict[str, float], shape=(6, 6)) -> dict[str, np.ndarray]:
    return {image: np.full(shape, value, dtype=np.float64) for image, value in spec.items()}


def _flat_aod(value: float | None = 0.2):
    return lambda days: dict.fromkeys(days, value)


def test_image_day_parses_the_earth_engine_id() -> None:
    assert image_day("20181120T074211_20181120T074211_T37PDL") == "2018-11-20"


def test_build_index_selects_candidates_and_composes_months() -> None:
    planes = _planes(
        {
            "20190610T000000_x_T00AAA": 0.95,
            "20190620T000000_x_T00AAA": 0.30,
            "20200610T000000_x_T00AAA": 0.90,
        }
    )
    payload = build_index(planes, _flat_aod())
    # The 0.30 day is below the June mean and must be dropped.
    assert payload["counts"]["days_available"] == 3
    assert payload["counts"]["candidate_days"] == 2
    assert set(payload["months"]) == {"2019-06", "2020-06"}
    assert payload["winners"].dtype == np.int16


def test_day_coverage_uses_the_best_contributing_image() -> None:
    # Two tiles on one day: a partial tile must not drag the day down.
    planes = _planes(
        {
            "20190610T000000_x_T00AAA": 0.95,
            "20190610T000000_x_T00BBB": 0.05,
            "20190620T000000_x_T00AAA": 0.50,
        }
    )
    payload = build_index(planes, _flat_aod())
    days = {row["day"] for row in payload["day_scalars"]}
    assert "2019-06-10" in days


def test_missing_aod_falls_back_to_the_locked_constant() -> None:
    planes = _planes({"20190610T000000_x_T00AAA": 0.9, "20200610T000000_x_T00AAA": 0.9})
    payload = build_index(planes, _flat_aod(None))
    assert {row["aod_source"] for row in payload["day_scalars"]} == {"locked_constant_0p1"}
    assert all(row["aod"] == pytest.approx(0.1) for row in payload["day_scalars"])


def test_empty_input_is_rejected() -> None:
    with pytest.raises(ValueError, match="no clear-score planes"):
        build_index({}, _flat_aod())


def test_all_zero_coverage_is_rejected_rather_than_silently_empty() -> None:
    planes = _planes({"20190610T000000_x_T00AAA": 0.1, "20190620T000000_x_T00AAA": 0.1})
    with pytest.raises(ValueError, match="selected no days"):
        build_index(planes, _flat_aod())


def test_written_index_matches_the_consumed_schema(tmp_path) -> None:
    planes = _planes({"20190610T000000_x_T00AAA": 0.95, "20200610T000000_x_T00AAA": 0.9})
    payload = build_index(planes, _flat_aod())
    out = tmp_path / "scene.npz"
    write_index(payload, out, "SITE__T00AAA_20200610T000000")
    archive = np.load(out, allow_pickle=False)
    assert set(archive.files) == {
        "months",
        "winners",
        "image_table",
        "day_scalars",
        "index_policy_json",
    }
    assert archive["winners"].shape[0] == len(archive["months"])
    assert json.loads(str(archive["image_table"][0]))["idx"].startswith("2019")
    policy = json.loads(str(archive["index_policy_json"]))
    # The archive must not claim Earth Engine composed its winners.
    assert policy["winner_source"] == "local_mosaic_from_edown_cloud_score_plus"
    receipt = json.loads(out.with_suffix(".json").read_text())
    assert receipt["status"] == "ok"
    assert receipt["candidate_days"] == 2
