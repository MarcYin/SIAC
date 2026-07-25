from __future__ import annotations

import json

import numpy as np
import pytest
import xarray as xr
from tools.aeronet_validation.target_tcwv import (
    l2a_wvp_grid_from_template,
    normalise_s2_l2a_wvp,
    normalise_s2_l2a_wvp_field,
    product_acquisition_time,
    resolve_target_tcwv,
)


@pytest.mark.parametrize(
    ("raw", "expected", "scale"),
    [
        (2345, 2.345, "dn_div_1000"),
        (2.345, 2.345, "physical_cm"),
        (None, None, None),
        (-1, None, None),
        (20000, None, None),
    ],
)
def test_normalise_s2_l2a_wvp(raw, expected, scale):
    value, actual_scale = normalise_s2_l2a_wvp(raw)
    assert value == pytest.approx(expected) if expected is not None else value is None
    assert actual_scale == scale


def test_resolve_target_tcwv_modes(tmp_path):
    matchup_id = "site__scene"
    wvp_dir = tmp_path / "wvp"
    wvp_dir.mkdir()
    record_path = wvp_dir / f"{matchup_id}.json"
    record_path.write_text(
        json.dumps(
            {
                "status": "OK",
                "source": "COPERNICUS/S2_SR_HARMONIZED.WVP",
                "tcwv_cm": 2.4,
            }
        ),
        encoding="utf-8",
    )

    assert resolve_target_tcwv(matchup_id, "prior", root=tmp_path, wvp_dir="wvp") is None
    fixed = resolve_target_tcwv(matchup_id, "fixed", root=tmp_path, fixed_cm=1.7, wvp_dir="wvp")
    assert fixed is not None and fixed.value_cm == pytest.approx(1.7)
    raw = resolve_target_tcwv(matchup_id, "l2a", root=tmp_path, wvp_dir="wvp")
    assert raw is not None and raw.value_cm == pytest.approx(2.4)
    scaled = resolve_target_tcwv(matchup_id, "l2a_085", root=tmp_path, wvp_dir="wvp")
    assert scaled is not None and scaled.value_cm == pytest.approx(2.04)
    assert scaled.input_value_cm == pytest.approx(2.4)
    assert scaled.record_path == str(record_path)


def test_resolve_target_tcwv_rejects_failed_record(tmp_path):
    matchup_id = "site__scene"
    wvp_dir = tmp_path / "wvp"
    wvp_dir.mkdir()
    (wvp_dir / f"{matchup_id}.json").write_text(json.dumps({"status": "FAILED"}), encoding="utf-8")
    with pytest.raises(RuntimeError, match="not OK"):
        resolve_target_tcwv(matchup_id, "l2a", root=tmp_path, wvp_dir="wvp")


def test_product_acquisition_time_matches_l1c_and_l2a():
    l1c = product_acquisition_time("S2A_MSIL1C_20240610T113321_N0510_R080_T28PCB_20240610T151027")
    l2a = product_acquisition_time("S2A_MSIL2A_20240610T113321_N0510_R080_T28PCB_20240610T164555")
    assert l1c == l2a


def test_normalise_spatial_l2a_wvp_preserves_valid_values_and_records_fill() -> None:
    field, info = normalise_s2_l2a_wvp_field(
        np.array([[1200.0, 0.0], [2.3, np.nan]], dtype=np.float32),
        fallback_cm=1.7,
    )

    assert np.allclose(field, [[1.2, 1.7], [2.3, 1.7]])
    assert info["valid_pixel_count"] == 2
    assert info["fallback_pixel_count"] == 2
    assert info["fallback_cm"] == pytest.approx(1.7)
    assert info["all_pixels_fallback"] is False


def test_l2a_wvp_grid_matches_atmospheric_template() -> None:
    template = xr.DataArray(
        np.zeros((3, 4), dtype=np.float32),
        dims=("y", "x"),
        coords={"x": [100.0, 220.0, 340.0, 460.0], "y": [940.0, 820.0, 700.0]},
    )

    grid = l2a_wvp_grid_from_template(template, crs="EPSG:32632")

    assert grid == {
        "crs": "EPSG:32632",
        "res": 120.0,
        "x0": 40.0,
        "y1": 1000.0,
        "W": 4,
        "H": 3,
    }
