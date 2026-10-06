"""Unit tests for the Landsat 8/9 OLI Collection-2 Level-1 preprocessor."""

from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np
import pytest

from siac.adapters.satellite import detect_sensor, get_preprocessor
from siac.adapters.satellite.landsat import (
    LandsatPreprocessor,
    flatten_mtl,
    parse_mtl_text,
    qa_pixel_to_cloud_classes,
)
from siac.adapters.satellite.observation import preprocess_to_observation_bundle
from siac.algorithms.rt.emulator.two_nn import TwoLayerNNEmulator
from siac.runtime.validation import validate_observation_bundle
from tests.landsat_scene import (
    ADD,
    CRS,
    HEIGHT,
    MULT,
    ORIGIN_X,
    PRODUCT_ID,
    QA_CIRRUS,
    QA_CLEAR,
    QA_CLOUD,
    QA_FILL,
    QA_SHADOW,
    SZA_DEG,
    WIDTH,
    _make_scene,
    _mtl_text,
    _qa_grid,
)

if TYPE_CHECKING:
    from pathlib import Path


def _expected_reflectance(n: int, sza_deg: float = SZA_DEG) -> float:
    return ((10000 + 1000 * n) * MULT + ADD) / np.cos(np.deg2rad(sza_deg))


# ---------------------------------------------------------------------------
# MTL parsing
# ---------------------------------------------------------------------------


def test_parse_mtl_text_nests_groups_and_strips_quotes() -> None:
    tree = parse_mtl_text(_mtl_text("LANDSAT_8", "L1TP", with_angle_names=False))
    contents = tree["LANDSAT_METADATA_FILE"]["PRODUCT_CONTENTS"]
    assert contents["LANDSAT_PRODUCT_ID"] == PRODUCT_ID
    assert contents["COLLECTION_NUMBER"] == "02"
    flat = flatten_mtl(tree)
    assert flat["SPACECRAFT_ID"] == "LANDSAT_8"
    assert float(flat["REFLECTANCE_MULT_BAND_4"]) == pytest.approx(MULT)


# ---------------------------------------------------------------------------
# Registration and detection
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("spacecraft", "expected"), [("LANDSAT_8", "l8"), ("LANDSAT_9", "l9")])
def test_detect_sensor_distinguishes_landsat_8_and_9(
    tmp_path: Path, spacecraft: str, expected: str
) -> None:
    scene = _make_scene(tmp_path, spacecraft=spacecraft)
    assert detect_sensor(scene) == expected
    assert detect_sensor(scene / f"{PRODUCT_ID}_MTL.txt") == expected


def test_detect_sensor_reads_json_only_scene(tmp_path: Path) -> None:
    scene = _make_scene(tmp_path, mtl_format="json")
    assert detect_sensor(scene) == "l9"


@pytest.mark.parametrize("sensor_id", ["l8", "l9", "landsat"])
def test_registry_resolves_landsat_preprocessor(sensor_id: str) -> None:
    assert isinstance(get_preprocessor(sensor_id), LandsatPreprocessor)


def test_platform_comes_from_mtl_not_registry_key(tmp_path: Path) -> None:
    scene = _make_scene(tmp_path, spacecraft="LANDSAT_9")
    pre = get_preprocessor("l8")
    pre.get_metadata(scene)
    assert pre.sensor_config.sensor_id == "OLI"
    assert pre.sensor_config.satellite_id == "L9"


@pytest.mark.parametrize(
    ("spacecraft", "level", "message"),
    [("LANDSAT_7", "L1TP", "Unsupported Landsat spacecraft"), ("LANDSAT_8", "L2SP", "Level-1")],
)
def test_rejects_unsupported_scenes(
    tmp_path: Path, spacecraft: str, level: str, message: str
) -> None:
    scene = _make_scene(tmp_path, spacecraft=spacecraft, level=level)
    with pytest.raises(ValueError, match=message):
        LandsatPreprocessor().get_metadata(scene)


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("mtl_format", ["txt", "json"])
def test_metadata_from_mtl(tmp_path: Path, mtl_format: str) -> None:
    scene = _make_scene(tmp_path, mtl_format=mtl_format)
    metadata = LandsatPreprocessor().get_metadata(scene)
    assert metadata["observation_time"] == datetime(2024, 6, 15, 10, 47, 12, 345678)
    assert metadata["satellite"] == "L9"
    assert metadata["tile_id"] == "199024"
    assert metadata["sun_elevation"] == pytest.approx(90.0 - SZA_DEG)
    assert metadata["reflectance_mult"]["B7"] == pytest.approx(MULT)
    assert set(metadata["reflectance_add"]) == {f"B{n}" for n in range(1, 8)}


# ---------------------------------------------------------------------------
# TOA reflectance
# ---------------------------------------------------------------------------


def test_load_toa_applies_rescaling_and_per_pixel_solar_zenith(tmp_path: Path) -> None:
    scene = _make_scene(tmp_path)
    toa = LandsatPreprocessor().load_toa(scene)
    assert set(toa.data_vars) == {f"B{n}" for n in range(1, 8)}
    for n in range(1, 8):
        band = toa[f"B{n}"]
        assert band.dtype == np.float32
        assert np.all(np.isnan(band.values[:5]))
        np.testing.assert_allclose(band.values[5:], _expected_reflectance(n), rtol=1e-5)
    assert toa.attrs["sensor"] == "OLI"
    assert toa.attrs["satellite"] == "L9"


def test_load_toa_falls_back_to_scene_sun_elevation(tmp_path: Path) -> None:
    scene = _make_scene(tmp_path, with_angles=False)
    toa = LandsatPreprocessor().load_toa(scene, band_names=("B4",))
    assert list(toa.data_vars) == ["B4"]
    np.testing.assert_allclose(toa["B4"].values[5:], _expected_reflectance(4), rtol=1e-5)


def test_band_loader_serves_bands_not_preloaded(tmp_path: Path) -> None:
    scene = _make_scene(tmp_path)
    toa = LandsatPreprocessor().load_toa(scene, band_names=("B4",))
    loader = toa.attrs["_siac_toa_band_loader"]
    b6 = loader("B6")
    np.testing.assert_allclose(b6.values[5:], _expected_reflectance(6), rtol=1e-5)
    with pytest.raises(KeyError):
        loader("B4")


def test_band_files_found_by_suffix_when_mtl_names_absent(tmp_path: Path) -> None:
    scene = _make_scene(tmp_path, with_angle_names=False)
    geometry = LandsatPreprocessor().extract_geometry(scene)
    assert float(geometry.sza.mean()) == pytest.approx(np.deg2rad(SZA_DEG), rel=1e-5)


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def test_extract_geometry_block_averages_angle_bands(tmp_path: Path) -> None:
    scene = _make_scene(tmp_path)
    pre = LandsatPreprocessor(config={"angle_block_pixels": 10})
    geometry = pre.extract_geometry(scene)
    assert geometry.sza.shape == (4, 5)
    # Fill rows are excluded from the block mean, and no block is left empty.
    np.testing.assert_allclose(geometry.sza.values, np.deg2rad(SZA_DEG), rtol=1e-5)
    np.testing.assert_allclose(geometry.saa.values, np.deg2rad(150.0), rtol=1e-5)
    assert np.all(np.isfinite(geometry.vza.values))
    # -179 deg and +179 deg average circularly to 180 deg (an arithmetic mean gives 0).
    vaa_deg = np.mod(np.rad2deg(geometry.vaa.values), 360.0)
    np.testing.assert_allclose(vaa_deg, 180.0, atol=1e-3)
    # Block centres sit on the native grid.
    assert float(geometry.sza.x[0]) == pytest.approx(ORIGIN_X + 15.0 + 4.5 * 30.0)
    assert str(geometry.sza.rio.crs) == CRS


def test_extract_geometry_requires_angle_bands(tmp_path: Path) -> None:
    scene = _make_scene(tmp_path, with_angles=False)
    with pytest.raises(FileNotFoundError, match="angle band"):
        LandsatPreprocessor().extract_geometry(scene)


# ---------------------------------------------------------------------------
# Cloud classes
# ---------------------------------------------------------------------------


def test_qa_pixel_bits_map_to_standard_classes() -> None:
    qa = np.array([QA_FILL, QA_CLEAR, QA_CLOUD, QA_SHADOW, QA_CIRRUS, QA_CLOUD | QA_SHADOW])
    np.testing.assert_array_equal(qa_pixel_to_cloud_classes(qa), [0, 1, 2, 3, 2, 2])


def test_product_qa_cloud_mask(tmp_path: Path) -> None:
    scene = _make_scene(tmp_path)
    pre = LandsatPreprocessor(config={"cloud_mask": {"provider": "product_qa"}})
    mask = pre.extract_cloud_mask(scene)
    classes = pre._last_cloud_classes
    assert classes is not None
    np.testing.assert_array_equal(classes.values, qa_pixel_to_cloud_classes(_qa_grid()))
    assert mask.dtype == bool
    assert bool(mask.values[0, 0])  # fill is not usable
    assert not bool(mask.values[7, 0])  # clear
    assert bool(mask.values[12, 12])  # cloud
    assert bool(mask.values[22, 12])  # shadow


# ---------------------------------------------------------------------------
# M1 contract
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cloud_mask", [{"provider": "product_qa"}, {"mode": "none"}], ids=["product_qa", "none"]
)
def test_preprocess_produces_valid_observation_bundle(tmp_path: Path, cloud_mask: dict) -> None:
    scene = _make_scene(tmp_path)
    pre = LandsatPreprocessor(config={"cloud_mask": cloud_mask, "angle_block_pixels": 10})
    bundle = preprocess_to_observation_bundle(pre, scene)
    validate_observation_bundle(bundle)
    assert bundle.sensor_config.satellite_id == "L9"
    assert bundle.crs == CRS
    assert bundle.cloud_mask.shape == (HEIGHT, WIDTH)
    assert set(bundle.toa.data_vars) == {f"B{n}" for n in range(1, 8)}


# ---------------------------------------------------------------------------
# RT emulator platform fallback
# ---------------------------------------------------------------------------


def test_landsat9_emulator_falls_back_to_landsat8_weights(tmp_path: Path) -> None:
    for n in range(1, 8):
        for part in ("xap", "xbp", "xcp"):
            (tmp_path / f"isotropic_OLI_L8_B{n}_{part}.npz").touch()
    emulator = TwoLayerNNEmulator(emulator_dir=tmp_path, sensor_id="OLI", satellite_id="L9")
    assert emulator.available_bands == [f"B{n}" for n in range(1, 8)]
    path = emulator._get_emulator_path("B4")
    assert path is not None
    assert "_L8_B4_" in path.name
