import json

import numpy as np
import pytest
import tools.aeronet_validation.build_harmonized_l2a_histories as harmonized_histories
import tools.aeronet_validation.build_paired_siac_histories as paired_histories
import xarray as xr
from tools.aeronet_validation.build_harmonized_l2a_histories import (
    _candidate,
    _correct_scene,
)
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import (
    CANONICAL_TO_S2_BAND,
    PHYSICAL_TARGET_RT,
    TARGET_RT,
    _cached_result_is_compatible,
    _has_only_sparse_land_errors,
    _physical_current_rt_surface,
    _selected_metadata,
    _sentinel_satellite_id,
)
from tools.aeronet_validation.build_per_acquisition_l2a_l1c_histories import (
    _correct_scene as _correct_pair_affine_scene,
)
from tools.aeronet_validation.build_scene_harmonized_l2a_histories import (
    _apply_scene_correction,
)
from tools.aeronet_validation.per_acquisition_l2a_l1c_mapping import (
    apply_band_map,
    build_scene_maps,
    fit_band_map,
    scene_key,
)
from tools.aeronet_validation.terrain_features import TerrainFields, local_solar_incidence
from tools.aeronet_validation.train_l2a_l1c_harmonizer import (
    BAND_NAMES,
    _acquisition_identity,
    _fit_model,
    _fixed_fold_splits,
    _predict_model,
    _scene_balanced_sample_weight,
    build_features,
    feature_names,
    load_pairs,
)
from tools.aeronet_validation.train_scene_l2a_harmonizer import (
    build_scene_features,
)
from tools.aeronet_validation.train_scene_l2a_harmonizer import (
    feature_names as scene_feature_names,
)

from siac.runtime import RTCoefficients


def test_selected_metadata_keeps_lowest_maiac_aod_fraction() -> None:
    rows = [
        {"l1c_id": "high", "maiac": 0.5},
        {"l1c_id": "low", "maiac": 0.1},
        {"l1c_id": "invalid", "maiac": None},
        {"l1c_id": "middle", "maiac": 0.3},
        {"l1c_id": "lower", "maiac": 0.2},
    ]

    selected = _selected_metadata(rows, 0.5)

    assert [row["l1c_id"] for row in selected] == ["low", "lower"]


def test_sparse_land_support_is_not_classified_as_a_processing_failure() -> None:
    sparse = [
        {
            "error_type": "RuntimeError",
            "reason": "only 52 paired clear-land pixels",
        }
    ]
    network = [{"error_type": "OSError", "reason": "connection reset"}]

    assert _has_only_sparse_land_errors(sparse) is True
    assert _has_only_sparse_land_errors([]) is True
    assert _has_only_sparse_land_errors(network) is False


def test_pair_cache_rejects_post_cutoff_scenes(tmp_path) -> None:
    archive = tmp_path / "pairs.npz"
    np.savez_compressed(
        archive,
        scenes_json=np.asarray(json.dumps([{"day": "2023-08-01"}, {"day": "2025-01-02"}])),
    )
    existing = {
        "status": "OK",
        "target_rt": PHYSICAL_TARGET_RT,
        "terrain_features": {"enabled": True},
    }

    assert not _cached_result_is_compatible(
        existing,
        archive,
        keep_fraction=0.6,
        max_samples_per_scene=1200,
        include_terrain=True,
        target_rt_mode="physical",
        scene_day_max="2023-12-31",
    )


def test_pair_cache_accepts_legacy_archive_that_already_satisfies_cutoff(tmp_path) -> None:
    archive = tmp_path / "pairs.npz"
    np.savez_compressed(
        archive,
        scenes_json=np.asarray(json.dumps([{"day": "2022-07-01"}, {"day": "2023-08-01"}])),
    )
    existing = {
        "status": "OK",
        "target_rt": PHYSICAL_TARGET_RT,
        "terrain_features": {"enabled": True},
    }

    assert _cached_result_is_compatible(
        existing,
        archive,
        keep_fraction=0.6,
        max_samples_per_scene=1200,
        include_terrain=True,
        target_rt_mode="physical",
        scene_day_max="2023-12-31",
    )
    assert not _cached_result_is_compatible(
        existing | {"scene_day_max": "2022-12-31"},
        archive,
        keep_fraction=0.6,
        max_samples_per_scene=1200,
        include_terrain=True,
        target_rt_mode="physical",
        scene_day_max="2023-12-31",
    )


def test_pair_cache_requires_explicit_config_for_unavailable_results(tmp_path) -> None:
    existing = {
        "status": "DATA_UNAVAILABLE",
        "keep_fraction": 0.6,
        "max_samples_per_scene": 1200,
        "scene_day_max": "2023-12-31",
    }
    kwargs = {
        "keep_fraction": 0.6,
        "max_samples_per_scene": 1200,
        "include_terrain": True,
        "target_rt_mode": "physical",
        "scene_day_max": "2023-12-31",
    }

    assert not _cached_result_is_compatible(existing, tmp_path / "missing.npz", **kwargs)
    assert _cached_result_is_compatible(
        existing | {"target_rt_mode": "physical", "include_terrain": True},
        tmp_path / "missing.npz",
        **kwargs,
    )


def test_physical_pair_target_uses_l2a_wvp_and_glo30_elevation() -> None:
    class FakeSensor:
        def get_band(self, name: str) -> str:
            return name

    class FakeRT:
        def __init__(self) -> None:
            self.geometry = None
            self.atmo = None
            self.bands = None

        def preload_scene_subset(self, geometry, atmo, bands) -> None:  # noqa: ANN001
            self.geometry = geometry
            self.atmo = atmo
            self.bands = bands

        def compute_coefficients(self, geometry, atmo, _band):  # noqa: ANN001
            return RTCoefficients(
                xap=xr.ones_like(geometry.sza),
                xbp=xr.zeros_like(geometry.sza),
                xcp=xr.zeros_like(geometry.sza),
            )

    terrain = TerrainFields(
        elevation_m=np.asarray([[100.0, 250.0], [300.0, 450.0]], dtype=np.float32),
        slope_deg=np.zeros((2, 2), dtype=np.float32),
        gradient_east=np.zeros((2, 2), dtype=np.float32),
        gradient_north=np.zeros((2, 2), dtype=np.float32),
        valid=np.ones((2, 2), dtype=bool),
    )
    l1_toa = np.linspace(0.01, 0.28, 28, dtype=np.float32).reshape(7, 2, 2)
    l2a_tcwv = np.asarray([[0.8, 1.1], [1.5, 1.9]], dtype=np.float32)
    backend = FakeRT()

    corrected = _physical_current_rt_surface(
        l1_toa=l1_toa,
        l2a_tcwv=l2a_tcwv,
        terrain=terrain,
        scene={"maiac": 0.2, "sza": 30.0, "saa": 120.0, "vza": 5.0, "vaa": 15.0},
        spacecraft="Sentinel-2A",
        fallback_tcwv=1.0,
        rt_backend=backend,
        sensor_cache={"S2A": FakeSensor()},
    )

    np.testing.assert_allclose(corrected, l1_toa)
    np.testing.assert_allclose(backend.atmo.tcwv.values, l2a_tcwv)
    np.testing.assert_allclose(backend.atmo.elevation.values, terrain.elevation_m / 1000.0)
    assert backend.bands == [CANONICAL_TO_S2_BAND[name] for name in CANONICAL_TO_S2_BAND]
    assert _sentinel_satellite_id("Sentinel-2A") == "S2A"


def test_compact_features_encode_atmospheric_differences() -> None:
    l2a = np.asarray(
        [[0.01, 0.02, 0.03, 0.04, 0.2, 0.1, 0.08]] * 2,
        dtype=np.float32,
    )
    features = build_features(
        l2a,
        band_index=1,
        feature_set="compact",
        l2a_aot=np.asarray([0.1, 0.3]),
        l2a_tcwv=np.asarray([1.0, 1.5]),
        maiac_aot=np.asarray([0.3, 0.2]),
        maiac_tcwv=np.asarray([1.4, 1.0]),
        sza_deg=np.asarray([30.0, 40.0]),
        vza_deg=np.asarray([5.0, 10.0]),
        raa_deg=np.asarray([0.0, 180.0]),
        elevation_km=np.asarray([0.1, 1.0]),
        month=np.asarray([3, 9]),
        sensor_is_s2b=np.asarray([0.0, 1.0]),
        sensor_is_s2c=np.asarray([0.0, 0.0]),
        processing_baseline=np.asarray([3.01, 5.11]),
    )

    assert features.shape == (2, len(feature_names("compact", 1)))
    np.testing.assert_allclose(features[:, 0], 0.02)
    np.testing.assert_allclose(features[:, 1], [0.2, -0.1])
    np.testing.assert_allclose(features[:, 2], [0.004, -0.002])
    np.testing.assert_allclose(features[:, 3], [0.4, -0.5])
    np.testing.assert_allclose(features[:, 7], [1.0, -1.0], atol=1.0e-12)


def test_exported_linear_model_prediction_is_deterministic() -> None:
    rng = np.random.default_rng(42)
    features = rng.normal(size=(80, 5))
    target = 0.003 + features @ np.asarray([0.2, -0.1, 0.04, 0.0, 0.08])

    model = _fit_model(features, target, alpha=10.0)
    first = _predict_model(features, model)
    second = _predict_model(features.copy(), json.loads(json.dumps(model)))

    np.testing.assert_allclose(first, second, rtol=0.0, atol=1.0e-12)
    assert np.sqrt(np.mean(np.square(first - target))) < 0.03


def test_scene_balancing_gives_each_acquisition_equal_total_weight() -> None:
    scene_code = np.asarray([0, 0, 0, 1, 2, 2], dtype=np.int16)

    weight = _scene_balanced_sample_weight(scene_code)

    np.testing.assert_allclose(
        [weight[scene_code == code].sum() for code in range(3)],
        [1.0, 1.0, 1.0],
    )


def test_acquisition_identity_groups_adjacent_tiles_from_one_overpass() -> None:
    first = {
        "scene_id": "S2A_48MXU_20200608_0_L1C",
        "l2a": {"system_index": "20200608T025551_20200608T031913_T48MXU"},
    }
    second = {
        "scene_id": "S2A_48MYU_20200608_0_L1C",
        "l2a": {"system_index": "20200608T025551_20200608T031913_T48MYU"},
    }

    assert _acquisition_identity(first, 0) == "20200608T025551"
    assert _acquisition_identity(second, 1) == _acquisition_identity(first, 0)


def test_locked_holdout_surface_model_trains_only_on_development_folds() -> None:
    raw_folds = np.repeat(np.arange(5, dtype=np.int16), 2)

    application_folds, splits, protocol = _fixed_fold_splits(
        raw_folds,
        holdout_folds={2, 3},
    )

    assert protocol["development_folds"] == [0, 1, 4]
    assert protocol["holdout_model_fold"] == 5
    np.testing.assert_array_equal(application_folds[4:8], 5)
    by_fold = {fold: (train, test) for fold, train, test in splits}
    np.testing.assert_array_equal(np.unique(raw_folds[by_fold[0][0]]), [1, 4])
    np.testing.assert_array_equal(np.unique(raw_folds[by_fold[5][0]]), [0, 1, 4])
    np.testing.assert_array_equal(np.unique(raw_folds[by_fold[5][1]]), [2, 3])


def test_cross_rt_features_retain_equal_aod_baseline_and_terrain_terms() -> None:
    l2a = np.asarray(
        [[0.02, 0.04, 0.05, 0.06, 0.2, 0.1, 0.08]] * 2,
        dtype=np.float32,
    )
    feature_set = "cross_rt_terrain"
    names = feature_names(feature_set, 1)
    features = build_features(
        l2a,
        band_index=1,
        feature_set=feature_set,
        l2a_aot=np.asarray([0.2, 0.4]),
        l2a_tcwv=np.asarray([1.0, 1.5]),
        maiac_aot=np.asarray([0.2, 0.4]),
        maiac_tcwv=np.asarray([1.0, 1.5]),
        sza_deg=np.asarray([30.0, 40.0]),
        vza_deg=np.asarray([5.0, 8.0]),
        raa_deg=np.asarray([20.0, 100.0]),
        elevation_km=np.asarray([0.1, 0.2]),
        month=np.asarray([3, 9]),
        sensor_is_s2b=np.asarray([0.0, 1.0]),
        sensor_is_s2c=np.asarray([0.0, 0.0]),
        processing_baseline=np.asarray([5.0, 5.11]),
        terrain_elevation_km=np.asarray([0.15, 0.35]),
        terrain_slope_deg=np.asarray([5.0, 20.0]),
        terrain_incidence_cos=np.asarray([0.8, 0.5]),
    )

    assert features.shape == (2, len(names))
    np.testing.assert_allclose(
        features[:, names.index("delta_aot_maiac_minus_sen2cor")],
        0.0,
    )
    np.testing.assert_allclose(features[:, names.index("mean_aot")], [0.2, 0.4])
    np.testing.assert_allclose(features[:, names.index("l2a_tcwv_cm")], [1.0, 1.5])
    assert "delta_tcwv_maiac_minus_sen2cor" not in names
    np.testing.assert_allclose(
        features[:, names.index("terrain_elevation_delta_km")],
        [0.05, 0.15],
    )
    assert np.all(np.isfinite(features))


def test_load_pairs_aligns_scene_metadata_by_sample_index(tmp_path) -> None:
    matchup_id = "Example_Site__T30ABC_20200101T000000"
    sample_count = 4
    l2a = (
        np.arange(sample_count * len(BAND_NAMES), dtype=np.float32).reshape(
            sample_count, len(BAND_NAMES)
        )
        / 100.0
    )
    scenes = [
        {
            "day": "2023-06-01",
            "maiac_aot": 0.2,
            "maiac_tcwv_cm": 1.0,
            "sza_deg": 30.0,
            "vza_deg": 5.0,
            "raa_deg": 20.0,
            "elevation_km": 0.1,
            "month": 2,
        },
        {
            "day": "2024-06-01",
            "maiac_aot": 0.4,
            "maiac_tcwv_cm": 2.0,
            "sza_deg": 40.0,
            "vza_deg": 10.0,
            "raa_deg": 60.0,
            "elevation_km": 0.5,
            "month": 8,
        },
    ]
    np.savez_compressed(
        tmp_path / f"{matchup_id}.npz",
        l2a=l2a,
        siac=l2a + 0.01,
        l2a_aot=np.full(sample_count, 0.15, dtype=np.float32),
        l2a_tcwv=np.full(sample_count, 1.2, dtype=np.float32),
        scene_index=np.asarray([1, 0, 1, 0], dtype=np.int16),
        band_names=np.asarray(BAND_NAMES),
        scenes_json=np.asarray(json.dumps(scenes)),
    )
    (tmp_path / f"{matchup_id}.json").write_text(
        json.dumps({"status": "OK", "uses_aeronet": False}), encoding="utf-8"
    )

    dataset = load_pairs(tmp_path, [matchup_id])

    np.testing.assert_allclose(dataset.maiac_aot, [0.4, 0.2, 0.4, 0.2])
    np.testing.assert_allclose(dataset.maiac_tcwv, [2.0, 1.0, 2.0, 1.0])
    np.testing.assert_array_equal(dataset.month, [8, 2, 8, 2])
    np.testing.assert_array_equal(dataset.scene_code, [1, 0, 1, 0])
    np.testing.assert_array_equal(dataset.acquisition_code, [1, 0, 1, 0])
    assert dataset.target_rt == TARGET_RT

    cutoff_dataset = load_pairs(tmp_path, [matchup_id], scene_day_max="2023-12-31")
    np.testing.assert_allclose(cutoff_dataset.maiac_aot, [0.2, 0.2])
    np.testing.assert_array_equal(cutoff_dataset.month, [2, 2])

    clean_dataset = load_pairs(tmp_path, [matchup_id], maiac_aot_max=0.3)
    np.testing.assert_allclose(clean_dataset.maiac_aot, [0.2, 0.2])
    np.testing.assert_array_equal(clean_dataset.month, [2, 2])

    unavailable_id = "No_Pairs__T00XXX_20240101T000000"
    (tmp_path / f"{unavailable_id}.json").write_text(
        json.dumps({"status": "DATA_UNAVAILABLE"}), encoding="utf-8"
    )
    partial_dataset = load_pairs(
        tmp_path,
        [matchup_id, unavailable_id],
        allow_missing_matchups=True,
    )
    assert partial_dataset.l2a.shape == dataset.l2a.shape
    with pytest.raises(FileNotFoundError, match="missing pair audit"):
        load_pairs(
            tmp_path,
            [matchup_id, "Missing_Audit"],
            allow_missing_matchups=True,
        )


def test_daily_harmonizer_caps_requested_solver_bands() -> None:
    band_models = {
        band: {
            "mean": [0.0] * len(feature_names("compact", index)),
            "scale": [1.0] * len(feature_names("compact", index)),
            "coef": [0.0] * len(feature_names("compact", index)),
            "intercept": 0.01,
            "feature_names": feature_names("compact", index),
        }
        for index, band in enumerate(BAND_NAMES)
    }
    artifact = {
        "models": {
            "test": {
                "feature_set": "compact",
                "full": band_models,
                "folds": {"0": band_models},
            }
        }
    }
    fetched = {
        "surface": np.full((len(BAND_NAMES), 2, 2), 0.1, dtype=np.float32),
        "l2a_aot": np.full((2, 2), 0.1, dtype=np.float32),
        "l2a_tcwv": np.full((2, 2), 1.0, dtype=np.float32),
        "valid": np.ones((2, 2), dtype=bool),
    }
    scene = {
        "maiac_aot": 0.3,
        "maiac_tcwv_cm": 1.5,
        "sza_deg": 30.0,
        "vza_deg": 5.0,
        "raa_deg": 90.0,
        "elevation_km": 0.2,
        "month": 6,
    }

    corrected, stats = _correct_scene(
        fetched,
        scene,
        model=artifact,
        candidate=("test", "solver", 0.005),
        model_scope="crossfit",
        fold=0,
    )

    np.testing.assert_allclose(corrected[[0, 4, 5, 6]], 0.1)
    np.testing.assert_allclose(corrected[[1, 2, 3]], 0.105)
    assert stats["bands"]["blue"]["cap_fraction"] == 1.0

    gated, gated_stats = _correct_scene(
        fetched,
        scene,
        model=artifact,
        candidate=("test", "solver", 0.005),
        model_scope="crossfit",
        fold=0,
        application_maiac_aot_max=0.15,
    )

    np.testing.assert_allclose(gated, 0.1)
    assert gated_stats["mapping_applied"] is False
    assert gated_stats["skip_reason"] == "outside_training_maiac_aot_domain"

    blue, blue_stats = _correct_scene(
        fetched,
        scene,
        model=artifact,
        candidate=("test", "blue", 0.005),
        model_scope="crossfit",
        fold=0,
    )

    np.testing.assert_allclose(blue[1], 0.105)
    np.testing.assert_allclose(blue[[0, 2, 3, 4, 5, 6]], 0.1)
    assert set(blue_stats["bands"]) == {"blue"}
    assert _candidate("test:blue:0.030") == ("test", "blue", 0.03)

    all_bands, all_stats = _correct_scene(
        fetched,
        scene,
        model=artifact,
        candidate=("test", "all", 0.005),
        model_scope="crossfit",
        fold=0,
    )

    np.testing.assert_allclose(all_bands, 0.105)
    assert set(all_stats["bands"]) == set(BAND_NAMES)


def test_history_uses_uncorrected_single_source_fallback_without_pair_metadata(
    tmp_path, monkeypatch
) -> None:
    matchup_id = "No_Pairs__T00XXX_20240101T000000"
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    comp = np.arange(3 * len(BAND_NAMES) * 2 * 2, dtype=np.float32).reshape(
        3, len(BAND_NAMES), 2, 2
    )
    np.savez_compressed(
        source_dir / f"{matchup_id}.npz",
        comp=comp,
        epsg=np.asarray(32630),
        transform=np.arange(6, dtype=float),
        realizations=np.asarray(["2023-06", "2025-06", "2026-06"]),
        scene_year=np.asarray(2024),
        month=np.asarray(6),
    )
    monkeypatch.setattr(harmonized_histories, "FALLBACK_HISTORY", source_dir)
    outputs = {
        "identity": tmp_path / "out" / "identity.npz",
        "candidate": tmp_path / "out" / "candidate.npz",
    }

    status = harmonized_histories._write_uncorrected_fallback(
        matchup_id,
        outputs=outputs,
        audit_path=tmp_path / "out" / "audit.json",
        model={
            "fold_by_matchup_id": {matchup_id: 5},
            "training_cutoff": "2023-12-31",
            "target": "physical target",
        },
        cutoff="2023-12-31",
        model_scope="crossfit",
        skip_reason="no_exact_pair_atmospheric_support",
        trigger_details={"harmonized_windows": ["2020-09"]},
    )

    assert status["status"] == "OK"
    assert status["mapping_applied"] is False
    assert status["fallback_trigger"] == {"harmonized_windows": ["2020-09"]}
    for path in outputs.values():
        with np.load(path, allow_pickle=False) as history:
            np.testing.assert_allclose(history["comp"], comp[:1])
            provenance = json.loads(str(history["harmonization_json"].item()))
        assert provenance["skip_reason"] == "no_exact_pair_atmospheric_support"
        assert provenance["low_temporal_support"] is True


def test_daily_terrain_harmonizer_uses_local_solar_incidence() -> None:
    terrain = TerrainFields(
        elevation_m=np.asarray([[120.0, 130.0], [140.0, 150.0]], dtype=np.float32),
        slope_deg=np.asarray([[0.0, 1.2], [2.3, 1.2]], dtype=np.float32),
        gradient_east=np.asarray([[0.0, 0.02], [0.04, -0.02]], dtype=np.float32),
        gradient_north=np.zeros((2, 2), dtype=np.float32),
        valid=np.ones((2, 2), dtype=bool),
    )
    feature_count = len(feature_names("terrain", 1))
    incidence_index = feature_names("terrain", 1).index("terrain_incidence_delta")
    band_models = {}
    for index, band in enumerate(BAND_NAMES):
        coefficient = np.zeros(len(feature_names("terrain", index)), dtype=float)
        if band == "blue":
            coefficient[incidence_index] = 1.0
        band_models[band] = {
            "mean": [0.0] * feature_count,
            "scale": [1.0] * feature_count,
            "coef": coefficient.tolist(),
            "intercept": 0.0,
            "feature_names": feature_names("terrain", index),
        }
    artifact = {
        "models": {
            "terrain": {
                "feature_set": "terrain",
                "full": band_models,
                "folds": {"0": band_models},
            }
        }
    }
    fetched = {
        "surface": np.full((len(BAND_NAMES), 2, 2), 0.1, dtype=np.float32),
        "l2a_aot": np.full((2, 2), 0.1, dtype=np.float32),
        "l2a_tcwv": np.full((2, 2), 1.0, dtype=np.float32),
        "valid": np.ones((2, 2), dtype=bool),
    }
    scene = {
        "maiac_aot": 0.3,
        "maiac_tcwv_cm": 1.5,
        "sza_deg": 30.0,
        "saa_deg": 90.0,
        "vza_deg": 5.0,
        "raa_deg": 90.0,
        "elevation_km": 0.2,
        "month": 6,
        "l2a": {"spacecraft": "Sentinel-2A", "processing_baseline": "05.00"},
    }

    corrected, stats = _correct_scene(
        fetched,
        scene,
        model=artifact,
        candidate=("terrain", "solver", 0.03),
        model_scope="crossfit",
        fold=0,
        terrain=terrain,
    )

    incidence = local_solar_incidence(terrain, sza_deg=30.0, saa_deg=90.0)
    expected_blue = 0.1 + np.clip(incidence - np.cos(np.radians(30.0)), -0.03, 0.03)
    np.testing.assert_allclose(corrected[1], expected_blue, rtol=0.0, atol=2.0e-7)
    np.testing.assert_allclose(corrected[[0, 2, 3, 4, 5, 6]], 0.1)
    assert stats["terrain"]["enabled"] is True
    assert stats["terrain"]["median_slope_deg"] == np.median(terrain.slope_deg)


def test_paired_history_fetch_retries_transient_errors(monkeypatch) -> None:
    calls = 0
    expected = (np.ones((7, 2, 2), dtype=np.float32), np.ones((2, 2), dtype=bool))

    def fake_fetch(_ee, _grid, _scene):
        nonlocal calls
        calls += 1
        if calls < 3:
            raise OSError("temporary")
        return expected

    monkeypatch.setattr(paired_histories, "_fetch_target", fake_fetch)
    monkeypatch.setattr(paired_histories.time, "sleep", lambda _seconds: None)

    actual = paired_histories._fetch_with_retries(None, {}, {}, attempts=3)

    assert calls == 3
    assert actual is expected


def test_scene_harmonizer_features_are_finite_and_complete() -> None:
    l2a = np.tile(
        np.asarray([0.02, 0.03, 0.04, 0.05, 0.2, 0.12, 0.08]),
        (120, 1),
    )
    features = build_scene_features(
        l2a,
        np.linspace(0.1, 0.2, 120),
        np.linspace(1.0, 1.5, 120),
        {
            "maiac_aot": 0.3,
            "maiac_tcwv_cm": 1.8,
            "sza_deg": 30.0,
            "vza_deg": 5.0,
            "raa_deg": 90.0,
            "elevation_km": 0.4,
            "month": 6,
            "l2a": {
                "spacecraft": "Sentinel-2B",
                "processing_baseline": "05.11",
            },
        },
    )

    assert features.shape == (len(scene_feature_names()),)
    assert np.all(np.isfinite(features))


def test_scene_harmonizer_solver_mode_preserves_anchor_bands() -> None:
    class Model:
        def predict(self, _features):
            return np.asarray([[0.01] * len(BAND_NAMES)])

    fetched = {
        "surface": np.full((len(BAND_NAMES), 2, 2), 0.1, dtype=np.float32),
        "l2a_aot": np.full((2, 2), 0.1, dtype=np.float32),
        "l2a_tcwv": np.full((2, 2), 1.0, dtype=np.float32),
        "valid": np.ones((2, 2), dtype=bool),
    }
    scene = {
        "maiac_aot": 0.2,
        "maiac_tcwv_cm": 1.5,
        "sza_deg": 30.0,
        "vza_deg": 5.0,
        "raa_deg": 90.0,
        "elevation_km": 0.2,
        "month": 6,
        "l2a": {"spacecraft": "Sentinel-2A", "processing_baseline": "05.00"},
    }

    corrected, stats = _apply_scene_correction(
        fetched, scene, model=Model(), candidate=("solver", 0.03)
    )

    np.testing.assert_allclose(corrected[[0, 4, 5, 6]], 0.1)
    np.testing.assert_allclose(corrected[[1, 2, 3]], 0.11)
    assert stats["valid_pixels"] == 4


def test_same_acquisition_affine_fit_is_robust_and_bounded() -> None:
    source = np.linspace(0.02, 0.22, 500)
    target = 0.01 + 0.85 * source
    target[0] = 0.7  # A bad paired pixel should be trimmed by the robust fit.

    band_map = fit_band_map(source, target, min_samples=100)

    assert band_map["valid"] is True
    assert band_map["method"] == "affine"
    assert band_map["slope"] == pytest.approx(0.85, abs=0.01)
    assert band_map["intercept"] == pytest.approx(0.01, abs=0.002)
    corrected, correction = apply_band_map(np.asarray([0.02, 0.2]), band_map, cap=0.01)
    np.testing.assert_allclose(correction, [0.007, -0.01])
    np.testing.assert_allclose(corrected, [0.027, 0.19])


def test_same_acquisition_low_spread_map_falls_back_to_offset() -> None:
    source = np.full(240, 0.1)
    target = source + 0.012

    band_map = fit_band_map(source, target, min_samples=100)

    assert band_map["valid"] is True
    assert band_map["method"] == "offset_fallback"
    assert band_map["slope"] == 1.0
    assert band_map["intercept"] == pytest.approx(0.012)


def test_same_acquisition_archive_maps_and_solver_mode(tmp_path) -> None:
    rng = np.random.default_rng(4)
    source = rng.uniform(0.03, 0.25, size=(300, len(BAND_NAMES))).astype(np.float32)
    target = (0.01 + 0.9 * source).astype(np.float32)
    matchup_id = "Example_Site__T30ABC_20240101T000000"
    np.savez_compressed(
        tmp_path / f"{matchup_id}.npz",
        l2a=source,
        siac=target,
        scene_index=np.zeros(source.shape[0], dtype=np.int16),
        band_names=np.asarray(BAND_NAMES),
        scenes_json=np.asarray(
            json.dumps(
                [
                    {
                        "scene_id": "S2A_30ABC_20200101_0_L1C",
                        "day": "2020-01-01",
                        "window": "2020-01",
                        "maiac_aot": 0.3,
                    }
                ]
            )
        ),
    )
    (tmp_path / f"{matchup_id}.json").write_text(
        json.dumps({"status": "OK", "uses_aeronet": False}), encoding="utf-8"
    )

    maps = build_scene_maps(tmp_path / f"{matchup_id}.npz", cutoff="2023-12-31")
    key = scene_key({"scene_id": "S2A_30ABC_20200101_0_L1C"})
    assert maps[key]["bands"]["blue"]["valid"] is True
    with pytest.raises(ValueError, match="no paired acquisitions"):
        build_scene_maps(
            tmp_path / f"{matchup_id}.npz",
            cutoff="2023-12-31",
            maiac_aot_max=0.15,
        )

    fetched = {
        "surface": np.full((len(BAND_NAMES), 2, 2), 0.1, dtype=np.float32),
        "valid": np.ones((2, 2), dtype=bool),
    }
    corrected, stats = _correct_pair_affine_scene(fetched, maps[key], candidate=("solver", 0.03))

    np.testing.assert_allclose(corrected[[0, 4, 5, 6]], 0.1)
    np.testing.assert_allclose(corrected[[1, 2, 3]], 0.1)
    assert stats["mapping_applied"] is True
    assert set(stats["bands"]) == {"blue", "green", "red"}


def test_daily_sklearn_harmonizer_applies_fold_models(tmp_path) -> None:
    import joblib
    from sklearn.ensemble import HistGradientBoostingRegressor
    from tools.aeronet_validation.train_l2a_l1c_nonlinear_harmonizer import (
        feature_names as nonlinear_feature_names,
    )

    names = nonlinear_feature_names()
    rng = np.random.default_rng(0)
    train_x = rng.normal(size=(64, len(names)))
    train_y = np.full(64, 0.01)
    model_dir = tmp_path / "models"
    model_dir.mkdir()
    model_files = {}
    for band in ("blue", "green", "red"):
        estimator = HistGradientBoostingRegressor(max_iter=1, min_samples_leaf=1, random_state=0)
        estimator.fit(train_x, train_y)
        path = model_dir / f"hgb_fold5_{band}.joblib"
        joblib.dump(estimator, path)
        model_files[band] = path.name
    artifact = {
        "model_type": "hist_gradient_boosting",
        "model_name": "hgb_full_state",
        "feature_names_by_band": dict.fromkeys(("blue", "green", "red"), names),
        "model_files": {"5": model_files},
    }
    terrain = TerrainFields(
        elevation_m=np.asarray([[120.0, 130.0], [140.0, 150.0]], dtype=np.float32),
        slope_deg=np.asarray([[0.0, 1.2], [2.3, 1.2]], dtype=np.float32),
        gradient_east=np.zeros((2, 2), dtype=np.float32),
        gradient_north=np.zeros((2, 2), dtype=np.float32),
        valid=np.ones((2, 2), dtype=bool),
    )
    fetched = {
        "surface": np.full((len(BAND_NAMES), 2, 2), 0.1, dtype=np.float32),
        "l2a_aot": np.full((2, 2), 0.1, dtype=np.float32),
        "l2a_tcwv": np.full((2, 2), 1.0, dtype=np.float32),
        "valid": np.ones((2, 2), dtype=bool),
    }
    scene = {
        "maiac_aot": 0.3,
        "maiac_tcwv_cm": 1.5,
        "sza_deg": 30.0,
        "saa_deg": 90.0,
        "vza_deg": 5.0,
        "raa_deg": 90.0,
        "elevation_km": 0.2,
        "month": 6,
        "l2a": {"spacecraft": "Sentinel-2A", "processing_baseline": "05.00"},
    }

    corrected, stats = _correct_scene(
        fetched,
        scene,
        model=artifact,
        candidate=("hgb_full_state", "solver", 0.005),
        model_scope="crossfit",
        fold=5,
        terrain=terrain,
        artifact_dir=tmp_path,
    )

    np.testing.assert_allclose(corrected[[1, 2, 3]], 0.105)
    np.testing.assert_allclose(corrected[[0, 4, 5, 6]], 0.1)
    assert stats["mapping_applied"] is True
    assert stats["bands"]["blue"]["cap_fraction"] == 1.0
    assert stats["terrain"]["enabled"] is True

    with pytest.raises(ValueError, match="does not contain model"):
        _correct_scene(
            fetched,
            scene,
            model=artifact,
            candidate=("other_model", "solver", 0.005),
            model_scope="crossfit",
            fold=5,
            terrain=terrain,
            artifact_dir=tmp_path,
        )


def test_daily_sklearn_harmonizer_supplies_scene_state_columns(tmp_path) -> None:
    import joblib
    from sklearn.ensemble import HistGradientBoostingRegressor
    from tools.aeronet_validation.train_l2a_l1c_nonlinear_harmonizer import (
        feature_names as nonlinear_feature_names,
    )

    names = [*nonlinear_feature_names("target", "blue"), "ozone_du_cams"]
    rng = np.random.default_rng(0)
    estimator = HistGradientBoostingRegressor(max_iter=1, min_samples_leaf=1, random_state=0)
    estimator.fit(rng.normal(size=(64, len(names))), np.full(64, 0.01))
    model_dir = tmp_path / "models"
    model_dir.mkdir()
    joblib.dump(estimator, model_dir / "hgb_fold5_blue.joblib")
    artifact = {
        "model_type": "hist_gradient_boosting",
        "model_name": "hgb_target_band_cams_o3",
        "feature_names_by_band": {"blue": names},
        "model_files": {"5": {"blue": "hgb_fold5_blue.joblib"}},
    }
    terrain = TerrainFields(
        elevation_m=np.full((2, 2), 100.0, dtype=np.float32),
        slope_deg=np.zeros((2, 2), dtype=np.float32),
        gradient_east=np.zeros((2, 2), dtype=np.float32),
        gradient_north=np.zeros((2, 2), dtype=np.float32),
        valid=np.ones((2, 2), dtype=bool),
    )
    fetched = {
        "surface": np.full((len(BAND_NAMES), 2, 2), 0.1, dtype=np.float32),
        "l2a_aot": np.full((2, 2), 0.1, dtype=np.float32),
        "l2a_tcwv": np.full((2, 2), 1.0, dtype=np.float32),
        "valid": np.ones((2, 2), dtype=bool),
    }
    scene = {
        "maiac_aot": 0.3,
        "maiac_tcwv_cm": 1.5,
        "sza_deg": 30.0,
        "saa_deg": 90.0,
        "vza_deg": 5.0,
        "raa_deg": 90.0,
        "elevation_km": 0.2,
        "month": 6,
        "l2a": {"spacecraft": "Sentinel-2A", "processing_baseline": "05.00"},
    }

    corrected, stats = _correct_scene(
        fetched,
        scene,
        model=artifact,
        candidate=("hgb_target_band_cams_o3", "blue", 0.005),
        model_scope="crossfit",
        fold=5,
        terrain=terrain,
        artifact_dir=tmp_path,
        scene_state={"ozone_du_cams": 312.0},
    )
    np.testing.assert_allclose(corrected[1], 0.105)
    assert stats["bands"]["blue"]["cap_fraction"] == 1.0

    with pytest.raises(ValueError, match="unknown feature column"):
        _correct_scene(
            fetched,
            scene,
            model=artifact,
            candidate=("hgb_target_band_cams_o3", "blue", 0.005),
            model_scope="crossfit",
            fold=5,
            terrain=terrain,
            artifact_dir=tmp_path,
        )
