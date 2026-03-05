"""Unit tests for cloud/cloud-shadow class module."""

from __future__ import annotations

from types import ModuleType
from typing import TYPE_CHECKING

import numpy as np
import pytest
import xarray as xr

from siac.cloud.mapping import (
    _ensure_expected_classes,
    apply_class_mapping,
    validate_class_mapping,
)
from siac.cloud.mask import (
    _call_user_callable,
    _external_classes,
    _extract_band,
    _group_band_names,
    _mean_group,
    _prepare_rgbnir,
    _resample_classes,
    _resample_continuous,
    build_cloud_classes,
    classes_to_bool_mask,
)
from siac.cloud.providers.omnicloudmask import OmniCloudMaskProvider
from siac.core.types import SensorBand, SensorConfig

if TYPE_CHECKING:
    from pathlib import Path


def _sensor_config_with_duplicate_red() -> SensorConfig:
    return SensorConfig(
        sensor_id="HS",
        satellite_id="TEST",
        bands=(
            SensorBand("G1", 545.0, 10.0, 5.0, 0),
            SensorBand("R1", 650.0, 10.0, 5.0, 1),
            SensorBand("R2", 670.0, 10.0, 5.0, 2),
            SensorBand("N1", 800.0, 20.0, 5.0, 3),
        ),
    )


def _toa_hs(shape: tuple[int, int] = (3, 3)) -> xr.Dataset:
    y, x = shape
    return xr.Dataset(
        {
            "G1": xr.DataArray(np.full((y, x), 0.2, dtype=np.float32), dims=["y", "x"]),
            "R1": xr.DataArray(np.full((y, x), 0.3, dtype=np.float32), dims=["y", "x"]),
            "R2": xr.DataArray(np.full((y, x), 0.5, dtype=np.float32), dims=["y", "x"]),
            "N1": xr.DataArray(np.full((y, x), 0.6, dtype=np.float32), dims=["y", "x"]),
        }
    )


def _spatial_da(value: float, shape: tuple[int, int] = (8, 8)) -> xr.DataArray:
    x = np.linspace(500000.0, 500000.0 + (shape[1] - 1) * 5.0, shape[1])
    y = np.linspace(4500000.0 + (shape[0] - 1) * 5.0, 4500000.0, shape[0])
    return xr.DataArray(
        np.full(shape, value, dtype=np.float32),
        dims=["y", "x"],
        coords={"x": x, "y": y},
    ).rio.write_crs("EPSG:32632")


def test_validate_class_mapping_allows_many_to_one_and_rejects_one_to_many():
    valid = validate_class_mapping({2: [8, 9, 10], 3: [3], 1: [4, 5, 6, 7]})
    assert valid[2] == [8, 9, 10]

    with pytest.raises(ValueError, match="maps to multiple"):
        validate_class_mapping({2: [8], 3: [8]})



def test_apply_class_mapping_identity_and_unmapped_behaviour():
    src = xr.DataArray(np.array([[0, 1, 2, 3, 9]], dtype=np.int16), dims=["y", "x"])

    out = apply_class_mapping(src, None, unmapped_to_missing=True)
    assert out.dtype == np.uint8
    np.testing.assert_array_equal(out.values, np.array([[0, 1, 2, 3, 0]], dtype=np.uint8))

    with pytest.raises(ValueError, match="unmapped"):
        apply_class_mapping(src, None, unmapped_to_missing=False)



def test_apply_class_mapping_explicit_mapping():
    src = xr.DataArray(np.array([[255, 1, 2, 10, 3]], dtype=np.int16), dims=["y", "x"])
    mapping = {0: [255], 1: [1], 2: [2, 10], 3: [3]}
    out = apply_class_mapping(src, mapping, unmapped_to_missing=False)
    np.testing.assert_array_equal(out.values, np.array([[0, 1, 2, 2, 3]], dtype=np.uint8))



def test_classes_to_bool_mask():
    classes = xr.DataArray(np.array([[0, 1, 2, 3]], dtype=np.uint8), dims=["y", "x"])
    mask = classes_to_bool_mask(classes)
    np.testing.assert_array_equal(mask.values, np.array([[True, False, True, True]], dtype=bool))



def test_prepare_rgbnir_averages_multiple_red_bands():
    cfg = _sensor_config_with_duplicate_red()
    toa = _toa_hs()

    red, green, nir = _prepare_rgbnir(
        toa,
        cfg,
        target_resolution_m=10.0,
        resolution_policy="auto",
        allow_upsample_to_target=False,
    )

    assert red.shape == green.shape == nir.shape
    assert float(red.mean()) == pytest.approx(0.4)



def test_group_band_names_raises_when_missing_color():
    cfg = _sensor_config_with_duplicate_red()
    toa = xr.Dataset({"R1": xr.DataArray(np.ones((2, 2), dtype=np.float32), dims=["y", "x"])})

    with pytest.raises(ValueError, match="green"):
        _group_band_names(toa, cfg, "green")



def test_resample_policy_branches(monkeypatch: pytest.MonkeyPatch):
    da = _spatial_da(0.3, shape=(10, 10))
    calls: list[tuple[str, float]] = []

    def _fake_resample(data, target_resolution, resampling="bilinear"):  # noqa: ANN001
        calls.append((resampling, float(target_resolution)))
        return data

    monkeypatch.setattr("siac.cloud.mask.resample", _fake_resample)

    # Finer than target: auto should downsample.
    _ = _resample_continuous(
        da,
        target_resolution_m=10.0,
        resolution_policy="auto",
        allow_upsample_to_target=False,
    )
    assert calls and calls[-1] == ("bilinear", 10.0)

    calls.clear()
    # Force mode always resamples when different.
    _ = _resample_classes(
        da,
        target_resolution_m=15.0,
        resolution_policy="force",
        allow_upsample_to_target=False,
    )
    assert calls and calls[-1] == ("nearest", 15.0)



def test_build_cloud_classes_none_and_user_callable_modes():
    toa = _toa_hs((2, 2))
    cfg = _sensor_config_with_duplicate_red()

    none_out = build_cloud_classes(toa, cfg, mode="none")
    assert set(np.unique(none_out.values)).issubset({1})

    def _simple_callable(_toa):
        return xr.DataArray(np.array([[8, 9], [3, 4]], dtype=np.int16), dims=["y", "x"])

    mapped = build_cloud_classes(
        toa,
        cfg,
        mode="user_callable",
        user_callable=_simple_callable,
        class_mapping={2: [8, 9], 3: [3], 1: [4]},
        unmapped_to_missing=False,
    )
    np.testing.assert_array_equal(mapped.values, np.array([[2, 2], [3, 1]], dtype=np.uint8))



def test_build_cloud_classes_external_file_mode(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    toa = _toa_hs((2, 2))
    cfg = _sensor_config_with_duplicate_red()

    external = tmp_path / "mask.tif"
    external.write_text("x")

    raw = xr.DataArray(np.array([[8, 3], [1, 255]], dtype=np.int16), dims=["y", "x"])
    monkeypatch.setattr("siac.cloud.mask.read_raster", lambda _path: raw)

    out = build_cloud_classes(
        toa,
        cfg,
        mode="external_file",
        external_mask_path=external,
        class_mapping={2: [8], 3: [3], 1: [1], 0: [255]},
        unmapped_to_missing=False,
    )
    np.testing.assert_array_equal(out.values, np.array([[2, 3], [1, 0]], dtype=np.uint8))



def test_build_cloud_classes_auto_mode_uses_provider(monkeypatch: pytest.MonkeyPatch):
    toa = _toa_hs((2, 2))
    cfg = _sensor_config_with_duplicate_red()

    class _FakeProvider:
        def predict(self, red, green, nir, class_mapping=None, unmapped_to_missing=True):  # noqa: ANN001
            assert red.shape == green.shape == nir.shape
            # returns already standardized classes
            return xr.DataArray(np.array([[1, 2], [3, 0]], dtype=np.uint8), dims=["y", "x"])

    monkeypatch.setattr("siac.cloud.mask.OmniCloudMaskProvider", _FakeProvider)

    out = build_cloud_classes(toa, cfg, mode="auto")
    np.testing.assert_array_equal(out.values, np.array([[1, 2], [3, 0]], dtype=np.uint8))



def test_omnicloud_provider_predictor_paths():
    red = xr.DataArray(np.array([[0.1, 0.5]], dtype=np.float32), dims=["y", "x"])
    green = xr.DataArray(np.array([[0.1, 0.5]], dtype=np.float32), dims=["y", "x"])
    nir = xr.DataArray(np.array([[0.1, 0.5]], dtype=np.float32), dims=["y", "x"])

    # Binary predictor: 1=cloud, 0=clear
    p_binary = OmniCloudMaskProvider(
        predictor=lambda arr: (arr[..., 0] > 0.2).astype(np.uint8)
    )
    out_binary = p_binary.predict(red, green, nir)
    np.testing.assert_array_equal(out_binary.values, np.array([[1, 2]], dtype=np.uint8))

    # Predictor with custom labels requiring mapping.
    p_custom = OmniCloudMaskProvider(
        predictor=lambda _arr: np.array([[4, 9]], dtype=np.int16)
    )
    out_custom = p_custom.predict(
        red,
        green,
        nir,
        class_mapping={1: [4], 2: [9]},
        unmapped_to_missing=False,
    )
    np.testing.assert_array_equal(out_custom.values, np.array([[1, 2]], dtype=np.uint8))



def test_build_cloud_classes_input_validation_errors():
    cfg = _sensor_config_with_duplicate_red()

    with pytest.raises(ValueError, match="at least one band"):
        build_cloud_classes(xr.Dataset(), cfg)

    toa = _toa_hs((2, 2))

    with pytest.raises(ValueError, match="external_mask_path"):
        build_cloud_classes(toa, cfg, mode="external_file")

    with pytest.raises(ValueError, match="user_callable"):
        build_cloud_classes(toa, cfg, mode="user_callable")

    with pytest.raises(ValueError, match="Unsupported cloud mode"):
        build_cloud_classes(toa, cfg, mode="bad")


def test_mapping_additional_validation_paths():
    normalized = validate_class_mapping({2: np.int16(8)})
    assert normalized[2] == [8]

    with pytest.raises(TypeError, match="dictionary"):
        validate_class_mapping("bad")  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="target class"):
        validate_class_mapping({9: [1]})

    with pytest.raises(TypeError, match="int or iterable"):
        validate_class_mapping({2: object()})  # type: ignore[arg-type]


def test_apply_mapping_empty_sources_and_strict_unmapped_error():
    src = xr.DataArray(np.array([[7, 8]], dtype=np.int16), dims=["y", "x"])
    out = apply_class_mapping(src, {2: []}, unmapped_to_missing=True)
    np.testing.assert_array_equal(out.values, np.array([[0, 0]], dtype=np.uint8))

    with pytest.raises(ValueError, match="unmapped values"):
        apply_class_mapping(src, {1: [7]}, unmapped_to_missing=False)


def test_ensure_expected_classes_empty_and_invalid():
    empty = xr.DataArray(np.array([], dtype=np.uint8), dims=["x"])
    assert _ensure_expected_classes(empty) is empty

    with pytest.raises(ValueError, match="must be in"):
        _ensure_expected_classes(xr.DataArray(np.array([4], dtype=np.int16), dims=["x"]))


def test_extract_band_and_resample_error_paths(monkeypatch: pytest.MonkeyPatch):
    stacked = xr.DataArray(np.arange(18, dtype=np.float32).reshape(2, 3, 3), dims=["band", "y", "x"])
    assert _extract_band(stacked).shape == (3, 3)
    assert _extract_band(stacked.isel(band=slice(0, 1))).shape == (3, 3)

    with pytest.raises(ValueError, match="must be > 0"):
        _resample_continuous(stacked.isel(band=0), target_resolution_m=0, resolution_policy="auto", allow_upsample_to_target=False)
    with pytest.raises(ValueError, match="must be > 0"):
        _resample_classes(stacked.isel(band=0), target_resolution_m=0, resolution_policy="auto", allow_upsample_to_target=False)

    calls: list[tuple[str, float]] = []

    def _fake_resample(data, target_resolution, resampling="bilinear"):  # noqa: ANN001
        calls.append((resampling, float(target_resolution)))
        return data

    monkeypatch.setattr("siac.cloud.mask.resample", _fake_resample)
    coarse = _spatial_da(0.2, shape=(5, 5)).assign_coords(
        x=np.linspace(0.0, 80.0, 5),
        y=np.linspace(80.0, 0.0, 5),
    )
    _ = _resample_continuous(
        coarse,
        target_resolution_m=10.0,
        resolution_policy="auto",
        allow_upsample_to_target=True,
    )
    assert calls and calls[-1] == ("bilinear", 10.0)

    calls.clear()
    keep = _resample_continuous(
        coarse,
        target_resolution_m=10.0,
        resolution_policy="auto",
        allow_upsample_to_target=False,
    )
    assert keep is coarse
    assert not calls

    _ = _resample_classes(
        coarse,
        target_resolution_m=10.0,
        resolution_policy="auto",
        allow_upsample_to_target=True,
    )
    assert calls and calls[-1] == ("nearest", 10.0)

    with pytest.raises(ValueError, match="resolution_policy"):
        _resample_continuous(coarse, target_resolution_m=10.0, resolution_policy="bad", allow_upsample_to_target=False)
    with pytest.raises(ValueError, match="resolution_policy"):
        _resample_classes(coarse, target_resolution_m=10.0, resolution_policy="bad", allow_upsample_to_target=False)


def test_mean_group_alignment_reproject_and_interp(monkeypatch: pytest.MonkeyPatch):
    calls = {"reproject": 0}

    def _fake_reproject(src, ref, resampling="bilinear"):  # noqa: ANN001
        calls["reproject"] += 1
        return xr.full_like(ref, float(np.nanmean(src.values)))

    monkeypatch.setattr("siac.cloud.mask.reproject_match", _fake_reproject)

    toa_geo = xr.Dataset(
        {
            "R1": _spatial_da(0.3, shape=(4, 4)),
            "R2": _spatial_da(0.5, shape=(4, 4)).assign_coords(
                x=np.linspace(100.0, 115.0, 4),
                y=np.linspace(215.0, 200.0, 4),
            ),
        }
    )
    out_geo = _mean_group(
        toa_geo,
        ["R1", "R2"],
        target_resolution_m=10.0,
        resolution_policy="auto",
        allow_upsample_to_target=False,
    )
    assert out_geo.ndim == 2
    assert calls["reproject"] >= 1

    toa_non_geo = xr.Dataset(
        {
            "A": xr.DataArray(np.full((4, 4), 0.3, dtype=np.float32), dims=["y", "x"]),
            "B": xr.DataArray(np.full((4, 4), 0.5, dtype=np.float32), dims=["y", "x"]),
        }
    )
    out_non_geo = _mean_group(
        toa_non_geo,
        ["A", "B"],
        target_resolution_m=10.0,
        resolution_policy="auto",
        allow_upsample_to_target=False,
    )
    assert out_non_geo.shape == (4, 4)


def test_prepare_rgbnir_alignment_branches(monkeypatch: pytest.MonkeyPatch):
    cfg = _sensor_config_with_duplicate_red()
    toa = _toa_hs((3, 3))

    monkeypatch.setattr("siac.cloud.mask._group_band_names", lambda *_args, **_kwargs: ["R1"])  # noqa: ARG005
    monkeypatch.setattr(
        "siac.cloud.mask._mean_group",
        lambda *_args, **_kwargs: _spatial_da(0.2, shape=(3, 3)).assign_coords(
            x=np.linspace(0.0, 20.0, 3),
            y=np.linspace(20.0, 0.0, 3),
        ),
    )
    rep_calls = {"n": 0}

    def _fake_reproject(src, ref, resampling="bilinear"):  # noqa: ANN001
        rep_calls["n"] += 1
        return xr.full_like(ref, float(np.nanmean(src.values)))

    monkeypatch.setattr("siac.cloud.mask.reproject_match", _fake_reproject)
    _ = _prepare_rgbnir(
        toa,
        cfg,
        target_resolution_m=10.0,
        resolution_policy="auto",
        allow_upsample_to_target=False,
    )
    assert rep_calls["n"] >= 2

    seq = iter(
        [
            xr.DataArray(
                np.full((4, 4), 0.2, dtype=np.float32),
                dims=["y", "x"],
                coords={"y": np.arange(4), "x": np.arange(4)},
            ),
            xr.DataArray(
                np.full((2, 2), 0.3, dtype=np.float32),
                dims=["y", "x"],
                coords={"y": np.array([0, 3]), "x": np.array([0, 3])},
            ),
            xr.DataArray(
                np.full((3, 3), 0.4, dtype=np.float32),
                dims=["y", "x"],
                coords={"y": np.array([0, 2, 3]), "x": np.array([0, 2, 3])},
            ),
        ]
    )
    monkeypatch.setattr("siac.cloud.mask._mean_group", lambda *_args, **_kwargs: next(seq))
    red, green, nir = _prepare_rgbnir(
        toa,
        cfg,
        target_resolution_m=10.0,
        resolution_policy="auto",
        allow_upsample_to_target=False,
    )
    assert red.shape == (4, 4)
    assert green.shape == red.shape
    assert nir.shape == red.shape


def test_external_classes_and_user_callable_fallbacks(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    raw = _spatial_da(8.0, shape=(2, 2)).astype(np.int16)
    ref = _spatial_da(0.2, shape=(2, 2))
    calls = {"reproject": 0}

    monkeypatch.setattr("siac.cloud.mask.read_raster", lambda _path: raw)

    def _fake_reproject(src, target, resampling="nearest"):  # noqa: ANN001
        calls["reproject"] += 1
        return xr.full_like(target, 8, dtype=np.int16)

    monkeypatch.setattr("siac.cloud.mask.reproject_match", _fake_reproject)
    out = _external_classes(
        tmp_path / "mask.tif",
        reference=ref,
        class_mapping={2: [8]},
        target_resolution_m=5.0,
        resolution_policy="auto",
        allow_upsample_to_target=False,
        unmapped_to_missing=False,
    )
    assert calls["reproject"] == 1
    np.testing.assert_array_equal(out.values, np.array([[2, 2], [2, 2]], dtype=np.uint8))

    def _kw_only(*, toa, sensor_config, input_path):  # noqa: ANN001
        _ = (toa, sensor_config, input_path)
        return np.array([[1]], dtype=np.int16)

    def _positional(toa, sensor_config, input_path):  # noqa: ANN001
        _ = (toa, sensor_config, input_path)
        return np.array([[1]], dtype=np.int16)

    def _single(toa):  # noqa: ANN001
        _ = toa
        return np.array([[1]], dtype=np.int16)

    sample = _toa_hs((2, 2))
    cfg = _sensor_config_with_duplicate_red()
    assert _call_user_callable(_kw_only, toa=sample, sensor_config=cfg, input_path=None).shape == (1, 1)
    assert _call_user_callable(_positional, toa=sample, sensor_config=cfg, input_path=None).shape == (1, 1)
    assert _call_user_callable(_single, toa=sample, sensor_config=cfg, input_path=None).shape == (1, 1)


def test_build_cloud_classes_user_callable_shape_reproject_and_provider_error(
    monkeypatch: pytest.MonkeyPatch,
):
    toa = _toa_hs((2, 2))
    cfg = _sensor_config_with_duplicate_red()

    monkeypatch.setattr(
        "siac.cloud.mask.reproject_match",
        lambda _src, target, **_kwargs: xr.full_like(target, 1, dtype=np.int16),
    )

    out = build_cloud_classes(
        toa,
        cfg,
        mode="user_callable",
        user_callable=lambda *_args, **_kwargs: np.array([[1]], dtype=np.int16),
        class_mapping={1: [1]},
        unmapped_to_missing=False,
    )
    np.testing.assert_array_equal(out.values, np.ones((2, 2), dtype=np.uint8))

    with pytest.raises(ValueError, match="Unsupported cloud provider"):
        build_cloud_classes(toa, cfg, mode="auto", provider="unsupported-provider")


def test_omnicloud_default_predictor_and_normalize_paths(monkeypatch: pytest.MonkeyPatch):
    # OmniCloudMask class + predict method
    module_predict = ModuleType("omnicloudmask")

    class _ModelPredict:
        def predict(self, arr):  # noqa: ANN001
            return (arr[..., 0] > 0).astype(np.uint8)

    module_predict.OmniCloudMask = _ModelPredict  # type: ignore[attr-defined]
    monkeypatch.setitem(__import__("sys").modules, "omnicloudmask", module_predict)
    assert callable(OmniCloudMaskProvider()._default_predictor())

    # OmniCloudMask class + __call__
    module_call = ModuleType("omnicloudmask")

    class _ModelCall:
        def __call__(self, arr):  # noqa: ANN001
            return (arr[..., 0] > 0).astype(np.uint8)

    module_call.OmniCloudMask = _ModelCall  # type: ignore[attr-defined]
    monkeypatch.setitem(__import__("sys").modules, "omnicloudmask", module_call)
    assert callable(OmniCloudMaskProvider()._default_predictor())

    # top-level predict function
    module_fn = ModuleType("omnicloudmask")
    module_fn.predict = lambda arr: (arr[..., 0] > 0).astype(np.uint8)  # type: ignore[attr-defined]
    monkeypatch.setitem(__import__("sys").modules, "omnicloudmask", module_fn)
    assert callable(OmniCloudMaskProvider()._default_predictor())

    # module with no known entrypoint -> None
    module_none = ModuleType("omnicloudmask")
    monkeypatch.setitem(__import__("sys").modules, "omnicloudmask", module_none)
    assert OmniCloudMaskProvider()._default_predictor() is None

    red = xr.DataArray(np.array([[0.1, 0.2]], dtype=np.float32), dims=["y", "x"])
    green = xr.DataArray(np.array([[0.1, 0.2]], dtype=np.float32), dims=["y", "x"])
    nir = xr.DataArray(np.array([[0.1, 0.2]], dtype=np.float32), dims=["y", "x"])

    # Heuristic fallback path
    provider = OmniCloudMaskProvider(predictor=None)
    monkeypatch.setattr(provider, "_default_predictor", lambda: None)
    out = provider.predict(red, green, nir)
    assert out.shape == (1, 2)

    # Shape mismatch in predict
    with pytest.raises(ValueError, match="identical shape"):
        provider.predict(red, green, xr.DataArray(np.array([[0.1], [0.2]], dtype=np.float32), dims=["y", "x"]))

    # 3D logits/probabilities branch and shape mismatch error
    norm = OmniCloudMaskProvider._normalize_raw_output(
        np.array([[[0.1, 0.9], [0.8, 0.2]]], dtype=np.float32),
        red,
    )
    assert norm.shape == red.shape
    with pytest.raises(ValueError, match="does not match input shape"):
        OmniCloudMaskProvider._normalize_raw_output(np.zeros((2, 2), dtype=np.float32), red)
