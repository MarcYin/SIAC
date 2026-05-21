"""Tests for wave 18f preview downsampling and ``include_previews`` opt-out."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pathlib import Path
import xarray as xr

from siac.storage.writers import (
    _downsample_for_preview,
    _field_to_uint8,
    write_cloud_mask_preview,
    write_false_colour_preview,
    write_field_preview,
)


def _band(shape: tuple[int, int], value: float = 0.15) -> xr.DataArray:
    arr = np.full(shape, value, dtype=np.float32)
    return xr.DataArray(arr, dims=("y", "x"), name="band")


class TestDownsample:
    def test_no_op_when_smaller_than_max(self) -> None:
        data = np.ones((200, 300), dtype=np.float32)
        out = _downsample_for_preview(data, max_size_px=2048)
        assert out.shape == data.shape
        assert out is data  # numpy returns the same array (no copy)

    def test_decimates_when_larger(self) -> None:
        data = np.arange(10980 * 10980, dtype=np.float32).reshape(10980, 10980)
        out = _downsample_for_preview(data, max_size_px=2048)
        # Stride was 10980 // 2048 = 5, so the larger axis becomes
        # ceil(10980 / 5) = 2196.
        assert out.shape[0] <= 10980 // 4
        assert out.shape[1] <= 10980 // 4
        assert max(out.shape) <= 10980 // 4

    def test_disabled_when_max_zero(self) -> None:
        data = np.ones((10980, 10980), dtype=np.float32)
        out = _downsample_for_preview(data, max_size_px=0)
        assert out.shape == data.shape

    def test_preserves_dtype_and_ordering(self) -> None:
        data = np.arange(2000 * 2000, dtype=np.float32).reshape(2000, 2000)
        out = _downsample_for_preview(data, max_size_px=500)
        # Stride is 2000 // 500 = 4, so out should be every 4th element.
        np.testing.assert_array_equal(out, data[::4, ::4])
        assert out.dtype == data.dtype


class TestFieldToUint8:
    def test_clip_to_range(self) -> None:
        data = np.array([[-0.5, 0.0], [0.5, 2.0]], dtype=np.float32)
        out = _field_to_uint8(data, 0.0, 1.0)
        # -0.5 → 0; 0.0 → 0; 0.5 → 127; 2.0 → 255.
        np.testing.assert_array_equal(out, np.array([[0, 0], [127, 255]], dtype=np.uint8))

    def test_nans_become_zero(self) -> None:
        data = np.array([[0.5, np.nan], [np.nan, 1.0]], dtype=np.float32)
        out = _field_to_uint8(data, 0.0, 1.0)
        np.testing.assert_array_equal(out, np.array([[127, 0], [0, 255]], dtype=np.uint8))

    def test_mask_overrides_pixels(self) -> None:
        data = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=np.float32)
        mask = np.array([[True, False], [False, True]])
        out = _field_to_uint8(data, 0.0, 1.0, mask=mask)
        np.testing.assert_array_equal(out, np.array([[0, 127], [127, 0]], dtype=np.uint8))

    def test_returns_uint8(self) -> None:
        data = np.ones((10, 10), dtype=np.float64)
        out = _field_to_uint8(data, 0.0, 1.0)
        assert out.dtype == np.uint8


class TestFalseColourPreview:
    def test_downsample_caps_output_size(self, tmp_path: Path) -> None:
        boa = xr.Dataset(
            {
                "B08": _band((4096, 4096), value=0.30),
                "B04": _band((4096, 4096), value=0.10),
                "B03": _band((4096, 4096), value=0.08),
            }
        )
        out_path = tmp_path / "false_colour.png"
        result = write_false_colour_preview(boa, out_path, max_size_px=1024)
        assert result == out_path

        from PIL import Image

        img = Image.open(out_path)
        # Stride is 4096 // 1024 = 4 → output is 1024×1024 (or a hair less).
        assert max(img.size) <= 4096 // 4

    def test_returns_none_when_bands_missing(self, tmp_path: Path) -> None:
        boa = xr.Dataset({"B02": _band((100, 100))})  # only one band
        result = write_false_colour_preview(boa, tmp_path / "missing.png")
        assert result is None


class TestCloudMaskPreview:
    def test_downsample_caps_output_size(self, tmp_path: Path) -> None:
        boa = xr.Dataset(
            {
                "B04": _band((4096, 4096), value=0.10),
                "B03": _band((4096, 4096), value=0.12),
                "B02": _band((4096, 4096), value=0.08),
            }
        )
        cloud = xr.DataArray(np.zeros((4096, 4096), dtype=bool), dims=("y", "x"), name="cloud_mask")
        # Stripe a few cloudy pixels so the overlay code path is exercised.
        cloud.values[100:200, 100:300] = True
        out_path = tmp_path / "cloud_mask.png"
        result = write_cloud_mask_preview(boa, cloud, out_path, max_size_px=1024)
        assert result == out_path
        from PIL import Image

        img = Image.open(out_path)
        assert max(img.size) <= 4096 // 4

    def test_cloud_mask_resizes_when_different_shape(self, tmp_path: Path) -> None:
        boa = xr.Dataset(
            {
                "B04": _band((2048, 2048), value=0.10),
                "B03": _band((2048, 2048), value=0.12),
                "B02": _band((2048, 2048), value=0.08),
            }
        )
        # Cloud mask at a coarser grid (1/8 the resolution) — must be
        # resized to match the (downsampled) BOA preview.
        cloud = xr.DataArray(np.zeros((256, 256), dtype=bool), dims=("y", "x"), name="cloud_mask")
        cloud.values[10:30, 10:50] = True
        out_path = tmp_path / "cloud_mask.png"
        result = write_cloud_mask_preview(boa, cloud, out_path, max_size_px=1024)
        assert result == out_path


class TestFieldPreview:
    def test_downsample_keeps_within_max_size(self, tmp_path: Path) -> None:
        field = xr.DataArray(
            np.random.RandomState(0).uniform(0.05, 0.4, (5000, 5000)).astype(np.float32),
            dims=("y", "x"),
            name="aot",
        )
        out_path = tmp_path / "aot.png"
        result = write_field_preview(
            field,
            out_path,
            vmin=0.0,
            vmax=0.5,
            palette="magma",
            max_size_px=1024,
        )
        assert result == out_path
        from PIL import Image

        img = Image.open(out_path)
        # Output includes a colour bar + label rows so height > main panel.
        # Just check the width is bounded.
        assert img.size[0] <= 5000 // 4


class TestIncludePreviewsOptOut:
    """The output writer's gating: ``include_previews=False`` must skip the
    expensive PNG renders entirely (wave 18f)."""

    def test_default_includes_previews(self) -> None:
        from siac.config.system import OutputDefaultsConfig

        cfg = OutputDefaultsConfig()
        assert cfg.include_previews is True
        assert cfg.preview_max_size_px == 2048

    def test_opt_out_disables_previews(self) -> None:
        from siac.config.system import OutputDefaultsConfig

        cfg = OutputDefaultsConfig(include_previews=False)
        assert cfg.include_previews is False
