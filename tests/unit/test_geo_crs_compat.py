"""Tests for siac.geo._crs_compat and the call sites in readers / reprojection.

Covers REVIEW.md §2.8 and §3.7: authority-based CRS comparison and
dtype-aware default resampling for ``check_rasters_aligned``,
``reproject_match``, ``reproject_dataset_match``, ``resample``,
``resample_to_shape``, and ``compute_common_bounds``.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
import rasterio
import xarray as xr
from pyproj import CRS as PyprojCRS
from rasterio.enums import Resampling

import siac.geo.reprojection as reproj_mod
import siac.storage.readers as readers_mod
from siac.geo._crs_compat import (
    crs_equivalent,
    default_resampling_for_dtype,
)

# =============================================================================
# crs_equivalent
# =============================================================================


def test_crs_equivalent_handles_none() -> None:
    assert crs_equivalent(None, None) is True
    assert crs_equivalent(None, "EPSG:4326") is False
    assert crs_equivalent("EPSG:4326", None) is False


def test_crs_equivalent_epsg_string_matches_canonical_wkt() -> None:
    """The whole point of the helper: EPSG:4326 == verbose WKT for 4326."""
    wkt_4326 = PyprojCRS.from_epsg(4326).to_wkt()
    assert wkt_4326 != "EPSG:4326"  # sanity: different strings
    assert crs_equivalent("EPSG:4326", wkt_4326) is True


def test_crs_equivalent_distinguishes_different_crs() -> None:
    assert crs_equivalent("EPSG:4326", "EPSG:32632") is False
    assert crs_equivalent("EPSG:4326", "EPSG:4269") is False


def test_crs_equivalent_accepts_pyproj_and_rasterio_crs() -> None:
    """Mixed CRS objects from pyproj/rasterio should compare equal."""
    pp = PyprojCRS.from_epsg(4326)
    rio = rasterio.CRS.from_epsg(4326)
    assert crs_equivalent(pp, rio) is True
    assert crs_equivalent(pp, "EPSG:4326") is True
    assert crs_equivalent(rio, "EPSG:4326") is True


def test_crs_equivalent_falls_back_to_string_match_on_garbage() -> None:
    """If both inputs are non-CRS strings the helper degrades to ==."""
    assert crs_equivalent("not a crs", "not a crs") is True
    assert crs_equivalent("not a crs", "different garbage") is False


# =============================================================================
# default_resampling_for_dtype
# =============================================================================


@pytest.mark.parametrize(
    "dtype",
    [np.uint8, np.uint16, np.int8, np.int16, np.int32, np.int64, np.bool_],
)
def test_default_resampling_integer_dtypes_use_nearest(dtype) -> None:
    assert default_resampling_for_dtype(dtype) is Resampling.nearest


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_default_resampling_float_dtypes_use_bilinear(dtype) -> None:
    assert default_resampling_for_dtype(dtype) is Resampling.bilinear


def test_default_resampling_unknown_dtype_returns_bilinear() -> None:
    """A garbage dtype falls through to bilinear (legacy behaviour)."""

    class _NotADtype:
        pass

    assert default_resampling_for_dtype(_NotADtype()) is Resampling.bilinear


# =============================================================================
# check_rasters_aligned uses authority comparison
# =============================================================================


def test_check_rasters_aligned_treats_epsg_string_and_wkt_as_equal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wkt_4326 = PyprojCRS.from_epsg(4326).to_wkt()

    def _info(path):  # noqa: ANN001
        # Same raster, same bounds/resolution — only the CRS encoding
        # differs between the two paths.
        crs_repr = "EPSG:4326" if str(path).endswith("a") else wkt_4326
        return {
            "crs": crs_repr,
            "resolution": (0.001, 0.001),
            "bounds": SimpleNamespace(left=0.0, bottom=0.0, right=1.0, top=1.0),
        }

    monkeypatch.setattr(readers_mod, "get_raster_info", _info)
    # Pre-fix this returned False because of string equality.
    assert readers_mod.check_rasters_aligned("a", "b") is True


def test_check_rasters_aligned_still_rejects_genuinely_different_crs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _info(path):  # noqa: ANN001
        return {
            "crs": "EPSG:32632" if str(path).endswith("a") else "EPSG:4326",
            "resolution": (10.0, 10.0),
            "bounds": SimpleNamespace(left=0.0, bottom=0.0, right=1.0, top=1.0),
        }

    monkeypatch.setattr(readers_mod, "get_raster_info", _info)
    assert readers_mod.check_rasters_aligned("a", "b") is False


# =============================================================================
# read_zarr_array remote-scheme detection (gs://, azure://, abfs://, file://)
# =============================================================================


@pytest.mark.parametrize(
    "url",
    [
        "http://example.com/data.zarr",
        "https://example.com/data.zarr",
        "s3://bucket/data.zarr",
        "gs://bucket/data.zarr",
        "azure://container/data.zarr",
        "abfs://container/data.zarr",
        "file:///tmp/data.zarr",
    ],
)
def test_read_zarr_array_remote_scheme_uses_mapper(
    monkeypatch: pytest.MonkeyPatch, url: str
) -> None:
    captured: dict[str, object] = {}

    def _fake_mapper(path: str) -> dict[str, str]:
        captured["mapper_called_with"] = path
        return {"sentinel": path}

    monkeypatch.setattr(readers_mod, "_get_remote_zarr_mapper", _fake_mapper)
    monkeypatch.setattr(
        readers_mod.xr,
        "open_zarr",
        lambda _obj, chunks=None: xr.Dataset({"v": xr.DataArray([1], dims=["x"])}),  # noqa: ARG005
    )

    readers_mod.read_zarr_array(url)
    assert captured["mapper_called_with"] == url


def test_read_zarr_array_local_path_does_not_use_mapper(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def _fake_mapper(path: str) -> dict[str, str]:
        captured["mapper_called"] = True
        return {}

    monkeypatch.setattr(readers_mod, "_get_remote_zarr_mapper", _fake_mapper)
    monkeypatch.setattr(
        readers_mod.xr,
        "open_zarr",
        lambda _obj, chunks=None: xr.Dataset({"v": xr.DataArray([1], dims=["x"])}),  # noqa: ARG005
    )

    readers_mod.read_zarr_array("/tmp/local.zarr")
    assert "mapper_called" not in captured


# =============================================================================
# reprojection helpers
# =============================================================================


def _spatial_da(
    shape: tuple[int, int] = (8, 8),
    crs: str = "EPSG:32632",
    dtype: np.dtype = np.float32,
) -> xr.DataArray:
    data = np.arange(shape[0] * shape[1], dtype=dtype).reshape(shape)
    x = np.linspace(500000.0, 500000.0 + (shape[1] - 1) * 10.0, shape[1])
    y = np.linspace(4500000.0 + (shape[0] - 1) * 10.0, 4500000.0, shape[0])
    return xr.DataArray(data, dims=["y", "x"], coords={"x": x, "y": y}).rio.write_crs(crs)


def test_resolve_resampling_unknown_string_raises() -> None:
    with pytest.raises(ValueError, match="Unknown resampling method"):
        reproj_mod._resolve_resampling("cubicc")


def test_resolve_resampling_unknown_string_lists_valid_methods() -> None:
    with pytest.raises(ValueError, match="bilinear"):
        # Error message advertises the valid set.
        reproj_mod._resolve_resampling("not-a-method")


def test_resolve_resampling_none_dispatches_on_dtype() -> None:
    assert reproj_mod._resolve_resampling(None, dtype=np.uint8) is Resampling.nearest
    assert reproj_mod._resolve_resampling(None, dtype=np.float32) is Resampling.bilinear


def test_resolve_resampling_passes_through_resampling_enum() -> None:
    assert reproj_mod._resolve_resampling(Resampling.cubic) is Resampling.cubic


def test_resolve_resampling_bad_type_raises_typeerror() -> None:
    with pytest.raises(TypeError, match="resampling must be"):
        reproj_mod._resolve_resampling(42)  # type: ignore[arg-type]


def test_reproject_match_default_uses_nearest_for_integer(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Integer mask must NOT be silently bilinear-resampled."""
    src = _spatial_da((6, 6), dtype=np.uint8)
    target = _spatial_da((4, 4), dtype=np.float32)

    captured: dict[str, object] = {}

    def _fake_match(self, target, resampling=None, **kwargs):  # noqa: ANN001
        captured["resampling"] = resampling
        return self

    monkeypatch.setattr(type(src.rio), "reproject_match", _fake_match, raising=True)

    reproj_mod.reproject_match(src, target)  # resampling defaults to None
    assert captured["resampling"] is Resampling.nearest


def test_reproject_match_default_uses_bilinear_for_float(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    src = _spatial_da((6, 6), dtype=np.float32)
    target = _spatial_da((4, 4), dtype=np.float32)

    captured: dict[str, object] = {}

    def _fake_match(self, target, resampling=None, **kwargs):  # noqa: ANN001
        captured["resampling"] = resampling
        return self

    monkeypatch.setattr(type(src.rio), "reproject_match", _fake_match, raising=True)

    reproj_mod.reproject_match(src, target)
    assert captured["resampling"] is Resampling.bilinear


def test_reproject_match_unknown_string_raises() -> None:
    src = _spatial_da((4, 4), dtype=np.float32)
    target = _spatial_da((4, 4), dtype=np.float32)
    with pytest.raises(ValueError, match="Unknown resampling method"):
        reproj_mod.reproject_match(src, target, resampling="bilinearr")


def test_reproject_dataset_match_picks_per_variable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Mask variable should get nearest, reflectance should get bilinear."""
    target = _spatial_da((4, 4), dtype=np.float32)

    mask = _spatial_da((6, 6), dtype=np.uint8).rename("mask")
    refl = _spatial_da((6, 6), dtype=np.float32).rename("refl")
    ds = xr.Dataset({"mask": mask, "refl": refl}).rio.write_crs("EPSG:32632")

    # Pair up dtype kind -> resampling so we can verify what each variable
    # got without relying on rio accessor having a ``.name`` attr.
    captured: list[tuple[str, Resampling]] = []

    def _fake_match(self, target, resampling=None, **kwargs):  # noqa: ANN001
        captured.append((self._obj.dtype.kind, resampling))
        return self._obj

    monkeypatch.setattr(type(target.rio), "reproject_match", _fake_match, raising=True)

    reproj_mod.reproject_dataset_match(ds, target)  # default None → per-var
    by_kind = dict(captured)
    assert by_kind["u"] is Resampling.nearest  # mask
    assert by_kind["f"] is Resampling.bilinear  # refl


def test_reproject_dataset_match_explicit_resampling_applies_to_all(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = _spatial_da((4, 4), dtype=np.float32)
    mask = _spatial_da((6, 6), dtype=np.uint8).rename("mask")
    refl = _spatial_da((6, 6), dtype=np.float32).rename("refl")
    ds = xr.Dataset({"mask": mask, "refl": refl}).rio.write_crs("EPSG:32632")

    captured: list[Resampling] = []

    def _fake_match(self, target, resampling=None, **kwargs):  # noqa: ANN001
        captured.append(resampling)
        return self._obj

    monkeypatch.setattr(type(target.rio), "reproject_match", _fake_match, raising=True)

    reproj_mod.reproject_dataset_match(ds, target, resampling="cubic")
    assert captured == [Resampling.cubic, Resampling.cubic]


def test_reproject_to_crs_default_uses_dtype_aware_resampling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    src = _spatial_da((6, 6), dtype=np.uint8)
    captured: dict[str, object] = {}

    def _fake_reproject(self, target_crs, **kwargs):  # noqa: ANN001
        captured.update(kwargs)
        return self

    monkeypatch.setattr(type(src.rio), "reproject", _fake_reproject, raising=True)

    reproj_mod.reproject_to_crs(src, "EPSG:4326")  # resampling=None default
    assert captured["resampling"] is Resampling.nearest


def test_resample_default_uses_dtype_aware_resampling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    src = _spatial_da((10, 10), dtype=np.int16)
    captured: dict[str, object] = {}

    def _fake_reproject(self, target_crs, **kwargs):  # noqa: ANN001
        captured.update(kwargs)
        return self

    monkeypatch.setattr(type(src.rio), "reproject", _fake_reproject, raising=True)

    reproj_mod.resample(src, target_resolution=20.0)  # default None
    assert captured["resampling"] is Resampling.nearest


def test_resample_uses_round_not_truncation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """``int(round(...))`` covers the original extent; ``int(...)`` truncates."""
    src = _spatial_da((10, 10), dtype=np.float32)
    captured_shapes: dict[str, tuple[int, int]] = {}

    def _fake_reproject(self, target_crs, *, shape, resampling, **kwargs):  # noqa: ANN001
        captured_shapes["shape"] = shape
        return self

    monkeypatch.setattr(type(src.rio), "reproject", _fake_reproject, raising=True)

    # Source is 10×10 at 10 m → 100 m extent.
    # Target resolution 30 m → scale = 10/30 ≈ 0.3333…
    # 10 * 0.3333 = 3.333 → truncating gives 3, rounding gives 3 too.
    # Use a scale that demonstrates rounding: 10 * 0.55 = 5.5 → trunc 5, round 6.
    # So pick target_resolution such that scale = 0.55: 10/x = 0.55 -> x ≈ 18.18
    reproj_mod.resample(src, target_resolution=18.18)
    h, w = captured_shapes["shape"]
    # scale = 10/18.18 ≈ 0.5500..., 10*0.55 = 5.5 → round → 6
    assert (h, w) == (6, 6)


def test_resample_to_shape_default_uses_dtype_aware_resampling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    src = _spatial_da((10, 10), dtype=np.uint16)
    captured: dict[str, object] = {}

    def _fake_reproject(self, target_crs, **kwargs):  # noqa: ANN001
        captured.update(kwargs)
        return self

    monkeypatch.setattr(type(src.rio), "reproject", _fake_reproject, raising=True)

    reproj_mod.resample_to_shape(src, target_shape=(5, 5))  # default None
    assert captured["resampling"] is Resampling.nearest


def test_compute_common_bounds_treats_equivalent_crs_as_same(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Two arrays whose CRS are authority-equivalent must not trigger a
    bounds transform.

    The case that the helper actually defends against: a pyproj ``CRS`` and
    a string encoding of the same authority. Rasterio's own ``__eq__`` is
    already authority-aware for objects from the same library, but
    cross-library comparisons (pyproj.CRS vs rasterio.CRS, or strings) used
    to fall through to plain string equality. ``_crs_equivalent`` makes the
    branch decision robust to the encoding.
    """
    a = _spatial_da((4, 4), crs="EPSG:32632", dtype=np.float32)
    b = _spatial_da((4, 4), crs="EPSG:32632", dtype=np.float32)

    calls: list[object] = []

    def _spy(bounds, src, dst):  # noqa: ANN001
        calls.append((bounds, src, dst))
        return bounds

    monkeypatch.setattr(reproj_mod, "transform_bounds", _spy)
    reproj_mod.compute_common_bounds(a, b, method="union")

    # Equivalent CRS → no bounds transform.
    assert calls == []


def test_compute_common_bounds_uses_helper_for_crs_check(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Direct evidence that the function delegates to ``_crs_equivalent``."""
    a = _spatial_da((4, 4), crs="EPSG:32632", dtype=np.float32)
    b = _spatial_da((4, 4), crs="EPSG:32632", dtype=np.float32)

    seen: list[tuple[object, object]] = []

    def _spy_eq(lhs, rhs):  # noqa: ANN001
        seen.append((lhs, rhs))
        return True

    monkeypatch.setattr(reproj_mod, "_crs_equivalent", _spy_eq)
    monkeypatch.setattr(
        reproj_mod, "transform_bounds", lambda *_a, **_k: pytest.fail("should not transform")
    )
    reproj_mod.compute_common_bounds(a, b, method="union")
    assert len(seen) >= 1  # at least one CRS comparison happened


def test_compute_common_bounds_still_transforms_for_real_crs_difference(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    a = _spatial_da((4, 4), crs="EPSG:32632", dtype=np.float32)
    b = _spatial_da((4, 4), crs="EPSG:4326", dtype=np.float32)

    calls: list[object] = []

    def _spy(bounds, src, dst):  # noqa: ANN001
        calls.append((bounds, src, dst))
        return bounds

    monkeypatch.setattr(reproj_mod, "transform_bounds", _spy)
    reproj_mod.compute_common_bounds(a, b, method="union")

    assert len(calls) == 1  # only b transforms; a is the reference
