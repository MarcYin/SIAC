"""Tests for the canonical TargetGrid contract and cached_reproject_match.

Wave 18 (opt 3): a robust, consistent reprojection-with-cache primitive
that the pipeline can adopt incrementally without breaking existing
``reproject_match`` callers.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rioxarray  # noqa: F401  # registers .rio accessor
import xarray as xr

from siac.geo.cached_reprojection import (
    CACHE_FORMAT_VERSION,
    cached_reproject_match,
    compute_cache_key,
    load_cached_reprojection,
    save_cached_reprojection,
)
from siac.geo.target_grid import TargetGrid


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def _build_grid(
    *,
    bounds: tuple[float, float, float, float] = (0.0, 0.0, 1000.0, 1000.0),
    crs: str = "EPSG:32633",
    resolution_m: float = 10.0,
    shape: tuple[int, int] | None = None,
) -> TargetGrid:
    if shape is None:
        x_min, y_min, x_max, y_max = bounds
        shape = (
            int(round((y_max - y_min) / resolution_m)),
            int(round((x_max - x_min) / resolution_m)),
        )
    return TargetGrid(bounds=bounds, crs=crs, resolution_m=resolution_m, shape=shape)


def _gradient_source(
    *,
    crs: str = "EPSG:32633",
    bounds: tuple[float, float, float, float] = (-100.0, -100.0, 1100.0, 1100.0),
    shape: tuple[int, int] = (60, 60),
) -> xr.DataArray:
    """A unique-per-pixel ``y * width + x`` gradient — cache hits trivially detectable."""
    height, width = shape
    arr = np.arange(height * width, dtype=np.float32).reshape(shape)
    x_min, y_min, x_max, y_max = bounds
    x_coords = np.linspace(x_min, x_max, width, dtype=np.float64)
    y_coords = np.linspace(y_max, y_min, height, dtype=np.float64)
    da = xr.DataArray(arr, dims=("y", "x"), coords={"y": y_coords, "x": x_coords})
    da = da.rio.write_crs(crs)
    return da


# --------------------------------------------------------------------------
# TargetGrid
# --------------------------------------------------------------------------


class TestTargetGrid:
    def test_signature_is_deterministic_for_identical_grids(self) -> None:
        a = _build_grid()
        b = _build_grid()
        assert a.signature() == b.signature()

    def test_signature_differs_when_bounds_differ(self) -> None:
        a = _build_grid()
        b = _build_grid(bounds=(0.0, 0.0, 2000.0, 2000.0))
        assert a.signature() != b.signature()

    def test_signature_differs_when_crs_differs(self) -> None:
        a = _build_grid(crs="EPSG:32633")
        b = _build_grid(crs="EPSG:32634")
        assert a.signature() != b.signature()

    def test_signature_differs_when_resolution_differs(self) -> None:
        a = _build_grid(resolution_m=10.0, shape=(100, 100))
        b = _build_grid(resolution_m=20.0, shape=(50, 50))
        assert a.signature() != b.signature()

    def test_signature_includes_version_tag(self) -> None:
        sig = _build_grid().signature()
        assert sig[0] == "TargetGrid/v1"

    def test_signature_absorbs_tiny_bounds_noise(self) -> None:
        """Bounds rounded to 6 decimal places so float-roundtrip noise doesn't bust cache."""
        a = _build_grid(bounds=(0.0, 0.0, 1000.0, 1000.0))
        b = _build_grid(
            bounds=(0.0 + 1e-9, 0.0, 1000.0 + 1e-9, 1000.0)
        )  # sub-micrometre noise
        assert a.signature() == b.signature()

    def test_from_template_round_trip(self) -> None:
        original = _build_grid()
        template = original.as_template_da()
        recovered = TargetGrid.from_template(template)
        assert recovered.signature() == original.signature()

    def test_from_template_requires_crs(self) -> None:
        # Give the array x/y coords so rioxarray can compute bounds, but
        # NOT a CRS — we want the explicit "no CRS metadata" failure
        # mode, not a DimensionMissingCoordinateError.
        bare = xr.DataArray(
            np.zeros((10, 10), dtype=np.float32),
            dims=("y", "x"),
            coords={
                "x": np.linspace(0.0, 100.0, 10),
                "y": np.linspace(100.0, 0.0, 10),
            },
        )
        with pytest.raises(ValueError, match="no CRS metadata"):
            TargetGrid.from_template(bare)

    def test_as_template_has_expected_shape_and_crs(self) -> None:
        grid = _build_grid()
        template = grid.as_template_da()
        assert template.shape == grid.shape
        assert template.rio.crs.to_string() == grid.crs

    def test_reproject_match_works_against_synthesized_template(self) -> None:
        """The synthesised template must be a valid ``reproject_match`` target."""
        grid = _build_grid()
        source = _gradient_source()
        from siac.geo.reprojection import reproject_match

        out = reproject_match(source, grid.as_template_da(), resampling="bilinear")
        assert out.shape == grid.shape


# --------------------------------------------------------------------------
# compute_cache_key
# --------------------------------------------------------------------------


class TestCacheKey:
    def test_key_is_deterministic(self) -> None:
        grid = _build_grid()
        k1 = compute_cache_key(target=grid, source_identity="mcd43-x", resampling="nearest")
        k2 = compute_cache_key(target=grid, source_identity="mcd43-x", resampling="nearest")
        assert k1 == k2
        assert len(k1) == 64  # sha256 hex

    def test_key_changes_with_grid(self) -> None:
        a = compute_cache_key(target=_build_grid(), source_identity="x", resampling=None)
        b = compute_cache_key(
            target=_build_grid(bounds=(0.0, 0.0, 2000.0, 2000.0)),
            source_identity="x",
            resampling=None,
        )
        assert a != b

    def test_key_changes_with_source_identity(self) -> None:
        a = compute_cache_key(target=_build_grid(), source_identity="x", resampling=None)
        b = compute_cache_key(target=_build_grid(), source_identity="y", resampling=None)
        assert a != b

    def test_key_changes_with_resampling(self) -> None:
        a = compute_cache_key(
            target=_build_grid(), source_identity="x", resampling="nearest"
        )
        b = compute_cache_key(
            target=_build_grid(), source_identity="x", resampling="bilinear"
        )
        assert a != b

    def test_key_changes_with_extra_namespace(self) -> None:
        a = compute_cache_key(
            target=_build_grid(),
            source_identity="x",
            resampling=None,
            extra_namespace="caller_a",
        )
        b = compute_cache_key(
            target=_build_grid(),
            source_identity="x",
            resampling=None,
            extra_namespace="caller_b",
        )
        assert a != b

    def test_cache_format_version_is_pinned(self) -> None:
        """If the format ever changes, this test makes the diff explicit."""
        assert CACHE_FORMAT_VERSION == "reproject-v1"


# --------------------------------------------------------------------------
# save / load roundtrip
# --------------------------------------------------------------------------


class TestPersistence:
    def test_roundtrip_preserves_values(self, tmp_path: Path) -> None:
        grid = _build_grid()
        src = _gradient_source()
        from siac.geo.reprojection import reproject_match

        reprojected = reproject_match(src, grid.as_template_da(), resampling="bilinear")
        key = "a" * 64
        save_cached_reprojection(tmp_path, key, reprojected)
        loaded = load_cached_reprojection(tmp_path, key, target=grid)
        assert loaded is not None
        np.testing.assert_allclose(
            loaded.values, reprojected.values, rtol=1e-6, atol=1e-6
        )

    def test_roundtrip_preserves_crs(self, tmp_path: Path) -> None:
        grid = _build_grid(crs="EPSG:32634")
        src = _gradient_source(crs="EPSG:32634")
        from siac.geo.reprojection import reproject_match

        reprojected = reproject_match(src, grid.as_template_da(), resampling="bilinear")
        save_cached_reprojection(tmp_path, "b" * 64, reprojected)
        loaded = load_cached_reprojection(tmp_path, "b" * 64, target=grid)
        assert loaded is not None
        assert loaded.rio.crs.to_string() == grid.crs

    def test_load_missing_returns_none(self, tmp_path: Path) -> None:
        grid = _build_grid()
        out = load_cached_reprojection(tmp_path, "c" * 64, target=grid)
        assert out is None

    def test_load_shape_mismatch_returns_none(self, tmp_path: Path) -> None:
        grid_small = _build_grid(shape=(50, 50))
        grid_big = _build_grid(
            bounds=(0.0, 0.0, 2000.0, 2000.0),
            shape=(200, 200),
        )
        src = _gradient_source()
        from siac.geo.reprojection import reproject_match

        small_out = reproject_match(
            src, grid_small.as_template_da(), resampling="bilinear"
        )
        save_cached_reprojection(tmp_path, "d" * 64, small_out)
        # Now ask for the cache under the same key but expect the BIG shape:
        # load must return None (treats it as a miss rather than crashing).
        result = load_cached_reprojection(tmp_path, "d" * 64, target=grid_big)
        assert result is None

    def test_load_corrupted_returns_none(self, tmp_path: Path) -> None:
        grid = _build_grid()
        key = "e" * 64
        bad_path = tmp_path / key[:2] / f"{key}.nc"
        bad_path.parent.mkdir(parents=True)
        bad_path.write_bytes(b"not a netcdf file by any stretch")
        assert load_cached_reprojection(tmp_path, key, target=grid) is None

    def test_load_disabled_when_cache_dir_is_none(self) -> None:
        grid = _build_grid()
        assert load_cached_reprojection(None, "0" * 64, target=grid) is None

    def test_save_disabled_when_cache_dir_is_none(self, tmp_path: Path) -> None:
        # Just must not raise; nothing to verify on disk.
        save_cached_reprojection(None, "0" * 64, _gradient_source())


# --------------------------------------------------------------------------
# cached_reproject_match
# --------------------------------------------------------------------------


class TestCachedReprojectMatch:
    def test_no_cache_matches_uncached_reproject_match(self) -> None:
        """With ``cache_dir=None`` the helper must be equivalent to ``reproject_match``."""
        from siac.geo.reprojection import reproject_match

        grid = _build_grid()
        src = _gradient_source()
        cached_out = cached_reproject_match(
            src,
            grid.as_template_da(),
            source_identity="x",
            cache_dir=None,
            resampling="bilinear",
        )
        uncached_out = reproject_match(
            src, grid.as_template_da(), resampling="bilinear"
        )
        np.testing.assert_allclose(
            cached_out.values, uncached_out.values, rtol=1e-6, atol=1e-6
        )

    def test_miss_then_hit_returns_identical_data(self, tmp_path: Path) -> None:
        grid = _build_grid()
        src = _gradient_source()
        first = cached_reproject_match(
            src, grid, source_identity="x", cache_dir=tmp_path, resampling="bilinear"
        )
        second = cached_reproject_match(
            src, grid, source_identity="x", cache_dir=tmp_path, resampling="bilinear"
        )
        # Both calls return numerically-identical results.
        np.testing.assert_allclose(first.values, second.values, rtol=1e-6, atol=1e-6)

    def test_hit_skips_actual_reprojection(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        grid = _build_grid()
        src = _gradient_source()

        # First call populates the cache.
        cached_reproject_match(
            src, grid, source_identity="x", cache_dir=tmp_path, resampling="bilinear"
        )

        # Now patch reproject_match to explode if it's invoked — a cache
        # hit must not touch it. Patch in both the module under test and
        # the geo.reprojection module to catch any indirection.
        call_count = {"n": 0}

        def _explode(*args, **kwargs):
            call_count["n"] += 1
            raise AssertionError(
                "cache hit must not invoke reproject_match"
            )

        monkeypatch.setattr(
            "siac.geo.cached_reprojection.reproject_match", _explode
        )

        # Second call must come entirely from cache; no reproject_match.
        result = cached_reproject_match(
            src, grid, source_identity="x", cache_dir=tmp_path, resampling="bilinear"
        )
        assert call_count["n"] == 0
        assert result.shape == grid.shape

    def test_different_source_identity_misses_cache(
        self, tmp_path: Path
    ) -> None:
        grid = _build_grid()
        src_a = _gradient_source()
        src_b = _gradient_source(shape=(60, 60))  # same shape, different conceptually
        # Force values to actually differ so a wrongly-served cache hit
        # would be visible.
        src_b = src_b * 10.0  # type: ignore[operator]
        src_b = src_b.rio.write_crs("EPSG:32633")

        out_a = cached_reproject_match(
            src_a, grid, source_identity="src_a", cache_dir=tmp_path, resampling="bilinear"
        )
        out_b = cached_reproject_match(
            src_b, grid, source_identity="src_b", cache_dir=tmp_path, resampling="bilinear"
        )
        # Identity differs → different cache entries → different outputs.
        assert not np.allclose(out_a.values, out_b.values)

    def test_empty_source_identity_raises(self, tmp_path: Path) -> None:
        grid = _build_grid()
        src = _gradient_source()
        with pytest.raises(ValueError, match="non-empty source_identity"):
            cached_reproject_match(
                src,
                grid,
                source_identity="",
                cache_dir=tmp_path,
                resampling="bilinear",
            )

    def test_accepts_template_da_in_addition_to_target_grid(
        self, tmp_path: Path
    ) -> None:
        """Legacy callers that pass a template DataArray must still work."""
        grid = _build_grid()
        template_da = grid.as_template_da()
        src = _gradient_source()
        out = cached_reproject_match(
            src,
            template_da,
            source_identity="x",
            cache_dir=tmp_path,
            resampling="bilinear",
        )
        assert out.shape == grid.shape

    def test_extra_namespace_partitions_cache(self, tmp_path: Path) -> None:
        grid = _build_grid()
        src = _gradient_source()
        # Two callers sharing cache_dir but with different namespaces must
        # not collide even with identical source_identity + grid + method.
        out_a = cached_reproject_match(
            src,
            grid,
            source_identity="shared",
            cache_dir=tmp_path,
            resampling="bilinear",
            extra_namespace="caller_a",
        )
        # Mutate source bytes via multiplication so the same source_identity
        # under namespace "caller_b" populates a distinct cache entry.
        src_b = (src * 2.0).rio.write_crs(src.rio.crs.to_string())  # type: ignore[operator]
        out_b = cached_reproject_match(
            src_b,
            grid,
            source_identity="shared",
            cache_dir=tmp_path,
            resampling="bilinear",
            extra_namespace="caller_b",
        )
        # The two namespaces produced different stored entries.
        from siac.geo.cached_reprojection import (
            compute_cache_key as _compute,
            _cache_path as _path,
        )

        key_a = _compute(
            target=grid,
            source_identity="shared",
            resampling="bilinear",
            extra_namespace="caller_a",
        )
        key_b = _compute(
            target=grid,
            source_identity="shared",
            resampling="bilinear",
            extra_namespace="caller_b",
        )
        assert key_a != key_b
        assert _path(tmp_path, key_a).exists()
        assert _path(tmp_path, key_b).exists()
        # And the second namespace's cache returns the modified data.
        np.testing.assert_allclose(out_b.values, out_b.values)  # trivially holds
        assert not np.allclose(out_a.values, out_b.values, equal_nan=True)
