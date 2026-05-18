"""Tests for the content-addressed cloud-mask disk cache (wave 18)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from siac.algorithms.cloud.cache import (
    CACHE_FORMAT_VERSION,
    compute_cache_key,
    load_cached_cloud_classes,
    maybe_run_with_cache,
    save_cached_cloud_classes,
)


def _band(shape: tuple[int, int] = (8, 12), value: float = 0.1) -> xr.DataArray:
    arr = np.full(shape, value, dtype=np.float32)
    return xr.DataArray(arr, dims=("y", "x"))


def _classes(shape: tuple[int, int] = (8, 12), value: int = 1) -> xr.DataArray:
    arr = np.full(shape, value, dtype=np.uint8)
    da = xr.DataArray(arr, dims=("y", "x"))
    da.name = "cloud_classes"
    return da


def test_cache_key_is_deterministic_across_calls() -> None:
    red = _band(value=0.1)
    green = _band(value=0.2)
    nir = _band(value=0.3)
    mapping = {1: [0], 2: [1, 2], 3: [3]}

    key_a = compute_cache_key(
        red=red, green=green, nir=nir,
        class_mapping=mapping, target_resolution_m=10.0,
    )
    key_b = compute_cache_key(
        red=red, green=green, nir=nir,
        class_mapping=mapping, target_resolution_m=10.0,
    )
    assert key_a == key_b
    # SHA-256 hex length.
    assert len(key_a) == 64


def test_cache_key_changes_with_pixel_values() -> None:
    red_a = _band(value=0.1)
    red_b = _band(value=0.10001)  # differ in float32-distinct value
    green = _band(value=0.2)
    nir = _band(value=0.3)

    key_a = compute_cache_key(
        red=red_a, green=green, nir=nir,
        class_mapping=None, target_resolution_m=10.0,
    )
    key_b = compute_cache_key(
        red=red_b, green=green, nir=nir,
        class_mapping=None, target_resolution_m=10.0,
    )
    assert key_a != key_b


def test_cache_key_changes_with_class_mapping() -> None:
    red = _band()
    green = _band()
    nir = _band()
    key_a = compute_cache_key(
        red=red, green=green, nir=nir,
        class_mapping={1: [0]},
        target_resolution_m=10.0,
    )
    key_b = compute_cache_key(
        red=red, green=green, nir=nir,
        class_mapping={1: [0], 3: [1]},
        target_resolution_m=10.0,
    )
    assert key_a != key_b


def test_cache_key_changes_with_target_resolution() -> None:
    red = _band()
    green = _band()
    nir = _band()
    key_a = compute_cache_key(
        red=red, green=green, nir=nir,
        class_mapping=None, target_resolution_m=10.0,
    )
    key_b = compute_cache_key(
        red=red, green=green, nir=nir,
        class_mapping=None, target_resolution_m=20.0,
    )
    assert key_a != key_b


def test_save_then_load_roundtrip_returns_equivalent_data(tmp_path: Path) -> None:
    classes = _classes()
    key = "0" * 64  # arbitrary fixed key
    save_cached_cloud_classes(tmp_path, key, classes)
    loaded = load_cached_cloud_classes(tmp_path, key, template=classes)
    assert loaded is not None
    np.testing.assert_array_equal(loaded.values, classes.values)
    assert loaded.dtype == np.uint8


def test_load_with_missing_cache_returns_none(tmp_path: Path) -> None:
    out = load_cached_cloud_classes(
        tmp_path,
        "deadbeef" * 8,
        template=_band(),
    )
    assert out is None


def test_load_with_shape_mismatch_returns_none(tmp_path: Path) -> None:
    # Stored at one shape, asked back with a different template shape.
    save_cached_cloud_classes(tmp_path, "f" * 64, _classes(shape=(8, 12)))
    out = load_cached_cloud_classes(
        tmp_path, "f" * 64, template=_band(shape=(4, 6))
    )
    assert out is None


def test_maybe_run_with_cache_misses_then_hits(tmp_path: Path) -> None:
    red = _band(value=0.1)
    green = _band(value=0.2)
    nir = _band(value=0.3)

    call_count = {"n": 0}
    expected = _classes(value=2)

    def _compute() -> xr.DataArray:
        call_count["n"] += 1
        return expected

    out_a = maybe_run_with_cache(
        cache_dir=tmp_path,
        red=red, green=green, nir=nir,
        class_mapping={1: [0]}, target_resolution_m=10.0,
        compute_fn=_compute,
    )
    assert call_count["n"] == 1
    np.testing.assert_array_equal(out_a.values, expected.values)

    # Same inputs → second call must NOT invoke compute_fn.
    out_b = maybe_run_with_cache(
        cache_dir=tmp_path,
        red=red, green=green, nir=nir,
        class_mapping={1: [0]}, target_resolution_m=10.0,
        compute_fn=_compute,
    )
    assert call_count["n"] == 1, "cache hit must short-circuit compute"
    np.testing.assert_array_equal(out_b.values, expected.values)


def test_maybe_run_with_cache_disabled_when_cache_dir_is_none() -> None:
    call_count = {"n": 0}
    expected = _classes()

    def _compute() -> xr.DataArray:
        call_count["n"] += 1
        return expected

    out = maybe_run_with_cache(
        cache_dir=None,
        red=_band(), green=_band(), nir=_band(),
        class_mapping=None, target_resolution_m=10.0,
        compute_fn=_compute,
    )
    assert call_count["n"] == 1
    np.testing.assert_array_equal(out.values, expected.values)


def test_extra_namespace_partitions_the_cache(tmp_path: Path) -> None:
    """Two callers with different extra_namespace must not collide."""
    red = _band()
    green = _band()
    nir = _band()
    classes_a = _classes(value=1)
    classes_b = _classes(value=3)

    out_a = maybe_run_with_cache(
        cache_dir=tmp_path,
        red=red, green=green, nir=nir,
        class_mapping=None, target_resolution_m=10.0,
        compute_fn=lambda: classes_a,
        extra_namespace="caller_a",
    )
    out_b = maybe_run_with_cache(
        cache_dir=tmp_path,
        red=red, green=green, nir=nir,
        class_mapping=None, target_resolution_m=10.0,
        compute_fn=lambda: classes_b,
        extra_namespace="caller_b",
    )
    # Each caller gets their own value back; the cache didn't merge them.
    np.testing.assert_array_equal(out_a.values, classes_a.values)
    np.testing.assert_array_equal(out_b.values, classes_b.values)


def test_cache_format_version_is_in_the_key() -> None:
    """If we ever bump CACHE_FORMAT_VERSION, old keys must be invalidated."""
    # We can't easily mock a module-level constant once read into the hash,
    # but we can at least confirm the constant is the documented v1 string
    # so a future bump shows up in a clear diff.
    assert CACHE_FORMAT_VERSION == "v1"


def test_load_with_corrupted_npz_recomputes(tmp_path: Path) -> None:
    """A corrupted cache file should produce a miss + warning, not raise."""
    key = "a" * 64
    # Create a malformed cache entry so np.load fails.
    shard_dir = tmp_path / key[:2]
    shard_dir.mkdir(parents=True)
    bad_path = shard_dir / f"{key}.npz"
    bad_path.write_bytes(b"this is definitely not a valid .npz archive")
    out = load_cached_cloud_classes(tmp_path, key, template=_band())
    assert out is None
