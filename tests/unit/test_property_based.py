"""Property-based tests using Hypothesis for core resampling and mask functions."""

from __future__ import annotations

import numpy as np
import xarray as xr
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from siac.geo.resample import (
    resample_field_to_template,
    resample_mask_to_template,
    shares_template_grid,
    should_resample_for_policy,
)


def _make_da(shape: tuple[int, int], fill: float = 0.5) -> xr.DataArray:
    """Create a simple (y, x) DataArray with coordinates."""
    h, w = shape
    return xr.DataArray(
        np.full(shape, fill, dtype=np.float32),
        dims=["y", "x"],
        coords={
            "y": np.arange(h, dtype=np.float64),
            "x": np.arange(w, dtype=np.float64),
        },
    )


@given(
    h=st.integers(min_value=1, max_value=64),
    w=st.integers(min_value=1, max_value=64),
)
@settings(max_examples=30, deadline=5000)
def test_shares_template_grid_reflexive(h: int, w: int):
    """Any DataArray shares a grid with itself."""
    da = _make_da((h, w))
    assert shares_template_grid(da, da)


@given(
    src_h=st.integers(min_value=2, max_value=32),
    src_w=st.integers(min_value=2, max_value=32),
    dst_h=st.integers(min_value=2, max_value=32),
    dst_w=st.integers(min_value=2, max_value=32),
)
@settings(max_examples=30, deadline=5000)
def test_resample_field_output_shape(src_h: int, src_w: int, dst_h: int, dst_w: int):
    """Resampled field always has the target shape."""
    source = _make_da((src_h, src_w), fill=0.3)
    template = _make_da((dst_h, dst_w))
    result = resample_field_to_template(source, template)
    assert result.shape == template.shape


@given(
    src_h=st.integers(min_value=2, max_value=32),
    src_w=st.integers(min_value=2, max_value=32),
    dst_h=st.integers(min_value=2, max_value=32),
    dst_w=st.integers(min_value=2, max_value=32),
)
@settings(max_examples=30, deadline=5000)
def test_resample_mask_output_shape(src_h: int, src_w: int, dst_h: int, dst_w: int):
    """Resampled mask always has the target shape."""
    source = _make_da((src_h, src_w), fill=1.0).astype(bool)
    template = _make_da((dst_h, dst_w))
    result = resample_mask_to_template(source, template)
    assert result.shape == template.shape
    assert result.dtype == bool


@given(
    h=st.integers(min_value=2, max_value=32),
    w=st.integers(min_value=2, max_value=32),
)
@settings(max_examples=20, deadline=5000)
def test_resample_field_preserves_finite_values(h: int, w: int):
    """Resampling a finite field should produce finite output."""
    source = _make_da((h, w), fill=0.5)
    template = _make_da((h, w))  # same size
    result = resample_field_to_template(source, template)
    assert np.all(np.isfinite(result.values))


@given(
    h=st.integers(min_value=2, max_value=32),
    w=st.integers(min_value=2, max_value=32),
)
@settings(max_examples=20, deadline=5000)
def test_resample_all_true_mask_stays_true(h: int, w: int):
    """An all-True mask should remain all-True after resampling."""
    source = xr.DataArray(
        np.ones((h, w), dtype=bool),
        dims=["y", "x"],
        coords={"y": np.arange(h, dtype=np.float64), "x": np.arange(w, dtype=np.float64)},
    )
    template = _make_da((h, w))
    result = resample_mask_to_template(source, template)
    assert np.all(result.values)


@given(
    current=st.floats(min_value=1.0, max_value=1000.0),
    target=st.floats(min_value=1.0, max_value=1000.0),
)
@settings(max_examples=50, deadline=5000)
def test_should_resample_force_symmetric(current: float, target: float):
    """Force policy: resample iff resolutions differ."""
    assume(abs(current - target) > 1e-6 or current == target)
    result = should_resample_for_policy(current, target, policy="force")
    expected = abs(current - target) > 1e-6
    assert result == expected
