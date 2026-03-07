from __future__ import annotations

import numpy as np
import xarray as xr

from siac.priors.surface.brdf_monthly_composite import (
    MonthlyBestPixelComposite,
    build_monthly_best_pixel_composite,
)


def _reflectance_cube() -> xr.DataArray:
    data = np.array(
        [
            [
                [[0.10, 0.20], [0.30, 0.40]],
                [[0.50, 0.60], [0.70, 0.80]],
            ],
            [
                [[0.11, 0.21], [0.31, 0.41]],
                [[0.51, 0.61], [0.71, 0.81]],
            ],
            [
                [[0.12, 0.22], [0.32, 0.42]],
                [[0.52, 0.62], [0.72, 0.82]],
            ],
        ],
        dtype=np.float32,
    )
    return xr.DataArray(
        data,
        dims=["time", "band", "y", "x"],
        coords={
            "time": [0, 1, 2],
            "band": ["B02", "B08"],
            "y": [0, 1],
            "x": [0, 1],
        },
    )


def test_monthly_best_pixel_composite_uses_lowest_quality_pixel() -> None:
    reflectance = _reflectance_cube()
    quality = xr.DataArray(
        np.array(
            [
                [[3, 2], [2, 1]],
                [[2, 3], [1, 2]],
                [[1, 1], [3, 3]],
            ],
            dtype=np.int16,
        ),
        dims=["time", "y", "x"],
        coords={"time": reflectance.time, "y": reflectance.y, "x": reflectance.x},
    )

    result = build_monthly_best_pixel_composite(
        reflectance,
        quality,
        year=2024,
        month=7,
    )

    assert isinstance(result, MonthlyBestPixelComposite)
    assert result.year == 2024
    assert result.month == 7

    expected_index = np.array([[2, 2], [1, 0]], dtype=np.int16)
    np.testing.assert_array_equal(result.sample_index.values, expected_index)

    expected_b02 = np.array([[0.12, 0.22], [0.31, 0.40]], dtype=np.float32)
    expected_b08 = np.array([[0.52, 0.62], [0.71, 0.80]], dtype=np.float32)
    np.testing.assert_allclose(result.reflectance.sel(band="B02").values, expected_b02)
    np.testing.assert_allclose(result.reflectance.sel(band="B08").values, expected_b08)
    np.testing.assert_array_equal(
        result.quality.values,
        np.array([[1, 1], [1, 1]], dtype=np.int16),
    )


def test_monthly_best_pixel_composite_skips_invalid_samples() -> None:
    reflectance = _reflectance_cube().copy()
    reflectance.loc[{"time": 2, "band": "B02", "y": 0, "x": 0}] = np.nan
    quality = xr.DataArray(
        np.array(
            [
                [[3, 2], [2, 1]],
                [[2, 3], [1, 2]],
                [[1, 1], [3, 3]],
            ],
            dtype=np.int16,
        ),
        dims=["time", "y", "x"],
        coords={"time": reflectance.time, "y": reflectance.y, "x": reflectance.x},
    )

    result = build_monthly_best_pixel_composite(reflectance, quality, year=2024, month=7)

    assert int(result.sample_index.sel(y=0, x=0)) == 1
    np.testing.assert_allclose(
        float(result.reflectance.sel(band="B02", y=0, x=0)),
        0.11,
    )
