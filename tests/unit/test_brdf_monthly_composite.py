from __future__ import annotations

import numpy as np
import rasterio
import xarray as xr

from siac.algorithms.surface.brdf_monthly_composite import (
    MonthlyBestPixelComposite,
    build_monthly_best_pixel_composite,
    build_monthly_best_pixel_kernel_composite,
)
from siac.runtime import BRDFKernelWeights


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


def _georeferenced_reflectance_cube() -> xr.DataArray:
    cube = _reflectance_cube()
    transform = rasterio.transform.from_bounds(399960.0, 4699040.0, 400960.0, 4700040.0, 2, 2)
    return (
        cube.rio.set_spatial_dims(x_dim="x", y_dim="y")
        .rio.write_crs("EPSG:32615")
        .rio.write_transform(transform)
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


def test_monthly_best_pixel_composite_preserves_spatial_metadata() -> None:
    reflectance = _georeferenced_reflectance_cube()
    quality = (
        xr.DataArray(
            np.array(
                [
                    [[3, 2], [2, 1]],
                    [[2, 3], [1, 2]],
                    [[1, 1], [3, 3]],
                ],
                dtype=np.float32,
            ),
            dims=["time", "y", "x"],
            coords={"time": reflectance.time, "y": reflectance.y, "x": reflectance.x},
        )
        .rio.set_spatial_dims(x_dim="x", y_dim="y")
        .rio.write_crs("EPSG:32615")
        .rio.write_transform(reflectance.isel(time=0, band=0, drop=True).rio.transform(recalc=True))
    )

    result = build_monthly_best_pixel_composite(reflectance, quality, year=2024, month=7)

    assert str(result.reflectance.rio.crs) == "EPSG:32615"
    assert result.reflectance.rio.transform(recalc=True) == reflectance.isel(
        time=0, drop=True
    ).rio.transform(recalc=True)
    assert str(result.quality.rio.crs) == "EPSG:32615"
    assert result.quality.rio.transform(recalc=True) == quality.isel(
        time=0, drop=True
    ).rio.transform(recalc=True)


def test_monthly_best_pixel_kernel_composite_preserves_spatial_metadata() -> None:
    reflectance = _georeferenced_reflectance_cube()
    quality = (
        xr.DataArray(
            np.array(
                [
                    [[0.3, 0.2], [0.2, 0.1]],
                    [[0.2, 0.3], [0.1, 0.2]],
                    [[0.1, 0.1], [0.3, 0.3]],
                ],
                dtype=np.float32,
            ),
            dims=["time", "y", "x"],
            coords={"time": reflectance.time, "y": reflectance.y, "x": reflectance.x},
        )
        .rio.set_spatial_dims(x_dim="x", y_dim="y")
        .rio.write_crs("EPSG:32615")
        .rio.write_transform(reflectance.isel(time=0, band=0, drop=True).rio.transform(recalc=True))
    )
    weights = BRDFKernelWeights(
        f0=reflectance,
        f1=reflectance + 1.0,
        f2=reflectance + 2.0,
        f0_unc=reflectance + 3.0,
        f1_unc=reflectance + 4.0,
        f2_unc=reflectance + 5.0,
    )

    result = build_monthly_best_pixel_kernel_composite(weights, quality, year=2024, month=7)

    assert str(result.kernels.f0.rio.crs) == "EPSG:32615"
    assert result.kernels.f0.rio.transform(recalc=True) == reflectance.isel(
        time=0, drop=True
    ).rio.transform(recalc=True)
    assert str(result.quality.rio.crs) == "EPSG:32615"
    assert result.quality.rio.transform(recalc=True) == quality.isel(
        time=0, drop=True
    ).rio.transform(recalc=True)
