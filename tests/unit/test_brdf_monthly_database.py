from __future__ import annotations

import numpy as np
import xarray as xr

from siac.algorithms.surface.brdf_monthly_composite import MonthlyBestPixelComposite
from siac.algorithms.surface.brdf_monthly_database import (
    MonthlyCompositeDatabase,
    build_monthly_composite_database,
)


def _make_composite(month_idx: int) -> MonthlyBestPixelComposite:
    bands = ["B02", "B03", "B08", "B11", "B12"]
    cube = np.zeros((len(bands), 2, 1), dtype=np.float32)

    # Query-band values are identical across the two pixels for every month.
    cube[bands.index("B08"), :, 0] = 0.40
    cube[bands.index("B11"), :, 0] = 0.30
    cube[bands.index("B12"), :, 0] = 0.20

    # Visible outputs differ by spectral family and should only be separable
    # through the median-summary context from the 15 monthly composites.
    cube[bands.index("B02"), 0, 0] = 0.10 + 0.001 * month_idx
    cube[bands.index("B03"), 0, 0] = 0.15 + 0.001 * month_idx
    cube[bands.index("B02"), 1, 0] = 0.40 + 0.001 * month_idx
    cube[bands.index("B03"), 1, 0] = 0.45 + 0.001 * month_idx

    # Median summaries differ strongly across the two pixels.
    cube[bands.index("B08"), 0, 0] += 0.01 * month_idx
    cube[bands.index("B11"), 0, 0] += 0.01 * month_idx
    cube[bands.index("B12"), 0, 0] += 0.01 * month_idx

    cube[bands.index("B08"), 1, 0] += 0.10 + 0.01 * month_idx
    cube[bands.index("B11"), 1, 0] += 0.10 + 0.01 * month_idx
    cube[bands.index("B12"), 1, 0] += 0.10 + 0.01 * month_idx

    reflectance = xr.DataArray(
        cube,
        dims=["band", "y", "x"],
        coords={"band": bands, "y": [0, 1], "x": [0]},
    )
    quality = xr.DataArray(
        np.zeros((2, 1), dtype=np.int16),
        dims=["y", "x"],
        coords={"y": [0, 1], "x": [0]},
    )
    sample_index = xr.DataArray(
        np.full((2, 1), month_idx, dtype=np.int16),
        dims=["y", "x"],
        coords={"y": [0, 1], "x": [0]},
    )
    return MonthlyBestPixelComposite(
        reflectance=reflectance,
        quality=quality,
        sample_index=sample_index,
        year=2020 + month_idx // 3,
        month=(month_idx % 3) + 6,
    )


def test_build_monthly_composite_database_requires_15_composites() -> None:
    composites = [_make_composite(i) for i in range(14)]
    try:
        build_monthly_composite_database(
            composites,
            query_bands=("B08", "B11", "B12"),
            visible_bands=("B02", "B03"),
        )
    except ValueError as exc:
        assert "15" in str(exc)
    else:
        raise AssertionError("Expected ValueError when composite count != 15")


def test_monthly_composite_database_uses_median_summary_in_query() -> None:
    composites = [_make_composite(i) for i in range(15)]
    database = build_monthly_composite_database(
        composites,
        query_bands=("B08", "B11", "B12"),
        visible_bands=("B02", "B03"),
    )

    assert isinstance(database, MonthlyCompositeDatabase)
    assert database.feature_names == (
        "nir",
        "swir1",
        "swir2",
        "median_nir",
        "median_swir1",
        "median_swir2",
    )

    corrected = xr.Dataset(
        {
            "B08": xr.DataArray(np.full((2, 1), 0.47, dtype=np.float32), dims=["y", "x"]),
            "B11": xr.DataArray(np.full((2, 1), 0.37, dtype=np.float32), dims=["y", "x"]),
            "B12": xr.DataArray(np.full((2, 1), 0.27, dtype=np.float32), dims=["y", "x"]),
        }
    )

    visible, visible_unc = database.predict_visible(corrected, k_neighbors=1)

    np.testing.assert_allclose(
        visible.sel(band="B02").values[:, 0],
        np.array([0.107, 0.400], dtype=np.float32),
        atol=1e-6,
    )
    np.testing.assert_allclose(
        visible.sel(band="B03").values[:, 0],
        np.array([0.157, 0.450], dtype=np.float32),
        atol=1e-6,
    )
    assert visible_unc.shape == visible.shape
    assert np.all(visible_unc.values >= 0.0)
