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


def test_build_monthly_composite_database_requires_non_empty_composites() -> None:
    composites: list[MonthlyBestPixelComposite] = []
    try:
        build_monthly_composite_database(
            composites,
            query_bands=("B08", "B11", "B12"),
            visible_bands=("B02", "B03"),
        )
    except ValueError as exc:
        assert "at least one" in str(exc)
    else:
        raise AssertionError("Expected ValueError when composite count is empty")


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

    visible, visible_unc, visible_quality = database.predict_visible(corrected, k_neighbors=1)

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
    assert visible_quality.shape == (2, 1)
    assert np.all(visible_quality.values >= 0.0)


def test_monthly_composite_database_all_median_key_appends_visible_fingerprint() -> None:
    composites = [_make_composite(i) for i in range(15)]
    default_db = build_monthly_composite_database(
        composites,
        query_bands=("B08", "B11", "B12"),
        visible_bands=("B02", "B03"),
    )
    all_db = build_monthly_composite_database(
        composites,
        query_bands=("B08", "B11", "B12"),
        visible_bands=("B02", "B03"),
        median_key="all",
    )

    # Default key: [query(3) | median_query(3)] = 6 features.
    assert default_db.entries_features.shape[1] == 6
    assert default_db.median_summary.sizes["feature"] == 3

    # "all" key appends the per-pixel visible temporal medians as a fingerprint:
    # [query(3) | median_query(3) | median_visible(2)] = 8 features.
    assert all_db.feature_names == (
        "nir",
        "swir1",
        "swir2",
        "median_nir",
        "median_swir1",
        "median_swir2",
        "median_B02",
        "median_B03",
    )
    assert all_db.entries_features.shape[1] == 8
    assert all_db.median_summary.sizes["feature"] == 5

    corrected = xr.Dataset(
        {
            "B08": xr.DataArray(np.full((2, 1), 0.47, dtype=np.float32), dims=["y", "x"]),
            "B11": xr.DataArray(np.full((2, 1), 0.37, dtype=np.float32), dims=["y", "x"]),
            "B12": xr.DataArray(np.full((2, 1), 0.27, dtype=np.float32), dims=["y", "x"]),
        }
    )
    visible, visible_unc, _ = all_db.predict_visible(corrected, k_neighbors=1)
    b02 = visible.sel(band="B02").values[:, 0]
    # The fingerprint preserves (and sharpens) per-pixel separation: the dark
    # pixel stays dark, the bright pixel stays bright.
    assert b02[0] < 0.20 < b02[1]
    assert np.all(np.isfinite(visible.values))
    assert visible_unc.shape == visible.shape


def test_monthly_composite_database_accepts_query_on_coarser_grid() -> None:
    composites = [_make_composite(i) for i in range(15)]
    database = build_monthly_composite_database(
        composites,
        query_bands=("B08", "B11", "B12"),
        visible_bands=("B02", "B03"),
    )

    corrected = xr.Dataset(
        {
            "B08": xr.DataArray(
                np.array([[0.47]], dtype=np.float32), dims=["y", "x"], coords={"y": [10], "x": [20]}
            ),
            "B11": xr.DataArray(
                np.array([[0.37]], dtype=np.float32), dims=["y", "x"], coords={"y": [10], "x": [20]}
            ),
            "B12": xr.DataArray(
                np.array([[0.27]], dtype=np.float32), dims=["y", "x"], coords={"y": [10], "x": [20]}
            ),
        }
    )

    visible, visible_unc, visible_quality = database.predict_visible(corrected, k_neighbors=1)

    assert visible.shape == (2, 1, 1)
    assert visible_unc.shape == visible.shape
    assert visible_quality.shape == (1, 1)
    assert list(visible.coords["y"].values) == [10]
    assert list(visible.coords["x"].values) == [20]
    assert np.isfinite(visible.values).all()


def test_build_monthly_composite_database_filters_high_source_fit_entries() -> None:
    composites = [_make_composite(i) for i in range(2)]
    source_fit_low = xr.DataArray(
        np.array([[0.01], [0.02]], dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [0, 1], "x": [0]},
    )
    source_fit_high = xr.DataArray(
        np.array([[0.01], [0.20]], dtype=np.float32),
        dims=["y", "x"],
        coords={"y": [0, 1], "x": [0]},
    )
    composites[0] = MonthlyBestPixelComposite(
        reflectance=composites[0].reflectance,
        quality=composites[0].quality,
        sample_index=composites[0].sample_index,
        year=composites[0].year,
        month=composites[0].month,
        source_fit_rmse=source_fit_low,
    )
    composites[1] = MonthlyBestPixelComposite(
        reflectance=composites[1].reflectance,
        quality=composites[1].quality,
        sample_index=composites[1].sample_index,
        year=composites[1].year,
        month=composites[1].month,
        source_fit_rmse=source_fit_high,
    )

    database = build_monthly_composite_database(
        composites,
        query_bands=("B08", "B11", "B12"),
        visible_bands=("B02", "B03"),
        max_source_fit_rmse=0.05,
    )

    assert database.entries_features.shape[0] == 3
    assert np.all(database.entries_source_fit_rmse <= 0.05)


def test_predict_visible_with_diagnostics_reports_query_feature_distance() -> None:
    composites = [_make_composite(i) for i in range(15)]
    database = build_monthly_composite_database(
        composites,
        query_bands=("B08", "B11", "B12"),
        visible_bands=("B02", "B03"),
    )

    corrected = xr.Dataset(
        {
            "B08": xr.DataArray(np.array([[0.47], [0.80]], dtype=np.float32), dims=["y", "x"]),
            "B11": xr.DataArray(np.array([[0.37], [0.70]], dtype=np.float32), dims=["y", "x"]),
            "B12": xr.DataArray(np.array([[0.27], [0.60]], dtype=np.float32), dims=["y", "x"]),
        }
    )

    prediction = database.predict_visible_with_diagnostics(corrected, k_neighbors=1)

    assert prediction.knn_feature_distance.shape == (2, 1)
    assert float(prediction.knn_feature_distance.values[0, 0]) == 0.0
    assert float(prediction.knn_feature_distance.values[1, 0]) > 0.0
