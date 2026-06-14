"""Route-B database built from monthly BRDF composites."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING, cast

import numpy as np
import xarray as xr
from scipy.spatial import cKDTree

from siac.algorithms.grid.assembler import _resample_da
from siac.runtime.models import copy_spatial_metadata_like

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.algorithms.surface.brdf_monthly_composite import MonthlyBestPixelComposite

logger = logging.getLogger(__name__)

# Systematic surface-prior uncertainty floor for the Route-B monthly database.
# The kNN neighbour spread captures only local sampling variability and omits the
# systematic MODIS->S2 mapping error, BRDF absolute uncertainty, and scene-vs-
# composite temporal mismatch. Left raw (~0.003-0.005) the solver over-trusts the
# prior and AOT overshoots (broad AERONET sweep: 40% within EE, +0.09 bias). The
# spread is combined in quadrature with this floor so the total lands near the
# kernel/whittaker ~0.007 sweet spot.
_MONTHLY_UNCERTAINTY_FLOOR = 0.006


@dataclass(frozen=True)
class _MonthlyPredictionDiagnostics:
    predicted: xr.DataArray
    uncertainty: xr.DataArray
    quality: xr.DataArray
    source_fit_rmse: xr.DataArray
    knn_feature_distance: xr.DataArray


def _query_feature_names(query_bands: Sequence[str]) -> tuple[str, ...]:
    """Names for the per-pixel current-observation query-band features."""
    if len(query_bands) == 3:
        return ("nir", "swir1", "swir2")
    return tuple(f"query_{idx}" for idx in range(len(query_bands)))


def _median_feature_names(
    query_bands: Sequence[str],
    visible_bands: Sequence[str],
    *,
    include_visible: bool,
) -> tuple[str, ...]:
    """Names for the temporal-median ("climatology") key features.

    Always includes the NIR/SWIR query medians; when ``include_visible`` is set
    (the ``monthly_database_median_key="all"`` option) it also appends the
    per-pixel temporal median of the visible bands as a full-spectrum location
    fingerprint.
    """
    query_median: tuple[str, ...]
    if len(query_bands) == 3:
        query_median = ("median_nir", "median_swir1", "median_swir2")
    else:
        query_median = tuple(f"median_query_{idx}" for idx in range(len(query_bands)))
    if not include_visible:
        return query_median
    return query_median + tuple(f"median_{name}" for name in visible_bands)


def _dataset_to_cube(
    dataset_or_array: xr.Dataset | xr.DataArray, band_names: Sequence[str]
) -> xr.DataArray:
    if isinstance(dataset_or_array, xr.DataArray):
        if "band" not in dataset_or_array.dims:
            raise ValueError("DataArray inputs must have a 'band' dimension")
        return dataset_or_array.sel(band=list(band_names))
    return xr.concat(
        [dataset_or_array[name] for name in band_names],
        dim=xr.IndexVariable("band", list(band_names)),
    )


def _resample_summary_to_query_grid(
    summary: xr.DataArray,
    query_cube: xr.DataArray,
) -> xr.DataArray:
    target_shape = (int(query_cube.sizes["y"]), int(query_cube.sizes["x"]))
    if summary.sizes.get("y") == target_shape[0] and summary.sizes.get("x") == target_shape[1]:
        return summary

    feature_coords = (
        summary.coords["feature"].values
        if "feature" in summary.coords
        else np.arange(summary.sizes["feature"])
    )
    resampled = xr.concat(
        [
            _resample_da(summary.sel(feature=feature, drop=True), target_shape, "bilinear")
            for feature in feature_coords
        ],
        dim=xr.IndexVariable("feature", feature_coords),
    )
    coords: dict[str, object] = {"feature": feature_coords}
    if "y" in query_cube.coords:
        coords["y"] = query_cube.coords["y"]
    if "x" in query_cube.coords:
        coords["x"] = query_cube.coords["x"]
    return resampled.assign_coords(**coords)


@dataclass(frozen=True)
class MonthlyCompositeDatabase:
    """Searchable database over monthly composites."""

    entries_features: np.ndarray
    entries_visible: np.ndarray
    entries_quality: np.ndarray
    entries_source_fit_rmse: np.ndarray
    median_summary: xr.DataArray
    visible_band_names: tuple[str, ...]
    query_band_names: tuple[str, ...]
    feature_names: tuple[str, ...]
    y_coords: xr.DataArray
    x_coords: xr.DataArray
    composites: tuple[MonthlyBestPixelComposite, ...] = ()

    @cached_property
    def _neighbor_index(self) -> cKDTree | None:
        if self.entries_features.size == 0:
            return None
        return cKDTree(self.entries_features, copy_data=False)

    def predict_visible(
        self,
        corrected_reflectance: xr.Dataset | xr.DataArray,
        *,
        k_neighbors: int = 3,
    ) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
        """Predict visible reflectance from corrected NIR/SWIR query bands."""
        diagnostics = self.predict_visible_with_diagnostics(
            corrected_reflectance,
            k_neighbors=k_neighbors,
        )
        return diagnostics.predicted, diagnostics.uncertainty, diagnostics.quality

    def predict_visible_with_diagnostics(
        self,
        corrected_reflectance: xr.Dataset | xr.DataArray,
        *,
        k_neighbors: int = 3,
    ) -> _MonthlyPredictionDiagnostics:
        """Predict visible reflectance plus query-quality diagnostics."""
        if k_neighbors < 1:
            raise ValueError("k_neighbors must be >= 1")

        query_cube = _dataset_to_cube(corrected_reflectance, self.query_band_names).astype(
            np.float32
        )
        query_summary = _resample_summary_to_query_grid(self.median_summary, query_cube)

        query_values = np.asarray(query_cube.values, dtype=np.float32)
        median_values = np.asarray(query_summary.values, dtype=np.float32)

        n_visible = len(self.visible_band_names)
        n_query = len(self.query_band_names)
        ny = query_cube.sizes["y"]
        nx = query_cube.sizes["x"]
        n_pixels = ny * nx
        predicted = np.full((n_visible, ny, nx), np.nan, dtype=np.float32)
        uncertainty = np.full((n_visible, ny, nx), np.nan, dtype=np.float32)
        quality = np.full((ny, nx), np.nan, dtype=np.float32)
        source_fit_rmse = np.full((ny, nx), np.nan, dtype=np.float32)
        knn_feature_distance = np.full((ny, nx), np.nan, dtype=np.float32)
        y_coords = query_cube.coords["y"] if "y" in query_cube.coords else np.arange(ny)
        x_coords = query_cube.coords["x"] if "x" in query_cube.coords else np.arange(nx)

        def _output_arrays() -> _MonthlyPredictionDiagnostics:
            coords = {"band": list(self.visible_band_names), "y": y_coords, "x": x_coords}
            predicted_da = xr.DataArray(predicted, dims=["band", "y", "x"], coords=coords)
            # Combine the neighbour spread with the systematic floor in quadrature;
            # NaN (no-prediction) pixels stay NaN. See _MONTHLY_UNCERTAINTY_FLOOR.
            floored_uncertainty = np.sqrt(
                np.asarray(uncertainty, dtype=np.float32) ** 2 + _MONTHLY_UNCERTAINTY_FLOOR**2
            ).astype(np.float32)
            uncertainty_da = xr.DataArray(floored_uncertainty, dims=["band", "y", "x"], coords=coords)
            quality_da = xr.DataArray(
                quality, dims=["y", "x"], coords={"y": y_coords, "x": x_coords}
            )
            source_fit_da = xr.DataArray(
                source_fit_rmse,
                dims=["y", "x"],
                coords={"y": y_coords, "x": x_coords},
            )
            distance_da = xr.DataArray(
                knn_feature_distance,
                dims=["y", "x"],
                coords={"y": y_coords, "x": x_coords},
            )
            return _MonthlyPredictionDiagnostics(
                copy_spatial_metadata_like(predicted_da, query_cube),
                copy_spatial_metadata_like(uncertainty_da, query_cube),
                copy_spatial_metadata_like(
                    quality_da, cast("xr.DataArray", query_cube.isel(band=0, drop=True))
                ),
                copy_spatial_metadata_like(
                    source_fit_da, cast("xr.DataArray", query_cube.isel(band=0, drop=True))
                ),
                copy_spatial_metadata_like(
                    distance_da, cast("xr.DataArray", query_cube.isel(band=0, drop=True))
                ),
            )

        features_flat = np.empty((n_pixels, n_query + median_values.shape[0]), dtype=np.float32)
        features_flat[:, :n_query] = query_values.reshape(n_query, n_pixels).T
        features_flat[:, n_query:] = median_values.reshape(median_values.shape[0], n_pixels).T
        valid_query_rows = np.flatnonzero(np.all(np.isfinite(features_flat), axis=1))
        if valid_query_rows.size == 0 or self.entries_features.size == 0:
            return _output_arrays()

        neighbor_index = self._neighbor_index
        if neighbor_index is None:
            return _output_arrays()

        neighbor_count = min(k_neighbors, self.entries_features.shape[0])
        predicted_flat = predicted.reshape(n_visible, n_pixels).T
        uncertainty_flat = uncertainty.reshape(n_visible, n_pixels).T
        quality_flat = quality.reshape(n_pixels)
        source_fit_flat = source_fit_rmse.reshape(n_pixels)
        distance_flat = knn_feature_distance.reshape(n_pixels)

        # Single batch KNN query with multi-threaded tree traversal.
        distances, neighbor_rows = neighbor_index.query(
            features_flat[valid_query_rows],
            k=neighbor_count,
            workers=-1,
        )
        if neighbor_count == 1:
            distances = np.asarray(distances, dtype=np.float32)[:, np.newaxis]
            neighbor_rows = np.asarray(neighbor_rows, dtype=np.intp)[:, np.newaxis]
        else:
            distances = np.asarray(distances, dtype=np.float32)
            neighbor_rows = np.asarray(neighbor_rows, dtype=np.intp)

        # Gather all neighbor data in bulk.
        selected_visible = self.entries_visible[neighbor_rows]
        selected_quality = self.entries_quality[neighbor_rows]
        selected_source_fit = self.entries_source_fit_rmse[neighbor_rows]
        selected_query_features = self.entries_features[neighbor_rows, :n_query]

        query_feature_values = features_flat[valid_query_rows, :n_query]
        feature_distance = np.sqrt(
            np.mean(
                np.square(
                    selected_query_features - query_feature_values[:, np.newaxis, :],
                    dtype=np.float32,
                ),
                axis=2,
                dtype=np.float32,
            )
        ).astype(np.float32, copy=False)

        # Vectorised inverse-distance weighting for all non-exact-match rows.
        zero_mask = distances == 0.0
        nonzero_rows = ~np.any(zero_mask, axis=1)

        if np.any(nonzero_rows):
            nz_idx = np.flatnonzero(nonzero_rows)
            nz_distances = distances[nz_idx]
            nz_visible = selected_visible[nz_idx]
            weights = 1.0 / np.maximum(nz_distances, 1e-6)
            weights = weights / np.sum(weights, axis=1, keepdims=True)
            estimate = np.sum(nz_visible * weights[..., np.newaxis], axis=1)
            spread = np.sqrt(
                np.sum(
                    ((nz_visible - estimate[:, np.newaxis, :]) ** 2) * weights[..., np.newaxis],
                    axis=1,
                )
            )
            out_idx = valid_query_rows[nz_idx]
            predicted_flat[out_idx] = estimate.astype(np.float32, copy=False)
            uncertainty_flat[out_idx] = spread.astype(np.float32, copy=False)
            quality_flat[out_idx] = np.sum(selected_quality[nz_idx] * weights, axis=1)
            source_fit_flat[out_idx] = np.sum(selected_source_fit[nz_idx] * weights, axis=1)
            distance_flat[out_idx] = np.sum(feature_distance[nz_idx] * weights, axis=1)

        # Handle exact-match rows (rare but possible).
        if np.any(~nonzero_rows):
            zero_rows = np.flatnonzero(~nonzero_rows)
            for local_index in zero_rows:
                matched = selected_visible[local_index][zero_mask[local_index]]
                matched_quality = selected_quality[local_index][zero_mask[local_index]]
                matched_source_fit = selected_source_fit[local_index][zero_mask[local_index]]
                estimate = matched.mean(axis=0)
                spread = (
                    matched.std(axis=0)
                    if matched.shape[0] > 1
                    else np.zeros(n_visible, dtype=np.float32)
                )
                flat_index = valid_query_rows[local_index]
                predicted_flat[flat_index] = estimate.astype(np.float32, copy=False)
                uncertainty_flat[flat_index] = spread.astype(np.float32, copy=False)
                quality_flat[flat_index] = np.float32(matched_quality.mean())
                source_fit_flat[flat_index] = np.float32(matched_source_fit.mean())
                distance_flat[flat_index] = 0.0

        return _output_arrays()


def build_monthly_composite_database(
    composites: Sequence[MonthlyBestPixelComposite],
    *,
    query_bands: Sequence[str],
    visible_bands: Sequence[str],
    max_source_fit_rmse: float | None = None,
    median_key: str = "query",
) -> MonthlyCompositeDatabase:
    """Build the Route-B database from one or more monthly composites.

    ``median_key`` controls the kNN lookup key's "climatology" block:
    ``"query"`` (default) appends only the temporal median of the NIR/SWIR query
    bands; ``"all"`` additionally appends the temporal median of the visible
    bands as a full-spectrum per-pixel fingerprint.
    """
    if len(composites) < 1:
        raise ValueError("Route-B monthly database requires at least one monthly composite")
    if len(query_bands) < 1:
        raise ValueError("query_bands must not be empty")
    if len(visible_bands) < 1:
        raise ValueError("visible_bands must not be empty")
    if median_key not in ("query", "all"):
        raise ValueError(f"median_key must be 'query' or 'all', got {median_key!r}")
    include_visible_median = median_key == "all"

    query_names = tuple(query_bands)
    visible_names = tuple(visible_bands)
    first = composites[0].reflectance
    ny = first.sizes["y"]
    nx = first.sizes["x"]
    n_composites = len(composites)
    n_query = len(query_names)
    n_visible = len(visible_names)
    n_pixels = ny * nx
    logger.info(
        "Monthly composite database build start: composites=%d query_bands=%d visible_bands=%d grid=%dx%d",
        n_composites,
        n_query,
        n_visible,
        int(ny),
        int(nx),
    )

    query_values = np.empty((n_composites, n_query, ny, nx), dtype=np.float32)
    # Per-composite visible stack is only needed when the visible temporal
    # median feeds the lookup key (median_key="all").
    visible_stack = (
        np.empty((n_composites, n_visible, ny, nx), dtype=np.float32)
        if include_visible_median
        else None
    )
    entries_visible = np.empty((n_composites * n_pixels, n_visible), dtype=np.float32)
    entries_quality = np.empty(n_composites * n_pixels, dtype=np.float32)
    entries_source_fit_rmse = np.zeros(n_composites * n_pixels, dtype=np.float32)
    for index, composite in enumerate(composites):
        if composite.reflectance.sizes.get("y") != ny or composite.reflectance.sizes.get("x") != nx:
            raise ValueError("All monthly composites must share the same spatial shape")
        query_values[index] = composite.reflectance.sel(band=list(query_names)).values.astype(
            np.float32
        )
        visible_values = composite.reflectance.sel(band=list(visible_names)).values.astype(
            np.float32
        )
        if visible_stack is not None:
            visible_stack[index] = visible_values
        quality_values = composite.quality.values.astype(np.float32)
        source_fit_rmse = getattr(composite, "source_fit_rmse", None)
        source_fit_values = (
            source_fit_rmse.values.astype(np.float32)
            if source_fit_rmse is not None
            else np.zeros((ny, nx), dtype=np.float32)
        )
        start = index * n_pixels
        end = start + n_pixels
        entries_visible[start:end] = visible_values.reshape(n_visible, n_pixels).T
        entries_quality[start:end] = quality_values.reshape(n_pixels)
        entries_source_fit_rmse[start:end] = source_fit_values.reshape(n_pixels)

    # Temporal-median ("climatology") key block: always the NIR/SWIR query
    # medians, optionally extended with the visible-band medians.
    query_median = np.nanmedian(query_values, axis=0)  # (n_query, ny, nx)
    if include_visible_median and visible_stack is not None:
        visible_median = np.nanmedian(visible_stack, axis=0)  # (n_visible, ny, nx)
        median_summary = np.concatenate([query_median, visible_median], axis=0)
    else:
        median_summary = query_median
    median_block = int(median_summary.shape[0])

    entries_features = np.empty((n_composites * n_pixels, n_query + median_block), dtype=np.float32)
    median_flat = median_summary.reshape(median_block, n_pixels).T
    for index in range(n_composites):
        start = index * n_pixels
        end = start + n_pixels
        entries_features[start:end, :n_query] = query_values[index].reshape(n_query, n_pixels).T
        entries_features[start:end, n_query:] = median_flat

    valid_entries = (
        np.all(np.isfinite(entries_features), axis=1)
        & np.all(np.isfinite(entries_visible), axis=1)
        & np.isfinite(entries_quality)
        & np.isfinite(entries_source_fit_rmse)
    )
    if max_source_fit_rmse is not None:
        valid_entries &= entries_source_fit_rmse <= float(max_source_fit_rmse)
    entries_features = entries_features[valid_entries]
    entries_visible = entries_visible[valid_entries]
    entries_quality = entries_quality[valid_entries]
    entries_source_fit_rmse = entries_source_fit_rmse[valid_entries]
    logger.info(
        "Monthly composite database build complete: valid_entries=%d total_entries=%d",
        int(entries_features.shape[0]),
        int(n_composites * n_pixels),
    )

    median_feature_names = _median_feature_names(
        query_names, visible_names, include_visible=include_visible_median
    )
    feature_names = _query_feature_names(query_names) + median_feature_names
    median_summary_da = copy_spatial_metadata_like(
        xr.DataArray(
            median_summary,
            dims=["feature", "y", "x"],
            coords={
                "feature": list(median_feature_names),
                "y": first.coords["y"],
                "x": first.coords["x"],
            },
        ),
        cast("xr.DataArray", first.isel(band=0, drop=True)),
    )

    return MonthlyCompositeDatabase(
        entries_features=entries_features,
        entries_visible=entries_visible,
        entries_quality=entries_quality,
        entries_source_fit_rmse=entries_source_fit_rmse,
        median_summary=median_summary_da,
        visible_band_names=visible_names,
        query_band_names=query_names,
        feature_names=feature_names,
        y_coords=first.coords["y"],
        x_coords=first.coords["x"],
        composites=tuple(composites),
    )
