"""Route-B database built from 15 monthly BRDF composites."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
import logging
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr
from scipy.spatial import cKDTree

from siac.algorithms.grid.assembler import _resample_da

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.algorithms.surface.brdf_monthly_composite import MonthlyBestPixelComposite

logger = logging.getLogger(__name__)


def _feature_names_for_query_bands(query_bands: Sequence[str]) -> tuple[str, ...]:
    if len(query_bands) == 3:
        return ("nir", "swir1", "swir2", "median_nir", "median_swir1", "median_swir2")
    query_names = tuple(f"query_{idx}" for idx in range(len(query_bands)))
    median_names = tuple(f"median_query_{idx}" for idx in range(len(query_bands)))
    return query_names + median_names


def _dataset_to_cube(dataset_or_array: xr.Dataset | xr.DataArray, band_names: Sequence[str]) -> xr.DataArray:
    if isinstance(dataset_or_array, xr.DataArray):
        if "band" not in dataset_or_array.dims:
            raise ValueError("DataArray inputs must have a 'band' dimension")
        return dataset_or_array.sel(band=list(band_names))
    return xr.concat([dataset_or_array[name] for name in band_names], dim=xr.IndexVariable("band", list(band_names)))


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
    """Searchable database over 15 monthly composites."""

    entries_features: np.ndarray
    entries_visible: np.ndarray
    median_summary: xr.DataArray
    visible_band_names: tuple[str, ...]
    query_band_names: tuple[str, ...]
    feature_names: tuple[str, ...]
    y_coords: xr.DataArray
    x_coords: xr.DataArray

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
    ) -> tuple[xr.DataArray, xr.DataArray]:
        """Predict visible reflectance from corrected NIR/SWIR query bands."""
        if k_neighbors < 1:
            raise ValueError("k_neighbors must be >= 1")

        query_cube = _dataset_to_cube(corrected_reflectance, self.query_band_names).astype(np.float32)
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
        y_coords = query_cube.coords["y"] if "y" in query_cube.coords else np.arange(ny)
        x_coords = query_cube.coords["x"] if "x" in query_cube.coords else np.arange(nx)

        features_flat = np.empty((n_pixels, n_query + median_values.shape[0]), dtype=np.float32)
        features_flat[:, :n_query] = query_values.reshape(n_query, n_pixels).T
        features_flat[:, n_query:] = median_values.reshape(median_values.shape[0], n_pixels).T
        valid_query_rows = np.flatnonzero(np.all(np.isfinite(features_flat), axis=1))
        if valid_query_rows.size == 0 or self.entries_features.size == 0:
            coords = {"band": list(self.visible_band_names), "y": y_coords, "x": x_coords}
            return (
                xr.DataArray(predicted, dims=["band", "y", "x"], coords=coords),
                xr.DataArray(uncertainty, dims=["band", "y", "x"], coords=coords),
            )

        neighbor_index = self._neighbor_index
        if neighbor_index is None:
            coords = {"band": list(self.visible_band_names), "y": y_coords, "x": x_coords}
            return (
                xr.DataArray(predicted, dims=["band", "y", "x"], coords=coords),
                xr.DataArray(uncertainty, dims=["band", "y", "x"], coords=coords),
            )

        neighbor_count = min(k_neighbors, self.entries_features.shape[0])
        predicted_flat = predicted.reshape(n_visible, n_pixels).T
        uncertainty_flat = uncertainty.reshape(n_visible, n_pixels).T
        chunk_size = 4096

        for start in range(0, valid_query_rows.size, chunk_size):
            query_rows = valid_query_rows[start : start + chunk_size]
            distances, neighbor_rows = neighbor_index.query(
                features_flat[query_rows],
                k=neighbor_count,
                workers=1,
            )
            if neighbor_count == 1:
                distances = np.asarray(distances, dtype=np.float32)[:, np.newaxis]
                neighbor_rows = np.asarray(neighbor_rows, dtype=np.intp)[:, np.newaxis]
            else:
                distances = np.asarray(distances, dtype=np.float32)
                neighbor_rows = np.asarray(neighbor_rows, dtype=np.intp)

            selected_visible = self.entries_visible[neighbor_rows]
            zero_mask = distances == 0.0
            nonzero_rows = ~np.any(zero_mask, axis=1)

            if np.any(nonzero_rows):
                selected_distances = distances[nonzero_rows]
                weighted_visible = selected_visible[nonzero_rows]
                weights = 1.0 / np.maximum(selected_distances, 1e-6)
                weights = weights / np.sum(weights, axis=1, keepdims=True)
                estimate = np.sum(weighted_visible * weights[..., np.newaxis], axis=1)
                spread = np.sqrt(
                    np.sum(
                        ((weighted_visible - estimate[:, np.newaxis, :]) ** 2)
                        * weights[..., np.newaxis],
                        axis=1,
                    )
                )
                predicted_flat[query_rows[nonzero_rows]] = estimate.astype(np.float32, copy=False)
                uncertainty_flat[query_rows[nonzero_rows]] = spread.astype(np.float32, copy=False)

            if np.any(~nonzero_rows):
                zero_rows = np.flatnonzero(~nonzero_rows)
                for local_index in zero_rows:
                    matched = selected_visible[local_index][zero_mask[local_index]]
                    estimate = matched.mean(axis=0)
                    spread = (
                        matched.std(axis=0)
                        if matched.shape[0] > 1
                        else np.zeros(n_visible, dtype=np.float32)
                    )
                    flat_index = query_rows[local_index]
                    predicted_flat[flat_index] = estimate.astype(np.float32, copy=False)
                    uncertainty_flat[flat_index] = spread.astype(np.float32, copy=False)

        coords = {"band": list(self.visible_band_names), "y": y_coords, "x": x_coords}
        return (
            xr.DataArray(predicted, dims=["band", "y", "x"], coords=coords),
            xr.DataArray(uncertainty, dims=["band", "y", "x"], coords=coords),
        )


def build_monthly_composite_database(
    composites: Sequence[MonthlyBestPixelComposite],
    *,
    query_bands: Sequence[str],
    visible_bands: Sequence[str],
) -> MonthlyCompositeDatabase:
    """Build the Route-B database from exactly 15 monthly composites."""
    if len(composites) != 15:
        raise ValueError("Route-B monthly database requires exactly 15 monthly composites")
    if len(query_bands) < 1:
        raise ValueError("query_bands must not be empty")
    if len(visible_bands) < 1:
        raise ValueError("visible_bands must not be empty")

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
    entries_visible = np.empty((n_composites * n_pixels, n_visible), dtype=np.float32)
    for index, composite in enumerate(composites):
        if composite.reflectance.sizes.get("y") != ny or composite.reflectance.sizes.get("x") != nx:
            raise ValueError("All monthly composites must share the same spatial shape")
        query_values[index] = composite.reflectance.sel(band=list(query_names)).values.astype(np.float32)
        visible_values = composite.reflectance.sel(band=list(visible_names)).values.astype(np.float32)
        start = index * n_pixels
        end = start + n_pixels
        entries_visible[start:end] = visible_values.reshape(n_visible, n_pixels).T

    median_summary = np.nanmedian(query_values, axis=0)

    entries_features = np.empty((n_composites * n_pixels, n_query * 2), dtype=np.float32)
    median_flat = median_summary.reshape(n_query, n_pixels).T
    for index in range(n_composites):
        start = index * n_pixels
        end = start + n_pixels
        entries_features[start:end, :n_query] = query_values[index].reshape(n_query, n_pixels).T
        entries_features[start:end, n_query:] = median_flat

    valid_entries = np.all(np.isfinite(entries_features), axis=1) & np.all(np.isfinite(entries_visible), axis=1)
    entries_features = entries_features[valid_entries]
    entries_visible = entries_visible[valid_entries]
    logger.info(
        "Monthly composite database build complete: valid_entries=%d total_entries=%d",
        int(entries_features.shape[0]),
        int(n_composites * n_pixels),
    )

    return MonthlyCompositeDatabase(
        entries_features=entries_features,
        entries_visible=entries_visible,
        median_summary=xr.DataArray(
            median_summary,
            dims=["feature", "y", "x"],
            coords={
                "feature": list(_feature_names_for_query_bands(query_names)[len(query_names):]),
                "y": first.coords["y"],
                "x": first.coords["x"],
            },
        ),
        visible_band_names=visible_names,
        query_band_names=query_names,
        feature_names=_feature_names_for_query_bands(query_names),
        y_coords=first.coords["y"],
        x_coords=first.coords["x"],
    )
