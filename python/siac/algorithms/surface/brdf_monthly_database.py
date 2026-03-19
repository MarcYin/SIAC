"""Route-B database built from 15 monthly BRDF composites."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.algorithms.surface.brdf_monthly_composite import MonthlyBestPixelComposite


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
        if query_cube.sizes.get("y") != self.median_summary.sizes.get("y") or query_cube.sizes.get("x") != self.median_summary.sizes.get("x"):
            raise ValueError("corrected_reflectance must share the database spatial grid")

        query_values = np.asarray(query_cube.values, dtype=np.float32)
        median_values = np.asarray(self.median_summary.values, dtype=np.float32)
        query_features = np.concatenate([query_values, median_values], axis=0)

        n_visible = len(self.visible_band_names)
        ny = query_cube.sizes["y"]
        nx = query_cube.sizes["x"]
        predicted = np.full((n_visible, ny, nx), np.nan, dtype=np.float32)
        uncertainty = np.full((n_visible, ny, nx), np.nan, dtype=np.float32)

        features_flat = query_features.reshape(query_features.shape[0], -1).T
        valid_entries = np.all(np.isfinite(self.entries_features), axis=1) & np.all(np.isfinite(self.entries_visible), axis=1)
        database_features = self.entries_features[valid_entries]
        database_visible = self.entries_visible[valid_entries]

        for flat_index, query_vector in enumerate(features_flat):
            if not np.all(np.isfinite(query_vector)):
                continue
            distances = np.linalg.norm(database_features - query_vector[np.newaxis, :], axis=1)
            if distances.size == 0:
                continue
            order = np.argsort(distances)[: min(k_neighbors, distances.size)]
            selected_distances = distances[order]
            selected_visible = database_visible[order]
            if np.any(selected_distances == 0.0):
                matched = selected_visible[selected_distances == 0.0]
                estimate = matched.mean(axis=0)
                spread = matched.std(axis=0) if matched.shape[0] > 1 else np.zeros(n_visible, dtype=np.float32)
            else:
                weights = 1.0 / np.maximum(selected_distances, 1e-6)
                weights = weights / np.sum(weights)
                estimate = np.sum(selected_visible * weights[:, np.newaxis], axis=0)
                spread = np.sqrt(
                    np.sum(((selected_visible - estimate[np.newaxis, :]) ** 2) * weights[:, np.newaxis], axis=0)
                )
            row = flat_index // nx
            col = flat_index % nx
            predicted[:, row, col] = estimate.astype(np.float32)
            uncertainty[:, row, col] = spread.astype(np.float32)

        coords = {"band": list(self.visible_band_names), "y": self.y_coords, "x": self.x_coords}
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

    query_stack = []
    visible_stack = []
    for composite in composites:
        if composite.reflectance.sizes.get("y") != ny or composite.reflectance.sizes.get("x") != nx:
            raise ValueError("All monthly composites must share the same spatial shape")
        query_stack.append(
            composite.reflectance.sel(band=list(query_names)).values.astype(np.float32)
        )
        visible_stack.append(
            composite.reflectance.sel(band=list(visible_names)).values.astype(np.float32)
        )

    query_values = np.stack(query_stack, axis=0)
    visible_values = np.stack(visible_stack, axis=0)
    median_summary = np.nanmedian(query_values, axis=0)

    median_repeated = np.broadcast_to(median_summary[np.newaxis, ...], query_values.shape)
    feature_stack = np.concatenate([query_values, median_repeated], axis=1)

    entries_features = feature_stack.transpose(0, 2, 3, 1).reshape(-1, feature_stack.shape[1])
    entries_visible = visible_values.transpose(0, 2, 3, 1).reshape(-1, visible_values.shape[1])

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
