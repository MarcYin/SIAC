"""Quality-aware spatial smoothing helpers for retrieved AOD fields."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numpy import typing as npt
from scipy.ndimage import (
    binary_dilation,
    convolve,
    distance_transform_edt,
    gaussian_filter,
    median_filter,
)
from scipy.spatial import cKDTree

from siac._rust_compat import whittaker_smooth_cube

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

BoolArray = npt.NDArray[np.bool_]
FloatArray = npt.NDArray[np.float32]

SEED_EXCLUSION_QA_NAMES = (
    "invalid_retrieval",
    "low_quality",
    "parameter_boundary",
    "zero_obs_support",
    "aot_lower_boundary",
    "aot_upper_boundary",
    "tcwv_lower_boundary",
    "tcwv_upper_boundary",
)


@dataclass(frozen=True)
class HoldoutMetrics:
    """Cross-validation metrics computed on held-out trusted pixels."""

    valid_count: int
    mae: float
    rmse: float
    bias: float
    corr: float

    def to_dict(self) -> dict[str, float | int]:
        return {
            "valid_count": self.valid_count,
            "mae": self.mae,
            "rmse": self.rmse,
            "bias": self.bias,
            "corr": self.corr,
        }


def _as_bool_mask(mask: np.ndarray, *, shape: tuple[int, int]) -> BoolArray:
    array = np.asarray(mask, dtype=bool)
    if array.shape != shape:
        raise ValueError(f"Expected mask shape {shape}, got {array.shape}")
    return array.astype(bool, copy=False)


def build_trusted_aod_seed_mask(
    aot: np.ndarray,
    *,
    cloud_mask: np.ndarray,
    qa_masks: Mapping[str, np.ndarray] | None = None,
    spectral_mapping_mask: np.ndarray | None = None,
    border_pixels: int = 1,
    cloud_buffer_pixels: int = 1,
    sharp_transition_buffer_pixels: int = 1,
) -> BoolArray:
    """Return a strict mask of trusted AOD pixels suitable as smoothing seeds."""

    values = np.asarray(aot, dtype=np.float32)
    if values.ndim != 2:
        raise ValueError(f"AOT must be 2-D, got {values.shape}")

    shape = values.shape
    trusted = np.isfinite(values) & (values > 0.0)
    exclusions = np.zeros(shape, dtype=bool)

    cloud = _as_bool_mask(cloud_mask, shape=shape)
    exclusions |= cloud

    qa_lookup = (
        {}
        if qa_masks is None
        else {name: _as_bool_mask(mask, shape=shape) for name, mask in qa_masks.items()}
    )
    for name in SEED_EXCLUSION_QA_NAMES:
        exclusions |= qa_lookup.get(name, np.zeros(shape, dtype=bool))

    if spectral_mapping_mask is not None:
        exclusions |= _as_bool_mask(spectral_mapping_mask, shape=shape)

    sharp_transition = qa_lookup.get("sharp_transition_excluded", np.zeros(shape, dtype=bool))
    exclusions |= sharp_transition
    if cloud_buffer_pixels > 0:
        exclusions |= binary_dilation(cloud, iterations=int(cloud_buffer_pixels))
    if sharp_transition_buffer_pixels > 0:
        exclusions |= binary_dilation(
            sharp_transition, iterations=int(sharp_transition_buffer_pixels)
        )

    if border_pixels > 0:
        buffer = int(border_pixels)
        exclusions[:buffer, :] = True
        exclusions[-buffer:, :] = True
        exclusions[:, :buffer] = True
        exclusions[:, -buffer:] = True

    trusted &= ~exclusions
    return trusted.astype(bool, copy=False)


def distance_to_seed_pixels(seed_mask: np.ndarray) -> FloatArray:
    """Return Euclidean pixel distance to the nearest trusted seed."""

    mask = np.asarray(seed_mask, dtype=bool)
    if mask.ndim != 2:
        raise ValueError(f"Seed mask must be 2-D, got {mask.shape}")
    return np.asarray(distance_transform_edt(~mask), dtype=np.float32)


def nearest_seed_fill(values: np.ndarray, seed_mask: np.ndarray) -> FloatArray:
    """Fill every pixel with the value from the nearest trusted seed."""

    source = np.asarray(values, dtype=np.float32)
    mask = np.asarray(seed_mask, dtype=bool)
    if source.shape != mask.shape:
        raise ValueError(f"Value/seed shape mismatch: {source.shape} vs {mask.shape}")
    if not np.any(mask):
        raise ValueError("At least one trusted seed pixel is required")

    _, indices = distance_transform_edt(~mask, return_indices=True)
    return np.asarray(source[tuple(indices)], dtype=np.float32)


def normalized_gaussian_surface(
    values: np.ndarray,
    seed_mask: np.ndarray,
    *,
    sigma: float,
) -> FloatArray:
    """Interpolate a smooth field using normalized Gaussian convolution of trusted seeds."""

    if sigma <= 0.0:
        raise ValueError("sigma must be positive")

    source = np.asarray(values, dtype=np.float32)
    mask = np.asarray(seed_mask, dtype=bool)
    if source.shape != mask.shape:
        raise ValueError(f"Value/seed shape mismatch: {source.shape} vs {mask.shape}")

    weights = mask.astype(np.float32)
    numerator = gaussian_filter(np.where(mask, source, 0.0), sigma=float(sigma), mode="nearest")
    denominator = gaussian_filter(weights, sigma=float(sigma), mode="nearest")
    surface = np.divide(
        numerator,
        denominator,
        out=np.full(source.shape, np.nan, dtype=np.float32),
        where=denominator > 1.0e-6,
    )
    missing = ~np.isfinite(surface)
    if np.any(missing):
        surface[missing] = nearest_seed_fill(source, mask)[missing]
    return np.asarray(surface, dtype=np.float32)


def multiscale_normalized_gaussian_surface(
    values: np.ndarray,
    seed_mask: np.ndarray,
    *,
    sigmas: Sequence[float] = (1.5, 3.0, 6.0),
    min_support_weight: float = 1.0e-2,
) -> FloatArray:
    """Fill gaps using the smallest Gaussian scale that has enough trusted support."""

    if min_support_weight <= 0.0:
        raise ValueError("min_support_weight must be positive")

    source = np.asarray(values, dtype=np.float32)
    mask = np.asarray(seed_mask, dtype=bool)
    if source.shape != mask.shape:
        raise ValueError(f"Value/seed shape mismatch: {source.shape} vs {mask.shape}")

    sigma_values = tuple(float(sigma) for sigma in sigmas)
    if not sigma_values or any(sigma <= 0.0 for sigma in sigma_values):
        raise ValueError("sigmas must contain one or more positive values")

    weights = mask.astype(np.float32)
    filled = np.full(source.shape, np.nan, dtype=np.float32)
    unresolved = np.ones(source.shape, dtype=bool)
    last_surface: FloatArray | None = None
    fallback = nearest_seed_fill(source, mask)

    for sigma in sigma_values:
        numerator = gaussian_filter(np.where(mask, source, 0.0), sigma=sigma, mode="nearest")
        denominator = gaussian_filter(weights, sigma=sigma, mode="nearest")
        surface = np.divide(
            numerator,
            denominator,
            out=np.full(source.shape, np.nan, dtype=np.float32),
            where=denominator > 1.0e-6,
        )
        supported = unresolved & (denominator >= float(min_support_weight)) & np.isfinite(surface)
        filled[supported] = surface[supported]
        unresolved &= ~supported
        last_surface = np.asarray(surface, dtype=np.float32)
        if not np.any(unresolved):
            break

    if np.any(unresolved):
        if last_surface is not None:
            usable = unresolved & np.isfinite(last_surface)
            filled[usable] = last_surface[usable]
            unresolved &= ~usable
        filled[unresolved] = fallback[unresolved]

    return np.asarray(filled, dtype=np.float32)


def idw_knn_surface(
    values: np.ndarray,
    seed_mask: np.ndarray,
    *,
    k: int = 16,
    power: float = 2.0,
) -> FloatArray:
    """Interpolate a surface from trusted seeds using inverse-distance weighting."""

    if k <= 0:
        raise ValueError("k must be positive")
    if power <= 0.0:
        raise ValueError("power must be positive")

    source = np.asarray(values, dtype=np.float32)
    mask = np.asarray(seed_mask, dtype=bool)
    if source.shape != mask.shape:
        raise ValueError(f"Value/seed shape mismatch: {source.shape} vs {mask.shape}")

    seed_rows, seed_cols = np.nonzero(mask)
    if seed_rows.size < 1:
        raise ValueError("At least one trusted seed pixel is required")

    seed_points = np.column_stack([seed_rows, seed_cols])
    tree = cKDTree(seed_points)
    grid_rows, grid_cols = np.indices(source.shape, dtype=np.int32)
    query_points = np.column_stack([grid_rows.ravel(), grid_cols.ravel()])

    distances, indices = tree.query(query_points, k=min(int(k), seed_points.shape[0]), workers=-1)
    if distances.ndim == 1:
        distances = distances[:, np.newaxis]
        indices = indices[:, np.newaxis]

    weights = 1.0 / np.maximum(distances, 1.0e-3) ** float(power)
    exact = distances[:, 0] < 1.0e-6
    if np.any(exact):
        weights[exact] = 0.0
        weights[exact, 0] = 1.0

    seed_values = source[seed_rows, seed_cols]
    predicted = np.sum(weights * seed_values[indices], axis=1) / np.sum(weights, axis=1)
    return predicted.reshape(source.shape).astype(np.float32)


def harmonic_surface(
    values: np.ndarray,
    seed_mask: np.ndarray,
    *,
    iterations: int = 250,
) -> FloatArray:
    """Interpolate a smooth harmonic field anchored on trusted seeds."""

    if iterations <= 0:
        raise ValueError("iterations must be positive")

    source = np.asarray(values, dtype=np.float32)
    mask = np.asarray(seed_mask, dtype=bool)
    if source.shape != mask.shape:
        raise ValueError(f"Value/seed shape mismatch: {source.shape} vs {mask.shape}")
    if not np.any(mask):
        raise ValueError("At least one trusted seed pixel is required")

    filled = nearest_seed_fill(source, mask)
    kernel = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]], dtype=np.float32) / 4.0
    target = ~mask
    for _ in range(int(iterations)):
        averaged = convolve(filled, kernel, mode="nearest")
        filled[target] = averaged[target]
        filled[mask] = source[mask]
    return np.asarray(filled, dtype=np.float32)


def preserve_seed_values(
    surface: np.ndarray, original: np.ndarray, seed_mask: np.ndarray
) -> FloatArray:
    """Overwrite trusted seed locations with the original retrieved values."""

    out = np.asarray(surface, dtype=np.float32).copy()
    source = np.asarray(original, dtype=np.float32)
    mask = np.asarray(seed_mask, dtype=bool)
    if out.shape != source.shape or out.shape != mask.shape:
        raise ValueError("Surface, original, and seed mask must share the same shape")
    out[mask] = source[mask]
    return out


def _whittaker_smooth_axis(
    values: np.ndarray,
    weights: np.ndarray,
    *,
    axis: int,
    lambda_: float,
) -> FloatArray:
    source = np.asarray(values, dtype=np.float32)
    weight_array = np.asarray(weights, dtype=np.float32)
    if source.shape != weight_array.shape:
        raise ValueError(f"Value/weight shape mismatch: {source.shape} vs {weight_array.shape}")
    if lambda_ < 0.0:
        raise ValueError("lambda_ must be non-negative")

    if axis == 1:
        cube_values = np.ascontiguousarray(source.T[:, np.newaxis, :, np.newaxis], dtype=np.float32)
        cube_weights = np.ascontiguousarray(
            weight_array.T[:, np.newaxis, :, np.newaxis], dtype=np.float32
        )
        smoothed = np.asarray(
            whittaker_smooth_cube(cube_values, cube_weights, float(lambda_)), dtype=np.float32
        )
        return np.asarray(smoothed[:, 0, :, 0].T, dtype=np.float32)
    if axis == 0:
        cube_values = np.ascontiguousarray(source[:, np.newaxis, np.newaxis, :], dtype=np.float32)
        cube_weights = np.ascontiguousarray(
            weight_array[:, np.newaxis, np.newaxis, :], dtype=np.float32
        )
        smoothed = np.asarray(
            whittaker_smooth_cube(cube_values, cube_weights, float(lambda_)), dtype=np.float32
        )
        return np.asarray(smoothed[:, 0, 0, :], dtype=np.float32)
    raise ValueError(f"axis must be 0 or 1, got {axis}")


def whittaker_xy_surface(
    values: np.ndarray,
    seed_mask: np.ndarray,
    *,
    lambda_x: float = 20.0,
    lambda_y: float = 20.0,
    seed_weight: float = 50.0,
    carry_weight: float = 1.0,
    passes: int = 3,
    init_surface: np.ndarray | None = None,
) -> FloatArray:
    """Apply a weighted Whittaker smoother along x then y, preserving trusted seeds."""

    if lambda_x < 0.0 or lambda_y < 0.0:
        raise ValueError("lambda_x and lambda_y must be non-negative")
    if seed_weight <= 0.0:
        raise ValueError("seed_weight must be positive")
    if carry_weight <= 0.0:
        raise ValueError("carry_weight must be positive")
    if passes <= 0:
        raise ValueError("passes must be positive")

    source = np.asarray(values, dtype=np.float32)
    mask = np.asarray(seed_mask, dtype=bool)
    if source.shape != mask.shape:
        raise ValueError(f"Value/seed shape mismatch: {source.shape} vs {mask.shape}")
    if not np.any(mask):
        raise ValueError("At least one trusted seed pixel is required")

    if init_surface is None:
        surface = multiscale_normalized_gaussian_surface(source, mask)
    else:
        surface = np.asarray(init_surface, dtype=np.float32).copy()
        if surface.shape != source.shape:
            raise ValueError(f"Initial surface shape mismatch: {surface.shape} vs {source.shape}")
        missing = ~np.isfinite(surface)
        if np.any(missing):
            surface[missing] = nearest_seed_fill(source, mask)[missing]
    surface[mask] = source[mask]

    fallback = nearest_seed_fill(source, mask)
    weights = np.where(mask, float(seed_weight), float(carry_weight)).astype(np.float32)

    for _ in range(int(passes)):
        surface = _whittaker_smooth_axis(surface, weights, axis=1, lambda_=float(lambda_x))
        bad = ~np.isfinite(surface)
        if np.any(bad):
            surface[bad] = fallback[bad]
        surface[mask] = source[mask]

        surface = _whittaker_smooth_axis(surface, weights, axis=0, lambda_=float(lambda_y))
        bad = ~np.isfinite(surface)
        if np.any(bad):
            surface[bad] = fallback[bad]
        surface[mask] = source[mask]

    return np.asarray(surface, dtype=np.float32)


def median_from_seed_surface(
    values: np.ndarray,
    seed_mask: np.ndarray,
    *,
    window_size: tuple[int, int] = (20, 20),
    init_surface: np.ndarray | None = None,
) -> FloatArray:
    """Apply a spatial median filter to a surface derived from trusted seeds only."""

    if len(window_size) != 2:
        raise ValueError("window_size must be a 2-tuple")
    if any(int(size) <= 0 for size in window_size):
        raise ValueError("window_size entries must be positive")

    source = np.asarray(values, dtype=np.float32)
    mask = np.asarray(seed_mask, dtype=bool)
    if source.shape != mask.shape:
        raise ValueError(f"Value/seed shape mismatch: {source.shape} vs {mask.shape}")
    if not np.any(mask):
        raise ValueError("At least one trusted seed pixel is required")

    if init_surface is None:
        base = harmonic_surface(source, mask, iterations=250)
    else:
        base = np.asarray(init_surface, dtype=np.float32)
        if base.shape != source.shape:
            raise ValueError(f"Initial surface shape mismatch: {base.shape} vs {source.shape}")

    smoothed = median_filter(
        base,
        size=tuple(int(size) for size in window_size),
        mode="nearest",
    )
    return np.asarray(smoothed, dtype=np.float32)


def sample_holdout_mask(
    seed_mask: np.ndarray,
    *,
    rng_seed: int = 42,
    holdout_fraction: float = 0.1,
    min_holdout: int = 5_000,
    max_holdout: int = 20_000,
) -> BoolArray:
    """Sample a validation subset from trusted seed pixels."""

    if not 0.0 < holdout_fraction < 1.0:
        raise ValueError("holdout_fraction must be between 0 and 1")

    mask = np.asarray(seed_mask, dtype=bool)
    seed_indices = np.flatnonzero(mask)
    if seed_indices.size == 0:
        raise ValueError("Cannot sample holdout pixels from an empty seed mask")

    requested = int(round(seed_indices.size * float(holdout_fraction)))
    holdout_size = min(int(max_holdout), max(int(min_holdout), requested))
    holdout_size = min(holdout_size, seed_indices.size)

    rng = np.random.default_rng(int(rng_seed))
    chosen = rng.choice(seed_indices, size=holdout_size, replace=False)
    holdout = np.zeros(mask.shape, dtype=bool)
    holdout.flat[chosen] = True
    return holdout.astype(bool, copy=False)


def score_holdout(
    surface: np.ndarray, truth: np.ndarray, holdout_mask: np.ndarray
) -> HoldoutMetrics:
    """Score predictions on a held-out subset of trusted pixels."""

    predicted = np.asarray(surface, dtype=np.float32)
    expected = np.asarray(truth, dtype=np.float32)
    holdout = np.asarray(holdout_mask, dtype=bool)
    if predicted.shape != expected.shape or predicted.shape != holdout.shape:
        raise ValueError("Surface, truth, and holdout mask must share the same shape")

    valid = holdout & np.isfinite(predicted) & np.isfinite(expected)
    count = int(np.count_nonzero(valid))
    if count == 0:
        return HoldoutMetrics(
            valid_count=0, mae=float("nan"), rmse=float("nan"), bias=float("nan"), corr=float("nan")
        )

    pred_values = np.asarray(predicted[valid], dtype=np.float64)
    truth_values = np.asarray(expected[valid], dtype=np.float64)
    delta = pred_values - truth_values

    if (
        pred_values.size > 1
        and not np.allclose(pred_values, pred_values[0])
        and not np.allclose(truth_values, truth_values[0])
    ):
        corr = float(np.corrcoef(pred_values, truth_values)[0, 1])
    else:
        corr = float("nan")

    return HoldoutMetrics(
        valid_count=count,
        mae=float(np.mean(np.abs(delta))),
        rmse=float(np.sqrt(np.mean(delta**2))),
        bias=float(np.mean(delta)),
        corr=corr,
    )
