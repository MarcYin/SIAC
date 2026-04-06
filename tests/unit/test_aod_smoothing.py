from __future__ import annotations

import numpy as np

import siac.algorithms.solver.aod_smoothing as aod_smoothing_mod
from siac.algorithms.solver.aod_smoothing import (
    build_trusted_aod_seed_mask,
    harmonic_surface,
    median_from_seed_surface,
    multiscale_normalized_gaussian_surface,
    nearest_seed_fill,
    normalized_gaussian_surface,
    preserve_seed_values,
    sample_holdout_mask,
    score_holdout,
    whittaker_xy_surface,
)


def test_build_trusted_aod_seed_mask_excludes_flags_buffers_and_border() -> None:
    aot = np.ones((5, 5), dtype=np.float32)
    aot[3, 3] = 0.0
    cloud_mask = np.zeros((5, 5), dtype=bool)
    cloud_mask[2, 2] = True
    qa_masks = {
        "low_quality": np.zeros((5, 5), dtype=bool),
        "sharp_transition_excluded": np.zeros((5, 5), dtype=bool),
    }
    qa_masks["low_quality"][1, 3] = True
    qa_masks["sharp_transition_excluded"][3, 1] = True

    trusted = build_trusted_aod_seed_mask(
        aot,
        cloud_mask=cloud_mask,
        qa_masks=qa_masks,
        border_pixels=1,
        cloud_buffer_pixels=1,
        sharp_transition_buffer_pixels=1,
    )

    expected = np.zeros((5, 5), dtype=bool)
    expected[1, 1] = True
    np.testing.assert_array_equal(trusted, expected)


def test_build_trusted_aod_seed_mask_excludes_sharp_transition_without_buffer() -> None:
    aot = np.ones((3, 3), dtype=np.float32)
    cloud_mask = np.zeros((3, 3), dtype=bool)
    qa_masks = {"sharp_transition_excluded": np.zeros((3, 3), dtype=bool)}
    qa_masks["sharp_transition_excluded"][1, 1] = True

    trusted = build_trusted_aod_seed_mask(
        aot,
        cloud_mask=cloud_mask,
        qa_masks=qa_masks,
        border_pixels=0,
        cloud_buffer_pixels=0,
        sharp_transition_buffer_pixels=0,
    )

    expected = np.ones((3, 3), dtype=bool)
    expected[1, 1] = False
    np.testing.assert_array_equal(trusted, expected)


def test_smoothers_fill_gaps_and_holdout_scoring_is_finite() -> None:
    aot = np.array(
        [
            [0.10, 0.12, 0.14],
            [0.11, 0.40, 0.15],
            [0.12, 0.13, 0.16],
        ],
        dtype=np.float32,
    )
    seed_mask = np.array(
        [
            [True, True, False],
            [True, False, False],
            [False, False, True],
        ],
        dtype=bool,
    )

    gaussian = normalized_gaussian_surface(aot, seed_mask, sigma=1.2)
    harmonic = harmonic_surface(aot, seed_mask, iterations=25)
    nearest = nearest_seed_fill(aot, seed_mask)
    preserved = preserve_seed_values(gaussian, aot, seed_mask)

    assert gaussian.shape == aot.shape
    assert harmonic.shape == aot.shape
    assert np.all(np.isfinite(gaussian))
    assert np.all(np.isfinite(harmonic))
    assert np.all(np.isfinite(nearest))
    np.testing.assert_array_equal(preserved[seed_mask], aot[seed_mask])

    holdout = sample_holdout_mask(seed_mask, rng_seed=3, holdout_fraction=0.25, min_holdout=1, max_holdout=2)
    metrics = score_holdout(harmonic, aot, holdout)
    assert metrics.valid_count >= 1
    assert np.isfinite(metrics.mae)
    assert np.isfinite(metrics.rmse)


def test_multiscale_gaussian_surface_fills_using_progressively_broader_support() -> None:
    aot = np.full((5, 5), 0.2, dtype=np.float32)
    aot[0, 0] = 0.1
    aot[4, 4] = 0.3
    seed_mask = np.zeros((5, 5), dtype=bool)
    seed_mask[0, 0] = True
    seed_mask[4, 4] = True

    surface = multiscale_normalized_gaussian_surface(
        aot,
        seed_mask,
        sigmas=(0.75, 1.5, 3.0),
        min_support_weight=5.0e-2,
    )

    assert np.all(np.isfinite(surface))
    assert float(surface[2, 2]) > 0.1
    assert float(surface[2, 2]) < 0.3


def test_whittaker_xy_surface_preserves_trusted_seeds(monkeypatch) -> None:
    def _fake_whittaker(cube: np.ndarray, weights: np.ndarray, lambda_: float) -> np.ndarray:
        assert cube.shape == weights.shape
        out = np.asarray(cube, dtype=np.float32).copy()
        positive = weights > 0
        mean_value = float(np.mean(cube[positive])) if np.any(positive) else 0.0
        out[~positive] = mean_value
        return out

    monkeypatch.setattr(aod_smoothing_mod, "whittaker_smooth_cube", _fake_whittaker)

    aot = np.array(
        [
            [0.10, 0.12, 0.14],
            [0.11, 0.40, 0.15],
            [0.12, 0.13, 0.16],
        ],
        dtype=np.float32,
    )
    seed_mask = np.array(
        [
            [True, False, False],
            [False, False, False],
            [False, False, True],
        ],
        dtype=bool,
    )

    surface = whittaker_xy_surface(
        aot,
        seed_mask,
        lambda_x=10.0,
        lambda_y=10.0,
        seed_weight=10.0,
        carry_weight=1.0,
        passes=2,
    )

    assert np.all(np.isfinite(surface))
    np.testing.assert_array_equal(surface[seed_mask], aot[seed_mask])


def test_median_from_seed_surface_returns_finite_field() -> None:
    aot = np.array(
        [
            [0.10, 0.12, 0.14, 0.15],
            [0.11, 0.40, 0.16, 0.17],
            [0.12, 0.13, 0.18, 0.19],
            [0.13, 0.14, 0.20, 0.21],
        ],
        dtype=np.float32,
    )
    seed_mask = np.array(
        [
            [True, True, False, False],
            [True, False, False, False],
            [False, False, False, True],
            [False, False, True, True],
        ],
        dtype=bool,
    )

    surface = median_from_seed_surface(aot, seed_mask, window_size=(2, 2))

    assert surface.shape == aot.shape
    assert np.all(np.isfinite(surface))
