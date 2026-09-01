"""Tests for global candidate-pool generation."""

from __future__ import annotations

import numpy as np
import pytest
from tools.aeronet_validation.global_candidate_pool import (
    EXCLUDED_CLASSES,
    candidates_from_class_grid,
    continent_of,
    global_tiles,
    pool_composition,
)


def _coords(n: int, start: float = 0.0, step: float = 0.5):
    return [start + step * i for i in range(n)]


def test_continent_boxes_cover_the_major_landmasses() -> None:
    assert continent_of(20.0, 5.0) == "africa"
    assert continent_of(100.0, 30.0) == "asia"
    assert continent_of(-100.0, 40.0) == "north_america"
    assert continent_of(-60.0, -20.0) == "south_america"
    assert continent_of(140.0, -25.0) == "oceania"
    assert continent_of(-30.0, -60.0) == "other"


def test_global_tiles_skip_the_polar_caps() -> None:
    tiles = global_tiles(30.0)
    assert min(t[1] for t in tiles) == -60.0
    assert max(t[3] for t in tiles) == 82.0
    assert all(t[0] >= -180.0 and t[2] <= 180.0 for t in tiles)


def test_draw_is_balanced_across_classes_not_uniform_over_pixels() -> None:
    # A tile that is 98% tree cover with a sliver of bare ground: a uniform
    # pixel draw would return almost no bare_sparse, which is the class the
    # catalogue most needs.
    codes = np.full((10, 10), 1, dtype=int)  # IGBP 1 = evergreen needleleaf
    codes[0, :2] = 16  # IGBP 16 = barren
    got = candidates_from_class_grid(
        codes, longitudes=_coords(10), latitudes=_coords(10), tile_id="t", per_tile=4, seed=1
    )
    composition = pool_composition(got)
    assert composition.get("bare_sparse", 0) == 2
    assert composition.get("tree_cover", 0) == 2


def test_water_and_permanent_ice_are_excluded() -> None:
    codes = np.array([[17, 15], [1, 16]])  # water, ice, forest, barren
    got = candidates_from_class_grid(
        codes, longitudes=_coords(2), latitudes=_coords(2), tile_id="t", per_tile=10
    )
    assert {c.land_cover for c in got} == {"tree_cover", "bare_sparse"}
    assert not ({c.land_cover for c in got} & EXCLUDED_CLASSES)


def test_unknown_codes_are_ignored() -> None:
    codes = np.array([[255, 0], [1, 1]])
    got = candidates_from_class_grid(
        codes, longitudes=_coords(2), latitudes=_coords(2), tile_id="t", per_tile=10
    )
    assert len(got) == 2


def test_candidates_carry_coordinates_and_continent() -> None:
    codes = np.array([[16]])
    (candidate,) = candidates_from_class_grid(
        codes, longitudes=[20.0], latitudes=[5.0], tile_id="t", per_tile=1
    )
    assert candidate.longitude == 20.0
    assert candidate.latitude == 5.0
    assert candidate.continent == "africa"
    assert candidate.land_cover == "bare_sparse"


def test_all_water_tile_yields_nothing_rather_than_failing() -> None:
    codes = np.full((4, 4), 17, dtype=int)
    assert (
        candidates_from_class_grid(
            codes, longitudes=_coords(4), latitudes=_coords(4), tile_id="t", per_tile=5
        )
        == ()
    )


def test_grid_shape_must_match_the_coordinates() -> None:
    with pytest.raises(ValueError, match="does not match"):
        candidates_from_class_grid(
            np.zeros((3, 4)), longitudes=_coords(3), latitudes=_coords(3), tile_id="t", per_tile=1
        )


def test_draw_is_deterministic_for_a_seed() -> None:
    codes = np.random.default_rng(0).choice([1, 6, 10, 16], size=(20, 20))
    kwargs = {"longitudes": _coords(20), "latitudes": _coords(20), "tile_id": "t", "per_tile": 12}
    a = candidates_from_class_grid(codes, seed=5, **kwargs)
    b = candidates_from_class_grid(codes, seed=5, **kwargs)
    assert [c.sample_id for c in a] == [c.sample_id for c in b]
