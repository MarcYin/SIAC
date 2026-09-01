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


def _grid_resolver(start: float = 0.0, step: float = 0.5):
    """Separable resolver standing in for an equirectangular grid."""

    def resolve(rows, cols):
        import numpy as np

        return (
            start + step * np.asarray(cols, dtype=float),
            start + step * np.asarray(rows, dtype=float),
        )

    return resolve


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
        codes, to_lonlat=_grid_resolver(), tile_id="t", per_tile=4, seed=1
    )
    composition = pool_composition(got)
    assert composition.get("bare_sparse", 0) == 2
    assert composition.get("tree_cover", 0) == 2


def test_water_and_permanent_ice_are_excluded() -> None:
    codes = np.array([[17, 15], [1, 16]])  # water, ice, forest, barren
    got = candidates_from_class_grid(codes, to_lonlat=_grid_resolver(), tile_id="t", per_tile=10)
    assert {c.land_cover for c in got} == {"tree_cover", "bare_sparse"}
    assert not ({c.land_cover for c in got} & EXCLUDED_CLASSES)


def test_unknown_codes_are_ignored() -> None:
    codes = np.array([[255, 0], [1, 1]])
    got = candidates_from_class_grid(codes, to_lonlat=_grid_resolver(), tile_id="t", per_tile=10)
    assert len(got) == 2


def test_candidates_carry_coordinates_and_continent() -> None:
    codes = np.array([[16]])
    (candidate,) = candidates_from_class_grid(
        codes, to_lonlat=lambda _rows, _cols: ([20.0], [5.0]), tile_id="t", per_tile=1
    )
    assert candidate.longitude == 20.0
    assert candidate.latitude == 5.0
    assert candidate.continent == "africa"
    assert candidate.land_cover == "bare_sparse"


def test_all_water_tile_yields_nothing_rather_than_failing() -> None:
    codes = np.full((4, 4), 17, dtype=int)
    assert (
        candidates_from_class_grid(codes, to_lonlat=_grid_resolver(), tile_id="t", per_tile=5) == ()
    )


def test_grid_must_be_two_dimensional() -> None:
    with pytest.raises(ValueError, match="must be 2-D"):
        candidates_from_class_grid(np.zeros(4), to_lonlat=_grid_resolver(), tile_id="t", per_tile=1)


def test_draw_is_deterministic_for_a_seed() -> None:
    codes = np.random.default_rng(0).choice([1, 6, 10, 16], size=(20, 20))
    kwargs = {"to_lonlat": _grid_resolver(), "tile_id": "t", "per_tile": 12}
    a = candidates_from_class_grid(codes, seed=5, **kwargs)
    b = candidates_from_class_grid(codes, seed=5, **kwargs)
    assert [c.sample_id for c in a] == [c.sample_id for c in b]


def test_sinusoidal_rows_do_not_share_a_longitude_axis() -> None:
    """Regression: MCD12Q1 is sinusoidal, where longitude depends on the row.

    The first catalogue treated the raster's projected metres as degrees, which
    produced 5,001 points at longitudes near -5,000,000 -- every one outside the
    valid range, every latitude band collapsing to "polar", and nothing in the
    row count to show for it.
    """

    import numpy as np
    from pyproj import Transformer

    sinusoidal = (
        'PROJCS["unknown",GEOGCS["unknown",DATUM["unknown",SPHEROID["unknown",6371007.181,0]],'
        'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Sinusoidal"],'
        'PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],'
        'PARAMETER["false_northing",0],UNIT["metre",1]]'
    )
    to_wgs84 = Transformer.from_crs(sinusoidal, "EPSG:4326", always_xy=True)
    scale, origin_x, origin_y = 463.3127165279165, -11238113.252095148, -5352188.501330493

    def resolve(rows, cols):
        x = origin_x + scale * (np.asarray(cols, dtype=float) + 0.5)
        y = origin_y - scale * (np.asarray(rows, dtype=float) + 0.5)
        return to_wgs84.transform(x, y)

    # Same easting, 4,000 rows apart: separable axes would give one longitude.
    (west_lon, _), (east_lon, _) = (
        tuple(v[0] for v in resolve(np.array([0]), np.array([0]))),
        tuple(v[0] for v in resolve(np.array([4000]), np.array([0]))),
    )
    assert abs(east_lon - west_lon) > 100.0

    codes = np.full((4200, 200), 16, dtype=int)  # barren everywhere
    got = candidates_from_class_grid(
        codes, to_lonlat=resolve, tile_id="t", per_tile=50, within=(-152.0, -68.0, -120.0, -48.0)
    )
    assert got, "clipping to the nominal tile removed every candidate"
    assert all(-180.0 <= c.longitude <= 180.0 for c in got)
    assert all(-90.0 <= c.latitude <= 90.0 for c in got)


def test_points_outside_the_nominal_tile_are_clipped() -> None:
    import numpy as np

    codes = np.full((4, 4), 16, dtype=int)

    # Resolver puts every point at lon 100, outside the requested box.
    def resolve(rows, _cols):
        return np.full(len(rows), 100.0), np.full(len(rows), 10.0)

    assert (
        candidates_from_class_grid(
            codes, to_lonlat=resolve, tile_id="t", per_tile=4, within=(-10.0, 0.0, 10.0, 20.0)
        )
        == ()
    )
