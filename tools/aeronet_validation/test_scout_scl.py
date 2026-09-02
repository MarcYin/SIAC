"""Tests for the free SCL survey that feeds the Cloud Score+ prefilter."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from tools.aeronet_validation.build_l2a_scl_index import _imports
from tools.aeronet_validation.scout_scl import (
    SCOUT_CLEAR_CLASSES,
    ScoutGrid,
    choose_overview,
    probe_overview_level,
    read_scl,
    scout,
    scout_grid,
    summarize_scl,
)


def _grid(
    width: int = 600, height: int = 600, scale: float = 20.0, resolution_m: float = 20.0
) -> ScoutGrid:
    imports = _imports()
    return scout_grid(
        "EPSG:32631",
        imports["Affine"](scale, 0.0, 400000.0, 0.0, -scale, 5000000.0),
        height,
        width,
        resolution_m=resolution_m,
        imports=imports,
    )


def test_choose_overview_takes_the_coarsest_at_or_below_target() -> None:
    # 20 m native, 240 m target -> decimation limit 12, so 8 wins and 16 does not.
    assert choose_overview([2, 4, 8, 16], target_resolution_m=240.0) == 2


def test_choose_overview_never_returns_coarser_than_the_target() -> None:
    # Warping up from coarser data would alias a categorical field, and the
    # remaining bytes are negligible, so the full band is read instead.
    assert choose_overview([16, 32], target_resolution_m=240.0) is None
    assert choose_overview([], target_resolution_m=240.0) is None


def test_scout_grid_preserves_extent_and_rounds_pixels_up() -> None:
    grid = _grid(width=601, height=600, resolution_m=240.0)
    # 601 * 20 m = 12020 m, which is 50.08 scout pixels: rounding down would
    # leave a strip of the AOI unmeasured.
    assert (grid.width, grid.height) == (51, 50)
    assert grid.bounds[0] == pytest.approx(400000.0)
    assert grid.bounds[3] == pytest.approx(5000000.0)
    ring = grid.wgs84_polygon["coordinates"][0]
    assert ring[0] == ring[-1]


def test_clear_fraction_is_against_the_aoi_not_the_observed_part() -> None:
    # Half nodata, half vegetation: perfectly clear where seen, but only half
    # the AOI can contribute a composite pixel.
    values = np.array([[0, 0], [4, 4]], dtype=np.uint8)
    clear, coverage = summarize_scl(values)
    assert clear == pytest.approx(0.5)
    assert coverage == pytest.approx(0.5)


def test_snow_counts_as_clear_and_cloud_does_not() -> None:
    assert 11 in SCOUT_CLEAR_CLASSES
    assert {3, 8, 9, 10}.isdisjoint(SCOUT_CLEAR_CLASSES)
    values = np.array([[11, 11], [9, 3]], dtype=np.uint8)
    clear, coverage = summarize_scl(values)
    assert clear == pytest.approx(0.5)
    assert coverage == pytest.approx(1.0)


def _write_scl(path, array: np.ndarray, *, scale: float = 20.0, overviews=(2, 4, 8, 16)) -> None:
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.transform import from_origin

    height, width = array.shape
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="uint8",
        crs="EPSG:32631",
        transform=from_origin(400000.0, 5000000.0, scale, scale),
        tiled=True,
        blockxsize=256,
        blockysize=256,
    ) as handle:
        handle.write(array, 1)
        handle.build_overviews(list(overviews), Resampling.nearest)


def _footprint(grid: ScoutGrid, *, x_share: float = 1.0, pad: float = 0.05) -> dict[str, Any]:
    """A WGS84 footprint covering ``x_share`` of the AOI from its western edge.

    Derived from the grid rather than hard-coded so the overlap the test claims
    is the overlap the gate actually measures.
    """

    from pyproj import Transformer

    west, south, east, north = grid.bounds
    span = (east - west) * float(x_share)
    to_wgs84 = Transformer.from_crs(grid.crs, "EPSG:4326", always_xy=True).transform
    corners = [
        (west - pad * (east - west), south - pad * (north - south)),
        (west + span, south - pad * (north - south)),
        (west + span, north + pad * (north - south)),
        (west - pad * (east - west), north + pad * (north - south)),
    ]
    ring = [list(to_wgs84(*point)) for point in corners]
    ring.append(ring[0])
    return {"type": "Polygon", "coordinates": [ring]}


def _item(href: str, geometry: dict[str, Any] | None = None) -> dict[str, Any]:
    return {
        "id": "S2A_MSIL2A_20200615T103021_N0500_R108_T31UDQ_20230101T000000",
        "assets": {"SCL": {"href": href}},
        "geometry": geometry if geometry is not None else _footprint(_grid()),
        "properties": {"s2:mgrs_tile": "31UDQ", "datetime": "2020-06-15T10:30:21Z"},
    }


def test_overview_read_reproduces_the_native_class_fractions(tmp_path) -> None:
    # The efficiency claim: an overview is a nearest subsample of a categorical
    # field, so fractions estimated from it match the native read.
    rng = np.random.default_rng(20260901)
    field = np.where(
        rng.random((1200, 1200)) < 0.3,
        rng.choice([3, 8, 9, 10], size=(1200, 1200)),
        rng.choice([4, 5, 6], size=(1200, 1200)),
    ).astype(np.uint8)
    # Give it spatial structure, or subsampling is trivially unbiased.
    field = np.repeat(np.repeat(field[::8, ::8], 8, axis=0), 8, axis=1).astype(np.uint8)
    path = tmp_path / "scl.tif"
    _write_scl(path, field)

    imports = _imports()
    grid = _grid(width=1200, height=1200)
    item = _item(str(path), _footprint(grid))
    level = probe_overview_level(item, sas="")
    assert level == 2  # decimation 8 -> 160 m, the coarsest at or below 240 m
    decimated = read_scl(item, sas="", grid=grid, imports=imports, overview_level=level)
    native = read_scl(item, sas="", grid=grid, imports=imports, overview_level=None)
    assert summarize_scl(decimated)[0] == pytest.approx(summarize_scl(native)[0], abs=0.02)


def test_footprint_gate_skips_the_cog_read_entirely() -> None:
    # A swath-edge acquisition is rejected from geometry alone, so it costs no
    # HTTP read at all -- the point of gating before the raster.
    imports = _imports()
    grid = _grid()
    sliver = _footprint(grid, x_share=0.4)
    records, diagnostics = scout(
        [_item("/nonexistent/should-not-be-opened.tif", sliver)],
        sas="",
        grid=grid,
        imports=imports,
        workers=1,
    )
    assert records == ()
    assert diagnostics["footprint_rejected"] == 1
    assert diagnostics["scl_reads"] == 0


def test_one_unreadable_cog_does_not_end_the_survey(tmp_path) -> None:
    imports = _imports()
    grid = _grid()
    good = tmp_path / "good.tif"
    _write_scl(good, np.full((600, 600), 4, dtype=np.uint8))
    items = [_item(str(good)), dict(_item("/nonexistent/missing.tif"), id="broken")]
    records, diagnostics = scout(items, sas="", grid=grid, imports=imports, workers=2)
    assert diagnostics["scl_reads"] == 2  # both passed the footprint gate
    assert len(records) == 1
    assert diagnostics["scl_read_failures"] == 1
    assert records[0].clear_fraction == pytest.approx(1.0)


def test_grid_comes_from_the_trailing_axes_and_an_epsg_code(tmp_path) -> None:
    # Teacher archives store comp as (realization, band, y, x) and declare
    # ``epsg`` rather than ``crs``; taking shape[:2] would size the grid from
    # the realization and band axes.
    from tools.aeronet_validation.scout_scl import grid_from_archive

    path = tmp_path / "teacher.npz"
    np.savez(
        path,
        comp=np.zeros((14, 7, 384, 384), dtype=np.float32),
        transform=np.array([20.0, 0.0, 787740.0, 0.0, -20.0, 315540.0]),
        epsg=np.asarray(32647),
    )
    grid = grid_from_archive(path)
    assert grid.crs == "EPSG:32647"
    assert (grid.width, grid.height) == (384, 384)


def test_the_survey_reads_native_by_default() -> None:
    """Pins a measured decision, not a preference.

    Decimating to 240 m saved nothing at 8 concurrent workers on real Planetary
    Computer reads and reported the AOI 0.016 more clear on average, which moved
    the shortlist. Both defaults must stay native for that result to hold.
    """

    import inspect

    from tools.aeronet_validation.scout_scl import SCL_NATIVE_RESOLUTION_M, SCOUT_RESOLUTION_M

    assert SCOUT_RESOLUTION_M == SCL_NATIVE_RESOLUTION_M
    assert inspect.signature(scout).parameters["overview_level"].default is None


def test_sas_token_renews_before_it_expires() -> None:
    """Regression: a token fetched once lapses mid-run and takes the rest with it.

    Measured as all 25 shards of a 4.5 hour array ceasing to produce output at
    the same wall-clock minute -- the token's expiry -- having looked healthy
    until that moment.
    """

    import datetime as dt

    from tools.aeronet_validation.scout_scl import SasToken

    issued = []

    class _Response:
        def __init__(self, expiry_seconds: float) -> None:
            self._expiry = expiry_seconds

        def raise_for_status(self) -> None:
            return None

        def json(self) -> dict[str, str]:
            issued.append(len(issued))
            moment = dt.datetime.now(dt.timezone.utc) + dt.timedelta(seconds=self._expiry)
            return {
                "token": f"token{len(issued)}",
                "msft:expiry": moment.strftime("%Y-%m-%dT%H:%M:%SZ"),
            }

    class _Session:
        def __init__(self, expiry_seconds: float) -> None:
            self._expiry = expiry_seconds

        def get(self, url, timeout=None):  # noqa: ARG002 - signature mirrors requests
            return _Response(self._expiry)

    # Comfortably valid: fetched once, then reused.
    fresh = SasToken(_Session(3600.0), margin_seconds=300.0)
    assert fresh.value == "token1"
    assert fresh.value == "token1"
    assert fresh.renewals == 1

    # Inside the safety margin: renewed on every use rather than served stale.
    expiring = SasToken(_Session(60.0), margin_seconds=300.0)
    first = expiring.value
    second = expiring.value
    assert first != second
    assert expiring.renewals == 2


def test_sas_token_without_an_advertised_expiry_still_renews() -> None:
    from tools.aeronet_validation.scout_scl import SasToken

    class _Session:
        def get(self, url, timeout=None):  # noqa: ARG002 - signature mirrors requests
            class _Response:
                def raise_for_status(self) -> None:
                    return None

                def json(self) -> dict[str, str]:
                    return {"token": "bare"}

            return _Response()

    token = SasToken(_Session())
    assert token.value == "bare"
    # Assumed lifetime, not treated as never-expiring.
    assert token._expires_at > 0.0  # noqa: SLF001 - asserting the fallback was applied
