"""Tests for turning catalogue points into current-date acquisitions."""

from __future__ import annotations

import datetime as dt

import pytest
from tools.aeronet_validation.build_global_scene_catalog import (
    SUPPORTED_SPACECRAFT,
    aoi_bbox,
    month_window,
    options_from_images,
    year_order,
)
from tools.aeronet_validation.global_scene_selection import ELIGIBLE_YEARS


class _Image:
    def __init__(self, **properties) -> None:
        self.properties = properties
        self.acquisition_time_utc = dt.datetime(2023, 7, 2, 10, 6, tzinfo=dt.timezone.utc)


def _image(**overrides) -> _Image:
    base = {
        "PRODUCT_ID": "S2A_MSIL1C_20230702T100601_N0509_R022_T32TQR_20230702T135151",
        "MGRS_TILE": "32TQR",
        "CLOUDY_PIXEL_PERCENTAGE": 39.3,
        "SPACECRAFT_NAME": "Sentinel-2A",
    }
    base.update(overrides)
    return _Image(**base)


def test_aoi_bbox_widens_in_longitude_towards_the_poles() -> None:
    # A fixed ground half-width spans more degrees of longitude at high
    # latitude; a constant degree box would shrink the AOI away from the equator.
    equator = aoi_bbox(0.0, 0.0)
    polar = aoi_bbox(0.0, 70.0)
    assert (polar[2] - polar[0]) > 2.5 * (equator[2] - equator[0])
    assert (polar[3] - polar[1]) == pytest.approx(equator[3] - equator[1])


def test_the_mgrs_tile_is_normalised_to_the_t_prefix() -> None:
    assert options_from_images([_image()])[0].mgrs_tile == "T32TQR"


def test_images_without_an_l1c_product_id_are_not_candidates() -> None:
    # PRODUCT_ID is the whole reason selection runs against Earth Engine rather
    # than a STAC catalogue, so an image lacking it is unusable downstream.
    assert options_from_images([_image(PRODUCT_ID=None)]) == ()
    assert options_from_images([_image(CLOUDY_PIXEL_PERCENTAGE=None)]) == ()


def test_unsupported_spacecraft_are_rejected() -> None:
    assert "Sentinel-2C" in SUPPORTED_SPACECRAFT
    assert options_from_images([_image(SPACECRAFT_NAME="Sentinel-2Z")]) == ()


def test_year_order_is_a_stable_per_aoi_shuffle() -> None:
    # Stable so a rerun selects the same scene; shuffled so scene years do not
    # pile up on the first eligible year wherever it already offers a clear one.
    first = year_order("sample_a", ELIGIBLE_YEARS)
    assert first == year_order("sample_a", ELIGIBLE_YEARS)
    assert sorted(first) == sorted(ELIGIBLE_YEARS)
    assert any(
        year_order(f"sample_{index}", ELIGIBLE_YEARS) != first for index in range(20)
    ), "every AOI drew the same order"


def test_month_window_straddles_the_month_midpoint() -> None:
    start, end = month_window(2021, 2)
    assert start == "2021-01-26"
    assert end == "2021-03-07"


def test_projected_coordinates_are_refused_before_reaching_earth_engine() -> None:
    """Regression: a metres-not-degrees catalogue must fail fast.

    Earth Engine accepts the resulting geometry and spends minutes on it before
    answering "User memory limit exceeded", which reads as a quota problem
    rather than as bad input.
    """

    with pytest.raises(ValueError, match="not a WGS84 coordinate"):
        aoi_bbox(-5262074.17, -5529868.92)
