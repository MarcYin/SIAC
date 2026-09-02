"""Tests for deriving an AOI grid from a catalogue point."""

from __future__ import annotations

import pytest
from tools.aeronet_validation.build_global_corpus_index import (
    GRID_RESOLUTION_M,
    TEMPLATE_RESOLUTION_M,
    TEMPLATE_SIZE,
    crs_for_mgrs_tile,
    grid_for_point,
    scene_day_of,
)
from tools.aeronet_validation.build_l2a_scl_index import _imports


def test_crs_follows_from_the_tile_id() -> None:
    # Saves opening a band per AOI purely to read a CRS the name already fixes.
    assert crs_for_mgrs_tile("T31UDQ") == "EPSG:32631"
    assert crs_for_mgrs_tile("31UDQ") == "EPSG:32631"
    assert crs_for_mgrs_tile("T47NQD") == "EPSG:32647"


def test_southern_hemisphere_tiles_get_a_southern_crs() -> None:
    # Band N is the first northern band, so M and below are southern.
    assert crs_for_mgrs_tile("T20HLK") == "EPSG:32720"
    assert crs_for_mgrs_tile("T21LXK") == "EPSG:32721"
    assert crs_for_mgrs_tile("T47NQD").startswith("EPSG:326")


def test_malformed_tile_ids_are_rejected_not_guessed() -> None:
    for bad in ("", "TXXUDQ", "T99UDQ", "T31IDQ", "T31"):
        with pytest.raises(ValueError):
            crs_for_mgrs_tile(bad)


def test_grid_snaps_to_the_60m_template_so_the_scales_nest() -> None:
    # context60 depends on exact 3:1 nesting; snapping to 20 m would let the
    # coarse grid straddle the fine one.
    imports = _imports()
    _, transform, span = grid_for_point(12.0, 45.7, "T32TQR", imports=imports)
    assert span == TEMPLATE_SIZE * TEMPLATE_RESOLUTION_M / GRID_RESOLUTION_M == 384
    assert transform.c % TEMPLATE_RESOLUTION_M == pytest.approx(0.0)
    assert transform.f % TEMPLATE_RESOLUTION_M == pytest.approx(0.0)
    assert transform.a == GRID_RESOLUTION_M
    assert transform.e == -GRID_RESOLUTION_M


def test_grid_is_centred_on_the_catalogue_point() -> None:
    imports = _imports()
    crs, transform, span = grid_for_point(12.0, 45.7, "T32TQR", imports=imports)
    to_utm = imports["Transformer"].from_crs("EPSG:4326", crs, always_xy=True)
    x, y = to_utm.transform(12.0, 45.7)
    centre_x = transform.c + transform.a * span / 2
    centre_y = transform.f + transform.e * span / 2
    # Within one template pixel of the point, which is all the snap allows.
    assert abs(centre_x - x) <= TEMPLATE_RESOLUTION_M
    assert abs(centre_y - y) <= TEMPLATE_RESOLUTION_M


def test_scene_day_is_taken_from_the_iso_timestamp() -> None:
    assert scene_day_of("2024-05-24T14:47:04.685000+00:00") == "2024-05-24"


def test_maiac_windows_are_chunked_to_survive_the_memory_limit() -> None:
    # A whole +/-45 day window in one getInfo has already produced "User memory
    # limit exceeded" on this project, so windows are split before being sent.
    from tools.aeronet_validation.maiac_gee_day_aod import CHUNK_DAYS, _chunks

    chunks = list(_chunks("2020-11-01", "2021-01-29"))
    assert len(chunks) > 1
    assert chunks[0][0] == "2020-11-01"
    assert chunks[-1][1] == "2021-01-29"
    # Contiguous and non-overlapping: no day counted twice or dropped.
    import datetime as dt

    for earlier, later in zip(chunks, chunks[1:]):
        assert dt.date.fromisoformat(later[0]) == dt.date.fromisoformat(earlier[1]) + dt.timedelta(
            days=1
        )
    for start, end in chunks:
        assert (dt.date.fromisoformat(end) - dt.date.fromisoformat(start)).days < CHUNK_DAYS


def test_a_single_day_window_is_one_chunk() -> None:
    from tools.aeronet_validation.maiac_gee_day_aod import _chunks

    assert list(_chunks("2020-11-01", "2020-11-01")) == [("2020-11-01", "2020-11-01")]
