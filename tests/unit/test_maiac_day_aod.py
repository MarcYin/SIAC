"""Tests for the per-day MAIAC AOD gate source (bestpixel surface prior).

earthaccess (network + HDF I/O) is mocked: a fake source supplies granule lists
and local paths, and the granule-read / tile-merge helpers are monkeypatched so
the per-day area-median reduction is exercised without touching the filesystem
or the network.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
import xarray as xr

import siac.adapters.atmo.maiac_day_aod as maiac_mod
from siac.adapters.atmo.maiac_day_aod import MAIACDayAODProvider, _month_temporal_range

_CRS = "EPSG:32633"
_BOUNDS = (300000.0, 8199400.0, 300600.0, 8200000.0)


class _FakeCatalog:
    def resolve_short_name(self, key: str) -> str:
        assert key == "mcd19_aod"
        return "MCD19A2"


class _FakeSource:
    def __init__(
        self,
        granules: list[Any],
        paths: list[Path],
        *,
        raise_search: bool = False,
    ) -> None:
        self._granules = granules
        self._paths = paths
        self._raise_search = raise_search
        self.temporal_calls: list[Any] = []

    def search_granules(self, *, temporal: Any = None, **_kwargs: Any) -> list[Any]:
        self.temporal_calls.append(temporal)
        if self._raise_search:
            raise RuntimeError("cmr unavailable")
        return self._granules

    def download_granules(self, granules: list[Any], dest: Any) -> list[Path]:
        return self._paths


def _provider(source: _FakeSource) -> MAIACDayAODProvider:
    return MAIACDayAODProvider(source=source, catalog=_FakeCatalog())


def _patch_tile_helpers(monkeypatch: pytest.MonkeyPatch, day_value: dict[str, float]) -> None:
    """Patch the HDF-read / grid / merge helpers with in-memory fakes.

    Each granule's AOD tile is a constant keyed off the granule's day so the
    per-day area-median is deterministic.
    """

    def _fake_read(path: str | Path, dataset: str) -> tuple[np.ndarray, dict[str, Any]]:
        assert dataset == "Optical_Depth_055"
        day = maiac_mod.parse_granule_date(path).date().isoformat()
        return np.full((2, 2), day_value[day], dtype=np.float32), {}

    def _fake_native(values: np.ndarray, *, granule_path: str | Path) -> xr.DataArray:
        return xr.DataArray(np.asarray(values), dims=("y", "x"))

    def _fake_merge(arrays: list[xr.DataArray], **_kwargs: Any) -> xr.DataArray:
        stacked = xr.concat([a.expand_dims("t") for a in arrays], dim="t")
        return stacked.mean("t", skipna=True)

    monkeypatch.setattr(maiac_mod, "read_hdf4_dataset", _fake_read)
    monkeypatch.setattr(maiac_mod, "make_native_grid_dataarray", _fake_native)
    monkeypatch.setattr(maiac_mod, "merge_reprojected_tiles", _fake_merge)


def test_month_temporal_range_spans_full_month() -> None:
    assert _month_temporal_range(2021, 2) == (
        "2021-02-01T00:00:00Z",
        "2021-02-28T23:59:59Z",
    )
    # Leap February.
    assert _month_temporal_range(2024, 2)[1] == "2024-02-29T23:59:59Z"


def test_day_aod_map_reduces_per_day(monkeypatch: pytest.MonkeyPatch) -> None:
    paths = [
        Path("MCD19A2.A2021060.h18v04.061.hdf"),  # 2021-03-01
        Path("MCD19A2.A2021060.h19v04.061.hdf"),  # 2021-03-01 (second tile)
        Path("MCD19A2.A2021061.h18v04.061.hdf"),  # 2021-03-02
    ]
    source = _FakeSource(granules=["g1", "g2", "g3"], paths=paths)
    _patch_tile_helpers(monkeypatch, {"2021-03-01": 0.3, "2021-03-02": 0.5})

    result = _provider(source).day_aod_map(_BOUNDS, _CRS, [(2021, 3)])

    assert result == {"2021-03-01": pytest.approx(0.3), "2021-03-02": pytest.approx(0.5)}
    # The whole month was searched.
    assert source.temporal_calls == [("2021-03-01T00:00:00Z", "2021-03-31T23:59:59Z")]


def test_day_aod_map_omits_window_without_granules(monkeypatch: pytest.MonkeyPatch) -> None:
    source = _FakeSource(granules=[], paths=[])
    _patch_tile_helpers(monkeypatch, {})
    assert _provider(source).day_aod_map(_BOUNDS, _CRS, [(2021, 3)]) == {}


def test_day_aod_map_tolerates_search_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    source = _FakeSource(granules=["g"], paths=[], raise_search=True)
    _patch_tile_helpers(monkeypatch, {})
    # A failed window is logged and skipped, never raised.
    assert _provider(source).day_aod_map(_BOUNDS, _CRS, [(2021, 3)]) == {}


def test_day_aod_map_skips_unparseable_granule(monkeypatch: pytest.MonkeyPatch) -> None:
    paths = [
        Path("not-a-modis-granule.hdf"),  # parse_granule_date raises -> skipped
        Path("MCD19A2.A2021060.h18v04.061.hdf"),  # 2021-03-01
    ]
    source = _FakeSource(granules=["g1", "g2"], paths=paths)
    _patch_tile_helpers(monkeypatch, {"2021-03-01": 0.2})
    result = _provider(source).day_aod_map(_BOUNDS, _CRS, [(2021, 3)])
    assert result == {"2021-03-01": pytest.approx(0.2)}


def test_day_aod_map_drops_all_nan_day(monkeypatch: pytest.MonkeyPatch) -> None:
    paths = [Path("MCD19A2.A2021060.h18v04.061.hdf")]  # 2021-03-01
    source = _FakeSource(granules=["g1"], paths=paths)

    def _fake_read(path: str | Path, dataset: str) -> tuple[np.ndarray, dict[str, Any]]:
        return np.full((2, 2), np.nan, dtype=np.float32), {}

    monkeypatch.setattr(maiac_mod, "read_hdf4_dataset", _fake_read)
    def _fake_native(values: np.ndarray, *, granule_path: str | Path) -> xr.DataArray:
        return xr.DataArray(np.asarray(values), dims=("y", "x"))

    monkeypatch.setattr(maiac_mod, "make_native_grid_dataarray", _fake_native)
    monkeypatch.setattr(
        maiac_mod,
        "merge_reprojected_tiles",
        lambda arrays, **_k: arrays[0],
    )
    # Day is all-NaN -> no AOD value emitted.
    assert source and _provider(source).day_aod_map(_BOUNDS, _CRS, [(2021, 3)]) == {}
