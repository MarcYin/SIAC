"""Unit tests for EarthAccessSource wrapper behavior."""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

import siac.adapters.earthdata as earthdata_mod
from siac.adapters.data.earthaccess_catalog import EarthAccessCatalog
from siac.adapters.data.earthaccess_source import EarthAccessSource


class _FakeEarthAccessModule:
    def __init__(self):
        self.login_calls = 0
        self.search_datasets_calls: list[dict] = []
        self.search_data_calls: list[dict] = []
        self.open_calls: list[list] = []
        self.download_calls: list[tuple[list, str, int | None]] = []

    def login(self, **kwargs):
        self.login_calls += 1
        self.last_login_kwargs = kwargs

    def search_datasets(self, **kwargs):
        self.search_datasets_calls.append(kwargs)
        if kwargs.get("short_name") == "MCD43A1":
            out = [{"short_name": "MCD43A1"}]
        else:
            out = [{"umm": {"ShortName": "MCD19A2"}}]
        count = kwargs.get("count")
        if count == 0:
            return []
        if isinstance(count, int) and count > 0:
            return out[:count]
        return out

    def search_data(self, **kwargs):
        self.search_data_calls.append(kwargs)
        out = [{"id": "g1"}, {"id": "g2"}]
        count = kwargs.get("count")
        if count == 0:
            return []
        if isinstance(count, int) and count > 0:
            return out[:count]
        return out

    def open(self, granules):
        self.open_calls.append(list(granules))
        return [f"opened:{len(granules)}"]

    def download(self, granules, dest_dir, threads=None, **_kwargs):
        self.download_calls.append((list(granules), dest_dir, threads))
        return [str(Path(dest_dir) / "g1.hdf")]


@pytest.fixture
def fake_earthaccess(monkeypatch):
    fake = _FakeEarthAccessModule()
    monkeypatch.setitem(sys.modules, "earthaccess", fake)
    return fake


class TestEarthAccessSource:
    def test_search_does_not_force_authentication(self, fake_earthaccess):
        src = EarthAccessSource(provider="LPDAAC_ECS")
        assert src.is_authenticated is False
        assert fake_earthaccess.login_calls == 0

        _ = src.search_datasets(keyword="MCD43")
        assert src.is_authenticated is False
        assert fake_earthaccess.login_calls == 0

        _ = src.search_datasets(keyword="MCD43")
        assert fake_earthaccess.login_calls == 0

    def test_search_datasets_forwards_query(self, fake_earthaccess):
        src = EarthAccessSource(provider="LPDAAC_ECS")
        out = src.search_datasets(short_name="MCD43A1", count=1)

        assert len(out) == 1
        kwargs = fake_earthaccess.search_datasets_calls[-1]
        assert kwargs["short_name"] == "MCD43A1"
        assert kwargs["provider"] == "LPDAAC_ECS"
        assert kwargs["count"] == 1

    def test_search_granules_formats_bbox_and_temporal(self, fake_earthaccess):
        src = EarthAccessSource()
        granules = src.search_granules(
            short_name="MCD43A1",
            bounds=(10.0, 20.0, 11.0, 21.0),
            crs="EPSG:4326",
            temporal=(datetime(2024, 1, 1), datetime(2024, 1, 3)),
            count=1,
        )

        assert len(granules) == 1
        kwargs = fake_earthaccess.search_data_calls[-1]
        assert kwargs["short_name"] == "MCD43A1"
        assert kwargs["bounding_box"] == (10.0, 20.0, 11.0, 21.0)
        assert kwargs["temporal"] == ("2024-01-01T00:00:00", "2024-01-03T00:00:00")
        assert kwargs["count"] == 1

    def test_search_count_zero_returns_early(self, fake_earthaccess):
        src = EarthAccessSource()
        assert src.search_datasets(short_name="MCD43A1", count=0) == []
        assert src.search_granules(short_name="MCD43A1", count=0) == []
        assert fake_earthaccess.search_datasets_calls == []
        assert fake_earthaccess.search_data_calls == []

    def test_open_and_download(self, fake_earthaccess, tmp_path: Path):
        src = EarthAccessSource()
        opened = src.open_granules([{"id": "a"}, {"id": "b"}])
        assert opened == ["opened:2"]
        assert src.is_authenticated is True
        assert fake_earthaccess.login_calls == 1
        assert fake_earthaccess.last_login_kwargs == {"strategy": "all", "persist": False}

        out = src.download_granules([{"id": "a"}], tmp_path / "cache")
        assert len(out) == 1
        assert out[0].name == "g1.hdf"
        assert (tmp_path / "cache").exists()
        assert fake_earthaccess.download_calls[-1][2] == 8

    def test_login_uses_temporary_environment_credentials(self, fake_earthaccess, monkeypatch):
        monkeypatch.delenv("EARTHDATA_USERNAME", raising=False)
        monkeypatch.delenv("EARTHDATA_PASSWORD", raising=False)

        def _login(**kwargs):
            fake_earthaccess.login_calls += 1
            fake_earthaccess.last_login_kwargs = kwargs
            fake_earthaccess.last_env = {
                "EARTHDATA_USERNAME": __import__("os").environ.get("EARTHDATA_USERNAME"),
                "EARTHDATA_PASSWORD": __import__("os").environ.get("EARTHDATA_PASSWORD"),
            }

        fake_earthaccess.login = _login

        src = EarthAccessSource(
            earthdata_username="user",
            earthdata_password="secret",
        )
        src.open_granules([{"id": "a"}])

        assert fake_earthaccess.login_calls == 1
        assert fake_earthaccess.last_login_kwargs["strategy"] == "environment"
        assert fake_earthaccess.last_env == {
            "EARTHDATA_USERNAME": "user",
            "EARTHDATA_PASSWORD": "secret",
        }

    def test_temporal_window_helper(self):
        t = datetime(2024, 1, 15, 12, 30)
        temporal = EarthAccessSource.temporal_window(t, 2)
        assert temporal[0].startswith("2024-01-13")
        assert temporal[1] == "2024-01-17T12:30:00Z"

    def test_bounds_normalization_identity(self):
        bounds = (1.0, 2.0, 3.0, 4.0)
        out = EarthAccessSource.normalize_bounds_to_wgs84(bounds, "EPSG:4326")
        assert out == bounds
        assert EarthAccessSource.to_cmr_bounding_box(out) == "1.0,2.0,3.0,4.0"


def test_earthaccess_source_from_auth_without_manager_uses_plain_source():
    source = earthdata_mod.earthaccess_source_from_auth(None, provider="LPDAAC_ECS")

    assert isinstance(source, EarthAccessSource)
    assert source.provider == "LPDAAC_ECS"


def test_earthaccess_source_from_auth_uses_manager_builder():
    expected = EarthAccessSource(provider="LPDAAC_ECS")

    class _FakeEarthdataAuth:
        def build_earthaccess_source(self, *, provider: str | None = None) -> EarthAccessSource:
            assert provider == "LPDAAC_ECS"
            return expected

    class _FakeManager:
        def earthdata(self) -> _FakeEarthdataAuth:
            return _FakeEarthdataAuth()

    assert (
        earthdata_mod.earthaccess_source_from_auth(_FakeManager(), provider="LPDAAC_ECS")
        is expected
    )


def test_build_earthaccess_runtime_normalizes_cache_dir_and_keeps_supplied_dependencies(
    tmp_path: Path,
):
    source = EarthAccessSource(provider="LPDAAC_ECS")
    catalog = EarthAccessCatalog(source=source)

    runtime = earthdata_mod.build_earthaccess_runtime(
        cache_dir=tmp_path / "cache",
        source=source,
        catalog=catalog,
        provider="IGNORED",
    )

    assert runtime.cache_dir == (tmp_path / "cache").expanduser()
    assert runtime.source is source
    assert runtime.catalog is catalog
    assert (
        earthdata_mod.earthaccess_cache_dir(runtime.cache_dir, "MCD43A1")
        == (tmp_path / "cache").expanduser()
    )


def test_select_candidate_paths_orders_tiles_and_filters_sample_days(monkeypatch, tmp_path: Path):
    p1 = tmp_path / "tile_a.hdf"
    p2 = tmp_path / "tile_b.hdf"
    p3 = tmp_path / "tile_c.hdf"

    timestamps = {
        p1: datetime(2024, 1, 8, 12, 0, 0),
        p2: datetime(2024, 1, 1, 12, 0, 0),
        p3: datetime(2024, 1, 8, 10, 0, 0),
    }
    tiles = {
        p1: (30, 7),
        p2: (29, 7),
        p3: (29, 7),
    }

    monkeypatch.setattr(earthdata_mod, "parse_granule_date", lambda path: timestamps[path])
    monkeypatch.setattr(earthdata_mod, "granule_intersects_bounds", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(earthdata_mod, "parse_tile_indices", lambda path: tiles[path])

    selected = earthdata_mod.select_candidate_paths(
        [p1, p2, p3],
        obs_time=datetime(2024, 1, 8, 11, 0, 0),
        bounds=(0.0, 0.0, 1.0, 1.0),
        crs="EPSG:4326",
        sample_dates=np.array(["2024-01-08"], dtype="datetime64[D]"),
    )

    assert selected == [p3, p1]
