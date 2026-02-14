"""Unit tests for EarthAccessSource wrapper behavior."""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import pytest

from siac.io.earthaccess_source import EarthAccessSource


class _FakeEarthAccessModule:
    def __init__(self):
        self.login_calls = 0
        self.search_datasets_calls: list[dict] = []
        self.search_data_calls: list[dict] = []
        self.open_calls: list[list] = []
        self.download_calls: list[tuple[list, str]] = []

    def login(self, **kwargs):
        self.login_calls += 1
        self.last_login_kwargs = kwargs

    def search_datasets(self, **kwargs):
        self.search_datasets_calls.append(kwargs)
        if kwargs.get("short_name") == "MCD43A1":
            return [{"short_name": "MCD43A1"}]
        return [{"umm": {"ShortName": "MCD19A2"}}]

    def search_data(self, **kwargs):
        self.search_data_calls.append(kwargs)
        return [{"id": "g1"}, {"id": "g2"}]

    def open(self, granules):
        self.open_calls.append(list(granules))
        return [f"opened:{len(granules)}"]

    def download(self, granules, dest_dir):
        self.download_calls.append((list(granules), dest_dir))
        return [str(Path(dest_dir) / "g1.hdf")]


@pytest.fixture
def fake_earthaccess(monkeypatch):
    fake = _FakeEarthAccessModule()
    monkeypatch.setitem(sys.modules, "earthaccess", fake)
    return fake


class TestEarthAccessSource:
    def test_lazy_authentication(self, fake_earthaccess):
        src = EarthAccessSource(provider="LPDAAC_ECS")
        assert src.is_authenticated is False
        assert fake_earthaccess.login_calls == 0

        _ = src.search_datasets(keyword="MCD43")
        assert src.is_authenticated is True
        assert fake_earthaccess.login_calls == 1

        _ = src.search_datasets(keyword="MCD43")
        assert fake_earthaccess.login_calls == 1

    def test_search_datasets_forwards_query(self, fake_earthaccess):
        src = EarthAccessSource(provider="LPDAAC_ECS")
        out = src.search_datasets(short_name="MCD43A1", count=1)

        assert len(out) == 1
        kwargs = fake_earthaccess.search_datasets_calls[-1]
        assert kwargs["short_name"] == "MCD43A1"
        assert kwargs["provider"] == "LPDAAC_ECS"

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
        assert kwargs["bounding_box"] == "10.0,20.0,11.0,21.0"
        assert "2024-01-01" in kwargs["temporal"]

    def test_open_and_download(self, fake_earthaccess, tmp_path: Path):
        src = EarthAccessSource()
        opened = src.open_granules([{"id": "a"}, {"id": "b"}])
        assert opened == ["opened:2"]

        out = src.download_granules([{"id": "a"}], tmp_path / "cache")
        assert len(out) == 1
        assert out[0].name == "g1.hdf"
        assert (tmp_path / "cache").exists()

    def test_temporal_window_helper(self):
        t = datetime(2024, 1, 15, 12, 30)
        temporal = EarthAccessSource.temporal_window(t, 2)
        assert temporal.startswith("2024-01-13")
        assert temporal.endswith("2024-01-17T12:30:00Z")

    def test_bounds_normalization_identity(self):
        bounds = (1.0, 2.0, 3.0, 4.0)
        out = EarthAccessSource.normalize_bounds_to_wgs84(bounds, "EPSG:4326")
        assert out == bounds
        assert EarthAccessSource.to_cmr_bounding_box(out) == "1.0,2.0,3.0,4.0"
