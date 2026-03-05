"""Unit tests for S2 convenience entrypoints in siac.siac."""

from __future__ import annotations

from datetime import date, datetime
from pathlib import Path

import pytest

from siac.core.config import SIACConfig
from siac.core.exceptions import DataNotFoundError
from siac.io.s2_data_source import S2Product, S2Query
from siac.siac import resolve_s2_input, search_sentinel2, siac_process_s2


def test_resolve_s2_input_returns_existing_local_path(tmp_path: Path):
    safe_dir = tmp_path / "S2A_MSIL1C_20240101T000000_N0500_R001_T31UDQ_20240101T000001.SAFE"
    safe_dir.mkdir(parents=True, exist_ok=True)
    cfg = SIACConfig(sensor="s2", s2_data={"backend": "local"})
    out = resolve_s2_input(safe_dir, cfg)
    assert out == safe_dir


def test_resolve_s2_input_local_backend_missing_path_raises():
    cfg = SIACConfig(sensor="s2", s2_data={"backend": "local"})
    with pytest.raises(DataNotFoundError, match="backend is 'local'"):
        resolve_s2_input("T31UDQ_20240101", cfg)


def test_resolve_s2_input_uses_remote_backend_and_cache(monkeypatch, tmp_path: Path):
    cfg = SIACConfig(
        sensor="s2",
        s2_data={"backend": "gcs", "cache_dir": tmp_path, "max_cloud_cover": 35.0},
    )
    captured: dict[str, object] = {}

    def _fake_get(self, query, dest_dir=None):  # noqa: ANN001
        captured["query"] = query
        captured["dest_dir"] = dest_dir
        safe = Path(dest_dir) / "fake.SAFE"
        safe.mkdir(parents=True, exist_ok=True)
        return safe

    monkeypatch.setattr("siac.io.s2_data_source.S2DataAccess.get", _fake_get)
    out = resolve_s2_input("T31UDQ_20240101", cfg)

    assert out == tmp_path / "fake.SAFE"
    assert captured["dest_dir"] == tmp_path
    q = captured["query"]
    assert isinstance(q, S2Query)
    assert q.mgrs_tile == "31UDQ"
    assert q.date == date(2024, 1, 1)
    assert q.max_cloud_cover == 35.0


def test_search_sentinel2_builds_query_and_calls_search(monkeypatch):
    fake_product = S2Product(
        product_id="S2A_MSIL1C_20260102T024121_N0511_R089_T50QLD_20260102T035433",
        mgrs_tile="50QLD",
        sensing_date=datetime(2026, 1, 2, 2, 41, 21),
        processing_baseline="N0511",
        cloud_cover=10.0,
        satellite="S2A",
        orbit_number=89,
        source_url="gs://gcp-public-data-sentinel-2/tiles/50/Q/LD/foo.SAFE/",
    )
    captured: dict[str, object] = {}

    def _fake_search(backend, query):  # noqa: ANN001
        captured["backend"] = backend
        captured["query"] = query
        return [fake_product]

    monkeypatch.setattr("siac.io.s2_data_source.search_s2", _fake_search)
    products = search_sentinel2(
        tile="50QLD",
        start_date="2026-01-02",
        end_date="2026-01-03",
        max_cloud_cover=20.0,
        backend="gcs",
    )

    assert len(products) == 1
    q = captured["query"]
    assert isinstance(q, S2Query)
    assert q.mgrs_tile == "50QLD"
    assert q.start_date == date(2026, 1, 2)
    assert q.end_date == date(2026, 1, 3)
    assert q.max_cloud_cover == 20.0


def test_siac_process_s2_resolves_and_runs_process(monkeypatch, tmp_path: Path):
    cfg = SIACConfig(sensor="s2", s2_data={"backend": "local"})
    safe_dir = tmp_path / "S2A_fake.SAFE"
    safe_dir.mkdir(parents=True, exist_ok=True)
    captured: dict[str, object] = {}

    monkeypatch.setattr("siac.siac.resolve_s2_input", lambda *_args, **_kwargs: safe_dir)

    def _fake_process(self, input_path, output_path=None):  # noqa: ANN001
        captured["input_path"] = Path(input_path)
        captured["output_path"] = output_path
        return "ok"

    monkeypatch.setattr("siac.siac.SIAC.process", _fake_process)
    out = siac_process_s2(cfg, "T31UDQ_20240101", output_path=tmp_path / "out")

    assert out == "ok"
    assert captured["input_path"] == safe_dir
    assert captured["output_path"] == tmp_path / "out"
