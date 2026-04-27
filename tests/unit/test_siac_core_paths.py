"""High-value helper coverage for the request-driven public API."""

from __future__ import annotations

from datetime import date, datetime
from pathlib import Path

import pytest

from siac.adapters.data.s2_data_source import S2Query
from siac.api.public import (
    apply_s2_query_defaults,
    coerce_date,
    coerce_s2_query,
    resolve_s2_input,
    search_sentinel2,
    siac_process_s2,
)
from siac.api.requests import (
    Sentinel2ProcessRequest,
    Sentinel2ResolveRequest,
    Sentinel2SearchRequest,
)
from siac.config import SIACConfig
from siac.errors import DataNotFoundError


def test_public_s2_helper_exports_and_local_errors():
    cfg = SIACConfig(sensor="s2", providers={"s2": {"backend": "local"}})

    assert coerce_date(None) is None
    assert coerce_date(date(2026, 1, 2)) == date(2026, 1, 2)
    assert coerce_date(datetime(2026, 1, 2, 8, 0, 0)) == date(2026, 1, 2)
    assert coerce_date("2026-01-03") == date(2026, 1, 3)

    q = S2Query.from_tile_date("T50QLD_20260102")
    q.max_cloud_cover = 100.0
    q.processing_level = "L1C"
    out_q = apply_s2_query_defaults(
        q,
        config=SIACConfig(
            sensor="s2",
            providers={"s2": {"max_cloud_cover": 33.0, "processing_level": "L2A"}},
        ),
    )
    assert out_q.max_cloud_cover == 33.0
    assert out_q.processing_level == "L2A"
    assert coerce_s2_query(out_q, config=cfg) is not out_q

    with pytest.raises(DataNotFoundError, match="backend is 'local'"):
        resolve_s2_input("missing.SAFE", cfg)

    with pytest.raises(ValueError, match="does not support backend='local'"):
        search_sentinel2(tile="50QLD", date_value="2026-01-02", backend="local")


def test_resolve_s2_input_builds_request(monkeypatch):
    cfg = SIACConfig(sensor="s2")
    captured: dict[str, object] = {}

    def _fake_resolve(request):  # noqa: ANN001
        captured["request"] = request
        return Path("/tmp/fake.SAFE")

    monkeypatch.setattr("siac.api.public.app_resolve_s2_input", _fake_resolve)

    out = resolve_s2_input("T31UDQ_20240101", cfg)
    assert out == Path("/tmp/fake.SAFE")
    request = captured["request"]
    assert isinstance(request, Sentinel2ResolveRequest)
    assert request.config is cfg
    assert request.query == "T31UDQ_20240101"


def test_search_sentinel2_builds_request(monkeypatch):
    captured: dict[str, object] = {}

    def _fake_search(request):  # noqa: ANN001
        captured["request"] = request
        return ["ok"]

    monkeypatch.setattr("siac.api.public.app_search_sentinel2", _fake_search)

    out = search_sentinel2(
        tile="50QLD", start_date="2026-01-02", end_date="2026-01-03", backend="gcs"
    )
    assert out == ["ok"]
    request = captured["request"]
    assert isinstance(request, Sentinel2SearchRequest)
    assert request.tile == "50QLD"
    assert request.start_date == "2026-01-02"
    assert request.end_date == "2026-01-03"
    assert request.backend == "gcs"


def test_siac_process_s2_builds_request(monkeypatch, tmp_path: Path):
    cfg = SIACConfig(sensor="s2", providers={"s2": {"backend": "local"}})
    captured: dict[str, object] = {}

    def _fake_process(request):  # noqa: ANN001
        captured["request"] = request
        return "ok"

    monkeypatch.setattr("siac.api.public.workflow_process_s2", _fake_process)

    out = siac_process_s2(
        cfg, "T31UDQ_20240101", output_path=tmp_path / "out", aoi="/tmp/aoi.geojson"
    )
    assert out == "ok"
    request = captured["request"]
    assert isinstance(request, Sentinel2ProcessRequest)
    assert request.config is cfg
    assert request.query == "T31UDQ_20240101"
    assert request.output_path == tmp_path / "out"
    assert request.aoi == "/tmp/aoi.geojson"
