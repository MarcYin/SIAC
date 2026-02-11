"""Unit tests for CDSE Sentinel-2 backend."""

from __future__ import annotations

from datetime import date, datetime
from io import BytesIO
from pathlib import Path
import zipfile

import pytest
import requests

import siac.io.copernicus_dataspace as cdse
from siac.core.exceptions import DataNotFoundError
from siac.io.copernicus_dataspace import CopernicusDataspaceBackend, download_cdse, search_cdse
from siac.io.s2_data_source import S2Product, S2Query


class _FakeResponse:
    def __init__(
        self,
        *,
        status_code: int = 200,
        json_data: dict | None = None,
        raw_bytes: bytes = b"",
    ):
        self.status_code = status_code
        self._json_data = json_data or {}
        self._raw_bytes = raw_bytes

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}", response=self)

    def json(self) -> dict:
        return self._json_data

    def iter_content(self, chunk_size: int = 1024 * 1024):
        for i in range(0, len(self._raw_bytes), chunk_size):
            yield self._raw_bytes[i : i + chunk_size]


def _stac_item(
    product_id: str,
    *,
    cloud: float = 10.0,
    mgrs_tile: str = "31UDQ",
    sensing_dt: str = "2024-01-01T10:31:00Z",
    product_href: str = "https://download.dataspace.copernicus.eu/odata/v1/Products(abc)/$value",
) -> dict:
    return {
        "id": product_id,
        "properties": {
            "datetime": sensing_dt,
            "eo:cloud_cover": cloud,
            "s2:mgrs_tile": mgrs_tile,
        },
        "assets": {
            "Product": {
                "href": product_href,
                "file:size": 100 * 1024 * 1024,
            }
        },
    }


def test_search_cdse_product_id_lookup(monkeypatch):
    product_id = "S2B_MSIL1C_20240102T103129_N0500_R108_T31UDQ_20240102T120000"

    def _fake_get(url: str, timeout: int = 60, **kwargs):  # noqa: ARG001
        assert product_id in url
        return _FakeResponse(json_data=_stac_item(product_id))

    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.get", _fake_get)

    q = S2Query.from_product_id(product_id)
    products = search_cdse(q)
    assert len(products) == 1
    p = products[0]
    assert p.product_id == product_id
    assert p.mgrs_tile == "31UDQ"
    assert p.processing_baseline == "N0500"
    assert p.orbit_number == 108
    assert p.satellite == "S2B"
    assert p.cloud_cover == pytest.approx(10.0)
    assert p.size_mb is not None and p.size_mb > 0


def test_search_cdse_filters_tile_cloud_and_pagination(monkeypatch):
    page_1 = {
        "features": [
            _stac_item(
                "S2A_MSIL1C_20240101T103111_N0500_R008_T31UDQ_20240101T120000",
                cloud=11.0,
            ),
        ],
        "links": [{"rel": "next", "href": "https://example.test/next"}],
    }
    page_2 = {
        "features": [
            _stac_item(
                "S2B_MSIL1C_20240101T103111_N0500_R108_T31TCJ_20240101T120000",
                cloud=5.0,
                mgrs_tile="31TCJ",
            ),
            _stac_item(
                "S2B_MSIL1C_20240101T103111_N0500_R108_T31UDQ_20240101T120000",
                cloud=55.0,
            ),
        ],
        "links": [],
    }

    def _fake_post(url: str, json: dict, timeout: int = 60, **kwargs):  # noqa: ARG001
        assert url.endswith("/search")
        assert json["collections"] == ["sentinel-2-l1c"]
        return _FakeResponse(json_data=page_1)

    def _fake_get(url: str, timeout: int = 60, **kwargs):  # noqa: ARG001
        assert url == "https://example.test/next"
        return _FakeResponse(json_data=page_2)

    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.post", _fake_post)
    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.get", _fake_get)

    q = S2Query(
        mgrs_tile="31UDQ",
        start_date=date(2024, 1, 1),
        end_date=date(2024, 1, 1),
        max_cloud_cover=20.0,
    )
    products = search_cdse(q)
    assert len(products) == 1
    assert products[0].mgrs_tile == "31UDQ"
    assert products[0].cloud_cover <= 20.0


def test_search_cdse_product_404_returns_empty(monkeypatch):
    def _fake_get(url: str, timeout: int = 60, **kwargs):  # noqa: ARG001
        return _FakeResponse(status_code=404)

    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.get", _fake_get)
    products = search_cdse(
        S2Query.from_product_id("S2B_MSIL1C_20240101T000000_N0500_R000_T31UDQ_20240101T000000")
    )
    assert products == []


def test_download_cdse_extracts_safe_zip(monkeypatch, tmp_path: Path):
    product_id = "S2A_MSIL1C_20240101T103101_N0500_R008_T31UDQ_20240101T120000"
    product = S2Product(
        product_id=product_id,
        mgrs_tile="31UDQ",
        sensing_date=datetime(2024, 1, 1, 10, 31, 1),
        processing_baseline="N0500",
        cloud_cover=12.0,
        satellite="S2A",
        orbit_number=8,
        source_url="https://download.example/products/1/$value",
    )

    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, mode="w") as zf:
        zf.writestr(f"{product_id}.SAFE/manifest.safe", "ok")
    payload = zip_buffer.getvalue()

    def _fake_get(url: str, headers: dict | None = None, **kwargs):  # noqa: ARG001
        assert url == product.source_url
        assert headers == {"Authorization": "Bearer token"}
        return _FakeResponse(raw_bytes=payload)

    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.get", _fake_get)

    safe_path = download_cdse(product, tmp_path, access_key="token")
    assert safe_path.exists()
    assert safe_path.name.endswith(".SAFE")
    assert (safe_path / "manifest.safe").exists()


def test_cdse_backend_delegates(monkeypatch, tmp_path: Path):
    product = S2Product(
        product_id="S2A_MSIL1C_20240101T103101_N0500_R008_T31UDQ_20240101T120000",
        mgrs_tile="31UDQ",
        sensing_date=datetime(2024, 1, 1, 10, 31, 1),
        processing_baseline="N0500",
        cloud_cover=12.0,
        satellite="S2A",
        orbit_number=8,
        source_url="https://example.test",
    )

    def _fake_search(query, access_key=None, secret_key=None):  # noqa: ARG001
        return [product]

    def _fake_download(prod, dest, access_key=None, secret_key=None):  # noqa: ARG001
        out = Path(dest) / f"{prod.product_id}.SAFE"
        out.mkdir(parents=True, exist_ok=True)
        return out

    monkeypatch.setattr("siac.io.copernicus_dataspace.search_cdse", _fake_search)
    monkeypatch.setattr("siac.io.copernicus_dataspace.download_cdse", _fake_download)

    backend = CopernicusDataspaceBackend(access_key="ak", secret_key="sk")
    found = backend.search(S2Query(mgrs_tile="31UDQ"))
    assert len(found) == 1
    safe = backend.download(found[0], tmp_path)
    assert safe.exists()


def test_datetime_range_builder_variants():
    q_day = S2Query(mgrs_tile="31UDQ", date=date(2024, 1, 2))
    q_range = S2Query(mgrs_tile="31UDQ", start_date=date(2024, 1, 2), end_date=date(2024, 1, 3))
    q_start_only = S2Query(mgrs_tile="31UDQ", start_date=date(2024, 1, 2))
    q_end_only = S2Query(mgrs_tile="31UDQ", end_date=date(2024, 1, 3))
    q_none = S2Query(mgrs_tile="31UDQ")

    assert cdse._to_datetime_range(q_day) == "2024-01-02T00:00:00Z/2024-01-02T23:59:59Z"
    assert cdse._to_datetime_range(q_range) == "2024-01-02T00:00:00Z/2024-01-03T23:59:59Z"
    assert cdse._to_datetime_range(q_start_only) == "2024-01-02T00:00:00Z/.."
    assert cdse._to_datetime_range(q_end_only) == "../2024-01-03T23:59:59Z"
    assert cdse._to_datetime_range(q_none) is None


def test_item_to_product_fallbacks_and_defaults():
    item = {
        "id": "X_MSIL1C_20240101T103101_ABC_T31UDQ_20240101T120000",
        "properties": {},
        "assets": {"Product": {"href": "https://example.test/value"}},
    }
    p = cdse._item_to_product(item)
    assert p.mgrs_tile == "31UDQ"
    assert p.sensing_date == datetime(2024, 1, 1, 10, 31, 1)
    assert p.processing_baseline == "N0000"
    assert p.orbit_number == 0
    assert p.satellite == "S2"
    assert p.cloud_cover == pytest.approx(100.0)
    assert p.size_mb is None


def test_item_to_product_error_paths():
    with pytest.raises(DataNotFoundError, match="has no 'Product' asset"):
        cdse._item_to_product({"id": "X_MSIL1C_20240101T103101", "properties": {}, "assets": {}})

    with pytest.raises(DataNotFoundError, match="has no href"):
        cdse._item_to_product(
            {"id": "X_MSIL1C_20240101T103101", "properties": {}, "assets": {"Product": {"title": "x"}}}
        )

    with pytest.raises(ValueError, match="Cannot parse sensing datetime"):
        cdse._parse_iso_datetime(None, "BAD")


def test_search_payload_and_next_link_helpers():
    q = S2Query(
        mgrs_tile="31UDQ",
        start_date=date(2024, 1, 1),
        end_date=date(2024, 1, 2),
        bbox=(0.0, 1.0, 2.0, 3.0),
        max_cloud_cover=35.0,
    )
    payload = cdse._search_payload(q, limit=17)
    assert payload["collections"] == ["sentinel-2-l1c"]
    assert payload["limit"] == 17
    assert payload["bbox"] == [0.0, 1.0, 2.0, 3.0]
    assert payload["datetime"] == "2024-01-01T00:00:00Z/2024-01-02T23:59:59Z"
    assert payload["filter"]["op"] == "<="

    q_no_filter = S2Query(mgrs_tile="31UDQ", max_cloud_cover=100.0)
    payload_no_filter = cdse._search_payload(q_no_filter)
    assert "filter" not in payload_no_filter
    assert "filter-lang" not in payload_no_filter

    assert cdse._next_link({"links": [{"rel": "self", "href": "x"}]}) is None
    assert cdse._next_link({"links": [{"rel": "next", "href": ""}]}) is None
    assert cdse._next_link({"links": [{"rel": "next", "href": "https://example.test/n"}]}) == "https://example.test/n"


def test_auth_header_resolution_modes(monkeypatch):
    assert cdse._resolve_auth_header("token", None) == {"Authorization": "Bearer token"}
    assert cdse._resolve_auth_header(None, None) == {}

    called = {"ok": False}

    def _fake_token(username: str, password: str, timeout: int = 60):  # noqa: ARG001
        called["ok"] = (username, password) == ("user", "pass")
        return "issued-token"

    monkeypatch.setattr("siac.io.copernicus_dataspace._token_from_credentials", _fake_token)
    assert cdse._resolve_auth_header("user", "pass") == {"Authorization": "Bearer issued-token"}
    assert called["ok"] is True


def test_token_from_credentials_missing_token_raises(monkeypatch):
    def _fake_post(url: str, **kwargs):  # noqa: ARG001
        return _FakeResponse(json_data={"not_access_token": "x"})

    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.post", _fake_post)
    with pytest.raises(DataNotFoundError, match="access_token"):
        cdse._token_from_credentials("user", "pass")


def test_search_cdse_processing_level_filter_and_product_lookup_error(monkeypatch):
    page = {
        "features": [_stac_item("S2A_MSIL1C_20240101T103101_N0500_R008_T31UDQ_20240101T120000")],
        "links": [],
    }

    def _fake_post(url: str, json: dict, timeout: int = 60, **kwargs):  # noqa: ARG001
        return _FakeResponse(json_data=page)

    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.post", _fake_post)
    q = S2Query(mgrs_tile="31UDQ", processing_level="L2A")
    assert search_cdse(q) == []

    def _fake_get_500(url: str, timeout: int = 60, **kwargs):  # noqa: ARG001
        return _FakeResponse(status_code=500)

    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.get", _fake_get_500)
    with pytest.raises(requests.HTTPError):
        search_cdse(
            S2Query.from_product_id(
                "S2B_MSIL1C_20240101T000000_N0500_R000_T31UDQ_20240101T000000"
            )
        )


def test_download_cdse_existing_safe_short_circuit(monkeypatch, tmp_path: Path):
    product = S2Product(
        product_id="S2A_MSIL1C_20240101T103101_N0500_R008_T31UDQ_20240101T120000",
        mgrs_tile="31UDQ",
        sensing_date=datetime(2024, 1, 1, 10, 31, 1),
        processing_baseline="N0500",
        cloud_cover=12.0,
        satellite="S2A",
        orbit_number=8,
        source_url="https://download.example/products/1/$value",
    )
    safe = tmp_path / f"{product.product_id}.SAFE"
    safe.mkdir(parents=True, exist_ok=True)

    called = {"n": 0}

    def _fake_get(url: str, **kwargs):  # noqa: ARG001
        called["n"] += 1
        raise AssertionError("requests.get should not be called when SAFE already exists")

    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.get", _fake_get)
    out = download_cdse(product, tmp_path)
    assert out == safe
    assert called["n"] == 0


def test_download_cdse_fallback_candidate_and_missing_safe(monkeypatch, tmp_path: Path):
    product = S2Product(
        product_id="S2A_MSIL1C_20240101T103101_N0500_R008_T31UDQ_20240101T120000",
        mgrs_tile="31UDQ",
        sensing_date=datetime(2024, 1, 1, 10, 31, 1),
        processing_baseline="N0500",
        cloud_cover=12.0,
        satellite="S2A",
        orbit_number=8,
        source_url="https://download.example/products/1/$value",
    )

    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, mode="w") as zf:
        zf.writestr(f"{product.product_id}_ALT.SAFE/manifest.safe", "ok")
    payload_with_alt = zip_buffer.getvalue()

    def _fake_get_alt(url: str, headers: dict | None = None, **kwargs):  # noqa: ARG001
        return _FakeResponse(raw_bytes=payload_with_alt)

    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.get", _fake_get_alt)
    out = download_cdse(product, tmp_path / "with_alt")
    assert out.name == f"{product.product_id}_ALT.SAFE"
    assert (out / "manifest.safe").exists()

    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, mode="w") as zf:
        zf.writestr("README.txt", "no safe")
    payload_no_safe = zip_buffer.getvalue()

    def _fake_get_no_safe(url: str, headers: dict | None = None, **kwargs):  # noqa: ARG001
        return _FakeResponse(raw_bytes=payload_no_safe)

    monkeypatch.setattr("siac.io.copernicus_dataspace.requests.get", _fake_get_no_safe)
    with pytest.raises(DataNotFoundError, match="SAFE directory was not found"):
        download_cdse(product, tmp_path / "missing_safe")
