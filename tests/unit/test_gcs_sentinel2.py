"""Unit tests for the GCS Sentinel-2 backend."""

from __future__ import annotations

from datetime import date, datetime
from pathlib import Path

import pytest

import siac.io.gcs_sentinel2 as gcs_mod
from siac.core.exceptions import DataNotFoundError
from siac.io.s2_data_source import S2Product, S2Query


def _product(product_id: str, source_url: str = "") -> S2Product:
    return S2Product(
        product_id=product_id,
        mgrs_tile="31UDQ",
        sensing_date=datetime(2024, 1, 3, 10, 30, 21),
        processing_baseline="N0500",
        cloud_cover=100.0,
        satellite="S2A",
        orbit_number=51,
        source_url=source_url,
        size_mb=None,
    )


def test_search_gcs_by_tile_filters_date_and_level(monkeypatch):
    prefixes = [
        "tiles/31/U/DQ/S2A_MSIL1C_20240101T103021_N0500_R051_T31UDQ_20240101T120000.SAFE/",
        "tiles/31/U/DQ/S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000.SAFE/",
        "tiles/31/U/DQ/S2A_MSIL2A_20240103T103021_N0500_R051_T31UDQ_20240103T120000.SAFE/",
    ]

    def _fake_list_api(*, prefix, delimiter=None, page_token=None, max_results=5000):  # noqa: ARG001
        assert prefix == "tiles/31/U/DQ/"
        assert delimiter == "/"
        return {"prefixes": prefixes}

    monkeypatch.setattr(gcs_mod, "_list_api", _fake_list_api)

    query = S2Query(
        mgrs_tile="31UDQ",
        date=date(2024, 1, 3),
        processing_level="L1C",
        max_cloud_cover=20.0,  # ignored by GCS backend
    )
    out = gcs_mod.search_gcs(query)
    assert len(out) == 1
    assert out[0].product_id.startswith("S2A_MSIL1C_20240103")
    assert out[0].mgrs_tile == "31UDQ"
    assert out[0].processing_baseline == "N0500"
    assert out[0].orbit_number == 51
    assert out[0].satellite == "S2A"


def test_search_gcs_by_product_id_uses_prefix_exists(monkeypatch):
    def _fake_list_api(*, prefix, delimiter=None, page_token=None, max_results=5000):  # noqa: ARG001
        assert prefix.endswith(".SAFE/")
        return {"items": [{"name": f"{prefix}manifest.safe", "size": "3"}]}

    monkeypatch.setattr(gcs_mod, "_list_api", _fake_list_api)

    query = S2Query(
        product_id="S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000.SAFE"
    )
    out = gcs_mod.search_gcs(query)
    assert len(out) == 1
    assert out[0].product_id.endswith("_T31UDQ_20240103T120000")


def test_search_gcs_bbox_only_raises():
    query = S2Query(bbox=(1.0, 2.0, 3.0, 4.0))
    with pytest.raises(ValueError, match="product_id or mgrs_tile"):
        gcs_mod.search_gcs(query)


def test_download_gcs_materializes_safe_tree(monkeypatch, tmp_path: Path):
    product_id = "S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000"
    safe_prefix = f"tiles/31/U/DQ/{product_id}.SAFE/"
    product = _product(product_id, source_url=f"gs://{gcs_mod.GCS_BUCKET}/{safe_prefix}")

    objects = [
        {"name": f"tiles/31/U/DQ/{product_id}.SAFE_$folder$", "size": "0"},
        {"name": f"{safe_prefix}manifest.safe", "size": "2"},
        {"name": f"{safe_prefix}GRANULE/L1C.xml", "size": "3"},
    ]

    downloaded: list[Path] = []

    def _fake_list_objects_under(prefix):
        assert prefix == safe_prefix
        return objects

    def _fake_download(url, target: Path, timeout=300):  # noqa: ARG001
        downloaded.append(target)
        payload = b"ok" if target.name == "manifest.safe" else b"xml"
        target.write_bytes(payload)

    monkeypatch.setattr(gcs_mod, "_list_objects_under", _fake_list_objects_under)
    monkeypatch.setattr(gcs_mod, "_download_url_to_file", _fake_download)

    safe_dir = gcs_mod.download_gcs(product, tmp_path)
    assert safe_dir == tmp_path / f"{product_id}.SAFE"
    assert (safe_dir / "manifest.safe").exists()
    assert (safe_dir / "GRANULE" / "L1C.xml").exists()
    assert len(downloaded) == 2


def test_download_gcs_skips_fully_downloaded_files(monkeypatch, tmp_path: Path):
    product_id = "S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000"
    safe_prefix = f"tiles/31/U/DQ/{product_id}.SAFE/"
    product = _product(product_id, source_url=f"gs://{gcs_mod.GCS_BUCKET}/{safe_prefix}")
    safe_dir = tmp_path / f"{product_id}.SAFE"
    safe_dir.mkdir(parents=True, exist_ok=True)
    (safe_dir / "manifest.safe").write_bytes(b"ok")  # size=2, should be skipped

    objects = [
        {"name": f"{safe_prefix}manifest.safe", "size": "2"},
        {"name": f"{safe_prefix}GRANULE/L1C.xml", "size": "3"},
    ]
    downloaded: list[Path] = []

    def _fake_list_objects_under(prefix):
        assert prefix == safe_prefix
        return objects

    def _fake_download(url, target: Path, timeout=300):  # noqa: ARG001
        downloaded.append(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(b"xml")

    monkeypatch.setattr(gcs_mod, "_list_objects_under", _fake_list_objects_under)
    monkeypatch.setattr(gcs_mod, "_download_url_to_file", _fake_download)

    out = gcs_mod.download_gcs(product, tmp_path)
    assert out == safe_dir
    assert len(downloaded) == 1
    assert downloaded[0] == safe_dir / "GRANULE" / "L1C.xml"


def test_download_gcs_uses_parallel_pool_for_multiple_files(monkeypatch, tmp_path: Path):
    product_id = "S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000"
    safe_prefix = f"tiles/31/U/DQ/{product_id}.SAFE/"
    product = _product(product_id, source_url=f"gs://{gcs_mod.GCS_BUCKET}/{safe_prefix}")

    objects = [
        {"name": f"{safe_prefix}manifest.safe", "size": "1"},
        {"name": f"{safe_prefix}GRANULE/L1C.xml", "size": "1"},
        {"name": f"{safe_prefix}AUX_DATA/aux.txt", "size": "1"},
    ]

    observed: dict[str, int] = {}

    class _FakePool:
        def __init__(self, max_workers: int):
            observed["max_workers"] = max_workers

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):  # noqa: ANN001
            return False

        def map(self, fn, iterable):
            jobs = list(iterable)
            observed["job_count"] = len(jobs)
            for job in jobs:
                fn(job)
            return [None] * len(jobs)

    def _fake_list_objects_under(prefix):
        assert prefix == safe_prefix
        return objects

    def _fake_download(url, target: Path, timeout=300):  # noqa: ARG001
        target.write_text("x")

    monkeypatch.setattr(gcs_mod, "_list_objects_under", _fake_list_objects_under)
    monkeypatch.setattr(gcs_mod, "_download_url_to_file", _fake_download)
    monkeypatch.setattr(gcs_mod, "ThreadPoolExecutor", _FakePool)

    safe_dir = gcs_mod.download_gcs(product, tmp_path)
    assert safe_dir.exists()
    assert observed["job_count"] == 3
    assert observed["max_workers"] > 1


def test_download_gcs_raises_when_no_objects(monkeypatch, tmp_path: Path):
    product_id = "S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000"
    product = _product(product_id)
    monkeypatch.setattr(gcs_mod, "_list_objects_under", lambda prefix: [])  # noqa: ARG005
    with pytest.raises(DataNotFoundError, match="No objects found"):
        gcs_mod.download_gcs(product, tmp_path)


def test_download_gcs_retries_failed_transfer(monkeypatch, tmp_path: Path):
    product_id = "S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000"
    safe_prefix = f"tiles/31/U/DQ/{product_id}.SAFE/"
    product = _product(product_id, source_url=f"gs://{gcs_mod.GCS_BUCKET}/{safe_prefix}")
    objects = [{"name": f"{safe_prefix}manifest.safe", "size": "2"}]

    attempts = {"count": 0}

    def _fake_list_objects_under(prefix):
        assert prefix == safe_prefix
        return objects

    def _fake_download(url, target: Path, timeout=300):  # noqa: ARG001
        attempts["count"] += 1
        if attempts["count"] == 1:
            raise OSError("transient network failure")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(b"ok")

    monkeypatch.setattr(gcs_mod, "_list_objects_under", _fake_list_objects_under)
    monkeypatch.setattr(gcs_mod, "_download_url_to_file", _fake_download)
    monkeypatch.setattr(gcs_mod, "GCS_DOWNLOAD_RETRY_BACKOFF_SEC", 0.0)

    safe_dir = gcs_mod.download_gcs(product, tmp_path)
    assert safe_dir.exists()
    assert (safe_dir / "manifest.safe").read_bytes() == b"ok"
    assert attempts["count"] == 2


def test_gcs_backend_delegates(monkeypatch, tmp_path: Path):
    product_id = "S2A_MSIL1C_20240103T103021_N0500_R051_T31UDQ_20240103T120000"
    product = _product(product_id)
    expected_safe = tmp_path / f"{product_id}.SAFE"
    expected_safe.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(gcs_mod, "search_gcs", lambda query: [product])  # noqa: ARG005
    monkeypatch.setattr(gcs_mod, "download_gcs", lambda p, d: expected_safe)  # noqa: ARG005

    backend = gcs_mod.GCSSentinel2Backend()
    q = S2Query(mgrs_tile="31UDQ", date=date(2024, 1, 3))
    out = backend.search(q)
    got = backend.download(product, tmp_path)

    assert len(out) == 1
    assert out[0].product_id == product_id
    assert got == expected_safe
