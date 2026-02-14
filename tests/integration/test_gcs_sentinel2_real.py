"""Integration test for real Sentinel-2 file download from public GCS."""

from __future__ import annotations

from datetime import date
from pathlib import Path

import pytest

from siac.io.gcs_sentinel2 import (
    _download_with_retry,
    _list_objects_under,
    _object_download_url,
    _resolve_safe_prefix,
    search_gcs,
)
from siac.io.s2_data_source import S2Query

TARGET_PRODUCT_ID = "S2C_MSIL1C_20260102T024121_N0511_R089_T50QLD_20260102T035433"


def _is_network_unavailable_error(exc: Exception) -> bool:
    message = str(exc).lower()
    markers = [
        "name or service not known",
        "nodename nor servname provided",
        "temporary failure in name resolution",
        "network is unreachable",
        "connection refused",
        "connection reset",
        "timed out",
        "timeout",
        "proxyerror",
        "ssl",
        "forbidden",
        "403",
        "429",
    ]
    return any(marker in message for marker in markers)


def _search_or_skip(query: S2Query):
    try:
        return search_gcs(query)
    except Exception as exc:  # pragma: no cover - env-dependent network behavior
        if _is_network_unavailable_error(exc):
            pytest.skip(f"GCS API not reachable in this environment: {exc}")
        raise


def _download_one_real_file_or_skip(product_id: str, tmp_path: Path) -> tuple[Path, int]:
    safe_prefix = _resolve_safe_prefix(
        next(iter(_search_or_skip(S2Query.from_product_id(f"{product_id}.SAFE"))))
    )

    try:
        objects = _list_objects_under(safe_prefix)
    except Exception as exc:  # pragma: no cover - env-dependent network behavior
        if _is_network_unavailable_error(exc):
            pytest.skip(f"GCS object listing not reachable in this environment: {exc}")
        raise

    candidates: list[tuple[int, str]] = []
    for item in objects:
        name = item.get("name")
        size = item.get("size")
        if (
            isinstance(name, str)
            and isinstance(size, str)
            and size.isdigit()
            and not name.endswith("_$folder$")
            and not name.endswith("/")
        ):
            candidates.append((int(size), name))

    assert candidates, f"No downloadable file objects found under {safe_prefix!r}."
    preferred = next(
        (
            pair
            for pair in candidates
            if pair[1].endswith("rep_info/S2_User_Product_Level-1C_Metadata.xsd")
        ),
        None,
    )
    expected_size, object_name = preferred if preferred is not None else min(candidates)

    out_path = tmp_path / object_name.split("/")[-1]
    try:
        _download_with_retry(
            _object_download_url(object_name),
            out_path,
            expected_size=expected_size,
            retries=2,
            backoff_sec=0.5,
        )
    except Exception as exc:  # pragma: no cover - env-dependent network behavior
        if _is_network_unavailable_error(exc):
            pytest.skip(f"GCS object download failed due to environment/network: {exc}")
        raise
    return out_path, expected_size


@pytest.mark.integration
def test_gcs_real_search_by_product_id_and_download_file(tmp_path):
    """Validate product-id search independently, then download one real file."""
    products = _search_or_skip(S2Query.from_product_id(f"{TARGET_PRODUCT_ID}.SAFE"))
    assert products, "Expected product-id search to return at least one product."
    assert any(p.product_id == TARGET_PRODUCT_ID for p in products)

    out_path, expected_size = _download_one_real_file_or_skip(TARGET_PRODUCT_ID, tmp_path)
    assert out_path.exists()
    assert out_path.stat().st_size == expected_size


@pytest.mark.integration
def test_gcs_real_search_by_tile_date_window_and_download_file(tmp_path):
    """Validate tile/date-range search independently, then download one real file."""
    products = _search_or_skip(
        S2Query(
            mgrs_tile="50QLD",
            start_date=date(2026, 1, 2),
            end_date=date(2026, 1, 3),
            processing_level="L1C",
        )
    )
    assert products, "Expected tile/date-range search to return at least one product."
    assert any(p.product_id == TARGET_PRODUCT_ID for p in products)

    target = next((p for p in products if p.product_id == TARGET_PRODUCT_ID), products[0])
    out_path, expected_size = _download_one_real_file_or_skip(target.product_id, tmp_path)
    assert out_path.exists()
    assert out_path.stat().st_size == expected_size
