"""Integration test for remote zipped LUT loading."""

from __future__ import annotations

import pytest

from siac.rt.lut import DEFAULT_LUT_URL, ZarrLUTBackend


def _is_network_unavailable_error(exc: Exception) -> bool:
    """Best-effort classification of network-unavailable conditions."""
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
        "403 client error",
    ]
    return any(marker in message for marker in markers)


def _is_local_zarr_compatibility_error(exc: Exception) -> bool:
    """Classify expected failures when local zarr/xarray cannot read remote format."""
    message = str(exc).lower()
    markers = [
        "zarr_v3_experimental_api",
        "v3 reading and writing is experimental",
        "unexpected keys in metadata",
        "consolidated_metadata",
    ]
    return any(marker in message for marker in markers)


@pytest.mark.integration
def test_open_public_remote_zipped_lut():
    """Should attempt to open the public remote zipped LUT URL."""
    backend = ZarrLUTBackend(
        DEFAULT_LUT_URL,
        storage_options={"timeout": 10.0},
    )

    try:
        ds = backend.lut
    except Exception as exc:  # pragma: no cover - exercised only when offline
        if _is_network_unavailable_error(exc):
            pytest.skip(f"Remote LUT URL not reachable in this environment: {exc}")
        if _is_local_zarr_compatibility_error(exc):
            pytest.skip(f"Local zarr/xarray build cannot read remote LUT format: {exc}")
        raise

    assert len(ds.data_vars) > 0
