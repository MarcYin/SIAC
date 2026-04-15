"""Integration test for retry logic with mocked earthaccess failures."""

from __future__ import annotations

import pytest

from siac.adapters._retry import retry_transient


def test_retry_with_simulated_earthaccess_503():
    """Simulate an earthaccess HTTP 503 then success."""
    call_count = [0]

    def mock_search(**kwargs):
        call_count[0] += 1
        if call_count[0] <= 2:
            raise ConnectionError("HTTP 503 Service Unavailable")
        return [{"granule": "test"}]

    result = retry_transient(
        mock_search,
        max_attempts=3,
        base_delay_s=0.01,
        label="earthaccess.search",
        short_name="MCD43A1",
    )
    assert result == [{"granule": "test"}]
    assert call_count[0] == 3


def test_retry_with_timeout_then_success():
    """Simulate a timeout on first attempt then success."""
    call_count = [0]

    def mock_download(granules, dest):
        call_count[0] += 1
        if call_count[0] == 1:
            raise TimeoutError("Connection timed out")
        return [f"{dest}/file.hdf"]

    result = retry_transient(
        mock_download,
        "granule1",
        "/tmp/cache",
        max_attempts=3,
        base_delay_s=0.01,
        label="earthaccess.download",
    )
    assert result == ["/tmp/cache/file.hdf"]
    assert call_count[0] == 2


def test_retry_permanent_failure_raises():
    """ValueError (non-transient) should not be retried."""
    def bad_query(**kwargs):
        raise ValueError("Invalid short_name")

    with pytest.raises(ValueError, match="Invalid short_name"):
        retry_transient(
            bad_query,
            max_attempts=3,
            base_delay_s=0.01,
            label="earthaccess.search",
        )
