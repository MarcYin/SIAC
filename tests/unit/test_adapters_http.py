"""Unit tests for the shared adapter HTTP helper (REVIEW.md §2.6)."""

from __future__ import annotations

import requests
from urllib3.util.retry import Retry

from siac.adapters import _http


def test_make_session_mounts_retry_adapter_on_http_and_https() -> None:
    session = _http.make_session()
    try:
        http_adapter = session.get_adapter("http://example.test/")
        https_adapter = session.get_adapter("https://example.test/")
        # Both prefixes share the same HTTPAdapter instance.
        assert http_adapter is https_adapter

        retry = http_adapter.max_retries
        assert isinstance(retry, Retry)
        # Configured retries cover the documented status codes.
        assert retry.total == _http.DEFAULT_TOTAL_RETRIES
        assert retry.backoff_factor == _http.DEFAULT_BACKOFF_FACTOR
        for code in (429, 500, 502, 503, 504):
            assert code in retry.status_forcelist
        # ``raise_on_status`` is False so callers see the final response and
        # can inspect ``status_code`` (matches the existing call sites).
        assert retry.raise_on_status is False
        assert retry.respect_retry_after_header is True
        # POST is included so token exchanges retry on 5xx.
        assert "POST" in retry.allowed_methods
    finally:
        session.close()


def test_make_session_stashes_default_timeout() -> None:
    session = _http.make_session(timeout_default=12.5)
    try:
        assert getattr(session, _http._TIMEOUT_ATTR) == 12.5
    finally:
        session.close()


def test_request_with_default_timeout_applies_session_default() -> None:
    captured: dict[str, object] = {}

    class _StubSession:
        def request(self, method: str, url: str, **kwargs):  # noqa: ANN003
            captured["method"] = method
            captured["url"] = url
            captured["kwargs"] = kwargs
            return "stub-response"

    stub = _StubSession()
    setattr(stub, _http._TIMEOUT_ATTR, 7.5)

    out = _http.request_with_default_timeout(stub, "GET", "https://example.test/")
    assert out == "stub-response"
    assert captured["kwargs"] == {"timeout": 7.5}

    captured.clear()
    _http.request_with_default_timeout(stub, "POST", "https://example.test/", timeout=2.0)
    # Caller-provided timeout wins over the session default.
    assert captured["kwargs"]["timeout"] == 2.0


def test_make_session_returns_requests_session() -> None:
    session = _http.make_session()
    try:
        assert isinstance(session, requests.Session)
    finally:
        session.close()


def test_make_session_accepts_custom_status_forcelist() -> None:
    session = _http.make_session(status_forcelist=(418,), total_retries=1)
    try:
        retry = session.get_adapter("https://example.test/").max_retries
        assert isinstance(retry, Retry)
        assert tuple(retry.status_forcelist) == (418,)
        assert retry.total == 1
    finally:
        session.close()
