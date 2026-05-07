"""Shared HTTP session factory for adapter-layer network usage.

This module centralizes the construction of ``requests.Session`` objects with
retry/backoff configured for transient HTTP errors. Adapter modules that need
to make outbound HTTP calls should obtain their session from
:func:`make_session` so that retry, backoff, and default timeout behaviour are
consistent across the codebase (REVIEW.md §2.6).
"""

from __future__ import annotations

import logging
from typing import Any

import requests  # type: ignore[import-untyped]
from requests.adapters import HTTPAdapter  # type: ignore[import-untyped]
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

DEFAULT_TIMEOUT_SECONDS: float = 60.0
DEFAULT_TOTAL_RETRIES: int = 3
DEFAULT_BACKOFF_FACTOR: float = 0.5
DEFAULT_STATUS_FORCELIST: tuple[int, ...] = (429, 500, 502, 503, 504)

# Stash key used to attach the default timeout to a Session instance. Callers
# that want timeout-aware requests can either pass ``timeout=`` explicitly or
# go through :func:`request_with_default_timeout`.
_TIMEOUT_ATTR = "_siac_default_timeout"


def make_session(
    *,
    total_retries: int = DEFAULT_TOTAL_RETRIES,
    backoff_factor: float = DEFAULT_BACKOFF_FACTOR,
    status_forcelist: tuple[int, ...] = DEFAULT_STATUS_FORCELIST,
    timeout_default: float = DEFAULT_TIMEOUT_SECONDS,
) -> requests.Session:
    """Return a ``requests.Session`` configured with retry/backoff.

    The returned session mounts an :class:`~requests.adapters.HTTPAdapter`
    with a :class:`~urllib3.util.retry.Retry` policy on both ``http://`` and
    ``https://``. Retries are applied for the standard transient failure HTTP
    statuses (``429``, ``500``, ``502``, ``503``, ``504``) on idempotent
    methods plus ``POST`` (CDSE token exchange is not idempotent at the
    application level but is safe to retry on 5xx because the server has not
    issued a token yet).

    The default timeout is stashed on the session under a private attribute;
    callers should either pass ``timeout=`` explicitly on every call or use
    :func:`request_with_default_timeout`.
    """
    retry = Retry(
        total=int(total_retries),
        backoff_factor=float(backoff_factor),
        status_forcelist=tuple(status_forcelist),
        allowed_methods=frozenset({"HEAD", "GET", "PUT", "DELETE", "OPTIONS", "POST"}),
        raise_on_status=False,
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session = requests.Session()
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    setattr(session, _TIMEOUT_ATTR, float(timeout_default))
    return session


def request_with_default_timeout(
    session: requests.Session,
    method: str,
    url: str,
    **kwargs: Any,
) -> requests.Response:
    """Issue a request via *session*, applying the session's default timeout.

    If the caller passed ``timeout=`` explicitly, that value wins. Otherwise
    the session's stashed default (set by :func:`make_session`) is used. This
    keeps timeout policy in one place without forcing every call site to
    repeat the value.
    """
    if "timeout" not in kwargs:
        default = getattr(session, _TIMEOUT_ATTR, DEFAULT_TIMEOUT_SECONDS)
        kwargs["timeout"] = default
    return session.request(method, url, **kwargs)
