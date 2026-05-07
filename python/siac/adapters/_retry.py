"""Lightweight retry helper for transient I/O failures."""

from __future__ import annotations

import logging
import random
import time
from typing import Any, TypeVar

logger = logging.getLogger(__name__)

T = TypeVar("T")

# REVIEW.md §3.3 _retry.py: ``OSError`` is broader than we'd like (it covers
# ENOSPC, EROFS, ENOENT) but several upstream stacks raise plain ``OSError``
# for transient conditions (rasterio with curl, fsspec); narrowing this would
# need a per-callsite audit. Kept here with a note for the next pass.
_TRANSIENT_EXCEPTIONS: tuple[type[BaseException], ...] = (
    OSError,
    ConnectionError,
    TimeoutError,
)

# Cap individual sleeps so a long retry chain stays bounded.
_MAX_BACKOFF_S = 30.0


def retry_transient(
    fn: Any,
    *args: Any,
    max_attempts: int = 3,
    base_delay_s: float = 2.0,
    transient_exceptions: tuple[type[BaseException], ...] = _TRANSIENT_EXCEPTIONS,
    label: str = "",
    **kwargs: Any,
) -> Any:
    """Call *fn* with jittered exponential backoff on transient exceptions.

    Backoff is ``base_delay_s * 2**(attempt - 1)`` multiplied by a uniform
    factor in ``[0.5, 1.5)`` (full jitter), capped at ``_MAX_BACKOFF_S``. The
    jitter avoids thundering-herd collisions when many workers retry in lock
    step (REVIEW.md §2.6).

    Args:
        fn: Callable to invoke.
        max_attempts: Total number of attempts (including the first).
        base_delay_s: Initial delay in seconds; doubled on each retry.
        transient_exceptions: Exception types considered transient.
        label: Human-readable label for log messages.

    Returns:
        The return value of *fn*.

    Raises:
        The last exception if all attempts fail.
    """
    last_exc: BaseException | None = None
    for attempt in range(1, max_attempts + 1):
        try:
            return fn(*args, **kwargs)
        except transient_exceptions as exc:
            last_exc = exc
            if attempt >= max_attempts:
                break
            raw_delay = base_delay_s * (2 ** (attempt - 1))
            jitter = 0.5 + random.random()  # [0.5, 1.5) full-jitter
            delay = min(raw_delay * jitter, _MAX_BACKOFF_S)
            logger.warning(
                "%s: attempt %d/%d failed (%s: %s); retrying in %.1fs",
                label or fn.__name__,
                attempt,
                max_attempts,
                type(exc).__name__,
                exc,
                delay,
            )
            time.sleep(delay)

    raise last_exc  # type: ignore[misc]
