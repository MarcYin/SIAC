"""Lightweight retry helper for transient I/O failures."""

from __future__ import annotations

import logging
import time
from typing import Any, TypeVar

logger = logging.getLogger(__name__)

T = TypeVar("T")

_TRANSIENT_EXCEPTIONS: tuple[type[BaseException], ...] = (
    OSError,
    ConnectionError,
    TimeoutError,
)


def retry_transient(
    fn: Any,
    *args: Any,
    max_attempts: int = 3,
    base_delay_s: float = 2.0,
    transient_exceptions: tuple[type[BaseException], ...] = _TRANSIENT_EXCEPTIONS,
    label: str = "",
    **kwargs: Any,
) -> Any:
    """Call *fn* with exponential backoff on transient exceptions.

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
    delay = base_delay_s
    for attempt in range(1, max_attempts + 1):
        try:
            return fn(*args, **kwargs)
        except transient_exceptions as exc:
            last_exc = exc
            if attempt >= max_attempts:
                break
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
            delay *= 2.0

    raise last_exc  # type: ignore[misc]
