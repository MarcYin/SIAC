"""Tests for the retry utility."""

import pytest

from siac.adapters._retry import retry_transient


def test_retry_succeeds_first_attempt():
    """Function that succeeds on first call should return immediately."""
    calls = [0]

    def succeed():
        calls[0] += 1
        return 42

    result = retry_transient(succeed, max_attempts=3, base_delay_s=0.01)
    assert result == 42
    assert calls[0] == 1


def test_retry_succeeds_after_transient_failure():
    """Function that fails once then succeeds should be retried."""
    calls = [0]

    def fail_once():
        calls[0] += 1
        if calls[0] == 1:
            raise ConnectionError("transient")
        return "ok"

    result = retry_transient(fail_once, max_attempts=3, base_delay_s=0.01)
    assert result == "ok"
    assert calls[0] == 2


def test_retry_exhausted_raises_last_exception():
    """All attempts failing should raise the last exception."""
    calls = [0]

    def always_fail():
        calls[0] += 1
        raise TimeoutError(f"attempt {calls[0]}")

    with pytest.raises(TimeoutError, match="attempt 3"):
        retry_transient(always_fail, max_attempts=3, base_delay_s=0.01)
    assert calls[0] == 3


def test_retry_does_not_catch_non_transient():
    """Non-transient exceptions (e.g. ValueError) should not be retried."""
    calls = [0]

    def bad_input():
        calls[0] += 1
        raise ValueError("permanent error")

    with pytest.raises(ValueError, match="permanent"):
        retry_transient(bad_input, max_attempts=3, base_delay_s=0.01)
    assert calls[0] == 1


def test_retry_with_custom_transient_exceptions():
    """Custom transient exception types should be retried."""
    calls = [0]

    def custom_fail():
        calls[0] += 1
        if calls[0] < 3:
            raise KeyError("missing")
        return "found"

    result = retry_transient(
        custom_fail,
        max_attempts=3,
        base_delay_s=0.01,
        transient_exceptions=(KeyError,),
    )
    assert result == "found"
    assert calls[0] == 3


def test_retry_passes_args_and_kwargs():
    """Arguments should be forwarded to the target function."""

    def add(a, b, c=0):
        return a + b + c

    result = retry_transient(add, 1, 2, c=3, max_attempts=1, base_delay_s=0.01)
    assert result == 6
