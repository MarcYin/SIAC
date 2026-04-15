"""Logging filter that redacts secrets from log records."""

from __future__ import annotations

import logging
import re

# Patterns that look like tokens, keys, or passwords in log messages.
_SECRET_PATTERNS = (
    re.compile(r"(Bearer\s+)[A-Za-z0-9._\-]+", re.IGNORECASE),
    re.compile(r"(access_key_id[=:]\s*)[A-Za-z0-9/+]+", re.IGNORECASE),
    re.compile(r"(secret_access_key[=:]\s*)[A-Za-z0-9/+=]+", re.IGNORECASE),
    re.compile(r"(password[=:]\s*)\S+", re.IGNORECASE),
    re.compile(r"(token[=:]\s*)[A-Za-z0-9._\-]+", re.IGNORECASE),
)


class SecretRedactionFilter(logging.Filter):
    """Scrub credential-like patterns from log messages before emission."""

    def filter(self, record: logging.LogRecord) -> bool:
        if isinstance(record.msg, str):
            record.msg = _redact(record.msg)
        if record.args:
            if isinstance(record.args, dict):
                record.args = {k: _redact(v) if isinstance(v, str) else v for k, v in record.args.items()}
            elif isinstance(record.args, tuple):
                record.args = tuple(_redact(a) if isinstance(a, str) else a for a in record.args)
        return True


def _redact(text: str) -> str:
    for pattern in _SECRET_PATTERNS:
        text = pattern.sub(lambda m: m.group(1) + "***REDACTED***", text)
    return text
