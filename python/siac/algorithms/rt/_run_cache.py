"""Shared RT run-cache plumbing (libRadtran / 6S / remote LUT).

These backends persist RT computations on disk so a re-run reuses the work
instead of recomputing it: libRadtran ``uvspec`` run outputs, the 6S native
scene-LUT grid batch, and the remote-LUT materialised scene subset. This module
holds the pieces they share:

- :func:`resolve_run_cache_dir` — one resolution policy for the cache
  directory, in priority order: per-backend config dir; ``None`` when disabled
  (per-backend flag or ``$SIAC_RT_RUN_CACHE_DISABLED``); otherwise
  ``$SIAC_RT_RUN_CACHE_ROOT`` (or ``~/.cache/siac``) joined with the
  backend-specific ``subpath``. The environment variables give operators a
  single global switch and let the test-suite redirect every backend's cache
  into a throwaway directory.
- :func:`load_cache_entry` / :func:`store_cache_entry` — one implementation of
  the cache-entry contract: a missing or corrupt entry is a miss (never an
  error), and writes publish atomically (unique temp file + rename) so a crash
  or concurrent writer can never leave a half-written entry that reads back as
  valid.
"""

from __future__ import annotations

import logging
import os
import threading
from pathlib import Path
from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    from collections.abc import Callable

logger = logging.getLogger(__name__)

_TRUTHY = {"1", "true", "yes", "on"}

_T = TypeVar("_T")


def resolve_run_cache_dir(
    configured: str | Path | None,
    *,
    subpath: str,
    enabled: bool = True,
) -> Path | None:
    """Resolve a backend run-cache directory, or ``None`` when disabled.

    ``configured`` is the per-backend config value (wins when set); ``subpath``
    is the backend-specific suffix under the cache root (e.g.
    ``"rt6s/run_cache"``); ``enabled=False`` (the per-backend config switch)
    disables the cache outright.
    """
    if not enabled:
        return None
    if configured is not None:
        return Path(configured).expanduser()
    if os.environ.get("SIAC_RT_RUN_CACHE_DISABLED", "").strip().lower() in _TRUTHY:
        return None
    root_env = os.environ.get("SIAC_RT_RUN_CACHE_ROOT")
    root = Path(root_env).expanduser() if root_env else Path.home().expanduser() / ".cache" / "siac"
    return root / subpath


def load_cache_entry(path: Path, reader: Callable[[Path], _T]) -> _T | None:
    """Read a cache entry via ``reader``; missing/corrupt entries are misses.

    A truncated or corrupt entry (crash mid-write, disk issue) must be
    recomputed rather than fail the run, so any reader exception is logged at
    debug level and treated as a miss.
    """
    if not path.exists():
        return None
    try:
        return reader(path)
    except Exception:
        logger.debug("Ignoring unreadable cache entry %s", path, exc_info=True)
        return None


def store_cache_entry(path: Path, writer: Callable[[Path], None]) -> None:
    """Publish a cache entry atomically; failures are logged and swallowed.

    ``writer`` receives a unique temp path next to ``path`` (pid + thread id in
    the name, so concurrent writers never collide) which is renamed over
    ``path`` only after a complete write. Cache writes are best-effort: an
    unwritable cache degrades to recomputation, never to a failed run.
    """
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f".{path.stem}.{os.getpid()}.{threading.get_ident()}{path.suffix}.tmp")
        writer(tmp)
        tmp.replace(path)
    except Exception:
        logger.debug("Failed to write cache entry %s", path, exc_info=True)
