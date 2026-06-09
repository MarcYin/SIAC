"""Shared resolution of RT run-cache directories (libRadtran / 6S / remote LUT).

These backends persist RT computations on disk so a re-run reuses the work
instead of recomputing it: libRadtran ``uvspec`` run outputs, the 6S native
scene-LUT grid batch, and the remote-LUT materialised scene subset. They all
resolve their cache directory the same way, in priority order:

1. an explicit per-backend config directory (``configured``);
2. ``None`` (cache disabled) when ``$SIAC_RT_RUN_CACHE_DISABLED`` is truthy;
3. otherwise ``$SIAC_RT_RUN_CACHE_ROOT`` (or ``~/.cache/siac``) joined with the
   backend-specific ``subpath``.

The environment variables give operators a single global switch and let the
test-suite redirect every backend's cache into a throwaway directory, so tests
never touch the developer's real cache or see stale cross-test cache hits.
"""

from __future__ import annotations

import os
from pathlib import Path

_TRUTHY = {"1", "true", "yes", "on"}


def resolve_run_cache_dir(configured: str | Path | None, *, subpath: str) -> Path | None:
    """Resolve a backend run-cache directory, or ``None`` when disabled.

    ``configured`` is the per-backend config value (wins when set); ``subpath``
    is the backend-specific suffix under the cache root (e.g.
    ``"rt6s/run_cache"``).
    """
    if configured is not None:
        return Path(configured).expanduser()
    if os.environ.get("SIAC_RT_RUN_CACHE_DISABLED", "").strip().lower() in _TRUTHY:
        return None
    root_env = os.environ.get("SIAC_RT_RUN_CACHE_ROOT")
    root = Path(root_env).expanduser() if root_env else Path.home().expanduser() / ".cache" / "siac"
    return root / subpath
