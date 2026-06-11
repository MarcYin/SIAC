"""Shared helpers for the from-source RT engine build harnesses (6S, libRadtran).

Both engines are fetched as upstream archives, integrity-checked, and compiled
against the active conda/pixi toolchain. The download/verify mechanics and the
conda library-path plumbing are identical concerns, so they live here once:

- :func:`fetch_archive` — atomic streaming download (pid-unique temp + rename,
  size cap) with SHA-256 verification. A checksum mismatch removes the
  poisoned/corrupt archive so the next run re-downloads, and the expected hash
  can be overridden per-run via an environment variable.
- :func:`prepend_conda_lib_path` — prepend ``$CONDA_PREFIX/lib`` to the
  platform's runtime library search path (``DYLD_FALLBACK_LIBRARY_PATH`` on
  macOS, ``LD_LIBRARY_PATH`` elsewhere) so engine binaries built against conda
  NetCDF/GSL resolve their shared libraries at run time as well as build time.
"""

from __future__ import annotations

import contextlib
import hashlib
import logging
import os
import sys
from pathlib import Path

import requests

logger = logging.getLogger(__name__)


def archive_sha256(path: Path) -> str:
    """SHA-256 of ``path`` in fixed-size chunks (never loads the whole archive)."""
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_archive(
    url: str,
    archive_path: Path,
    *,
    expected_sha256: str | None,
    sha_env_var: str,
    max_bytes: int,
    timeout_s: float = 900.0,
    what: str = "archive",
) -> None:
    """Download ``url`` to ``archive_path`` atomically and verify its SHA-256.

    ``sha_env_var`` names an environment variable that overrides
    ``expected_sha256`` (for intentional upstream changes); when neither is
    set the check is skipped with a warning. ``what`` flavours log/error text
    (e.g. ``"libRadtran"``, ``"6SV2.1"``).
    """
    if not archive_path.exists():
        logger.info("Downloading %s archive %s -> %s", what, url, archive_path.name)
        # Unique temp name so concurrent builds don't clobber each other's
        # partial download; rename publishes only a complete file.
        tmp_path = archive_path.with_suffix(archive_path.suffix + f".{os.getpid()}.tmp")
        try:
            written = 0
            with requests.get(url, timeout=timeout_s, stream=True) as response:
                response.raise_for_status()
                with tmp_path.open("wb") as fh:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if not chunk:
                            continue
                        written += len(chunk)
                        if written > max_bytes:
                            raise RuntimeError(
                                f"{what} download {url} exceeded the "
                                f"{max_bytes // 1024**3} GB cap; aborting."
                            )
                        fh.write(chunk)
            tmp_path.replace(archive_path)
        except BaseException:
            with contextlib.suppress(OSError):
                tmp_path.unlink()
            raise

    expected = os.environ.get(sha_env_var) or expected_sha256
    if expected:
        actual = archive_sha256(archive_path)
        if actual.lower() != expected.lower():
            # Remove the poisoned/corrupt archive so the next run re-downloads it.
            with contextlib.suppress(OSError):
                archive_path.unlink()
            raise RuntimeError(
                f"{what} archive from {url} has SHA-256 {actual!r} but expected "
                f"{expected!r}; the cached file has been removed. Re-run to re-download, "
                f"or set {sha_env_var} to the new hash if the change is intentional."
            )
    else:
        logger.warning(
            "%s archive checksum not configured (%s); proceeding without integrity "
            "verification of %s.",
            what,
            sha_env_var,
            archive_path.name,
        )


def conda_prefix() -> str:
    """Active conda/pixi environment prefix (falls back to ``sys.prefix``)."""
    return os.environ.get("CONDA_PREFIX") or sys.prefix


def prepend_conda_lib_path(env: dict[str, str]) -> dict[str, str]:
    """Prepend ``<conda prefix>/lib`` to the runtime library search path in ``env``."""
    prefix = env.get("CONDA_PREFIX") or os.environ.get("CONDA_PREFIX")
    if prefix:
        lib = str(Path(prefix) / "lib")
        var = "DYLD_FALLBACK_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
        env[var] = lib + os.pathsep + env.get(var, "")
    return env
