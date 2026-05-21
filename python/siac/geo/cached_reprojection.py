"""Content-addressed disk cache for ``reproject_match`` outputs.

The wave-17 cProfile breakdown (``/tmp/wave17_profile.txt``) put 35 s of
T33KWP's 528 s wall-clock inside ``warp.py:reproject`` — mostly the MCD43
monthly-composite reprojection from the BRDF cache's native grid onto
the scene's solver grid. Source data (raw MCD43 HDF tiles, CAMS NetCDFs)
is already disk-cached by the upstream providers, but the **reprojected
output** is recomputed every run.

This module persists the reprojected ``xr.DataArray`` under a SHA-256
content-hash key. The key folds in:

- The :class:`siac.geo.target_grid.TargetGrid` signature (bounds, CRS,
  resolution, shape) — guarantees that a cache built for one scene
  grid can't be served to a different one.
- A caller-supplied ``source_identity`` string (e.g.
  ``f"mcd43-{tile_id}-{date}-{kernel_name}"``). The caller knows the
  provenance better than we can deduce it from the array.
- The resampling method.
- A short ``CACHE_FORMAT_VERSION`` so we can invalidate cleanly when
  the on-disk layout or our reproject contract changes.

Persistence uses NetCDF (engine=netcdf4 if available, else h5netcdf,
else scipy). NetCDF preserves dims/coords + CRS via rioxarray's
``spatial_ref`` coordinate variable. Atomic writes via tempfile + rename
keep the cache safe under interrupted runs.
"""

from __future__ import annotations

import hashlib
import logging
import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import xarray as xr

from siac.geo.reprojection import reproject_match
from siac.geo.target_grid import TargetGrid

if TYPE_CHECKING:
    from rasterio.enums import Resampling

logger = logging.getLogger(__name__)

#: Bumped on any change to the on-disk format or the cache-key
#: derivation so old entries become unreadable rather than silently
#: misread. The cloud cache uses ``cloud-v1``; we use a different
#: namespace so the two caches can share a directory without colliding.
CACHE_FORMAT_VERSION: str = "reproject-v1"


def _stable_resampling_repr(resampling: str | Resampling | None) -> str:
    """Deterministic string for the resampling enum/string/None."""
    if resampling is None:
        return "auto"
    if isinstance(resampling, str):
        return resampling.lower()
    # rasterio.enums.Resampling has a .name attribute.
    name = getattr(resampling, "name", None)
    if name is not None:
        return str(name).lower()
    return repr(resampling)


def compute_cache_key(
    *,
    target: TargetGrid,
    source_identity: str,
    resampling: str | Resampling | None,
    extra_namespace: str = "",
) -> str:
    """Compute the deterministic cache key for this reprojection request."""
    h = hashlib.sha256()
    h.update(CACHE_FORMAT_VERSION.encode("utf-8"))
    h.update(b"\0")
    h.update(extra_namespace.encode("utf-8"))
    h.update(b"\0")
    for chunk in target.signature():
        h.update(repr(chunk).encode("utf-8"))
        h.update(b"\0")
    h.update(source_identity.encode("utf-8"))
    h.update(b"\0")
    h.update(_stable_resampling_repr(resampling).encode("utf-8"))
    return h.hexdigest()


def _cache_path(cache_dir: Path, key: str) -> Path:
    # 2-char shard, same convention as the cloud cache.
    return Path(cache_dir) / key[:2] / f"{key}.nc"


def _restore_crs(da: xr.DataArray, target: TargetGrid) -> xr.DataArray:
    """Ensure rioxarray CRS metadata is set on a freshly-loaded DataArray."""
    import rioxarray  # noqa: F401  # registers .rio accessor

    try:
        if da.rio.crs is None:
            da = da.rio.write_crs(target.crs)
    except Exception:
        # Some xarray builds raise on missing crs; just write it.
        da = da.rio.write_crs(target.crs)
    return da


def load_cached_reprojection(
    cache_dir: Path | str | None,
    key: str,
    *,
    target: TargetGrid,
) -> xr.DataArray | None:
    """Return a cached reprojection if present + shape-compatible, else ``None``.

    Returns None on any load error (missing, corrupted, shape mismatch)
    rather than raising — the caller will fall through to the live
    ``reproject_match`` and overwrite the cache on success.
    """
    if cache_dir is None:
        return None
    path = _cache_path(Path(cache_dir), key)
    if not path.exists():
        return None
    try:
        # ``decode_coords="all"`` recovers rioxarray's spatial_ref coord
        # so .rio.crs round-trips through the NetCDF.
        ds = xr.open_dataset(path, decode_coords="all", engine=None)
    except (OSError, ValueError, KeyError) as exc:
        logger.warning(
            "Reproject cache hit at %s failed to load (%s: %s); recomputing.",
            path,
            type(exc).__name__,
            exc,
        )
        return None
    try:
        # We store the array under the var name "data" — see the save
        # path below. If the file uses an older layout fall through.
        if "data" not in ds.data_vars:
            logger.warning(
                "Reproject cache file %s has unexpected layout (vars=%s); recomputing.",
                path,
                list(ds.data_vars),
            )
            return None
        out = ds["data"].load()
        # ``ds.close`` after loading; the array is now in-memory.
        ds.close()
        if out.shape[-2:] != target.shape:
            logger.info(
                "Reproject cache shape mismatch at %s (%s vs expected %s); recomputing.",
                path,
                out.shape[-2:],
                target.shape,
            )
            return None
        out = _restore_crs(out, target)
        logger.info("Reproject cache HIT: %s", path)
        return out
    finally:
        import contextlib

        with contextlib.suppress(Exception):
            ds.close()


def save_cached_reprojection(
    cache_dir: Path | str | None,
    key: str,
    array: xr.DataArray,
) -> None:
    """Persist ``array`` under ``key``. Best-effort: log + continue on failure."""
    if cache_dir is None:
        return
    path = _cache_path(Path(cache_dir), key)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{key}.", suffix=".nc", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        named = array.rename("data") if array.name != "data" else array
        # zlib level 4 is the sweet-spot for atmospheric-correction
        # floats: ~5× compression, ~25% CPU overhead vs uncompressed.
        ds = named.to_dataset(name="data")
        encoding = {
            "data": {
                "zlib": True,
                "complevel": 4,
                "shuffle": True,
            }
        }
        ds.to_netcdf(tmp_path, encoding=encoding)
        tmp_path.replace(path)
        logger.info(
            "Reproject cache MISS → wrote %s (%.1f MB)",
            path,
            path.stat().st_size / 1e6,
        )
    except (OSError, ValueError) as exc:
        logger.warning(
            "Failed to write reproject cache %s (%s: %s); proceeding.",
            path,
            type(exc).__name__,
            exc,
        )
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass


def cached_reproject_match(
    source: xr.DataArray,
    target: xr.DataArray | TargetGrid,
    *,
    source_identity: str,
    cache_dir: Path | str | None = None,
    resampling: str | Resampling | None = None,
    extra_namespace: str = "",
    nodata: float | None = None,
) -> xr.DataArray:
    """Reproject ``source`` to match ``target`` with disk caching.

    Functionally equivalent to :func:`siac.geo.reprojection.reproject_match`
    when ``cache_dir`` is ``None`` — same inputs and outputs. When
    ``cache_dir`` is set, the output is persisted under a content-hash
    key and reused on subsequent calls with matching inputs.

    Args:
        source: The DataArray to reproject.
        target: Either a template :class:`xr.DataArray` (legacy entry
            point) or an explicit :class:`TargetGrid`. The latter is
            preferred for new code — it makes the cache key derivation
            explicit instead of inferring it from a transient template.
        source_identity: A short, stable string that uniquely identifies
            ``source``'s provenance (e.g. ``"mcd43-T33KWP-2026-03-29-fiso"``).
            Two callers using different source data MUST pass different
            identities; otherwise the cache will serve one's data to the
            other. **No fallback to content-hashing** — the helper raises
            ``ValueError`` if you pass an empty identity, because silent
            content-hashing would defeat the cache's whole purpose for
            very large priors (gigabytes of source data hashed per call).
        cache_dir: Directory to persist the cache under. ``None`` disables
            caching entirely (falls through to plain ``reproject_match``).
        resampling: Resampling method. ``None`` lets ``reproject_match``
            pick based on dtype.
        extra_namespace: Optional cache key salt — set to a non-empty
            string to partition the cache between unrelated consumers
            sharing the same ``cache_dir``.
        nodata: Optional nodata value written into the source before
            reprojection (forwarded to ``reproject_match``).

    Returns:
        The reprojected DataArray, semantically identical to
        ``reproject_match(source, target, resampling=resampling, nodata=nodata)``.
    """
    if not source_identity:
        raise ValueError(
            "cached_reproject_match requires a non-empty source_identity "
            "to keep the cache key meaningful. Pass a string that uniquely "
            "identifies the provenance of the source data."
        )

    if isinstance(target, TargetGrid):
        target_grid = target
        template_for_reproject = target.as_template_da()
    else:
        target_grid = TargetGrid.from_template(target)
        template_for_reproject = target

    if cache_dir is None:
        return reproject_match(
            source,
            template_for_reproject,
            resampling=resampling,
            nodata=nodata,
        )

    key = compute_cache_key(
        target=target_grid,
        source_identity=source_identity,
        resampling=resampling,
        extra_namespace=extra_namespace,
    )
    cached = load_cached_reprojection(cache_dir, key, target=target_grid)
    if cached is not None:
        return cached

    result = reproject_match(
        source,
        template_for_reproject,
        resampling=resampling,
        nodata=nodata,
    )
    save_cached_reprojection(cache_dir, key, result)
    return result
