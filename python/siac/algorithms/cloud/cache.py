"""Content-addressed disk cache for cloud-classification model outputs.

The OmniCloudMask PyTorch backend takes ~20-25 s per Sentinel-2 scene and
is non-deterministic across runs (small set of edge-of-cloud pixels flip
between adjacent classes — this is what wave 17 had to absorb via a
loosened AOT regression tolerance until we tightened it). Caching the
model output on disk under a content-addressed key gives us both
properties at once:

- Subsequent runs on the same TOA inputs skip the PyTorch inference
  entirely (single-digit-ms cache hit vs ~25 s miss).
- Run-to-run output becomes deterministic because every consumer (the
  regression suite, repeated production runs, dev iteration) reads the
  same persisted bytes.

The cache key is a SHA-256 digest of:

- The 3 raw input rasters (red/green/nir) coerced to float32 and
  serialised in C order.
- The class-mapping dict (sorted, JSON-encoded).
- The target resolution (the rasters are already at that resolution
  when this is called, but the parameter influences upstream
  resampling behaviour so it's included for safety).
- A short ``CACHE_FORMAT_VERSION`` constant so we can invalidate the
  cache when the OCM model or our normalisation logic changes.

The cache itself is a NumPy ``.npz`` file holding the standardised
``cloud_classes`` raster (``uint8``, shape == red/green/nir shape) plus
its dims/coords as JSON metadata. We don't use NetCDF because xarray's
write path involves dask/netcdf4 wiring that adds overhead larger than
the cache save itself; ``.npz`` is plain numpy and reads in <10 ms.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr

if TYPE_CHECKING:
    from collections.abc import Iterable

logger = logging.getLogger(__name__)

#: Bumped whenever the on-disk format changes OR the wrapped model
#: contract changes such that historic cache entries must not be served.
CACHE_FORMAT_VERSION: str = "v1"


def _stable_class_mapping_repr(
    class_mapping: dict[int, Iterable[int]] | None,
) -> str:
    """Return a deterministic string representation of the class mapping."""
    if class_mapping is None:
        return "null"
    items = sorted(
        (int(k), sorted(int(x) for x in v)) for k, v in class_mapping.items()
    )
    return json.dumps(items, separators=(",", ":"))


def compute_cache_key(
    *,
    red: xr.DataArray,
    green: xr.DataArray,
    nir: xr.DataArray,
    class_mapping: dict[int, Iterable[int]] | None,
    target_resolution_m: float,
    extra_namespace: str = "",
) -> str:
    """Compute the deterministic content-hash cache key for these inputs."""
    h = hashlib.sha256()
    h.update(CACHE_FORMAT_VERSION.encode("utf-8"))
    h.update(b"\0")
    h.update(extra_namespace.encode("utf-8"))
    h.update(b"\0")
    h.update(f"target_resolution_m={float(target_resolution_m):.6f}".encode("utf-8"))
    h.update(b"\0")
    h.update(_stable_class_mapping_repr(class_mapping).encode("utf-8"))
    for label, arr in (("red", red), ("green", green), ("nir", nir)):
        h.update(b"\0")
        h.update(label.encode("utf-8"))
        h.update(b"\0")
        h.update(repr(tuple(arr.shape)).encode("utf-8"))
        # ``np.ascontiguousarray`` guarantees byte order is well-defined.
        values = np.ascontiguousarray(np.asarray(arr.values, dtype=np.float32))
        h.update(values.tobytes())
    return h.hexdigest()


def _cache_path(cache_dir: Path, key: str) -> Path:
    # Shard by first 2 hex chars so a single directory never grows past
    # filesystem-unfriendly sizes (256-way fan-out, suffix is unique).
    return Path(cache_dir) / key[:2] / f"{key}.npz"


def load_cached_cloud_classes(
    cache_dir: Path | str | None,
    key: str,
    *,
    template: xr.DataArray,
) -> xr.DataArray | None:
    """Return cached classes if present and shape-compatible, else ``None``.

    A shape mismatch returns ``None`` rather than raising — it means an
    upstream input changed in a way the key didn't catch, so we should
    just miss the cache and recompute.
    """
    if cache_dir is None:
        return None
    path = _cache_path(Path(cache_dir), key)
    if not path.exists():
        return None
    try:
        with np.load(path, allow_pickle=False) as data:
            values = np.asarray(data["values"], dtype=np.uint8)
    except (OSError, ValueError, KeyError) as exc:
        logger.warning(
            "Cloud-mask cache hit at %s failed to load (%s: %s); recomputing.",
            path,
            type(exc).__name__,
            exc,
        )
        return None
    if values.shape != template.shape:
        logger.info(
            "Cloud-mask cache shape mismatch at %s (%s vs expected %s); "
            "recomputing.",
            path,
            values.shape,
            template.shape,
        )
        return None
    out = xr.DataArray(values, dims=template.dims, coords=template.coords)
    out.name = "cloud_classes"
    logger.info("Cloud-mask cache HIT: %s", path)
    return out


def save_cached_cloud_classes(
    cache_dir: Path | str | None,
    key: str,
    classes: xr.DataArray,
) -> None:
    """Persist classes under the given key. Best-effort: failures are logged."""
    if cache_dir is None:
        return
    path = _cache_path(Path(cache_dir), key)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Stage to a sibling tempfile alongside ``path`` and ``os.replace`` on
    # success — same atomic-rename pattern used for COG outputs in
    # ``storage/writers.py``. ``np.savez`` auto-appends ``.npz`` if the
    # destination doesn't already end with that suffix, so we build the
    # tempfile name to already include it and capture the exact resulting
    # path from numpy via tempfile's deterministic name.
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{key}.", suffix=".npz", dir=str(path.parent)
    )
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        values = np.ascontiguousarray(np.asarray(classes.values, dtype=np.uint8))
        np.savez(tmp_path, values=values)
        os.replace(tmp_path, path)
        logger.info("Cloud-mask cache MISS → wrote %s", path)
    except (OSError, ValueError) as exc:
        logger.warning(
            "Failed to write cloud-mask cache %s (%s: %s); proceeding.",
            path,
            type(exc).__name__,
            exc,
        )
        # Clean up the staged file if the save failed mid-stream.
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass


def maybe_run_with_cache(
    *,
    cache_dir: Path | str | None,
    red: xr.DataArray,
    green: xr.DataArray,
    nir: xr.DataArray,
    class_mapping: dict[int, Iterable[int]] | None,
    target_resolution_m: float,
    compute_fn: Any,
    extra_namespace: str = "",
) -> xr.DataArray:
    """Run ``compute_fn`` with disk-cache short-circuit when ``cache_dir`` set.

    ``compute_fn`` is invoked with no arguments and must return an
    ``xr.DataArray`` of standardised cloud classes (dtype uint8).
    """
    if cache_dir is None:
        return compute_fn()
    key = compute_cache_key(
        red=red,
        green=green,
        nir=nir,
        class_mapping=class_mapping,
        target_resolution_m=target_resolution_m,
        extra_namespace=extra_namespace,
    )
    cached = load_cached_cloud_classes(cache_dir, key, template=red)
    if cached is not None:
        return cached
    classes = compute_fn()
    save_cached_cloud_classes(cache_dir, key, classes)
    return classes
