#!/usr/bin/env python3
"""Re-encode the Cloud Score+ raster cache from scratch onto the group workspace.

The cache is what allows the winner index to be re-derived under a different
clear threshold or ordering without re-spending the Earth Engine requests that
produced it. It lives on scratch, which is purged, so keeping that option means
moving it -- and at 663 GB of float64 it is worth re-encoding on the way.

Quantising to uint8 leaves the winner field bit-identical: the per-day weight
term separates candidates far more coarsely than 1/254, so the argmax never
sees the difference. Measured over eight scenes, every pixel agreed.

One caveat on what this preserves. The cache holds only the ~111 shortlisted
acquisitions per scene, not all ~463 surveyed, so it supports retuning the
threshold or the ordering *within* that shortlist. A retune that widens the
shortlist needs a fresh fetch regardless.
"""

from __future__ import annotations

import argparse
import json
import shutil
import time
from pathlib import Path
from typing import Any

import numpy as np
from tools.aeronet_validation.build_cloudscore_index_local import (
    CLEAR_SCORE_NODATA,
    CLEAR_SCORE_SCALE,
    FILE_PREFIX,
)


def encode(values: np.ndarray) -> np.ndarray:
    """Clear scores in [0, 1] to uint8, reserving 255 for nodata."""

    source = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(source)
    out = np.full(source.shape, CLEAR_SCORE_NODATA, dtype=np.uint8)
    scaled = np.round(np.clip(source[finite], 0.0, 1.0) * CLEAR_SCORE_SCALE)
    out[finite] = scaled.astype(np.uint8)
    return out


def decode(values: np.ndarray) -> np.ndarray:
    """Inverse of :func:`encode`, for verification."""

    source = np.asarray(values)
    out = source.astype(np.float64)
    return np.where(source == CLEAR_SCORE_NODATA, np.nan, out / CLEAR_SCORE_SCALE)


def reencode_raster(source: Path, target: Path) -> tuple[int, int, float]:
    """Re-encode one raster, returning ``(source bytes, target bytes, max error)``."""

    import rasterio

    with rasterio.open(source) as handle:
        values = handle.read(1)
        profile = handle.profile
    encoded = encode(values)
    profile.update(
        dtype="uint8",
        nodata=CLEAR_SCORE_NODATA,
        compress="deflate",
        zlevel=6,
        predictor=2,
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(".tif.tmp")
    with rasterio.open(temporary, "w", **profile) as handle:
        handle.write(encoded, 1)
    temporary.replace(target)

    finite = np.isfinite(values)
    error = float(np.max(np.abs(decode(encoded)[finite] - values[finite]))) if finite.any() else 0.0
    return source.stat().st_size, target.stat().st_size, error


def archive_scene(source: Path, target: Path, *, overwrite: bool = False) -> dict[str, Any]:
    """Re-encode one scene's rasters and copy its manifests verbatim."""

    rasters = sorted(source.glob(f"images/*/{FILE_PREFIX}_*.tif"))
    source_bytes = target_bytes = 0
    worst = 0.0
    written = skipped = 0
    for raster in rasters:
        destination = target / raster.relative_to(source)
        if destination.is_file() and not overwrite:
            skipped += 1
            target_bytes += destination.stat().st_size
            source_bytes += raster.stat().st_size
            continue
        before, after, error = reencode_raster(raster, destination)
        source_bytes += before
        target_bytes += after
        worst = max(worst, error)
        written += 1
    # Manifests carry the pixel_api provenance and the exact request config.
    manifests = source / "manifests"
    if manifests.is_dir():
        (target / "manifests").mkdir(parents=True, exist_ok=True)
        for manifest in manifests.glob("*.json"):
            destination = target / "manifests" / manifest.name
            if overwrite or not destination.is_file():
                shutil.copy2(manifest, destination)
    return {
        "rasters": len(rasters),
        "written": written,
        "skipped": skipped,
        "source_bytes": source_bytes,
        "target_bytes": target_bytes,
        "max_error": worst,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    source_root = Path(args.source)
    target_root = Path(args.target)
    scenes = sorted(path for path in source_root.iterdir() if path.is_dir())
    if args.shard_count > 1:
        scenes = [s for i, s in enumerate(scenes) if i % args.shard_count == args.shard]
    if args.limit:
        scenes = scenes[: int(args.limit)]

    totals = {"rasters": 0, "written": 0, "skipped": 0, "source_bytes": 0, "target_bytes": 0}
    worst = 0.0
    failures: dict[str, str] = {}
    started = time.perf_counter()
    for position, scene in enumerate(scenes):
        try:
            result = archive_scene(scene, target_root / scene.name, overwrite=args.overwrite)
        except Exception as error:  # noqa: BLE001 - one bad scene must not end the shard
            failures[scene.name] = f"{type(error).__name__}: {error}"
            continue
        for key in totals:
            totals[key] += result[key]
        worst = max(worst, result["max_error"])
        if (position + 1) % 50 == 0:
            print(
                f"  {position + 1}/{len(scenes)} scenes, "
                f"{totals['target_bytes'] / 1e9:.1f} GB written",
                flush=True,
            )

    summary = {
        "scenes": len(scenes),
        "failed": len(failures),
        "failures": dict(list(failures.items())[:10]),
        **totals,
        "compression": (
            totals["source_bytes"] / totals["target_bytes"] if totals["target_bytes"] else 0.0
        ),
        # Quantisation error must not exceed half a step, or the encoding is wrong.
        "max_error": worst,
        "max_error_budget": 0.5 / CLEAR_SCORE_SCALE,
        "wall_seconds": round(time.perf_counter() - started, 1),
    }
    if args.summary:
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        Path(args.summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--source", required=True)
    value.add_argument("--target", required=True)
    value.add_argument("--overwrite", action="store_true")
    value.add_argument("--limit", type=int, default=0)
    value.add_argument("--shard", type=int, default=0)
    value.add_argument("--shard-count", type=int, default=1)
    value.add_argument("--summary")
    return value


def main() -> None:
    print(json.dumps(run(parser().parse_args()), indent=2))


if __name__ == "__main__":
    main()
