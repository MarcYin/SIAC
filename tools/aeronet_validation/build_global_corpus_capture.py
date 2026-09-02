#!/usr/bin/env python3
"""Capture L1C tensors for every catalogue scene, sharded across a Slurm array.

``capture_mixed_resolution_l1c`` is written per matchup: one invocation, one
scene. At 4,849 scenes that would be as many process launches, each paying GDAL
and PROJ import costs to read a few windows, so this drives its ``run`` in
process and keeps the per-scene locking and resume behaviour it already has.

The scene catalogue doubles as capture's manifest -- it already carries
``matchup_id``, ``product_id``, ``longitude`` and ``latitude`` -- so there is no
conversion step to drift out of sync with the catalogue it came from.

No calibration templates exist for a global catalogue, so the grid is derived
from the catalogue point. That derivation was verified to agree exactly with
``build_global_corpus_index.grid_for_point`` -- same CRS, same origin, same
extent, across six UTM zones and both hemispheres -- because a disagreement of
one pixel would silently misalign the winner index against the reflectance it
indexes.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Any

from experiments.current_date_toa_prior.capture_mixed_resolution_l1c import (
    GRID_SCALE,
    _existing_quality_approved,
)
from experiments.current_date_toa_prior.capture_mixed_resolution_l1c import (
    build_parser as build_capture_parser,
)
from experiments.current_date_toa_prior.capture_mixed_resolution_l1c import run as capture_run

#: Matches the index grid: a power of two on the 60 m aerosol grid, 384 at 20 m.
TEMPLATE_SIZE = 128
TEMPLATE_RESOLUTION_M = 60.0


def capture_args(
    manifest: Path,
    matchup_id: str,
    output_root: Path,
    *,
    minimum_finite_fraction: float,
    scale: int,
) -> argparse.Namespace:
    """Build the namespace ``capture_mixed_resolution_l1c.run`` expects."""

    # Parsed through capture's own parser rather than hand-built, so every
    # option it declares arrives with its declared default. Hand-listing them
    # silently omitted minimum_source_fraction on the first attempt, which the
    # capture would have read as missing rather than as its 0.50 default.
    namespace = build_capture_parser().parse_args(
        [
            "--manifest",
            str(manifest),
            "--matchup-id",
            str(matchup_id),
            # Unused once the template is derived, but the parser requires it.
            "--calibration-root",
            str(manifest.parent),
            "--output-root",
            str(output_root),
            "--scale",
            str(int(scale)),
            "--minimum-finite-fraction",
            str(float(minimum_finite_fraction)),
            "--allow-derived-template",
            "--template-size",
            str(TEMPLATE_SIZE),
            "--template-resolution",
            str(TEMPLATE_RESOLUTION_M),
        ]
    )
    return namespace


def run(args: argparse.Namespace) -> dict[str, Any]:
    manifest = Path(args.scene_catalog)
    with manifest.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if args.shard_count > 1:
        rows = [row for index, row in enumerate(rows) if index % args.shard_count == args.shard]
    if args.limit:
        rows = rows[: int(args.limit)]

    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    captured = skipped = 0
    failures: dict[str, str] = {}
    started = time.perf_counter()
    for position, row in enumerate(rows):
        matchup = row["matchup_id"]
        namespace = capture_args(
            manifest,
            matchup,
            output_root,
            minimum_finite_fraction=args.minimum_finite_fraction,
            scale=args.scale,
        )
        try:
            # Reuses capture's own quality gate, so a rerun neither redoes work
            # nor silently accepts an archive that failed its coverage check.
            if _existing_quality_approved(namespace) is not None:
                skipped += 1
                continue
            capture_run(namespace)
            captured += 1
        except Exception as error:  # noqa: BLE001 - one bad scene must not end the shard
            failures[matchup] = f"{type(error).__name__}: {error}"
        if (position + 1) % 25 == 0:
            print(
                f"  {position + 1}/{len(rows)} scenes, {captured} captured, "
                f"{skipped} skipped, {len(failures)} failed",
                flush=True,
            )
    elapsed = time.perf_counter() - started

    summary = {
        "scenes": len(rows),
        "captured": captured,
        "skipped_existing": skipped,
        "failed": len(failures),
        "failures": dict(list(failures.items())[:10]),
        "wall_seconds": round(elapsed, 1),
        "seconds_per_capture": round(elapsed / captured, 2) if captured else 0.0,
    }
    if args.summary:
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        Path(args.summary).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--scene-catalog", required=True)
    value.add_argument("--output-root", required=True)
    # capture's own constant. It is the ratio between grids -- context is
    # ``scale`` times wider and detail is 60/``scale`` metres -- not a
    # resolution, and it must be a positive odd integer.
    value.add_argument("--scale", type=int, default=GRID_SCALE)
    value.add_argument("--minimum-finite-fraction", type=float, default=0.98)
    value.add_argument("--limit", type=int, default=0)
    value.add_argument("--shard", type=int, default=0)
    value.add_argument("--shard-count", type=int, default=1)
    value.add_argument("--summary")
    return value


def main() -> None:
    print(json.dumps(run(parser().parse_args()), indent=2))


if __name__ == "__main__":
    main()
