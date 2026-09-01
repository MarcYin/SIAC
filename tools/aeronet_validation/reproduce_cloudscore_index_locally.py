#!/usr/bin/env python3
"""Reproduce a server-side Cloud Score+ winner index with a local mosaic.

The committed index runs ``qualityMosaic`` inside Earth Engine and downloads
only the winner plane. This driver holds the candidate selection and the daily
scalars fixed -- both are read from the reference index -- and varies only
*where the argmax happens*. If the winners agree, the local composer is a
faithful replica and the weighting becomes an offline, editable choice.

edown names each raster by the Earth Engine image id, which matches the ``idx``
field of the reference ``image_table`` exactly, so candidates join per image
rather than per date.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from tools.aeronet_validation.cloudscore_local_mosaic import (
    Candidate,
    agreement,
    compose_monthly_winners,
    erode_valid,
)

COLLECTION = "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED"
BAND = "cs_cdf"
FILE_PREFIX = "GOOGLE_CLOUD_SCORE_PLUS_V1_S2_HARMONIZED"


def read_reference(path: Path) -> dict[str, Any]:
    archive = np.load(path, allow_pickle=False)
    images = [json.loads(str(row)) for row in archive["image_table"]]
    scalars = {json.loads(str(row))["day"]: json.loads(str(row)) for row in archive["day_scalars"]}
    policy = json.loads(str(archive["index_policy_json"]))
    return {
        "months": [str(m) for m in archive["months"]],
        "winners": np.asarray(archive["winners"]),
        "images": images,
        "day_scalars": scalars,
        "policy": policy,
    }


def fetch_scores(
    bbox: tuple[float, float, float, float], days: list[str], output_root: Path, edown: str
) -> None:
    """Fetch cs_cdf for the given days. edown resumes, so this is idempotent."""
    output_root.mkdir(parents=True, exist_ok=True)
    for day in sorted(set(days)):
        command = [
            edown,
            "download",
            "--collection-id",
            COLLECTION,
            "--band",
            BAND,
            "--bbox",
            *[str(c) for c in bbox],
            "--start-date",
            day,
            "--end-date",
            day,
            "--output-root",
            str(output_root),
        ]
        result = subprocess.run(command, capture_output=True, text=True, timeout=900)
        if result.returncode != 0:
            raise RuntimeError(f"edown failed for {day}: {result.stderr[-400:]}")


def load_plane(
    output_root: Path, image_id: str, template: Any, resampling: str = "nearest"
) -> np.ndarray | None:
    """Read one image's cs_cdf onto the scene grid, or None if absent."""
    import rioxarray  # noqa: F401
    from rasterio.enums import Resampling

    matches = list(output_root.glob(f"images/*/{FILE_PREFIX}_{image_id}.tif"))
    if not matches:
        return None
    import rioxarray as rxr

    raster = rxr.open_rasterio(matches[0]).squeeze("band", drop=True)
    method = getattr(Resampling, resampling)
    return raster.rio.reproject_match(template, resampling=method).values.astype(np.float64)


def build_candidates(
    reference: dict[str, Any],
    output_root: Path,
    template: Any,
    erosion_px: int,
    gate_threshold: bool = False,
    resampling: str = "nearest",
) -> dict[str, list[Candidate]]:
    by_month: dict[str, list[Candidate]] = defaultdict(list)
    for record in reference["images"]:
        day = str(record["day"])
        plane = load_plane(output_root, str(record["idx"]), template, resampling)
        if plane is None:
            continue
        # ``(cs_cdf - threshold) / span`` is a continuous ordering, not a mask:
        # a pixel below the clear threshold scores negative but still competes,
        # and qualityMosaic gives it a winner if nothing better exists. Gating
        # on the threshold here left most pixels unassigned while the reference
        # had full coverage. The declared 500 m erosion likewise belongs to the
        # daily clean-coverage statistic that selects candidate days, not to the
        # per-pixel mosaic. Both are retained as opt-in diagnostics.
        valid = np.isfinite(plane)
        if gate_threshold:
            valid &= plane >= float(reference["policy"]["cloud_score_clear_threshold"])
        if erosion_px > 0:
            valid = erode_valid(valid, radius_px=erosion_px)
        masked = np.where(valid, plane, np.nan)
        scalar = reference["day_scalars"].get(day, {})
        by_month[day[:7]].append(
            Candidate(
                day=day,
                score=masked,
                day_weight=float(scalar.get("weight", 0.0)),
                coverage_ratio=float(record.get("ratio", 0.0)),
            )
        )
    return dict(by_month)


def run(args: argparse.Namespace) -> dict[str, Any]:
    import numpy as np
    import rioxarray  # noqa: F401
    import xarray as xr

    reference = read_reference(Path(args.reference))
    teacher = np.load(Path(args.teacher), allow_pickle=False)
    transform = np.asarray(teacher["transform"], dtype=np.float64)
    height, width = teacher["surface"].shape[:2]
    x = transform[2] + transform[0] * (np.arange(width) + 0.5)
    y = transform[5] + transform[4] * (np.arange(height) + 0.5)
    template = xr.DataArray(
        np.zeros((height, width), dtype=np.float32), dims=("y", "x"), coords={"y": y, "x": x}
    ).rio.write_crs(str(teacher["crs"]))

    months = sorted({str(r["day"])[:7] for r in reference["images"]})
    if args.max_months > 0:
        months = months[: args.max_months]
    wanted = [r for r in reference["images"] if str(r["day"])[:7] in set(months)]
    if not args.skip_fetch:
        fetch_scores(
            tuple(args.bbox), [str(r["day"]) for r in wanted], Path(args.cache), args.edown
        )

    erosion_px = (
        int(round(float(reference["policy"].get("cloud_score_erosion_radius_m", 0.0)) / 20.0))
        if args.erode_mosaic
        else 0
    )
    subset = dict(reference)
    subset["images"] = wanted
    candidates = build_candidates(
        subset, Path(args.cache), template, erosion_px, args.gate_threshold, args.resampling
    )

    local_months, local_winners, _ = compose_monthly_winners(candidates)
    reference_index = {month: i for i, month in enumerate(reference["months"])}
    rows = []
    for position, month in enumerate(local_months):
        if month not in reference_index:
            continue
        stats = agreement(reference["winners"][reference_index[month]], local_winners[position])
        stats["month"] = month
        stats["candidates"] = len(candidates[month])
        rows.append(stats)
    return {
        "reference": str(args.reference),
        "erosion_px": erosion_px,
        "months_compared": len(rows),
        "per_month": rows,
        "mean_identical_fraction": float(np.mean([r["identical_fraction"] for r in rows]))
        if rows
        else float("nan"),
    }


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--reference", required=True)
    value.add_argument("--teacher", required=True)
    value.add_argument("--cache", required=True)
    value.add_argument("--bbox", nargs=4, type=float, required=True)
    value.add_argument("--edown", default="/home/users/marcyin/SIAC/tools/edown_runtime")
    value.add_argument("--max-months", type=int, default=2)
    value.add_argument("--skip-fetch", action="store_true")
    # Measured against the server-side index: nearest 0.965, average 0.950,
    # bilinear 0.891 mean identical fraction. Earth Engine samples the 10 m
    # cs_cdf at 20 m without an averaging reduction.
    value.add_argument(
        "--resampling",
        default="nearest",
        choices=("average", "nearest", "bilinear"),
        help="10 m cs_cdf -> 20 m grid reduction.",
    )
    value.add_argument(
        "--gate-threshold",
        action="store_true",
        help="Diagnostic: exclude pixels below the clear threshold.",
    )
    value.add_argument(
        "--erode-mosaic",
        action="store_true",
        help="Diagnostic: erode the mosaic input by the declared radius.",
    )
    value.add_argument("--output")
    return value


def main() -> None:
    args = parser().parse_args()
    payload = run(args)
    text = json.dumps(payload, indent=2)
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
