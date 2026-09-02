#!/usr/bin/env python3
"""Build a Cloud Score+ winner index with a locally composed quality mosaic.

Earth Engine is used only as a data source: ``edown`` fetches per-image
``cs`` rasters and everything else -- clean coverage, candidate selection,
daily weighting and the per-pixel argmax -- runs here. That removes the L2A
dependency of the SCL proxy index and makes the weighting an offline choice
instead of one baked into a server-side ``qualityMosaic``.

Output matches ``siac_l1c_cloudscore_winner_index_v1`` so existing consumers
read it unchanged, except that ``winner_source`` truthfully records the local
composition.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
from collections import defaultdict
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
from tools.aeronet_validation.cloudscore_index_policy import (
    CLEAR_THRESHOLD,
    CLOUD_SCORE_BAND,
    LIBRARY_YEARS,
    clean_coverage,
    daily_weights,
    index_policy,
    seasonal_windows,
    select_candidate_days,
)
from tools.aeronet_validation.cloudscore_local_mosaic import Candidate, compose_monthly_winners
from tools.aeronet_validation.scout_scl import archive_grid_spec

COLLECTION = "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED"
BAND = CLOUD_SCORE_BAND
FILE_PREFIX = "GOOGLE_CLOUD_SCORE_PLUS_V1_S2_HARMONIZED"
SCHEMA = "siac_l1c_cloudscore_winner_index_v1"

#: Chunk budget for a Cloud Score+ fetch, below edown's 48 MiB default.
#: edown sizes chunks against this limit but the server measured one at
#: 50,466,816 bytes against its own 50,331,648 ceiling -- a 0.3% overshoot, so
#: the estimate omits per-request overhead. Rejection is total for the image, so
#: the margin is worth more than the extra chunks it costs: without it, roughly
#: 1.1% of AOIs failed outright at the fetch.
REQUEST_BYTE_LIMIT = 32 * 1024 * 1024

#: ``day -> AOD`` lookup. Injectable so the policy can be validated against a
#: reference index's own AOD values before the MAIAC path is trusted.
AODLookup = Callable[[Sequence[str]], "dict[str, float | None]"]


def grid_template(
    crs: str, transform_values: Sequence[float], height: int, width: int
) -> tuple[Any, tuple[float, float, float, float]]:
    """Build a resampling template and its WGS84 bbox from a grid spec."""
    import rioxarray  # noqa: F401
    import xarray as xr
    from pyproj import Transformer

    transform = np.asarray(transform_values, dtype=np.float64)
    x = transform[2] + transform[0] * (np.arange(width) + 0.5)
    y = transform[5] + transform[4] * (np.arange(height) + 0.5)
    template = xr.DataArray(
        np.zeros((height, width), dtype=np.float32), dims=("y", "x"), coords={"y": y, "x": x}
    ).rio.write_crs(crs)
    xs = (transform[2], transform[2] + transform[0] * width)
    ys = (transform[5], transform[5] + transform[4] * height)
    lon, lat = Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform(
        [min(xs), max(xs)], [min(ys), max(ys)]
    )
    return template, (min(lon), min(lat), max(lon), max(lat))


def scene_grid(teacher_path: Path) -> tuple[Any, tuple[float, float, float, float]]:
    """Build the 20 m template and its WGS84 bbox from a teacher archive."""

    # Shared with the scout: campaign archives disagree on both the raster key
    # and how they name the CRS, so the spec is resolved in one place.
    return grid_template(*archive_grid_spec(teacher_path))


def fetch_windows(
    bbox: tuple[float, float, float, float],
    windows: Sequence[tuple[str, str]],
    cache: Path,
    edown: str,
    *,
    timeout_s: int = 1800,
) -> None:
    """Fetch the clear-score band for each seasonal window. edown resumes, so this is idempotent."""
    cache.mkdir(parents=True, exist_ok=True)
    for start, end in windows:
        result = subprocess.run(
            [
                edown,
                "download",
                "--collection-id",
                COLLECTION,
                "--band",
                BAND,
                "--bbox",
                *[str(c) for c in bbox],
                "--start-date",
                start,
                "--end-date",
                end,
                "--output-root",
                str(cache),
            ],
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
        if result.returncode != 0:
            raise RuntimeError(f"edown failed for {start}..{end}: {result.stderr[-400:]}")


def fetch_selected(
    bbox: tuple[float, float, float, float],
    image_ids: Sequence[str],
    cache: Path,
    edown: str,
    *,
    timeout_s: int = 3600,
) -> None:
    """Fetch exactly the shortlisted acquisitions, in one call.

    The seasonal windows only ever existed to bound the search; once the scout
    has named the acquisitions, the id list is the selection and a single span
    from the first to the last shortlisted day is enough. edown filters
    server-side on ``system:index``, so the span being wide costs nothing --
    the collection is already reduced to the shortlist before any count or
    metadata transfer. Requires edown >= 0.2.2.
    """

    if not image_ids:
        raise ValueError("no image ids selected; refusing to fetch a whole window by accident")
    days = sorted({image_day(image) for image in image_ids})
    cache.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [
            edown,
            "download",
            "--collection-id",
            COLLECTION,
            "--band",
            BAND,
            "--bbox",
            *[str(coordinate) for coordinate in bbox],
            "--start-date",
            days[0],
            "--end-date",
            days[-1],
            "--image-id",
            ",".join(sorted(image_ids)),
            "--request-byte-limit",
            str(REQUEST_BYTE_LIMIT),
            "--output-root",
            str(cache),
        ],
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"edown failed for {len(image_ids)} selected images: {result.stderr[-400:]}"
        )
    assert_direct_pixel_reads(cache)


def assert_direct_pixel_reads(cache: Path) -> dict[str, int]:
    """Fail if any image was served by ``computePixels`` rather than ``getPixels``.

    A plain band read is a direct asset read; only a scale map or transform
    plugin makes Earth Engine evaluate an expression graph, which is markedly
    slower. Nothing in this pipeline needs one, so a ``computePixels`` result
    means an option leaked in -- a change that costs a lot of wall clock while
    producing identical pixels, and so would otherwise go unnoticed.

    Requires edown >= 0.2.2, which records the API per image. Older manifests
    do not carry the field and are reported as unknown rather than failing.
    """

    counts: dict[str, int] = defaultdict(int)
    for manifest in sorted(cache.glob("manifests/*.json")):
        document = json.loads(manifest.read_text(encoding="utf-8"))
        for record in (document.get("download") or {}).get("results", []) or []:
            if record.get("status") != "downloaded":
                continue
            counts[str(record.get("pixel_api") or "unknown")] += 1
    if counts.get("computePixels"):
        raise RuntimeError(
            f"{counts['computePixels']} images were fetched with computePixels; "
            "a plain band read must use getPixels. Check for a scale map or "
            "transform plugin on the download config."
        )
    return dict(counts)


def read_planes(cache: Path, template: Any) -> dict[str, np.ndarray]:
    """Read every cached clear-score raster onto the scene grid, keyed by image id.

    Nearest resampling is deliberate: it reproduced the server-side winners
    best (0.965 versus 0.950 for average and 0.891 for bilinear), because Earth
    Engine samples the 10 m band at 20 m without an averaging reduction.
    """
    import rioxarray as rxr
    from rasterio.enums import Resampling

    planes: dict[str, np.ndarray] = {}
    for path in sorted(cache.glob(f"images/*/{FILE_PREFIX}_*.tif")):
        image_id = path.stem[len(FILE_PREFIX) + 1 :]
        raster = rxr.open_rasterio(path).squeeze("band", drop=True)
        planes[image_id] = raster.rio.reproject_match(
            template, resampling=Resampling.nearest
        ).values.astype(np.float64)
    return planes


def image_day(image_id: str) -> str:
    """Acquisition day from an Earth Engine S2 image id (``YYYYMMDDTHHMMSS_...``)."""
    stamp = image_id.split("_")[0]
    return dt.datetime.strptime(stamp[:8], "%Y%m%d").date().isoformat()


def build_index(
    planes: dict[str, np.ndarray],
    aod_lookup: AODLookup,
    *,
    clear_threshold: float = CLEAR_THRESHOLD,
    library_years: Sequence[int] = LIBRARY_YEARS,
) -> dict[str, Any]:
    """Compose the index payload from per-image clear-score planes."""
    if not planes:
        raise ValueError("no clear-score planes available")
    coverage = {
        image: clean_coverage(p, clear_threshold=clear_threshold) for image, p in planes.items()
    }
    ratio = {image: float(np.mean(np.isfinite(p))) for image, p in planes.items()}

    by_day: dict[str, list[str]] = defaultdict(list)
    for image in planes:
        by_day[image_day(image)].append(image)
    # A day's clean coverage is its best contributing image: a partial tile must
    # not drag down a day that is otherwise clear over the AOI.
    day_coverage = {day: max(coverage[i] for i in images) for day, images in by_day.items()}

    candidate_days = select_candidate_days(day_coverage)
    if not candidate_days:
        raise ValueError("candidate rule selected no days")
    weights = daily_weights(
        aod_lookup(candidate_days), {d: day_coverage[d] for d in candidate_days}
    )

    by_month: dict[str, list[Candidate]] = defaultdict(list)
    image_rows: list[dict[str, Any]] = []
    for day in candidate_days:
        for image in sorted(by_day[day]):
            by_month[day[:7]].append(
                Candidate(
                    day=day,
                    score=planes[image],
                    day_weight=float(weights[day]["weight"]),
                    coverage_ratio=ratio[image],
                )
            )
            image_rows.append({"day": day, "idx": image, "ratio": ratio[image]})

    months, winners, ordering = compose_monthly_winners(by_month)
    return {
        "months": months,
        "winners": winners,
        "image_table": image_rows,
        "day_scalars": [weights[day] for day in candidate_days],
        "ordering": ordering,
        "policy": index_policy(library_years=library_years),
        "counts": {
            "images_available": len(planes),
            "days_available": len(day_coverage),
            "candidate_days": len(candidate_days),
            "months": len(months),
        },
    }


def write_index(payload: dict[str, Any], output: Path, matchup_id: str) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        months=np.asarray(payload["months"], dtype="<U7"),
        winners=payload["winners"].astype(np.int16),
        image_table=np.asarray([json.dumps(r, sort_keys=True) for r in payload["image_table"]]),
        day_scalars=np.asarray([json.dumps(r, sort_keys=True) for r in payload["day_scalars"]]),
        index_policy_json=np.asarray(json.dumps(payload["policy"], sort_keys=True)),
    )
    receipt = {
        "schema_version": SCHEMA,
        "matchup_id": matchup_id,
        "status": "ok",
        **payload["counts"],
        "index_policy": payload["policy"],
    }
    output.with_suffix(".json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")


def maiac_lookup(bbox: tuple[float, float, float, float], crs: str = "EPSG:4326") -> AODLookup:
    """Per-day MAIAC AOD over the AOI; days without a retrieval return None."""

    def lookup(days: Sequence[str]) -> dict[str, float | None]:
        from siac.adapters.atmo.maiac_day_aod import MAIACDayAODProvider

        periods = {(int(d[:4]), int(d[5:7])) for d in days}
        resolved = MAIACDayAODProvider().day_aod_map(bbox, crs, sorted(periods))
        return {day: resolved.get(day) for day in days}

    return lookup


def run(args: argparse.Namespace) -> dict[str, Any]:
    template, bbox = scene_grid(Path(args.teacher))
    scene_day = (
        dt.datetime.strptime(args.matchup_id.split("_")[-1][:8], "%Y%m%d").date().isoformat()
    )
    windows = seasonal_windows(scene_day, library_years=tuple(args.library_years))
    cache = Path(args.cache)
    if not args.skip_fetch:
        fetch_windows(bbox, windows, cache, args.edown)
    planes = read_planes(cache, template)

    if args.reference_aod:
        reference = np.load(args.reference_aod, allow_pickle=False)
        stored = {
            json.loads(str(r))["day"]: json.loads(str(r))["aod"] for r in reference["day_scalars"]
        }
        lookup: AODLookup = lambda days: {d: stored.get(d) for d in days}  # noqa: E731
    else:
        lookup = maiac_lookup(bbox)

    payload = build_index(planes, lookup, library_years=tuple(args.library_years))
    if args.output:
        write_index(payload, Path(args.output), args.matchup_id)
    return {"matchup_id": args.matchup_id, "windows": len(windows), **payload["counts"]}


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--matchup-id", required=True)
    value.add_argument("--teacher", required=True, help="Teacher npz supplying the 20 m grid.")
    value.add_argument("--cache", required=True)
    value.add_argument("--output")
    value.add_argument("--edown", default="/home/users/marcyin/SIAC/tools/edown_runtime")
    value.add_argument("--library-years", type=int, nargs="+", default=list(LIBRARY_YEARS))
    value.add_argument("--skip-fetch", action="store_true")
    value.add_argument(
        "--reference-aod",
        help="Reference index whose day AOD values replace the MAIAC lookup, "
        "isolating the mosaic and candidate rule from MAIAC availability.",
    )
    return value


def main() -> None:
    args = parser().parse_args()
    print(json.dumps(run(args), indent=2))


if __name__ == "__main__":
    main()
