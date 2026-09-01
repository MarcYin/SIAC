#!/usr/bin/env python3
"""Does scouting first change the winner index it was meant to make cheaper?

The prefilter's whole claim is that Cloud Score+ only needs to score a
shortlist, because the acquisitions it drops could not have won a pixel anyway.
That is a testable claim, and it is testable without spending any quota: a
matchup whose full window is already cached can be composed twice from the same
planes -- once from all of them, once from the scouted shortlist -- and the two
winner fields compared.

Winners are compared by *acquisition identity*, not by index. Removing
candidates renumbers a month's winner values, so equal indices would mean
nothing; each index is mapped back through the month's ordering first.

Both arms are also compared against the committed Earth Engine index, which
separates two different questions: whether the local composer reproduces the
server (already measured) and whether shortlisting perturbs it (measured here).
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from tools.aeronet_validation.build_cloudscore_index_local import (
    build_index,
    image_day,
    read_planes,
    scene_grid,
)
from tools.aeronet_validation.build_l2a_scl_index import SAS_URL, _imports, _session
from tools.aeronet_validation.cloudscore_index_policy import seasonal_windows
from tools.aeronet_validation.scout_prefilter import shortlist
from tools.aeronet_validation.scout_scl import grid_from_archive, query_items, scout

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence


def winner_days(payload: dict[str, Any]) -> dict[str, np.ndarray]:
    """Per-month plane of winning *day*, resolving indices through the ordering."""

    resolved: dict[str, np.ndarray] = {}
    for position, month in enumerate(payload["months"]):
        order = list(payload["ordering"][month])
        labels = np.asarray(order, dtype=object)
        index = payload["winners"][position]
        flat = np.clip(index.astype(np.int64), 0, max(0, len(order) - 1))
        plane = labels[flat] if order else np.full(index.shape, "", dtype=object)
        resolved[month] = np.where(index >= 0, plane, "")
    return resolved


def compare(left: dict[str, np.ndarray], right: dict[str, np.ndarray]) -> dict[str, Any]:
    """Pixel agreement on the winning day, month by month."""

    shared = sorted(set(left) & set(right))
    rows = []
    for month in shared:
        same = np.asarray(left[month] == right[month])
        rows.append({"month": month, "identical_fraction": float(np.mean(same))})
    return {
        "months_compared": len(shared),
        "months_only_left": sorted(set(left) - set(right)),
        "months_only_right": sorted(set(right) - set(left)),
        "per_month": rows,
        "mean_identical_fraction": (
            float(np.mean([r["identical_fraction"] for r in rows])) if rows else 0.0
        ),
        "months_identical": sum(1 for r in rows if r["identical_fraction"] == 1.0),
    }


def winner_quality(
    payload: dict[str, Any], planes: Mapping[str, np.ndarray]
) -> dict[str, dict[str, np.ndarray]]:
    """Per-pixel clearness and day AOD of whichever acquisition won.

    "How many pixels changed hands" cannot say whether a shortlist is worse; a
    different winner is only a loss if it is a *dimmer* one. This recovers what
    each arm actually achieved at every pixel so the two can be compared on
    quality rather than on identity.
    """

    weights = {str(row["day"]): row for row in payload["day_scalars"]}
    by_month: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in payload["image_table"]:
        by_month[str(row["day"])[:7]].append(row)
    achieved: dict[str, dict[str, np.ndarray]] = {}
    for position, month in enumerate(payload["months"]):
        rows = by_month.get(month, [])
        if not rows:
            continue
        scores = np.stack([np.asarray(planes[str(row["idx"])], dtype=np.float64) for row in rows])
        aod = np.asarray(
            [float(weights[str(row["day"])].get("aod") or np.nan) for row in rows],
            dtype=np.float64,
        )
        index = np.clip(payload["winners"][position].astype(np.int64), 0, len(rows) - 1)
        picked = np.take_along_axis(scores, index[None, ...], axis=0)[0]
        achieved[month] = {"score": picked, "aod": aod[index]}
    return achieved


def quality_delta(
    left: Mapping[str, dict[str, np.ndarray]],
    right: Mapping[str, dict[str, np.ndarray]],
    left_days: Mapping[str, np.ndarray],
    right_days: Mapping[str, np.ndarray],
) -> dict[str, Any]:
    """Compare achieved quality, overall and where the winner changed."""

    months = sorted(set(left) & set(right) & set(left_days) & set(right_days))
    changed_left: list[np.ndarray] = []
    changed_right: list[np.ndarray] = []
    aod_left: list[np.ndarray] = []
    aod_right: list[np.ndarray] = []
    all_left: list[np.ndarray] = []
    all_right: list[np.ndarray] = []
    for month in months:
        differs = np.asarray(left_days[month] != right_days[month])
        all_left.append(left[month]["score"].ravel())
        all_right.append(right[month]["score"].ravel())
        if differs.any():
            changed_left.append(left[month]["score"][differs])
            changed_right.append(right[month]["score"][differs])
            aod_left.append(left[month]["aod"][differs])
            aod_right.append(right[month]["aod"][differs])

    def summarize(values: list[np.ndarray]) -> float:
        if not values:
            return float("nan")
        stacked = np.concatenate(values)
        finite = stacked[np.isfinite(stacked)]
        return float(np.mean(finite)) if finite.size else float("nan")

    changed_l, changed_r = summarize(changed_left), summarize(changed_right)
    swapped = (
        float(np.mean(np.concatenate(changed_right) >= np.concatenate(changed_left)))
        if changed_left
        else float("nan")
    )
    return {
        "mean_winner_score_left": summarize(all_left),
        "mean_winner_score_right": summarize(all_right),
        "changed_pixels": int(sum(v.size for v in changed_left)),
        "changed_mean_winner_score_left": changed_l,
        "changed_mean_winner_score_right": changed_r,
        "changed_score_delta": changed_r - changed_l,
        "changed_fraction_right_at_least_left": swapped,
        "changed_mean_winner_aod_left": summarize(aod_left),
        "changed_mean_winner_aod_right": summarize(aod_right),
    }


def reference_winner_days(path: Path) -> dict[str, np.ndarray]:
    """Winning day per month from a committed index archive."""

    archive = np.load(path, allow_pickle=False)
    months = [str(m) for m in archive["months"]]
    winners = np.asarray(archive["winners"])
    table = [json.loads(str(row)) for row in archive["image_table"]]
    # A month's winner index addresses its candidate list, which is one entry
    # per *image* in image_table order -- not one per distinct day. Collapsing
    # to unique days silently renumbers every month that has two acquisitions on
    # one day, which is most of them.
    by_month: dict[str, list[str]] = defaultdict(list)
    for row in table:
        by_month[str(row["day"])[:7]].append(str(row["day"]))
    resolved: dict[str, np.ndarray] = {}
    for position, month in enumerate(months):
        order = by_month.get(month, [])
        if not order:
            continue
        labels = np.asarray(order, dtype=object)
        index = winners[position]
        flat = np.clip(index.astype(np.int64), 0, len(order) - 1)
        resolved[month] = np.where(index >= 0, labels[flat], "")
    return resolved


def scouted_image_ids(
    grid_archive: Path,
    scene_day: str,
    *,
    library_years: Sequence[int],
    top_k: int,
    workers: int,
) -> tuple[set[str], dict[str, Any]]:
    """Shortlist the window from free signals, returning the surviving ids."""

    imports = _imports()
    grid = grid_from_archive(grid_archive, imports=imports)
    session = _session(imports["requests"])
    sas = str(session.get(SAS_URL, timeout=60).json()["token"])
    items = query_items(session, grid, seasonal_windows(scene_day, library_years=library_years))
    records, diagnostics = scout(items, sas=sas, grid=grid, imports=imports, workers=workers)
    chosen = shortlist(records, top_k=top_k)
    ids = {record.image_id for entries in chosen.values() for record in entries}
    return ids, {**diagnostics, "shortlisted": len(ids), "scouted": len(records)}


def restrict(planes: Mapping[str, np.ndarray], keep: set[str]) -> dict[str, np.ndarray]:
    """Keep only shortlisted acquisitions, matching on the sensing/datastrip id."""

    return {image: plane for image, plane in planes.items() if image in keep}


def run(args: argparse.Namespace) -> dict[str, Any]:
    template, bbox = scene_grid(Path(args.teacher))
    planes = read_planes(Path(args.cache), template)
    if not planes:
        raise SystemExit(f"no cached Cloud Score+ planes under {args.cache}")
    years = tuple(args.library_years)
    days = sorted({image_day(image) for image in planes})
    scene_day = args.scene_day or (
        f"{args.matchup_id.split('_')[-1][:4]}-"
        f"{args.matchup_id.split('_')[-1][4:6]}-{args.matchup_id.split('_')[-1][6:8]}"
    )

    keep, diagnostics = scouted_image_ids(
        Path(args.grid_archive or args.teacher),
        scene_day,
        library_years=years,
        top_k=args.top_k,
        workers=args.workers,
    )
    matched = set(planes) & keep
    shortlisted = restrict(planes, keep)

    # The per-day AOD weight is an additive term in the quality ordering, so
    # feeding a constant where the committed index used measured MAIAC reorders
    # candidates across days and would be scored as a mosaic disagreement. Both
    # arms therefore reuse the reference's own day AOD when one is supplied,
    # which is what isolates the shortlist from MAIAC availability.
    stored: dict[str, float | None] = {}
    if args.reference:
        reference_archive = np.load(Path(args.reference), allow_pickle=False)
        for row in reference_archive["day_scalars"]:
            record = json.loads(str(row))
            stored[str(record["day"])] = record.get("aod")

    def lookup(wanted: Sequence[str]) -> dict[str, float | None]:
        return {day: stored.get(day) for day in wanted}

    full_index = build_index(dict(planes), lookup, library_years=years)
    short_index = build_index(shortlisted, lookup, library_years=years)

    result: dict[str, Any] = {
        "matchup_id": args.matchup_id,
        "scene_day": scene_day,
        "library_years": list(years),
        "cached_planes": len(planes),
        "cached_days": len(days),
        "scout": diagnostics,
        "shortlist_matched_cache": len(matched),
        "planes_dropped": len(planes) - len(shortlisted),
        "request_saving": 1.0 - len(shortlisted) / len(planes),
        "full_counts": full_index["counts"],
        "shortlist_counts": short_index["counts"],
        "shortlist_vs_full": compare(winner_days(full_index), winner_days(short_index)),
        "shortlist_quality_vs_full": quality_delta(
            winner_quality(full_index, planes),
            winner_quality(short_index, shortlisted),
            winner_days(full_index),
            winner_days(short_index),
        ),
        "day_aod_source": "reference_day_scalars" if args.reference else "locked_constant",
        "days_with_reference_aod": sum(1 for value in stored.values() if value is not None),
    }
    if args.reference:
        reference = reference_winner_days(Path(args.reference))
        result["full_vs_reference"] = compare(reference, winner_days(full_index))
        result["shortlist_vs_reference"] = compare(reference, winner_days(short_index))
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--matchup-id", required=True)
    value.add_argument("--teacher", required=True, help="npz supplying the 20 m grid.")
    value.add_argument("--grid-archive", help="Grid for scouting; defaults to --teacher.")
    value.add_argument("--cache", required=True, help="edown cache of Cloud Score+ rasters.")
    value.add_argument("--reference", help="Committed index archive to compare both arms to.")
    value.add_argument("--scene-day", help="Override the day parsed from the matchup id.")
    value.add_argument("--library-years", type=int, nargs="+", default=list(range(2018, 2024)))
    value.add_argument("--top-k", type=int, default=5)
    value.add_argument("--workers", type=int, default=8)
    value.add_argument("--output")
    return value


def main() -> None:
    result = run(parser().parse_args())
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
