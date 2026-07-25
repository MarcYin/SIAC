"""Build matched-history controls directly in the current SIAC RT space.

This is an AERONET-free upper-bound control for the learned L2A harmonizer.
It uses the exact historical dates retained by the pair campaign, fetches the
matching L1C pixels, and applies their saved MAIAC-conditioned libRadtran
sidecars before the same daily and monthly medians are formed.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage
from tools.aeronet_validation.build_harmonized_l2a_histories import (
    DEFAULT_CASES,
    DEFAULT_OUTPUT,
    HALF_SIZE_DEGREES,
    MATCHUPS,
    _load_scenes,
    _nanmedian,
    _write_history,
)
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import (
    CANONICAL_BANDS,
    LAND_SCL,
    REFLECTANCE_SCALE,
)


def _fetch_target(
    ee: Any,
    grid: dict[str, Any],
    scene: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    import bestpixel as bp
    from bestpixel._gee import get_patch
    from bestpixel.atmosphere import AtmoSidecar
    from bestpixel.l1c import SIDECAR_TO_GEE

    sidecar = AtmoSidecar.load(str(scene["sidecar"]))
    scene_id = str(scene["scene_id"])
    if scene_id not in sidecar.scenes:
        raise KeyError(f"{scene_id} is absent from {scene['sidecar']}")
    band_ids = [SIDECAR_TO_GEE[name] for name in CANONICAL_BANDS]
    l1_raw = get_patch(ee, str(scene["l1c_asset"]), band_ids, grid=grid)
    cloud_score = get_patch(ee, str(scene["cloud_score_asset"]), ["cs"], grid=grid)[0]
    l2a_asset = f"COPERNICUS/S2_SR_HARMONIZED/{scene['l2a']['system_index']}"
    scl = get_patch(ee, l2a_asset, ["SCL"], grid=grid)[0]

    atmospheric = sidecar.scenes[scene_id]
    triples = [
        atmospheric.coeffs(sidecar.band_index(name), atmospheric.wvp) for name in CANONICAL_BANDS
    ]
    l1_toa = np.asarray(l1_raw, dtype=np.float32) * np.float32(REFLECTANCE_SCALE)
    target = np.asarray(
        bp.correct_toa(
            np.ascontiguousarray(l1_toa),
            [triple[0] for triple in triples],
            [triple[1] for triple in triples],
            [triple[2] for triple in triples],
        ),
        dtype=np.float32,
    )
    valid = (
        (np.asarray(cloud_score, dtype=np.float32) > 0.60)
        & np.isin(np.asarray(scl, dtype=np.int16), LAND_SCL)
        & np.all(np.isfinite(target), axis=0)
        & np.all((target > 0.001) & (target < 0.8), axis=0)
    )
    valid = ndimage.binary_erosion(valid, structure=np.ones((3, 3), dtype=bool), iterations=1)
    target[:, ~valid] = np.nan
    return target, valid


def _fetch_with_retries(
    ee: Any,
    grid: dict[str, Any],
    scene: dict[str, Any],
    *,
    attempts: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    for attempt in range(1, attempts + 1):
        try:
            return _fetch_target(ee, grid, scene)
        except Exception:  # noqa: BLE001 - bounded retry is recorded by caller
            if attempt == attempts:
                raise
            time.sleep(5.0 * attempt)
    raise RuntimeError("unreachable retry state")


def build_one(
    matchup_id: str,
    row: dict[str, str],
    *,
    pair_dir: Path,
    output_dir: Path,
    cutoff: str,
    force: bool,
) -> dict[str, Any]:
    from bestpixel._gee import init_ee, utm_epsg_from_bbox, utm_grid

    destination = output_dir / "paired_siac_daily" / f"{matchup_id}.npz"
    audit_path = output_dir / "paired_siac_audits" / f"{matchup_id}.json"
    if destination.exists() and audit_path.exists() and not force:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if audit.get("status") == "OK":
            return audit

    lon = float(row["longitude"])
    lat = float(row["latitude"])
    bbox = (
        lon - HALF_SIZE_DEGREES,
        lat - HALF_SIZE_DEGREES,
        lon + HALF_SIZE_DEGREES,
        lat + HALF_SIZE_DEGREES,
    )
    grid = utm_grid(bbox, utm_epsg_from_bbox(bbox), 60.0)
    scenes = _load_scenes(pair_dir / f"{matchup_id}.npz", cutoff)
    grouped: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for scene in scenes:
        grouped[str(scene["window"])][str(scene["day"])].append(scene)

    composites: list[np.ndarray] = []
    windows: list[str] = []
    used_scenes: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    ee = init_ee()
    for window in sorted(grouped):
        day_surfaces: list[np.ndarray] = []
        for day in sorted(grouped[window]):
            tile_surfaces: list[np.ndarray] = []
            for scene in grouped[window][day]:
                try:
                    target, valid = _fetch_with_retries(ee, grid, scene)
                except Exception as exc:  # noqa: BLE001 - audit every unavailable scene
                    errors.append(
                        {
                            "window": window,
                            "day": day,
                            "scene_id": str(scene.get("scene_id", "")),
                            "error_type": type(exc).__name__,
                            "reason": str(exc)[:500],
                        }
                    )
                    continue
                valid_count = int(np.count_nonzero(valid))
                if valid_count < 100:
                    errors.append(
                        {
                            "window": window,
                            "day": day,
                            "scene_id": str(scene.get("scene_id", "")),
                            "error_type": "SparseScene",
                            "reason": f"only {valid_count} clear-land pixels",
                        }
                    )
                    continue
                tile_surfaces.append(target)
                used_scenes.append(
                    {
                        "window": window,
                        "day": day,
                        "scene_id": scene.get("scene_id"),
                        "maiac_aot": scene.get("maiac_aot"),
                        "valid_pixels": valid_count,
                    }
                )
            if tile_surfaces:
                day_surfaces.append(_nanmedian(tile_surfaces))
        if day_surfaces:
            windows.append(window)
            composites.append(_nanmedian(day_surfaces))

    status = "OK" if len(windows) >= 2 else "FAILED"
    if status == "OK":
        sensing = row["sensing_time_utc"]
        _write_history(
            destination,
            composites,
            windows,
            grid=grid,
            scene_year=int(sensing[:4]),
            scene_month=int(sensing[5:7]),
            provenance={
                "uses_aeronet": False,
                "candidate": "paired_siac_daily",
                "history_cutoff": cutoff,
                "surface_source": "same-day Sentinel-2 L1C corrected by saved MAIAC/current-RT sidecar",
                "application": "per acquisition before tile mosaic and temporal median",
                "role": "mapping upper-bound control",
            },
        )
    audit = {
        "status": status,
        "matchup_id": matchup_id,
        "uses_aeronet": False,
        "window_count": len(windows),
        "windows": windows,
        "scene_count": len(used_scenes),
        "scenes": used_scenes,
        "errors": errors,
        "output": str(destination),
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matchup_id", nargs="+")
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--history-cutoff", default="2023-12-31")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    requested = set(args.matchup_id)
    case_ids = {row["matchup_id"] for row in csv.DictReader(args.cases.open())}
    if not requested <= case_ids:
        raise SystemExit("requested matchup is outside the locked development case list")
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    exit_code = 0
    for matchup_id in args.matchup_id:
        audit = build_one(
            matchup_id,
            rows[matchup_id],
            pair_dir=args.pairs,
            output_dir=args.output,
            cutoff=args.history_cutoff,
            force=args.force,
        )
        print("PAIRED_HISTORY_DONE " + json.dumps(audit), flush=True)
        if audit["status"] != "OK":
            exit_code = 1
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
