"""Build one-L2A-source histories with same-acquisition physical harmonization.

Each retained historical L2A acquisition is mapped to the current L1C/RT BOA
frame using only clear pixels from that same historical acquisition.  The map
is applied before tile mosaicking and temporal compositing.  AERONET is never
read; the output remains an operational L2A surface prior.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from tools.aeronet_validation.build_harmonized_l2a_histories import (
    DEFAULT_CASES,
    HALF_SIZE_DEGREES,
    MATCHUPS,
    _fetch_scene,
    _load_scenes,
    _nanmedian,
    _write_history,
)
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import TARGET_RT
from tools.aeronet_validation.per_acquisition_l2a_l1c_mapping import (
    apply_band_map,
    build_scene_maps,
    scene_key,
)
from tools.aeronet_validation.train_l2a_l1c_harmonizer import BAND_NAMES, VISIBLE_INDICES

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_PAIRS = ROOT / "analysis/l2a_l1c_exact_pairs_mediumdev_20260713"
DEFAULT_OUTPUT = ROOT / "analysis/l2a_l1c_per_acquisition_affine_20260715/daily_histories"


def _candidate(value: str) -> tuple[str, float]:
    parts = value.split(":")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("candidate must be MODE:CAP")
    mode, cap_text = parts
    if mode not in {"all", "visible", "solver", "blue"}:
        raise argparse.ArgumentTypeError("MODE must be all, visible, solver, or blue")
    try:
        cap = float(cap_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("CAP must be numeric") from exc
    if not 0.0 < cap <= 0.2:
        raise argparse.ArgumentTypeError("CAP must be in (0, 0.2]")
    return mode, cap


def _tag(candidate: tuple[str, float]) -> str:
    mode, cap = candidate
    return f"pair_affine_{mode}_cap{cap:.3f}".replace(".", "p")


def _band_indices(mode: str) -> tuple[int, ...]:
    if mode == "all":
        return tuple(range(len(BAND_NAMES)))
    if mode == "visible":
        return VISIBLE_INDICES
    if mode == "solver":
        return (1, 2, 3)
    if mode == "blue":
        return (1,)
    raise ValueError(f"unsupported mapping mode {mode!r}")


def _correct_scene(
    fetched: dict[str, np.ndarray],
    scene_map: dict[str, Any] | None,
    *,
    candidate: tuple[str, float],
) -> tuple[np.ndarray, dict[str, Any]]:
    mode, cap = candidate
    surface = np.asarray(fetched["surface"], dtype=np.float64).copy()
    valid = np.asarray(fetched["valid"], dtype=bool)
    selected = np.flatnonzero(valid.ravel())
    if selected.size == 0:
        return surface.astype(np.float32), {"valid_pixels": 0, "bands": {}}
    if scene_map is None:
        return surface.astype(np.float32), {
            "valid_pixels": int(selected.size),
            "mapping_applied": False,
            "skip_reason": "no_same_acquisition_pair",
            "bands": {},
        }

    flat = surface.reshape(len(BAND_NAMES), -1).T
    band_stats: dict[str, Any] = {}
    mapped_bands = 0
    for band_index in _band_indices(mode):
        band_name = BAND_NAMES[band_index]
        band_map = (scene_map.get("bands") or {}).get(band_name)
        if not isinstance(band_map, dict) or not bool(band_map.get("valid")):
            band_stats[band_name] = {
                "mapping_applied": False,
                "skip_reason": None if not isinstance(band_map, dict) else band_map.get("reason"),
            }
            continue
        corrected, correction = apply_band_map(flat[selected, band_index], band_map, cap=cap)
        flat[selected, band_index] = corrected
        mapped_bands += 1
        band_stats[band_name] = {
            "mapping_applied": True,
            "method": band_map["method"],
            "intercept": band_map["intercept"],
            "slope": band_map["slope"],
            "paired_sample_count": band_map["sample_count"],
            "inlier_count": band_map["inlier_count"],
            "source_iqr": band_map["source_iqr"],
            "correlation": band_map["correlation"],
            "pair_residual_rmse": band_map["residual_rmse"],
            "median_correction": float(np.nanmedian(correction)),
            "p95_abs_correction": float(np.nanpercentile(np.abs(correction), 95.0)),
            "cap_fraction": float(np.mean(np.abs(correction) >= cap - 1.0e-12)),
        }
    return surface.astype(np.float32), {
        "valid_pixels": int(selected.size),
        "mapping_applied": bool(mapped_bands),
        "mode": mode,
        "cap": cap,
        "pair_scene": {
            "scene_id": scene_map.get("scene_id"),
            "day": scene_map.get("day"),
            "paired_sample_count": scene_map.get("paired_sample_count"),
        },
        "bands": band_stats,
    }


def build_one(
    matchup_id: str,
    row: dict[str, str],
    *,
    pair_dir: Path,
    output_dir: Path,
    candidates: list[tuple[str, float]],
    cutoff: str,
    min_samples: int,
    min_source_iqr: float,
    pair_maiac_aot_max: float | None,
    force: bool,
) -> dict[str, Any]:
    from bestpixel._gee import init_ee, utm_epsg_from_bbox, utm_grid

    outputs = {
        "identity": output_dir / "identity_daily" / f"{matchup_id}.npz",
        **{
            _tag(candidate): output_dir / _tag(candidate) / f"{matchup_id}.npz"
            for candidate in candidates
        },
    }
    audit_path = output_dir / "audits" / f"{matchup_id}.json"
    if not force and audit_path.exists() and all(path.exists() for path in outputs.values()):
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if audit.get("status") == "OK":
            return audit

    pair_path = pair_dir / f"{matchup_id}.npz"
    pair_audit = json.loads(pair_path.with_suffix(".json").read_text(encoding="utf-8"))
    scene_maps = build_scene_maps(
        pair_path,
        cutoff=cutoff,
        min_samples=min_samples,
        min_source_iqr=min_source_iqr,
        maiac_aot_max=pair_maiac_aot_max,
    )
    lon = float(row["longitude"])
    lat = float(row["latitude"])
    bbox = (
        lon - HALF_SIZE_DEGREES,
        lat - HALF_SIZE_DEGREES,
        lon + HALF_SIZE_DEGREES,
        lat + HALF_SIZE_DEGREES,
    )
    grid = utm_grid(bbox, utm_epsg_from_bbox(bbox), 60.0)
    scenes = _load_scenes(pair_path, cutoff)
    by_window: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for scene in scenes:
        by_window[str(scene["window"])][str(scene["day"])].append(scene)

    products: dict[str, list[np.ndarray]] = {name: [] for name in outputs}
    windows: list[str] = []
    scene_audits: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    ee = init_ee()
    for window in sorted(by_window):
        window_days: dict[str, list[np.ndarray]] = {name: [] for name in outputs}
        for day in sorted(by_window[window]):
            day_tiles: dict[str, list[np.ndarray]] = {name: [] for name in outputs}
            for scene in by_window[window][day]:
                try:
                    fetched = _fetch_scene(ee, grid, scene)
                except Exception as exc:  # noqa: BLE001 - retain usable acquisitions
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
                if int(np.count_nonzero(fetched["valid"])) < 100:
                    errors.append(
                        {
                            "window": window,
                            "day": day,
                            "scene_id": str(scene.get("scene_id", "")),
                            "error_type": "SparseScene",
                            "reason": "fewer than 100 clear-land pixels",
                        }
                    )
                    continue
                day_tiles["identity"].append(fetched["surface"])
                audit: dict[str, Any] = {
                    "window": window,
                    "day": day,
                    "scene_id": scene.get("scene_id"),
                    "pair_scene_key": scene_key(scene),
                    "candidates": {},
                }
                for candidate in candidates:
                    tag = _tag(candidate)
                    corrected, stats = _correct_scene(
                        fetched,
                        scene_maps.get(scene_key(scene)),
                        candidate=candidate,
                    )
                    day_tiles[tag].append(corrected)
                    audit["candidates"][tag] = stats
                scene_audits.append(audit)
            if not day_tiles["identity"]:
                continue
            for name, tiles in day_tiles.items():
                window_days[name].append(_nanmedian(tiles))
        if not window_days["identity"]:
            continue
        windows.append(window)
        for name, days in window_days.items():
            products[name].append(_nanmedian(days))

    common_provenance = {
        "uses_aeronet": False,
        "application": "same acquisition map before tile mosaic and temporal median",
        "history_cutoff": cutoff,
        "pair_archive": str(pair_path),
        "pair_target": pair_audit.get("target"),
        # Older pair audits predate explicit target_rt provenance.  Their
        # archive is still tied to the builder's immutable target settings.
        "pair_target_rt": pair_audit.get("target_rt") or TARGET_RT,
        "mapping": "robust same-acquisition L2A to current-RT affine frame conversion",
        "minimum_paired_samples": int(min_samples),
        "minimum_source_iqr": float(min_source_iqr),
        "pair_maiac_aot_max": pair_maiac_aot_max,
        "available_pair_scene_count": len(scene_maps),
    }
    if len(windows) < 2:
        status = {
            "status": "FAILED",
            "matchup_id": matchup_id,
            "windows": windows,
            "errors": errors,
            **common_provenance,
        }
    else:
        sensing = row["sensing_time_utc"]
        for name, destination in outputs.items():
            _write_history(
                destination,
                products[name],
                windows,
                grid=grid,
                scene_year=int(sensing[:4]),
                scene_month=int(sensing[5:7]),
                provenance={
                    **common_provenance,
                    "candidate": None if name == "identity" else name,
                },
            )
        status = {
            "status": "OK",
            "matchup_id": matchup_id,
            "uses_aeronet": False,
            "windows": windows,
            "window_count": len(windows),
            "scene_count": len(scene_audits),
            "errors": errors,
            "outputs": {name: str(path) for name, path in outputs.items()},
            "scenes": scene_audits,
            **common_provenance,
        }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return status


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matchup_id", nargs="+")
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--history-cutoff", default="2023-12-31")
    parser.add_argument("--candidate", action="append", type=_candidate, default=[])
    parser.add_argument("--min-samples", type=int, default=200)
    parser.add_argument("--min-source-iqr", type=float, default=0.01)
    parser.add_argument(
        "--pair-maiac-aot-max",
        type=float,
        help="use maps only from historical MAIAC AOD values at or below this clean-scene limit",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.min_samples < 100:
        raise SystemExit("--min-samples must be at least 100")
    if not 0.0 < args.min_source_iqr < 0.2:
        raise SystemExit("--min-source-iqr must be in (0, 0.2)")
    if args.pair_maiac_aot_max is not None and not 0.0 < args.pair_maiac_aot_max <= 2.0:
        raise SystemExit("--pair-maiac-aot-max must be in (0, 2]")
    candidates = args.candidate or [("visible", 0.015), ("visible", 0.03), ("all", 0.03)]
    requested = set(args.matchup_id)
    case_ids = {row["matchup_id"] for row in csv.DictReader(args.cases.open())}
    if not requested <= case_ids:
        raise SystemExit("requested matchup is outside the locked development case list")
    rows = {row["matchup_id"]: row for row in csv.DictReader(MATCHUPS.open())}
    exit_code = 0
    for matchup_id in args.matchup_id:
        status = build_one(
            matchup_id,
            rows[matchup_id],
            pair_dir=args.pairs,
            output_dir=args.output,
            candidates=candidates,
            cutoff=args.history_cutoff,
            min_samples=args.min_samples,
            min_source_iqr=args.min_source_iqr,
            pair_maiac_aot_max=args.pair_maiac_aot_max,
            force=args.force,
        )
        print(
            "PAIR_AFFINE_HISTORY_DONE "
            + json.dumps(
                {
                    "status": status["status"],
                    "matchup_id": matchup_id,
                    "window_count": status.get("window_count"),
                    "scene_count": status.get("scene_count"),
                    "error_count": len(status.get("errors", [])),
                }
            ),
            flush=True,
        )
        if status["status"] != "OK":
            exit_code = 1
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
