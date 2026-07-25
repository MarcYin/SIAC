"""Build daily L2A histories with one learned correction per acquisition."""

from __future__ import annotations

import argparse
import csv
import json
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import joblib
import numpy as np
from tools.aeronet_validation.build_harmonized_l2a_histories import (
    DEFAULT_CASES,
    DEFAULT_OUTPUT,
    DEFAULT_PAIRS,
    HALF_SIZE_DEGREES,
    MATCHUPS,
    _fetch_scene,
    _load_scenes,
    _nanmedian,
    _write_history,
)
from tools.aeronet_validation.train_l2a_l1c_harmonizer import BAND_NAMES
from tools.aeronet_validation.train_scene_l2a_harmonizer import (
    build_scene_features,
    feature_names,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_SCENE_MODEL = (
    ROOT / "analysis/l2a_l1c_harmonizer_mediumdev_20260713/scene_model/scene_harmonizer.joblib"
)


def _candidate(value: str) -> tuple[str, float]:
    mode, separator, cap_text = value.partition(":")
    if not separator or mode not in {"all", "solver"}:
        raise argparse.ArgumentTypeError("candidate must be all:CAP or solver:CAP")
    try:
        cap = float(cap_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("candidate cap must be numeric") from exc
    if not 0.0 < cap <= 0.2:
        raise argparse.ArgumentTypeError("candidate cap must be in (0, 0.2]")
    return mode, cap


def _tag(candidate: tuple[str, float]) -> str:
    mode, cap = candidate
    return f"scene_et_{mode}_cap{cap:.3f}".replace(".", "p")


def _apply_scene_correction(
    fetched: dict[str, np.ndarray],
    scene: dict[str, Any],
    *,
    model: Any,
    candidate: tuple[str, float],
) -> tuple[np.ndarray, dict[str, Any]]:
    mode, cap = candidate
    surface = np.asarray(fetched["surface"], dtype=np.float64).copy()
    valid = np.asarray(fetched["valid"], dtype=bool)
    selected = np.flatnonzero(valid.ravel())
    if selected.size == 0:
        return surface.astype(np.float32), {"valid_pixels": 0, "correction": {}}
    flat = surface.reshape(len(BAND_NAMES), -1).T
    features = build_scene_features(
        flat[selected],
        fetched["l2a_aot"].ravel()[selected],
        fetched["l2a_tcwv"].ravel()[selected],
        scene,
    )
    correction = np.clip(np.asarray(model.predict(features[None])[0]), -cap, cap)
    indices = range(len(BAND_NAMES)) if mode == "all" else (1, 2, 3)
    for band_index in indices:
        flat[selected, band_index] = np.clip(
            flat[selected, band_index] + correction[band_index], 0.001, 0.8
        )
    return surface.astype(np.float32), {
        "valid_pixels": int(selected.size),
        "correction": {BAND_NAMES[index]: float(correction[index]) for index in indices},
    }


def _fetch_l2a_with_retries(
    ee: Any,
    grid: dict[str, Any],
    scene: dict[str, Any],
) -> dict[str, np.ndarray]:
    for attempt in range(1, 4):
        try:
            return _fetch_scene(ee, grid, scene)
        except Exception:  # noqa: BLE001 - bounded acquisition retry
            if attempt == 3:
                raise
            time.sleep(5.0 * attempt)
    raise RuntimeError("unreachable retry state")


def build_one(
    matchup_id: str,
    row: dict[str, str],
    *,
    pair_dir: Path,
    output_dir: Path,
    artifact: dict[str, Any],
    candidates: list[tuple[str, float]],
    cutoff: str,
    force: bool,
) -> dict[str, Any]:
    from bestpixel._gee import init_ee, utm_epsg_from_bbox, utm_grid

    outputs = {
        _tag(candidate): output_dir / _tag(candidate) / f"{matchup_id}.npz"
        for candidate in candidates
    }
    audit_path = output_dir / "scene_et_audits" / f"{matchup_id}.json"
    if audit_path.exists() and all(path.exists() for path in outputs.values()) and not force:
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

    fold = int(artifact["fold_by_matchup_id"][matchup_id])
    model = artifact["fold_models"][str(fold)]
    products: dict[str, list[np.ndarray]] = {tag: [] for tag in outputs}
    windows: list[str] = []
    scene_audits: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    ee = init_ee()
    for window in sorted(grouped):
        window_days: dict[str, list[np.ndarray]] = {tag: [] for tag in outputs}
        for day in sorted(grouped[window]):
            day_tiles: dict[str, list[np.ndarray]] = {tag: [] for tag in outputs}
            for scene in grouped[window][day]:
                try:
                    fetched = _fetch_l2a_with_retries(ee, grid, scene)
                except Exception as exc:  # noqa: BLE001 - retain other exact pairs
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
                    continue
                audit = {
                    "window": window,
                    "day": day,
                    "scene_id": scene.get("scene_id"),
                    "candidates": {},
                }
                for candidate in candidates:
                    tag = _tag(candidate)
                    corrected, stats = _apply_scene_correction(
                        fetched, scene, model=model, candidate=candidate
                    )
                    day_tiles[tag].append(corrected)
                    audit["candidates"][tag] = stats
                scene_audits.append(audit)
            if any(day_tiles.values()):
                for tag, tiles in day_tiles.items():
                    if tiles:
                        window_days[tag].append(_nanmedian(tiles))
        if all(window_days.values()):
            windows.append(window)
            for tag, days in window_days.items():
                products[tag].append(_nanmedian(days))

    status = "OK" if len(windows) >= 2 else "FAILED"
    if status == "OK":
        sensing = row["sensing_time_utc"]
        for tag, destination in outputs.items():
            _write_history(
                destination,
                products[tag],
                windows,
                grid=grid,
                scene_year=int(sensing[:4]),
                scene_month=int(sensing[5:7]),
                provenance={
                    "uses_aeronet": False,
                    "candidate": tag,
                    "history_cutoff": cutoff,
                    "model": artifact["best_model"],
                    "model_scope": "site-group crossfit",
                    "application": "one correction per acquisition before tile mosaic",
                },
            )
    audit = {
        "status": status,
        "matchup_id": matchup_id,
        "uses_aeronet": False,
        "fold": fold,
        "window_count": len(windows),
        "windows": windows,
        "scene_count": len(scene_audits),
        "scenes": scene_audits,
        "errors": errors,
        "outputs": {tag: str(path) for tag, path in outputs.items()},
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matchup_id", nargs="+")
    parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
    parser.add_argument("--model", type=Path, default=DEFAULT_SCENE_MODEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--history-cutoff", default="2023-12-31")
    parser.add_argument("--candidate", action="append", type=_candidate, default=[])
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    artifact = joblib.load(args.model)
    if artifact.get("uses_aeronet") is not False:
        raise SystemExit("refusing a scene model without uses_aeronet=false provenance")
    if artifact.get("feature_names") != feature_names():
        raise SystemExit("scene harmonizer feature schema mismatch")
    candidates = args.candidate or [("all", 0.03), ("solver", 0.03)]
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
            artifact=artifact,
            candidates=candidates,
            cutoff=args.history_cutoff,
            force=args.force,
        )
        print("SCENE_HISTORY_DONE " + json.dumps(audit), flush=True)
        if audit["status"] != "OK":
            exit_code = 1
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
