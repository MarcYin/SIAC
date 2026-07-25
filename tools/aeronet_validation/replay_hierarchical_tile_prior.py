"""Replay saved SIAC cost cubes with a physical tile-level AOD hierarchy.

The tile posterior is derived only from the surface likelihood and the supplied
atmospheric prior.  A second local solve keeps the original prior uncertainty
but recentres it on that tile posterior.  AERONET truth is copied into the
output records for later scoring; it is never used by either retrieval.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from siac._rust_compat import surface_driven_pool_argmin

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_BASELINE = ROOT / "phaseD_results_lowcloud20_native_maiac_adaptive_b03_chi2_20260713"
DEFAULT_CUBES = ROOT / "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
DEFAULT_OUTPUT = ROOT / "analysis/medium_aod_hierarchical_tile_prior_development_20260713"


def _finite_mean(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.mean(finite)) if finite.size else None


def _within_ee(retrieved: float | None, truth: float) -> bool | None:
    if retrieved is None or not math.isfinite(retrieved):
        return None
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


def _tile_posterior(
    cube: np.ndarray,
    axis: np.ndarray,
    prior: np.ndarray,
    prior_unc: np.ndarray,
    valid: np.ndarray,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    usable = np.asarray(valid, dtype=bool) & np.isfinite(prior) & np.isfinite(prior_unc)
    if not np.any(usable):
        raise ValueError("cost cube has no pixels with a finite atmospheric prior")

    with np.errstate(invalid="ignore"):
        surface_curve = np.nanmedian(np.asarray(cube, dtype=np.float64)[:, usable], axis=1)
    prior_center = float(np.nanmedian(np.asarray(prior, dtype=np.float64)[usable]))
    prior_width = max(
        float(np.nanmedian(np.abs(np.asarray(prior_unc, dtype=np.float64)[usable]))),
        1.0e-6,
    )
    total_curve = surface_curve + np.square((axis - prior_center) / prior_width)
    finite = np.flatnonzero(np.isfinite(total_curve))
    if finite.size == 0:
        raise ValueError("tile posterior has no finite AOD node")
    best = int(finite[np.argmin(total_curve[finite])])

    # Delta-cost width is diagnostic only. The local replay deliberately keeps
    # the operational atmospheric-prior width and changes just its centre.
    basin = axis[np.isfinite(total_curve) & (total_curve <= total_curve[best] + 0.5)]
    width = (
        max(0.5 * float(np.max(basin) - np.min(basin)), 0.02)
        if basin.size
        else 0.02
    )
    return float(axis[best]), width, surface_curve, total_curve


def _local_replay(
    cube: np.ndarray,
    axis: np.ndarray,
    prior_unc: np.ndarray,
    valid: np.ndarray,
    *,
    tile_aod: float,
    pool_window: int,
    min_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    recentered = np.full(np.asarray(prior_unc).shape, tile_aod, dtype=np.float64)
    return surface_driven_pool_argmin(
        np.ascontiguousarray(cube, dtype=np.float64),
        np.ascontiguousarray(axis, dtype=np.float64),
        recentered,
        np.ascontiguousarray(prior_unc, dtype=np.float64),
        np.ascontiguousarray(valid, dtype=bool),
        int(pool_window),
        int(min_count),
    )


def _result_record(
    baseline: dict[str, Any],
    *,
    method: str,
    retrieved: float | None,
    tile_aod: float,
    tile_width: float,
    local_aod: np.ndarray | None,
    local_unc: np.ndarray | None,
    local_cost: np.ndarray | None,
    surface_curve: np.ndarray,
    total_curve: np.ndarray,
    axis: np.ndarray,
) -> dict[str, Any]:
    record = dict(baseline)
    truth = float(record["truth"])
    record["status"] = "OK" if retrieved is not None else "NO_RETRIEVAL"
    record["retrieved"] = retrieved
    record["scene_mean"] = retrieved
    record["err"] = retrieved - truth if retrieved is not None else None
    record["within_ee"] = _within_ee(retrieved, truth)
    record["retrieval_extraction"] = {
        "mode": "scene_mean",
        "selected": "scene_mean",
        "scene_mean": retrieved,
    }
    if local_aod is not None:
        record["valid_aot_count"] = int(np.count_nonzero(np.isfinite(local_aod)))
        record["valid_aot_fraction"] = float(np.mean(np.isfinite(local_aod)))

    finite_surface = np.flatnonzero(np.isfinite(surface_curve))
    surface_min = (
        float(axis[int(finite_surface[np.argmin(surface_curve[finite_surface])])])
        if finite_surface.size
        else None
    )
    record["hierarchical_tile_prior"] = {
        "method": method,
        "uses_aeronet_in_retrieval": False,
        "tile_posterior_aot": tile_aod,
        "tile_posterior_delta_cost_width": tile_width,
        "tile_surface_curve_min_aot": surface_min,
        "tile_surface_curve": surface_curve.tolist(),
        "tile_total_curve": total_curve.tolist(),
        "aot_axis": axis.tolist(),
        "local_aot_mean": _finite_mean(local_aod) if local_aod is not None else None,
        "local_aot_unc_mean": _finite_mean(local_unc) if local_unc is not None else None,
        "local_observation_cost_mean": _finite_mean(local_cost) if local_cost is not None else None,
    }
    return record


def replay_case(
    cube_path: Path,
    baseline_path: Path,
    output_root: Path,
) -> None:
    baseline = json.loads(baseline_path.read_text())
    with np.load(cube_path) as data:
        cube = np.asarray(data["cube"], dtype=np.float64)
        axis = np.asarray(data["aot_axis"], dtype=np.float64)
        prior = np.asarray(data["aot_prior"], dtype=np.float64)
        prior_unc = np.asarray(data["aot_prior_unc"], dtype=np.float64)
        valid = np.asarray(data["solve_valid"], dtype=bool)
        pool_window = int(np.asarray(data["pool_window"]).reshape(-1)[0])
        min_count = int(np.asarray(data["min_count"]).reshape(-1)[0])

    tile_aod, tile_width, surface_curve, total_curve = _tile_posterior(
        cube, axis, prior, prior_unc, valid
    )
    local_aod, local_unc, local_cost = _local_replay(
        cube,
        axis,
        prior_unc,
        valid,
        tile_aod=tile_aod,
        pool_window=pool_window,
        min_count=min_count,
    )

    records = {
        "tile_map": _result_record(
            baseline,
            method="tile_map",
            retrieved=tile_aod,
            tile_aod=tile_aod,
            tile_width=tile_width,
            local_aod=None,
            local_unc=None,
            local_cost=None,
            surface_curve=surface_curve,
            total_curve=total_curve,
            axis=axis,
        ),
        "hierarchical_recenter": _result_record(
            baseline,
            method="hierarchical_recenter",
            retrieved=_finite_mean(local_aod),
            tile_aod=tile_aod,
            tile_width=tile_width,
            local_aod=local_aod,
            local_unc=local_unc,
            local_cost=local_cost,
            surface_curve=surface_curve,
            total_curve=total_curve,
            axis=axis,
        ),
    }
    for method, record in records.items():
        out_dir = output_root / method
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / baseline_path.name).write_text(json.dumps(record, indent=2) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--cube-dir", type=Path, default=DEFAULT_CUBES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    completed = 0
    skipped = 0
    for cube_path in sorted(args.cube_dir.glob("*.npz")):
        baseline_path = args.baseline_dir / f"{cube_path.stem}.json"
        if not baseline_path.exists():
            skipped += 1
            continue
        replay_case(cube_path, baseline_path, args.output_dir)
        completed += 1
    print(f"replayed={completed} skipped_without_baseline={skipped} output={args.output_dir}")


if __name__ == "__main__":
    main()
