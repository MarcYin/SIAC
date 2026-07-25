"""Re-solve saved surface cost cubes under sparse-support pooling variants."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from siac._rust_compat import surface_driven_pool_argmin
from siac.algorithms.solver.surface_driven import _select_auto2_solution

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
REPORT = Path("notes/reports/masked-r2-regression-20260711/comparison.json")
ROWS = REPORT.with_name("comparison-rows.csv")
CUBES = ROOT / "phaseD_cost_cubes_campaign250_masked_r2_l2awvp_6s_20260710"

POOL_VARIANTS = (
    ("current_20x20_min80", 20, 80),
    ("same_window_min20", 20, 20),
    ("same_window_min5", 20, 5),
    ("half_window_min20", 10, 20),
    ("half_window_min5", 10, 5),
    ("per_pixel", 1, 1),
)


def _native_pool_argmin(
    cost_cube: np.ndarray,
    aot_axis: np.ndarray,
    aot_prior: np.ndarray,
    aot_prior_unc: np.ndarray,
    valid_mask: np.ndarray,
    pool_window: int,
    pool_min_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return surface_driven_pool_argmin(
        np.ascontiguousarray(cost_cube, dtype=np.float64),
        np.ascontiguousarray(aot_axis, dtype=np.float64),
        np.ascontiguousarray(aot_prior, dtype=np.float64),
        np.ascontiguousarray(aot_prior_unc, dtype=np.float64),
        np.ascontiguousarray(valid_mask, dtype=bool),
        int(pool_window),
        int(pool_min_count),
    )


def _numpy_pool_argmin(
    cost_cube: np.ndarray,
    aot_axis: np.ndarray,
    aot_prior: np.ndarray,
    aot_prior_unc: np.ndarray,
    valid_mask: np.ndarray,
    pool_window: int,
    pool_min_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Serial equivalent of the Rust surface-driven pooling kernel.

    This script is a diagnostic replay over saved cost cubes. On busy shared
    nodes, the native Rayon pool can fail to create a worker thread even with
    RAYON_NUM_THREADS=1. This slower path is available for process-light
    one-off diagnostics when needed.
    """
    cost = np.asarray(cost_cube, dtype=np.float64)
    axis = np.asarray(aot_axis, dtype=np.float64)
    prior = np.asarray(aot_prior, dtype=np.float64)
    prior_unc = np.abs(np.asarray(aot_prior_unc, dtype=np.float64))
    valid = np.asarray(valid_mask, dtype=bool)
    n_aot, ny, nx = cost.shape
    w = max(1, int(pool_window))
    lo = (w - 1) // 2
    hi = w // 2
    min_count = max(1, int(pool_min_count))

    pooled = np.full((n_aot, ny, nx), np.nan, dtype=np.float64)
    try:
        from scipy import ndimage
    except ImportError:  # pragma: no cover - scipy is present in the pixi env
        for k in range(n_aot):
            slab = cost[k]
            for iy in range(ny):
                y0 = max(0, iy - lo)
                y1 = min(ny, iy + hi + 1)
                for ix in range(nx):
                    x0 = max(0, ix - lo)
                    x1 = min(nx, ix + hi + 1)
                    window = slab[y0:y1, x0:x1]
                    finite = window[np.isfinite(window)]
                    if finite.size >= min_count:
                        pooled[k, iy, ix] = float(np.median(finite))
    else:
        origin = -1 if w % 2 == 0 else 0

        def window_median(window: np.ndarray, *, axis: int | tuple[int, ...]) -> np.ndarray:
            count = np.sum(np.isfinite(window), axis=axis)
            with np.errstate(all="ignore"):
                median = np.nanmedian(window, axis=axis)
            return np.where(count >= min_count, median, np.nan)

        for k in range(n_aot):
            pooled[k] = ndimage.vectorized_filter(
                cost[k],
                window_median,
                size=(w, w),
                mode="constant",
                cval=np.nan,
                origin=origin,
                batch_memory=256 * 1024 * 1024,
            )

    aod_out = np.full((ny, nx), np.nan, dtype=np.float64)
    unc_out = np.full((ny, nx), np.nan, dtype=np.float64)
    jmin_out = np.full((ny, nx), np.nan, dtype=np.float64)
    all_finite = np.isfinite(pooled).all(axis=0)
    for iy in range(ny):
        for ix in range(nx):
            if not valid[iy, ix] or not all_finite[iy, ix]:
                continue
            pooled_curve = pooled[:, iy, ix]
            total = pooled_curve.copy()
            if (
                math.isfinite(float(prior[iy, ix]))
                and math.isfinite(float(prior_unc[iy, ix]))
                and prior_unc[iy, ix] > 0.0
            ):
                total += np.square(axis - prior[iy, ix]) / max(prior_unc[iy, ix] ** 2, 1e-12)
            finite = np.isfinite(total)
            if not np.any(finite):
                continue
            best_k = int(np.nanargmin(total))
            best = float(total[best_k])
            aod_out[iy, ix] = float(axis[best_k])
            jmin_out[iy, ix] = float(np.nanmin(pooled_curve))
            basin = axis[np.isfinite(total) & (total <= best + 0.5)]
            unc_out[iy, ix] = max(0.5 * (float(np.max(basin)) - float(np.min(basin))), 0.02)
    return aod_out, unc_out, jmin_out


def _finite_mean(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.mean(finite)) if finite.size else None


def _within_ee(retrieved: float | None, truth: float) -> bool | None:
    if retrieved is None or not math.isfinite(retrieved):
        return None
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


def _pool_solve(
    *,
    cube: np.ndarray,
    cube_abs: np.ndarray,
    aot_axis: np.ndarray,
    aot_prior: np.ndarray,
    aot_prior_unc: np.ndarray,
    solve_valid: np.ndarray,
    window: int,
    min_count: int,
    engine: str,
) -> tuple[np.ndarray, np.ndarray]:
    pool_argmin = _native_pool_argmin if engine == "native" else _numpy_pool_argmin
    common = (
        aot_axis,
        aot_prior,
        aot_prior_unc,
        solve_valid,
        int(window),
        int(min_count),
    )
    main_aod, main_unc, main_cost = pool_argmin(cube, *common)
    abs_aod, abs_unc, abs_cost = pool_argmin(cube_abs, *common)
    selected, _, cost, _ = _select_auto2_solution(
        aot_main=np.asarray(main_aod),
        unc_main=np.asarray(main_unc),
        cost_main=np.asarray(main_cost),
        aot_abs=np.asarray(abs_aod),
        unc_abs=np.asarray(abs_unc),
        cost_abs=np.asarray(abs_cost),
        clean_threshold=0.15,
        high_threshold=0.6,
        cost_gain=0.2,
    )
    return np.asarray(selected), np.asarray(cost)


def _scene_curve_solve(
    *,
    cube: np.ndarray,
    cube_abs: np.ndarray,
    aot_axis: np.ndarray,
    aot_prior: np.ndarray,
    aot_prior_unc: np.ndarray,
    solve_valid: np.ndarray,
) -> tuple[float | None, float | None]:
    valid = np.asarray(solve_valid, dtype=bool)
    if not np.any(valid):
        return None, None
    prior = float(np.nanmedian(np.asarray(aot_prior, dtype=np.float64)[valid]))
    prior_unc = float(np.nanmedian(np.asarray(aot_prior_unc, dtype=np.float64)[valid]))
    prior_unc = max(prior_unc, 1.0e-6)
    axis = np.asarray(aot_axis, dtype=np.float64)

    def resolve(values: np.ndarray) -> tuple[float, float]:
        curve = np.nanmedian(np.asarray(values, dtype=np.float64)[:, valid], axis=1)
        total = curve + np.square((axis - prior) / prior_unc)
        finite = np.flatnonzero(np.isfinite(total))
        if not finite.size:
            return math.nan, math.nan
        index = int(finite[np.argmin(total[finite])])
        return float(axis[index]), float(total[index])

    main_aod, main_cost = resolve(cube)
    abs_aod, abs_cost = resolve(cube_abs)
    selected, _, selected_cost, _ = _select_auto2_solution(
        aot_main=np.asarray([main_aod]),
        unc_main=np.asarray([prior_unc]),
        cost_main=np.asarray([main_cost]),
        aot_abs=np.asarray([abs_aod]),
        unc_abs=np.asarray([prior_unc]),
        cost_abs=np.asarray([abs_cost]),
        clean_threshold=0.15,
        high_threshold=0.6,
        cost_gain=0.2,
    )
    return _finite_mean(selected), _finite_mean(selected_cost)


def analyze_site(row: dict[str, Any], cube_dir: Path, *, engine: str) -> dict[str, Any]:
    matchup_id = str(row["matchup_id"])
    truth = float(row["truth"])
    cube_path = cube_dir / f"{matchup_id}.npz"
    if not cube_path.exists():
        return {
            "matchup_id": matchup_id,
            "site": row["site"],
            "truth": truth,
            "status": "MISSING_CUBE",
            "variants": [],
        }
    with np.load(cube_path) as data:
        cube = np.asarray(data["cube"], dtype=np.float32)
        cube_abs = np.asarray(data["cube_abs"], dtype=np.float32)
        aot_axis = np.asarray(data["aot_axis"], dtype=np.float64)
        aot_prior = np.asarray(data["aot_prior"], dtype=np.float64)
        aot_prior_unc = np.asarray(data["aot_prior_unc"], dtype=np.float64)
        solve_valid = np.asarray(data["solve_valid"], dtype=bool)
        saved_window = int(np.asarray(data["pool_window"]).reshape(-1)[0])
        saved_min_count = int(np.asarray(data["min_count"]).reshape(-1)[0])

    variants = []
    for name, window, min_count in POOL_VARIANTS:
        aod, cost = _pool_solve(
            cube=cube,
            cube_abs=cube_abs,
            aot_axis=aot_axis,
            aot_prior=aot_prior,
            aot_prior_unc=aot_prior_unc,
            solve_valid=solve_valid,
            window=window,
            min_count=min_count,
            engine=engine,
        )
        retrieved = _finite_mean(aod)
        variants.append(
            {
                "name": name,
                "window": window,
                "min_count": min_count,
                "retrieved": retrieved,
                "within_ee": _within_ee(retrieved, truth),
                "solved_pixels": int(np.count_nonzero(np.isfinite(aod))),
                "median_cost": (
                    float(np.nanmedian(cost[np.isfinite(cost)]))
                    if np.any(np.isfinite(cost))
                    else None
                ),
            }
        )

    scene_aod, scene_cost = _scene_curve_solve(
        cube=cube,
        cube_abs=cube_abs,
        aot_axis=aot_axis,
        aot_prior=aot_prior,
        aot_prior_unc=aot_prior_unc,
        solve_valid=solve_valid,
    )
    variants.append(
        {
            "name": "scene_curve",
            "window": None,
            "min_count": int(np.count_nonzero(solve_valid)),
            "retrieved": scene_aod,
            "within_ee": _within_ee(scene_aod, truth),
            "solved_pixels": int(np.count_nonzero(solve_valid)) if scene_aod is not None else 0,
            "median_cost": scene_cost,
        }
    )
    return {
        "matchup_id": matchup_id,
        "site": row["site"],
        "truth": truth,
        "status": "OK",
        "cloud_fraction": row["cloud_frac"],
        "valid_support_count": int(np.count_nonzero(solve_valid)),
        "saved_pool_window": saved_window,
        "saved_min_count": saved_min_count,
        "variants": variants,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, default=REPORT)
    parser.add_argument("--rows-csv", type=Path, default=ROWS)
    parser.add_argument("--cube-dir", type=Path, default=CUBES)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cohort", choices=("sparse", "all"), default="sparse")
    parser.add_argument("--index", type=int, required=True, help="One-based cohort index")
    parser.add_argument(
        "--engine",
        choices=("numpy", "native"),
        default="native",
        help="Pooling backend. native is fast; numpy is serial and avoids Rayon thread creation.",
    )
    args = parser.parse_args()

    comparison = json.loads(args.comparison.read_text(encoding="utf-8"))
    if args.cohort == "sparse":
        rows = comparison["new_abstention_diagnostics"]["rows"]
        cohort = [
            row for row in rows if row["abstention_class"] == "valid_support_but_pool_unsolved"
        ]
    else:
        with args.rows_csv.open(encoding="utf-8", newline="") as handle:
            cohort = list(csv.DictReader(handle))
    if not 1 <= args.index <= len(cohort):
        raise SystemExit(f"index {args.index} outside 1..{len(cohort)}")
    result = analyze_site(cohort[args.index - 1], args.cube_dir, engine=args.engine)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    path = args.output_dir / f"{result['matchup_id']}.json"
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
