"""Replay saved surface likelihoods with a robust staged-MAIAC AOD prior.

This diagnostic changes only the probability model used for the atmospheric
backstop.  It keeps the staged-MAIAC centre, surface prior, RT coefficients, masks,
pooling, and AOT grid from the saved end-to-end run.  AERONET truth is copied
to the output for scoring but is never used by the solve.

The default is a Student-t prior with three degrees of freedom and the legacy
50 percent relative MAIAC scale.  Unlike the tuned AOD-regime backstop, this is
a single heavy-tailed error model: it behaves quadratically near MAIAC but lets
strong surface evidence escape a bad atmospheric forecast in either direction.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Literal

import numpy as np

from siac._rust_compat import surface_driven_pool_argmin

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_BASELINE = ROOT / (
    "phaseD_results_lowcloud20_mediumphysics_baseline_costcube_mediumdev_20260713"
)
DEFAULT_CUBES = ROOT / (
    "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
)
DEFAULT_OUTPUT = ROOT / (
    "analysis/medium_aod_robust_cams_student_t3_development_20260713"
)

PriorLoss = Literal["gaussian_linear", "student_t_linear", "gaussian_log", "student_t_log"]


def _finite_mean(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.mean(finite)) if finite.size else None


def _within_ee(retrieved: float | None, truth: float) -> bool | None:
    if retrieved is None or not math.isfinite(retrieved):
        return None
    return abs(retrieved - truth) <= 0.05 + 0.15 * truth


def prior_penalty(
    aot_axis: np.ndarray,
    *,
    centre: float,
    linear_scale: float | None,
    loss: PriorLoss,
    relative_scale: float,
    absolute_floor: float,
    log_factor: float,
    degrees_of_freedom: float,
) -> np.ndarray:
    """Return minus-two-log-prior values, up to an additive constant."""
    axis = np.asarray(aot_axis, dtype=np.float64)
    if loss.endswith("_linear"):
        scale = (
            max(float(linear_scale), absolute_floor)
            if linear_scale is not None
            else max(relative_scale * centre, absolute_floor)
        )
        standardized = (axis - centre) / scale
    else:
        positive_axis = np.maximum(axis, absolute_floor)
        positive_centre = max(centre, absolute_floor)
        standardized = np.log(positive_axis / positive_centre) / math.log(log_factor)

    if loss.startswith("gaussian_"):
        return np.square(standardized)
    return (degrees_of_freedom + 1.0) * np.log1p(
        np.square(standardized) / degrees_of_freedom
    )


def _build_record(
    baseline: dict[str, Any],
    *,
    solution: tuple[np.ndarray, np.ndarray, np.ndarray],
    metadata: dict[str, Any],
) -> dict[str, Any]:
    record = dict(baseline)
    retrieved = _finite_mean(solution[0])
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
    record["robust_atmo_prior_replay"] = {
        **metadata,
        "aot_unc_mean": _finite_mean(solution[1]),
        "penalized_cost_mean": _finite_mean(solution[2]),
    }
    return record


def replay_case(
    cube_path: Path,
    baseline_path: Path,
    output_dir: Path,
    *,
    loss: PriorLoss,
    relative_scale: float,
    absolute_floor: float,
    log_factor: float,
    degrees_of_freedom: float,
    linear_scale_source: Literal["relative", "saved"],
) -> None:
    baseline: dict[str, Any] = json.loads(baseline_path.read_text())
    with np.load(cube_path, allow_pickle=False) as data:
        cube = np.asarray(data["cube"], dtype=np.float64)
        axis = np.asarray(data["aot_axis"], dtype=np.float64)
        prior = np.asarray(data["aot_prior"], dtype=np.float64)
        saved_prior_unc = np.asarray(data["aot_prior_unc"], dtype=np.float64)
        valid = np.asarray(data["solve_valid"], dtype=bool)
        pool_window = int(np.asarray(data["pool_window"]).reshape(-1)[0])
        min_count = int(np.asarray(data["min_count"]).reshape(-1)[0])

    finite_prior = prior[np.isfinite(prior)]
    if finite_prior.size == 0:
        raise ValueError(f"No finite atmospheric prior in {cube_path}")
    centre = float(np.mean(finite_prior))
    if float(np.ptp(finite_prior)) > 1e-6:
        raise ValueError(
            f"{cube_path.name}: robust replay requires a tile-constant prior; "
            f"range={float(np.ptp(finite_prior)):.6g}"
        )

    saved_scale = _finite_mean(saved_prior_unc)
    if linear_scale_source == "saved" and saved_scale is None:
        raise ValueError(f"No finite saved atmospheric-prior uncertainty in {cube_path}")
    linear_scale = saved_scale if linear_scale_source == "saved" else None
    penalty = prior_penalty(
        axis,
        centre=centre,
        linear_scale=linear_scale,
        loss=loss,
        relative_scale=relative_scale,
        absolute_floor=absolute_floor,
        log_factor=log_factor,
        degrees_of_freedom=degrees_of_freedom,
    )
    penalized_cube = cube + penalty[:, None, None]
    disabled_backstop = np.full_like(prior, np.inf, dtype=np.float64)
    aod, aod_unc, cost = surface_driven_pool_argmin(
        np.ascontiguousarray(penalized_cube, dtype=np.float64),
        np.ascontiguousarray(axis, dtype=np.float64),
        np.ascontiguousarray(prior, dtype=np.float64),
        np.ascontiguousarray(disabled_backstop, dtype=np.float64),
        np.ascontiguousarray(valid, dtype=bool),
        pool_window,
        min_count,
    )
    metadata = {
        "uses_aeronet_in_retrieval": False,
        "shared_surface_prior_and_rt_cube": True,
        "atmospheric_prior_source": "CAMS",
        "atmospheric_prior_aot": centre,
        "loss": loss,
        "linear_scale_source": linear_scale_source if loss.endswith("_linear") else None,
        "linear_scale": linear_scale if loss.endswith("_linear") else None,
        "relative_scale": relative_scale if loss.endswith("_linear") else None,
        "absolute_floor": absolute_floor,
        "log_factor": log_factor if loss.endswith("_log") else None,
        "degrees_of_freedom": (
            degrees_of_freedom if loss.startswith("student_t_") else None
        ),
    }
    record = _build_record(
        baseline,
        solution=(
            np.asarray(aod, dtype=np.float64),
            np.asarray(aod_unc, dtype=np.float64),
            np.asarray(cost, dtype=np.float64),
        ),
        metadata=metadata,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / baseline_path.name).write_text(json.dumps(record, indent=2) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--cube-dir", type=Path, default=DEFAULT_CUBES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--loss",
        choices=("gaussian_linear", "student_t_linear", "gaussian_log", "student_t_log"),
        default="student_t_linear",
    )
    parser.add_argument("--relative-scale", type=float, default=0.5)
    parser.add_argument("--absolute-floor", type=float, default=0.02)
    parser.add_argument("--log-factor", type=float, default=2.0)
    parser.add_argument("--degrees-of-freedom", type=float, default=3.0)
    parser.add_argument(
        "--linear-scale-source",
        choices=("relative", "saved"),
        default="relative",
        help="Use the relative scale or the uncertainty saved with each cost cube.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.relative_scale <= 0.0:
        raise SystemExit("--relative-scale must be positive")
    if args.absolute_floor <= 0.0:
        raise SystemExit("--absolute-floor must be positive")
    if args.log_factor <= 1.0:
        raise SystemExit("--log-factor must exceed one")
    if args.degrees_of_freedom <= 0.0:
        raise SystemExit("--degrees-of-freedom must be positive")

    completed = 0
    skipped = 0
    for cube_path in sorted(args.cube_dir.glob("*.npz")):
        baseline_path = args.baseline_dir / f"{cube_path.stem}.json"
        if not baseline_path.exists():
            skipped += 1
            continue
        replay_case(
            cube_path,
            baseline_path,
            args.output_dir,
            loss=args.loss,
            relative_scale=args.relative_scale,
            absolute_floor=args.absolute_floor,
            log_factor=args.log_factor,
            degrees_of_freedom=args.degrees_of_freedom,
            linear_scale_source=args.linear_scale_source,
        )
        completed += 1
    summary = {
        "uses_aeronet_in_retrieval": False,
        "loss": args.loss,
        "relative_scale": args.relative_scale,
        "absolute_floor": args.absolute_floor,
        "log_factor": args.log_factor,
        "degrees_of_freedom": args.degrees_of_freedom,
        "linear_scale_source": args.linear_scale_source,
        "completed": completed,
        "skipped": skipped,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "replay_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
