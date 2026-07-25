"""Replay one S2 surface likelihood with local historical OOF error maps.

The replay keeps the saved surface source, RT coefficients, masks, AOD grid,
and atmospheric prior.  It derives a local visible-surface bias and uncertainty
map only from leave-one-realization-out errors in the historical S2 composite
stack.  AERONET truth is copied to receipts for later scoring and is never used
to build the maps or choose a result.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from tools.aeronet_validation.replay_cost_cube_prior_uncertainty import (
    _build_record,
    _finite_mean,
    _solve,
    _surface_curve_min_index,
)
from tools.aeronet_validation.replay_surface_oof_uncertainty import (
    estimate_temporal_oof_error_maps,
)

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
DEFAULT_BASELINE = ROOT / (
    "phaseD_results_lowcloud20_mediumphysics_baseline_costcube_mediumdev_20260713"
)
DEFAULT_CUBES = ROOT / (
    "phaseD_cost_cubes_lowcloud20_mediumphysics_baseline_mediumdev_20260713"
)
DEFAULT_DUMPS = ROOT / "calib_dumps_c250"
DEFAULT_FIELDS = ROOT / (
    "analysis/medium_aod_surface_oof_map_et20_fields_development_20260713"
)
VISIBLE_DEBIAS = {
    "B02": (-0.0003, 0.0243),
    "B03": (-0.0006, 0.0235),
    "B04": (-0.0011, 0.0223),
}


def _default_output(mode: str, prior: str) -> Path:
    return ROOT / (
        f"analysis/medium_aod_surface_oof_map_et20_{mode}_{prior}_development_20260713"
    )


def _solutions(
    cube: np.ndarray,
    axis: np.ndarray,
    prior: np.ndarray,
    calibrated_unc: np.ndarray,
    valid: np.ndarray,
    pool_window: int,
    min_count: int,
    *,
    z_threshold: float,
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], dict[str, Any]]:
    product_unc = np.maximum(0.5 * prior, 0.02)
    calibrated = _solve(
        cube, axis, prior, calibrated_unc, valid, pool_window, min_count
    )
    product = _solve(cube, axis, prior, product_unc, valid, pool_window, min_count)
    surface_index = _surface_curve_min_index(cube)
    surface_min = float(axis[surface_index]) if surface_index is not None else None
    atmo_aot = _finite_mean(prior)
    calibrated_sigma = _finite_mean(calibrated_unc)
    conflict = None
    if (
        surface_min is not None
        and atmo_aot is not None
        and calibrated_sigma is not None
        and calibrated_sigma > 0.0
    ):
        conflict = (surface_min - atmo_aot) / calibrated_sigma
    use_product = bool(conflict is not None and conflict > z_threshold)
    return (
        {
            "calibrated": calibrated,
            "product": product,
            "adaptive": product if use_product else calibrated,
        },
        {
            "surface_curve_min_aot": surface_min,
            "atmospheric_prior_aot_mean": atmo_aot,
            "calibrated_prior_sigma_mean": calibrated_sigma,
            "standardized_positive_conflict": conflict,
            "z_threshold": z_threshold,
            "adaptive_decision": "product" if use_product else "calibrated",
        },
    )


def replay_case(
    cube_path: Path,
    dump_path: Path,
    baseline_path: Path,
    output_dirs: dict[tuple[str, str], Path],
    fields_dir: Path,
    *,
    z_threshold: float,
) -> None:
    baseline: dict[str, Any] = json.loads(baseline_path.read_text())
    with np.load(cube_path, allow_pickle=False) as data:
        signed_residual = np.asarray(data["band_signed_residual_cube"], dtype=np.float64)
        base_uncertainty = np.asarray(data["boa_unc"], dtype=np.float64)
        axis = np.asarray(data["aot_axis"], dtype=np.float64)
        prior = np.asarray(data["aot_prior"], dtype=np.float64)
        calibrated_unc = np.asarray(data["aot_prior_unc"], dtype=np.float64)
        valid = np.asarray(data["solve_valid"], dtype=bool)
        band_names = [str(value) for value in np.asarray(data["band_names"]).tolist()]
        pool_window = int(np.asarray(data["pool_window"]).reshape(-1)[0])
        min_count = int(np.asarray(data["min_count"]).reshape(-1)[0])
    with np.load(dump_path, allow_pickle=False) as data:
        composites = np.asarray(data["comp"], dtype=np.float64)

    field_path = fields_dir / cube_path.name
    if field_path.exists():
        with np.load(field_path, allow_pickle=False) as fields:
            cached_names = [
                str(value) for value in np.asarray(fields["band_names"]).tolist()
            ]
            if cached_names != band_names:
                raise ValueError(
                    f"cached OOF bands {cached_names} do not match cost bands {band_names}"
                )
            bias = np.asarray(fields["bias"], dtype=np.float64)
            oof_sigma = np.asarray(fields["sigma"], dtype=np.float64)
            support = np.asarray(fields["support"], dtype=np.int64)
        oof_metadata = {
            "estimator": "cached leave-one-realization-out ET20 local signed error and MAD",
            "band_bias_median": {
                name: float(np.nanmedian(bias[index]))
                for index, name in enumerate(band_names)
            },
            "band_sigma_median": {
                name: float(np.nanmedian(oof_sigma[index]))
                for index, name in enumerate(band_names)
            },
            "band_sigma_p90": {
                name: float(np.nanpercentile(oof_sigma[index], 90.0))
                for index, name in enumerate(band_names)
            },
            "band_support_median": {
                name: float(np.nanmedian(support[index]))
                for index, name in enumerate(band_names)
            },
        }
    else:
        bias, oof_sigma, support, oof_metadata = estimate_temporal_oof_error_maps(
            composites, band_names
        )
    if bias.shape != base_uncertainty.shape:
        raise ValueError(
            f"OOF map shape {bias.shape} does not match cost field {base_uncertainty.shape}"
        )
    effective_uncertainty = np.sqrt(
        np.square(base_uncertainty) + np.square(oof_sigma)
    )

    anchor_iteration = baseline.get("anchor_iterate") or {}
    final_surface_anchor_aot = anchor_iteration.get("pass1_scene_mean")
    if final_surface_anchor_aot is None:
        final_surface_anchor_aot = baseline.get("atmo_prior", {}).get("aot_median")
    final_surface_anchor_aot = float(final_surface_anchor_aot)
    debias_offsets = np.asarray(
        [
            VISIBLE_DEBIAS.get(name, (0.0, 0.0))[0]
            + VISIBLE_DEBIAS.get(name, (0.0, 0.0))[1] * final_surface_anchor_aot
            for name in band_names
        ],
        dtype=np.float64,
    )[:, np.newaxis, np.newaxis]

    residual_by_mode = {
        "unc": signed_residual,
        "bias_unc": signed_residual + bias[:, np.newaxis],
        # The saved prior contains the AERONET-site surface debias. Adding its
        # offset back to BOA-prior residuals reconstructs the coefficient-free
        # historical-S2 prediction before the OOF correction is evaluated.
        "nodebias_unc": signed_residual + debias_offsets[:, np.newaxis],
        "oof_corrected_unc": (
            signed_residual
            + debias_offsets[:, np.newaxis]
            + bias[:, np.newaxis]
        ),
    }
    if not field_path.exists():
        fields_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            field_path,
            bias=bias.astype(np.float32),
            sigma=oof_sigma.astype(np.float32),
            support=support.astype(np.int16),
            band_names=np.asarray(band_names),
        )

    for mode, residual in residual_by_mode.items():
        with np.errstate(invalid="ignore", divide="ignore"):
            cube = np.sum(
                np.square(residual)
                / np.square(effective_uncertainty[:, np.newaxis]),
                axis=0,
            )
        solutions, solve_metadata = _solutions(
            cube,
            axis,
            prior,
            calibrated_unc,
            valid,
            pool_window,
            min_count,
            z_threshold=z_threshold,
        )
        for prior_mode, solution in solutions.items():
            record = _build_record(
                baseline,
                solution=solution,
                metadata={
                    "uses_aeronet_in_retrieval": False,
                    "surface_source": "historical_s2_bestpixel_only",
                    "shared_surface_prior_and_rt": True,
                    "map_mode": mode,
                    "prior_mode": prior_mode,
                    "removed_aeronet_surface_debias": mode
                    in {"nodebias_unc", "oof_corrected_unc"},
                    "surface_anchor_aot": final_surface_anchor_aot,
                    "removed_surface_debias": {
                        name: float(debias_offsets[index, 0, 0])
                        for index, name in enumerate(band_names)
                    },
                    "oof": oof_metadata,
                    "base_uncertainty_mean": {
                        name: _finite_mean(base_uncertainty[index])
                        for index, name in enumerate(band_names)
                    },
                    "effective_uncertainty_mean": {
                        name: _finite_mean(effective_uncertainty[index])
                        for index, name in enumerate(band_names)
                    },
                    **solve_metadata,
                },
            )
            record["surface_oof_map"] = record.pop("same_cube_prior_uncertainty")
            destination = output_dirs[(mode, prior_mode)]
            destination.mkdir(parents=True, exist_ok=True)
            (destination / baseline_path.name).write_text(
                json.dumps(record, indent=2) + "\n"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--cube-dir", type=Path, default=DEFAULT_CUBES)
    parser.add_argument("--dump-dir", type=Path, default=DEFAULT_DUMPS)
    parser.add_argument("--fields-dir", type=Path, default=DEFAULT_FIELDS)
    parser.add_argument("--z-threshold", type=float, default=2.576)
    parser.add_argument("--matchup-id", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.z_threshold <= 0.0:
        raise SystemExit("--z-threshold must be positive")
    output_dirs = {
        (mode, prior): _default_output(mode, prior)
        for mode in ("unc", "bias_unc", "nodebias_unc", "oof_corrected_unc")
        for prior in ("calibrated", "product", "adaptive")
    }
    selected = set(args.matchup_id)
    completed = skipped = 0
    for cube_path in sorted(args.cube_dir.glob("*.npz")):
        if selected and cube_path.stem not in selected:
            continue
        baseline_path = args.baseline_dir / f"{cube_path.stem}.json"
        dump_path = args.dump_dir / f"{cube_path.stem}.npz"
        if not baseline_path.exists() or not dump_path.exists():
            skipped += 1
            continue
        replay_case(
            cube_path,
            dump_path,
            baseline_path,
            output_dirs,
            args.fields_dir,
            z_threshold=args.z_threshold,
        )
        completed += 1
    print(f"replayed={completed} skipped={skipped} fields={args.fields_dir}")


if __name__ == "__main__":
    main()
