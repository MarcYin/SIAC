"""Audit operational feature dependence of the fixed AOD calibrator."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from tools.aeronet_validation.aod_residual_calibration import (
    CalibrationSample,
    feature_matrix,
    feature_schema,
    load_samples,
    site_fold,
)
from tools.aeronet_validation.evaluate_global_aod_offset import (
    apply_shifted_log_offset,
)
from tools.aeronet_validation.evaluate_selected_aod_seed_robustness import (
    evaluate_slice,
)
from tools.aeronet_validation.select_generic_aod_calibrator import (
    AOD_ERROR_OFFSET,
    DEFAULT_MAX_SCENE_CLOUD_COVER,
    DEFAULT_ROOT,
    RANDOM_STATE,
    _ids,
    _model,
    candidate_specs,
    filter_scene_cloud_cover,
)

TREE_COUNT = 300


def ablation_schemas(names: tuple[str, ...]) -> dict[str, tuple[str, ...]]:
    """Return fixed feature-removal audits without consulting target labels."""

    def keep(predicate: Any) -> tuple[str, ...]:
        return tuple(name for name in names if predicate(name))

    def is_cams_species(name: str) -> bool:
        return (
            name.startswith("context_cams_")
            and "total_aerosol" not in name
            and "angstrom" not in name
            and "ratio_" not in name
        )

    return {
        "full_operational_schema": names,
        "without_utc_phase": keep(lambda name: not name.startswith("time_utc_")),
        "without_all_time": keep(lambda name: not name.startswith("time_")),
        "without_cams_species_components": keep(lambda name: not is_cams_species(name)),
        "without_cams_context": keep(
            lambda name: not name.startswith("context_cams_")
            and not name.startswith("consistency_cams_")
        ),
        "without_maiac_context": keep(
            lambda name: not name.startswith("atmo_")
            and "_atmo" not in name
            and "atmo_" not in name
        ),
        "retrieval_surface_solver_only": keep(
            lambda name: not name.startswith(
                ("context_", "atmo_", "consistency_", "time_", "platform_")
            )
        ),
    }


def _predict(
    spec: Any,
    train: list[CalibrationSample],
    test: list[CalibrationSample],
    names: tuple[str, ...],
    offset: float,
) -> np.ndarray:
    train_x = feature_matrix(train, names)
    test_x = feature_matrix(test, names)
    baseline_train = np.asarray([sample.retrieved for sample in train], dtype=float)
    truth_train = np.asarray([sample.truth for sample in train], dtype=float)
    baseline_test = np.asarray([sample.retrieved for sample in test], dtype=float)
    target = np.log(
        (truth_train + AOD_ERROR_OFFSET) / (baseline_train + AOD_ERROR_OFFSET)
    )
    model = _model(
        spec,
        len(names),
        random_state=RANDOM_STATE,
        tree_estimators=TREE_COUNT,
    )
    model.fit(train_x, target)
    raw = np.clip(
        (baseline_test + AOD_ERROR_OFFSET) * np.exp(model.predict(test_x))
        - AOD_ERROR_OFFSET,
        0.0,
        4.0,
    )
    return apply_shifted_log_offset(raw, offset)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--selection-analysis", type=Path, required=True)
    parser.add_argument("--offset-analysis", type=Path, required=True)
    parser.add_argument("--train-result-dir", type=Path, required=True)
    parser.add_argument("--train-context-dir", type=Path, required=True)
    parser.add_argument("--train-atmo-context-dir", type=Path, required=True)
    parser.add_argument("--train-mids", type=Path, required=True)
    parser.add_argument("--target-result-dir", type=Path, required=True)
    parser.add_argument("--target-context-dir", type=Path, required=True)
    parser.add_argument("--target-atmo-context-dir", type=Path, required=True)
    parser.add_argument("--target-mids", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    selection = json.loads(args.selection_analysis.read_text(encoding="utf-8"))
    offset = float(
        json.loads(args.offset_analysis.read_text(encoding="utf-8"))["selected_offset"]
    )
    spec = next(
        spec
        for spec in candidate_specs()
        if spec.name == selection["selected_recipe"]["candidate"]
    )
    matchup_path = args.root / "matchups" / "matchups.csv"
    common: dict[str, Any] = {"include_geography": False}
    train = filter_scene_cloud_cover(
        load_samples(
            args.train_result_dir,
            args.train_context_dir,
            matchup_path,
            _ids(args.train_mids),
            atmo_context_dir=args.train_atmo_context_dir,
            **common,
        ),
        DEFAULT_MAX_SCENE_CLOUD_COVER,
    )
    target = load_samples(
        args.target_result_dir,
        args.target_context_dir,
        matchup_path,
        _ids(args.target_mids),
        atmo_context_dir=args.target_atmo_context_dir,
        **common,
    )
    external_development = [sample for sample in train if site_fold(sample.site) != 4]
    external_holdout = [sample for sample in train if site_fold(sample.site) == 4]
    external_sites = {sample.site for sample in train}
    folds = np.asarray([site_fold(sample.site) for sample in target], dtype=int)
    seen = np.asarray([sample.site in external_sites for sample in target], dtype=bool)
    variants: dict[str, Any] = {}
    for name, names in ablation_schemas(feature_schema(train)).items():
        target_prediction = _predict(spec, train, target, names, offset)
        holdout_prediction = _predict(
            spec, external_development, external_holdout, names, offset
        )
        variants[name] = {
            "feature_count": len(names),
            "target_all": evaluate_slice(
                target, target_prediction, np.arange(len(target))
            ),
            "target_development": evaluate_slice(
                target, target_prediction, np.flatnonzero(folds != 4)
            ),
            "target_confirmation": evaluate_slice(
                target, target_prediction, np.flatnonzero(folds == 4)
            ),
            "target_seen_sites": evaluate_slice(
                target, target_prediction, np.flatnonzero(seen)
            ),
            "target_unseen_sites": evaluate_slice(
                target, target_prediction, np.flatnonzero(~seen)
            ),
            "external_holdout": evaluate_slice(
                external_holdout,
                holdout_prediction,
                np.arange(len(external_holdout)),
            ),
        }
    analysis = {
        "schema_version": "siac-generic-aod-feature-ablation-v1",
        "policy": {
            "diagnostic_only": True,
            "target_labels_used_to_select_ablation": False,
            "base_recipe": selection["selected_recipe"]["name"],
            "global_log_offset": offset,
            "tree_count": TREE_COUNT,
            "random_state": RANDOM_STATE,
            "geography_features": False,
        },
        "variants": variants,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                name: {
                    "features": row["feature_count"],
                    "target": row["target_all"]["candidate"],
                    "confirmation": row["target_confirmation"]["candidate"],
                    "unseen_sites": row["target_unseen_sites"]["candidate"],
                    "external_holdout": row["external_holdout"]["candidate"],
                }
                for name, row in variants.items()
            },
            indent=2,
        )
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
