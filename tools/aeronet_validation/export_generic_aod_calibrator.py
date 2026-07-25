"""Fit, validate, and export the fixed low-cloud AOD calibration artifact."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import sklearn
from tools.aeronet_validation.aod_residual_calibration import (
    expected_error,
    feature_matrix,
    feature_schema,
    load_samples,
    metrics,
    site_fold,
)
from tools.aeronet_validation.select_generic_aod_calibrator import (
    AOD_ERROR_OFFSET,
    DEFAULT_MAX_SCENE_CLOUD_COVER,
    DEFAULT_ROOT,
    _ids,
    _model,
    candidate_specs,
    filter_scene_cloud_cover,
)

from siac.algorithms.aod_calibration import GenericAODCalibrator


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _variant_parameters(name: str) -> tuple[int, int]:
    match = re.fullmatch(r"seed(\d+)_n(\d+)", name)
    if match is None:
        raise ValueError(f"Primary variant is not a fixed seed/tree model: {name}")
    return int(match.group(1)), int(match.group(2))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--selection-analysis", type=Path, required=True)
    parser.add_argument("--offset-analysis", type=Path, required=True)
    parser.add_argument("--seed-analysis", type=Path, required=True)
    parser.add_argument("--train-result-dir", type=Path, required=True)
    parser.add_argument("--train-context-dir", type=Path, required=True)
    parser.add_argument("--train-atmo-context-dir", type=Path, required=True)
    parser.add_argument("--train-mids", type=Path, required=True)
    parser.add_argument("--target-result-dir", type=Path, required=True)
    parser.add_argument("--target-context-dir", type=Path, required=True)
    parser.add_argument("--target-atmo-context-dir", type=Path, required=True)
    parser.add_argument("--target-mids", type=Path, required=True)
    parser.add_argument("--model-output", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    args = parser.parse_args()

    selection = json.loads(args.selection_analysis.read_text(encoding="utf-8"))
    offset_analysis = json.loads(args.offset_analysis.read_text(encoding="utf-8"))
    seed_analysis = json.loads(args.seed_analysis.read_text(encoding="utf-8"))
    recipe_record = selection["selected_recipe"]
    spec = next(
        spec for spec in candidate_specs() if spec.name == recipe_record["candidate"]
    )
    if spec.target != "log_ratio" or float(recipe_record["shrinkage"]) != 1.0:
        raise ValueError("Exporter requires the selected unshrunk log-ratio recipe")
    primary_variant = str(seed_analysis["primary_stability_variant"])
    random_state, tree_count = _variant_parameters(primary_variant)

    matchup_path = args.root / "matchups" / "matchups.csv"
    common: dict[str, Any] = {"include_geography": False}
    train = filter_scene_cloud_cover(
        load_samples(
            args.train_result_dir,
            args.train_context_dir,
            matchup_path,
            _ids(args.train_mids),
            atmo_context_dir=args.train_atmo_context_dir,
            require_complete=True,
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
        require_complete=True,
        **common,
    )
    names = feature_schema(train)
    train_x = feature_matrix(train, names)
    baseline_train = np.asarray([sample.retrieved for sample in train], dtype=float)
    truth_train = np.asarray([sample.truth for sample in train], dtype=float)
    train_y = np.log(
        (truth_train + AOD_ERROR_OFFSET) / (baseline_train + AOD_ERROR_OFFSET)
    )
    estimator = _model(
        spec,
        len(names),
        random_state=random_state,
        tree_estimators=tree_count,
    )
    estimator.fit(train_x, train_y)
    selected_mask = estimator.named_steps["select"].get_support()
    selected_names = [
        name for name, selected in zip(names, selected_mask, strict=True) if selected
    ]
    importances = np.asarray(
        estimator.named_steps["model"].feature_importances_, dtype=float
    )
    feature_importance = sorted(
        (
            {"feature": name, "importance": float(importance)}
            for name, importance in zip(selected_names, importances, strict=True)
        ),
        key=lambda row: row["importance"],
        reverse=True,
    )
    calibrator = GenericAODCalibrator(
        feature_names=names,
        estimator=estimator,
        global_log_offset=float(offset_analysis["selected_offset"]),
        aod_shift=AOD_ERROR_OFFSET,
        metadata={
            "schema_version": "siac-generic-aod-calibrator-artifact-v1",
            "recipe": recipe_record["name"],
            "primary_variant": primary_variant,
            "cloud_domain": "metadata scene cloud cover < 20%",
            "geography_features": False,
        },
    )

    target_prediction = calibrator.predict(
        [sample.retrieved for sample in target],
        [sample.features for sample in target],
    )
    expected_ids = [row["matchup_id"] for row in offset_analysis["case_receipts"]]
    if expected_ids != [sample.matchup_id for sample in target]:
        raise ValueError("Target matchup IDs or order do not agree with validation receipts")
    expected_prediction = np.asarray(
        [row["adjusted"] for row in offset_analysis["case_receipts"]], dtype=float
    )
    max_abs_delta = float(np.max(np.abs(target_prediction - expected_prediction)))
    truth_target = np.asarray([sample.truth for sample in target], dtype=float)
    expected_hit = np.abs(expected_prediction - truth_target) <= expected_error(truth_target)
    artifact_hit = np.abs(target_prediction - truth_target) <= expected_error(truth_target)
    classification_disagreements = int(np.sum(expected_hit != artifact_hit))
    if max_abs_delta > 0.005 or classification_disagreements:
        raise ValueError(
            "Exported model does not reproduce validation within the cross-host "
            f"tolerance: max_delta={max_abs_delta} hit_disagreements="
            f"{classification_disagreements}"
        )
    target_folds = np.asarray([site_fold(sample.site) for sample in target], dtype=int)
    artifact_metrics = metrics(truth_target, target_prediction)
    artifact_confirmation_metrics = metrics(
        truth_target[target_folds == 4], target_prediction[target_folds == 4]
    )

    args.model_output.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(calibrator, args.model_output, compress=3)
    training_ids = "\n".join(sample.matchup_id for sample in train).encode("utf-8")
    manifest = {
        "schema_version": "siac-generic-aod-calibrator-manifest-v1",
        "model_path": str(args.model_output),
        "model_sha256": _sha256(args.model_output),
        "recipe": recipe_record,
        "global_log_offset": calibrator.global_log_offset,
        "formula": (
            "raw = clip((retrieved + 1/3) * exp(model_log_ratio) - 1/3, 0, 4); "
            "final = clip((raw + 1/3) * exp(-0.0125) - 1/3, 0, 4)"
        ),
        "random_state": random_state,
        "tree_count": tree_count,
        "training": {
            "count": len(train),
            "site_count": len({sample.site for sample in train}),
            "matchup_id_sha256": hashlib.sha256(training_ids).hexdigest(),
            "max_scene_cloud_cover_percent": DEFAULT_MAX_SCENE_CLOUD_COVER,
            "geography_features": False,
            "target_labels_used_for_fit": False,
        },
        "features": {
            "schema_count": len(names),
            "selected_count": len(selected_names),
            "schema": list(names),
            "selected": selected_names,
            "importance": feature_importance,
        },
        "validation_reproduction": {
            "target_count": len(target),
            "max_abs_prediction_delta": max_abs_delta,
            "prediction_tolerance": 0.005,
            "within_ee_classification_disagreements": classification_disagreements,
            "target": artifact_metrics,
            "target_confirmation": artifact_confirmation_metrics,
            "target_nested_offset_selection": offset_analysis[
                "target_nested_offset_selection"
            ]["aggregate"]["candidate"],
            "external_holdout": offset_analysis["external_holdout"]["candidate"],
            "all_replayed_variants_meet_target": offset_analysis[
                "all_replayed_variants_meet_target"
            ],
        },
        "software": {
            "numpy": np.__version__,
            "scikit_learn": sklearn.__version__,
            "joblib": joblib.__version__,
        },
        "source_reports": {
            "selection": str(args.selection_analysis),
            "offset": str(args.offset_analysis),
            "seed_robustness": str(args.seed_analysis),
        },
    }
    args.manifest_output.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_output.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "model": str(args.model_output),
                "model_sha256": manifest["model_sha256"],
                "train_count": len(train),
                "feature_count": len(names),
                "selected_feature_count": len(selected_names),
                "target_reproduction_max_abs_delta": max_abs_delta,
            },
            indent=2,
        )
    )
    print(args.manifest_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
