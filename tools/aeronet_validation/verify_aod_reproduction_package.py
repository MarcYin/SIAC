"""Verify and replay the frozen SIAC AOD production reproduction package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping

import joblib
import numpy as np


DEFAULT_PACKAGE = Path(
    "/gws/ssde/j25a/nceo_isp/public/siac_refactor/reports/"
    "aod-production-reproduction-spec-20260713/downloads"
)
DEFAULT_REPO = Path("/home/users/marcyin/SIAC")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_jsonl(path: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip():
                continue
            record = json.loads(line)
            matchup_id = str(record["matchup_id"])
            if matchup_id in records:
                raise ValueError(f"Duplicate matchup_id in {path}: {matchup_id}")
            records[matchup_id] = record
    return records


def _ids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _matchups(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return {str(row["matchup_id"]): row for row in csv.DictReader(stream)}


def _references(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return {str(row["matchup_id"]): row for row in csv.DictReader(stream)}


def _validate_inventory(package: Path) -> tuple[int, str]:
    report_root = package.parent
    manifest_path = package / "reproduction-manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    failures = []
    for relative, expected in manifest["files"].items():
        path = report_root / relative
        if not path.exists():
            failures.append(f"missing: {relative}")
            continue
        if path.stat().st_size != int(expected["bytes"]):
            failures.append(f"size: {relative}")
            continue
        if _sha256(path) != str(expected["sha256"]):
            failures.append(f"sha256: {relative}")
    if failures:
        raise ValueError("Package inventory failed: " + ", ".join(failures))
    return len(manifest["files"]), str(manifest["model_sha256"])


def _ordered(records: Mapping[str, Any], ids: Iterable[str], label: str) -> list[Any]:
    missing = [matchup_id for matchup_id in ids if matchup_id not in records]
    if missing:
        raise ValueError(f"{label} is missing {len(missing)} target records")
    return [records[matchup_id] for matchup_id in ids]


def verify(package: Path, repo_root: Path) -> dict[str, Any]:
    inventory_count, expected_model_hash = _validate_inventory(package)
    model_path = package / "aod-calibrator.joblib"
    if _sha256(model_path) != expected_model_hash:
        raise ValueError("Model hash does not match reproduction manifest")

    for path in (repo_root, repo_root / "python"):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))
    from tools.aeronet_validation.aod_residual_calibration import (
        extract_operational_features,
    )

    target_ids = _ids(package / "campaign250-lowcloud20-mids.txt")
    if len(target_ids) != 152 or len(set(target_ids)) != 152:
        raise ValueError("Target denominator is not 152 unique matchup IDs")
    result_by_id = _load_jsonl(package / "target-physical-results.jsonl")
    maiac_by_id = _load_jsonl(package / "target-maiac-inputs.jsonl")
    cams_by_id = _load_jsonl(package / "target-cams-context.jsonl")
    matchup_by_id = _matchups(package / "target-matchups.csv")
    reference_by_id = _references(package / "reference-predictions.csv")
    results = _ordered(result_by_id, target_ids, "physical results")
    maiac = _ordered(maiac_by_id, target_ids, "MAIAC inputs")
    cams = _ordered(cams_by_id, target_ids, "CAMS context")
    matchups = _ordered(matchup_by_id, target_ids, "matchups")
    references = _ordered(reference_by_id, target_ids, "reference predictions")

    if any(result.get("status") != "OK" for result in results):
        raise ValueError("Not every physical result is OK")
    feature_maps = [
        extract_operational_features(
            result,
            context,
            matchup,
            atmo_context=atmo,
            include_geography=False,
        )
        for result, context, matchup, atmo in zip(
            results, cams, matchups, maiac, strict=True
        )
    ]
    baseline = np.asarray([float(result["retrieved"]) for result in results], dtype=float)
    truth = np.asarray([float(result["truth"]) for result in results], dtype=float)
    expected = np.asarray([float(row["adjusted"]) for row in references], dtype=float)
    calibrator = joblib.load(model_path)
    prediction = np.asarray(calibrator.predict(baseline, feature_maps), dtype=float)
    if prediction.shape != (152,) or not np.all(np.isfinite(prediction)):
        raise ValueError("Calibrator output is not 152 finite scalar predictions")

    max_abs_delta = float(np.max(np.abs(prediction - expected)))
    threshold = 0.05 + 0.15 * truth
    expected_hit = np.abs(expected - truth) <= threshold + 1e-12
    replay_hit = np.abs(prediction - truth) <= threshold + 1e-12
    disagreements = int(np.sum(expected_hit != replay_hit))
    hits = int(np.sum(replay_hit))
    if max_abs_delta > 0.005:
        raise ValueError(f"Prediction delta {max_abs_delta:.9f} exceeds 0.005")
    if disagreements != 0 or hits != 134:
        raise ValueError(
            f"Replay result mismatch: hits={hits}, classification_disagreements={disagreements}"
        )

    error = prediction - truth
    training_ids = _ids(package / "external-training-ids.txt")
    training_digest = hashlib.sha256("\n".join(training_ids).encode("utf-8")).hexdigest()
    model_manifest = json.loads((package / "model-manifest.json").read_text(encoding="utf-8"))
    if training_digest != model_manifest["training"]["matchup_id_sha256"]:
        raise ValueError("Training ID order does not match model manifest")
    return {
        "status": "PASS",
        "inventory_files_verified": inventory_count,
        "model_sha256": expected_model_hash,
        "target_count": len(target_ids),
        "training_id_count": len(training_ids),
        "feature_schema_count": len(calibrator.feature_names),
        "max_abs_prediction_delta": max_abs_delta,
        "prediction_tolerance": 0.005,
        "within_ee_classification_disagreements": disagreements,
        "metrics": {
            "hits": hits,
            "within_ee_percent": 100.0 * hits / len(target_ids),
            "rmse": float(np.sqrt(np.mean(np.square(error)))),
            "mae": float(np.mean(np.abs(error))),
            "bias": float(np.mean(error)),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", type=Path, default=DEFAULT_PACKAGE)
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO)
    args = parser.parse_args()
    result = verify(args.package.resolve(), args.repo_root.resolve())
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
