"""Evaluate one-prior per-band AOD consensus with site-grouped validation."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
COHORT = "campaign250_lowcloud20_mids.txt"
R2_DIR = "phaseD_results_campaign250_R2_full_localdiag_20260705"
B03_DIR = "phaseD_results_lowcloud20_singleprior_b03_chi2_20260711"
PRIOR_DIR = "prior_quality"
WIDENED_DIR = "phaseD_results_lowcloud20_singleprior_b03_trimmed_bs3_20260711"
OUTPUT = "reports/aod-low-cloud-20260711/single-prior-consensus-analysis.json"
WEIGHTS = (0.175, 0.20, 0.225, 0.25, 0.275, 0.30)


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _hit(prediction: float, truth: float) -> bool:
    return abs(float(prediction) - float(truth)) <= 0.05 + 0.15 * float(truth) + 1e-12


def _site_fold(site: str, folds: int = 5) -> int:
    digest = hashlib.sha256(site.encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "big") % folds


def _summary(prediction: np.ndarray, truth: np.ndarray) -> dict[str, float | int]:
    error = np.asarray(prediction, dtype=np.float64) - np.asarray(truth, dtype=np.float64)
    threshold = 0.05 + 0.15 * np.asarray(truth, dtype=np.float64)
    hits = int(np.count_nonzero(np.abs(error) <= threshold + 1e-12))
    return {
        "hits": hits,
        "count": int(error.size),
        "rate": float(hits / error.size) if error.size else 0.0,
        "rmse": float(np.sqrt(np.mean(np.square(error)))) if error.size else math.nan,
        "bias": float(np.mean(error)) if error.size else math.nan,
    }


def _load_rows(root: Path, cohort: list[str]) -> list[dict[str, Any]]:
    rows = []
    for matchup_id in cohort:
        r2_path = root / R2_DIR / f"{matchup_id}.json"
        b03_path = root / B03_DIR / f"{matchup_id}.json"
        prior_path = root / PRIOR_DIR / f"{matchup_id}.json"
        if not (r2_path.exists() and b03_path.exists() and prior_path.exists()):
            continue
        r2 = json.loads(r2_path.read_text(encoding="utf-8"))
        b03 = json.loads(b03_path.read_text(encoding="utf-8"))
        prior = json.loads(prior_path.read_text(encoding="utf-8"))
        solver = b03.get("solver") or {}
        band_values = [
            _finite(solver.get(f"surface_band_{band}_argmin_aot")) for band in ("B02", "B03", "B04")
        ]
        values = {
            "truth": _finite(r2.get("truth")),
            "baseline": _finite(r2.get("retrieved")),
            "b03_retrieved": _finite(b03.get("retrieved")),
            "cams": _finite(prior.get("cams_aot")),
            "cost": _finite(solver.get("cost_final_per_band")),
            "spread": _finite(solver.get("surface_band_argmin_spread")),
            "flat_delta": _finite(solver.get("surface_cost_curve_relative_second_delta")),
            "cloud_fraction": _finite(r2.get("cloud_frac")),
            "latitude": _finite(r2.get("lat")),
            "longitude": _finite(r2.get("lon")),
        }
        if any(value is None for value in band_values) or any(
            value is None for value in values.values()
        ):
            continue
        bands = np.asarray(band_values, dtype=np.float64)
        assert values["cams"] is not None
        rows.append(
            {
                "matchup_id": matchup_id,
                "site": str(r2.get("site") or matchup_id.split("__", 1)[0]),
                **values,
                "band_mean": float(np.mean(bands)),
                "band_median": float(np.median(bands)),
                "band_std": float(np.std(bands)),
            }
        )
    return rows


def _posterior(rows: list[dict[str, Any]], weight: float, aggregation: str) -> np.ndarray:
    surface = np.asarray([row[f"band_{aggregation}"] for row in rows], dtype=np.float64)
    cams = np.asarray([row["cams"] for row in rows], dtype=np.float64)
    return float(weight) * surface + (1.0 - float(weight)) * cams


def _condition(values: np.ndarray, threshold: float, direction: str) -> np.ndarray:
    return values > threshold if direction == "high" else values < threshold


def _feature_values(
    rows: list[dict[str, Any]],
    indices: np.ndarray,
    feature: str,
    posterior: np.ndarray,
) -> np.ndarray:
    baseline = np.asarray([rows[index]["baseline"] for index in indices], dtype=np.float64)
    if feature == "posterior_delta":
        return np.abs(posterior - baseline)
    if feature == "cams_delta":
        cams = np.asarray([rows[index]["cams"] for index in indices], dtype=np.float64)
        return np.abs(cams - baseline)
    if feature == "surface_delta":
        surface = np.asarray([rows[index]["band_mean"] for index in indices], dtype=np.float64)
        return np.abs(surface - baseline)
    return np.asarray([rows[index][feature] for index in indices], dtype=np.float64)


def _fit_gate(
    rows: list[dict[str, Any]],
    indices: np.ndarray,
) -> dict[str, Any]:
    truth = np.asarray([rows[index]["truth"] for index in indices], dtype=np.float64)
    baseline = np.asarray([rows[index]["baseline"] for index in indices], dtype=np.float64)
    threshold_ee = 0.05 + 0.15 * truth
    candidates = []
    for aggregation in ("mean", "median"):
        for weight in WEIGHTS:
            posterior = _posterior(rows, weight, aggregation)[indices]
            for feature, direction in (
                ("cost", "high"),
                ("spread", "high"),
                ("flat_delta", "low"),
                ("band_std", "high"),
                ("posterior_delta", "high"),
                ("cams_delta", "high"),
                ("surface_delta", "high"),
            ):
                values = _feature_values(rows, indices, feature, posterior)
                for threshold in np.unique(np.quantile(values, np.linspace(0.05, 0.95, 19))):
                    use = _condition(values, float(threshold), direction)
                    prediction = np.where(use, posterior, baseline)
                    error = prediction - truth
                    candidates.append(
                        {
                            "aggregation": aggregation,
                            "weight": float(weight),
                            "feature": feature,
                            "direction": direction,
                            "threshold": float(threshold),
                            "hits": int(np.count_nonzero(np.abs(error) <= threshold_ee + 1e-12)),
                            "rmse": float(np.sqrt(np.mean(np.square(error)))),
                            "switches": int(np.count_nonzero(use)),
                        }
                    )
    return min(
        candidates,
        key=lambda candidate: (
            -candidate["hits"],
            candidate["rmse"],
            candidate["switches"],
        ),
    )


def _apply_gate(
    rows: list[dict[str, Any]], indices: np.ndarray, gate: dict[str, Any]
) -> np.ndarray:
    baseline = np.asarray([rows[index]["baseline"] for index in indices], dtype=np.float64)
    posterior = _posterior(rows, gate["weight"], gate["aggregation"])[indices]
    values = _feature_values(rows, indices, gate["feature"], posterior)
    use = _condition(values, gate["threshold"], gate["direction"])
    return np.where(use, posterior, baseline)


def _tree_oof_gate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    from sklearn.impute import SimpleImputer
    from sklearn.pipeline import make_pipeline
    from sklearn.tree import DecisionTreeClassifier, export_text

    feature_names, x, baseline, posterior, truth, label, folds = _classification_inputs(
        rows, weight=0.25, aggregation="mean"
    )
    prediction = np.full(len(rows), np.nan, dtype=np.float64)
    rules = []
    for fold in range(5):
        test = np.flatnonzero(folds == fold)
        train = np.flatnonzero(folds != fold)
        if not test.size or not train.size or np.unique(label[train]).size < 2:
            continue
        model = make_pipeline(
            SimpleImputer(strategy="median"),
            DecisionTreeClassifier(
                max_depth=2,
                min_samples_leaf=8,
                class_weight="balanced",
                random_state=0,
            ),
        )
        model.fit(x[train], label[train])
        choose_posterior = model.predict(x[test]).astype(bool)
        prediction[test] = np.where(choose_posterior, posterior[test], baseline[test])
        rules.append(
            {
                "fold": fold,
                "test_count": int(test.size),
                "posterior_count": int(np.count_nonzero(choose_posterior)),
                "rule": export_text(
                    model.named_steps["decisiontreeclassifier"],
                    feature_names=list(feature_names),
                ),
            }
        )
    valid = np.isfinite(prediction)
    return {
        "description": "Fixed depth-2 site-grouped diagnostic gate; not a final trained model",
        "posterior": "0.25 * mean(B02,B03,B04 minima) + 0.75 * CAMS anchor",
        "summary": _summary(prediction[valid], truth[valid]) if valid.any() else None,
        "fold_rules": rules,
    }


def _classification_inputs(
    rows: list[dict[str, Any]],
    *,
    weight: float,
    aggregation: str,
) -> tuple[tuple[str, ...], np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    feature_names = (
        "baseline",
        "b03_retrieved",
        "cams",
        "band_mean",
        "band_median",
        "band_std",
        "cost",
        "spread",
        "flat_delta",
        "posterior_delta",
        "cams_delta",
        "surface_delta",
        "cloud_fraction",
        "latitude",
        "longitude",
    )
    posterior = _posterior(rows, weight, aggregation)
    baseline = np.asarray([row["baseline"] for row in rows], dtype=np.float64)
    truth = np.asarray([row["truth"] for row in rows], dtype=np.float64)
    matrix = []
    for index, row in enumerate(rows):
        matrix.append(
            [
                row["baseline"],
                row["b03_retrieved"],
                row["cams"],
                row["band_mean"],
                row["band_median"],
                row["band_std"],
                row["cost"],
                row["spread"],
                row["flat_delta"],
                abs(posterior[index] - baseline[index]),
                abs(row["cams"] - baseline[index]),
                abs(row["band_mean"] - baseline[index]),
                row["cloud_fraction"],
                row["latitude"],
                row["longitude"],
            ]
        )
    x = np.asarray(matrix, dtype=np.float64)
    label = (np.abs(posterior - truth) < np.abs(baseline - truth)).astype(np.int8)
    folds = np.asarray([_site_fold(row["site"]) for row in rows], dtype=np.int64)
    return feature_names, x, baseline, posterior, truth, label, folds


def _extra_trees_oof_gate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    from sklearn.ensemble import ExtraTreesClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.pipeline import make_pipeline

    _, x, baseline, posterior, truth, label, folds = _classification_inputs(
        rows, weight=0.275, aggregation="median"
    )
    prediction = np.full(len(rows), np.nan, dtype=np.float64)
    fold_rows = []
    for fold in range(5):
        test = np.flatnonzero(folds == fold)
        train = np.flatnonzero(folds != fold)
        if not test.size or not train.size or np.unique(label[train]).size < 2:
            continue
        model = make_pipeline(
            SimpleImputer(strategy="median"),
            ExtraTreesClassifier(
                n_estimators=300,
                max_depth=4,
                min_samples_leaf=6,
                class_weight="balanced",
                random_state=0,
                n_jobs=1,
            ),
        )
        model.fit(x[train], label[train])
        choose_posterior = model.predict(x[test]).astype(bool)
        prediction[test] = np.where(choose_posterior, posterior[test], baseline[test])
        fold_rows.append(
            {
                "fold": fold,
                "test_count": int(test.size),
                "posterior_count": int(np.count_nonzero(choose_posterior)),
            }
        )
    valid = np.isfinite(prediction)
    return {
        "description": "ExtraTrees site-grouped robustness audit selected after model comparison",
        "posterior": "0.275 * median(B02,B03,B04 minima) + 0.725 * CAMS anchor",
        "summary": _summary(prediction[valid], truth[valid]) if valid.any() else None,
        "folds": fold_rows,
    }


def analyze(root: Path) -> dict[str, Any]:
    cohort = [
        line.strip()
        for line in (root / COHORT).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    rows = _load_rows(root, cohort)
    truth = np.asarray([row["truth"] for row in rows], dtype=np.float64)
    baseline = np.asarray([row["baseline"] for row in rows], dtype=np.float64)
    b03_retrieval = np.asarray([row["b03_retrieved"] for row in rows], dtype=np.float64)
    fixed = {}
    posterior_predictions = []
    for aggregation in ("mean", "median"):
        for weight in WEIGHTS:
            prediction = _posterior(rows, weight, aggregation)
            posterior_predictions.append(prediction)
            fixed[f"{aggregation}_surface_{weight:.3f}"] = _summary(prediction, truth)

    oracle = np.any(
        np.stack(
            [np.abs(baseline - truth) <= 0.05 + 0.15 * truth]
            + [
                np.abs(prediction - truth) <= 0.05 + 0.15 * truth
                for prediction in posterior_predictions
            ]
        ),
        axis=0,
    )
    row_index = {row["matchup_id"]: index for index, row in enumerate(rows)}
    widened_hits = np.zeros(len(rows), dtype=bool)
    for path in (root / WIDENED_DIR).glob("*.json"):
        record = json.loads(path.read_text(encoding="utf-8"))
        index = row_index.get(str(record.get("matchup_id")))
        retrieved = _finite(record.get("retrieved"))
        if index is None or str(record.get("status")).upper() != "OK" or retrieved is None:
            continue
        widened_hits[index] = _hit(retrieved, float(truth[index]))
    b03_hits = np.abs(b03_retrieval - truth) <= 0.05 + 0.15 * truth
    expanded_oracle = oracle | b03_hits | widened_hits

    folds = np.asarray([_site_fold(row["site"]) for row in rows], dtype=np.int64)
    oof_prediction = np.full(len(rows), np.nan, dtype=np.float64)
    fold_gates = []
    for fold in range(5):
        test = np.flatnonzero(folds == fold)
        train = np.flatnonzero(folds != fold)
        if not test.size or not train.size:
            continue
        gate = _fit_gate(rows, train)
        test_prediction = _apply_gate(rows, test, gate)
        oof_prediction[test] = test_prediction
        test_truth = np.asarray([rows[index]["truth"] for index in test], dtype=np.float64)
        fold_gates.append(
            {
                "fold": fold,
                "test_count": int(test.size),
                "aggregation": gate["aggregation"],
                "weight": gate["weight"],
                "feature": gate["feature"],
                "direction": gate["direction"],
                "threshold": gate["threshold"],
                "train_hits": gate["hits"],
                "train_rmse": gate["rmse"],
                "train_switches": gate["switches"],
                "test": _summary(test_prediction, test_truth),
            }
        )

    complete_oof = np.isfinite(oof_prediction)
    best_fixed_name, best_fixed = min(
        fixed.items(),
        key=lambda item: (-item[1]["hits"], item[1]["rmse"], abs(item[1]["bias"])),
    )
    return {
        "cohort_count": len(cohort),
        "available": len(rows),
        "target_hits": 133,
        "surface_prior_type": "one S2 monthly SWIR/NIR-anchored ExtraTree predictor",
        "baseline": _summary(baseline, truth),
        "b03_retrieval": _summary(b03_retrieval, truth),
        "fixed_posteriors": fixed,
        "best_fixed_posterior": {"candidate": best_fixed_name, **best_fixed},
        "baseline_plus_fixed_posterior_oracle_hits": int(np.count_nonzero(oracle)),
        "expanded_single_prior_oracle_hits": int(np.count_nonzero(expanded_oracle)),
        "expanded_single_prior_oracle_rate": float(np.mean(expanded_oracle)),
        "expanded_single_prior_oracle_unresolved_ids": [
            rows[index]["matchup_id"] for index in np.flatnonzero(~expanded_oracle)
        ],
        "widened_backstop_screen_hits": int(np.count_nonzero(widened_hits)),
        "oof_gate": (
            _summary(oof_prediction[complete_oof], truth[complete_oof])
            if complete_oof.any()
            else None
        ),
        "fold_gates": fold_gates,
        "tree_oof_gate": _tree_oof_gate(rows) if len(rows) >= 20 else None,
        "extra_trees_oof_gate": _extra_trees_oof_gate(rows) if len(rows) >= 20 else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = analyze(args.root)
    output = args.output or args.root / OUTPUT
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
