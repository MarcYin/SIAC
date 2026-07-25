"""Evaluate a uniform, CAMS-aware residual calibration for SIAC AOD.

The calibration consumes one SIAC retrieval, its solver diagnostics, the one
surface prior used by that retrieval, and CAMS aerosol context.  AERONET values
are used only to form the supervised target and evaluation metrics; they are
never exposed to the feature extractor.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

DEFAULT_ROOT = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")
SCHEMA_VERSION = "siac-cams-residual-v1"
TARGET_HITS = 133
N_SELECT = 25
GUARD_AOD = 0.75
UNIFORM_SVR_CONFIGS = (
    (20, 0.6, 0.03, 0.002),
    (35, 0.6, 0.01, 0.01),
    (35, 0.6, 0.01, 0.002),
)
UNIFORM_SHRINKAGE = 0.75
UNIFORM_GUARD_AOD = 0.8
HIGH_CAMS_AOD = 0.5
HIGH_BAND_SPREAD = 0.4
HIGH_CAMS_GAIN = 1.8

# These fields are evaluation products.  Keeping the denylist close to the
# extractor makes accidental label leakage fail loudly during future changes.
FORBIDDEN_FEATURE_TOKENS = (
    "aeronet",
    "truth",
    "within_ee",
    "ee_threshold",
    "error",
    "err",
    "regime",
    "runtime",
    "total_s",
)

TOP_LEVEL_FEATURES = (
    "retrieved",
    "retrieved_winmed",
    "scene_mean",
    "cloud_frac",
    "invalid_frac",
    "valid_aot_count",
    "valid_aot_fraction",
    "lookback_years",
    "aod_quality_score",
)


@dataclass(frozen=True)
class CalibrationSample:
    """One retrieval with operational features and a held-back target."""

    matchup_id: str
    site: str
    retrieved: float
    truth: float
    features: Mapping[str, float]


@dataclass
class ResidualEnsemble:
    """The fixed two-SVR residual ensemble used by the candidate algorithm."""

    feature_names: tuple[str, ...]
    model_a: Pipeline
    model_b: Pipeline
    model_b_offset: float

    def predict(self, samples: Sequence[CalibrationSample]) -> np.ndarray:
        matrix = feature_matrix(samples, self.feature_names)
        residual_a = np.asarray(self.model_a.predict(matrix), dtype=float)
        residual_b = np.asarray(self.model_b.predict(matrix), dtype=float) + self.model_b_offset
        residual = 0.25 * residual_a + 0.75 * residual_b
        baseline = np.asarray([sample.retrieved for sample in samples], dtype=float)
        return np.clip(baseline + residual, 0.0, 4.0)


@dataclass
class UniformResidualCalibrator:
    """One fixed residual ensemble plus a physical high-AOD safeguard."""

    feature_names: tuple[str, ...]
    models: tuple[Pipeline, ...]

    def predict(self, samples: Sequence[CalibrationSample]) -> tuple[np.ndarray, np.ndarray]:
        matrix = feature_matrix(samples, self.feature_names)
        baseline = np.asarray([sample.retrieved for sample in samples], dtype=float)
        candidates = np.stack(
            [baseline + np.asarray(model.predict(matrix), dtype=float) for model in self.models]
        )
        prediction = baseline + UNIFORM_SHRINKAGE * (np.mean(candidates, axis=0) - baseline)
        prediction = np.where(prediction > UNIFORM_GUARD_AOD, baseline, prediction)
        cams = np.asarray(
            [
                sample.features.get(
                    "context_cams_total_aerosol_optical_depth_at_550nm_surface",
                    np.nan,
                )
                for sample in samples
            ],
            dtype=float,
        )
        spread = np.asarray(
            [sample.features.get("solver_surface_band_argmin_spread", np.nan) for sample in samples],
            dtype=float,
        )
        high_disagreement = (cams >= HIGH_CAMS_AOD) & (spread >= HIGH_BAND_SPREAD)
        prediction = np.where(
            high_disagreement,
            np.maximum(prediction, HIGH_CAMS_GAIN * cams),
            prediction,
        )
        return np.clip(prediction, 0.0, 4.0), high_disagreement


def _finite(value: Any) -> float | None:
    if isinstance(value, bool):
        return float(value)
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _add(features: dict[str, float], name: str, value: Any) -> None:
    if any(token in name.lower() for token in FORBIDDEN_FEATURE_TOKENS):
        raise ValueError(f"Forbidden calibration feature name: {name}")
    finite = _finite(value)
    if finite is not None:
        features[name] = finite


def _ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator is None or abs(denominator) < 1e-8:
        return None
    return numerator / denominator


def site_fold(site: str, n_folds: int = 5) -> int:
    """Return the deterministic site-group fold used by the benchmark."""
    digest = hashlib.sha256(site.encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % n_folds


def extract_operational_features(
    result: Mapping[str, Any],
    aerosol_context: Mapping[str, Any],
    matchup: Mapping[str, Any],
    *,
    spectral_context: Mapping[str, Any] | None = None,
    atmo_context: Mapping[str, Any] | None = None,
    include_geography: bool = True,
) -> dict[str, float]:
    """Extract features available at retrieval time, without target fields."""
    features: dict[str, float] = {}

    for key in TOP_LEVEL_FEATURES:
        _add(features, f"siac_{key}", result.get(key))

    for group_name in ("prior_boa", "prior_unc"):
        group = result.get(group_name)
        if isinstance(group, Mapping):
            for band in ("B02", "B03", "B04"):
                _add(features, f"surface_{group_name}_{band}", group.get(band))

    anchor = result.get("anchor_iterate")
    if isinstance(anchor, Mapping):
        pass1 = _finite(anchor.get("pass1_scene_mean"))
        pass2 = _finite(anchor.get("pass2_scene_mean"))
        _add(features, "anchor_pass1_scene_mean", pass1)
        _add(features, "anchor_pass2_scene_mean", pass2)
        if pass1 is not None and pass2 is not None:
            _add(features, "anchor_pass_delta", pass2 - pass1)
            _add(features, "anchor_pass_abs_delta", abs(pass2 - pass1))

    solver = result.get("solver")
    if isinstance(solver, Mapping):
        for key, value in sorted(solver.items()):
            if isinstance(value, (bool, int, float)):
                _add(features, f"solver_{key}", value)

    atmo_aot = None
    atmo_prior = result.get("atmo_prior")
    if isinstance(atmo_prior, Mapping):
        for key, value in sorted(atmo_prior.items()):
            _add(features, f"atmo_{key}", value)
        atmo_aot = _finite(atmo_prior.get("aot_mean") or atmo_prior.get("aot_median"))
    if atmo_context is not None:
        for key in ("aot", "aot_unc", "tcwv", "tcwv_unc", "tco3", "tco3_unc", "elevation"):
            value = atmo_context.get(key)
            if value is not None:
                _add(features, f"atmo_{key}_mean", value)
                _add(features, f"atmo_{key}_median", value)
        if atmo_aot is None:
            atmo_aot = _finite(atmo_context.get("aot"))

    values = aerosol_context.get("values")
    if isinstance(values, Mapping):
        for key, value in sorted(values.items()):
            # MERRA is deliberately excluded: this candidate has one external
            # atmospheric context source and never routes between products.
            if key.startswith("cams_"):
                _add(features, f"context_{key}", value)

    sensing_time = str(
        matchup.get("sensing_time_utc")
        or aerosol_context.get("sensing_time_utc")
        or ""
    )
    try:
        when = datetime.fromisoformat(sensing_time.replace("Z", "+00:00"))
    except ValueError:
        when = None
    if when is not None:
        year_phase = 2.0 * math.pi * (when.timetuple().tm_yday - 1) / 365.25
        hour = when.hour + when.minute / 60.0 + when.second / 3600.0
        hour_phase = 2.0 * math.pi * hour / 24.0
        _add(features, "time_year_sin", math.sin(year_phase))
        _add(features, "time_year_cos", math.cos(year_phase))
        _add(features, "time_utc_sin", math.sin(hour_phase))
        _add(features, "time_utc_cos", math.cos(hour_phase))

    satellite = str(matchup.get("satellite") or "")
    _add(features, "platform_is_s2a", satellite == "S2A")
    _add(features, "platform_is_s2b", satellite == "S2B")
    _add(features, "metadata_scene_cloud_cover", matchup.get("scene_cloud_cover"))
    _add(features, "metadata_elevation_m", matchup.get("elevation_m"))

    if spectral_context is not None:
        for angle in ("sza", "vza", "raa"):
            _add(features, f"geometry_{angle}", spectral_context.get(angle))
        toa = spectral_context.get("toa")
        if isinstance(toa, Mapping):
            for band in ("B02", "B03", "B04", "B8A", "B11", "B12"):
                _add(features, f"toa_{band}", toa.get(band))
            for left, right in (
                ("B02", "B03"),
                ("B03", "B04"),
                ("B02", "B04"),
                ("B8A", "B11"),
                ("B11", "B12"),
            ):
                left_value = _finite(toa.get(left))
                right_value = _finite(toa.get(right))
                _add(
                    features,
                    f"toa_delta_{left}_{right}",
                    left_value - right_value
                    if left_value is not None and right_value is not None
                    else None,
                )
                _add(features, f"toa_ratio_{left}_{right}", _ratio(left_value, right_value))
            red = _finite(toa.get("B04"))
            nir = _finite(toa.get("B8A"))
            swir1 = _finite(toa.get("B11"))
            swir2 = _finite(toa.get("B12"))
            _add(
                features,
                "toa_ndvi_B8A_B04",
                _ratio(nir - red, nir + red)
                if nir is not None and red is not None
                else None,
            )
            _add(
                features,
                "toa_ndmi_B8A_B11",
                _ratio(nir - swir1, nir + swir1)
                if nir is not None and swir1 is not None
                else None,
            )
            _add(
                features,
                "toa_ndsi_B11_B12",
                _ratio(swir1 - swir2, swir1 + swir2)
                if swir1 is not None and swir2 is not None
                else None,
            )

    lat = _finite(result.get("lat") if result.get("lat") is not None else matchup.get("latitude"))
    lon = _finite(result.get("lon") if result.get("lon") is not None else matchup.get("longitude"))
    if include_geography:
        _add(features, "geo_latitude", lat)
        _add(features, "geo_longitude", lon)
        _add(features, "geo_abs_latitude", abs(lat) if lat is not None else None)
        if lat is not None:
            _add(features, "geo_latitude_sin", math.sin(math.radians(lat)))
            _add(features, "geo_latitude_cos", math.cos(math.radians(lat)))
        if lon is not None:
            _add(features, "geo_longitude_sin", math.sin(math.radians(lon)))
            _add(features, "geo_longitude_cos", math.cos(math.radians(lon)))

    # Compact physical contrasts make the fixed learner less dependent on the
    # raw feature scales and expose surface/atmosphere consistency directly.
    retrieved = _finite(result.get("retrieved"))
    prior_boa = result.get("prior_boa") if isinstance(result.get("prior_boa"), Mapping) else {}
    prior_unc = result.get("prior_unc") if isinstance(result.get("prior_unc"), Mapping) else {}
    for left, right in (("B02", "B03"), ("B03", "B04"), ("B02", "B04")):
        left_boa = _finite(prior_boa.get(left))
        right_boa = _finite(prior_boa.get(right))
        left_unc = _finite(prior_unc.get(left))
        right_unc = _finite(prior_unc.get(right))
        _add(
            features,
            f"surface_boa_delta_{left}_{right}",
            left_boa - right_boa if left_boa is not None and right_boa is not None else None,
        )
        _add(features, f"surface_boa_ratio_{left}_{right}", _ratio(left_boa, right_boa))
        _add(features, f"surface_unc_ratio_{left}_{right}", _ratio(left_unc, right_unc))

    cams_total = None
    cams_469 = None
    cams_670 = None
    cams_865 = None
    if isinstance(values, Mapping):
        cams_total = _finite(values.get("cams_total_aerosol_optical_depth_at_550nm_surface"))
        cams_469 = _finite(values.get("cams_total_aerosol_optical_depth_at_469nm_surface"))
        cams_670 = _finite(values.get("cams_total_aerosol_optical_depth_at_670nm_surface"))
        cams_865 = _finite(values.get("cams_total_aerosol_optical_depth_at_865nm_surface"))
    _add(features, "consistency_cams_minus_siac", cams_total - retrieved if cams_total is not None and retrieved is not None else None)
    _add(features, "consistency_cams_to_siac", _ratio(cams_total, retrieved))
    _add(
        features,
        "consistency_atmo_minus_siac",
        atmo_aot - retrieved if atmo_aot is not None and retrieved is not None else None,
    )
    _add(features, "consistency_atmo_to_siac", _ratio(atmo_aot, retrieved))
    _add(
        features,
        "consistency_cams_minus_atmo",
        cams_total - atmo_aot if cams_total is not None and atmo_aot is not None else None,
    )
    _add(features, "consistency_cams_to_atmo", _ratio(cams_total, atmo_aot))
    _add(features, "context_cams_ratio_469_550", _ratio(cams_469, cams_total))
    _add(features, "context_cams_ratio_670_550", _ratio(cams_670, cams_total))
    _add(features, "context_cams_ratio_865_550", _ratio(cams_865, cams_total))
    if cams_469 is not None and cams_865 is not None and cams_469 > 0 and cams_865 > 0:
        angstrom = -math.log(cams_469 / cams_865) / math.log(469.0 / 865.0)
        _add(features, "context_cams_angstrom_469_865", angstrom)

    return features


def load_matchups(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return {row["matchup_id"]: row for row in csv.DictReader(stream)}


def load_samples(
    result_dir: Path,
    context_dir: Path,
    matchups_path: Path,
    matchup_ids: Iterable[str],
    *,
    spectral_dir: Path | None = None,
    atmo_context_dir: Path | None = None,
    include_geography: bool = True,
    require_complete: bool = True,
) -> list[CalibrationSample]:
    matchups = load_matchups(matchups_path)
    samples: list[CalibrationSample] = []
    missing: list[str] = []
    failed: list[str] = []
    for matchup_id in matchup_ids:
        result_path = result_dir / f"{matchup_id}.json"
        context_path = context_dir / f"{matchup_id}.json"
        spectral_path = spectral_dir / f"{matchup_id}.json" if spectral_dir is not None else None
        atmo_context_path = (
            atmo_context_dir / f"{matchup_id}.json" if atmo_context_dir is not None else None
        )
        if (
            not result_path.exists()
            or not context_path.exists()
            or matchup_id not in matchups
            or (spectral_path is not None and not spectral_path.exists())
        ):
            missing.append(matchup_id)
            continue
        result = json.loads(result_path.read_text(encoding="utf-8"))
        context = json.loads(context_path.read_text(encoding="utf-8"))
        spectral = (
            json.loads(spectral_path.read_text(encoding="utf-8"))
            if spectral_path is not None
            else None
        )
        atmo_context = (
            json.loads(atmo_context_path.read_text(encoding="utf-8"))
            if atmo_context_path is not None and atmo_context_path.exists()
            else None
        )
        if result.get("status") != "OK" or context.get("status") != "OK":
            failed.append(matchup_id)
            continue
        retrieved = _finite(result.get("retrieved"))
        truth = _finite(result.get("truth"))
        if retrieved is None or truth is None:
            failed.append(matchup_id)
            continue
        samples.append(
            CalibrationSample(
                matchup_id=matchup_id,
                site=str(result.get("site") or matchups[matchup_id].get("site") or ""),
                retrieved=retrieved,
                truth=truth,
                features=extract_operational_features(
                    result,
                    context,
                    matchups[matchup_id],
                    spectral_context=spectral,
                    atmo_context=atmo_context,
                    include_geography=include_geography,
                ),
            )
        )
    if require_complete and (missing or failed):
        raise ValueError(f"Incomplete samples: missing={len(missing)} failed={len(failed)}")
    return samples


def feature_schema(samples: Sequence[CalibrationSample]) -> tuple[str, ...]:
    return tuple(sorted({name for sample in samples for name in sample.features}))


def feature_matrix(samples: Sequence[CalibrationSample], names: Sequence[str]) -> np.ndarray:
    return np.asarray(
        [[sample.features.get(name, np.nan) for name in names] for sample in samples],
        dtype=float,
    )


def _pipeline(*, c: float, gamma: float, epsilon: float, n_features: int) -> Pipeline:
    return Pipeline(
        steps=(
            ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("scale", StandardScaler()),
            ("select", SelectKBest(f_regression, k=min(N_SELECT, n_features))),
            ("svr", SVR(C=c, gamma=gamma, epsilon=epsilon)),
        )
    )


def fit_ensemble(
    samples: Sequence[CalibrationSample],
    feature_names: Sequence[str] | None = None,
) -> ResidualEnsemble:
    names = tuple(feature_names or feature_schema(samples))
    matrix = feature_matrix(samples, names)
    residual = np.asarray([sample.truth - sample.retrieved for sample in samples], dtype=float)
    model_a = _pipeline(c=0.6, gamma=0.03, epsilon=0.002, n_features=len(names))
    model_b = _pipeline(c=1.0, gamma=0.001, epsilon=0.02, n_features=len(names))
    model_a.fit(matrix, residual)
    model_b.fit(matrix, residual)
    offset = float(np.median(residual - model_b.predict(matrix)))
    return ResidualEnsemble(names, model_a, model_b, offset)


def fit_uniform_calibrator(
    samples: Sequence[CalibrationSample],
    feature_names: Sequence[str] | None = None,
) -> UniformResidualCalibrator:
    """Fit the frozen uniform candidate on an independent development cohort."""
    names = tuple(feature_names or feature_schema(samples))
    matrix = feature_matrix(samples, names)
    residual = np.asarray([sample.truth - sample.retrieved for sample in samples], dtype=float)
    models: list[Pipeline] = []
    for n_select, c, gamma, epsilon in UNIFORM_SVR_CONFIGS:
        model = Pipeline(
            steps=(
                ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
                ("scale", StandardScaler()),
                ("select", SelectKBest(f_regression, k=min(n_select, len(names)))),
                ("svr", SVR(C=c, gamma=gamma, epsilon=epsilon)),
            )
        )
        model.fit(matrix, residual)
        models.append(model)
    return UniformResidualCalibrator(names, tuple(models))


def apply_guard(baseline: np.ndarray, candidate: np.ndarray, threshold: float = GUARD_AOD) -> np.ndarray:
    """Retain the uncalibrated retrieval when the candidate is high-AOD."""
    return np.where(candidate > threshold, baseline, candidate)


def expected_error(truth: np.ndarray) -> np.ndarray:
    return 0.05 + 0.15 * truth


def metrics(truth: np.ndarray, prediction: np.ndarray) -> dict[str, float | int]:
    error = prediction - truth
    within = np.abs(error) <= expected_error(truth)
    return {
        "count": int(truth.size),
        "hits": int(within.sum()),
        "within_ee_percent": float(100.0 * within.mean()),
        "rmse": float(np.sqrt(np.mean(error**2))),
        "mae": float(np.mean(np.abs(error))),
        "bias": float(np.mean(error)),
    }


def selected_features(model: Pipeline, names: Sequence[str]) -> list[str]:
    support = model.named_steps["select"].get_support()
    return [name for name, selected in zip(names, support, strict=True) if selected]


def evaluate_site_oof(samples: Sequence[CalibrationSample]) -> dict[str, Any]:
    names = feature_schema(samples)
    truth = np.asarray([sample.truth for sample in samples], dtype=float)
    baseline = np.asarray([sample.retrieved for sample in samples], dtype=float)
    candidate = np.full(len(samples), np.nan, dtype=float)
    fold_receipts: list[dict[str, Any]] = []
    folds = np.asarray([site_fold(sample.site) for sample in samples], dtype=int)
    for fold in range(5):
        train = [sample for sample, sample_fold in zip(samples, folds, strict=True) if sample_fold != fold]
        test_indices = np.flatnonzero(folds == fold)
        test = [samples[index] for index in test_indices]
        model = fit_ensemble(train, names)
        candidate[test_indices] = model.predict(test)
        fold_receipts.append(
            {
                "fold": fold,
                "train_count": len(train),
                "test_count": len(test),
                "model_a_features": selected_features(model.model_a, names),
                "model_b_features": selected_features(model.model_b, names),
                "model_b_offset": model.model_b_offset,
            }
        )
    guarded = apply_guard(baseline, candidate)
    baseline_hit = np.abs(baseline - truth) <= expected_error(truth)
    candidate_hit = np.abs(candidate - truth) <= expected_error(truth)
    guarded_hit = np.abs(guarded - truth) <= expected_error(truth)
    rows = []
    for index, sample in enumerate(samples):
        rows.append(
            {
                "matchup_id": sample.matchup_id,
                "site": sample.site,
                "fold": int(folds[index]),
                "truth": sample.truth,
                "baseline": sample.retrieved,
                "candidate": float(candidate[index]),
                "guarded": float(guarded[index]),
                "baseline_within_ee": bool(baseline_hit[index]),
                "candidate_within_ee": bool(candidate_hit[index]),
                "guarded_within_ee": bool(guarded_hit[index]),
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "feature_count": len(names),
        "feature_names": list(names),
        "guard_aod": GUARD_AOD,
        "target_hits": TARGET_HITS,
        "baseline": metrics(truth, baseline),
        "candidate": metrics(truth, candidate),
        "guarded": metrics(truth, guarded),
        "gains": int((guarded_hit & ~baseline_hit).sum()),
        "losses": int((~guarded_hit & baseline_hit).sum()),
        "folds": fold_receipts,
        "predictions": rows,
    }


def evaluate_train_test(
    train: Sequence[CalibrationSample],
    test: Sequence[CalibrationSample],
) -> dict[str, Any]:
    names = feature_schema(train)
    model = fit_uniform_calibrator(train, names)
    truth = np.asarray([sample.truth for sample in test], dtype=float)
    baseline = np.asarray([sample.retrieved for sample in test], dtype=float)
    candidate, high_disagreement = model.predict(test)
    baseline_hit = np.abs(baseline - truth) <= expected_error(truth)
    candidate_hit = np.abs(candidate - truth) <= expected_error(truth)
    return {
        "schema_version": SCHEMA_VERSION,
        "train_count": len(train),
        "test_count": len(test),
        "feature_count": len(names),
        "guard_aod": UNIFORM_GUARD_AOD,
        "model_features": [selected_features(fitted, names) for fitted in model.models],
        "uniform_recipe": {
            "svr_configs": [list(config) for config in UNIFORM_SVR_CONFIGS],
            "shrinkage": UNIFORM_SHRINKAGE,
            "guard_aod": UNIFORM_GUARD_AOD,
            "high_cams_aod": HIGH_CAMS_AOD,
            "high_band_spread": HIGH_BAND_SPREAD,
            "high_cams_gain": HIGH_CAMS_GAIN,
        },
        "baseline": metrics(truth, baseline),
        "candidate": metrics(truth, candidate),
        "gains": int((candidate_hit & ~baseline_hit).sum()),
        "losses": int((~candidate_hit & baseline_hit).sum()),
        "high_disagreement_count": int(high_disagreement.sum()),
        "predictions": [
            {
                "matchup_id": sample.matchup_id,
                "site": sample.site,
                "truth": sample.truth,
                "baseline": sample.retrieved,
                "candidate": float(candidate[index]),
                "high_disagreement": bool(high_disagreement[index]),
                "baseline_within_ee": bool(baseline_hit[index]),
                "candidate_within_ee": bool(candidate_hit[index]),
                "transition": (
                    "gain"
                    if candidate_hit[index] and not baseline_hit[index]
                    else "loss"
                    if baseline_hit[index] and not candidate_hit[index]
                    else "retained_hit"
                    if candidate_hit[index]
                    else "retained_miss"
                ),
            }
            for index, sample in enumerate(test)
        ],
    }


def _ids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--context-dir", type=Path, required=True)
    parser.add_argument("--spectral-dir", type=Path)
    parser.add_argument("--atmo-context-dir", type=Path)
    parser.add_argument("--mids", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--train-result-dir", type=Path)
    parser.add_argument("--train-context-dir", type=Path)
    parser.add_argument("--train-mids", type=Path)
    parser.add_argument("--train-atmo-context-dir", type=Path)
    parser.add_argument("--no-geography", action="store_true")
    args = parser.parse_args()
    samples = load_samples(
        args.result_dir,
        args.context_dir,
        args.root / "matchups" / "matchups.csv",
        _ids(args.mids),
        spectral_dir=args.spectral_dir,
        atmo_context_dir=args.atmo_context_dir,
        include_geography=not args.no_geography,
    )
    train_options = (args.train_result_dir, args.train_context_dir, args.train_mids)
    if any(train_options) and not all(train_options):
        parser.error("--train-result-dir, --train-context-dir, and --train-mids are required together")
    if all(train_options):
        train = load_samples(
            args.train_result_dir,
            args.train_context_dir,
            args.root / "matchups" / "matchups.csv",
            _ids(args.train_mids),
            atmo_context_dir=args.train_atmo_context_dir,
            include_geography=not args.no_geography,
            require_complete=False,
        )
        analysis = evaluate_train_test(train, samples)
    else:
        analysis = evaluate_site_oof(samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    summary_keys = ("baseline", "candidate", "guarded", "gains", "losses")
    print(json.dumps({key: analysis[key] for key in summary_keys if key in analysis}, indent=2))
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
