"""Fit guarded, same-acquisition L2A-to-current-RT surface mappings.

The mapping is deliberately local to one historical Sentinel-2 acquisition:
clear operational L2A pixels are paired with L1C TOA corrected by the current
RT configuration for that *same* image.  It therefore contains no AERONET
values, no target-scene information, and no cross-site learned offset.

The output is retained as an operational-L2A surface prior.  The paired L1C
surface is used only to estimate a robust per-band affine frame conversion
before the L2A acquisitions are mosaicked and temporally composited.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from typing import TYPE_CHECKING, Any

import numpy as np
from tools.aeronet_validation.build_l2a_l1c_harmonization_pairs import CANONICAL_BANDS

if TYPE_CHECKING:
    from pathlib import Path

_MAD_TO_STD = 1.4826
_MIN_REFLECTANCE = 0.001
_MAX_REFLECTANCE = 0.8


def _invalid_map(reason: str, *, samples: int = 0) -> dict[str, Any]:
    return {
        "valid": False,
        "method": "identity",
        "reason": reason,
        "sample_count": int(samples),
    }


def _line_fit(source: np.ndarray, target: np.ndarray) -> tuple[float, float]:
    design = np.column_stack([np.ones(source.size, dtype=np.float64), source])
    intercept, slope = np.linalg.lstsq(design, target, rcond=None)[0]
    return float(intercept), float(slope)


def _correlation(source: np.ndarray, target: np.ndarray) -> float:
    if source.size < 2 or np.std(source) <= 1.0e-12 or np.std(target) <= 1.0e-12:
        return float("nan")
    return float(np.corrcoef(source, target)[0, 1])


def fit_band_map(
    source: np.ndarray,
    target: np.ndarray,
    *,
    min_samples: int = 200,
    min_source_iqr: float = 0.01,
    min_slope: float = 0.4,
    max_slope: float = 1.6,
    min_correlation: float = 0.7,
    max_iterations: int = 3,
) -> dict[str, Any]:
    """Fit one robust target = intercept + slope * L2A map.

    A low-spread or weakly correlated acquisition cannot constrain a slope
    safely.  In that case the map falls back to a robust scene offset instead
    of extrapolating a noisy affine fit.  The caller still clips the resulting
    pixel correction, making the map a bounded frame conversion rather than a
    free surface model.
    """
    x = np.asarray(source, dtype=np.float64).ravel()
    y = np.asarray(target, dtype=np.float64).ravel()
    usable = (
        np.isfinite(x)
        & np.isfinite(y)
        & (x > _MIN_REFLECTANCE)
        & (x < _MAX_REFLECTANCE)
        & (y > _MIN_REFLECTANCE)
        & (y < _MAX_REFLECTANCE)
    )
    x = x[usable]
    y = y[usable]
    if x.size < int(min_samples):
        return _invalid_map("insufficient_paired_pixels", samples=x.size)

    inlier = np.ones(x.size, dtype=bool)
    for _ in range(max(1, int(max_iterations))):
        intercept, slope = _line_fit(x[inlier], y[inlier])
        residual = y - (intercept + slope * x)
        center = float(np.median(residual[inlier]))
        scale = float(np.median(np.abs(residual[inlier] - center)) * _MAD_TO_STD)
        # A small absolute floor avoids over-trimming quantised dark surfaces.
        threshold = max(3.0 * scale, 0.003)
        candidate = np.abs(residual - center) <= threshold
        if int(np.count_nonzero(candidate)) < int(min_samples) or np.array_equal(candidate, inlier):
            break
        inlier = candidate

    intercept, slope = _line_fit(x[inlier], y[inlier])
    correlation = _correlation(x[inlier], y[inlier])
    source_iqr = float(np.subtract(*np.percentile(x[inlier], [75.0, 25.0])))
    method = "affine"
    if (
        source_iqr < float(min_source_iqr)
        or not math.isfinite(correlation)
        or correlation < float(min_correlation)
        or slope < float(min_slope)
        or slope > float(max_slope)
    ):
        # A robust intercept-only correction is identifiable even when a
        # scene has little spectral spread.  Keeping a unit slope avoids a
        # poorly conditioned bright/dark extrapolation.
        intercept = float(np.median(y[inlier] - x[inlier]))
        slope = 1.0
        method = "offset_fallback"

    residual = y[inlier] - (intercept + slope * x[inlier])
    residual_mad = float(np.median(np.abs(residual - np.median(residual))) * _MAD_TO_STD)
    return {
        "valid": True,
        "method": method,
        "reason": None,
        "sample_count": int(x.size),
        "inlier_count": int(np.count_nonzero(inlier)),
        "intercept": float(intercept),
        "slope": float(slope),
        "source_iqr": source_iqr,
        "correlation": correlation,
        "residual_rmse": float(np.sqrt(np.mean(np.square(residual)))),
        "residual_mad": residual_mad,
        "source_p01": float(np.percentile(x[inlier], 1.0)),
        "source_p99": float(np.percentile(x[inlier], 99.0)),
    }


def apply_band_map(
    source: np.ndarray,
    band_map: dict[str, Any],
    *,
    cap: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply a valid affine map with a bounded reflectance correction."""
    values = np.asarray(source, dtype=np.float64)
    if not bool(band_map.get("valid")):
        return values.copy(), np.zeros(values.shape, dtype=np.float64)
    intercept = float(band_map["intercept"])
    slope = float(band_map["slope"])
    mapped = intercept + slope * values
    correction = np.clip(mapped - values, -float(cap), float(cap))
    corrected = np.clip(values + correction, _MIN_REFLECTANCE, _MAX_REFLECTANCE)
    corrected = np.where(np.isfinite(values), corrected, values)
    correction = np.where(np.isfinite(values), correction, np.nan)
    return corrected, correction


def scene_key(scene: dict[str, Any]) -> str:
    """Return a stable same-acquisition identifier shared by pairs and history."""
    value = str(scene.get("scene_id") or "")
    if value:
        return value
    l2a = scene.get("l2a") or {}
    value = str(l2a.get("system_index") or "")
    if value:
        return value
    return ":".join(str(scene.get(name) or "") for name in ("day", "tile", "window"))


def build_scene_maps(
    archive_path: Path,
    *,
    cutoff: str,
    min_samples: int = 200,
    min_source_iqr: float = 0.01,
    maiac_aot_max: float | None = None,
) -> dict[str, dict[str, Any]]:
    """Load a pair archive and fit one map per eligible historical acquisition."""
    audit_path = archive_path.with_suffix(".json")
    if not archive_path.exists() or not audit_path.exists():
        raise FileNotFoundError(f"missing exact-pair archive or audit for {archive_path.stem}")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("status") != "OK" or audit.get("uses_aeronet") is not False:
        raise RuntimeError(f"pair archive is not a clean terminal result: {audit_path}")
    with np.load(archive_path, allow_pickle=False) as archive:
        source = np.asarray(archive["l2a"], dtype=np.float64)
        target = np.asarray(archive["siac"], dtype=np.float64)
        local_scene = np.asarray(archive["scene_index"], dtype=np.int64)
        names = tuple(str(value) for value in archive["band_names"].tolist())
        scenes = json.loads(str(archive["scenes_json"].item()))
    if source.shape != target.shape or source.ndim != 2 or names != CANONICAL_BANDS:
        raise ValueError(f"invalid pair archive schema in {archive_path}")
    if (
        local_scene.shape != (source.shape[0],)
        or np.any(local_scene < 0)
        or np.any(local_scene >= len(scenes))
    ):
        raise ValueError(f"invalid scene indices in {archive_path}")

    grouped: dict[str, list[np.ndarray]] = defaultdict(list)
    grouped_target: dict[str, list[np.ndarray]] = defaultdict(list)
    grouped_scene: dict[str, dict[str, Any]] = {}
    for index, scene in enumerate(scenes):
        if str(scene.get("day", "")) > cutoff:
            continue
        scene_aot = scene.get("maiac_aot", scene.get("maiac"))
        try:
            scene_aot_value = float(scene_aot)
        except (TypeError, ValueError):
            scene_aot_value = float("nan")
        if maiac_aot_max is not None and (
            not np.isfinite(scene_aot_value) or scene_aot_value > maiac_aot_max
        ):
            continue
        selected = local_scene == index
        if not np.any(selected):
            continue
        key = scene_key(scene)
        grouped[key].append(source[selected])
        grouped_target[key].append(target[selected])
        grouped_scene.setdefault(key, scene)

    output: dict[str, dict[str, Any]] = {}
    for key in sorted(grouped):
        x = np.concatenate(grouped[key], axis=0)
        y = np.concatenate(grouped_target[key], axis=0)
        output[key] = {
            "scene_id": key,
            "day": str(grouped_scene[key].get("day", "")),
            "window": str(grouped_scene[key].get("window", "")),
            "paired_sample_count": int(x.shape[0]),
            "bands": {
                band: fit_band_map(
                    x[:, band_index],
                    y[:, band_index],
                    min_samples=min_samples,
                    min_source_iqr=min_source_iqr,
                )
                for band_index, band in enumerate(CANONICAL_BANDS)
            },
        }
    if not output:
        raise ValueError(f"no paired acquisitions on or before {cutoff} in {archive_path}")
    return output
