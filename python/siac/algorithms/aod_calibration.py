"""Operational inference for a fixed supervised AOD calibration model."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence


@dataclass
class GenericAODCalibrator:
    """Apply one fixed log-ratio model uniformly to SIAC retrievals.

    The estimator predicts ``log((truth + shift) / (retrieved + shift))`` from
    operational features. ``global_log_offset`` is one scalar correction shared
    by every retrieval. The feature maps must not contain evaluation labels.
    """

    feature_names: tuple[str, ...]
    estimator: Any
    global_log_offset: float = 0.0
    aod_shift: float = 1.0 / 3.0
    aod_min: float = 0.0
    aod_max: float = 4.0
    metadata: dict[str, Any] = field(default_factory=dict)

    def feature_matrix(self, feature_maps: Sequence[Mapping[str, float]]) -> np.ndarray:
        """Create the estimator matrix in the frozen training-schema order."""
        return np.asarray(
            [
                [features.get(name, np.nan) for name in self.feature_names]
                for features in feature_maps
            ],
            dtype=float,
        )

    def predict(
        self,
        retrieved_aod: Sequence[float] | np.ndarray,
        feature_maps: Sequence[Mapping[str, float]],
    ) -> np.ndarray:
        """Return uniformly calibrated AOD without using site or target labels."""
        baseline = np.asarray(retrieved_aod, dtype=float)
        if baseline.ndim != 1:
            raise ValueError("retrieved_aod must be one-dimensional")
        if len(feature_maps) != baseline.size:
            raise ValueError("Feature-map and retrieved-AOD counts do not agree")
        if not np.all(np.isfinite(baseline)):
            raise ValueError("retrieved_aod contains non-finite values")
        log_ratio = np.asarray(
            self.estimator.predict(self.feature_matrix(feature_maps)), dtype=float
        )
        if log_ratio.shape != baseline.shape or not np.all(np.isfinite(log_ratio)):
            raise ValueError("AOD calibrator produced invalid log-ratio predictions")
        model_prediction = np.clip(
            (baseline + self.aod_shift) * np.exp(log_ratio) - self.aod_shift,
            self.aod_min,
            self.aod_max,
        )
        prediction = (model_prediction + self.aod_shift) * np.exp(
            self.global_log_offset
        ) - self.aod_shift
        return np.clip(prediction, self.aod_min, self.aod_max)
