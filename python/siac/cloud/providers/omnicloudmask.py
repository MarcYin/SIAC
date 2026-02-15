"""Default cloud/cloud-shadow provider using OmniCloudMask when available."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any

import numpy as np
import xarray as xr

from siac.cloud.mapping import apply_class_mapping


class OmniCloudMaskProvider:
    """Run OmniCloudMask (or fallback heuristic) and return standardized classes."""

    def __init__(self, predictor: Callable[[np.ndarray], np.ndarray] | None = None):
        self._predictor = predictor

    def _default_predictor(self) -> Callable[[np.ndarray], np.ndarray] | None:
        try:
            import omnicloudmask  # type: ignore[import-not-found]
        except Exception:
            return None

        # Adapter for likely OmniCloudMask APIs.
        if hasattr(omnicloudmask, "OmniCloudMask"):
            model = omnicloudmask.OmniCloudMask()
            if hasattr(model, "predict"):
                return model.predict
            if hasattr(model, "__call__"):
                return model

        if hasattr(omnicloudmask, "predict"):
            return omnicloudmask.predict

        return None

    @staticmethod
    def _heuristic(red: np.ndarray, green: np.ndarray, nir: np.ndarray) -> np.ndarray:
        missing = ~np.isfinite(red) | ~np.isfinite(green) | ~np.isfinite(nir)

        classes = np.full(red.shape, 1, dtype=np.uint8)  # clear
        cloud = ((red > 0.28) & (green > 0.30)) | (nir > 0.45)
        shadow = (~cloud) & (nir < 0.12) & (red < 0.14) & (green < 0.14)

        classes[cloud] = 2
        classes[shadow] = 3
        classes[missing] = 0
        return classes

    @staticmethod
    def _normalize_raw_output(raw: np.ndarray, template: xr.DataArray) -> xr.DataArray:
        if raw.ndim == 3 and raw.shape[-1] > 1:
            # probabilities/logits in last dimension
            raw = np.argmax(raw, axis=-1)

        if raw.shape != template.shape:
            raise ValueError(
                "OmniCloudMask output shape does not match input shape: "
                f"{raw.shape} vs {template.shape}"
            )

        arr = xr.DataArray(raw, dims=template.dims, coords=template.coords)

        # Common binary output convention: 0 clear, 1 cloud.
        values = np.unique(np.asarray(arr.values[np.isfinite(arr.values)]))
        if set(int(v) for v in values).issubset({0, 1}):
            mapped = xr.where(arr.astype(np.int16) == 1, 2, 1).astype(np.uint8)
            mapped.name = "cloud_classes"
            return mapped

        return arr

    def predict(
        self,
        red: xr.DataArray,
        green: xr.DataArray,
        nir: xr.DataArray,
        *,
        class_mapping: dict[int, Iterable[int]] | None = None,
        unmapped_to_missing: bool = True,
    ) -> xr.DataArray:
        """Predict cloud classes and map to SIAC standard classes."""
        if red.shape != green.shape or red.shape != nir.shape:
            raise ValueError("red, green, and nir arrays must have identical shape")

        predictor = self._predictor or self._default_predictor()
        if predictor is None:
            classes = self._heuristic(
                np.asarray(red.values, dtype=np.float32),
                np.asarray(green.values, dtype=np.float32),
                np.asarray(nir.values, dtype=np.float32),
            )
            out = xr.DataArray(classes, dims=red.dims, coords=red.coords, name="cloud_classes")
            # Heuristic already matches the standardized classes.
            return apply_class_mapping(out, None, unmapped_to_missing=True)

        rgbnir = np.stack(
            [
                np.asarray(red.values, dtype=np.float32),
                np.asarray(green.values, dtype=np.float32),
                np.asarray(nir.values, dtype=np.float32),
            ],
            axis=-1,
        )
        raw = predictor(rgbnir)
        normalized = self._normalize_raw_output(np.asarray(raw), red)

        # If the caller provides mapping, apply it. Otherwise keep standardized labels.
        if class_mapping:
            return apply_class_mapping(
                normalized,
                class_mapping,
                unmapped_to_missing=unmapped_to_missing,
            )

        return apply_class_mapping(
            normalized,
            None,
            unmapped_to_missing=unmapped_to_missing,
        )
