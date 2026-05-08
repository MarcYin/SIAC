"""Default cloud/cloud-shadow provider using required OmniCloudMask inference."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, cast

import numpy as np
import xarray as xr

from siac.algorithms.cloud.mapping import apply_class_mapping

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable


class OmniCloudMaskProvider:
    """Run OmniCloudMask and return SIAC-standardized cloud classes."""

    def __init__(self, predictor: Callable[[np.ndarray], np.ndarray] | None = None):
        self._predictor = predictor

    @staticmethod
    def _default_predictor() -> Callable[[np.ndarray], np.ndarray]:
        try:
            omnicloudmask = import_module("omnicloudmask")
        except ImportError as exc:
            # Narrowed from ``except Exception`` (REVIEW.md §2.1):
            # only swallow import-time failures here. Other exception
            # classes (a partial install raising AttributeError during
            # module init, for example) should propagate so the user can
            # see the real cause.
            raise ImportError(
                "omnicloudmask is required for SIAC cloud masking. "
                "Install the 'omnicloudmask' package."
            ) from exc

        predict_from_array = getattr(omnicloudmask, "predict_from_array", None)
        if not callable(predict_from_array):
            raise RuntimeError(
                "omnicloudmask does not expose a callable predict_from_array() entrypoint"
            )

        return cast("Callable[[np.ndarray], np.ndarray]", predict_from_array)

    @staticmethod
    def _normalize_raw_output(raw: np.ndarray, template: xr.DataArray) -> xr.DataArray:
        if raw.ndim == 3:
            if raw.shape[0] == 1 and raw.shape[1:] == template.shape:
                raw = raw[0]
            elif raw.shape[0] > 1 and raw.shape[1:] == template.shape:
                # Channel-first probabilities/logits from OCM confidence output.
                raw = np.argmax(raw, axis=0)
            elif raw.shape[-1] > 1 and raw.shape[:-1] == template.shape:
                raw = np.argmax(raw, axis=-1)

        if raw.shape != template.shape:
            raise ValueError(
                "OmniCloudMask output shape does not match input shape: "
                f"{raw.shape} vs {template.shape}"
            )

        arr = xr.DataArray(raw, dims=template.dims, coords=template.coords)

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
        red_np = np.asarray(red.values, dtype=np.float32)
        green_np = np.asarray(green.values, dtype=np.float32)
        nir_np = np.asarray(nir.values, dtype=np.float32)
        missing = ~np.isfinite(red_np) | ~np.isfinite(green_np) | ~np.isfinite(nir_np)
        rgbnir = np.stack(
            [
                np.where(np.isfinite(red_np), red_np, 0.0),
                np.where(np.isfinite(green_np), green_np, 0.0),
                np.where(np.isfinite(nir_np), nir_np, 0.0),
            ],
            axis=0,
        )
        raw = predictor(rgbnir)
        normalized = self._normalize_raw_output(np.asarray(raw), red)

        if class_mapping:
            out = apply_class_mapping(
                normalized,
                class_mapping,
                unmapped_to_missing=unmapped_to_missing,
            )
        else:
            # OmniCloudMask classes: 0 clear, 1 thick cloud, 2 thin cloud, 3 shadow.
            out = apply_class_mapping(
                normalized,
                {1: [0], 2: [1, 2], 3: [3]},
                unmapped_to_missing=unmapped_to_missing,
            )

        out = out.where(~missing, other=np.uint8(0)).astype(np.uint8)
        out.name = "cloud_classes"
        return out
