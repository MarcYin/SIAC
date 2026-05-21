"""Default cloud/cloud-shadow provider using required OmniCloudMask inference."""

from __future__ import annotations

import logging
from functools import partial
from importlib import import_module
from typing import TYPE_CHECKING, cast

import numpy as np
import xarray as xr

from siac.algorithms.cloud.mapping import apply_class_mapping

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

logger = logging.getLogger(__name__)


#: Default OmniCloudMask inference device.
#:
#: PyTorch's MPS (Apple GPU) and CUDA backends are bit-deterministic
#: *within* a Python process but **not across processes** — Apple's
#: Metal driver and NVIDIA's cuDNN can pick different shader kernels or
#: algorithms on different process launches, producing ULP-level softmax
#: drift that flips argmax at edge-of-cloud pixels. Wave 18g of
#: ``REVIEW_FIXES.md`` traces how 79 such pixel flips (0.0024 % on a
#: typical S2 scene) propagate through the M5 solver's Voronoi-fill
#: amplifier to produce 100 % of AOT pixels drifting with a ~4 %
#: multiplicative bias. CPU inference is bit-deterministic across
#: processes, so we default to it.
#:
#: The performance cost on T33KWP-sized scenes is ~5-10 s per scene
#: (CPU ~25 s vs MPS ~15 s for the OCM inference); negligible against
#: the ~8 min total pipeline wall-clock. Override via the
#: ``OmniCloudMaskProvider(inference_device=...)`` constructor or the
#: ``cloud_mask.inference_device`` config field if cross-process
#: determinism isn't required.
DEFAULT_INFERENCE_DEVICE: str = "cpu"


class OmniCloudMaskProvider:
    """Run OmniCloudMask and return SIAC-standardized cloud classes."""

    def __init__(
        self,
        predictor: Callable[[np.ndarray], np.ndarray] | None = None,
        *,
        inference_device: str | None = None,
    ) -> None:
        """Construct the provider.

        Args:
            predictor: Optional pre-built predictor callable. When provided,
                bypasses the lazy ``omnicloudmask.predict_from_array`` resolve
                and the ``inference_device`` argument is ignored. Used by
                tests to inject deterministic stubs.
            inference_device: Torch device name passed to
                ``omnicloudmask.predict_from_array`` (e.g. ``"cpu"``,
                ``"cuda"``, ``"mps"``). ``None`` selects
                :data:`DEFAULT_INFERENCE_DEVICE` (CPU). Pass the string
                ``"auto"`` to let OmniCloudMask pick its own default
                (``default_device()`` — typically the fastest available
                GPU). See the module docstring on
                :data:`DEFAULT_INFERENCE_DEVICE` for why CPU is the
                conservative default.
        """
        self._predictor = predictor
        if inference_device is None:
            self._inference_device: str | None = DEFAULT_INFERENCE_DEVICE
        elif str(inference_device).lower() == "auto":
            self._inference_device = None  # let OCM pick default_device()
        else:
            self._inference_device = str(inference_device)

    @staticmethod
    def _default_predictor(
        inference_device: str | None,
    ) -> Callable[[np.ndarray], np.ndarray]:
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

        if inference_device is None:
            # Let OmniCloudMask pick — non-deterministic on GPU backends.
            logger.info(
                "OmniCloudMask: using library default device (non-deterministic "
                "across processes on GPU backends). Set "
                "cloud_mask.inference_device='cpu' for bit-deterministic results."
            )
            return cast(
                "Callable[[np.ndarray], np.ndarray]", predict_from_array
            )
        logger.info(
            "OmniCloudMask: pinning inference_device=%r for cross-process "
            "determinism.",
            inference_device,
        )
        return cast(
            "Callable[[np.ndarray], np.ndarray]",
            partial(predict_from_array, inference_device=inference_device),
        )

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

        predictor = self._predictor or self._default_predictor(self._inference_device)
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
