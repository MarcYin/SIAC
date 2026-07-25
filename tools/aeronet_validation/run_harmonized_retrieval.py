"""Run one Phase-D matchup with opt-in harmonized-history overrides.

The shared Phase-D harness predates the B03 seasonal experiment and its
seasonal wrapper does not forward the tau-predictor payload.  This entry point
keeps those campaign-only adaptations explicit without changing the shared
harness or the production defaults.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
HARNESS = Path("/gws/ssde/j25a/nceo_isp/public/siac_refactor")


def _band_names(value: str | None) -> tuple[str, ...]:
    bands = tuple(dict.fromkeys(part.strip() for part in (value or "").split(",") if part.strip()))
    if not bands:
        return ("B01", "B02", "B03", "B04")
    unknown = set(bands) - {"B01", "B02", "B03", "B04", "B8A", "B11", "B12"}
    if unknown:
        raise ValueError(f"unsupported direct-history bands: {sorted(unknown)}")
    return bands


def _install_overrides() -> None:
    python_dir = str(REPO / "python")
    harness_dir = str(HARNESS)
    if python_dir not in sys.path:
        sys.path.insert(0, python_dir)
    if harness_dir not in sys.path:
        sys.path.insert(0, harness_dir)

    import phaseD_prebuilt_prior

    original_prebuilt = phaseD_prebuilt_prior.make_prebuilt_surface_prior_fn
    direct_bands = _band_names(os.environ.get("PHASE_D_PREBUILT_BANDS"))

    def _prebuilt_with_requested_bands(
        matchup_id: str,
        *,
        target_band_names: tuple[str, ...] | None = None,
        **kwargs: Any,
    ) -> Any:
        return original_prebuilt(
            matchup_id,
            target_band_names=target_band_names or direct_bands,
            **kwargs,
        )

    phaseD_prebuilt_prior.make_prebuilt_surface_prior_fn = _prebuilt_with_requested_bands

    if os.environ.get("PHASE_D_SEASONAL_TAU_PREDICTOR", "0") == "1":
        from siac.algorithms.surface import seasonal_predictor

        original_predictor = seasonal_predictor.seasonal_extra_tree_prior

        def _predictor_with_tau(*args: Any, **kwargs: Any) -> Any:
            kwargs["attach_tau_predictor"] = True
            return original_predictor(*args, **kwargs)

        seasonal_predictor.seasonal_extra_tree_prior = _predictor_with_tau


def main() -> None:
    _install_overrides()
    import phaseD_array_runner

    phaseD_array_runner.main()


if __name__ == "__main__":
    main()
