"""Prepared monthly-composite build workflow."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from siac.adapters.auth import CredentialManager
from siac.algorithms.surface.brdf_monthly_composite import MonthlyKernelWeightComposite
from siac.algorithms.surface.monthly_composite_store import (
    filter_materialized_monthly_composite_collection,
    write_monthly_composite_collection,
)
from siac.algorithms.surface.swir_refine import build_monthly_composites_from_brdf
from siac.app.assembly import resolve_brdf_provider
from siac.app.planning import coerce_aoi_spec, resolve_run_config
from siac.public_models import PreparedMonthlyCompositeBuildResult

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path

    from siac.app.requests import AOISpec
    from siac.config import SIACConfig

logger = logging.getLogger(__name__)


def prepare_monthly_composites(
    config: SIACConfig,
    *,
    aoi: AOISpec,
    year_months: Sequence[tuple[int, int]],
    resolution: float | None,
    output_path: str | Path,
    auth: CredentialManager | None = None,
) -> PreparedMonthlyCompositeBuildResult:
    """Build and persist prepared monthly composites for an AOI."""
    runtime_aoi = coerce_aoi_spec(aoi)
    if runtime_aoi is None:
        raise ValueError("prepare_monthly_composites requires an AOI")
    resolved_config = resolve_run_config(
        config,
        sensor=getattr(config, "sensor", "auto"),
    )
    auth_obj = auth or CredentialManager.from_config(resolved_config)
    brdf_provider = resolve_brdf_provider(resolved_config, auth=auth_obj)
    effective_resolution = _resolve_requested_resolution(brdf_provider, resolution)
    collection = build_monthly_composites_from_brdf(
        brdf_provider=brdf_provider,
        bounds=runtime_aoi.get_bounds(),
        crs=str(runtime_aoi.crs),
        resolution=effective_resolution,
        year_months=year_months,
    )
    materialized_collection = filter_materialized_monthly_composite_collection(collection)
    requested_periods = {
        (int(year), int(month))
        for year, month in sorted((int(y), int(m)) for y, m in year_months)
    }
    written_periods = {
        (int(composite.year), int(composite.month))
        for composite in materialized_collection.composites
    }
    skipped_periods = tuple(
        f"{year:04d}-{month:02d}"
        for year, month in sorted(requested_periods - written_periods)
    )
    if skipped_periods:
        logger.info(
            "Skipping monthly composites with no valid samples: %s",
            ", ".join(skipped_periods),
        )
    if not materialized_collection.composites:
        raise ValueError(
            "No monthly composites with valid BRDF data were produced for the requested AOI and periods."
        )
    store_path = write_monthly_composite_collection(
        materialized_collection,
        output_path,
    )
    periods = tuple(
        f"{int(composite.year):04d}-{int(composite.month):02d}"
        for composite in sorted(
            materialized_collection.composites,
            key=lambda composite: (int(composite.year), int(composite.month)),
        )
    )
    representation = (
        "kernel_weights"
        if all(isinstance(composite, MonthlyKernelWeightComposite) for composite in materialized_collection.composites)
        else "reflectance"
    )
    return PreparedMonthlyCompositeBuildResult(
        store_path=store_path,
        period_count=len(periods),
        periods=periods,
        source_name=materialized_collection.source_name,
        source_band_names=tuple(band.name for band in materialized_collection.source_bands),
        representation=representation,
    )


def _resolve_requested_resolution(
    brdf_provider: object,
    resolution: float | None,
) -> float:
    if resolution is not None:
        if resolution <= 0:
            raise ValueError("resolution must be > 0")
        return float(resolution)

    source_bands = tuple(getattr(brdf_provider, "source_bands", ()) or ())
    for band in source_bands:
        band_resolution = float(getattr(band, "resolution", 0.0))
        if band_resolution > 0:
            logger.info(
                "Monthly composite generation: resolution not provided, defaulting to source dataset resolution %.3f",
                band_resolution,
            )
            return band_resolution
    raise ValueError("resolution must be provided when the BRDF provider does not expose a positive source resolution")


__all__ = [
    "PreparedMonthlyCompositeBuildResult",
    "prepare_monthly_composites",
]
