"""Prepared monthly-composite build workflow."""

from __future__ import annotations

from typing import TYPE_CHECKING

from siac.adapters.auth import CredentialManager
from siac.algorithms.surface.brdf_monthly_composite import MonthlyKernelWeightComposite
from siac.algorithms.surface.monthly_composite_store import (
    MonthlyCompositeStoreGridSpec,
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


def prepare_monthly_composites(
    config: SIACConfig,
    *,
    aoi: AOISpec,
    year_months: Sequence[tuple[int, int]],
    resolution: float,
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
    collection = build_monthly_composites_from_brdf(
        brdf_provider=brdf_provider,
        bounds=runtime_aoi.get_bounds(),
        crs=str(runtime_aoi.crs),
        resolution=resolution,
        year_months=year_months,
    )
    store_path = write_monthly_composite_collection(
        collection,
        output_path,
        grid=MonthlyCompositeStoreGridSpec.from_bounds(
            runtime_aoi.get_bounds(),
            crs=str(runtime_aoi.crs),
            resolution=resolution,
        ),
    )
    periods = tuple(f"{year:04d}-{month:02d}" for year, month in sorted((int(y), int(m)) for y, m in year_months))
    representation = (
        "kernel_weights"
        if all(isinstance(composite, MonthlyKernelWeightComposite) for composite in collection.composites)
        else "reflectance"
    )
    return PreparedMonthlyCompositeBuildResult(
        store_path=store_path,
        period_count=len(periods),
        periods=periods,
        source_name=collection.source_name,
        source_band_names=tuple(band.name for band in collection.source_bands),
        representation=representation,
    )


__all__ = [
    "PreparedMonthlyCompositeBuildResult",
    "prepare_monthly_composites",
]
