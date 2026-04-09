"""Route-B monthly-database construction and runtime query helpers."""

from __future__ import annotations

import calendar
import collections.abc
import logging
import time
from dataclasses import dataclass
from datetime import datetime
from typing import TYPE_CHECKING, Any, Protocol, TypeAlias, cast

import numpy as np
import xarray as xr
from numpy import typing as npt

from siac.algorithms.brdf.kernels import BRDFKernels, compute_reflectance
from siac.algorithms.correction.atmospheric import AtmosphericCorrector
from siac.algorithms.grid.assembler import (
    _build_target_template,
    _resample_cloud_mask,
    _resample_da,
)
from siac.algorithms.surface.brdf_monthly_composite import (
    MonthlyBestPixelComposite,
    MonthlyCompositeCollection,
    MonthlyKernelWeightComposite,
    build_monthly_best_pixel_composite,
    build_monthly_best_pixel_kernel_composite,
)
from siac.algorithms.surface.brdf_monthly_database import (
    MonthlyCompositeDatabase,
    build_monthly_composite_database,
)
from siac.algorithms.surface.spectral_mapping import (
    SpectralMapper,
    SpectralMappingConfig,
    _write_distance_metric_diagnostics,
    needs_spectral_mapping,
)
from siac.domain import SensorBand
from siac.runtime import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    ObservationBundle,
    SurfacePrior,
)
from siac.runtime.models import copy_spatial_metadata_like

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

logger = logging.getLogger(__name__)

BoolArray: TypeAlias = npt.NDArray[np.bool_]


_HISTORY_YEARS = 5
_HISTORY_MONTH_OFFSETS = (-1, 0, 1)
_WEEKLY_STEP_DAYS = 7
_BRDF_BATCH_MONTHS = 3


class _TemporalBRDFProvider(Protocol):
    source_bands: collections.abc.Sequence[SensorBand]

    def get_temporal_brdf_parameters_batch(
        self,
        *,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_times: collections.abc.Sequence[datetime],
        target_resolution: float,
        bands: collections.abc.Sequence[SensorBand],
        temporal_windows: collections.abc.Sequence[int],
        sample_date_sets: collections.abc.Sequence[collections.abc.Sequence[datetime]],
    ) -> collections.abc.Sequence[BRDFKernelWeights]: ...


@dataclass(frozen=True)
class _MonthlyModeledReflectance:
    year: int
    month: int
    reflectance: xr.DataArray
    quality: xr.DataArray
    reflectance_unc: xr.DataArray
    source_fit_rmse: xr.DataArray | None = None


@dataclass(frozen=True)
class _PreparedMonthlyComposite:
    year: int
    month: int
    reflectance: xr.DataArray
    quality: xr.DataArray
    sample_index: xr.DataArray
    reflectance_unc: xr.DataArray | None


@dataclass(frozen=True)
class _MonthSpec:
    year: int
    month: int
    center_time: datetime
    temporal_window: int
    sample_dates: tuple[datetime, ...]


def build_monthly_surface_prior_database(
    *,
    monthly_composites: MonthlyCompositeCollection | None = None,
    observation: ObservationBundle | None = None,
    brdf_provider: object | None = None,
    resolution: float | None = None,
    geometry: GeometryAngles,
    visible_bands: collections.abc.Sequence[SensorBand],
    query_bands: collections.abc.Sequence[SensorBand],
    spectral_library: SpectralMappingConfig | None = None,
    spectral_k_neighbors: int = 5,
    max_source_fit_rmse: float | None = None,
) -> MonthlyCompositeDatabase:
    """Build a Route-B database from prepared monthly composites."""
    target_bands = _deduplicate_bands([*visible_bands, *query_bands])
    if monthly_composites is None:
        if observation is None or brdf_provider is None or resolution is None:
            raise TypeError(
                "build_monthly_surface_prior_database requires monthly_composites=... "
                "or the legacy observation=..., brdf_provider=..., resolution=... inputs"
            )
        monthly_composites = generate_monthly_composites_from_brdf(
            observation=observation,
            brdf_provider=brdf_provider,
            resolution=resolution,
            fallback_source_bands=target_bands,
        )

    source_bands = tuple(monthly_composites.source_bands) or tuple(target_bands)
    if not monthly_composites.source_bands:
        logger.info(
            "Monthly surface-prior database: source band metadata missing; "
            "falling back to requested target bands"
        )

    logger.info(
        "Monthly surface-prior database build: prepared_composites=%d target_bands=%d source_bands=%d",
        len(monthly_composites.composites),
        len(target_bands),
        len(source_bands),
    )
    spectral_mapper = (
        SpectralMapper(
            source_bands,
            target_bands,
            spectral_library=spectral_library,
            k_neighbors=spectral_k_neighbors,
        )
        if needs_spectral_mapping(source_bands, target_bands)
        else None
    )
    composites = _normalize_monthly_composites_to_target_basis(
        monthly_composites.composites,
        geometry=geometry,
        target_bands=target_bands,
        spectral_mapper=spectral_mapper,
    )
    logger.info("Monthly surface-prior database: assembling final database from %d composite(s)", len(composites))
    database = build_monthly_composite_database(
        composites,
        query_bands=tuple(band.name for band in query_bands),
        visible_bands=tuple(band.name for band in visible_bands),
        max_source_fit_rmse=max_source_fit_rmse,
    )
    logger.info(
        "Monthly surface-prior database complete: entries=%d visible_bands=%d query_bands=%d",
        int(database.entries_features.shape[0]),
        len(database.visible_band_names),
        len(database.query_band_names),
    )
    return database


def generate_monthly_composites_from_brdf(
    *,
    observation: ObservationBundle,
    brdf_provider: object,
    resolution: float,
    fallback_source_bands: collections.abc.Sequence[SensorBand] = (),
) -> MonthlyCompositeCollection:
    """Generate default monthly kernel-weight composites from a temporal BRDF provider."""
    obs_time = cast("datetime", observation.metadata["observation_time"])
    return build_monthly_composites_from_brdf(
        brdf_provider=brdf_provider,
        bounds=observation.bounds,
        crs=observation.crs,
        resolution=resolution,
        year_months=_iter_history_months(obs_time),
        fallback_source_bands=fallback_source_bands,
        template_time=obs_time,
    )


def build_monthly_composites_from_brdf(
    *,
    brdf_provider: object,
    bounds: tuple[float, float, float, float],
    crs: str,
    resolution: float,
    year_months: collections.abc.Sequence[tuple[int, int]],
    fallback_source_bands: collections.abc.Sequence[SensorBand] = (),
    template_time: datetime | None = None,
) -> MonthlyCompositeCollection:
    """Build monthly kernel-weight composites for explicit year/month selections."""
    if resolution <= 0:
        raise ValueError("resolution must be > 0")

    source_bands = _resolve_provider_source_bands(brdf_provider, fallback_source_bands)
    month_specs = _build_explicit_month_specs(year_months, template_time=template_time)
    logger.info(
        "Monthly composite generation from BRDF: first_period=%s months=%d batch_size=%d source_bands=%d",
        f"{month_specs[0].year:04d}-{month_specs[0].month:02d}" if month_specs else "n/a",
        len(month_specs),
        _BRDF_BATCH_MONTHS,
        len(source_bands),
    )

    composites: list[MonthlyKernelWeightComposite] = []
    month_batches = [
        (batch_start, month_specs[batch_start: batch_start + _BRDF_BATCH_MONTHS])
        for batch_start in range(0, len(month_specs), _BRDF_BATCH_MONTHS)
    ]
    for batch_start, month_batch in month_batches:
        logger.info(
            "Monthly composite generation batch %d-%d/%d: fetching temporal BRDF for months=%s",
            batch_start + 1,
            batch_start + len(month_batch),
            len(month_specs),
            [f"{spec.year:04d}-{spec.month:02d}" for spec in month_batch],
        )
        temporal_weights_batch = _get_temporal_brdf_parameters_batch(
            brdf_provider,
            bounds=bounds,
            crs=crs,
            obs_times=[spec.center_time for spec in month_batch],
            target_resolution=resolution,
            bands=source_bands,
            temporal_windows=[spec.temporal_window for spec in month_batch],
            sample_date_sets=[spec.sample_dates for spec in month_batch],
        )
        for spec, temporal_weights in zip(month_batch, temporal_weights_batch, strict=True):
            monthly_weights = _select_temporal_weights_for_month(
                temporal_weights,
                year=spec.year,
                month=spec.month,
            )
            quality = _quality_from_temporal_weights(monthly_weights)
            monthly_weights, quality, collapsed = _collapse_repeated_temporal_weights(monthly_weights, quality)
            if collapsed:
                logger.info(
                    "Monthly composite generation month %04d-%02d: collapsed identical sampled days to %d sample",
                    spec.year,
                    spec.month,
                    int(monthly_weights.f0.sizes["time"]),
                )
            composites.append(
                build_monthly_best_pixel_kernel_composite(
                    monthly_weights,
                    quality,
                    year=spec.year,
                    month=spec.month,
                )
            )

    return MonthlyCompositeCollection(
        composites=tuple(composites),
        source_bands=tuple(source_bands),
        source_name=getattr(brdf_provider, "source_name", None),
    )


def _fetch_monthly_batch_inputs(
    brdf_provider: object,
    *,
    observation: ObservationBundle,
    resolution: float,
    geometry: GeometryAngles,
    source_bands: collections.abc.Sequence[SensorBand],
    month_batch: collections.abc.Sequence[_MonthSpec],
    total_months: int,
    batch_start: int,
) -> list[_MonthlyModeledReflectance]:
    logger.info(
        "Monthly surface-prior database batch %d-%d/%d: fetching temporal BRDF for months=%s",
        batch_start + 1,
        batch_start + len(month_batch),
        total_months,
        [f"{spec.year:04d}-{spec.month:02d}" for spec in month_batch],
    )
    temporal_weights_batch = _get_temporal_brdf_parameters_batch(
        brdf_provider,
        bounds=observation.bounds,
        crs=observation.crs,
        obs_times=[spec.center_time for spec in month_batch],
        target_resolution=resolution,
        bands=source_bands,
        temporal_windows=[spec.temporal_window for spec in month_batch],
        sample_date_sets=[spec.sample_dates for spec in month_batch],
    )
    logger.info(
        "Monthly surface-prior database batch %d-%d/%d: received %d temporal BRDF stack(s)",
        batch_start + 1,
        batch_start + len(month_batch),
        total_months,
        len(temporal_weights_batch),
    )

    monthly_inputs: list[_MonthlyModeledReflectance] = []
    for spec, temporal_weights in zip(month_batch, temporal_weights_batch, strict=True):
        logger.info(
            "Monthly surface-prior database month %04d-%02d: forward modelling %d sampled day(s)",
            spec.year,
            spec.month,
            int(temporal_weights.f0.sizes["time"]),
        )
        reflectance, quality, reflectance_unc = _forward_model_monthly_reflectance(
            temporal_weights,
            geometry=geometry,
            year=spec.year,
            month=spec.month,
        )
        reflectance, quality, reflectance_unc, collapsed = _collapse_repeated_temporal_samples(
            reflectance,
            quality,
            reflectance_unc,
        )
        if collapsed:
            logger.info(
                "Monthly surface-prior database month %04d-%02d: collapsed identical sampled days to %d sample",
                spec.year,
                spec.month,
                int(reflectance.sizes["time"]),
            )
        monthly_inputs.append(
            _MonthlyModeledReflectance(
                year=spec.year,
                month=spec.month,
                reflectance=reflectance,
                quality=quality,
                reflectance_unc=reflectance_unc,
            )
        )
    return monthly_inputs


def _build_monthly_composites(
    monthly_inputs: collections.abc.Sequence[_MonthlyModeledReflectance],
) -> list[MonthlyBestPixelComposite]:
    composites: list[MonthlyBestPixelComposite] = []
    for monthly_input in monthly_inputs:
        logger.info(
            "Monthly surface-prior database month %04d-%02d: building best-pixel composite from %d sample(s)",
            monthly_input.year,
            monthly_input.month,
            int(monthly_input.reflectance.sizes["time"]),
        )
        composites.append(
            build_monthly_best_pixel_composite(
                monthly_input.reflectance,
                monthly_input.quality,
                year=monthly_input.year,
                month=monthly_input.month,
            )
        )
        logger.info(
            "Monthly surface-prior database month %04d-%02d: composite complete",
            monthly_input.year,
            monthly_input.month,
        )
    return composites


def query_surface_prior_from_monthly_database(
    *,
    observation: ObservationBundle,
    atmo_prior: AtmosphericState,
    rt_model: object,
    database: MonthlyCompositeDatabase,
    query_band_names: tuple[str, ...] | None = None,
    visible_band_names: tuple[str, ...] | None = None,
    target_resolution: float | None = None,
    k_neighbors: int = 3,
    max_prediction_uncertainty: float | None = 0.05,
    max_composite_quality: float | None = 0.05,
    max_knn_feature_distance: float | None = 0.05,
    diagnostic_cache_dir: Path | str | None = None,
) -> SurfacePrior:
    """Build a visible-band surface prior from first-pass corrected NIR/SWIR."""
    expected_query = tuple(database.query_band_names)
    if query_band_names is not None and tuple(query_band_names) != expected_query:
        raise ValueError("query_band_names must match the database query-band ordering")

    expected_visible = tuple(database.visible_band_names)
    if visible_band_names is not None and tuple(visible_band_names) != expected_visible:
        raise ValueError("visible_band_names must match the database visible-band ordering")

    query_template = _database_query_template(database)
    if target_resolution is not None:
        resolution = float(target_resolution)
        if not np.isfinite(resolution) or resolution <= 0.0:
            raise ValueError("target_resolution must be a finite positive value")
        query_template = _build_target_template(observation.bounds, observation.crs, resolution)
    target_shape = (
        int(query_template.sizes["y"]),
        int(query_template.sizes["x"]),
    )
    native_query_valid = _query_observation_valid_mask(
        observation.toa,
        band_names=expected_query,
        cloud_mask=observation.cloud_mask,
    )
    # Pre-warm the KDTree in a background thread so it's ready when needed.
    import concurrent.futures

    def _prewarm_tree(db: object) -> None:
        idx = getattr(db, "_neighbor_index", None)
        if idx is None and hasattr(db, "predict_visible_with_diagnostics"):
            # Trigger any lazy property on the underlying database.
            inner = getattr(db, "_inner", None) or getattr(db, "_database", None)
            if inner is not None:
                getattr(inner, "_neighbor_index", None)

    _tree_pool = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    try:
        _tree_future = _tree_pool.submit(_prewarm_tree, database)
    except Exception:
        _tree_pool.shutdown(wait=False)
        raise

    _t0 = time.perf_counter()
    coarse_query_toa, coarse_query_valid = _resample_dataset_with_validity(
        observation.toa,
        band_names=expected_query,
        valid_mask=native_query_valid,
        target_shape=target_shape,
        template=query_template,
        cloud_mask=observation.cloud_mask,
    )
    logger.info("M3 timing: _resample_dataset_with_validity %.3f s", time.perf_counter() - _t0)
    _t0 = time.perf_counter()
    coarse_geometry = _resample_geometry_to_target_shape(observation.geometry, target_shape)
    coarse_atmo = _resample_atmo_to_target_shape(atmo_prior, target_shape)
    coarse_cloud_mask = _resample_cloud_mask_to_target_shape(observation.cloud_mask, target_shape)
    logger.info("M3 timing: resample geometry/atmo/cloud %.3f s", time.perf_counter() - _t0)
    coarse_invalid = coarse_cloud_mask | (~coarse_query_valid)
    corrector = AtmosphericCorrector(rt_model, observation.sensor_config)
    _t0 = time.perf_counter()
    correction = corrector.correct(
        coarse_query_toa,
        coarse_geometry,
        coarse_atmo,
        cloud_mask=coarse_invalid,
    )
    logger.info("M3 timing: AtmosphericCorrector.correct (first-pass) %.3f s", time.perf_counter() - _t0)
    corrected_query_mask = _resample_cloud_mask_to_target_shape(correction.cloud_mask, target_shape)
    corrected_query = _apply_invalid_mask_to_dataset(
        _resample_dataset(
            correction.boa,
            band_names=expected_query,
            target_shape=target_shape,
            template=query_template,
        ),
        band_names=expected_query,
        invalid_mask=corrected_query_mask,
    )
    if hasattr(database, "predict_visible_with_diagnostics"):
        # Ensure KDTree is built before entering the prediction loop.
        _tree_future.result()
        _t0 = time.perf_counter()
        prediction = database.predict_visible_with_diagnostics(
            corrected_query,
            k_neighbors=k_neighbors,
        )
        logger.info("M3 timing: predict_visible_with_diagnostics %.3f s", time.perf_counter() - _t0)
        predicted_visible = prediction.predicted
        predicted_unc = prediction.uncertainty
        predicted_quality = prediction.quality
        predicted_source_fit = prediction.source_fit_rmse
        predicted_distance = prediction.knn_feature_distance
    else:
        _tree_future.result()
        predicted_visible, predicted_unc, predicted_quality = database.predict_visible(
            corrected_query,
            k_neighbors=k_neighbors,
        )
        predicted_source_fit = xr.zeros_like(
            cast("xr.DataArray", predicted_quality),
            dtype=np.float32,
        )
        predicted_distance = xr.zeros_like(
            cast("xr.DataArray", predicted_quality),
            dtype=np.float32,
        )

    _tree_pool.shutdown(wait=True)
    spatial_reference = cast("xr.DataArray", predicted_visible.isel(band=0, drop=True))
    predicted_source_fit = copy_spatial_metadata_like(predicted_source_fit.astype(np.float32), spatial_reference)
    predicted_distance = copy_spatial_metadata_like(predicted_distance.astype(np.float32), spatial_reference)

    _write_distance_metric_diagnostics(
        diagnostic_cache_dir,
        prefix="swir_refine_distances",
        metrics={
            "source_fit_rmse": predicted_source_fit,
            "knn_feature_distance": predicted_distance,
        },
        metadata={
            "query_band_names": list(expected_query),
            "visible_band_names": list(expected_visible),
            "k_neighbors": int(k_neighbors),
            "target_shape": [int(size) for size in target_shape],
        },
    )

    cloud_mask = corrected_query_mask
    uncertainty_ok: BoolArray = np.ones(target_shape, dtype=bool)
    if max_prediction_uncertainty is not None:
        uncertainty_ok = cast(
            "BoolArray",
            np.all(predicted_unc.values <= float(max_prediction_uncertainty), axis=0),
        )
    quality_ok: BoolArray = np.ones(target_shape, dtype=bool)
    if max_composite_quality is not None:
        quality_ok = cast("BoolArray", predicted_quality.values <= float(max_composite_quality))
    distance_ok: BoolArray = np.ones(target_shape, dtype=bool)
    if max_knn_feature_distance is not None:
        distance_ok = cast("BoolArray", predicted_distance.values <= float(max_knn_feature_distance))
    valid = (
        np.all(np.isfinite(predicted_visible.values), axis=0)
        & np.all(np.isfinite(predicted_unc.values), axis=0)
        & np.isfinite(predicted_quality.values)
        & np.isfinite(predicted_distance.values)
        & coarse_query_valid.values
        & uncertainty_ok
        & quality_ok
        & distance_ok
    )
    mask = xr.DataArray(
        valid & (~cloud_mask.values.astype(bool)),
        dims=["y", "x"],
        coords={"y": predicted_visible.coords["y"], "x": predicted_visible.coords["x"]},
    )

    return SurfacePrior(
        boa=predicted_visible,
        boa_unc=predicted_unc,
        kernels=None,
        mask=mask,
    )


def resample_geometry_for_surface_prior(
    observation: ObservationBundle,
    *,
    resolution: float,
) -> GeometryAngles:
    """Resample geometry to the Route-B target grid."""
    target_template = _build_target_template(observation.bounds, observation.crs, resolution)
    target_shape: tuple[int, int] = (
        int(target_template.shape[0]),
        int(target_template.shape[1]),
    )
    return GeometryAngles(
        sza=_resample_da(observation.geometry.sza, target_shape, "bilinear", template=target_template),
        saa=_resample_da(observation.geometry.saa, target_shape, "bilinear", template=target_template),
        vza=_resample_da(observation.geometry.vza, target_shape, "bilinear", template=target_template),
        vaa=_resample_da(observation.geometry.vaa, target_shape, "bilinear", template=target_template),
    )


def _iter_history_months(obs_time: datetime) -> list[tuple[int, int]]:
    months: list[tuple[int, int]] = []
    for year_offset in range(1, _HISTORY_YEARS + 1):
        base_year = obs_time.year - year_offset
        for month_offset in _HISTORY_MONTH_OFFSETS:
            month = obs_time.month + month_offset
            year = base_year
            while month < 1:
                month += 12
                year -= 1
            while month > 12:
                month -= 12
                year += 1
            months.append((year, month))
    return months


def _build_explicit_month_specs(
    year_months: collections.abc.Sequence[tuple[int, int]],
    *,
    template_time: datetime | None,
) -> list[_MonthSpec]:
    if not year_months:
        raise ValueError("year_months must not be empty")

    specs: list[_MonthSpec] = []
    seen: set[tuple[int, int]] = set()
    for year, month in sorted((int(year), int(month)) for year, month in year_months):
        if not 1 <= month <= 12:
            raise ValueError(f"month must be between 1 and 12, got {month}")
        key = (year, month)
        if key in seen:
            raise ValueError(f"Duplicate monthly composite selection: {year:04d}-{month:02d}")
        seen.add(key)
        center_template = template_time or datetime(year, month, 15, 12, 0, 0)
        specs.append(
            _MonthSpec(
                year=year,
                month=month,
                center_time=_month_center_datetime(year, month, center_template),
                temporal_window=max(1, calendar.monthrange(year, month)[1] // 2 + 1),
                sample_dates=_weekly_sample_dates(year, month),
            )
        )
    return specs


def _month_center_datetime(year: int, month: int, template_time: datetime) -> datetime:
    n_days = calendar.monthrange(year, month)[1]
    center_day = int(np.ceil(n_days / 2.0))
    return template_time.replace(year=year, month=month, day=center_day)


def _weekly_sample_dates(year: int, month: int) -> tuple[datetime, ...]:
    n_days = calendar.monthrange(year, month)[1]
    days = list(range(1, n_days + 1, _WEEKLY_STEP_DAYS))
    if days[-1] != n_days and (n_days - days[-1]) >= (_WEEKLY_STEP_DAYS // 2):
        days.append(n_days)
    return tuple(datetime(year, month, day) for day in days)


def _forward_model_monthly_reflectance(
    temporal_weights: BRDFKernelWeights,
    *,
    geometry: GeometryAngles,
    year: int,
    month: int,
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    logger.info(
        "Monthly forward model start: month=%04d-%02d time=%d bands=%d grid=%dx%d",
        year,
        month,
        int(temporal_weights.f0.sizes["time"]),
        int(temporal_weights.f0.sizes["band"]),
        int(temporal_weights.f0.sizes["y"]),
        int(temporal_weights.f0.sizes["x"]),
    )
    kernels = BRDFKernels(hb=2.0, br=1.0)
    k_vol, k_geo = kernels.compute(geometry.vza, geometry.sza, geometry.raa)
    reflectance = cast(
        "xr.DataArray",
        compute_reflectance(
            temporal_weights.f0,
            temporal_weights.f1,
            temporal_weights.f2,
            cast("xr.DataArray", k_vol),
            cast("xr.DataArray", k_geo),
        ),
    ).transpose("time", "band", "y", "x")
    reflectance_unc = temporal_weights.compute_reflectance_uncertainty(
        cast("xr.DataArray", k_vol),
        cast("xr.DataArray", k_geo),
    ).transpose("time", "band", "y", "x")

    month_mask = _select_month_mask(temporal_weights.f0.coords["time"].values, year=year, month=month)
    if month_mask.any():
        reflectance = reflectance.isel(time=month_mask)
        reflectance_unc = reflectance_unc.isel(time=month_mask)

    quality = reflectance_unc.mean(dim="band", skipna=True).astype(np.float32)
    logger.info(
        "Monthly forward model complete: month=%04d-%02d retained_samples=%d",
        year,
        month,
        int(reflectance.sizes["time"]),
    )
    return reflectance.astype(np.float32), quality, reflectance_unc.astype(np.float32)


def _select_month_mask(time_values: np.ndarray, *, year: int, month: int) -> np.ndarray:
    month_strings = np.asarray(time_values, dtype="datetime64[D]").astype("datetime64[M]")
    target = np.datetime64(f"{year:04d}-{month:02d}", "M")
    month_mask = cast(
        "np.ndarray[Any, np.dtype[np.bool_]]",
        np.asarray(month_strings == target, dtype=np.bool_),
    )
    return month_mask


def _collapse_repeated_temporal_samples(
    reflectance: xr.DataArray,
    quality: xr.DataArray,
    reflectance_unc: xr.DataArray,
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray, bool]:
    if int(reflectance.sizes.get("time", 1)) <= 1:
        return reflectance, quality, reflectance_unc, False
    if not _all_time_slices_identical(reflectance):
        return reflectance, quality, reflectance_unc, False
    if not _all_time_slices_identical(quality):
        return reflectance, quality, reflectance_unc, False
    if not _all_time_slices_identical(reflectance_unc):
        return reflectance, quality, reflectance_unc, False
    first = {"time": slice(0, 1)}
    return (
        reflectance.isel(first),
        quality.isel(first),
        reflectance_unc.isel(first),
        True,
    )


def _all_time_slices_identical(data: xr.DataArray) -> bool:
    if int(data.sizes.get("time", 1)) <= 1:
        return True
    values = np.asarray(data.values)
    template = np.broadcast_to(values[0:1], values.shape)
    return bool(np.allclose(values, template, rtol=0.0, atol=0.0, equal_nan=True))


def _map_monthly_inputs_to_target_basis(
    monthly_inputs: collections.abc.Sequence[_MonthlyModeledReflectance],
    spectral_mapper: SpectralMapper,
) -> list[_MonthlyModeledReflectance]:
    if not monthly_inputs:
        return []

    stacked_reflectance = xr.concat(
        [monthly_input.reflectance for monthly_input in monthly_inputs],
        dim="time",
    ).transpose("time", "band", "y", "x")
    stacked_uncertainty = xr.concat(
        [monthly_input.reflectance_unc for monthly_input in monthly_inputs],
        dim="time",
    ).transpose("time", "band", "y", "x")
    mapped_reflectance, mapped_uncertainty, mapped_source_fit_rmse = spectral_mapper.map(
        stacked_reflectance,
        source_uncertainty=stacked_uncertainty,
    )

    remapped: list[_MonthlyModeledReflectance] = []
    start = 0
    for monthly_input in monthly_inputs:
        stop = start + int(monthly_input.reflectance.sizes["time"])
        mapped_month = mapped_reflectance.isel(time=slice(start, stop))
        mapped_unc_month = mapped_uncertainty.isel(time=slice(start, stop))
        mapped_fit_month = mapped_source_fit_rmse.isel(time=slice(start, stop))
        remapped.append(
            _MonthlyModeledReflectance(
                year=monthly_input.year,
                month=monthly_input.month,
                reflectance=mapped_month,
                quality=_combine_quality_with_mapping_uncertainty(
                    monthly_input.quality,
                    mapped_unc_month,
                    mapped_fit_month,
                ),
                reflectance_unc=mapped_unc_month,
                source_fit_rmse=mapped_fit_month.mean(dim="time", skipna=True).astype(np.float32),
            )
        )
        start = stop
    return remapped


def _combine_quality_with_mapping_uncertainty(
    quality: xr.DataArray,
    mapped_uncertainty: xr.DataArray,
    mapped_source_fit_rmse: xr.DataArray,
) -> xr.DataArray:
    combined = np.sqrt(
        np.square(quality.values, dtype=np.float32)
        + np.square(mapped_uncertainty.mean(dim="band", skipna=True).values, dtype=np.float32)
        + np.square(mapped_source_fit_rmse.values, dtype=np.float32)
    ).astype(np.float32)
    return xr.DataArray(
        combined,
        dims=["time", "y", "x"],
        coords={
            "time": quality.coords["time"],
            "y": quality.coords["y"],
            "x": quality.coords["x"],
        },
    )


def _combine_composite_quality_with_mapping_uncertainty(
    quality: xr.DataArray,
    mapped_uncertainty: xr.DataArray,
    mapped_source_fit_rmse: xr.DataArray,
) -> xr.DataArray:
    combined = np.sqrt(
        np.square(quality.values, dtype=np.float32)
        + np.square(mapped_uncertainty.mean(dim="band", skipna=True).values, dtype=np.float32)
        + np.square(mapped_source_fit_rmse.values, dtype=np.float32)
    ).astype(np.float32)
    return xr.DataArray(
        combined,
        dims=["y", "x"],
        coords={
            "y": quality.coords["y"],
            "x": quality.coords["x"],
        },
    )


def _normalize_monthly_composite_to_target_basis(
    composite: MonthlyBestPixelComposite | MonthlyKernelWeightComposite,
    *,
    geometry: GeometryAngles,
) -> _PreparedMonthlyComposite:
    template = geometry.sza
    downsample_method = _monthly_composite_downsample_method(composite, template)
    quality = _resample_spatial_field_to_template(
        composite.quality.astype(np.float32),
        template,
        downsample_method,
    )
    sample_index = _resample_spatial_field_to_template(
        composite.sample_index.astype(np.float32),
        template,
        "nearest",
    ).round().astype(np.int16)

    if isinstance(composite, MonthlyKernelWeightComposite):
        weights = _resample_brdf_weights_to_template(composite.kernels, template, method=downsample_method)
        reflectance, reflectance_unc = _reflectance_from_kernel_weights(weights, geometry)
    else:
        reflectance = _resample_band_cube_to_template(
            composite.reflectance.astype(np.float32),
            template,
            downsample_method,
        )
        reflectance_unc = None

    if reflectance_unc is not None:
        reflectance_unc = reflectance_unc.astype(np.float32)

    return _PreparedMonthlyComposite(
        reflectance=reflectance.astype(np.float32),
        quality=quality.astype(np.float32),
        sample_index=sample_index.astype(np.int16),
        year=composite.year,
        month=composite.month,
        reflectance_unc=reflectance_unc,
    )


def _normalize_monthly_composites_to_target_basis(
    composites: collections.abc.Sequence[MonthlyBestPixelComposite | MonthlyKernelWeightComposite],
    *,
    geometry: GeometryAngles,
    target_bands: collections.abc.Sequence[SensorBand],
    spectral_mapper: SpectralMapper | None,
) -> list[MonthlyBestPixelComposite]:
    prepared = [
        _normalize_monthly_composite_to_target_basis(
            composite,
            geometry=geometry,
        )
        for composite in composites
    ]
    if not prepared:
        return []

    if spectral_mapper is None:
        return [
            MonthlyBestPixelComposite(
                reflectance=prepared_composite.reflectance.sel(band=[band.name for band in target_bands]).astype(np.float32),
                quality=prepared_composite.quality.astype(np.float32),
                sample_index=prepared_composite.sample_index.astype(np.int16),
                year=prepared_composite.year,
                month=prepared_composite.month,
            )
            for prepared_composite in prepared
        ]

    remapped: list[MonthlyBestPixelComposite | None] = [None] * len(prepared)
    groups: dict[bool, list[tuple[int, _PreparedMonthlyComposite]]] = {False: [], True: []}
    for index, prepared_composite in enumerate(prepared):
        groups[prepared_composite.reflectance_unc is not None].append((index, prepared_composite))

    for has_uncertainty, group in groups.items():
        if not group:
            continue
        stacked_reflectance = xr.concat(
            [entry.reflectance for _, entry in group],
            dim=xr.IndexVariable("time", np.arange(len(group), dtype=np.int32)),
        ).transpose("time", "band", "y", "x")
        stacked_uncertainty = None
        if has_uncertainty:
            stacked_uncertainty = xr.concat(
                [cast("xr.DataArray", entry.reflectance_unc) for _, entry in group],
                dim=xr.IndexVariable("time", np.arange(len(group), dtype=np.int32)),
            ).transpose("time", "band", "y", "x")

        mapped_reflectance, mapped_uncertainty, mapped_source_fit_rmse = spectral_mapper.map(
            stacked_reflectance,
            source_uncertainty=stacked_uncertainty,
        )
        for time_index, (result_index, entry) in enumerate(group):
            mapped_reflectance_entry = mapped_reflectance.isel(time=time_index, drop=True).astype(np.float32)
            mapped_uncertainty_entry = mapped_uncertainty.isel(time=time_index, drop=True).astype(np.float32)
            mapped_source_fit_entry = mapped_source_fit_rmse.isel(time=time_index, drop=True).astype(np.float32)
            remapped[result_index] = MonthlyBestPixelComposite(
                reflectance=mapped_reflectance_entry,
                quality=_combine_composite_quality_with_mapping_uncertainty(
                    entry.quality,
                    mapped_uncertainty_entry,
                    mapped_source_fit_entry,
                ).astype(np.float32),
                source_fit_rmse=mapped_source_fit_entry,
                sample_index=entry.sample_index.astype(np.int16),
                year=entry.year,
                month=entry.month,
            )

    return [cast("MonthlyBestPixelComposite", composite) for composite in remapped]


def _reflectance_from_kernel_weights(
    weights: BRDFKernelWeights,
    geometry: GeometryAngles,
) -> tuple[xr.DataArray, xr.DataArray]:
    kernels = BRDFKernels(hb=2.0, br=1.0)
    k_vol, k_geo = kernels.compute(geometry.vza, geometry.sza, geometry.raa)
    reflectance = weights.compute_reflectance(
        cast("xr.DataArray", k_vol),
        cast("xr.DataArray", k_geo),
    ).transpose("band", "y", "x")
    reflectance_unc = weights.compute_reflectance_uncertainty(
        cast("xr.DataArray", k_vol),
        cast("xr.DataArray", k_geo),
    ).transpose("band", "y", "x")
    return reflectance.astype(np.float32), reflectance_unc.astype(np.float32)


def _resample_brdf_weights_to_template(
    weights: BRDFKernelWeights,
    template: xr.DataArray,
    *,
    method: str,
) -> BRDFKernelWeights:
    return BRDFKernelWeights(
        f0=_resample_band_cube_to_template(weights.f0, template, method),
        f1=_resample_band_cube_to_template(weights.f1, template, method),
        f2=_resample_band_cube_to_template(weights.f2, template, method),
        f0_unc=_resample_band_cube_to_template(weights.f0_unc, template, method),
        f1_unc=_resample_band_cube_to_template(weights.f1_unc, template, method),
        f2_unc=_resample_band_cube_to_template(weights.f2_unc, template, method),
        reflectance_unc=(
            _resample_band_cube_to_template(weights.reflectance_unc, template, method)
            if weights.reflectance_unc is not None
            else None
        ),
    )


def _monthly_composite_downsample_method(
    composite: MonthlyBestPixelComposite | MonthlyKernelWeightComposite,
    template: xr.DataArray,
) -> str:
    if isinstance(composite, MonthlyKernelWeightComposite):
        source = composite.kernels.f0
    else:
        source = composite.reflectance
    if _is_coarser_target_grid(source, template):
        return "area"
    return "bilinear"


def _is_coarser_target_grid(
    source: xr.DataArray,
    template: xr.DataArray,
) -> bool:
    source_height = int(source.sizes.get("y", 0))
    source_width = int(source.sizes.get("x", 0))
    target_height = int(template.sizes.get("y", 0))
    target_width = int(template.sizes.get("x", 0))
    return source_height > target_height or source_width > target_width


def _resample_band_cube_to_template(
    data: xr.DataArray,
    template: xr.DataArray,
    method: str,
) -> xr.DataArray:
    if _shares_spatial_grid(data, template):
        return copy_spatial_metadata_like(data.astype(np.float32), template)

    target_shape = (int(template.sizes["y"]), int(template.sizes["x"]))
    band_coords = data.coords["band"].values if "band" in data.coords else np.arange(data.sizes["band"])
    resampled = xr.concat(
        [
            _resample_da(data.sel(band=band, drop=True), target_shape, method, template=template)
            for band in band_coords
        ],
        dim=xr.IndexVariable("band", band_coords),
    )
    coords: dict[str, object] = {"band": band_coords}
    if "y" in template.coords:
        coords["y"] = template.coords["y"]
    if "x" in template.coords:
        coords["x"] = template.coords["x"]
    return copy_spatial_metadata_like(resampled.assign_coords(**coords).astype(np.float32), template)


def _resample_spatial_field_to_template(
    data: xr.DataArray,
    template: xr.DataArray,
    method: str,
) -> xr.DataArray:
    if _shares_spatial_grid(data, template):
        return copy_spatial_metadata_like(data.astype(np.float32), template)
    target_shape = (int(template.sizes["y"]), int(template.sizes["x"]))
    resampled = _resample_da(data, target_shape, method, template=template)
    coords: dict[str, object] = {}
    if "y" in template.coords:
        coords["y"] = template.coords["y"]
    if "x" in template.coords:
        coords["x"] = template.coords["x"]
    return copy_spatial_metadata_like(resampled.assign_coords(**coords).astype(np.float32), template)


def _shares_spatial_grid(data: xr.DataArray, template: xr.DataArray) -> bool:
    if data.sizes.get("y") != template.sizes.get("y") or data.sizes.get("x") != template.sizes.get("x"):
        return False
    if "y" in data.coords and "y" in template.coords and not np.array_equal(
        np.asarray(data.coords["y"].values),
        np.asarray(template.coords["y"].values),
    ):
        return False
    if "x" in data.coords and "x" in template.coords and not np.array_equal(
        np.asarray(data.coords["x"].values),
        np.asarray(template.coords["x"].values),
    ):
        return False
    return not (
        "x" in data.coords
        and "x" in template.coords
        and not np.array_equal(
            np.asarray(data.coords["x"].values),
            np.asarray(template.coords["x"].values),
        )
    )


def _select_temporal_weights_for_month(
    temporal_weights: BRDFKernelWeights,
    *,
    year: int,
    month: int,
) -> BRDFKernelWeights:
    month_mask = _select_month_mask(temporal_weights.f0.coords["time"].values, year=year, month=month)
    if not month_mask.any():
        return temporal_weights
    return _select_temporal_weight_indexer(temporal_weights, month_mask)


def _quality_from_temporal_weights(
    temporal_weights: BRDFKernelWeights,
) -> xr.DataArray:
    if temporal_weights.reflectance_unc is not None:
        return temporal_weights.reflectance_unc.mean(dim="band", skipna=True).astype(np.float32)
    stacked = xr.concat(
        [temporal_weights.f0_unc, temporal_weights.f1_unc, temporal_weights.f2_unc],
        dim=xr.IndexVariable("parameter", ["f0", "f1", "f2"]),
    )
    return stacked.mean(dim=("parameter", "band"), skipna=True).astype(np.float32)


def _collapse_repeated_temporal_weights(
    temporal_weights: BRDFKernelWeights,
    quality: xr.DataArray,
) -> tuple[BRDFKernelWeights, xr.DataArray, bool]:
    if int(temporal_weights.f0.sizes.get("time", 1)) <= 1:
        return temporal_weights, quality, False
    if not _all_time_slices_identical(temporal_weights.f0):
        return temporal_weights, quality, False
    if not _all_time_slices_identical(temporal_weights.f1):
        return temporal_weights, quality, False
    if not _all_time_slices_identical(temporal_weights.f2):
        return temporal_weights, quality, False
    if not _all_time_slices_identical(temporal_weights.f0_unc):
        return temporal_weights, quality, False
    if not _all_time_slices_identical(temporal_weights.f1_unc):
        return temporal_weights, quality, False
    if not _all_time_slices_identical(temporal_weights.f2_unc):
        return temporal_weights, quality, False
    if temporal_weights.reflectance_unc is not None and not _all_time_slices_identical(temporal_weights.reflectance_unc):
        return temporal_weights, quality, False
    if not _all_time_slices_identical(quality):
        return temporal_weights, quality, False
    first = slice(0, 1)
    return (
        _select_temporal_weight_indexer(temporal_weights, first),
        quality.isel(time=first),
        True,
    )


def _select_temporal_weight_indexer(
    temporal_weights: BRDFKernelWeights,
    indexer: Any,
) -> BRDFKernelWeights:
    return BRDFKernelWeights(
        f0=temporal_weights.f0.isel(time=indexer),
        f1=temporal_weights.f1.isel(time=indexer),
        f2=temporal_weights.f2.isel(time=indexer),
        f0_unc=temporal_weights.f0_unc.isel(time=indexer),
        f1_unc=temporal_weights.f1_unc.isel(time=indexer),
        f2_unc=temporal_weights.f2_unc.isel(time=indexer),
        reflectance_unc=(
            temporal_weights.reflectance_unc.isel(time=indexer)
            if temporal_weights.reflectance_unc is not None
            else None
        ),
    )


def _deduplicate_bands(bands: collections.abc.Sequence[SensorBand]) -> list[SensorBand]:
    seen: set[str] = set()
    ordered: list[SensorBand] = []
    for band in bands:
        if band.name in seen:
            continue
        seen.add(band.name)
        ordered.append(band)
    return ordered


def _native_observation_resolution(observation: ObservationBundle) -> float:
    available = {
        band.name: band.resolution
        for band in observation.sensor_config.bands
        if band.name in observation.toa.data_vars
    }
    if not available:
        return 10.0
    return float(min(available.values()))


def _observation_shape(observation: ObservationBundle) -> tuple[int, int]:
    first_var = next(iter(observation.toa.data_vars))
    shape = observation.toa[first_var].shape
    if len(shape) != 2:
        raise ValueError("Observation TOA variables must be 2-D over y/x")
    height, width = shape
    return int(height), int(width)


def _resample_atmo_to_observation_grid(
    observation: ObservationBundle,
    atmo_prior: AtmosphericState,
) -> AtmosphericState:
    return _resample_atmo_to_target_shape(atmo_prior, _observation_shape(observation))


def _resample_geometry_to_target_shape(
    geometry: GeometryAngles,
    target_shape: tuple[int, int],
) -> GeometryAngles:
    from concurrent.futures import ThreadPoolExecutor

    if geometry.sza.shape == target_shape:
        return geometry
    fields = {"sza": geometry.sza, "saa": geometry.saa, "vza": geometry.vza, "vaa": geometry.vaa}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {name: pool.submit(_resample_da, da, target_shape, "bilinear") for name, da in fields.items()}
        results = {name: fut.result() for name, fut in futures.items()}
    return GeometryAngles(**results)


def _resample_atmo_to_target_shape(
    atmo_prior: AtmosphericState,
    target_shape: tuple[int, int],
) -> AtmosphericState:
    from concurrent.futures import ThreadPoolExecutor

    if atmo_prior.aot.shape == target_shape:
        return atmo_prior
    fields = {
        "aot": atmo_prior.aot, "tcwv": atmo_prior.tcwv, "tco3": atmo_prior.tco3,
        "aot_unc": atmo_prior.aot_unc, "tcwv_unc": atmo_prior.tcwv_unc,
        "tco3_unc": atmo_prior.tco3_unc, "elevation": atmo_prior.elevation,
    }
    with ThreadPoolExecutor(max_workers=7) as pool:
        futures = {name: pool.submit(_resample_da, da, target_shape, "bilinear") for name, da in fields.items()}
        results = {name: fut.result() for name, fut in futures.items()}
    return AtmosphericState(**results)


def _resample_cloud_mask_to_target_shape(
    cloud_mask: xr.DataArray,
    target_shape: tuple[int, int],
) -> xr.DataArray:
    if "band" in cloud_mask.dims:
        cloud_mask = cloud_mask.any(dim="band")
    if cloud_mask.shape == target_shape:
        return cloud_mask
    return _resample_cloud_mask(cloud_mask, target_shape)


def _query_observation_valid_mask(
    dataset: xr.Dataset,
    *,
    band_names: collections.abc.Sequence[str],
    cloud_mask: xr.DataArray,
) -> xr.DataArray:
    """Compute a boolean valid-pixel mask at the cloud-mask resolution.

    Only bands that are already present in *dataset.data_vars* (i.e. at the
    same spatial grid as the cloud mask) are incorporated.  Bands that must be
    loaded lazily at a different native resolution are handled separately
    inside :func:`_resample_dataset_with_validity`.
    """
    valid = ~cloud_mask.values.astype(bool)
    for name in band_names:
        if name not in dataset.data_vars:
            continue
        band = dataset[name].values
        valid &= np.isfinite(band) & (band > 0.0) & (band < 1.0)
    return xr.DataArray(valid, dims=cloud_mask.dims, coords=cloud_mask.coords)


_TOA_BAND_LOADER_ATTR = "_siac_toa_band_loader"


def _resample_dataset_with_validity(
    dataset: xr.Dataset,
    *,
    band_names: collections.abc.Sequence[str],
    valid_mask: xr.DataArray,
    target_shape: tuple[int, int],
    template: xr.DataArray | None = None,
    cloud_mask: xr.DataArray | None = None,
) -> tuple[xr.Dataset, xr.DataArray]:
    """Resample bands to *target_shape* honouring a per-pixel validity mask.

    Bands already in *dataset.data_vars* (same grid as *valid_mask*) are
    processed using the precomputed mask.  Bands **not** in data_vars are
    fetched via the TOA band loader at their native resolution; a per-band
    validity mask is computed independently at that resolution so that
    mixed-resolution inputs are handled correctly without an expensive
    intermediate reproject.
    """
    from concurrent.futures import Future, ThreadPoolExecutor

    valid_fraction = _resample_da(valid_mask.astype(np.float32), target_shape, "area", template=template)
    valid_support = valid_fraction > 0.0
    valid_denominator = valid_fraction.where(valid_support)

    combined_valid_support = valid_support
    data_vars: dict[str, xr.DataArray] = {}
    band_loader = dataset.attrs.get(_TOA_BAND_LOADER_ATTR)

    # ── Separate bands by resolution path ────────────────────────────
    in_dataset_names: list[str] = []
    lazy_names: list[str] = []
    for name in band_names:
        if name in dataset.data_vars:
            in_dataset_names.append(name)
        elif callable(band_loader):
            lazy_names.append(name)

    # ── Resample same-resolution bands in parallel ───────────────────
    def _resample_in_dataset_band(name: str) -> tuple[str, xr.DataArray]:
        masked_band = dataset[name].where(valid_mask, 0.0)
        masked_mean = _resample_da(masked_band, target_shape, "area", template=template)
        return name, masked_mean / valid_denominator

    def _resample_lazy_band(name: str) -> tuple[str, xr.DataArray, xr.DataArray] | None:
        try:
            band_da = band_loader(name, native=True)
        except (KeyError, RuntimeError):
            return None
        band_vals = band_da.values
        band_finite = np.isfinite(band_vals) & (band_vals > 0.0) & (band_vals < 1.0)
        if cloud_mask is not None:
            cm_at_band_res = _resample_da(
                cloud_mask.astype(np.float32), band_da.shape, "nearest",
            ).values > 0.5
            band_valid = band_finite & ~cm_at_band_res
        else:
            band_valid = band_finite
        band_valid_da = xr.DataArray(
            band_valid.astype(np.float32), dims=band_da.dims, coords=band_da.coords,
        )
        band_valid_frac = _resample_da(band_valid_da, target_shape, "area", template=template)
        band_valid_sup = band_valid_frac > 0.0
        band_valid_denom = band_valid_frac.where(band_valid_sup)
        masked_band = band_da.where(
            xr.DataArray(band_valid, dims=band_da.dims, coords=band_da.coords), 0.0,
        )
        masked_mean = _resample_da(masked_band, target_shape, "area", template=template)
        return name, masked_mean / band_valid_denom, band_valid_sup

    n_workers = max(len(in_dataset_names) + len(lazy_names), 1)
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        ds_futures: list[Future[tuple[str, xr.DataArray]]] = [
            pool.submit(_resample_in_dataset_band, n) for n in in_dataset_names
        ]
        lazy_futures: list[Future[tuple[str, xr.DataArray, xr.DataArray] | None]] = [
            pool.submit(_resample_lazy_band, n) for n in lazy_names
        ]

        for fut in ds_futures:
            name, resampled = fut.result()
            data_vars[name] = resampled

        for lazy_fut in lazy_futures:
            result = lazy_fut.result()
            if result is None:
                continue
            name, resampled, band_valid_sup = result
            data_vars[name] = resampled
            combined_valid_support = combined_valid_support & band_valid_sup

    if not data_vars:
        raise ValueError("No query bands were available in the corrected reflectance dataset")
    return xr.Dataset(data_vars), combined_valid_support.astype(bool)


def _resample_dataset(
    dataset: xr.Dataset,
    *,
    band_names: collections.abc.Sequence[str],
    target_shape: tuple[int, int],
    template: xr.DataArray | None = None,
) -> xr.Dataset:
    data_vars = {
        name: _resample_da(dataset[name], target_shape, "area", template=template)
        for name in band_names
        if name in dataset.data_vars
    }
    if not data_vars:
        raise ValueError("No query bands were available in the corrected reflectance dataset")
    return xr.Dataset(data_vars)


def _apply_invalid_mask_to_dataset(
    dataset: xr.Dataset,
    *,
    band_names: collections.abc.Sequence[str],
    invalid_mask: xr.DataArray,
) -> xr.Dataset:
    invalid = invalid_mask.astype(bool)
    data_vars = {
        name: dataset[name].where(~invalid)
        for name in band_names
        if name in dataset.data_vars
    }
    if not data_vars:
        raise ValueError("No query bands were available in the corrected reflectance dataset")
    return xr.Dataset(data_vars)


def _database_query_template(database: MonthlyCompositeDatabase) -> xr.DataArray:
    summary = database.median_summary
    if "feature" in summary.dims:
        return cast("xr.DataArray", summary.isel(feature=0, drop=True))
    return cast("xr.DataArray", summary)


def _resolve_provider_source_bands(
    brdf_provider: object,
    fallback: collections.abc.Sequence[SensorBand],
) -> list[SensorBand]:
    raw_source_bands = getattr(brdf_provider, "source_bands", None)
    if raw_source_bands is None:
        return _deduplicate_bands(fallback)
    if isinstance(raw_source_bands, (str, bytes, collections.abc.Mapping)):
        raise TypeError("brdf_provider.source_bands must be a sequence of SensorBand objects")
    try:
        source_bands = list(raw_source_bands)
    except TypeError as exc:
        raise TypeError("brdf_provider.source_bands must be a sequence of SensorBand objects") from exc
    if not all(isinstance(band, SensorBand) for band in source_bands):
        raise TypeError("brdf_provider.source_bands must contain only SensorBand objects")
    return _deduplicate_bands(source_bands)


def _get_temporal_brdf_parameters_batch(
    brdf_provider: object,
    *,
    bounds: tuple[float, float, float, float],
    crs: str,
    obs_times: collections.abc.Sequence[datetime],
    target_resolution: float,
    bands: collections.abc.Sequence[SensorBand],
    temporal_windows: collections.abc.Sequence[int],
    sample_date_sets: collections.abc.Sequence[collections.abc.Sequence[datetime]],
) -> collections.abc.Sequence[BRDFKernelWeights]:
    method = getattr(brdf_provider, "get_temporal_brdf_parameters_batch", None)
    if method is None or not callable(method):
        raise TypeError("brdf_provider must define get_temporal_brdf_parameters_batch(...)")
    batch_fetcher = cast("Callable[..., collections.abc.Sequence[BRDFKernelWeights]]", method)
    return batch_fetcher(
        bounds=bounds,
        crs=crs,
        obs_times=obs_times,
        target_resolution=target_resolution,
        bands=bands,
        temporal_windows=temporal_windows,
        sample_date_sets=sample_date_sets,
    )
