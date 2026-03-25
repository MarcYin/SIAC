"""Route-B monthly-database construction and runtime query helpers."""

from __future__ import annotations

import calendar
import collections.abc
import logging
from dataclasses import dataclass
from datetime import datetime
from typing import TYPE_CHECKING, Any, Protocol, cast

import numpy as np
import xarray as xr

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
    HyperspectralLibrary,
    SpectralMapper,
    SpectralMappingConfig,
)
from siac.domain import SensorBand
from siac.runtime import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    ObservationBundle,
    SurfacePrior,
)

if TYPE_CHECKING:
    from collections.abc import Callable

logger = logging.getLogger(__name__)


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
    spectral_library: HyperspectralLibrary | SpectralMappingConfig | None = None,
    spectral_k_neighbors: int = 5,
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
        if tuple(band.name for band in source_bands) != tuple(band.name for band in target_bands)
        else None
    )
    composites = [
        _normalize_monthly_composite_to_target_basis(
            composite,
            geometry=geometry,
            target_bands=target_bands,
            spectral_mapper=spectral_mapper,
        )
        for composite in monthly_composites.composites
    ]
    logger.info("Monthly surface-prior database: assembling final database from %d composite(s)", len(composites))
    database = build_monthly_composite_database(
        composites,
        query_bands=tuple(band.name for band in query_bands),
        visible_bands=tuple(band.name for band in visible_bands),
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
    if resolution <= 0:
        raise ValueError("resolution must be > 0")

    source_bands = _resolve_provider_source_bands(brdf_provider, fallback_source_bands)
    obs_time = cast("datetime", observation.metadata["observation_time"])
    month_specs = [
        _MonthSpec(
            year=year,
            month=month,
            center_time=_month_center_datetime(year, month, obs_time),
            temporal_window=max(1, calendar.monthrange(year, month)[1] // 2 + 1),
            sample_dates=_weekly_sample_dates(year, month),
        )
        for year, month in _iter_history_months(obs_time)
    ]
    logger.info(
        "Monthly composite generation from BRDF: obs_time=%s months=%d batch_size=%d source_bands=%d",
        obs_time.isoformat(),
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
            bounds=observation.bounds,
            crs=observation.crs,
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
    k_neighbors: int = 3,
) -> SurfacePrior:
    """Build a visible-band surface prior from first-pass corrected NIR/SWIR."""
    expected_query = tuple(database.query_band_names)
    if query_band_names is not None and tuple(query_band_names) != expected_query:
        raise ValueError("query_band_names must match the database query-band ordering")

    expected_visible = tuple(database.visible_band_names)
    if visible_band_names is not None and tuple(visible_band_names) != expected_visible:
        raise ValueError("visible_band_names must match the database visible-band ordering")

    target_shape = (
        int(database.median_summary.sizes["y"]),
        int(database.median_summary.sizes["x"]),
    )
    coarse_query_toa = _resample_dataset(
        observation.toa,
        band_names=expected_query,
        target_shape=target_shape,
    )
    coarse_geometry = _resample_geometry_to_target_shape(observation.geometry, target_shape)
    coarse_atmo = _resample_atmo_to_target_shape(atmo_prior, target_shape)
    coarse_cloud_mask = _resample_cloud_mask_to_target_shape(observation.cloud_mask, target_shape)
    corrector = AtmosphericCorrector(rt_model, observation.sensor_config)
    correction = corrector.correct(
        coarse_query_toa,
        coarse_geometry,
        coarse_atmo,
        cloud_mask=coarse_cloud_mask,
    )
    corrected_query = _resample_dataset(
        correction.boa,
        band_names=expected_query,
        target_shape=target_shape,
    )
    predicted_visible, predicted_unc = database.predict_visible(
        corrected_query,
        k_neighbors=k_neighbors,
    )

    cloud_mask = _resample_cloud_mask_to_target_shape(correction.cloud_mask, target_shape)
    valid = np.all(np.isfinite(predicted_visible.values), axis=0) & np.all(np.isfinite(predicted_unc.values), axis=0)
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
    mapped_reflectance, mapped_uncertainty = spectral_mapper.map(
        stacked_reflectance,
        source_uncertainty=stacked_uncertainty,
    )

    remapped: list[_MonthlyModeledReflectance] = []
    start = 0
    for monthly_input in monthly_inputs:
        stop = start + int(monthly_input.reflectance.sizes["time"])
        mapped_month = mapped_reflectance.isel(time=slice(start, stop))
        mapped_unc_month = mapped_uncertainty.isel(time=slice(start, stop))
        remapped.append(
            _MonthlyModeledReflectance(
                year=monthly_input.year,
                month=monthly_input.month,
                reflectance=mapped_month,
                quality=_combine_quality_with_mapping_uncertainty(monthly_input.quality, mapped_unc_month),
                reflectance_unc=mapped_unc_month,
            )
        )
        start = stop
    return remapped


def _combine_quality_with_mapping_uncertainty(
    quality: xr.DataArray,
    mapped_uncertainty: xr.DataArray,
) -> xr.DataArray:
    combined = np.sqrt(
        np.square(quality.values, dtype=np.float32)
        + np.square(mapped_uncertainty.mean(dim="band", skipna=True).values, dtype=np.float32)
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
) -> xr.DataArray:
    combined = np.sqrt(
        np.square(quality.values, dtype=np.float32)
        + np.square(mapped_uncertainty.mean(dim="band", skipna=True).values, dtype=np.float32)
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
    target_bands: collections.abc.Sequence[SensorBand],
    spectral_mapper: SpectralMapper | None,
) -> MonthlyBestPixelComposite:
    template = geometry.sza
    quality = _resample_spatial_field_to_template(composite.quality.astype(np.float32), template, "bilinear")
    sample_index = _resample_spatial_field_to_template(
        composite.sample_index.astype(np.float32),
        template,
        "nearest",
    ).round().astype(np.int16)

    if isinstance(composite, MonthlyKernelWeightComposite):
        weights = _resample_brdf_weights_to_template(composite.kernels, template)
        reflectance, reflectance_unc = _reflectance_from_kernel_weights(weights, geometry)
    else:
        reflectance = _resample_band_cube_to_template(composite.reflectance.astype(np.float32), template, "bilinear")
        reflectance_unc = None

    if spectral_mapper is not None:
        reflectance, mapped_unc = spectral_mapper.map(
            reflectance,
            source_uncertainty=reflectance_unc,
        )
        quality = _combine_composite_quality_with_mapping_uncertainty(quality, mapped_unc)
    else:
        reflectance = reflectance.sel(band=[band.name for band in target_bands]).astype(np.float32)

    return MonthlyBestPixelComposite(
        reflectance=reflectance.astype(np.float32),
        quality=quality.astype(np.float32),
        sample_index=sample_index.astype(np.int16),
        year=composite.year,
        month=composite.month,
    )


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
) -> BRDFKernelWeights:
    return BRDFKernelWeights(
        f0=_resample_band_cube_to_template(weights.f0, template, "bilinear"),
        f1=_resample_band_cube_to_template(weights.f1, template, "bilinear"),
        f2=_resample_band_cube_to_template(weights.f2, template, "bilinear"),
        f0_unc=_resample_band_cube_to_template(weights.f0_unc, template, "bilinear"),
        f1_unc=_resample_band_cube_to_template(weights.f1_unc, template, "bilinear"),
        f2_unc=_resample_band_cube_to_template(weights.f2_unc, template, "bilinear"),
        reflectance_unc=(
            _resample_band_cube_to_template(weights.reflectance_unc, template, "bilinear")
            if weights.reflectance_unc is not None
            else None
        ),
    )


def _resample_band_cube_to_template(
    data: xr.DataArray,
    template: xr.DataArray,
    method: str,
) -> xr.DataArray:
    if _shares_spatial_grid(data, template):
        return data.astype(np.float32)

    target_shape = (int(template.sizes["y"]), int(template.sizes["x"]))
    band_coords = data.coords["band"].values if "band" in data.coords else np.arange(data.sizes["band"])
    resampled = xr.concat(
        [
            _resample_da(data.sel(band=band, drop=True), target_shape, method)
            for band in band_coords
        ],
        dim=xr.IndexVariable("band", band_coords),
    )
    coords: dict[str, object] = {"band": band_coords}
    if "y" in template.coords:
        coords["y"] = template.coords["y"]
    if "x" in template.coords:
        coords["x"] = template.coords["x"]
    return resampled.assign_coords(**coords).astype(np.float32)


def _resample_spatial_field_to_template(
    data: xr.DataArray,
    template: xr.DataArray,
    method: str,
) -> xr.DataArray:
    if _shares_spatial_grid(data, template):
        return data.astype(np.float32)
    target_shape = (int(template.sizes["y"]), int(template.sizes["x"]))
    resampled = _resample_da(data, target_shape, method)
    coords: dict[str, object] = {}
    if "y" in template.coords:
        coords["y"] = template.coords["y"]
    if "x" in template.coords:
        coords["x"] = template.coords["x"]
    return resampled.assign_coords(**coords).astype(np.float32)


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
    if geometry.sza.shape == target_shape:
        return geometry
    return GeometryAngles(
        sza=_resample_da(geometry.sza, target_shape, "bilinear"),
        saa=_resample_da(geometry.saa, target_shape, "bilinear"),
        vza=_resample_da(geometry.vza, target_shape, "bilinear"),
        vaa=_resample_da(geometry.vaa, target_shape, "bilinear"),
    )


def _resample_atmo_to_target_shape(
    atmo_prior: AtmosphericState,
    target_shape: tuple[int, int],
) -> AtmosphericState:
    if atmo_prior.aot.shape == target_shape:
        return atmo_prior
    return AtmosphericState(
        aot=_resample_da(atmo_prior.aot, target_shape, "bilinear"),
        tcwv=_resample_da(atmo_prior.tcwv, target_shape, "bilinear"),
        tco3=_resample_da(atmo_prior.tco3, target_shape, "bilinear"),
        aot_unc=_resample_da(atmo_prior.aot_unc, target_shape, "bilinear"),
        tcwv_unc=_resample_da(atmo_prior.tcwv_unc, target_shape, "bilinear"),
        tco3_unc=_resample_da(atmo_prior.tco3_unc, target_shape, "bilinear"),
        elevation=_resample_da(atmo_prior.elevation, target_shape, "bilinear"),
    )


def _resample_cloud_mask_to_target_shape(
    cloud_mask: xr.DataArray,
    target_shape: tuple[int, int],
) -> xr.DataArray:
    if "band" in cloud_mask.dims:
        cloud_mask = cloud_mask.any(dim="band")
    if cloud_mask.shape == target_shape:
        return cloud_mask
    return _resample_cloud_mask(cloud_mask, target_shape)


def _resample_dataset(
    dataset: xr.Dataset,
    *,
    band_names: collections.abc.Sequence[str],
    target_shape: tuple[int, int],
) -> xr.Dataset:
    data_vars = {
        name: _resample_da(dataset[name], target_shape, "area")
        for name in band_names
        if name in dataset.data_vars
    }
    if not data_vars:
        raise ValueError("No query bands were available in the corrected reflectance dataset")
    return xr.Dataset(data_vars)


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
