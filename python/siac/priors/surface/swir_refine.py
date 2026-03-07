"""Route-B monthly-database construction and runtime query helpers."""

from __future__ import annotations

import calendar
from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from siac.core.types import (
    AtmosphericState,
    BRDFKernelWeights,
    GeometryAngles,
    ObservationBundle,
    SensorBand,
    SurfacePrior,
)
from siac.correction.atmospheric import AtmosphericCorrector
from siac.grid.assembler import _compute_target_shape, _resample_cloud_mask, _resample_da
from siac.priors.brdf.kernels import BRDFKernels, compute_reflectance
from siac.priors.surface.brdf_monthly_composite import build_monthly_best_pixel_composite
from siac.priors.surface.brdf_monthly_database import (
    MonthlyCompositeDatabase,
    build_monthly_composite_database,
)
from siac.priors.surface.spectral_mapping import map_multispectral_reflectance

if TYPE_CHECKING:
    from collections.abc import Sequence


_HISTORY_YEARS = 5
_HISTORY_MONTH_OFFSETS = (-1, 0, 1)
_WEEKLY_STEP_DAYS = 7


def build_monthly_surface_prior_database(
    *,
    observation: ObservationBundle,
    brdf_provider: object,
    resolution: float,
    geometry: GeometryAngles,
    visible_bands: Sequence[SensorBand],
    query_bands: Sequence[SensorBand],
) -> MonthlyCompositeDatabase:
    """Build the 15-month Route-B database for the current scene."""
    if resolution <= 0:
        raise ValueError("resolution must be > 0")

    required_bands = _deduplicate_bands([*visible_bands, *query_bands])
    source_bands = _deduplicate_bands(getattr(brdf_provider, "source_bands", required_bands))
    obs_time = observation.metadata["observation_time"]
    month_specs = [
        (
            year,
            month,
            _month_center_datetime(year, month, obs_time),
            max(1, calendar.monthrange(year, month)[1] // 2 + 1),
            _weekly_sample_dates(year, month),
        )
        for year, month in _iter_history_months(obs_time)
    ]
    temporal_weights_list = brdf_provider.get_temporal_brdf_parameters_batch(
        bounds=observation.bounds,
        crs=observation.crs,
        obs_times=[spec[2] for spec in month_specs],
        target_resolution=resolution,
        bands=source_bands,
        temporal_windows=[spec[3] for spec in month_specs],
        sample_date_sets=[spec[4] for spec in month_specs],
    )
    composites = []
    for (year, month, _center_time, _temporal_window, _sample_dates), temporal_weights in zip(
        month_specs,
        temporal_weights_list,
        strict=True,
    ):
        reflectance, quality, reflectance_unc = _forward_model_monthly_reflectance(
            temporal_weights,
            geometry=geometry,
            year=year,
            month=month,
        )
        if tuple(band.name for band in source_bands) != tuple(band.name for band in required_bands):
            reflectance, mapped_unc = map_multispectral_reflectance(
                reflectance,
                source_bands=source_bands,
                target_bands=required_bands,
                source_uncertainty=reflectance_unc,
            )
            quality = np.sqrt(
                np.square(quality.values, dtype=np.float32)
                + np.square(mapped_unc.mean(dim="band", skipna=True).values, dtype=np.float32)
            ).astype(np.float32)
            quality = xr.DataArray(
                quality,
                dims=["time", "y", "x"],
                coords={
                    "time": reflectance.coords["time"],
                    "y": reflectance.coords["y"],
                    "x": reflectance.coords["x"],
                },
            )
        composites.append(
            build_monthly_best_pixel_composite(
                reflectance,
                quality,
                year=year,
                month=month,
            )
        )

    return build_monthly_composite_database(
        composites,
        query_bands=tuple(band.name for band in query_bands),
        visible_bands=tuple(band.name for band in visible_bands),
    )


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

    aligned_atmo = _resample_atmo_to_observation_grid(observation, atmo_prior)
    corrector = AtmosphericCorrector(rt_model, observation.sensor_config)
    correction = corrector.correct(
        observation.toa,
        observation.geometry,
        aligned_atmo,
        cloud_mask=observation.cloud_mask,
    )

    target_shape = (database.median_summary.sizes["y"], database.median_summary.sizes["x"])
    corrected_query = _resample_dataset(
        correction.boa,
        band_names=expected_query,
        target_shape=target_shape,
    )
    predicted_visible, predicted_unc = database.predict_visible(
        corrected_query,
        k_neighbors=k_neighbors,
    )

    cloud_mask = correction.cloud_mask
    if "band" in cloud_mask.dims:
        cloud_mask = cloud_mask.any(dim="band")
    if cloud_mask.shape != target_shape:
        cloud_mask = _resample_cloud_mask(cloud_mask, target_shape)
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
    native_resolution = _native_observation_resolution(observation)
    first_var = next(iter(observation.toa.data_vars))
    target_shape = _compute_target_shape(observation.toa[first_var].shape, native_resolution, resolution)
    return GeometryAngles(
        sza=_resample_da(observation.geometry.sza, target_shape, "bilinear"),
        saa=_resample_da(observation.geometry.saa, target_shape, "bilinear"),
        vza=_resample_da(observation.geometry.vza, target_shape, "bilinear"),
        vaa=_resample_da(observation.geometry.vaa, target_shape, "bilinear"),
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
    kernels = BRDFKernels(hb=2.0, br=1.0)
    k_vol, k_geo = kernels.compute(geometry.vza, geometry.sza, geometry.raa)
    reflectance = compute_reflectance(
        temporal_weights.f0,
        temporal_weights.f1,
        temporal_weights.f2,
        k_vol,
        k_geo,
    ).transpose("time", "band", "y", "x")
    reflectance_unc = temporal_weights.compute_reflectance_uncertainty(k_vol, k_geo).transpose("time", "band", "y", "x")

    month_mask = _select_month_mask(temporal_weights.f0.coords["time"].values, year=year, month=month)
    if month_mask.any():
        reflectance = reflectance.isel(time=month_mask)
        reflectance_unc = reflectance_unc.isel(time=month_mask)

    quality = reflectance_unc.mean(dim="band", skipna=True).astype(np.float32)
    return reflectance.astype(np.float32), quality, reflectance_unc.astype(np.float32)


def _select_month_mask(time_values: np.ndarray, *, year: int, month: int) -> np.ndarray:
    month_strings = np.asarray(time_values, dtype="datetime64[D]").astype("datetime64[M]")
    target = np.datetime64(f"{year:04d}-{month:02d}", "M")
    return month_strings == target


def _deduplicate_bands(bands: Sequence[SensorBand]) -> list[SensorBand]:
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
    return tuple(int(size) for size in observation.toa[first_var].shape)


def _resample_atmo_to_observation_grid(
    observation: ObservationBundle,
    atmo_prior: AtmosphericState,
) -> AtmosphericState:
    target_shape = _observation_shape(observation)
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


def _resample_dataset(
    dataset: xr.Dataset,
    *,
    band_names: Sequence[str],
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
