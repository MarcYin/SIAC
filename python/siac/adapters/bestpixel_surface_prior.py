"""Bestpixel-backed surface prior for the surface-driven aerosol solver.

This builds an opt-in surface-reflectance prior (``algorithms.surface_prior
.method = "bestpixel"``) directly from the external ``bestpixel`` package's
cloud-free monthly best-pixel composites, feeding the
:class:`~siac.algorithms.solver.surface_driven.SurfaceDrivenSolver`.

The validated recipe (from the experimental ``diag_seasonal_retrieval.py``
``comp_ref`` harness): build one composite per LOYO ``(year, month)`` window,
spectrally map each onto the scene's sensor bands, then reduce *across* the N
realizations per band with a robust temporal statistic —

    boa      = nanmedian(realizations after a robust MAD sigma-clip)
    boa_unc  = sqrt(nanmean((x - median)**2))   (RMSE about the median), floored

The temporal RMSE-across-realizations is precisely the per-band ``boa_unc`` the
surface-driven solver wants (self-calibrating temporal spread), matching the
``comp_ref`` harness. When only one realization survives, the per-composite
bestpixel uncertainty is used instead.

A per-acquisition-day MAIAC AOD gate (see
:mod:`siac.adapters.atmo.maiac_day_aod`) optionally drops the high-aerosol days
per window before compositing so the prior is built from the cleanest scenes.

``bestpixel`` and ``earthaccess`` are optional dependencies imported lazily;
importing this module never requires them.
"""

from __future__ import annotations

import logging
import math
from collections import defaultdict
from datetime import date, datetime, timedelta
from typing import TYPE_CHECKING, Any, Protocol, cast

import numpy as np
import xarray as xr

from siac.adapters.bestpixel import (
    _BESTPIXEL_DN_SCALE,
    _QUALITY_NODATA,
    _coords_from_grid,
    _make_source_da,
)
from siac.algorithms.grid.assembler import _build_target_template
from siac.algorithms.surface.spectral_mapping import SpectralMapper
from siac.geo._spatial import copy_spatial_metadata_like
from siac.geo.reprojection import reproject_match, transform_bounds
from siac.runtime.models import SurfacePrior

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.domain import SensorBand
    from siac.runtime import ObservationBundle

logger = logging.getLogger(__name__)

#: BOA-uncertainty floor (matches the surface-driven solver's ``_MIN_BOA_UNC``).
_BOA_UNC_FLOOR = 0.006
#: MAD -> robust standard-deviation scale for a normal distribution.
_MAD_TO_STD = 1.4826


class MAIACDayAODCallable(Protocol):
    """Signature of the injectable per-day MAIAC AOD gate source."""

    def __call__(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        periods: Sequence[tuple[int, int]],
    ) -> dict[str, float]: ...


def build_aod_by_day_gate(
    raw_day_aod: dict[str, float],
    *,
    aod_max: float | None,
    low_aod_frac: float,
) -> tuple[dict[str, float] | None, float | None]:
    """Reduce a raw per-day AOD map to the days bestpixel should keep.

    Returns ``(gate_map, effective_aod_max)``. ``gate_map`` is ``None`` when no
    gate should be applied (the raw map was empty or nothing survived), in which
    case the caller keeps every scene day. Otherwise it contains only the kept
    days; passed with ``reject_unknown=True`` it both drops the high-AOD days and
    rejects days bestpixel sees that the gate never scored.

    Per ``(year, month)`` window: keep all days with ``AOD <= aod_max`` when
    ``aod_max`` is set, else keep the lowest ``low_aod_frac`` fraction of days.
    """
    if not raw_day_aod:
        return None, None
    by_window: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for day, aod in raw_day_aod.items():
        by_window[day[:7]].append((day, float(aod)))

    kept: dict[str, float] = {}
    for items in by_window.values():
        ordered = sorted(items, key=lambda kv: kv[1])
        if aod_max is not None:
            chosen = [(day, aod) for day, aod in ordered if aod <= float(aod_max)]
        else:
            n_keep = max(1, math.ceil(float(low_aod_frac) * len(ordered)))
            chosen = ordered[:n_keep]
        for day, aod in chosen:
            kept[day] = aod

    if not kept:
        return None, None
    effective_aod_max = float(aod_max) if aod_max is not None else max(kept.values()) + 1.0
    return kept, effective_aod_max


def reproject_period_to_template(
    period: dict[str, Any],
    target_template: xr.DataArray,
    bands: Sequence[str],
) -> tuple[xr.DataArray, xr.DataArray] | None:
    """Reproject one bestpixel period's reflectance + uncertainty to the grid.

    Returns ``(reflectance, boa_unc)`` as ``(band, y, x)`` DataArrays in
    reflectance units on ``target_template``'s grid (NaN at nodata), or ``None``
    when the period carries no valid pixels over the AOI.
    """
    grid = period["grid"]
    source_crs = str(grid.get("crs") or f"EPSG:{int(grid['epsg'])}")
    x_coords, y_coords = _coords_from_grid(grid)

    bands_dict = period["bands"]
    unc_dict = period.get("boa_unc") or {}
    reflectance_scale = float(period.get("reflectance_scale", _BESTPIXEL_DN_SCALE))
    quality = np.asarray(period["quality"])
    height, width = quality.shape
    nodata_mask = quality == _QUALITY_NODATA

    reflectance = np.empty((len(bands), height, width), dtype=np.float32)
    uncertainty = np.full((len(bands), height, width), np.nan, dtype=np.float32)
    for index, name in enumerate(bands):
        if name not in bands_dict:
            raise KeyError(
                f"bestpixel composite is missing requested band {name!r}; "
                f"returned bands: {sorted(bands_dict)}"
            )
        reflectance[index] = np.asarray(bands_dict[name], dtype=np.float32) * reflectance_scale
        if name in unc_dict:
            uncertainty[index] = np.asarray(unc_dict[name], dtype=np.float32) * reflectance_scale
    reflectance[:, nodata_mask] = np.nan
    uncertainty[:, nodata_mask] = np.nan

    if not np.isfinite(reflectance).any():
        return None

    reflectance_da = _make_source_da(reflectance, x_coords, y_coords, source_crs, band_names=bands)
    uncertainty_da = _make_source_da(uncertainty, x_coords, y_coords, source_crs, band_names=bands)
    reflectance_aligned = reproject_match(
        reflectance_da, target_template, resampling="bilinear", nodata=np.nan
    ).astype(np.float32)
    uncertainty_aligned = reproject_match(
        uncertainty_da, target_template, resampling="bilinear", nodata=np.nan
    ).astype(np.float32)
    return reflectance_aligned, uncertainty_aligned


def _native_grid_template(period: dict[str, Any]) -> xr.DataArray:
    """Build a 2-D georeferenced template on a bestpixel period's NATIVE grid.

    Used by the surface-driven ``resolve_on_prior_grid`` path: the prior (and
    every realization) is built on the composite's own tile-aligned grid instead
    of a fresh observation-bounds grid, so the solver resolves on the native grid
    and the prior is never smeared by a sub-pixel reprojection onto obs bounds.
    """
    grid = period["grid"]
    source_crs = str(grid.get("crs") or f"EPSG:{int(grid['epsg'])}")
    x_coords, y_coords = _coords_from_grid(grid)
    template = xr.DataArray(
        np.full((y_coords.size, x_coords.size), np.nan, dtype=np.float32),
        dims=("y", "x"),
        coords={"y": y_coords, "x": x_coords},
    )
    template = template.rio.set_spatial_dims(x_dim="x", y_dim="y")
    return template.rio.write_crs(source_crs)


def reduce_realizations(
    reflectance_stack: np.ndarray,
    uncertainty_first: np.ndarray | None,
    *,
    robust_clip: float = 0.0,
    floor: float = _BOA_UNC_FLOOR,
) -> tuple[np.ndarray, np.ndarray]:
    """Reduce ``(N, band, y, x)`` realizations to ``(band, y, x)`` boa + boa_unc.

    ``boa`` is the per-band temporal nanmedian; ``boa_unc`` is the temporal RMSE
    about that median (``sqrt(nanmean((x - median)**2))``), floored — matching the
    ``comp_ref`` harness (``diag_seasonal_retrieval.py``). With a single
    realization the per-composite ``uncertainty_first`` is used. An optional
    ``robust_clip`` drops realizations more than ``robust_clip * MAD`` from the
    median before the final reduction (so bad-AC outlier days don't inflate the
    RMSE).
    """
    stack = np.asarray(reflectance_stack, dtype=np.float32)
    if stack.ndim != 4:
        raise ValueError(f"reflectance_stack must be (N, band, y, x), got shape {stack.shape}")
    n_realizations = stack.shape[0]

    with np.errstate(invalid="ignore"):
        median = np.nanmedian(stack, axis=0)
        if n_realizations == 1:
            if uncertainty_first is not None:
                boa_unc = np.asarray(uncertainty_first, dtype=np.float32)
            else:
                boa_unc = np.full_like(median, floor)
        else:
            if robust_clip > 0.0:
                mad0 = np.nanmedian(np.abs(stack - median[np.newaxis]), axis=0) * _MAD_TO_STD
                deviation = np.abs(stack - median[np.newaxis])
                keep = (deviation <= (float(robust_clip) * mad0[np.newaxis])) | (
                    mad0[np.newaxis] == 0.0
                )
                stack = np.where(keep, stack, np.nan)
                median = np.nanmedian(stack, axis=0)
            # RMSE about the (clipped) median — matches the comp_ref harness
            # (diag_seasonal_retrieval.py:959). The robust_clip above already
            # removed outlier days, so the mean-square spread is not inflated.
            boa_unc = np.sqrt(np.nanmean((stack - median[np.newaxis]) ** 2, axis=0))

    boa = np.asarray(median, dtype=np.float32)
    boa_unc = np.maximum(np.asarray(boa_unc, dtype=np.float32), float(floor))
    boa_unc = np.where(np.isfinite(boa), boa_unc, np.nan).astype(np.float32)
    return boa, boa_unc


def _resolve_years(obs_time: datetime, lookback_years: int) -> tuple[int, ...]:
    if lookback_years <= 0:
        raise ValueError(f"bestpixel_lookback_years must be > 0, got {lookback_years}")
    end_year = obs_time.year - 1
    start_year = end_year - int(lookback_years) + 1
    return tuple(range(start_year, end_year + 1))


def _resolve_months(
    obs_time: datetime, months: Sequence[int] | None, seasonal_window: int = 0
) -> tuple[int, ...]:
    if months:
        return tuple(int(m) for m in months)
    center = int(obs_time.month)
    window = max(0, int(seasonal_window))
    if window == 0:
        return (center,)
    # Scene month +/- window, wrapping into 1..12 (May +/-1 -> Apr,May,Jun;
    # Jan +/-1 -> Dec,Jan,Feb), matching the comp_ref harness seasonal window.
    return tuple(sorted({((center - 1 + off) % 12) + 1 for off in range(-window, window + 1)}))


def _fetch_periods(
    config: Any,
    observation: ObservationBundle,
    resolution: float,
    bands: Sequence[str],
    *,
    maiac_day_aod: MAIACDayAODCallable | None,
) -> list[dict[str, Any]]:
    """Run the gated ``bestpixel.build_monthly_composites`` call once."""
    obs_time = observation.metadata.get("observation_time")
    if not isinstance(obs_time, datetime):
        raise ValueError(
            "bestpixel surface prior requires observation.metadata['observation_time'] "
            "to be a datetime."
        )
    mc = config.providers.monthly_composites
    reduction = str(
        getattr(config.algorithms.surface_prior, "bestpixel_window_reduction", "window")
    )
    if reduction == "daily_median":
        return _fetch_daily_median_periods(
            config,
            observation,
            resolution,
            bands,
            maiac_day_aod=maiac_day_aod,
        )
    if reduction != "window":
        raise ValueError(
            "surface_prior.bestpixel_window_reduction must be 'window' or 'daily_median', "
            f"got {reduction!r}"
        )

    years = _resolve_years(obs_time, int(mc.bestpixel_lookback_years))
    months = _resolve_months(
        obs_time, mc.bestpixel_months, int(getattr(mc, "bestpixel_seasonal_window_months", 0))
    )
    periods = [(year, month) for year in years for month in months]

    gate_map, effective_aod_max = _build_gate(
        config, observation, periods, maiac_day_aod=maiac_day_aod
    )

    west, south, east, north = transform_bounds(observation.bounds, observation.crs, "EPSG:4326")
    bbox = [float(west), float(south), float(east), float(north)]
    disk_cache = getattr(mc, "bestpixel_disk_cache", None)

    import bestpixel as bp  # lazy optional dependency

    kwargs: dict[str, Any] = {
        "bbox": bbox,
        "years": list(years),
        "months": list(months),
        "resolution": float(resolution),
        "top_k": int(mc.bestpixel_top_k),
        "endpoint": str(mc.bestpixel_endpoint),
        "bands": list(bands),
        "disk_cache": str(disk_cache) if disk_cache is not None else None,
        "max_cloud_cover": float(mc.bestpixel_max_cloud_cover),
        "output_crs": "utm",
        "emit_uncertainty": True,
    }
    if gate_map:
        kwargs["aod_by_day"] = gate_map
        kwargs["aod_max"] = effective_aod_max
        kwargs["reject_unknown"] = True
    logger.info(
        "bestpixel surface prior: endpoint=%s years=%s months=%s bands=%s gated_days=%d",
        mc.bestpixel_endpoint,
        list(years),
        list(months),
        list(bands),
        len(gate_map) if gate_map else 0,
    )
    return cast("list[dict[str, Any]]", list(bp.build_monthly_composites(**kwargs)))


def _fetch_daily_median_periods(
    config: Any,
    observation: ObservationBundle,
    resolution: float,
    bands: Sequence[str],
    *,
    maiac_day_aod: MAIACDayAODCallable | None,
) -> list[dict[str, Any]]:
    """Build harness-style period composites from clean daily top-1 surfaces.

    For each (year, month) window, select the low-AOD days first, fetch each day
    independently with ``top_k=1``, reproject onto the scene template, then median
    those daily reflectance surfaces. This mirrors ``build_l2a_seasonal.py`` more
    closely than bestpixel's single per-window composite.
    """
    obs_time = observation.metadata.get("observation_time")
    if not isinstance(obs_time, datetime):
        raise ValueError(
            "bestpixel surface prior requires observation.metadata['observation_time'] "
            "to be a datetime."
        )
    mc = config.providers.monthly_composites
    years = _resolve_years(obs_time, int(mc.bestpixel_lookback_years))
    months = _resolve_months(
        obs_time, mc.bestpixel_months, int(getattr(mc, "bestpixel_seasonal_window_months", 0))
    )
    periods = [(year, month) for year in years for month in months]
    gate_map, _effective_aod_max = _build_gate(
        config, observation, periods, maiac_day_aod=maiac_day_aod
    )
    if not gate_map:
        raise RuntimeError(
            "bestpixel_window_reduction='daily_median' requires a MAIAC day-AOD gate; "
            "no clean days were available."
        )

    west, south, east, north = transform_bounds(observation.bounds, observation.crs, "EPSG:4326")
    bbox = [float(west), float(south), float(east), float(north)]
    disk_cache = getattr(mc, "bestpixel_disk_cache", None)
    target_template = _build_target_template(
        observation.bounds, str(observation.crs), float(resolution)
    )

    import bestpixel as bp  # lazy optional dependency

    out_periods: list[dict[str, Any]] = []
    for year, month in periods:
        month_key = f"{year:04d}-{month:02d}"
        days = sorted(day for day in gate_map if day.startswith(month_key))
        if not days:
            continue
        daily_reflectance: list[np.ndarray] = []
        daily_uncertainty: list[np.ndarray] = []
        source_ids: list[str] = []
        for day in days:
            start = date.fromisoformat(day)
            stop = start + timedelta(days=1)
            try:
                period = bp.build_composite(
                    bbox=bbox,
                    datetime=f"{start.isoformat()}/{stop.isoformat()}",
                    resolution=float(resolution),
                    top_k=1,
                    endpoint=str(mc.bestpixel_endpoint),
                    bands=list(bands),
                    disk_cache=str(disk_cache) if disk_cache is not None else None,
                    max_cloud_cover=float(mc.bestpixel_max_cloud_cover),
                    output_crs="utm",
                    emit_uncertainty=True,
                )
            except Exception as exc:  # noqa: BLE001 - one bad day should not kill the window
                logger.info(
                    "bestpixel daily_median: %s failed (%s: %s); skipping",
                    day,
                    type(exc).__name__,
                    exc,
                )
                continue
            reprojected = reproject_period_to_template(period, target_template, bands)
            if reprojected is None:
                continue
            reflectance_da, uncertainty_da = reprojected
            reflectance = np.asarray(reflectance_da.values, dtype=np.float32)
            if not _has_enough_valid_surface(reflectance, bands):
                continue
            daily_reflectance.append(reflectance)
            daily_uncertainty.append(np.asarray(uncertainty_da.values, dtype=np.float32))
            source_ids.extend(str(src) for src in period.get("source_ids", ()))

        if not daily_reflectance:
            logger.info("bestpixel daily_median: %s had no usable clean days", month_key)
            continue
        stack = np.stack(daily_reflectance, axis=0)
        with np.errstate(invalid="ignore"):
            median = np.nanmedian(stack, axis=0).astype(np.float32)
            if stack.shape[0] == 1 and daily_uncertainty:
                uncertainty = daily_uncertainty[0].astype(np.float32)
            else:
                uncertainty = np.sqrt(np.nanmean((stack - median[np.newaxis]) ** 2, axis=0))
                uncertainty = np.maximum(uncertainty, _BOA_UNC_FLOOR).astype(np.float32)
        out_periods.append(
            _period_from_reflectance_template(
                target_template,
                bands,
                median,
                uncertainty,
                source_ids=source_ids,
                year=year,
                month=month,
            )
        )
        logger.info(
            "bestpixel daily_median: %s kept %d/%d clean day(s)",
            month_key,
            stack.shape[0],
            len(days),
        )
    logger.info(
        "bestpixel surface prior: endpoint=%s daily_median periods=%d/%d bands=%s gated_days=%d",
        mc.bestpixel_endpoint,
        len(out_periods),
        len(periods),
        list(bands),
        len(gate_map),
    )
    return out_periods


def _has_enough_valid_surface(reflectance: np.ndarray, bands: Sequence[str]) -> bool:
    if reflectance.size == 0:
        return False
    try:
        blue_index = list(bands).index("blue")
    except ValueError:
        blue_index = 0
    return bool(np.isfinite(reflectance[blue_index]).mean() > 0.2)


def _period_from_reflectance_template(
    template: xr.DataArray,
    bands: Sequence[str],
    reflectance: np.ndarray,
    uncertainty: np.ndarray,
    *,
    source_ids: Sequence[str],
    year: int,
    month: int,
) -> dict[str, Any]:
    x = np.asarray(template.coords["x"].values, dtype=float)
    y = np.asarray(template.coords["y"].values, dtype=float)
    x_res = float(abs(x[1] - x[0])) if x.size > 1 else 1.0
    y_res = float(abs(y[0] - y[1])) if y.size > 1 else x_res
    x0 = float(x[0] - 0.5 * x_res)
    y1 = float(y[0] + 0.5 * y_res)
    crs = str(template.rio.crs)
    epsg = template.rio.crs.to_epsg() if template.rio.crs is not None else None
    quality = np.where(
        np.all(np.isfinite(reflectance), axis=0),
        0,
        _QUALITY_NODATA,
    ).astype(np.uint16)
    return {
        "bands": {
            name: np.asarray(reflectance[index], dtype=np.float32)
            for index, name in enumerate(bands)
        },
        "boa_unc": {
            name: np.asarray(uncertainty[index], dtype=np.float32)
            for index, name in enumerate(bands)
        },
        "quality": quality,
        "reflectance_scale": 1.0,
        "grid": {
            "bounds": [
                x0,
                float(y[-1] - 0.5 * y_res),
                float(x[-1] + 0.5 * x_res),
                y1,
            ],
            "epsg": int(epsg) if epsg is not None else 0,
            "crs": crs,
            "resolution": float(max(x_res, y_res)),
            "width": int(x.size),
            "height": int(y.size),
            "transform": [x_res, 0.0, x0, 0.0, -y_res, y1],
        },
        "band_names": list(bands),
        "source_ids": list(source_ids),
        "year": int(year),
        "month": int(month),
    }


def _template_epsg_transform(template: xr.DataArray) -> tuple[int, np.ndarray]:
    x = np.asarray(template.coords["x"].values, dtype=float)
    y = np.asarray(template.coords["y"].values, dtype=float)
    x_res = float(abs(x[1] - x[0])) if x.size > 1 else 1.0
    y_res = float(abs(y[0] - y[1])) if y.size > 1 else x_res
    x0 = float(x[0] - 0.5 * x_res)
    y1 = float(y[0] + 0.5 * y_res)
    crs = template.rio.crs
    epsg = crs.to_epsg() if crs is not None else None
    if epsg is None:
        raise ValueError("bestpixel visible prediction requires an EPSG-resolvable prior grid CRS.")
    return int(epsg), np.asarray([x_res, 0.0, x0, 0.0, -y_res, y1], dtype=float)


def _build_gate(
    config: Any,
    observation: ObservationBundle,
    periods: Sequence[tuple[int, int]],
    *,
    maiac_day_aod: MAIACDayAODCallable | None,
) -> tuple[dict[str, float] | None, float | None]:
    sp = config.algorithms.surface_prior
    day_source = maiac_day_aod if maiac_day_aod is not None else _default_maiac_day_aod()
    try:
        raw_map = day_source(observation.bounds, str(observation.crs), list(periods))
    except Exception:  # noqa: BLE001 - a failed gate must not abort the prior
        logger.warning(
            "bestpixel surface prior: MAIAC day-AOD gate unavailable; keeping all scene days",
            exc_info=True,
        )
        raw_map = {}
    return build_aod_by_day_gate(
        raw_map,
        aod_max=sp.bestpixel_aod_max,
        low_aod_frac=float(sp.bestpixel_low_aod_frac),
    )


def _default_maiac_day_aod() -> MAIACDayAODCallable:
    """Lazily build the default earthaccess-backed MAIAC per-day AOD source."""
    from siac.adapters.atmo.maiac_day_aod import MAIACDayAODProvider

    return MAIACDayAODProvider().day_aod_map


def direct_source_indices(
    source_bands: Sequence[SensorBand], target_bands: Sequence[SensorBand]
) -> list[int] | None:
    """Per-target source-band index when every target has a direct counterpart.

    The bestpixel L2A canonical bands ARE Sentinel-2 bands (coastal→B01,
    blue→B02, red→B04, …), so for an S2 scene the bestpixel source basis already
    contains the solve bands and a KNN spectral remap is spurious — it biases the
    prior's blue/red ratio (dark-blue). Bands are matched by their canonical
    RSRF band id (the bestpixel source carries ``rsrf_band_id="B01"`` for the
    ``coastal`` source band, which equals the S2 target band's ``name``); a
    generous wavelength guard rejects accidental id collisions. When every target
    band has such a counterpart this returns the source-index list for an
    identity (direct band-column) passthrough — matching the harness ``comp_ref``
    direct-band reference. Returns ``None`` when any target lacks a direct
    counterpart (genuine cross-sensor case, e.g. mcd43a4 MODIS NBAR → S2), so the
    caller falls back to the KNN mapper.
    """

    def _identity_key(band: SensorBand) -> str:
        return str(getattr(band, "rsrf_band_id", None) or band.name)

    by_key: dict[str, int] = {}
    for source_index, band in enumerate(source_bands):
        by_key.setdefault(_identity_key(band), source_index)
    indices: list[int] = []
    for target in target_bands:
        match = by_key.get(_identity_key(target))
        if match is None:
            return None
        # Wavelength guard: identity only when the band centres are close
        # (nominal bestpixel vs catalog centres differ by at most a few nm).
        if (
            abs(float(source_bands[match].center_wavelength) - float(target.center_wavelength))
            > 30.0
        ):
            return None
        indices.append(match)
    return indices


def _resolve_surface_library(config: Any) -> Any | None:
    """Return the configured prepared surface library, or ``None`` for live build.

    A prepared library moves acquisition and atmospheric correction offline, so
    the library can be corrected in the same RT model the solver uses (which the
    live composite path cannot guarantee) and the slowest pipeline stage is paid
    once rather than per scene.
    """

    surface_cfg = getattr(getattr(config, "algorithms", None), "surface_prior", None)
    library_path = getattr(surface_cfg, "prepared_library_path", None)
    if not library_path:
        return None

    from siac.adapters.surface_library import PreparedSurfaceLibrary

    return PreparedSurfaceLibrary(
        library_path,
        band_names=getattr(surface_cfg, "prepared_library_bands", None),
        scene_key=getattr(surface_cfg, "prepared_library_scene_key", None),
    )


def build_bestpixel_surface_prior(
    config: Any,
    observation: ObservationBundle,
    resolution: float,
    *,
    bands: Sequence[str],
    source_bands: Sequence[SensorBand],
    target_bands: Sequence[SensorBand],
    spectral_library: Any | None,
    k_neighbors: int,
    maiac_day_aod: MAIACDayAODCallable | None = None,
    atmo_prior: Any | None = None,
    anchor_atmo_prior: Any | None = None,
    rt_model: Any | None = None,
) -> SurfacePrior:
    """Build the bestpixel-backed :class:`SurfacePrior` for the solver."""
    library = _resolve_surface_library(config)
    if library is None:
        periods = _fetch_periods(
            config, observation, resolution, bands, maiac_day_aod=maiac_day_aod
        )
        library_rt_space = None
    else:
        from siac.adapters.surface_library import realization_to_period

        periods = [
            realization_to_period(realization)
            for realization in library.realizations(observation, resolution, bands)
        ]
        library_rt_space = library.rt_space
    # Surface-driven "resolve on the prior's native grid": build the prior (and
    # every realization) on the composite's own tile-aligned grid instead of a
    # fresh observation-bounds grid. This removes the sub-pixel composite-vs-obs
    # resampling smear; M4 then adopts ``SurfacePrior.solver_grid`` as the solver
    # target. Default (flag off) keeps the legacy observation-bounds grid.
    solver_cfg = getattr(getattr(config, "algorithms", None), "solver", None)
    resolve_on_prior_grid = bool(getattr(solver_cfg, "surface_driven_resolve_on_prior_grid", False))
    used_native_grid = resolve_on_prior_grid and bool(periods)
    if used_native_grid:
        target_template = _native_grid_template(periods[0])
        logger.info(
            "bestpixel surface prior: resolving on the composite NATIVE grid "
            "(%dx%d, %s) — realizations co-registered onto it, no obs-bounds smear",
            int(target_template.sizes["x"]),
            int(target_template.sizes["y"]),
            target_template.rio.crs,
        )
    else:
        target_template = _build_target_template(
            observation.bounds, str(observation.crs), float(resolution)
        )
    # Identity (direct band-column) passthrough when the bestpixel source basis
    # already contains every solve band (S2 L2A → S2); only build the KNN mapper
    # for a genuine cross-sensor source basis (e.g. mcd43a4 MODIS NBAR).
    passthrough_indices = direct_source_indices(source_bands, target_bands)
    target_band_names = [str(band.name) for band in target_bands]
    mapper: SpectralMapper | None = None
    if passthrough_indices is None:
        mapper = SpectralMapper(
            source_bands,
            target_bands,
            spectral_library=spectral_library,
            k_neighbors=int(k_neighbors),
        )
        logger.info("bestpixel surface prior: KNN spectral mapping (cross-sensor source basis)")
    else:
        logger.info(
            "bestpixel surface prior: identity band passthrough %s (no KNN remap)",
            target_band_names,
        )

    mapped_reflectance: list[np.ndarray] = []
    mapped_uncertainty_first: np.ndarray | None = None
    predictor_reflectance: list[np.ndarray] = []
    template_da: xr.DataArray | None = None
    surface_cfg = config.algorithms.surface_prior
    predict_visible = bool(getattr(surface_cfg, "bestpixel_predict_visible", False))
    predictor_band_order = ("coastal", "blue", "green", "red", "nir", "swir16", "swir22")
    predictor_indices: list[int] | None = None
    if predict_visible:
        band_lookup = {str(name): index for index, name in enumerate(bands)}
        missing = [name for name in predictor_band_order if name not in band_lookup]
        if missing:
            raise ValueError(
                "bestpixel_predict_visible requires the default visible-to-SWIR bestpixel bands; "
                f"missing: {', '.join(missing)}"
            )
        predictor_indices = [band_lookup[name] for name in predictor_band_order]

    for period in periods:
        reprojected = reproject_period_to_template(period, target_template, bands)
        if reprojected is None:
            continue
        reflectance_da, uncertainty_da = reprojected
        if predictor_indices is not None:
            predictor_reflectance.append(
                np.asarray(reflectance_da.isel(band=predictor_indices).values, dtype=np.float32)
            )
        if passthrough_indices is not None:
            mapped_refl = (
                reflectance_da.isel(band=passthrough_indices)
                .assign_coords(band=target_band_names)
                .astype(np.float32)
            )
            mapped_unc = (
                uncertainty_da.isel(band=passthrough_indices)
                .assign_coords(band=target_band_names)
                .astype(np.float32)
            )
        else:
            assert mapper is not None
            mapped_refl, mapped_unc, _fit = mapper.map(
                reflectance_da, source_uncertainty=uncertainty_da
            )
        mapped_reflectance.append(np.asarray(mapped_refl.values, dtype=np.float32))
        if mapped_uncertainty_first is None:
            mapped_uncertainty_first = np.asarray(mapped_unc.values, dtype=np.float32)
            template_da = mapped_refl

    if not mapped_reflectance or template_da is None:
        raise RuntimeError(
            "bestpixel surface prior produced no usable composites over the scene AOI."
        )

    robust_clip = float(config.algorithms.surface_prior.bestpixel_robust_clip)
    boa_values, boa_unc_values = reduce_realizations(
        np.stack(mapped_reflectance, axis=0),
        mapped_uncertainty_first,
        robust_clip=robust_clip,
    )

    boa = xr.DataArray(boa_values, dims=template_da.dims, coords=template_da.coords, name="boa")
    boa_unc = xr.DataArray(
        boa_unc_values, dims=template_da.dims, coords=template_da.coords, name="boa_unc"
    )
    mask_values = np.all(np.isfinite(boa_values), axis=0)
    spatial_reference = cast("xr.DataArray", template_da.isel(band=0, drop=True))
    mask = copy_spatial_metadata_like(
        xr.DataArray(
            mask_values,
            dims=spatial_reference.dims,
            coords=spatial_reference.coords,
            name="mask",
        ),
        spatial_reference,
    )
    prior = SurfacePrior(
        boa=boa,
        boa_unc=boa_unc,
        kernels=None,
        mask=mask,
        monthly_composites=(),
        solver_grid=target_template if used_native_grid else None,
        rt_space=library_rt_space,
    )
    if not predict_visible:
        return prior

    predictor_anchor_prior = anchor_atmo_prior if anchor_atmo_prior is not None else atmo_prior
    if predictor_anchor_prior is None:
        raise ValueError("bestpixel_predict_visible requires an atmospheric prior.")
    anchor_missing = [name for name in ("B8A", "B11", "B12") if name not in observation.toa]
    if anchor_missing:
        raise ValueError(
            "bestpixel_predict_visible requires scene TOA anchor bands: "
            f"{', '.join(anchor_missing)}"
        )
    if not predictor_reflectance:
        raise RuntimeError("bestpixel visible prediction had no usable monthly composite stack.")
    aot_values = np.asarray(predictor_anchor_prior.aot.values, dtype=np.float64)
    finite_aot = aot_values[np.isfinite(aot_values)]
    if finite_aot.size == 0:
        raise ValueError("bestpixel_predict_visible requires a finite atmospheric-prior AOD.")
    anchor_aot = float(np.nanmedian(finite_aot))
    logger.info(
        "bestpixel visible predictor: anchor_source=%s anchor_aot=%.3f",
        "secondary_atmo" if anchor_atmo_prior is not None else "atmo_prior",
        anchor_aot,
    )
    from siac.algorithms.surface.seasonal_predictor import (
        DEFAULT_BAND_COLUMNS,
        seasonal_extra_tree_prior,
    )

    target_band_columns = {
        str(name): DEFAULT_BAND_COLUMNS[str(name)]
        for name in getattr(surface_cfg, "bestpixel_predict_visible_bands", ("B02", "B04"))
        if str(name) in DEFAULT_BAND_COLUMNS
    }
    epsg, transform = _template_epsg_transform(spatial_reference)
    debias_cfg = getattr(surface_cfg, "bestpixel_predict_visible_debias", None)
    debias = (
        {str(band): (float(pair[0]), float(pair[1])) for band, pair in debias_cfg.items()}
        if debias_cfg
        else None
    )
    return seasonal_extra_tree_prior(
        prior,
        observation,
        seasonal_composites=np.stack(predictor_reflectance, axis=0),
        epsg=epsg,
        transform=transform,
        anchor_aot=anchor_aot,
        anchor_aot_field=predictor_anchor_prior.aot,
        target_band_columns=target_band_columns,
        debias=debias,
        predictor_model=str(getattr(surface_cfg, "bestpixel_predict_visible_model", "extra_tree")),
        uncertainty_floor=float(
            getattr(surface_cfg, "bestpixel_predict_visible_uncertainty_floor", 0.006)
        ),
        b01_uncertainty_floor=float(
            getattr(surface_cfg, "bestpixel_predict_visible_uncertainty_floor", 0.006)
        ),
        atmo_prior=atmo_prior,
        rt_model=rt_model,
        attach_tau_predictor=(
            bool(
                getattr(
                    getattr(getattr(config, "algorithms", None), "solver", None),
                    "surface_driven_tau_dependent_prior",
                    False,
                )
            )
        ),
        ensemble_aggregation=str(
            getattr(surface_cfg, "bestpixel_predict_visible_aggregation", "median")
        ),
        anchor_match_scale=float(
            getattr(surface_cfg, "bestpixel_predict_visible_anchor_match_scale", 0.05)
        ),
        scene_mean_geometry=bool(
            getattr(
                getattr(getattr(config, "algorithms", None), "solver", None),
                "surface_driven_scene_mean_geometry",
                False,
            )
        ),
    )
