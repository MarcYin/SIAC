"""Thin adapter over the external ``bestpixel`` best-pixel composite package.

``bestpixel`` (PyPI, Rust-backed) builds cloud-free monthly best-pixel
composites directly from STAC/COG sources — Planetary Computer / Earth Search
Sentinel-2 L2A, HLS, or MCD43A4 — and returns them as in-memory ``numpy``
arrays given an area of interest and a set of composite periods. This adapter
drives :func:`bestpixel.build_composite` once per month behind SIAC's
:class:`~siac.domain.protocols.MonthlyCompositeProvider` protocol so those
composites can serve as the Route-B monthly surface-reflectance prior, as an
alternative to the MODIS-BRDF–generated composites.

The ``bestpixel`` package is an *optional* dependency: it is imported lazily
inside :meth:`BestPixelMonthlyCompositeProvider.get_monthly_composites` so that
importing SIAC (and running the rest of the test suite) never requires it.
"""

from __future__ import annotations

import calendar
import logging
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import xarray as xr

from siac.algorithms.grid.assembler import _build_target_template
from siac.algorithms.surface.brdf_monthly_composite import (
    MonthlyBestPixelComposite,
    MonthlyCompositeCollection,
)
from siac.domain import SensorBand
from siac.geo.reprojection import reproject_match, transform_bounds
from siac.runtime.models import copy_spatial_metadata_like

if TYPE_CHECKING:
    from collections.abc import Sequence

    from siac.runtime import ObservationBundle

logger = logging.getLogger(__name__)

#: bestpixel returns integer DN; reflectance = DN * this factor.
_BESTPIXEL_DN_SCALE = 1.0e-4

#: Max concurrent per-period bestpixel fetches. Each ``build_composite`` is
#: itself internally concurrent, so this caps the outer fan-out (a 5-year,
#: single-month prior runs 5 at once; a 12-month config is bounded here so we
#: don't launch 60 simultaneous STAC sessions).
_MAX_FETCH_WORKERS = 8

#: Total concurrent COG connections to spread across the parallel period
#: fetches. bestpixel's per-call default (192) is fine for ONE fetch, but with
#: N periods in flight that becomes N*192 simultaneous S3 connections — e.g. a
#: 5-period S2 prior opened ~960, which earth-search throttles hard (the fetch
#: that should take seconds then crawls for an hour). We divide this budget by
#: the worker count so total in-flight connections stay ~constant regardless of
#: how many periods run at once.
_FETCH_CONCURRENCY_BUDGET = 192

# bestpixel quality codes (uint16), per its ``build_composite`` contract.
_QUALITY_CLEAR = 0
_QUALITY_MARGINAL = 1
_QUALITY_DARK = 2
_QUALITY_NODATA = 65535

#: Map bestpixel quality codes to SIAC's composite "quality cost" (lower is
#: better — the Route-B database filter drops pixels whose cost exceeds
#: ``surface_prior.monthly_database_filter.max_composite_quality``, default
#: 0.05). Clear/marginal/dark all stay under that default so they are
#: retained; nodata becomes NaN and is dropped.
_QUALITY_COST: dict[int, float] = {
    _QUALITY_CLEAR: 0.0,
    _QUALITY_MARGINAL: 0.02,
    _QUALITY_DARK: 0.04,
}

#: Canonical bestpixel band name -> (center wavelength nm, bandwidth nm).
#: Nominal Sentinel-2 values. bestpixel uses these names across endpoints and
#: SIAC's spectral mapper projects them onto the target sensor's bands by
#: wavelength, so these only need to be accurate enough to drive the mapping.
_BESTPIXEL_BAND_SPECS: dict[str, tuple[float, float]] = {
    "coastal": (443.0, 20.0),
    "blue": (492.0, 65.0),
    "green": (560.0, 35.0),
    "red": (665.0, 30.0),
    "rededge1": (704.0, 15.0),
    "rededge2": (740.0, 15.0),
    "rededge3": (783.0, 20.0),
    "nir": (833.0, 105.0),
    "nir08": (865.0, 20.0),
    "nir09": (945.0, 20.0),
    "swir16": (1614.0, 90.0),
    "swir22": (2202.0, 175.0),
}

#: Canonical bestpixel band name -> Sentinel-2 ``rsrf_band_id`` for the
#: ``sentinel-2a_msi`` published source basis. Only B01/B02/B03/B04/B08/B11/B12
#: have a published source entry; the narrow red-edge / water-vapour bands are
#: ignored by the mapper if requested.
_BESTPIXEL_TO_S2_BAND_ID: dict[str, str] = {
    "coastal": "B01",
    "blue": "B02",
    "green": "B03",
    "red": "B04",
    "rededge1": "B05",
    "rededge2": "B06",
    "rededge3": "B07",
    "nir": "B08",
    "nir08": "B8A",
    "nir09": "B09",
    "swir16": "B11",
    "swir22": "B12",
}

#: Canonical bestpixel band name -> MODIS ``rsrf_band_id`` for the
#: ``terra_modis`` published source basis (the mcd43a4 endpoint returns MODIS
#: NBAR). MODIS has no 443 nm coastal band, so ``coastal`` is unmapped.
_BESTPIXEL_TO_MODIS_BAND_ID: dict[str, str] = {
    "blue": "B3",
    "green": "B4",
    "red": "B1",
    "nir": "B2",
    "swir16": "B6",
    "swir22": "B7",
}

#: Endpoint -> (published source sensor unit, canonical-band -> rsrf_band_id).
#: SIAC's spectral mapper is keyed by ``rsrf_sensor_unit_id``; tagging the
#: source bands correctly is what lets it project bestpixel reflectance onto
#: the scene's sensor. mcd43a4 returns MODIS NBAR (``terra_modis``); the S2 and
#: HLS endpoints return Sentinel-2 / S2-harmonised reflectance, mapped via the
#: ``sentinel-2a_msi`` basis. HLS is bandpass-adjusted toward S2, so S2A is the
#: right published basis for it too.
_SOURCE_BASIS: dict[str, tuple[str, dict[str, str]]] = {
    "mcd43a4": ("terra_modis", _BESTPIXEL_TO_MODIS_BAND_ID),
    "pc": ("sentinel-2a_msi", _BESTPIXEL_TO_S2_BAND_ID),
    "earth-search": ("sentinel-2a_msi", _BESTPIXEL_TO_S2_BAND_ID),
    "hls": ("sentinel-2a_msi", _BESTPIXEL_TO_S2_BAND_ID),
    "hls-s30": ("sentinel-2a_msi", _BESTPIXEL_TO_S2_BAND_ID),
    "hls2-s30": ("sentinel-2a_msi", _BESTPIXEL_TO_S2_BAND_ID),
    "auto": ("sentinel-2a_msi", _BESTPIXEL_TO_S2_BAND_ID),
}


def _source_basis(endpoint: str) -> tuple[str, dict[str, str]]:
    basis = _SOURCE_BASIS.get(endpoint)
    if basis is None:
        # Unknown/new endpoint: default to the S2 basis (the common case).
        return "sentinel-2a_msi", _BESTPIXEL_TO_S2_BAND_ID
    return basis


#: Default band set: spans visible→SWIR and covers SIAC's Route-B visible
#: (B01/B02/B04) + query (B08/B11/B12) target bands for spectral mapping.
_DEFAULT_BANDS: tuple[str, ...] = (
    "coastal",
    "blue",
    "green",
    "red",
    "nir",
    "swir16",
    "swir22",
)

#: Public alias for the default band set, shared by the bestpixel surface-prior
#: builder so it does not duplicate the visible→SWIR selection.
DEFAULT_BESTPIXEL_BANDS: tuple[str, ...] = _DEFAULT_BANDS


def bestpixel_source_bands(
    bands: Sequence[str],
    *,
    endpoint: str = "pc",
    native_resolution_m: float = 60.0,
) -> tuple[SensorBand, ...]:
    """Build the :class:`SensorBand` spectral basis for ``bands``.

    The bands are tagged with the published-library source sensor that matches
    ``endpoint`` (``terra_modis`` for mcd43a4, ``sentinel-2a_msi`` otherwise) so
    SIAC's spectral mapper projects them onto the scene's sensor correctly.

    Raises:
        ValueError: if a band name is unknown, or has no rsrf entry for the
            endpoint's source basis (e.g. ``coastal`` with mcd43a4).
    """
    sensor_unit_id, band_id_map = _source_basis(endpoint)
    resolved: list[SensorBand] = []
    for index, name in enumerate(bands):
        spec = _BESTPIXEL_BAND_SPECS.get(name)
        if spec is None:
            raise ValueError(
                f"Unknown bestpixel band {name!r}; known bands: {sorted(_BESTPIXEL_BAND_SPECS)}"
            )
        rsrf_band_id = band_id_map.get(name)
        if rsrf_band_id is None:
            raise ValueError(
                f"bestpixel band {name!r} has no {sensor_unit_id!r} source-basis entry "
                f"(endpoint={endpoint!r}); available: {sorted(band_id_map)}"
            )
        center, bandwidth = spec
        resolved.append(
            SensorBand(
                name=name,
                center_wavelength=center,
                bandwidth=bandwidth,
                resolution=float(native_resolution_m),
                band_index=index,
                rsrf_sensor_unit_id=sensor_unit_id,
                rsrf_band_id=rsrf_band_id,
            )
        )
    return tuple(resolved)


class BestPixelMonthlyCompositeProvider:
    """Serve monthly best-pixel composites from the ``bestpixel`` package.

    Implements the :class:`~siac.domain.protocols.MonthlyCompositeProvider`
    protocol. For each requested period, ``bestpixel`` is asked for a
    cloud-free best-pixel composite over the observation's AOI; the result is
    reprojected onto SIAC's Route-B target grid and wrapped as a
    :class:`MonthlyBestPixelComposite`.
    """

    def __init__(
        self,
        *,
        endpoint: str = "pc",
        bands: Sequence[str] | None = None,
        lookback_years: int = 5,
        months: Sequence[int] | None = None,
        output_crs: str = "utm",
        top_k: int = 3,
        max_cloud_cover: float = 90.0,
        resolution_m: float | None = None,
        fetch_resolution_m: float | None = None,
        disk_cache: str | None = None,
    ) -> None:
        if lookback_years <= 0:
            raise ValueError(f"lookback_years must be > 0, got {lookback_years}")
        self._endpoint = str(endpoint)
        self._bands: tuple[str, ...] = tuple(bands) if bands else _DEFAULT_BANDS
        self._lookback_years = int(lookback_years)
        # ``None`` defers month selection to the scene: a surface prior should
        # be seasonally matched, so by default we composite only the
        # observation's own calendar month across the lookback years (e.g. a
        # March scene -> March of each prior year). Set ``months`` explicitly
        # to widen the seasonal window (e.g. (2, 3, 4)) or composite a full year.
        self._months: tuple[int, ...] | None = tuple(int(m) for m in months) if months else None
        if self._months is not None:
            for month in self._months:
                if not 1 <= month <= 12:
                    raise ValueError(f"months must be in 1..12, got {month}")
        self._output_crs = str(output_crs)
        self._top_k = int(top_k)
        self._max_cloud_cover = float(max_cloud_cover)
        self._resolution_m = float(resolution_m) if resolution_m is not None else None
        self._fetch_resolution_m = (
            float(fetch_resolution_m) if fetch_resolution_m is not None else None
        )
        self._disk_cache = disk_cache
        # Validate band names + tag with the endpoint's source basis eagerly so
        # config errors (e.g. coastal requested with mcd43a4) surface at build
        # time rather than mid-run.
        self._source_bands = bestpixel_source_bands(
            self._bands,
            endpoint=self._endpoint,
            native_resolution_m=self._resolution_m or 60.0,
        )

    @property
    def source_name(self) -> str:
        return f"bestpixel:{self._endpoint} monthly composites"

    @property
    def source_bands(self) -> tuple[SensorBand, ...]:
        return self._source_bands

    @property
    def resolution_m(self) -> float | None:
        return self._resolution_m

    def get_monthly_composites(
        self,
        observation: ObservationBundle,
        resolution: float,
    ) -> MonthlyCompositeCollection:
        obs_time = observation.metadata.get("observation_time")
        if not isinstance(obs_time, datetime):
            raise ValueError(
                "bestpixel composites require observation.metadata['observation_time'] "
                "to be a datetime."
            )
        years = self._resolve_years(obs_time)
        months = self._resolve_months(obs_time)
        # bestpixel takes a WGS84 (lon/lat) bbox; the observation grid is the
        # scene's (typically UTM) CRS, so reproject the bounds first.
        west, south, east, north = transform_bounds(
            observation.bounds, observation.crs, "EPSG:4326"
        )
        bbox = [float(west), float(south), float(east), float(north)]

        import bestpixel as bp  # lazy optional dependency

        database_res = float(resolution)
        # Optionally fetch finer than the prior grid and area-average down, to
        # trade resolution for a smoother (lower-noise) prior.
        fetch_res = (
            self._fetch_resolution_m
            if (self._fetch_resolution_m is not None and self._fetch_resolution_m < database_res)
            else database_res
        )
        downsample = fetch_res < database_res

        target_template = _build_target_template(observation.bounds, observation.crs, database_res)
        logger.info(
            "bestpixel: building composites endpoint=%s bbox=%s years=%s months=%s "
            "fetch_res=%.1f prior_res=%.1f%s bands=%s grid=%dx%d",
            self._endpoint,
            [round(v, 4) for v in bbox],
            list(years),
            list(months),
            float(fetch_res),
            database_res,
            " (area-averaged)" if downsample else "",
            list(self._bands),
            int(target_template.sizes["y"]),
            int(target_template.sizes["x"]),
        )

        # NOTE: We call ``build_composite`` per (year, month) with a
        # calendar-correct end date rather than ``build_monthly_composites``.
        # bestpixel 0.1.1's batch helper builds every period's date range as
        # ``YYYY-MM-01/YYYY-MM-31``, so any month shorter than 31 days (Feb +
        # the 30-day months) is rejected by the STAC backend with a 400.
        # Per-period requests with the real last-of-month sidestep the bug and
        # let us skip individual periods that have no scene coverage.
        #
        # The periods are independent network fetches, so we run them
        # concurrently in a thread pool (bestpixel is thread-safe and releases
        # the GIL during its Rust-side fetch; measured ~11x faster than serial
        # for a 5-period prior). ``ThreadPoolExecutor.map`` preserves input
        # order, so the resulting composite order is deterministic.
        period_specs = [(year, month) for year in years for month in months]
        n_workers = min(len(period_specs), _MAX_FETCH_WORKERS)
        # Divide the connection budget across the in-flight periods so total
        # concurrent S3 connections stay bounded (avoids the earth-search
        # throttling that turned a seconds-long fetch into an hour).
        per_fetch_concurrency = max(8, _FETCH_CONCURRENCY_BUDGET // n_workers)
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            fetched = list(
                executor.map(
                    lambda spec: self._fetch_period(
                        bp, bbox, spec[0], spec[1], fetch_res, per_fetch_concurrency
                    ),
                    period_specs,
                )
            )

        # Area-average when fetching finer than the prior grid; bilinear when
        # the grids match (no resolution change).
        reflectance_resampling = "average" if downsample else "bilinear"
        composites: list[MonthlyBestPixelComposite] = []
        skipped: list[tuple[int, int]] = []
        for (year, month), period in zip(period_specs, fetched, strict=True):
            composite = (
                self._period_to_composite(
                    period,
                    target_template,
                    year=year,
                    month=month,
                    reflectance_resampling=reflectance_resampling,
                )
                if period is not None
                else None
            )
            if composite is None:
                skipped.append((year, month))
                continue
            composites.append(composite)

        logger.info(
            "bestpixel: built %d monthly composite(s); skipped %d period(s) with no "
            "coverage or fetch errors",
            len(composites),
            len(skipped),
        )
        if not composites:
            raise RuntimeError(
                f"bestpixel returned no usable composites for endpoint={self._endpoint!r}, "
                f"years={list(years)}, months={list(months)} over the scene AOI. "
                "Try a different endpoint (e.g. 'earth-search') or composite periods."
            )
        return MonthlyCompositeCollection(
            composites=tuple(composites),
            source_bands=self._source_bands,
            source_name=self.source_name,
        )

    def _resolve_years(self, obs_time: datetime) -> tuple[int, ...]:
        """The ``lookback_years`` full calendar years before the scene year."""
        end_year = obs_time.year - 1
        start_year = end_year - self._lookback_years + 1
        return tuple(range(start_year, end_year + 1))

    def _resolve_months(self, obs_time: datetime) -> tuple[int, ...]:
        """Composite months: the configured set, else the scene's own month.

        A surface prior should be seasonally matched to the scene, so the
        default (when ``months`` is unset) is just the observation's calendar
        month — e.g. a March scene composites March of each lookback year.
        """
        if self._months is not None:
            return self._months
        return (int(obs_time.month),)

    def _fetch_period(
        self,
        bp: Any,
        bbox: list[float],
        year: int,
        month: int,
        resolution: float,
        concurrency: int,
    ) -> dict[str, Any] | None:
        """Fetch one monthly composite, returning ``None`` on any failure."""
        last_day = calendar.monthrange(year, month)[1]
        date_range = f"{year:04d}-{month:02d}-01/{year:04d}-{month:02d}-{last_day:02d}"
        try:
            return cast(
                "dict[str, Any]",
                bp.build_composite(
                    bbox=bbox,
                    datetime=date_range,
                    resolution=resolution,
                    top_k=self._top_k,
                    max_cloud_cover=self._max_cloud_cover,
                    concurrency=int(concurrency),
                    endpoint=self._endpoint,
                    disk_cache=self._disk_cache,
                    bands=list(self._bands),
                    output_crs=self._output_crs,
                ),
            )
        except Exception as exc:  # noqa: BLE001 - network/STAC failures are skippable
            logger.warning(
                "bestpixel: period %04d-%02d failed (%s: %s); skipping",
                year,
                month,
                type(exc).__name__,
                exc,
            )
            return None

    def _period_to_composite(
        self,
        period: dict[str, Any],
        target_template: xr.DataArray,
        *,
        year: int,
        month: int,
        reflectance_resampling: str = "bilinear",
    ) -> MonthlyBestPixelComposite | None:
        grid = period["grid"]
        source_crs = str(grid.get("crs") or f"EPSG:{int(grid['epsg'])}")
        x_coords, y_coords = _coords_from_grid(grid)

        bands_dict = period["bands"]
        quality = np.asarray(period["quality"])
        height, width = quality.shape
        nodata_mask = quality == _QUALITY_NODATA

        reflectance = np.empty((len(self._bands), height, width), dtype=np.float32)
        for index, name in enumerate(self._bands):
            if name not in bands_dict:
                raise KeyError(
                    f"bestpixel composite is missing requested band {name!r}; "
                    f"returned bands: {sorted(bands_dict)}"
                )
            reflectance[index] = (
                np.asarray(bands_dict[name], dtype=np.float32) * _BESTPIXEL_DN_SCALE
            )
        reflectance[:, nodata_mask] = np.nan

        # A period with no scene coverage comes back entirely nodata; skip it
        # so the Route-B database isn't built from empty months.
        if not np.isfinite(reflectance).any():
            logger.info("bestpixel: period %04d-%02d has no valid pixels; skipping", year, month)
            return None

        reflectance_da = _make_source_da(
            reflectance, x_coords, y_coords, source_crs, band_names=self._bands
        )
        quality_cost = _quality_to_cost(quality)
        quality_da = _make_source_da(quality_cost, x_coords, y_coords, source_crs)

        # Average quality cost too when downsampling so it reflects the cell's
        # mean rather than a single sub-pixel sample.
        quality_resampling = "average" if reflectance_resampling == "average" else "nearest"
        reflectance_aligned = reproject_match(
            reflectance_da, target_template, resampling=reflectance_resampling, nodata=np.nan
        ).astype(np.float32)
        quality_aligned = reproject_match(
            quality_da, target_template, resampling=quality_resampling, nodata=np.nan
        ).astype(np.float32)
        sample_index = _sample_index_from(reflectance_aligned, quality_aligned)

        return MonthlyBestPixelComposite(
            reflectance=reflectance_aligned,
            quality=quality_aligned,
            sample_index=sample_index,
            year=int(year),
            month=int(month),
        )


def _coords_from_grid(grid: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    """Pixel-center x/y coordinate vectors from a bestpixel ``grid`` dict."""
    transform = [float(v) for v in grid["transform"]]
    x_step, _, x_origin, _, y_step, y_origin = transform
    width = int(grid["width"])
    height = int(grid["height"])
    x = x_origin + (np.arange(width, dtype=np.float64) + 0.5) * x_step
    y = y_origin + (np.arange(height, dtype=np.float64) + 0.5) * y_step
    return x, y


def _make_source_da(
    values: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    crs: str,
    *,
    band_names: Sequence[str] | None = None,
) -> xr.DataArray:
    if band_names is None:
        data = xr.DataArray(values, dims=("y", "x"), coords={"y": y, "x": x})
    else:
        data = xr.DataArray(
            values,
            dims=("band", "y", "x"),
            coords={"band": np.asarray(list(band_names), dtype=object), "y": y, "x": x},
        )
    data = data.rio.set_spatial_dims(x_dim="x", y_dim="y")
    return data.rio.write_crs(crs)


def _quality_to_cost(quality: np.ndarray) -> np.ndarray:
    """Map bestpixel quality codes to a float cost array (NaN at nodata)."""
    cost = np.full(quality.shape, np.nan, dtype=np.float32)
    mapped = np.zeros(quality.shape, dtype=bool)
    for code, value in _QUALITY_COST.items():
        match = quality == code
        cost[match] = value
        mapped |= match
    # Any non-nodata code we didn't enumerate is treated as the worst
    # retained ("dark") tier rather than silently dropped.
    unknown = ~mapped & (quality != _QUALITY_NODATA)
    cost[unknown] = _QUALITY_COST[_QUALITY_DARK]
    return cast("np.ndarray", cost)


def _sample_index_from(
    reflectance: xr.DataArray,
    quality: xr.DataArray,
) -> xr.DataArray:
    """0 where a pixel carries a finite composite payload, -1 otherwise."""
    finite_reflectance = np.all(np.isfinite(np.asarray(reflectance.values)), axis=0)
    finite_quality = np.isfinite(np.asarray(quality.values))
    valid = finite_reflectance & finite_quality
    index = np.where(valid, 0, -1).astype(np.int16)
    return copy_spatial_metadata_like(
        xr.DataArray(
            index,
            dims=("y", "x"),
            coords={"y": quality.coords["y"], "x": quality.coords["x"]},
        ),
        quality,
    )
