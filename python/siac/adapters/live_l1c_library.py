"""Surface library built live from Sentinel-2 L1C top-of-atmosphere data.

This is the in-package port of the validated research builder that produced the
faithful ``acix3`` surface dictionary. It replaces the offline three-script chain

``acix3_index_prototype.py`` -> ``build_index_l1c_toa_gcs.py`` ->
``convert_acix3_faithful_dict_6scci.py``

with one :class:`~siac.adapters.surface_library.SurfaceLibrary` implementation
that runs the last two stages for an arbitrary scene at retrieval time.

Why the shape of it matters:

* A monthly composite is not one acquisition. A *winner index* records, per
  pixel, which of the month's acquisitions won the cloud-score mosaic; only two
  or three acquisitions typically win any month, so a faithful mosaic costs a
  handful of windowed reads rather than a full month of imagery.
* The library is corrected here, in the solver's own RT model, at each winning
  day's own aerosol loading. A library corrected in one RT space and solved
  against in another injects a systematic surface offset that the solver turns
  into biased AOT — see :mod:`siac.domain.rt_space`.

The winner index itself is an input (``months``/``winners``/``image_table``/
``day_scalars`` in one ``.npz`` per scene), because building it needs a
cloud-score source outside this package. Everything downstream of it is
GEE-free: L1C top-of-atmosphere reflectance is read straight from the anonymous
``gs://gcp-public-data-sentinel-2`` bucket.

The native 6S backend evaluates one atmospheric state per invocation, so the
correction builds coefficient *curves* over an AOD node axis at scalar median
geometry and interpolates them per pixel — the same construction the validated
converter used, and the same one the solver's coefficient provider uses.
"""

from __future__ import annotations

import json
import logging
import re
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol, cast

import numpy as np

from siac.adapters.data import gcs_sentinel2 as gcs
from siac.adapters.data.gcs_sentinel2 import GCSSentinel2Backend
from siac.adapters.data.s2_data_source import S2Query
from siac.adapters.surface_library import (
    PREDICTOR_BAND_ORDER,
    SurfaceLibraryRealization,
    observation_scene_key,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    import xarray as xr

    from siac.adapters.data.s2_data_source import S2Product
    from siac.domain import SensorBand
    from siac.domain.rt_space import RTSpace
    from siac.runtime.models import ObservationBundle

logger = logging.getLogger(__name__)

__all__ = [
    "LIBRARY_BAND_NAMES",
    "GCSL1CTOAReader",
    "IndexImage",
    "L1CTOAReader",
    "LiveL1CSurfaceLibrary",
    "MosaicIndex",
    "correct_toa_realization",
    "median_scalar",
    "read_mosaic_index",
    "resolve_mgrs_tile",
]

#: Sentinel-2 band names read from L1C, ordered to match
#: :data:`~siac.adapters.surface_library.PREDICTOR_BAND_ORDER`
#: (coastal, blue, green, red, nir, swir16, swir22).
LIBRARY_BAND_NAMES: tuple[str, ...] = ("B01", "B02", "B03", "B04", "B8A", "B11", "B12")

#: L1C digital-number scale factor (DN -> reflectance).
_L1C_DN_SCALE = 10000.0
#: Processing baselines from this number on carry a radiometric DN offset.
_OFFSET_BASELINE = 400
#: The DN offset those baselines apply.
_BASELINE_OFFSET_DN = 1000.0
#: AOD used for a winning day the index never scored.
_DEFAULT_DAY_AOD = 0.1
#: Number of AOD nodes the 6S coefficient curves are evaluated on.
_AOD_NODE_COUNT = 8
#: AOD node-axis bounds (the native 6S aerosol axis is only valid inside these).
_AOD_NODE_MIN = 0.01
_AOD_NODE_MAX = 4.0
#: Minimum finite fraction a corrected realization must reach to be kept.
_KEEP_FRACTION = 0.6
#: Realizations always kept even if none reaches ``_KEEP_FRACTION``.
_MIN_REALIZATIONS = 4
#: Fallback correction state when no atmospheric prior is supplied.
_DEFAULT_TCWV_CM = 2.0
_DEFAULT_TCO3_ATM_CM = 0.3
#: Default parallelism for the L1C reads (only the winning days are read).
_DEFAULT_MAX_WORKERS = 8

#: GDAL settings for windowed ``/vsicurl`` reads of the public JP2 objects.
_GDAL_VSICURL_ENV = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".jp2",
    "GDAL_HTTP_MAX_RETRY": "3",
    "GDAL_HTTP_RETRY_DELAY": "2",
}

_MGRS_IN_TEXT = re.compile(r"T(\d{2}[A-Z]{3})")
_MGRS_BARE = re.compile(r"^T?(\d{2}[A-Z]{3})$")


@dataclass(frozen=True)
class IndexImage:
    """One acquisition the winner index can select, with its mean geometry."""

    image_id: str
    day: str
    sza: float
    saa: float
    vza: float
    vaa: float

    @property
    def sensing_token(self) -> str:
        """``YYYYMMDDTHHMMSS`` sensing token used to find the L1C SAFE product."""

        return self.image_id.split("_")[0]


@dataclass(frozen=True)
class MosaicIndex:
    """Per-pixel winning-acquisition index for one scene's monthly mosaics.

    ``winners`` is ``(month, y, x)``; each value indexes ``images_by_month`` for
    that month, with ``-1`` meaning no acquisition was clear at that pixel.
    """

    months: tuple[str, ...]
    winners: np.ndarray
    images_by_month: dict[str, tuple[IndexImage, ...]]
    day_aod: dict[str, float]


def read_mosaic_index(path: str | Path) -> MosaicIndex:
    """Read a winner-index ``.npz`` written by the mosaic-index builder.

    The stored ``image_table`` is one flat, month-ordered list of JSON records;
    it is split back into per-month blocks by walking it in order, exactly as the
    research loader did, because the winner values index those blocks.
    """

    index_path = Path(path)
    with np.load(index_path, allow_pickle=False) as payload:
        months = tuple(str(m) for m in np.asarray(payload["months"]).ravel().tolist())
        winners = np.asarray(payload["winners"], dtype=np.int16)
        table = [json.loads(str(s)) for s in np.asarray(payload["image_table"]).ravel().tolist()]
        scalars = [json.loads(str(s)) for s in np.asarray(payload["day_scalars"]).ravel().tolist()]

    if winners.ndim != 3 or winners.shape[0] != len(months):
        raise ValueError(
            f"Mosaic index {index_path} must hold a (month, y, x) 'winners' array matching its "
            f"{len(months)} months; got shape {winners.shape}."
        )

    images_by_month: dict[str, tuple[IndexImage, ...]] = {}
    cursor = 0
    for month in months:
        block: list[IndexImage] = []
        for record in table[cursor:]:
            if str(record["day"])[:7] != month:
                break
            block.append(
                IndexImage(
                    image_id=str(record["idx"]),
                    day=str(record["day"]),
                    sza=float(record["sza"]),
                    saa=float(record["saa"]),
                    vza=float(record["vza"]),
                    vaa=float(record["vaa"]),
                )
            )
        images_by_month[month] = tuple(block)
        cursor += len(block)

    day_aod = {str(item["day"]): float(item["aod"]) for item in scalars}
    return MosaicIndex(
        months=months,
        winners=winners,
        images_by_month=images_by_month,
        day_aod=day_aod,
    )


def resolve_mgrs_tile(observation: ObservationBundle, scene_key: str | None = None) -> str:
    """Find the scene's MGRS tile from ``scene_key`` or observation metadata."""

    metadata = observation.metadata or {}
    candidates = [
        scene_key,
        metadata.get("mgrs_tile"),
        metadata.get("tile_id"),
        metadata.get("product_id"),
        metadata.get("scene_key"),
        metadata.get("scene_id"),
    ]
    for candidate in candidates:
        if not candidate:
            continue
        text = str(candidate).upper()
        bare = _MGRS_BARE.match(text)
        if bare is not None:
            return bare.group(1)
        found = _MGRS_IN_TEXT.search(text)
        if found is not None:
            return found.group(1)
    raise ValueError(
        "Live L1C surface library needs the scene's MGRS tile: set observation metadata "
        "'mgrs_tile' (or 'tile_id'/'product_id'), or use a scene key containing it "
        "(e.g. '..._T31UDQ_20240101T104321')."
    )


def median_scalar(field: Any) -> float | None:
    """Finite median of an array-like / DataArray field, or ``None`` if empty."""

    if field is None:
        return None
    values = np.asarray(getattr(field, "values", field), dtype=np.float64)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return None
    return float(np.median(finite))


class L1CTOAReader(Protocol):
    """Reads one acquisition's L1C top-of-atmosphere reflectance onto a grid."""

    def read(
        self,
        *,
        mgrs_tile: str,
        sensing_token: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
    ) -> np.ndarray | None:
        """Return ``(band, y, x)`` reflectance in :data:`LIBRARY_BAND_NAMES` order.

        Returns ``None`` when the acquisition cannot be located or read; a single
        missing acquisition degrades that month's mosaic rather than the run.
        """


def _baseline_offset(processing_baseline: str | None) -> float:
    """Radiometric DN offset implied by a SAFE processing baseline."""

    text = str(processing_baseline or "").lstrip("Nn")
    if not text.isdigit():
        return 0.0
    return _BASELINE_OFFSET_DN if int(text) >= _OFFSET_BASELINE else 0.0


def _band_jp2_urls(product: S2Product) -> dict[str, str]:
    """Map each library band to its ``/vsicurl`` JP2 URL inside a SAFE product."""

    prefix = gcs._resolve_safe_prefix(product)
    objects = gcs._list_objects_under(prefix)
    urls: dict[str, str] = {}
    for band in LIBRARY_BAND_NAMES:
        for item in objects:
            name = str(item.get("name", ""))
            if name.endswith(f"_{band}.jp2") and "/IMG_DATA/" in name:
                urls[band] = "/vsicurl/" + gcs._object_download_url(name)
                break
    return urls


def _read_toa_window(
    urls: dict[str, str],
    *,
    offset: float,
    crs: str,
    transform: tuple[float, float, float, float, float, float],
    width: int,
    height: int,
    resampling: str = "average",
) -> np.ndarray:
    """Windowed-read + warp each band's JP2 straight onto the target grid.

    A :class:`~rasterio.vrt.WarpedVRT` is used rather than
    :mod:`siac.geo.reprojection` because the source is a remote 110 km JP2 and
    only the AOI window is wanted; the VRT reads solely the overview blocks the
    target grid touches.
    """

    import rasterio
    from rasterio.transform import Affine
    from rasterio.vrt import WarpedVRT

    from siac.geo.reprojection import RESAMPLING_METHODS

    method = RESAMPLING_METHODS[resampling]
    planes: np.ndarray = np.full((len(LIBRARY_BAND_NAMES), height, width), np.nan, dtype=np.float32)
    with rasterio.Env(**_GDAL_VSICURL_ENV):
        for index, band in enumerate(LIBRARY_BAND_NAMES):
            url = urls.get(band)
            if url is None:
                continue
            with (
                rasterio.open(url) as src,
                WarpedVRT(
                    src,
                    crs=crs,
                    transform=Affine(*transform),
                    width=width,
                    height=height,
                    resampling=method,
                ) as vrt,
            ):
                dn = vrt.read(1).astype(np.float32)
            dn[dn == 0] = np.nan
            planes[index] = (dn - offset) / _L1C_DN_SCALE
    return planes


class GCSL1CTOAReader:
    """Reads L1C TOA from the anonymous ``gcp-public-data-sentinel-2`` bucket.

    The tile listing is fetched once per MGRS tile and cached, then acquisitions
    are matched by their ``YYYYMMDDTHHMMSS`` sensing token. Where a tile holds
    several SAFE products for one sensing time (reprocessing baselines), the
    first by product id is used, reproducing the research builder's choice.
    """

    def __init__(
        self,
        backend: Any | None = None,
        *,
        resampling: str = "average",
    ) -> None:
        self._backend = backend if backend is not None else GCSSentinel2Backend()
        self._resampling = str(resampling)
        self._by_tile: dict[str, dict[str, S2Product]] = {}
        self._lock = threading.Lock()

    def _products_for_tile(self, mgrs_tile: str) -> dict[str, S2Product]:
        with self._lock:
            cached = self._by_tile.get(mgrs_tile)
        if cached is not None:
            return cached
        products = self._backend.search(S2Query(mgrs_tile=mgrs_tile, processing_level="L1C"))
        by_token: dict[str, S2Product] = {}
        for product in sorted(products, key=lambda item: str(item.product_id)):
            by_token.setdefault(product.sensing_date.strftime("%Y%m%dT%H%M%S"), product)
        with self._lock:
            self._by_tile[mgrs_tile] = by_token
        logger.info(
            "Live L1C library: %d L1C products listed for tile %s", len(by_token), mgrs_tile
        )
        return by_token

    def read(
        self,
        *,
        mgrs_tile: str,
        sensing_token: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
    ) -> np.ndarray | None:
        product = self._products_for_tile(mgrs_tile).get(sensing_token)
        if product is None:
            logger.info(
                "Live L1C library: no L1C product on GCS for %s at %s", mgrs_tile, sensing_token
            )
            return None
        urls = _band_jp2_urls(product)
        if not urls:
            logger.warning("Live L1C library: no JP2 bands under %s", product.product_id)
            return None
        return _read_toa_window(
            urls,
            offset=_baseline_offset(product.processing_baseline),
            crs=crs,
            transform=transform,
            width=width,
            height=height,
            resampling=self._resampling,
        )


def _scalar_field(value: float) -> xr.DataArray:
    """One-pixel ``(y, x)`` field: the native 6S backend runs per state."""

    import xarray as xr

    return xr.DataArray(np.full((1, 1), float(value), dtype=np.float32), dims=("y", "x"))


def _aod_node_axis(aod: np.ndarray, count: int = _AOD_NODE_COUNT) -> np.ndarray:
    low = max(_AOD_NODE_MIN, float(aod.min()) * 0.8)
    high = min(_AOD_NODE_MAX, max(float(aod.max()) * 1.2, low + 0.05))
    nodes = np.unique(np.geomspace(low, high, int(count)).astype(np.float32))
    return cast("np.ndarray", nodes)


def correct_toa_realization(
    *,
    rt_model: Any,
    sensor_bands: Sequence[SensorBand],
    toa: np.ndarray,
    aod: np.ndarray,
    sza_deg: float,
    saa_deg: float,
    vza_deg: float,
    vaa_deg: float,
    tcwv: float,
    tco3: float,
    elevation_km: float,
    node_count: int = _AOD_NODE_COUNT,
) -> np.ndarray:
    """Atmospherically correct one ``(band, y, x)`` TOA mosaic through ``rt_model``.

    ``aod`` is the per-pixel aerosol loading of whichever acquisition won each
    pixel, so the correction is evaluated on an AOD node axis spanning the
    mosaic's range and interpolated per pixel. Geometry, water vapour and ozone
    are scalar medians: the native 6S backend runs one invocation per state, and
    discretising the smoothly varying fields there is what the validated
    converter did.
    """

    from siac.runtime.models import AtmosphericState, GeometryAngles

    toa_values = np.asarray(toa, dtype=np.float64)
    aod_values = np.asarray(aod, dtype=np.float64)
    if toa_values.ndim != 3 or toa_values.shape[0] != len(sensor_bands):
        raise ValueError(
            f"TOA must be (band, y, x) with {len(sensor_bands)} bands; got {toa_values.shape}."
        )

    finite_aod = np.isfinite(aod_values) & (aod_values > 0.0)
    fill = float(np.median(aod_values[finite_aod])) if finite_aod.any() else _DEFAULT_DAY_AOD
    aod_values = np.where(finite_aod, aod_values, fill)
    nodes = _aod_node_axis(aod_values, node_count)

    geometry = GeometryAngles.from_degrees(
        _scalar_field(sza_deg),
        _scalar_field(saa_deg),
        _scalar_field(vza_deg),
        _scalar_field(vaa_deg),
    )
    curves: np.ndarray = np.full((len(sensor_bands), nodes.size, 3), np.nan, dtype=np.float64)
    for node_index, node in enumerate(nodes):
        state = AtmosphericState(
            aot=_scalar_field(float(node)),
            tcwv=_scalar_field(tcwv),
            tco3=_scalar_field(tco3),
            aot_unc=_scalar_field(0.1),
            tcwv_unc=_scalar_field(0.5),
            tco3_unc=_scalar_field(0.05),
            elevation=_scalar_field(elevation_km),
        )
        coefficients = rt_model.compute_coefficients_multi(geometry, state, list(sensor_bands))
        for band_index, coefficient in enumerate(coefficients):
            curves[band_index, node_index] = (
                float(np.ravel(np.asarray(coefficient.xap))[0]),
                float(np.ravel(np.asarray(coefficient.xbp))[0]),
                float(np.ravel(np.asarray(coefficient.xcp))[0]),
            )

    corrected: np.ndarray = np.full(toa_values.shape, np.nan, dtype=np.float32)
    for band_index in range(len(sensor_bands)):
        xap = np.interp(aod_values, nodes, curves[band_index, :, 0])
        xbp = np.interp(aod_values, nodes, curves[band_index, :, 1])
        xcp = np.interp(aod_values, nodes, curves[band_index, :, 2])
        intermediate = xap * toa_values[band_index] - xbp
        with np.errstate(invalid="ignore", divide="ignore"):
            corrected[band_index] = (intermediate / (1.0 + xcp * intermediate)).astype(np.float32)
    return corrected


def _resample_winners(winners: np.ndarray, height: int, width: int) -> np.ndarray:
    """Nearest-neighbour the winner planes onto a ``(height, width)`` grid.

    The index file carries no CRS or transform, so a shape mismatch can only be
    resolved by assuming it covers the same footprint as the observation. That
    assumption is logged rather than silently applied.
    """

    if winners.shape[1:] == (height, width):
        return winners
    source_height, source_width = winners.shape[1:]
    logger.warning(
        "Live L1C library: winner index is %dx%d but the solver grid is %dx%d; "
        "resampling the index nearest-neighbour and assuming a shared AOI footprint.",
        source_height,
        source_width,
        height,
        width,
    )
    rows = np.clip(
        ((np.arange(height) + 0.5) * source_height / height).astype(int), 0, source_height - 1
    )
    columns = np.clip(
        ((np.arange(width) + 0.5) * source_width / width).astype(int), 0, source_width - 1
    )
    resampled = np.asarray(winners[:, rows[:, None], columns[None, :]], dtype=np.int16)
    return cast("np.ndarray", resampled)


def _template_grid(
    template: xr.DataArray,
) -> tuple[str, tuple[float, float, float, float, float, float], int, int]:
    x = np.asarray(template.coords["x"].values, dtype=float)
    y = np.asarray(template.coords["y"].values, dtype=float)
    x_step = float(abs(x[1] - x[0])) if x.size > 1 else 1.0
    y_step = float(abs(y[0] - y[1])) if y.size > 1 else x_step
    transform = (
        x_step,
        0.0,
        float(x[0] - 0.5 * x_step),
        0.0,
        -y_step,
        float(y[0] + 0.5 * y_step),
    )
    crs = template.rio.crs
    if crs is None:
        raise ValueError("Live L1C surface library requires a CRS on the solver grid template.")
    return str(crs), transform, int(x.size), int(y.size)


class LiveL1CSurfaceLibrary:
    """Surface library assembled live from L1C, corrected in the solver's RT space.

    For each month in the scene's winner index this fetches only the winning
    acquisitions' L1C top-of-atmosphere reflectance from the anonymous Google
    Cloud Sentinel-2 bucket, mosaics them through the index, and corrects the
    mosaic with ``rt_model`` at each winning day's own aerosol loading. One
    :class:`~siac.adapters.surface_library.SurfaceLibraryRealization` per month
    comes out, matching the prepared store the same construction writes offline.
    """

    def __init__(
        self,
        index_root: str | Path,
        *,
        rt_model: Any,
        solver_config: Any | None = None,
        scene_key: str | None = None,
        tcwv: float | None = None,
        tco3: float | None = None,
        elevation_km: float | None = None,
        maiac_day_aod: Callable[
            [tuple[float, float, float, float], str, Sequence[tuple[int, int]]], dict[str, float]
        ]
        | None = None,
        toa_reader: L1CTOAReader | None = None,
        keep_fraction: float = _KEEP_FRACTION,
        max_workers: int = _DEFAULT_MAX_WORKERS,
    ) -> None:
        self._index_root = Path(index_root).expanduser()
        self._rt_model = rt_model
        self._solver_config = solver_config
        self._scene_key = scene_key
        self._tcwv = _DEFAULT_TCWV_CM if tcwv is None else float(tcwv)
        self._tco3 = _DEFAULT_TCO3_ATM_CM if tco3 is None else float(tco3)
        self._elevation_km = 0.0 if elevation_km is None else float(elevation_km)
        self._maiac_day_aod = maiac_day_aod
        self._toa_reader = toa_reader
        self._keep_fraction = float(keep_fraction)
        self._max_workers = max(1, int(max_workers))

    @property
    def rt_space(self) -> RTSpace | None:
        """The RT space the library is corrected in — the solver's own.

        Reported through :meth:`~siac.domain.rt_space.RTSpace.for_solver` so a
        scene-adaptive aerosol mode is identified by its rule rather than by one
        scene's resolved mixture, exactly as the solve is.
        """

        from siac.domain.rt_space import RTSpace

        return RTSpace.for_solver(self._rt_model, self._solver_config)

    def _index_path(self, observation: ObservationBundle) -> Path:
        key = self._scene_key or observation_scene_key(observation)
        candidate = self._index_root / f"{key}.npz"
        if candidate.is_file():
            return candidate
        raise FileNotFoundError(
            f"Live L1C surface library has no mosaic index for scene {key!r} under "
            f"{self._index_root}. Build the index for this scene, or point "
            "surface_prior.live_l1c_index_path at the correct store."
        )

    def _reader(self) -> L1CTOAReader:
        if self._toa_reader is None:
            self._toa_reader = GCSL1CTOAReader()
        return self._toa_reader

    def _day_aod(self, index: MosaicIndex, observation: ObservationBundle) -> dict[str, float]:
        """Per-day AOD for the correction: the index's own scalars, else MAIAC."""

        if index.day_aod:
            return dict(index.day_aod)
        source = self._maiac_day_aod
        if source is None:
            from siac.adapters.atmo.maiac_day_aod import MAIACDayAODProvider

            source = MAIACDayAODProvider().day_aod_map
        periods = sorted({(int(m[:4]), int(m[5:7])) for m in index.months})
        try:
            return dict(source(observation.bounds, str(observation.crs), periods))
        except Exception:  # noqa: BLE001 - a failed gate must not abort the library
            logger.warning(
                "Live L1C library: per-day MAIAC AOD unavailable; correcting at AOD=%.2f",
                _DEFAULT_DAY_AOD,
                exc_info=True,
            )
            return {}

    def _effective_rt_model(self, observation: ObservationBundle) -> Any:
        """Apply the solver's scene-adaptive aerosol species to the RT backend.

        The backend handed to the surface prior still carries the base aerosol
        setup; the solver resolves the scene's CCI climatology mixture itself.
        The library must resolve the same one, or it is corrected in a different
        RT space than the solve it feeds.
        """

        mode = str(getattr(self._solver_config, "surface_driven_aerosol_species", "none") or "none")
        if mode == "none":
            return self._rt_model
        if mode != "cci_climatology_exact":
            raise ValueError(f"Unsupported surface-driven aerosol species mode {mode!r}.")

        base_setup = getattr(self._rt_model, "rt_setup", None)
        clone = getattr(self._rt_model, "with_rt_setup", None)
        if base_setup is None or clone is None or not hasattr(base_setup, "model_copy"):
            raise ValueError(
                "The live L1C surface library needs an RT backend exposing rt_setup and "
                "with_rt_setup() to match the solver's species mode."
            )
        from siac.algorithms.rt.aerosol_species import climatology_cci_aerosol_setup
        from siac.geo.reprojection import transform_bounds

        west, south, east, north = transform_bounds(
            observation.bounds, str(observation.crs), "EPSG:4326"
        )
        lon = 0.5 * (float(west) + float(east))
        lat = 0.5 * (float(south) + float(north))
        month = self._scene_month(observation)
        setup = climatology_cci_aerosol_setup(lon, lat, month)
        logger.info(
            "Live L1C library: correcting in the solver's %s aerosol space "
            "(lon=%.3f lat=%.3f month=%d)",
            mode,
            lon,
            lat,
            month,
        )
        return clone(base_setup.model_copy(update={"aerosol": setup}, deep=True))

    def _scene_month(self, observation: ObservationBundle) -> int:
        observation_time = getattr(self._rt_model, "observation_time", None)
        if observation_time is None:
            observation_time = (observation.metadata or {}).get("observation_time")
        month = getattr(observation_time, "month", None)
        if month is None:
            raise ValueError(
                "The live L1C surface library needs the scene's observation time to resolve the "
                "solver's aerosol species; set observation metadata 'observation_time'."
            )
        return int(month)

    def _fetch_winning_toa(
        self,
        needed: dict[str, IndexImage],
        *,
        mgrs_tile: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
    ) -> dict[str, np.ndarray]:
        reader = self._reader()

        def _read(image: IndexImage) -> tuple[str, np.ndarray | None]:
            try:
                planes = reader.read(
                    mgrs_tile=mgrs_tile,
                    sensing_token=image.sensing_token,
                    crs=crs,
                    transform=transform,
                    width=width,
                    height=height,
                )
            except Exception:  # noqa: BLE001 - one bad acquisition must not kill the library
                logger.warning(
                    "Live L1C library: failed reading %s; its pixels stay unfilled",
                    image.image_id,
                    exc_info=True,
                )
                return image.image_id, None
            return image.image_id, planes

        fetched: dict[str, np.ndarray] = {}
        images = list(needed.values())
        with ThreadPoolExecutor(max_workers=min(self._max_workers, max(1, len(images)))) as pool:
            for image_id, planes in pool.map(_read, images):
                if planes is not None:
                    fetched[image_id] = planes
        return fetched

    def realizations(
        self,
        observation: ObservationBundle,
        resolution: float,
        bands: Sequence[str],  # noqa: ARG002 - the library declares its own bands
    ) -> list[SurfaceLibraryRealization]:
        from siac.algorithms.grid.assembler import _build_target_template

        index = read_mosaic_index(self._index_path(observation))
        template = _build_target_template(
            observation.bounds, str(observation.crs), float(resolution)
        )
        crs, transform, width, height = _template_grid(template)
        winners = _resample_winners(index.winners, height, width)
        mgrs_tile = resolve_mgrs_tile(observation, self._scene_key)
        day_aod = self._day_aod(index, observation)
        rt_model = self._effective_rt_model(observation)
        sensor_bands = [observation.sensor_config.get_band(name) for name in LIBRARY_BAND_NAMES]

        needed: dict[str, IndexImage] = {}
        winning_by_month: dict[str, list[int]] = {}
        for month_index, month in enumerate(index.months):
            plane = winners[month_index]
            images = index.images_by_month.get(month, ())
            selected = [int(j) for j in np.unique(plane[plane >= 0]) if int(j) < len(images)]
            winning_by_month[month] = selected
            for position in selected:
                needed[images[position].image_id] = images[position]

        if not needed:
            raise RuntimeError(
                f"Mosaic index for {self._index_root} selected no acquisitions over the scene AOI."
            )
        logger.info(
            "Live L1C library: %d month(s), %d winning acquisition(s) to read from GCS tile %s",
            len(index.months),
            len(needed),
            mgrs_tile,
        )
        fetched = self._fetch_winning_toa(
            needed,
            mgrs_tile=mgrs_tile,
            crs=crs,
            transform=transform,
            width=width,
            height=height,
        )

        corrected: list[np.ndarray] = []
        labels: list[str] = []
        for month_index, month in enumerate(index.months):
            selected = winning_by_month[month]
            if not selected:
                continue
            mosaic = self._assemble_month(
                plane=winners[month_index],
                images=index.images_by_month[month],
                selected=selected,
                fetched=fetched,
                day_aod=day_aod,
                shape=(height, width),
            )
            if mosaic is None:
                continue
            toa, aod, angles = mosaic
            corrected.append(
                correct_toa_realization(
                    rt_model=rt_model,
                    sensor_bands=sensor_bands,
                    toa=toa,
                    aod=aod,
                    sza_deg=angles["sza"],
                    saa_deg=angles["saa"],
                    vza_deg=angles["vza"],
                    vaa_deg=angles["vaa"],
                    tcwv=self._tcwv,
                    tco3=self._tco3,
                    elevation_km=self._elevation_km,
                )
            )
            labels.append(month)

        if not corrected:
            raise RuntimeError(
                "Live L1C surface library produced no usable realizations: none of the winning "
                "acquisitions could be read from the public Sentinel-2 bucket."
            )
        kept = _keep_realizations(corrected, self._keep_fraction)
        logger.info(
            "Live L1C library: kept %d/%d realizations (%s), rt_space=%s",
            len(kept),
            len(corrected),
            ",".join(labels[index] for index in kept),
            self.rt_space if self.rt_space is not None else "unmanaged",
        )
        return [
            SurfaceLibraryRealization(
                reflectance=corrected[index],
                band_names=PREDICTOR_BAND_ORDER,
                crs=crs,
                transform=transform,
                label=labels[index],
            )
            for index in kept
        ]

    @staticmethod
    def _assemble_month(
        *,
        plane: np.ndarray,
        images: Sequence[IndexImage],
        selected: Sequence[int],
        fetched: dict[str, np.ndarray],
        day_aod: dict[str, float],
        shape: tuple[int, int],
    ) -> tuple[np.ndarray, np.ndarray, dict[str, float]] | None:
        """Mosaic one month's winning acquisitions through the winner index."""

        height, width = shape
        toa: np.ndarray = np.full(
            (len(LIBRARY_BAND_NAMES), height, width), np.nan, dtype=np.float32
        )
        aod: np.ndarray = np.full((height, width), np.nan, dtype=np.float32)
        angle_planes: dict[str, np.ndarray] = {
            name: np.full((height, width), np.nan, dtype=np.float32)
            for name in ("sza", "saa", "vza", "vaa")
        }
        filled = False
        for position in selected:
            image = images[position]
            mask = plane == position
            planes = fetched.get(image.image_id)
            if planes is not None:
                toa[:, mask] = planes[:, mask]
                filled = True
            aod[mask] = day_aod.get(image.day, _DEFAULT_DAY_AOD)
            for name in angle_planes:
                angle_planes[name][mask] = float(getattr(image, name))
        if not filled or not np.isfinite(toa).any():
            return None
        angles = {name: _plane_median(values) for name, values in angle_planes.items()}
        return toa, aod, angles

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"LiveL1CSurfaceLibrary(index_root={self._index_root!s})"


def _plane_median(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    return float(np.median(finite)) if finite.size else 0.0


def _keep_realizations(realizations: Sequence[np.ndarray], keep_fraction: float) -> list[int]:
    """Drop realizations too sparse to be useful, keeping a workable minimum.

    A mosaic that is mostly cloud-masked contributes noise to the temporal
    median, but dropping down to one or two realizations destroys the temporal
    spread that becomes the prior uncertainty — so the filter yields only when
    enough realizations survive.
    """

    finite = np.array(
        [float(np.isfinite(realization).mean()) for realization in realizations], dtype=np.float64
    )
    keep = [index for index, value in enumerate(finite) if value >= float(keep_fraction)]
    if len(keep) < min(_MIN_REALIZATIONS, len(realizations)):
        return list(range(len(realizations)))
    return keep
