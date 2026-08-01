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
  day's own atmospheric state. A library corrected in one RT space and solved
  against in another injects a systematic surface offset that the solver turns
  into biased AOT — see :mod:`siac.domain.rt_space`.

Every atmospheric and terrain input the correction needs comes from a real
measured source, per winning day:

==============  ========================================================
AOD 550 nm      MAIAC, from the index's own per-day scalars (or the
                :class:`~siac.adapters.atmo.maiac_day_aod.MAIACDayAODProvider`)
Water vapour    retrieved from the acquisition's OWN L1C B8A/B09 band
                ratio (:mod:`siac.algorithms.water_vapour`), so no
                published Level-2A product has to exist or be reachable;
                ``water_vapour_source="l2a"`` reads Sen2Cor's instead
Ozone           CAMS ``gtco3`` for that day at the scene's overpass hour,
                via :class:`~siac.adapters.atmo.cams.CAMSProvider`
Terrain         the configured DEM (Copernicus GLO-30), sampled onto the
                library's own grid by :func:`siac.geo.dem.read_elevation_km`
==============  ========================================================

**None of them has a default.** A missing source raises and names what is
missing. Substituting a plausible constant for a measured atmospheric input is
the failure mode that made a broken aerosol prior invisible for a whole
campaign: runs still "succeed", and the damage only shows up in validation.

The winner index itself is an input (``months``/``winners``/``image_table``/
``day_scalars`` in one ``.npz`` per scene), because building it needs a
cloud-score source outside this package. Reflectance downstream of it is
GEE-free: L1C top-of-atmosphere reflectance is read straight from the anonymous
``gs://gcp-public-data-sentinel-2`` bucket.

The native 6S backend evaluates one atmospheric state per invocation, so the
correction builds coefficient *curves* over an AOD node axis at the mosaic's
median geometry, water vapour and ozone, then interpolates them per pixel — the
same construction the validated converter used, and the same one the solver's
coefficient provider uses.
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
    from datetime import datetime

    import xarray as xr

    from siac.adapters.data.s2_data_source import S2Product
    from siac.domain import SensorBand
    from siac.domain.rt_space import RTSpace
    from siac.runtime.models import ObservationBundle

logger = logging.getLogger(__name__)

__all__ = [
    "LIBRARY_BAND_NAMES",
    "WATER_VAPOUR_BAND_NAMES",
    "CAMSTCO3Source",
    "CIBRWaterVapourReader",
    "GCSL1CTOAReader",
    "IndexImage",
    "L1CBandReader",
    "L1CTOAReader",
    "LiveL1CSurfaceLibrary",
    "MissingLibraryInputError",
    "MosaicIndex",
    "PlanetaryComputerWVPReader",
    "TCO3Source",
    "WVPReader",
    "correct_toa_realization",
    "read_mosaic_index",
    "resolve_mgrs_tile",
]

#: Sentinel-2 band names read from L1C, ordered to match
#: :data:`~siac.adapters.surface_library.PREDICTOR_BAND_ORDER`
#: (coastal, blue, green, red, nir, swir16, swir22).
LIBRARY_BAND_NAMES: tuple[str, ...] = ("B01", "B02", "B03", "B04", "B8A", "B11", "B12")

#: Bands the CIBR water-vapour retrieval needs: the 940 nm absorption band and
#: its continuum neighbour.
WATER_VAPOUR_BAND_NAMES: tuple[str, ...] = ("B8A", "B09")

#: L1C digital-number scale factor (DN -> reflectance).
_L1C_DN_SCALE = 10000.0
#: Sentinel-2 L2A ``WVP`` digital-number scale factor (DN -> cm).
_WVP_DN_SCALE = 1000.0
#: Processing baselines from this number on carry a radiometric DN offset.
_OFFSET_BASELINE = 400
#: The DN offset those baselines apply.
_BASELINE_OFFSET_DN = 1000.0
#: Number of AOD nodes the 6S coefficient curves are evaluated on.
_AOD_NODE_COUNT = 8
#: AOD node-axis bounds (the native 6S aerosol axis is only valid inside these).
_AOD_NODE_MIN = 0.01
_AOD_NODE_MAX = 4.0
#: Physically plausible total ozone column (atm-cm); outside this the CAMS read
#: is wrong rather than merely unusual, so it is rejected instead of accepted.
_TCO3_VALID_RANGE = (0.15, 0.55)
#: Physically plausible total column water vapour (cm).
_TCWV_VALID_RANGE = (0.05, 15.0)
#: Minimum finite fraction a corrected realization must reach to be kept.
_KEEP_FRACTION = 0.6
#: Realizations always kept even if none reaches ``_KEEP_FRACTION``.
_MIN_REALIZATIONS = 4
#: Default parallelism for the per-acquisition reads.
_DEFAULT_MAX_WORKERS = 8
#: Local solar time of the Sentinel-2 descending-node overpass, in hours.
_S2_OVERPASS_LOCAL_HOUR = 10.5

#: Planetary Computer STAC endpoints for the Sentinel-2 L2A water-vapour band.
_PC_STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
_PC_SAS_URL = "https://planetarycomputer.microsoft.com/api/sas/v1/token/sentinel-2-l2a"

#: GDAL settings for windowed ``/vsicurl`` reads of remote imagery.
_GDAL_REMOTE_ENV = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "GDAL_HTTP_MAX_RETRY": "3",
    "GDAL_HTTP_RETRY_DELAY": "2",
}

_MGRS_IN_TEXT = re.compile(r"T(\d{2}[A-Z]{3})")
_MGRS_BARE = re.compile(r"^T?(\d{2}[A-Z]{3})$")


class MissingLibraryInputError(RuntimeError):
    """A measured input the correction needs is unavailable for this scene.

    Raised instead of substituting a plausible constant: a library corrected at
    an invented water vapour, ozone or terrain height is silently wrong, and the
    error only surfaces as biased AOT much later.
    """


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

    @property
    def raa(self) -> float:
        """Relative azimuth folded into ``[0, 180]`` (the LUT's convention)."""

        from siac.algorithms.water_vapour import relative_azimuth_deg

        return float(relative_azimuth_deg(self.vaa, self.saa))


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


def _finite_median(values: np.ndarray) -> float | None:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None
    return float(np.median(finite))


# --------------------------------------------------------------------------- #
# remote reads
# --------------------------------------------------------------------------- #
def _warp_read(
    url: str,
    *,
    crs: str,
    transform: tuple[float, float, float, float, float, float],
    width: int,
    height: int,
    resampling: str = "average",
) -> np.ndarray:
    """Windowed-read + warp one remote raster straight onto the target grid.

    A :class:`~rasterio.vrt.WarpedVRT` is used rather than
    :mod:`siac.geo.reprojection` because the source is a remote 110 km scene and
    only the AOI window is wanted; the VRT reads solely the overview blocks the
    target grid touches. ``average`` matches the mean-aggregated pyramids the
    reference mosaics were built from.
    """

    import rasterio
    from rasterio.transform import Affine
    from rasterio.vrt import WarpedVRT

    from siac.geo.reprojection import RESAMPLING_METHODS

    with (
        rasterio.open(url) as src,
        WarpedVRT(
            src,
            crs=crs,
            transform=Affine(*transform),
            width=width,
            height=height,
            resampling=RESAMPLING_METHODS[resampling],
        ) as vrt,
    ):
        dn = vrt.read(1).astype(np.float32)
    dn[dn == 0] = np.nan
    return cast("np.ndarray", dn)


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


class L1CBandReader(Protocol):
    """Reads named L1C bands of one acquisition onto a grid."""

    def read_bands(
        self,
        bands: Sequence[str],
        *,
        mgrs_tile: str,
        sensing_token: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
    ) -> dict[str, np.ndarray] | None:
        """Return TOA reflectance per requested band, or ``None`` if unavailable."""


class WVPReader(Protocol):
    """Resolves one acquisition's total column water vapour over the AOI."""

    def read(
        self,
        *,
        mgrs_tile: str,
        sensing_token: str,
        day: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
        sza_deg: float,
        vza_deg: float,
        raa_deg: float,
        elevation_km: np.ndarray | float,
    ) -> np.ndarray | None:
        """Return a ``(y, x)`` total column water vapour field in cm, or ``None``.

        The geometry and terrain arguments let a reader retrieve the column from
        the scene itself; a reader that reads a published product ignores them.
        """


class TCO3Source(Protocol):
    """Resolves one day's total ozone column over the scene, in atm-cm."""

    def __call__(self, day: str, *, lon: float, lat: float) -> float: ...


def _baseline_offset(processing_baseline: str | None) -> float:
    """Radiometric DN offset implied by a SAFE processing baseline."""

    text = str(processing_baseline or "").lstrip("Nn")
    if not text.isdigit():
        return 0.0
    return _BASELINE_OFFSET_DN if int(text) >= _OFFSET_BASELINE else 0.0


def _band_jp2_urls(product: S2Product, bands: Sequence[str]) -> dict[str, str]:
    """Map each requested band to its ``/vsicurl`` JP2 URL inside a SAFE product."""

    prefix = gcs._resolve_safe_prefix(product)
    objects = gcs._list_objects_under(prefix)
    urls: dict[str, str] = {}
    for band in bands:
        for item in objects:
            name = str(item.get("name", ""))
            if name.endswith(f"_{band}.jp2") and "/IMG_DATA/" in name:
                urls[band] = "/vsicurl/" + gcs._object_download_url(name)
                break
    return urls


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
        # The listing walks every SAFE prefix on the tile (thousands of objects),
        # so the lock is held across the fetch: releasing it first lets all the
        # reader threads issue the same listing concurrently.
        with self._lock:
            cached = self._by_tile.get(mgrs_tile)
            if cached is not None:
                return cached
            products = self._backend.search(S2Query(mgrs_tile=mgrs_tile, processing_level="L1C"))
            by_token: dict[str, S2Product] = {}
            for product in sorted(products, key=lambda item: str(item.product_id)):
                by_token.setdefault(product.sensing_date.strftime("%Y%m%dT%H%M%S"), product)
            self._by_tile[mgrs_tile] = by_token
        logger.info(
            "Live L1C library: %d L1C products listed for tile %s", len(by_token), mgrs_tile
        )
        return by_token

    def read_bands(
        self,
        bands: Sequence[str],
        *,
        mgrs_tile: str,
        sensing_token: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
    ) -> dict[str, np.ndarray] | None:
        """Read named L1C bands onto the target grid, in reflectance.

        Every band gets the same radiometric handling, including the
        ``RADIO_ADD_OFFSET`` that processing baselines N0400 and later apply.
        Missing that offset on one band silently biases anything derived from a
        band *ratio* — which is why v1's own TCWV call ended up commented out.
        """

        import rasterio

        product = self._products_for_tile(mgrs_tile).get(sensing_token)
        if product is None:
            logger.info(
                "Live L1C library: no L1C product on GCS for %s at %s", mgrs_tile, sensing_token
            )
            return None
        urls = _band_jp2_urls(product, bands)
        if not urls:
            logger.warning("Live L1C library: no JP2 bands under %s", product.product_id)
            return None
        offset = _baseline_offset(product.processing_baseline)
        planes: dict[str, np.ndarray] = {}
        with rasterio.Env(**_GDAL_REMOTE_ENV, CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".jp2"):
            for band in bands:
                url = urls.get(band)
                if url is None:
                    continue
                dn = _warp_read(
                    url,
                    crs=crs,
                    transform=transform,
                    width=width,
                    height=height,
                    resampling=self._resampling,
                )
                planes[band] = ((dn - offset) / _L1C_DN_SCALE).astype(np.float32)
        return planes

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
        planes = self.read_bands(
            LIBRARY_BAND_NAMES,
            mgrs_tile=mgrs_tile,
            sensing_token=sensing_token,
            crs=crs,
            transform=transform,
            width=width,
            height=height,
        )
        if planes is None:
            return None
        stack: np.ndarray = np.full(
            (len(LIBRARY_BAND_NAMES), height, width), np.nan, dtype=np.float32
        )
        for index, band in enumerate(LIBRARY_BAND_NAMES):
            if band in planes:
                stack[index] = planes[band]
        return stack


class PlanetaryComputerWVPReader:
    """Reads the Sentinel-2 L2A ``WVP`` band from the Planetary Computer.

    Water vapour is a *measured* input to the correction, and the matching L2A
    product carries Sen2Cor's own retrieval of it for the very acquisition being
    corrected — which is why the reference library used it rather than a model
    field. Items are matched to the index's acquisitions by sensing token, and
    the highest processing baseline for that acquisition wins.
    """

    def __init__(self, session: Any | None = None, *, resampling: str = "average") -> None:
        self._session = session
        self._resampling = str(resampling)
        self._token: str | None = None
        self._items: dict[tuple[str, str], dict[str, Any] | None] = {}
        self._lock = threading.Lock()

    def _http(self) -> Any:
        if self._session is None:
            from siac.adapters._http import make_session

            self._session = make_session()
        return self._session

    def _sas_token(self) -> str:
        with self._lock:
            if self._token is not None:
                return self._token
        response = self._http().get(_PC_SAS_URL, timeout=30)
        response.raise_for_status()
        token = str(response.json()["token"])
        with self._lock:
            self._token = token
        return token

    def _find_item(self, mgrs_tile: str, sensing_token: str, day: str) -> dict[str, Any] | None:
        key = (mgrs_tile, sensing_token)
        with self._lock:
            if key in self._items:
                return self._items[key]
        payload = {
            "collections": ["sentinel-2-l2a"],
            "datetime": f"{day}T00:00:00Z/{day}T23:59:59Z",
            "query": {"s2:mgrs_tile": {"eq": mgrs_tile}},
            "limit": 100,
        }
        response = self._http().post(f"{_PC_STAC_URL}/search", json=payload, timeout=60)
        response.raise_for_status()
        matches = [
            feature
            for feature in response.json().get("features", [])
            if sensing_token in str(feature.get("id", ""))
        ]
        item = (
            max(
                matches,
                key=lambda feature: str(feature["properties"].get("s2:processing_baseline", "0")),
            )
            if matches
            else None
        )
        with self._lock:
            self._items[key] = item
        return item

    def read(
        self,
        *,
        mgrs_tile: str,
        sensing_token: str,
        day: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
        sza_deg: float = 0.0,  # noqa: ARG002 - a published product needs no geometry
        vza_deg: float = 0.0,  # noqa: ARG002
        raa_deg: float = 0.0,  # noqa: ARG002
        elevation_km: np.ndarray | float = 0.0,  # noqa: ARG002
    ) -> np.ndarray | None:
        import rasterio

        item = self._find_item(mgrs_tile, sensing_token, day)
        if item is None or "WVP" not in item.get("assets", {}):
            logger.info("Live L1C library: no L2A WVP asset for %s at %s", mgrs_tile, sensing_token)
            return None
        href = f"{item['assets']['WVP']['href']}?{self._sas_token()}"
        with rasterio.Env(**_GDAL_REMOTE_ENV):
            dn = _warp_read(
                href,
                crs=crs,
                transform=transform,
                width=width,
                height=height,
                resampling=self._resampling,
            )
        return cast("np.ndarray", dn / np.float32(_WVP_DN_SCALE))


class CIBRWaterVapourReader:
    """Retrieves water vapour from the acquisition's own L1C TOA.

    The Continuum Interpolated Band Ratio (B8A continuum / B09 940 nm
    absorption) sounds the column directly from the imagery already being
    fetched, so no published Level-2A product has to exist, be reachable, or
    have been processed with a baseline whose radiometry we understand. See
    :mod:`siac.algorithms.water_vapour` for the retrieval itself.

    Both bands are read through the same :class:`GCSL1CTOAReader`, which applies
    one radiometric offset to all of them — a band *ratio* is exactly where a
    per-band offset mismatch would do its damage.
    """

    def __init__(self, band_reader: L1CBandReader) -> None:
        self._band_reader = band_reader

    def read(
        self,
        *,
        mgrs_tile: str,
        sensing_token: str,
        day: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
        sza_deg: float,
        vza_deg: float,
        raa_deg: float,
        elevation_km: np.ndarray | float,
    ) -> np.ndarray | None:
        from siac.algorithms.water_vapour import retrieve_water_vapour

        planes = self._band_reader.read_bands(
            WATER_VAPOUR_BAND_NAMES,
            mgrs_tile=mgrs_tile,
            sensing_token=sensing_token,
            crs=crs,
            transform=transform,
            width=width,
            height=height,
        )
        missing = (
            list(WATER_VAPOUR_BAND_NAMES)
            if planes is None
            else [band for band in WATER_VAPOUR_BAND_NAMES if band not in planes]
        )
        if planes is None or missing:
            logger.info(
                "Live L1C library: %s at %s has no %s for the CIBR water-vapour retrieval",
                mgrs_tile,
                sensing_token,
                ",".join(missing),
            )
            return None
        result = retrieve_water_vapour(
            toa_b09=planes["B09"],
            toa_b8a=planes["B8A"],
            sza_deg=sza_deg,
            vza_deg=vza_deg,
            raa_deg=raa_deg,
            elevation_km=elevation_km,
        )
        if not result.valid.any():
            logger.warning(
                "Live L1C library: CIBR retrieved no valid water vapour for %s (%s)",
                sensing_token,
                day,
            )
            return None
        logger.info(
            "Live L1C library: CIBR water vapour for %s (%s): median %.3f cm, %.1f%% filled "
            "from this acquisition's own median",
            sensing_token,
            day,
            float(np.nanmedian(result.water_vapour_cm[result.valid])),
            100.0 * result.masked_fraction,
        )
        return result.water_vapour_cm


class CAMSTCO3Source:
    """Per-day total ozone column at the scene, from the CAMS archive.

    :class:`~siac.adapters.atmo.cams.CAMSProvider` supplies both halves of the
    read that must not be reinvented: it resolves the ``{base}/YYYY-MM-DD.nc``
    archive layout (including the mirrored naming variants) and publishes the
    ``gtco3`` kg m^-2 -> atm-cm scale. The selection on top of it is the scene's
    nearest grid point at the overpass hour, which is what the reference library
    sampled.

    The day's file is required. ``get_prior`` is deliberately *not* used here:
    it degrades to a default 0.3 atm-cm state when the day is missing, which is
    exactly the invisible substitution this library refuses to make.
    """

    #: CAMS variable holding the total ozone column.
    _VARIABLE = "gtco3"

    def __init__(self, data_path: str | Path, *, cache_dir: str | Path | None = None) -> None:
        from siac.adapters.atmo.cams import CAMSProvider

        self._provider = CAMSProvider(
            data_path,
            temporal_interp=False,
            download_missing=False,
            cache_dir=cache_dir,
        )
        self._data_path = data_path
        self._cache: dict[str, float] = {}
        self._lock = threading.Lock()

    def __call__(self, day: str, *, lon: float, lat: float) -> float:
        with self._lock:
            cached = self._cache.get(day)
        if cached is not None:
            return cached

        from datetime import datetime

        from siac.adapters.atmo.cams import CAMSProvider

        # Sentinel-2 crosses the descending node at ~10:30 local solar time; the
        # reference library sampled CAMS at that hour, not at midnight.
        hour = int(round(_S2_OVERPASS_LOCAL_HOUR - float(lon) / 15.0)) % 24
        when = datetime.strptime(f"{day}T{hour:02d}", "%Y-%m-%dT%H")
        dataset = self._provider._load_cams_data(when)
        if dataset is None:
            raise MissingLibraryInputError(
                f"No CAMS ozone for {day} under {self._data_path!s}. The live L1C surface "
                "library corrects each winning day at its own measured ozone column and will "
                "not substitute a default; stage the day's CAMS file, or point "
                "providers.atmo.data_path at an archive that covers the library's months."
            )
        try:
            raw = self._select_point(dataset, day=day, hour=hour, lon=lon, lat=lat)
        finally:
            close = getattr(dataset, "close", None)
            if close is not None:
                close()
        value = float(raw) * float(CAMSProvider._RAW_CAMS_TO_SIAC_SCALE[self._VARIABLE])
        low, high = _TCO3_VALID_RANGE
        if not np.isfinite(value) or not low <= value <= high:
            raise MissingLibraryInputError(
                f"CAMS ozone for {day} at the scene is {value!r} atm-cm, outside the plausible "
                f"range {low}-{high}. The library will not correct at an implausible ozone "
                "column; check the CAMS archive for that day."
            )
        with self._lock:
            self._cache[day] = value
        return value

    @classmethod
    def _select_point(cls, dataset: Any, *, day: str, hour: int, lon: float, lat: float) -> float:
        """Nearest CAMS grid point at ``hour``, in the file's raw kg m^-2."""

        if cls._VARIABLE not in dataset:
            raise MissingLibraryInputError(
                f"The CAMS file for {day} carries no {cls._VARIABLE!r} variable, so the scene's "
                "ozone column cannot be read."
            )
        field = dataset[cls._VARIABLE]
        if "forecast_reference_time" in field.dims:
            field = field.squeeze("forecast_reference_time", drop=True)
        if "forecast_period" in field.dims:
            selector: dict[str, Any] = {"forecast_period": float(hour)}
        else:
            selector = {"time": np.datetime64(f"{day}T{hour:02d}:00")}
        # CAMS grids are published on both -180..180 and 0..360 longitudes.
        longitudes = np.asarray(field["longitude"].values, dtype=np.float64)
        point_lon = float(lon)
        if longitudes.max() > 180.0 and point_lon < 0.0:
            point_lon += 360.0
        return float(
            field.sel(latitude=float(lat), longitude=point_lon, **selector, method="nearest").item()
        )


# --------------------------------------------------------------------------- #
# correction
# --------------------------------------------------------------------------- #
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
    mosaic's range and interpolated per pixel. Geometry, water vapour, ozone and
    terrain are scalar: the native 6S backend runs one invocation per state, and
    discretising those smoothly varying fields there is what the validated
    converter did.
    """

    from siac.runtime.models import AtmosphericState, GeometryAngles

    toa_values = np.asarray(toa, dtype=np.float64)
    aod_values = np.asarray(aod, dtype=np.float64)
    if toa_values.ndim != 3 or toa_values.shape[0] != len(sensor_bands):
        raise ValueError(
            f"TOA must be (band, y, x) with {len(sensor_bands)} bands; got {toa_values.shape}."
        )

    # MAIAC reports 0.000 on very clean days, and such a day can win every pixel
    # of a realization (La Parguera 2018-06). A measured zero is an aerosol
    # state — the interpolation clamps it to the lowest node — so only NaN,
    # the absence of a measurement, refuses.
    finite_aod = np.isfinite(aod_values) & (aod_values >= 0.0)
    if not finite_aod.any():
        raise MissingLibraryInputError(
            "No finite per-pixel AOD for this realization; the correction has no aerosol state."
        )
    fill = float(np.median(aod_values[finite_aod]))
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


def _require_matching_winner_grid(
    winners: np.ndarray, height: int, width: int, index_root: Path
) -> None:
    """Reject a winner index that was not built on the library's own grid.

    The index file carries no CRS or transform, so nothing can verify that a
    differently shaped index covers the same ground. Resampling it on the
    assumption that it does would silently mis-attribute pixels to the wrong
    acquisition — and therefore to the wrong day's atmospheric state.
    """

    if winners.shape[1:] == (height, width):
        return
    source_height, source_width = winners.shape[1:]
    raise ValueError(
        f"Mosaic index under {index_root} was built on a {source_height}x{source_width} grid but "
        f"the library grid is {height}x{width}. The index carries no CRS or transform, so the two "
        "cannot be co-registered; rebuild the index on this scene's AOI and resolution, or run "
        "the retrieval on the AOI and resolution the index was built for."
    )


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


@dataclass(frozen=True)
class _MonthMosaic:
    """One month's mosaicked TOA plus the measured state it is corrected at."""

    toa: np.ndarray
    aod: np.ndarray
    tcwv: float
    tco3: float
    angles: dict[str, float]


class LiveL1CSurfaceLibrary:
    """Surface library assembled live from L1C, corrected in the solver's RT space.

    For each month in the scene's winner index this fetches only the winning
    acquisitions' L1C top-of-atmosphere reflectance from the anonymous Google
    Cloud Sentinel-2 bucket, mosaics them through the index, and corrects the
    mosaic with ``rt_model`` at that month's measured aerosol, water-vapour,
    ozone and terrain state. One
    :class:`~siac.adapters.surface_library.SurfaceLibraryRealization` per month
    comes out, matching the prepared store the same construction writes offline.
    """

    def __init__(
        self,
        index_root: str | Path,
        *,
        rt_model: Any,
        dem_path: str | Path | None,
        cams_data_path: str | Path | None,
        solver_config: Any | None = None,
        scene_key: str | None = None,
        cams_cache_dir: str | Path | None = None,
        maiac_day_aod: Callable[
            [tuple[float, float, float, float], str, Sequence[tuple[int, int]]], dict[str, float]
        ]
        | None = None,
        toa_reader: L1CTOAReader | None = None,
        wvp_reader: WVPReader | None = None,
        water_vapour_source: str = "cibr",
        tco3_source: TCO3Source | None = None,
        keep_fraction: float = _KEEP_FRACTION,
        max_workers: int = _DEFAULT_MAX_WORKERS,
    ) -> None:
        self._index_root = Path(index_root).expanduser()
        self._rt_model = rt_model
        self._dem_path = dem_path
        self._cams_data_path = cams_data_path
        self._cams_cache_dir = cams_cache_dir
        self._solver_config = solver_config
        self._scene_key = scene_key
        self._maiac_day_aod = maiac_day_aod
        self._toa_reader = toa_reader
        self._wvp_reader = wvp_reader
        self._water_vapour_source = str(water_vapour_source)
        self._tco3_source = tco3_source
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

    # -- inputs ------------------------------------------------------------- #
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

    def _water_vapour_reader(self) -> WVPReader:
        if self._wvp_reader is not None:
            return self._wvp_reader
        source = self._water_vapour_source
        if source == "cibr":
            reader = self._reader()
            if not hasattr(reader, "read_bands"):
                raise MissingLibraryInputError(
                    "water_vapour_source='cibr' retrieves the column from the scene's own L1C "
                    "B8A/B09, so the TOA reader must be able to read named bands."
                )
            self._wvp_reader = CIBRWaterVapourReader(cast("L1CBandReader", reader))
        elif source == "l2a":
            self._wvp_reader = PlanetaryComputerWVPReader()
        else:
            raise ValueError(f"Unknown water_vapour_source {source!r}; expected 'cibr' or 'l2a'.")
        return self._wvp_reader

    def _ozone_source(self) -> TCO3Source:
        if self._tco3_source is None:
            if not self._cams_data_path:
                raise MissingLibraryInputError(
                    "The live L1C surface library corrects each winning day at its own measured "
                    "ozone column, and no CAMS source is configured. Set providers.atmo.data_path "
                    "to a CAMS archive (e.g. the group mirror laid out as "
                    "'<base>/YYYY-MM-DD.nc') covering the library's months."
                )
            self._tco3_source = CAMSTCO3Source(self._cams_data_path, cache_dir=self._cams_cache_dir)
        return self._tco3_source

    def _elevation_field(self, template: xr.DataArray) -> np.ndarray:
        """Per-pixel terrain height (km) sampled from the configured DEM.

        The field, not just its median, because the water-vapour retrieval keys
        on altitude per pixel; the 6S correction takes the median from it.
        """

        from siac.geo.dem import read_elevation_km, use_sea_level_elevation

        dem = None if self._dem_path is None else str(self._dem_path)
        if use_sea_level_elevation(dem):
            raise MissingLibraryInputError(
                "The live L1C surface library corrects at the scene's real terrain height and no "
                "DEM is configured (paths.dem). Set paths.dem to a DEM such as "
                "siac.config.public.COPERNICUS_GLO30_DEM_VRT; a sea-level placeholder would "
                "under-attribute AOD over high terrain."
            )
        elevation: np.ndarray = np.asarray(
            read_elevation_km(template, dem).values, dtype=np.float64
        )
        value = _finite_median(elevation)
        if value is None:
            raise MissingLibraryInputError(
                f"The DEM at {dem} yielded no finite elevation over the scene AOI."
            )
        logger.info("Live L1C library: terrain elevation %.3f km from %s", value, dem)
        return elevation

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
            resolved = dict(source(observation.bounds, str(observation.crs), periods))
        except Exception as exc:
            raise MissingLibraryInputError(
                "The live L1C surface library corrects each winning day at its own measured AOD, "
                "and the mosaic index carries no per-day scalars, so MAIAC was queried and "
                f"failed ({type(exc).__name__}: {exc}). Rebuild the index with day scalars, or "
                "make the MAIAC source reachable."
            ) from exc
        if not resolved:
            raise MissingLibraryInputError(
                "Neither the mosaic index nor MAIAC provided a per-day AOD for this scene; the "
                "library will not correct at an assumed aerosol loading."
            )
        return resolved

    def _effective_rt_model(self, observation: ObservationBundle, lon: float, lat: float) -> Any:
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
        observation_time: datetime | None = getattr(self._rt_model, "observation_time", None)
        if observation_time is None:
            observation_time = (observation.metadata or {}).get("observation_time")
        month = getattr(observation_time, "month", None)
        if month is None:
            raise ValueError(
                "The live L1C surface library needs the scene's observation time to resolve the "
                "solver's aerosol species; set observation metadata 'observation_time'."
            )
        return int(month)

    # -- acquisition -------------------------------------------------------- #
    def _fetch_winning_inputs(
        self,
        needed: dict[str, IndexImage],
        *,
        mgrs_tile: str,
        crs: str,
        transform: tuple[float, float, float, float, float, float],
        width: int,
        height: int,
        elevation_km: np.ndarray,
    ) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        """Read each winning acquisition's L1C TOA and its water vapour."""

        toa_reader = self._reader()
        wvp_reader = self._water_vapour_reader()
        grid: dict[str, Any] = {
            "crs": crs,
            "transform": transform,
            "width": width,
            "height": height,
        }

        def _read(image: IndexImage) -> tuple[str, np.ndarray | None, np.ndarray | None]:
            planes: np.ndarray | None = None
            water: np.ndarray | None = None
            try:
                planes = toa_reader.read(
                    mgrs_tile=mgrs_tile, sensing_token=image.sensing_token, **grid
                )
            except Exception:  # noqa: BLE001 - one bad acquisition must not kill the library
                logger.warning(
                    "Live L1C library: failed reading L1C TOA for %s; its pixels stay unfilled",
                    image.image_id,
                    exc_info=True,
                )
            if planes is not None:
                try:
                    water = wvp_reader.read(
                        mgrs_tile=mgrs_tile,
                        sensing_token=image.sensing_token,
                        day=image.day,
                        sza_deg=image.sza,
                        vza_deg=image.vza,
                        raa_deg=image.raa,
                        elevation_km=elevation_km,
                        **grid,
                    )
                except Exception:  # noqa: BLE001 - reported as a missing input below
                    logger.warning(
                        "Live L1C library: failed resolving water vapour for %s",
                        image.image_id,
                        exc_info=True,
                    )
            return image.image_id, planes, water

        toa_by_id: dict[str, np.ndarray] = {}
        wvp_by_id: dict[str, np.ndarray] = {}
        images = list(needed.values())
        with ThreadPoolExecutor(max_workers=min(self._max_workers, max(1, len(images)))) as pool:
            for image_id, planes, water in pool.map(_read, images):
                if planes is not None:
                    toa_by_id[image_id] = planes
                if water is not None:
                    wvp_by_id[image_id] = water
        return toa_by_id, wvp_by_id

    # -- assembly ----------------------------------------------------------- #
    def realizations(
        self,
        observation: ObservationBundle,
        resolution: float,
        bands: Sequence[str],  # noqa: ARG002 - the library declares its own bands
    ) -> list[SurfaceLibraryRealization]:
        from siac.algorithms.grid.assembler import _build_target_template
        from siac.geo.reprojection import transform_bounds

        index = read_mosaic_index(self._index_path(observation))
        template = _build_target_template(
            observation.bounds, str(observation.crs), float(resolution)
        )
        crs, transform, width, height = _template_grid(template)
        _require_matching_winner_grid(index.winners, height, width, self._index_root)
        winners = index.winners

        west, south, east, north = transform_bounds(
            observation.bounds, str(observation.crs), "EPSG:4326"
        )
        lon = 0.5 * (float(west) + float(east))
        lat = 0.5 * (float(south) + float(north))

        mgrs_tile = resolve_mgrs_tile(observation, self._scene_key)
        day_aod = self._day_aod(index, observation)
        elevation_field = self._elevation_field(template)
        elevation_km = float(np.median(elevation_field[np.isfinite(elevation_field)]))
        ozone = self._ozone_source()
        rt_model = self._effective_rt_model(observation, lon, lat)
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
            raise MissingLibraryInputError(
                f"Mosaic index for {self._index_root} selected no acquisitions over the scene AOI."
            )
        logger.info(
            "Live L1C library: %d month(s), %d winning acquisition(s) to read for tile %s",
            len(index.months),
            len(needed),
            mgrs_tile,
        )
        toa_by_id, wvp_by_id = self._fetch_winning_inputs(
            needed,
            mgrs_tile=mgrs_tile,
            crs=crs,
            transform=transform,
            width=width,
            height=height,
            elevation_km=elevation_field,
        )

        corrected: list[np.ndarray] = []
        labels: list[str] = []
        for month_index, month in enumerate(index.months):
            selected = winning_by_month[month]
            if not selected:
                continue
            mosaic = self._assemble_month(
                month=month,
                plane=winners[month_index],
                images=index.images_by_month[month],
                selected=selected,
                toa_by_id=toa_by_id,
                wvp_by_id=wvp_by_id,
                day_aod=day_aod,
                ozone=ozone,
                lon=lon,
                lat=lat,
                shape=(height, width),
            )
            if mosaic is None:
                continue
            corrected.append(
                correct_toa_realization(
                    rt_model=rt_model,
                    sensor_bands=sensor_bands,
                    toa=mosaic.toa,
                    aod=mosaic.aod,
                    sza_deg=mosaic.angles["sza"],
                    saa_deg=mosaic.angles["saa"],
                    vza_deg=mosaic.angles["vza"],
                    vaa_deg=mosaic.angles["vaa"],
                    tcwv=mosaic.tcwv,
                    tco3=mosaic.tco3,
                    elevation_km=elevation_km,
                )
            )
            labels.append(month)

        if not corrected:
            raise MissingLibraryInputError(
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
        month: str,
        plane: np.ndarray,
        images: Sequence[IndexImage],
        selected: Sequence[int],
        toa_by_id: dict[str, np.ndarray],
        wvp_by_id: dict[str, np.ndarray],
        day_aod: dict[str, float],
        ozone: TCO3Source,
        lon: float,
        lat: float,
        shape: tuple[int, int],
    ) -> _MonthMosaic | None:
        """Mosaic one month's winning acquisitions through the winner index.

        Every provenance plane is filled from the acquisition that won each
        pixel, so the month's aerosol, water-vapour and ozone state is that of
        the acquisitions actually in the mosaic.
        """

        height, width = shape
        toa: np.ndarray = np.full(
            (len(LIBRARY_BAND_NAMES), height, width), np.nan, dtype=np.float32
        )
        aod: np.ndarray = np.full((height, width), np.nan, dtype=np.float32)
        tcwv: np.ndarray = np.full((height, width), np.nan, dtype=np.float32)
        tco3: np.ndarray = np.full((height, width), np.nan, dtype=np.float32)
        angle_planes: dict[str, np.ndarray] = {
            name: np.full((height, width), np.nan, dtype=np.float32)
            for name in ("sza", "saa", "vza", "vaa")
        }
        filled = False
        for position in selected:
            image = images[position]
            mask = plane == position
            planes = toa_by_id.get(image.image_id)
            if planes is None:
                continue
            toa[:, mask] = planes[:, mask]
            filled = True
            if image.day not in day_aod:
                raise MissingLibraryInputError(
                    f"No measured AOD for {image.day}, which wins pixels in {month}. The mosaic "
                    "index scored no scalar for that day; rebuild the index so every winning day "
                    "carries its MAIAC AOD."
                )
            aod[mask] = day_aod[image.day]
            water = wvp_by_id.get(image.image_id)
            if water is None:
                # The month's other acquisitions still carry measured water
                # vapour, and the correction consumes the month's median, so a
                # single L2A gap leaves those pixels unscored rather than
                # inventing a column for them. A month with no measured water
                # vapour at all is rejected below.
                logger.warning(
                    "Live L1C library: no Sentinel-2 L2A water vapour for %s (%s) in %s; its "
                    "pixels do not contribute to the month's measured column.",
                    image.image_id,
                    image.day,
                    month,
                )
            else:
                tcwv[mask] = water[mask]
            tco3[mask] = ozone(image.day, lon=lon, lat=lat)
            for name in angle_planes:
                angle_planes[name][mask] = float(getattr(image, name))
        if not filled or not np.isfinite(toa).any():
            return None

        tcwv_value = _finite_median(tcwv)
        low, high = _TCWV_VALID_RANGE
        if tcwv_value is None:
            raise MissingLibraryInputError(
                f"No Sentinel-2 L2A water vapour for any winning acquisition in {month}. The "
                "library corrects at the acquisitions' own measured water-vapour column and will "
                "not substitute a default; check the Planetary Computer L2A coverage for those "
                "acquisitions."
            )
        if not low <= tcwv_value <= high:
            raise MissingLibraryInputError(
                f"Sentinel-2 L2A water vapour over the scene for {month} reduced to "
                f"{tcwv_value:.3f} cm, outside the plausible range {low}-{high}. The library will "
                "not correct at an implausible water-vapour column."
            )
        tco3_value = _finite_median(tco3)
        if tco3_value is None:
            raise MissingLibraryInputError(
                f"No CAMS ozone resolved for any winning day in {month}."
            )
        angles = {name: _plane_median(values, month, name) for name, values in angle_planes.items()}
        return _MonthMosaic(toa=toa, aod=aod, tcwv=tcwv_value, tco3=tco3_value, angles=angles)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"LiveL1CSurfaceLibrary(index_root={self._index_root!s})"


def _plane_median(values: np.ndarray, month: str, name: str) -> float:
    value = _finite_median(values)
    if value is None:
        raise MissingLibraryInputError(
            f"No {name} geometry for any winning acquisition in {month}."
        )
    return value


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
