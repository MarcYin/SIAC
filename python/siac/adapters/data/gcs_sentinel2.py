"""
Google Cloud Storage backend for Sentinel-2.

Uses the public ``gs://gcp-public-data-sentinel-2`` bucket (no auth).
Supports MGRS-based listing only (no full catalog search).
See PLANS_S2.md §2, Phase 3.
"""

from __future__ import annotations

import json
import logging
import re
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from typing import Any, cast

from siac.adapters.data.s2_data_source import S2Product, S2Query
from siac.errors import DataNotFoundError

logger = logging.getLogger(__name__)

# Constants
GCS_BUCKET = "gcp-public-data-sentinel-2"
GCS_API_BASE = f"https://storage.googleapis.com/storage/v1/b/{GCS_BUCKET}/o"
GCS_DOWNLOAD_BASE = f"https://storage.googleapis.com/{GCS_BUCKET}/"
GCS_MAX_DOWNLOAD_WORKERS = 8
GCS_DOWNLOAD_MAX_RETRIES = 3
GCS_DOWNLOAD_RETRY_BACKOFF_SEC = 1.0

_RE_TILE = re.compile(r"_T(\d{2}[A-Z]{3})_")
_RE_SENSING = re.compile(r"_([0-9]{8})T([0-9]{6})_")
_RE_BASELINE = re.compile(r"_N(\d{4})_")
_RE_ORBIT = re.compile(r"_R(\d{3})_")
_RE_SAT = re.compile(r"^(S2[A-Z])_")


def _normalize_product_id(product_id: str) -> str:
    product = product_id.strip()
    if product.endswith(".SAFE"):
        product = product[:-5]
    return product


def _normalize_mgrs_tile(tile: str) -> str:
    value = tile.strip().upper().lstrip("T")
    if not re.match(r"^\d{2}[A-Z]{3}$", value):
        raise ValueError(f"Invalid MGRS tile: {tile!r}. Expected like '31UDQ' or 'T31UDQ'.")
    return value


def _mgrs_tile_to_prefix(tile: str) -> str:
    return f"tiles/{tile[:2]}/{tile[2:3]}/{tile[3:]}/"


def _safe_prefix_from_product_id(product_id: str) -> str:
    match = _RE_TILE.search(product_id)
    if match is None:
        raise ValueError(f"Cannot infer MGRS tile from product_id {product_id!r}.")
    tile = _normalize_mgrs_tile(match.group(1))
    return f"{_mgrs_tile_to_prefix(tile)}{product_id}.SAFE/"


def _http_json(url: str, timeout: int = 60) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"User-Agent": "siac-gcs-s2/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = json.load(resp)
    if not isinstance(payload, dict):
        raise ValueError(f"Unexpected JSON payload from {url!r}: expected an object")
    return cast("dict[str, Any]", payload)


def _list_api(
    *,
    prefix: str,
    delimiter: str | None = None,
    page_token: str | None = None,
    max_results: int = 5000,
) -> dict[str, Any]:
    params: dict[str, str] = {"prefix": prefix, "maxResults": str(max_results)}
    if delimiter is not None:
        params["delimiter"] = delimiter
    if page_token is not None:
        params["pageToken"] = page_token
    url = f"{GCS_API_BASE}?{urllib.parse.urlencode(params)}"
    return _http_json(url)


def _list_safe_prefixes(tile_prefix: str) -> list[str]:
    out: list[str] = []
    page_token: str | None = None
    while True:
        payload = _list_api(prefix=tile_prefix, delimiter="/", page_token=page_token)
        for value in payload.get("prefixes", []) or []:
            if isinstance(value, str) and value.endswith(".SAFE/"):
                out.append(value)
        page_token = payload.get("nextPageToken")
        if page_token is None:
            break
    return out


def _list_objects_under(prefix: str) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    page_token: str | None = None
    while True:
        payload = _list_api(prefix=prefix, page_token=page_token)
        for item in payload.get("items", []) or []:
            if isinstance(item, dict) and isinstance(item.get("name"), str):
                out.append(item)
        page_token = payload.get("nextPageToken")
        if page_token is None:
            break
    return out


def _prefix_exists(prefix: str) -> bool:
    page_token: str | None = None
    while True:
        payload = _list_api(prefix=prefix, page_token=page_token, max_results=1)
        if (payload.get("items") or []):
            return True
        page_token = payload.get("nextPageToken")
        if page_token is None:
            return False


def _parse_sensing_datetime(product_id: str) -> datetime:
    match = _RE_SENSING.search(product_id)
    if match is None:
        raise ValueError(f"Cannot parse sensing datetime from product_id {product_id!r}.")
    return datetime.strptime(match.group(1) + match.group(2), "%Y%m%d%H%M%S")


def _parse_baseline(product_id: str) -> str:
    match = _RE_BASELINE.search(product_id)
    if match is None:
        return "N0000"
    return f"N{match.group(1)}"


def _parse_orbit(product_id: str) -> int:
    match = _RE_ORBIT.search(product_id)
    if match is None:
        return 0
    return int(match.group(1))


def _parse_satellite(product_id: str) -> str:
    match = _RE_SAT.search(product_id)
    if match is None:
        return "S2"
    return match.group(1)


def _product_from_safe_prefix(safe_prefix: str) -> S2Product:
    clean = safe_prefix.rstrip("/")
    safe_name = clean.split("/")[-1]
    if not safe_name.endswith(".SAFE"):
        raise ValueError(f"Invalid SAFE prefix: {safe_prefix!r}")

    product_id = safe_name[:-5]
    tile_match = _RE_TILE.search(product_id)
    if tile_match is None:
        raise ValueError(f"Cannot parse MGRS tile from SAFE prefix: {safe_prefix!r}")

    return S2Product(
        product_id=product_id,
        mgrs_tile=_normalize_mgrs_tile(tile_match.group(1)),
        sensing_date=_parse_sensing_datetime(product_id),
        processing_baseline=_parse_baseline(product_id),
        cloud_cover=100.0,
        satellite=_parse_satellite(product_id),
        orbit_number=_parse_orbit(product_id),
        source_url=f"gs://{GCS_BUCKET}/{clean}/",
        size_mb=None,
    )


def _matches_query(product: S2Product, query: S2Query) -> bool:
    if query.mgrs_tile and product.mgrs_tile != _normalize_mgrs_tile(query.mgrs_tile):
        return False

    if query.processing_level:
        level = query.processing_level.upper()
        if f"MSI{level}" not in product.product_id:
            return False

    sensing_day = product.sensing_date.date()
    if query.date is not None and sensing_day != query.date:
        return False
    if query.start_date is not None and sensing_day < query.start_date:
        return False
    return not (query.end_date is not None and sensing_day > query.end_date)


def _safe_prefix_from_source_url(source_url: str) -> str:
    if source_url.startswith(f"gs://{GCS_BUCKET}/"):
        prefix = source_url[len(f"gs://{GCS_BUCKET}/") :]
    elif source_url.startswith(GCS_DOWNLOAD_BASE):
        prefix = urllib.parse.unquote(source_url[len(GCS_DOWNLOAD_BASE) :])
    else:
        raise ValueError(f"Unsupported GCS source URL: {source_url!r}")

    clean = prefix.rstrip("/")
    if clean.endswith(".SAFE"):
        clean = f"{clean}/"
    elif ".SAFE/" in clean:
        clean = f"{clean.split('.SAFE/')[0]}.SAFE/"
    else:
        raise ValueError(f"Cannot infer SAFE prefix from source URL: {source_url!r}")
    return clean


def _resolve_safe_prefix(product: S2Product) -> str:
    if product.source_url:
        try:
            return _safe_prefix_from_source_url(product.source_url)
        except ValueError:
            logger.debug(
                "Could not parse SAFE prefix from source_url=%s; falling back to product_id.",
                product.source_url,
            )
    return _safe_prefix_from_product_id(_normalize_product_id(product.product_id))


def _object_download_url(object_name: str) -> str:
    return GCS_DOWNLOAD_BASE + urllib.parse.quote(object_name, safe="/")


def _download_url_to_file(url: str, target: Path, timeout: int = 300) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "siac-gcs-s2/1.0"})
    tmp = target.with_suffix(target.suffix + ".part")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp, tmp.open("wb") as fh:
            while True:
                chunk = resp.read(1024 * 1024)
                if not chunk:
                    break
                fh.write(chunk)
        tmp.replace(target)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise


def _target_path_for_object(object_name: str, safe_prefix: str, safe_dir: Path) -> Path | None:
    if object_name.endswith("_$folder$") or object_name.endswith("/"):
        return None
    if not object_name.startswith(safe_prefix):
        return None

    rel = object_name[len(safe_prefix) :].lstrip("/")
    if not rel:
        return None

    # Guard against malformed object keys with path traversal.
    candidate = (safe_dir / rel).resolve()
    safe_root = safe_dir.resolve()
    if candidate != safe_root and safe_root not in candidate.parents:
        raise DataNotFoundError(f"Refusing to write object outside SAFE root: {object_name!r}")
    return candidate


def _is_fully_downloaded(target: Path, expected_size: int | None) -> bool:
    if not target.exists():
        return False
    if expected_size is None:
        return True
    return target.stat().st_size == expected_size


def _download_with_retry(
    url: str,
    target: Path,
    *,
    expected_size: int | None,
    retries: int = GCS_DOWNLOAD_MAX_RETRIES,
    backoff_sec: float = GCS_DOWNLOAD_RETRY_BACKOFF_SEC,
) -> None:
    attempts = max(0, int(retries)) + 1
    for attempt in range(1, attempts + 1):
        try:
            _download_url_to_file(url, target)
            if expected_size is not None and target.stat().st_size != expected_size:
                raise OSError(
                    f"Downloaded size mismatch for {target.name}: "
                    f"expected={expected_size} got={target.stat().st_size}"
                )
            return
        except Exception as exc:
            target.unlink(missing_ok=True)
            if attempt >= attempts:
                raise RuntimeError(f"Failed downloading {url} after {attempts} attempts") from exc

            delay = max(0.0, float(backoff_sec)) * (2 ** (attempt - 1))
            logger.warning(
                "Retrying download (%d/%d) for %s after error: %s",
                attempt,
                attempts - 1,
                url,
                exc,
            )
            if delay > 0:
                time.sleep(delay)


def _download_jobs_parallel(
    jobs: list[tuple[str, Path, int | None]],
    *,
    max_workers: int = GCS_MAX_DOWNLOAD_WORKERS,
    retries: int = GCS_DOWNLOAD_MAX_RETRIES,
    backoff_sec: float = GCS_DOWNLOAD_RETRY_BACKOFF_SEC,
) -> None:
    if not jobs:
        return

    workers = max(1, min(int(max_workers), len(jobs)))
    if workers == 1:
        for url, target, expected_size in jobs:
            _download_with_retry(
                url,
                target,
                expected_size=expected_size,
                retries=retries,
                backoff_sec=backoff_sec,
            )
        return

    def _run(job: tuple[str, Path, int | None]) -> None:
        _download_with_retry(
            job[0],
            job[1],
            expected_size=job[2],
            retries=retries,
            backoff_sec=backoff_sec,
        )

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in pool.map(_run, jobs):
            pass


def search_gcs(query: S2Query) -> list[S2Product]:
    """List S2 products via GCS prefix listing (MGRS tile or product ID required)."""
    query.validate()

    if query.product_id:
        product_id = _normalize_product_id(query.product_id)
        try:
            safe_prefix = _safe_prefix_from_product_id(product_id)
        except ValueError:
            return []
        if not _prefix_exists(safe_prefix):
            return []
        product = _product_from_safe_prefix(safe_prefix)
        return [product] if _matches_query(product, query) else []

    if query.mgrs_tile is None:
        raise ValueError(
            "GCS backend supports product_id or mgrs_tile queries; bbox-only search is not supported."
        )

    tile = _normalize_mgrs_tile(query.mgrs_tile)
    tile_prefix = _mgrs_tile_to_prefix(tile)
    safe_prefixes = _list_safe_prefixes(tile_prefix)

    products: list[S2Product] = []
    for safe_prefix in safe_prefixes:
        try:
            product = _product_from_safe_prefix(safe_prefix)
        except ValueError:
            continue
        if _matches_query(product, query):
            products.append(product)

    if query.max_cloud_cover < 100.0:
        logger.warning(
            "GCS backend does not expose cloud-cover metadata via prefix listing; "
            "max_cloud_cover filter is ignored."
        )

    return sorted(products, key=lambda p: p.sensing_date, reverse=True)


def download_gcs(product: S2Product, dest_dir: Path) -> Path:
    """Download S2 SAFE directory from GCS public bucket."""
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)

    product_id = _normalize_product_id(product.product_id)
    safe_dir = dest / f"{product_id}.SAFE"

    safe_prefix = _resolve_safe_prefix(product)
    objects = _list_objects_under(safe_prefix)
    if not objects:
        raise DataNotFoundError(f"No objects found under GCS SAFE prefix: {safe_prefix!r}")

    safe_dir.mkdir(parents=True, exist_ok=True)
    jobs: list[tuple[str, Path, int | None]] = []
    for item in objects:
        name = item.get("name")
        if not isinstance(name, str):
            continue

        target = _target_path_for_object(name, safe_prefix=safe_prefix, safe_dir=safe_dir)
        if target is None:
            continue
        target.parent.mkdir(parents=True, exist_ok=True)

        size_raw = item.get("size")
        size_bytes = int(size_raw) if isinstance(size_raw, str) and size_raw.isdigit() else None
        if _is_fully_downloaded(target, size_bytes):
            continue

        jobs.append((_object_download_url(name), target, size_bytes))

    _download_jobs_parallel(jobs)
    downloaded = len(jobs)

    if not any(safe_dir.rglob("*")):
        raise DataNotFoundError(f"Downloaded {downloaded} objects but SAFE directory is empty: {safe_dir}")

    logger.info("Downloaded %d objects for %s to %s", downloaded, product_id, safe_dir)
    return safe_dir


class GCSSentinel2Backend:
    """GCS public bucket backend implementing the S2DataBackend protocol."""

    def search(self, query: S2Query) -> list[S2Product]:
        return search_gcs(query)

    def download(self, product: S2Product, dest_dir: Path) -> Path:
        return download_gcs(product, dest_dir)
