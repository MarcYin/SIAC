"""
Copernicus Data Space Ecosystem (CDSE) backend for Sentinel-2.

This module implements ``search_cdse()`` and ``download_cdse()`` for the S2 data
access layer described in ``docs/PLANS_S2.md`` (section 2, phase 2).
"""

from __future__ import annotations

import logging
import re
import tempfile
import zipfile
from datetime import datetime, time
from pathlib import Path
from typing import TYPE_CHECKING, Any

import requests

from siac.errors import DataNotFoundError
from siac.io.s2_data_source import S2Product, S2Query

if TYPE_CHECKING:
    from siac.adapters.auth import CredentialManager

logger = logging.getLogger(__name__)

# Constants
CDSE_ENDPOINT = "https://eodata.dataspace.copernicus.eu"
CDSE_BUCKET = "eodata"
CDSE_PREFIX = "Sentinel-2/MSI/L1C"
CDSE_ODATA_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1"
CDSE_STAC_URL = "https://stac.dataspace.copernicus.eu/v1"
CDSE_S2_COLLECTION = "sentinel-2-l1c"
CDSE_TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
    "protocol/openid-connect/token"
)

_RE_TILE = re.compile(r"_T(\d{2}[A-Z]{3})_")
_RE_BASELINE = re.compile(r"_N(\d{4})_")
_RE_ORBIT = re.compile(r"_R(\d{3})_")
_RE_SAT = re.compile(r"^(S2[A-Z])_")


def _iso_datetime(dt: datetime) -> str:
    if dt.tzinfo is None:
        return dt.isoformat(timespec="seconds") + "Z"
    return dt.isoformat(timespec="seconds").replace("+00:00", "Z")


def _to_datetime_range(query: S2Query) -> str | None:
    if query.date is not None:
        start = datetime.combine(query.date, time.min)
        end = datetime.combine(query.date, time.max).replace(microsecond=0)
        return f"{_iso_datetime(start)}/{_iso_datetime(end)}"

    if query.start_date is None and query.end_date is None:
        return None

    if query.start_date is not None:
        start_dt = datetime.combine(query.start_date, time.min)
    else:
        return f"../{_iso_datetime(datetime.combine(query.end_date, time.max).replace(microsecond=0))}"

    if query.end_date is not None:
        end_dt = datetime.combine(query.end_date, time.max).replace(microsecond=0)
        return f"{_iso_datetime(start_dt)}/{_iso_datetime(end_dt)}"
    return f"{_iso_datetime(start_dt)}/.."


def _parse_iso_datetime(value: str | None, fallback_product_id: str) -> datetime:
    if value:
        try:
            return datetime.fromisoformat(value.replace("Z", "+00:00")).replace(tzinfo=None)
        except ValueError:
            pass

    # Fallback from product ID token: S2*_MSIL1C_YYYYMMDDTHHMMSS_...
    parts = fallback_product_id.split("_")
    if len(parts) >= 3:
        try:
            return datetime.strptime(parts[2], "%Y%m%dT%H%M%S")
        except ValueError:
            pass
    raise ValueError(f"Cannot parse sensing datetime from item {fallback_product_id!r}")


def _get_tile(product_id: str, props: dict[str, Any]) -> str:
    tile = props.get("s2:mgrs_tile")
    if isinstance(tile, str) and tile:
        return tile.lstrip("T")
    match = _RE_TILE.search(product_id)
    if match:
        return match.group(1)
    return ""


def _get_baseline(product_id: str) -> str:
    match = _RE_BASELINE.search(product_id)
    if match:
        return f"N{match.group(1)}"
    return "N0000"


def _get_orbit(product_id: str) -> int:
    match = _RE_ORBIT.search(product_id)
    if match:
        return int(match.group(1))
    return 0


def _get_satellite(product_id: str) -> str:
    match = _RE_SAT.search(product_id)
    if match:
        return match.group(1)
    return "S2"


def _pick_product_href(item: dict[str, Any]) -> tuple[str, float | None]:
    assets = item.get("assets", {})
    product_asset = assets.get("Product", {})
    if not product_asset:
        raise DataNotFoundError(
            f"CDSE STAC item {item.get('id', '<unknown>')} has no 'Product' asset."
        )

    href = product_asset.get("href")
    if not isinstance(href, str) or not href:
        raise DataNotFoundError(
            f"CDSE STAC item {item.get('id', '<unknown>')} Product asset has no href."
        )

    size_bytes = product_asset.get("file:size")
    size_mb = None
    if isinstance(size_bytes, (int, float)) and size_bytes > 0:
        size_mb = float(size_bytes) / (1024.0 * 1024.0)
    return href, size_mb


def _item_to_product(item: dict[str, Any]) -> S2Product:
    product_id = str(item["id"]).replace(".SAFE", "")
    props = item.get("properties", {})
    sensing_date = _parse_iso_datetime(props.get("datetime"), product_id)
    tile = _get_tile(product_id, props)
    cloud_cover = float(props.get("eo:cloud_cover", 100.0) or 100.0)
    source_url, size_mb = _pick_product_href(item)
    return S2Product(
        product_id=product_id,
        mgrs_tile=tile,
        sensing_date=sensing_date,
        processing_baseline=_get_baseline(product_id),
        cloud_cover=cloud_cover,
        satellite=_get_satellite(product_id),
        orbit_number=_get_orbit(product_id),
        source_url=source_url,
        size_mb=size_mb,
    )


def _next_link(page: dict[str, Any]) -> dict[str, Any] | None:
    for link in page.get("links", []):
        if link.get("rel") == "next":
            href = link.get("href")
            if isinstance(href, str) and href:
                return link
    return None


def _search_payload(query: S2Query, limit: int = 100) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "collections": [CDSE_S2_COLLECTION],
        "limit": limit,
    }
    if query.bbox is not None:
        payload["bbox"] = list(query.bbox)

    dt_range = _to_datetime_range(query)
    if dt_range is not None:
        payload["datetime"] = dt_range

    filters: list[dict[str, Any]] = []
    if query.mgrs_tile:
        filters.append({
            "op": "=",
            "args": [{"property": "grid:code"}, f"MGRS-{query.mgrs_tile.lstrip('T').upper()}"],
        })

    if query.max_cloud_cover < 100.0:
        filters.append({
            "op": "<=",
            "args": [{"property": "eo:cloud_cover"}, float(query.max_cloud_cover)],
        })

    if filters:
        payload["filter-lang"] = "cql2-json"
        payload["filter"] = filters[0] if len(filters) == 1 else {"op": "and", "args": filters}
    return payload


def _post_json(url: str, payload: dict[str, Any], timeout: int = 60) -> dict[str, Any]:
    resp = requests.post(url, json=payload, timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def _get_json(url: str, timeout: int = 60) -> dict[str, Any]:
    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def _load_page_from_link(link: dict[str, Any], timeout: int = 60) -> dict[str, Any]:
    href = link.get("href")
    if not isinstance(href, str) or not href:
        raise DataNotFoundError("CDSE STAC next link is missing href.")

    method = str(link.get("method", "GET")).upper()
    if method == "POST":
        body = link.get("body")
        if body is not None and not isinstance(body, dict):
            raise DataNotFoundError("CDSE STAC POST next link body is not a JSON object.")
        return _post_json(href, body or {}, timeout=timeout)
    return _get_json(href, timeout=timeout)


def _token_from_credentials(username: str, password: str, timeout: int = 60) -> str:
    token, _ = _token_exchange(username, password, timeout)
    return token


def _token_exchange(username: str, password: str, timeout: int = 60) -> tuple[str, int]:
    """Exchange credentials for an access token. Returns (token, expires_in)."""
    resp = requests.post(
        CDSE_TOKEN_URL,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data={
            "client_id": "cdse-public",
            "username": username,
            "password": password,
            "grant_type": "password",
        },
        timeout=timeout,
    )
    resp.raise_for_status()
    body = resp.json()
    token = body.get("access_token")
    if not isinstance(token, str) or not token:
        raise DataNotFoundError("CDSE token response does not include access_token.")
    expires_in = int(body.get("expires_in", 300))
    return token, expires_in


def _resolve_auth_header(access_key: str | None, secret_key: str | None) -> dict[str, str]:
    # If caller passes a bearer token directly in access_key, use it as-is.
    if access_key and not secret_key:
        return {"Authorization": f"Bearer {access_key}"}

    # Username/password mode.
    if access_key and secret_key:
        token = _token_from_credentials(access_key, secret_key)
        return {"Authorization": f"Bearer {token}"}
    return {}


def search_cdse(
    query: S2Query,
    access_key: str | None = None,
    secret_key: str | None = None,
) -> list[S2Product]:
    """Search CDSE STAC catalogue for Sentinel-2 L1C products."""
    del access_key, secret_key  # currently unused for search
    query.validate()

    # Direct lookup by product ID.
    if query.product_id:
        product_id = query.product_id.replace(".SAFE", "")
        url = f"{CDSE_STAC_URL}/collections/{CDSE_S2_COLLECTION}/items/{product_id}"
        try:
            item = _get_json(url)
        except requests.HTTPError as exc:
            status = exc.response.status_code if exc.response is not None else None
            if status == 404:
                return []
            raise
        return [_item_to_product(item)]

    payload = _search_payload(query=query)
    page = _post_json(f"{CDSE_STAC_URL}/search", payload)

    items: list[dict[str, Any]] = []
    page_count = 0
    while True:
        page_count += 1
        items.extend(page.get("features", []))
        link = _next_link(page)
        if link is None:
            break
        if page_count >= 50:
            logger.warning("CDSE search pagination limit reached (50 pages).")
            break
        page = _load_page_from_link(link)

    products = [_item_to_product(item) for item in items]

    # Client-side filters for fields that may not always be queryable server-side.
    if query.mgrs_tile:
        wanted = query.mgrs_tile.lstrip("T").upper()
        products = [p for p in products if p.mgrs_tile.upper() == wanted]

    if query.date:
        products = [p for p in products if p.sensing_date.date() == query.date]
    elif query.start_date or query.end_date:
        if query.start_date:
            products = [p for p in products if p.sensing_date.date() >= query.start_date]
        if query.end_date:
            products = [p for p in products if p.sensing_date.date() <= query.end_date]

    if query.max_cloud_cover < 100.0:
        products = [p for p in products if p.cloud_cover <= query.max_cloud_cover]

    if query.processing_level:
        level = query.processing_level.upper()
        products = [p for p in products if f"MSI{level}" in p.product_id]

    return products


def download_cdse(
    product: S2Product,
    dest_dir: Path,
    access_key: str | None = None,
    secret_key: str | None = None,
) -> Path:
    """Download and extract a Sentinel-2 SAFE ZIP from CDSE."""
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)

    safe_dir = dest / f"{product.product_id}.SAFE"
    if safe_dir.exists():
        return safe_dir

    headers = _resolve_auth_header(access_key=access_key, secret_key=secret_key)
    resp = requests.get(product.source_url, headers=headers, stream=True, timeout=300)
    resp.raise_for_status()

    with tempfile.NamedTemporaryFile(prefix="cdse_", suffix=".zip", delete=False, dir=dest) as tmp:
        tmp_path = Path(tmp.name)
        for chunk in resp.iter_content(chunk_size=1024 * 1024):
            if chunk:
                tmp.write(chunk)

    try:
        with zipfile.ZipFile(tmp_path) as zf:
            zf.extractall(dest)
    finally:
        tmp_path.unlink(missing_ok=True)

    if safe_dir.exists():
        return safe_dir

    # Fallback if product_id was slightly different to extracted folder name.
    candidates = sorted(dest.glob(f"{product.product_id}*.SAFE"))
    if candidates:
        return candidates[0]
    raise DataNotFoundError(
        f"Downloaded archive for {product.product_id!r} but SAFE directory was not found in {dest}."
    )


class CopernicusDataspaceBackend:
    """CDSE backend implementing the ``S2DataBackend`` protocol."""

    def __init__(
        self,
        access_key: str | None = None,
        secret_key: str | None = None,
        auth: CredentialManager | None = None,
    ):
        self._access_key = access_key
        self._secret_key = secret_key
        self._auth = auth

    def _get_auth_header(self) -> dict[str, str]:
        """Resolve auth header, preferring explicit keys over auth manager."""
        if self._access_key or self._secret_key:
            return _resolve_auth_header(self._access_key, self._secret_key)
        if self._auth is not None and self._auth.cdse().has_credentials():
            return self._auth.cdse().authorization_header()
        return {}

    def search(self, query: S2Query) -> list[S2Product]:
        return search_cdse(query, self._access_key, self._secret_key)

    def download(self, product: S2Product, dest_dir: Path) -> Path:
        if (
            not self._access_key
            and not self._secret_key
            and self._auth is not None
            and self._auth.cdse().has_credentials()
        ):
            token = self._auth.cdse().get_token()
            return download_cdse(product, dest_dir, access_key=token)
        return download_cdse(product, dest_dir, self._access_key, self._secret_key)
