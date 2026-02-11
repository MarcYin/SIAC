import os
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
from osgeo import gdal, ogr, osr

from SIAC.create_logger import create_logger
from SIAC.raster_boundary import get_boundary


logger = create_logger()

CDSE_STAC_ROOT = "https://stac.dataspace.copernicus.eu/v1"
CDSE_TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
    "protocol/openid-connect/token"
)
CDSE_S3_HTTP_ENDPOINT = "https://eodata.dataspace.copernicus.eu"

CDSE_S2_L1C_COLLECTION = "sentinel-2-l1c"
CDSE_DEM_COLLECTION_30 = "cop-dem-glo-30-dged-cog"
CDSE_DEM_COLLECTION_90 = "cop-dem-glo-90-dged-cog"
COPDEM_OPEN_S3_30 = "https://copernicus-dem-30m.s3.amazonaws.com"
COPDEM_OPEN_S3_90 = "https://copernicus-dem-90m.s3.amazonaws.com"

S2_BANDS = [
    "B01",
    "B02",
    "B03",
    "B04",
    "B05",
    "B06",
    "B07",
    "B08",
    "B8A",
    "B09",
    "B10",
    "B11",
    "B12",
]
S2_REQUIRED_METADATA_ASSETS = ["product_metadata", "granule_metadata"]


def _make_session(session=None):
    if session is None:
        return requests.Session()
    return session


def _raise_for_status_with_text(response):
    if response.ok:
        return
    txt = response.text[:4000]
    raise requests.HTTPError(
        "HTTP {} calling {}: {}".format(response.status_code, response.url, txt),
        response=response,
    )


def get_cdse_access_token(
    username,
    password,
    totp=None,
    client_id="cdse-public",
    timeout=60,
    session=None,
):
    """
    Get an OIDC access token for CDSE APIs (STAC HTTPS/OData downloads).
    """
    data = {
        "client_id": client_id,
        "username": username,
        "password": password,
        "grant_type": "password",
    }
    if totp is not None:
        data["totp"] = str(totp)
    ses = _make_session(session=session)
    resp = ses.post(
        CDSE_TOKEN_URL,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data=data,
        timeout=timeout,
    )
    _raise_for_status_with_text(resp)
    token = resp.json().get("access_token")
    if token is None:
        raise RuntimeError("CDSE token response did not include access_token.")
    return token


def _to_wgs84_geometry(geometry, source_srs):
    if geometry is None:
        return None
    geom = geometry.Clone()
    if source_srs is None:
        return geom
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(4326)
    if int(gdal.VersionInfo()) > 2500000:
        source_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        target_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transform = osr.CoordinateTransformation(source_srs, target_srs)
    geom.Transform(transform)
    return geom


def _geometry_from_feature_collection(obj):
    features = obj.get("features", [])
    if len(features) == 0:
        raise ValueError("FeatureCollection AOI has no features.")
    geom = features[0].get("geometry")
    if geom is None:
        raise ValueError("First feature in AOI FeatureCollection has no geometry.")
    return geom


def _geometry_from_path(aoi_path):
    # Raster AOI path.
    raster = gdal.Open(aoi_path)
    if raster is not None:
        geojson_wgs84 = get_boundary(aoi_path, to_wgs84=True)[0]
        return _geometry_from_feature_collection(json.loads(geojson_wgs84))

    # Vector AOI path.
    vector = ogr.Open(aoi_path)
    if vector is not None:
        layer = vector.GetLayer(0)
        source_srs = layer.GetSpatialRef()
        geom = None
        for feat in layer:
            g = feat.GetGeometryRef()
            if g is None:
                continue
            if geom is None:
                geom = g.Clone()
            else:
                geom = geom.Union(g)
        if geom is None:
            raise ValueError("AOI vector has no geometry: {}".format(aoi_path))
        geom = _to_wgs84_geometry(geom, source_srs)
        return json.loads(geom.ExportToJson())

    raise ValueError("Unsupported AOI path: {}".format(aoi_path))


def _normalize_aoi(aoi=None, bbox=None):
    """
    Returns (geometry_wgs84, bbox_wgs84).
    """
    if bbox is not None:
        if len(bbox) != 4:
            raise ValueError("bbox must be [min_lon, min_lat, max_lon, max_lat].")
        return None, [float(v) for v in bbox]

    if aoi is None:
        return None, None

    geometry = None
    if isinstance(aoi, dict):
        if aoi.get("type") == "FeatureCollection":
            geometry = _geometry_from_feature_collection(aoi)
        elif aoi.get("type") == "Feature":
            geometry = aoi.get("geometry")
        else:
            geometry = aoi
    elif isinstance(aoi, (list, tuple)) and len(aoi) == 4:
        return None, [float(v) for v in aoi]
    elif isinstance(aoi, str):
        if os.path.exists(aoi):
            geometry = _geometry_from_path(aoi)
        else:
            # Try GeoJSON string first, then WKT.
            try:
                obj = json.loads(aoi)
                if obj.get("type") == "FeatureCollection":
                    geometry = _geometry_from_feature_collection(obj)
                elif obj.get("type") == "Feature":
                    geometry = obj.get("geometry")
                else:
                    geometry = obj
            except Exception:
                wkt_geom = ogr.CreateGeometryFromWkt(aoi)
                if wkt_geom is None:
                    raise ValueError("Could not parse AOI input.")
                geometry = json.loads(wkt_geom.ExportToJson())
    else:
        raise ValueError("Unsupported AOI input type: {}".format(type(aoi)))

    if geometry is None:
        raise ValueError("Could not parse AOI geometry.")
    ogr_geom = ogr.CreateGeometryFromJson(json.dumps(geometry))
    if ogr_geom is None:
        raise ValueError("AOI geometry is invalid.")
    min_x, max_x, min_y, max_y = ogr_geom.GetEnvelope()
    return geometry, [min_x, min_y, max_x, max_y]


def _datetime_range(start_time=None, end_time=None):
    if start_time is None and end_time is None:
        return None
    start = ".." if start_time is None else str(start_time)
    end = ".." if end_time is None else str(end_time)
    return "{}/{}".format(start, end)


def _get_next_href(stac_page):
    for link in stac_page.get("links", []):
        if link.get("rel") == "next":
            return link.get("href"), (link.get("method", "GET") or "GET").upper()
    return None, None


def search_cdse_s2_items(
    aoi=None,
    bbox=None,
    start_time=None,
    end_time=None,
    datetime_range=None,
    max_cloud_cover=None,
    limit=50,
    max_items=None,
    collection=CDSE_S2_L1C_COLLECTION,
    stac_root=CDSE_STAC_ROOT,
    timeout=60,
    session=None,
):
    """
    Search Sentinel-2 L1C scenes from CDSE STAC.

    Returns a list of STAC Item dicts.
    """
    geometry, bbox = _normalize_aoi(aoi=aoi, bbox=bbox)
    if datetime_range is None:
        datetime_range = _datetime_range(start_time=start_time, end_time=end_time)

    payload = {"collections": [collection], "limit": int(limit)}
    if payload["limit"] < 1:
        payload["limit"] = 1
    if payload["limit"] > 1000:
        payload["limit"] = 1000
    if geometry is not None:
        payload["intersects"] = geometry
    elif bbox is not None:
        payload["bbox"] = bbox
    if datetime_range is not None:
        payload["datetime"] = datetime_range
    if max_cloud_cover is not None:
        payload["filter-lang"] = "cql2-json"
        payload["filter"] = {
            "op": "<=",
            "args": [{"property": "eo:cloud_cover"}, float(max_cloud_cover)],
        }

    ses = _make_session(session=session)
    url = stac_root.rstrip("/") + "/search"
    resp = ses.post(url, json=payload, timeout=timeout)
    _raise_for_status_with_text(resp)

    items = []
    page = resp.json()
    page_count = 0
    while True:
        page_count += 1
        feats = page.get("features", [])
        items.extend(feats)

        if max_items is not None and len(items) >= int(max_items):
            items = items[: int(max_items)]
            break

        next_href, next_method = _get_next_href(page)
        if next_href is None:
            break
        if page_count > 200:
            logger.warning("Stopping STAC pagination after 200 pages.")
            break

        if next_method == "POST":
            nxt = ses.post(next_href, timeout=timeout)
        else:
            nxt = ses.get(next_href, timeout=timeout)
        _raise_for_status_with_text(nxt)
        page = nxt.json()

    return items


def _normalize_s2_asset_keys(
    bands=None,
    include_product_metadata=True,
    include_granule_metadata=True,
    include_safe_manifest=False,
):
    keys = []
    if bands is None:
        keys.extend(S2_BANDS)
    else:
        for band in bands:
            b = str(band).upper()
            if b not in S2_BANDS:
                raise ValueError("Unsupported S2 band key: {}".format(band))
            keys.append(b)
    if include_product_metadata:
        keys.append("product_metadata")
    if include_granule_metadata:
        keys.append("granule_metadata")
    if include_safe_manifest:
        keys.append("safe_manifest")
    # Preserve order and remove duplicates.
    seen = set()
    normalized = []
    for key in keys:
        if key not in seen:
            normalized.append(key)
            seen.add(key)
    return normalized


def cdse_s3_to_https_url(
    s3_url,
    endpoint=CDSE_S3_HTTP_ENDPOINT,
    drop_bucket=True,
):
    """
    Convert s3://eodata/... to https://eodata.dataspace.copernicus.eu/...
    """
    if not str(s3_url).startswith("s3://"):
        return s3_url

    path = str(s3_url)[5:]
    if "/" not in path:
        raise ValueError("Invalid s3 URL: {}".format(s3_url))
    bucket, key = path.split("/", 1)
    if drop_bucket and bucket == "eodata":
        return endpoint.rstrip("/") + "/" + key
    return endpoint.rstrip("/") + "/" + bucket + "/" + key


def _copdem_open_s3_root(item):
    collection = item.get("collection")
    if collection == CDSE_DEM_COLLECTION_30:
        return COPDEM_OPEN_S3_30
    if collection == CDSE_DEM_COLLECTION_90:
        return COPDEM_OPEN_S3_90

    item_id = str(item.get("id", ""))
    if "_COG_10_" in item_id:
        return COPDEM_OPEN_S3_30
    if "_COG_30_" in item_id:
        return COPDEM_OPEN_S3_90
    return None


def _compose_copdem_open_url(item):
    """
    Compose public-open CopDEM URL:
      - 30m: https://copernicus-dem-30m.s3.amazonaws.com/{item_id}/{item_id}.tif
      - 90m: https://copernicus-dem-90m.s3.amazonaws.com/{item_id}/{item_id}.tif
    """
    item_id = item.get("id")
    if item_id is None:
        return None
    root = _copdem_open_s3_root(item)
    if root is None:
        return None
    return "{}/{}/{}.tif".format(root.rstrip("/"), item_id, item_id)


def _resolve_asset_url(asset):
    alt_https = (
        asset.get("alternate", {})
        .get("https", {})
        .get("href")
    )
    if alt_https is not None:
        auth_refs = asset.get("alternate", {}).get("https", {}).get("auth:refs", [])
        return alt_https, auth_refs

    href = asset.get("href")
    if href is None:
        return None, []
    auth_refs = asset.get("auth:refs", [])
    if str(href).startswith("s3://"):
        return cdse_s3_to_https_url(href), auth_refs
    return href, auth_refs


def _safe_root_from_item(item):
    assets = item.get("assets", {})
    for asset in assets.values():
        local_path = asset.get("file:local_path")
        if local_path is None:
            continue
        marker = ".SAFE"
        idx = local_path.find(marker)
        if idx >= 0:
            return local_path[: idx + len(marker)]
    return item["id"] + ".SAFE"


def _download_http_file(url, destination, headers=None, timeout=180, retries=3):
    if headers is None:
        headers = {}
    attempt = 0
    while True:
        attempt += 1
        try:
            req = requests.get(url, headers=headers, stream=True, timeout=timeout)
            _raise_for_status_with_text(req)
            with open(destination, "wb") as fh:
                for chunk in req.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        fh.write(chunk)
            return destination
        except Exception:
            if attempt >= retries:
                raise
            time.sleep(2 ** (attempt - 1))


def download_cdse_s2_item_assets(
    item,
    output_dir,
    access_token=None,
    bands=None,
    include_product_metadata=True,
    include_granule_metadata=True,
    include_safe_manifest=False,
    overwrite=False,
    max_workers=4,
    timeout=180,
):
    """
    Download selected assets from a Sentinel-2 STAC item into a .SAFE-like tree.
    """
    assets = item.get("assets", {})
    asset_keys = _normalize_s2_asset_keys(
        bands=bands,
        include_product_metadata=include_product_metadata,
        include_granule_metadata=include_granule_metadata,
        include_safe_manifest=include_safe_manifest,
    )
    missing = [k for k in asset_keys if k not in assets]
    if len(missing) > 0:
        raise KeyError(
            "Missing required S2 assets in item {}: {}".format(item["id"], missing)
        )

    out_dir = os.path.abspath(output_dir)
    os.makedirs(out_dir, exist_ok=True)

    jobs = []
    for key in asset_keys:
        asset = assets[key]
        url, auth_refs = _resolve_asset_url(asset)
        if url is None:
            raise ValueError("Asset {} has no download URL.".format(key))

        local_rel = asset.get("file:local_path")
        if local_rel is None:
            local_rel = os.path.join(item["id"] + ".SAFE", key)
        target = os.path.join(out_dir, local_rel)
        os.makedirs(os.path.dirname(target), exist_ok=True)

        if (not overwrite) and os.path.exists(target):
            continue

        headers = {}
        if "oidc" in auth_refs:
            if access_token is None:
                raise ValueError(
                    "Asset {} requires OIDC token, but no access_token was provided.".format(
                        key
                    )
                )
            headers["Authorization"] = "Bearer {}".format(access_token)

        jobs.append((key, url, target, headers))

    if len(jobs) > 0:
        workers = max(1, int(max_workers))
        with ThreadPoolExecutor(max_workers=workers) as exe:
            futs = {
                exe.submit(
                    _download_http_file,
                    url,
                    target,
                    headers=headers,
                    timeout=timeout,
                ): (key, target)
                for key, url, target, headers in jobs
            }
            for fut in as_completed(futs):
                key, target = futs[fut]
                fut.result()
                logger.info("Downloaded %s -> %s", key, target)

    safe_root = os.path.join(out_dir, _safe_root_from_item(item))
    return safe_root


def search_and_download_cdse_s2(
    output_dir,
    aoi=None,
    bbox=None,
    start_time=None,
    end_time=None,
    datetime_range=None,
    max_cloud_cover=None,
    limit=50,
    max_items=1,
    access_token=None,
    bands=None,
    include_safe_manifest=False,
    overwrite=False,
):
    """
    Convenience wrapper: STAC search + selected asset download.
    """
    items = search_cdse_s2_items(
        aoi=aoi,
        bbox=bbox,
        start_time=start_time,
        end_time=end_time,
        datetime_range=datetime_range,
        max_cloud_cover=max_cloud_cover,
        limit=limit,
        max_items=max_items,
    )
    safe_dirs = []
    for item in items:
        safe_dir = download_cdse_s2_item_assets(
            item=item,
            output_dir=output_dir,
            access_token=access_token,
            bands=bands,
            include_product_metadata=True,
            include_granule_metadata=True,
            include_safe_manifest=include_safe_manifest,
            overwrite=overwrite,
        )
        safe_dirs.append(safe_dir)
    return items, safe_dirs


def search_cdse_dem_items(
    aoi=None,
    bbox=None,
    prefer_30m=True,
    limit=500,
    max_items=None,
    timeout=60,
    session=None,
):
    """
    Search DEM tiles intersecting AOI from CDSE CopDEM collections.
    """
    geometry, bbox = _normalize_aoi(aoi=aoi, bbox=bbox)
    collections = [CDSE_DEM_COLLECTION_30, CDSE_DEM_COLLECTION_90]
    if not prefer_30m:
        collections = [CDSE_DEM_COLLECTION_90, CDSE_DEM_COLLECTION_30]

    ses = _make_session(session=session)
    all_items = []
    seen_ids = set()
    for coll in collections:
        payload = {"collections": [coll], "limit": int(limit)}
        if payload["limit"] < 1:
            payload["limit"] = 1
        if payload["limit"] > 1000:
            payload["limit"] = 1000
        if geometry is not None:
            payload["intersects"] = geometry
        elif bbox is not None:
            payload["bbox"] = bbox

        resp = ses.post(
            CDSE_STAC_ROOT.rstrip("/") + "/search",
            json=payload,
            timeout=timeout,
        )
        _raise_for_status_with_text(resp)
        page = resp.json()
        page_count = 0
        while True:
            page_count += 1
            feats = page.get("features", [])
            for feat in feats:
                fid = feat.get("id")
                if fid in seen_ids:
                    continue
                all_items.append(feat)
                seen_ids.add(fid)
                if max_items is not None and len(all_items) >= int(max_items):
                    return all_items

            next_href, next_method = _get_next_href(page)
            if next_href is None:
                break
            if page_count > 200:
                logger.warning("Stopping DEM STAC pagination after 200 pages.")
                break
            if next_method == "POST":
                nxt = ses.post(next_href, timeout=timeout)
            else:
                nxt = ses.get(next_href, timeout=timeout)
            _raise_for_status_with_text(nxt)
            page = nxt.json()
    return all_items


def dem_items_to_http_urls(dem_items, use_open_access_s3=True, drop_bucket=True):
    urls = []
    for item in dem_items:
        if use_open_access_s3:
            open_url = _compose_copdem_open_url(item)
            if open_url is not None:
                urls.append(open_url)
                continue

        asset = item.get("assets", {}).get("data", {})
        href = asset.get("href")
        if href is None:
            continue
        urls.append(cdse_s3_to_https_url(href, drop_bucket=drop_bucket))
    # stable order, no duplicates
    seen = set()
    deduped = []
    for url in urls:
        if url not in seen:
            deduped.append(url)
            seen.add(url)
    return deduped


def search_cdse_dem_http_urls(
    aoi=None,
    bbox=None,
    prefer_30m=True,
    limit=500,
    max_items=None,
    use_open_access_s3=True,
    drop_bucket=True,
):
    dem_items = search_cdse_dem_items(
        aoi=aoi,
        bbox=bbox,
        prefer_30m=prefer_30m,
        limit=limit,
        max_items=max_items,
    )
    return dem_items_to_http_urls(
        dem_items,
        use_open_access_s3=use_open_access_s3,
        drop_bucket=drop_bucket,
    )


def _as_vsicurl_url(url):
    if str(url).startswith("/vsicurl/"):
        return url
    return "/vsicurl/{}".format(url)


def open_cdse_dem_crop(
    dem_s2_http_urls,
    aoi,
    xRes=30,
    yRes=30,
    crs=None,
    access_token=None,
    chunks=None,
    masked=True,
):
    """
    Open and crop DEM tiles to AOI using rioxarray over /vsicurl/.

    Usage pattern follows:
        rioxarray.open_rasterio([f"/vsicurl/{i}" for i in dem_s2_http_urls], ...)
    """
    try:
        import rasterio
        import rioxarray
        from rioxarray.merge import merge_arrays
    except ImportError:
        raise ImportError(
            "open_cdse_dem_crop requires rasterio+rioxarray. "
            "Please install rioxarray to use DEM crop by /vsicurl/."
        )

    geometry, _ = _normalize_aoi(aoi=aoi, bbox=None)
    if geometry is None:
        raise ValueError("AOI is required for DEM crop.")

    vsi_urls = [_as_vsicurl_url(url) for url in dem_s2_http_urls]
    env_kwargs = {}
    if access_token is not None:
        env_kwargs["GDAL_HTTP_HEADERS"] = "Authorization: Bearer {}".format(access_token)

    arrays = []
    with rasterio.Env(**env_kwargs):
        for url in vsi_urls:
            arr = rioxarray.open_rasterio(url, chunks=chunks, masked=masked)
            if "band" in arr.dims and arr.sizes.get("band", 1) == 1:
                arr = arr.squeeze("band", drop=True)
            arrays.append(arr)

        merged = merge_arrays(arrays, res=(float(xRes), float(yRes)), crs=crs)
        cropped = merged.rio.clip([geometry], crs="EPSG:4326", drop=True)
    return cropped
