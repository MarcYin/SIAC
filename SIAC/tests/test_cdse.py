from SIAC.cdse import _normalize_s2_asset_keys
from SIAC.cdse import _compose_copdem_open_url
from SIAC.cdse import _resolve_asset_url
from SIAC.cdse import _safe_root_from_item
from SIAC.cdse import dem_items_to_http_urls
from SIAC.cdse import cdse_s3_to_https_url


def test_cdse_s3_to_https_url_drop_bucket():
    s3 = "s3://eodata/auxdata/CopDEM_COG/foo.tif"
    got = cdse_s3_to_https_url(s3, endpoint="https://eodata.dataspace.copernicus.eu")
    assert got == "https://eodata.dataspace.copernicus.eu/auxdata/CopDEM_COG/foo.tif"


def test_cdse_s3_to_https_url_keep_bucket():
    s3 = "s3://eodata/auxdata/CopDEM_COG/foo.tif"
    got = cdse_s3_to_https_url(
        s3,
        endpoint="https://eodata.dataspace.copernicus.eu",
        drop_bucket=False,
    )
    assert got == "https://eodata.dataspace.copernicus.eu/eodata/auxdata/CopDEM_COG/foo.tif"


def test_resolve_asset_url_prefers_https_alternate():
    asset = {
        "href": "s3://eodata/path/B02.jp2",
        "alternate": {
            "https": {
                "href": "https://zipper.dataspace.copernicus.eu/odata/v1/.../$value",
                "auth:refs": ["oidc"],
            }
        },
    }
    url, auth_refs = _resolve_asset_url(asset)
    assert url.startswith("https://zipper.dataspace.copernicus.eu/")
    assert auth_refs == ["oidc"]


def test_normalize_s2_asset_keys_defaults_and_subset():
    all_keys = _normalize_s2_asset_keys()
    assert "B02" in all_keys
    assert "product_metadata" in all_keys
    assert "granule_metadata" in all_keys

    subset = _normalize_s2_asset_keys(
        bands=["b02", "B11"],
        include_product_metadata=False,
        include_granule_metadata=True,
    )
    assert subset == ["B02", "B11", "granule_metadata"]


def test_safe_root_from_item_prefers_file_local_path():
    item = {
        "id": "S2A_TEST",
        "assets": {
            "B02": {
                "file:local_path": (
                    "S2A_MSIL1C_20260101T000000_N0512_R000_T00XXX_20260101T000000.SAFE/"
                    "GRANULE/L1C_T00XXX_A000000_20260101T000000/IMG_DATA/T00XXX_B02.jp2"
                )
            }
        },
    }
    safe_root = _safe_root_from_item(item)
    assert safe_root.endswith(".SAFE")
    assert safe_root.startswith("S2A_MSIL1C_")


def test_compose_copdem_open_url_30m():
    item = {
        "id": "Copernicus_DSM_COG_10_N00_00_E009_00_DEM",
        "collection": "cop-dem-glo-30-dged-cog",
    }
    got = _compose_copdem_open_url(item)
    assert (
        got
        == "https://copernicus-dem-30m.s3.amazonaws.com/"
        "Copernicus_DSM_COG_10_N00_00_E009_00_DEM/"
        "Copernicus_DSM_COG_10_N00_00_E009_00_DEM.tif"
    )


def test_compose_copdem_open_url_90m():
    item = {
        "id": "Copernicus_DSM_COG_30_N00_00_E006_00_DEM",
        "collection": "cop-dem-glo-90-dged-cog",
    }
    got = _compose_copdem_open_url(item)
    assert (
        got
        == "https://copernicus-dem-90m.s3.amazonaws.com/"
        "Copernicus_DSM_COG_30_N00_00_E006_00_DEM/"
        "Copernicus_DSM_COG_30_N00_00_E006_00_DEM.tif"
    )


def test_dem_items_to_http_urls_prefers_open_access_buckets():
    items = [
        {
            "id": "Copernicus_DSM_COG_10_N00_00_E009_00_DEM",
            "collection": "cop-dem-glo-30-dged-cog",
            "assets": {"data": {"href": "s3://eodata/auxdata/CopDEM_COG/x.tif"}},
        },
        {
            "id": "Copernicus_DSM_COG_30_N00_00_E006_00_DEM",
            "collection": "cop-dem-glo-90-dged-cog",
            "assets": {"data": {"href": "s3://eodata/auxdata/CopDEM_COG/y.tif"}},
        },
    ]
    got = dem_items_to_http_urls(items, use_open_access_s3=True)
    assert got == [
        "https://copernicus-dem-30m.s3.amazonaws.com/"
        "Copernicus_DSM_COG_10_N00_00_E009_00_DEM/"
        "Copernicus_DSM_COG_10_N00_00_E009_00_DEM.tif",
        "https://copernicus-dem-90m.s3.amazonaws.com/"
        "Copernicus_DSM_COG_30_N00_00_E006_00_DEM/"
        "Copernicus_DSM_COG_30_N00_00_E006_00_DEM.tif",
    ]
