"""
Copernicus Data Space Ecosystem (CDSE) backend for Sentinel-2.

Uses CDSE's S3-compatible API for download and OData for search.
See PLANS_S2.md §2, Phase 2.

Status: Stub — search/download functions defined but not yet implemented.
"""

from __future__ import annotations

import logging
from pathlib import Path

from siac.io.s2_data_source import S2Product, S2Query

logger = logging.getLogger(__name__)

# Constants
CDSE_ENDPOINT = "https://eodata.dataspace.copernicus.eu"
CDSE_BUCKET = "eodata"
CDSE_PREFIX = "Sentinel-2/MSI/L1C"
CDSE_ODATA_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1"


def search_cdse(
    query: S2Query,
    access_key: str | None = None,
    secret_key: str | None = None,
) -> list[S2Product]:
    """Search CDSE OData catalogue for S2 products.

    Not yet implemented — raises NotImplementedError.
    """
    raise NotImplementedError("CDSE search not yet implemented")


def download_cdse(
    product: S2Product,
    dest_dir: Path,
    access_key: str | None = None,
    secret_key: str | None = None,
) -> Path:
    """Download S2 SAFE directory from CDSE S3 bucket.

    Not yet implemented — raises NotImplementedError.
    """
    raise NotImplementedError("CDSE download not yet implemented")


class CopernicusDataspaceBackend:
    """CDSE S3 backend implementing the S2DataBackend protocol."""

    def __init__(
        self,
        access_key: str | None = None,
        secret_key: str | None = None,
    ):
        self._access_key = access_key
        self._secret_key = secret_key

    def search(self, query: S2Query) -> list[S2Product]:
        return search_cdse(query, self._access_key, self._secret_key)

    def download(self, product: S2Product, dest_dir: Path) -> Path:
        return download_cdse(product, dest_dir, self._access_key, self._secret_key)
