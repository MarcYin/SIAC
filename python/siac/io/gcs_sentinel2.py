"""
Google Cloud Storage backend for Sentinel-2.

Uses the public ``gs://gcp-public-data-sentinel-2`` bucket (no auth).
Supports MGRS-based listing only (no full catalogue search).
See PLANS_S2.md §2, Phase 3.

Status: Stub — search/download functions defined but not yet implemented.
"""

from __future__ import annotations

import logging
from pathlib import Path

from siac.io.s2_data_source import S2Product, S2Query

logger = logging.getLogger(__name__)

# Constants
GCS_BUCKET = "gcp-public-data-sentinel-2"


def search_gcs(query: S2Query) -> list[S2Product]:
    """List S2 products via GCS prefix listing (MGRS tile or product ID required).

    Not yet implemented — raises NotImplementedError.
    """
    raise NotImplementedError("GCS search not yet implemented")


def download_gcs(product: S2Product, dest_dir: Path) -> Path:
    """Download S2 SAFE directory from GCS public bucket.

    Not yet implemented — raises NotImplementedError.
    """
    raise NotImplementedError("GCS download not yet implemented")


class GCSSentinel2Backend:
    """GCS public bucket backend implementing the S2DataBackend protocol."""

    def search(self, query: S2Query) -> list[S2Product]:
        return search_gcs(query)

    def download(self, product: S2Product, dest_dir: Path) -> Path:
        return download_gcs(product, dest_dir)
