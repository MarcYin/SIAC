"""
Sentinel-2 data access — search, select, and download.

This module sits **before** M1 in the pipeline.  It resolves a user query
(product ID, tile+date shorthand, spatial/temporal search) into a local
SAFE directory path that the S2 preprocessor can consume.

See PLANS_S2.md §2 for the full specification.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Protocol

logger = logging.getLogger(__name__)


# ── Data types ─────────────────────────────────────────────────────────

@dataclass
class S2Query:
    """Flexible query for Sentinel-2 products."""

    product_id: str | None = None
    mgrs_tile: str | None = None
    date: date | None = None
    start_date: date | None = None
    end_date: date | None = None
    bbox: tuple[float, float, float, float] | None = None  # (W, S, E, N) WGS84
    max_cloud_cover: float = 100.0
    processing_level: str = "L1C"

    @classmethod
    def from_product_id(cls, product_id: str) -> S2Query:
        """Create a query from a full SAFE product ID."""
        return cls(product_id=product_id)

    @classmethod
    def from_tile_date(cls, tile_date: str) -> S2Query:
        """Parse ``'T31UDQ_20210801'`` or ``'31UDQ_20210801'``.

        Raises:
            ValueError: If the string cannot be parsed.
        """
        # Strip leading T if present
        s = tile_date.strip()
        m = re.match(r"^T?(\d{2}[A-Z]{3})_(\d{8})$", s)
        if not m:
            raise ValueError(
                f"Cannot parse tile+date shorthand: {tile_date!r}. "
                "Expected format: T31UDQ_20210801 or 31UDQ_20210801"
            )
        tile = m.group(1)
        d = datetime.strptime(m.group(2), "%Y%m%d").date()
        return cls(mgrs_tile=tile, date=d)

    def validate(self) -> None:
        """Ensure at least one spatial constraint is set.

        Raises:
            ValueError: If the query has no spatial constraint.
        """
        if not any([
            self.product_id,
            self.mgrs_tile,
            self.bbox,
        ]):
            raise ValueError(
                "S2Query must have at least one spatial constraint "
                "(product_id, mgrs_tile, or bbox)"
            )


@dataclass
class S2Product:
    """Metadata for a discovered Sentinel-2 product."""

    product_id: str
    mgrs_tile: str
    sensing_date: datetime
    processing_baseline: str  # e.g. "N0500"
    cloud_cover: float  # 0–100
    satellite: str  # "S2A" or "S2B"
    orbit_number: int
    source_url: str  # backend-specific URI
    size_mb: float | None = None

    @property
    def baseline_number(self) -> int:
        """Extract numeric baseline for comparison: ``'N0500'`` → ``500``."""
        return int(self.processing_baseline.replace("N", ""))


# ── Backend protocol ──────────────────────────────────────────────────

class S2DataBackend(Protocol):
    """Backend for searching and fetching S2 SAFE directories."""

    def search(self, query: S2Query) -> list[S2Product]: ...

    def download(self, product: S2Product, dest_dir: Path) -> Path: ...


# ── Plain helper functions (the real logic) ────────────────────────────

def deduplicate_products(products: list[S2Product]) -> list[S2Product]:
    """Group by (mgrs_tile, sensing_date), keep highest baseline_number.

    When multiple processing baselines exist for the same tile+date
    (e.g. N0301, N0400, N0500), only the newest baseline is kept.
    """
    best: dict[tuple[str, date], S2Product] = {}
    for p in products:
        key = (p.mgrs_tile, p.sensing_date.date() if isinstance(p.sensing_date, datetime) else p.sensing_date)
        existing = best.get(key)
        if existing is None or p.baseline_number > existing.baseline_number:
            best[key] = p
    return sorted(best.values(), key=lambda p: p.sensing_date, reverse=True)


def select_best_product(products: list[S2Product]) -> S2Product:
    """Pick the single best product: newest baseline, then newest sensing date.

    Raises:
        ValueError: If the product list is empty.
    """
    if not products:
        raise ValueError("No products to select from")
    deduped = deduplicate_products(products)
    return deduped[0]


def _parse_query(query: S2Query | str) -> S2Query:
    """Normalise a query argument to an S2Query instance."""
    if isinstance(query, S2Query):
        return query
    s = str(query).strip()
    # Try product ID first (contains .SAFE or MSIL1C)
    if "MSIL1C" in s or ".SAFE" in s:
        return S2Query.from_product_id(s.replace(".SAFE", ""))
    # Try tile+date shorthand
    try:
        return S2Query.from_tile_date(s)
    except ValueError:
        pass
    # Treat as product ID
    return S2Query.from_product_id(s)


def search_s2(
    backend: S2DataBackend,
    query: S2Query | str,
) -> list[S2Product]:
    """Search and deduplicate.  No download."""
    q = _parse_query(query)
    q.validate()
    raw = backend.search(q)
    return deduplicate_products(raw)


def fetch_s2(
    backend: S2DataBackend,
    query: S2Query | str,
    dest_dir: Path,
) -> Path:
    """Search, select best product, download.  Returns local SAFE path."""
    products = search_s2(backend, query)
    best = select_best_product(products)
    logger.info(f"Selected: {best.product_id} (baseline={best.processing_baseline})")
    return backend.download(best, dest_dir)


# ── Orchestrator class ─────────────────────────────────────────────────

class S2DataAccess:
    """Unified S2 data access — search, select, download.

    Holds a configured backend + cache directory.  Public methods are
    thin wrappers around the plain functions above.
    """

    def __init__(
        self,
        backend: S2DataBackend,
        cache_dir: Path | None = None,
    ):
        self._backend = backend
        self._cache_dir = cache_dir or Path.home() / ".cache" / "siac" / "s2"

    def get(self, query: S2Query | str, dest_dir: Path | None = None) -> Path:
        """Main entry point.  Returns path to local SAFE directory (input to M1)."""
        return fetch_s2(self._backend, query, dest_dir or self._cache_dir)

    def search(self, query: S2Query | str) -> list[S2Product]:
        """Search only (no download).  Returns deduplicated list."""
        return search_s2(self._backend, query)
