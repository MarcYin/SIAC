"""Earthaccess wrapper for NASA Earthdata access.

This module centralizes earthaccess login/search/open/download calls and keeps
all direct earthaccess interactions behind one boundary.
"""

from __future__ import annotations

import importlib
import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

from pyproj import Transformer

logger = logging.getLogger(__name__)


class EarthAccessSource:
    """Thin wrapper around the ``earthaccess`` package with lazy authentication."""

    def __init__(
        self,
        *,
        provider: str | None = None,
        login_strategy: str | None = None,
        persist: bool | None = None,
        login_kwargs: dict[str, Any] | None = None,
    ) -> None:
        self.provider = provider
        self.login_strategy = login_strategy
        self.persist = persist
        self.login_kwargs = dict(login_kwargs or {})

        self._earthaccess: Any | None = None
        self._authenticated = False

    @property
    def is_authenticated(self) -> bool:
        """Return whether ``login()`` has been executed in this instance."""
        return self._authenticated

    def _get_module(self) -> Any:
        if self._earthaccess is None:
            try:
                self._earthaccess = importlib.import_module("earthaccess")
            except ImportError as exc:  # pragma: no cover - env-dependent
                raise ImportError(
                    "earthaccess is required for NASA Earthdata access. "
                    "Install with `pip install earthaccess`."
                ) from exc
        return self._earthaccess

    def _ensure_auth(self) -> None:
        if self._authenticated:
            return

        ea = self._get_module()
        kwargs = dict(self.login_kwargs)
        if self.login_strategy is not None:
            kwargs.setdefault("strategy", self.login_strategy)
        if self.persist is not None:
            kwargs.setdefault("persist", self.persist)

        try:
            if kwargs:
                ea.login(**kwargs)
            else:
                ea.login()
        except TypeError:
            # Some earthaccess versions may not expose all kwargs.
            ea.login()

        self._authenticated = True

    @staticmethod
    def normalize_bounds_to_wgs84(
        bounds: tuple[float, float, float, float],
        crs: str,
    ) -> tuple[float, float, float, float]:
        """Convert bounds from ``crs`` to WGS84 longitude/latitude bounds."""
        xmin, ymin, xmax, ymax = bounds

        if crs.upper() in {"EPSG:4326", "CRS84", "OGC:CRS84"}:
            return (xmin, ymin, xmax, ymax)

        transformer = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        corners = [
            transformer.transform(xmin, ymin),
            transformer.transform(xmin, ymax),
            transformer.transform(xmax, ymin),
            transformer.transform(xmax, ymax),
        ]
        xs = [pt[0] for pt in corners]
        ys = [pt[1] for pt in corners]
        return (min(xs), min(ys), max(xs), max(ys))

    @staticmethod
    def to_cmr_bounding_box(bounds_wgs84: tuple[float, float, float, float]) -> str:
        """Format WGS84 bounds into CMR ``bounding_box`` string."""
        xmin, ymin, xmax, ymax = bounds_wgs84
        return f"{xmin},{ymin},{xmax},{ymax}"

    @staticmethod
    def temporal_window(obs_time: datetime, window_days: int) -> str:
        """Create CMR temporal range string centered on ``obs_time``."""
        delta = timedelta(days=max(0, int(window_days)))
        start = (obs_time - delta).strftime("%Y-%m-%dT%H:%M:%SZ")
        end = (obs_time + delta).strftime("%Y-%m-%dT%H:%M:%SZ")
        return f"{start},{end}"

    def search_datasets(
        self,
        *,
        keyword: str | None = None,
        short_name: str | None = None,
        provider: str | None = None,
        count: int | None = None,
        **kwargs: Any,
    ) -> list[Any]:
        """Search CMR datasets through earthaccess."""
        self._ensure_auth()
        ea = self._get_module()

        query = dict(kwargs)
        if keyword is not None:
            query["keyword"] = keyword
        if short_name is not None:
            query["short_name"] = short_name

        provider_name = provider if provider is not None else self.provider
        if provider_name is not None:
            query["provider"] = provider_name

        datasets = list(ea.search_datasets(**query))
        if count is None:
            return datasets
        return datasets[: max(0, int(count))]

    def search_granules(
        self,
        *,
        short_name: str,
        bounds: tuple[float, float, float, float] | None = None,
        crs: str = "EPSG:4326",
        temporal: str | tuple[Any, Any] | None = None,
        provider: str | None = None,
        count: int | None = None,
        **kwargs: Any,
    ) -> list[Any]:
        """Search CMR granules by short name with optional AOI/time filters."""
        self._ensure_auth()
        ea = self._get_module()

        query = dict(kwargs)
        query["short_name"] = short_name

        if bounds is not None:
            bbox_wgs84 = self.normalize_bounds_to_wgs84(bounds, crs)
            query["bounding_box"] = self.to_cmr_bounding_box(bbox_wgs84)

        if temporal is not None:
            if isinstance(temporal, tuple) and len(temporal) == 2:
                t0 = temporal[0].isoformat() if hasattr(temporal[0], "isoformat") else str(temporal[0])
                t1 = temporal[1].isoformat() if hasattr(temporal[1], "isoformat") else str(temporal[1])
                query["temporal"] = f"{t0},{t1}"
            else:
                query["temporal"] = str(temporal)

        provider_name = provider if provider is not None else self.provider
        if provider_name is not None:
            query["provider"] = provider_name

        granules = list(ea.search_data(**query))
        if count is None:
            return granules
        return granules[: max(0, int(count))]

    def open_granules(self, granules: list[Any]) -> Any:
        """Open granules through earthaccess data access helper."""
        self._ensure_auth()
        ea = self._get_module()
        return ea.open(granules)

    def download_granules(self, granules: list[Any], dest_dir: str | Path) -> list[Path]:
        """Download granules to ``dest_dir`` and return local paths."""
        self._ensure_auth()
        ea = self._get_module()

        dest = Path(dest_dir)
        dest.mkdir(parents=True, exist_ok=True)
        out = ea.download(granules, str(dest))

        if out is None:
            return []
        if isinstance(out, (str, Path)):
            return [Path(out)]
        if isinstance(out, list):
            return [Path(p) for p in out]
        return []
