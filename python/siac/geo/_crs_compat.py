"""
Internal helpers for CRS comparison and dtype-aware resampling defaults.

These small utilities back :mod:`siac.storage.readers` and
:mod:`siac.geo.reprojection` so we can:

* compare CRS values by authority/WKT semantics rather than string match
  (REVIEW.md §2.8, §3.7) — ``"EPSG:4326"`` and the equivalent verbose WKT
  must compare equal;
* pick a sensible default resampling kernel from the input dtype
  (REVIEW.md §3.7) — bilinear corrupts integer class IDs and boolean masks.

This module is private; nothing here is part of the public ``siac.geo`` API.
"""

from __future__ import annotations

from typing import Any

from rasterio.enums import Resampling


def crs_equivalent(a: Any, b: Any) -> bool:
    """Compare two CRS values by authority / WKT semantics.

    Accepts anything :func:`pyproj.CRS.from_user_input` handles — EPSG strings
    (``"EPSG:4326"``), pyproj ``CRS`` instances, rasterio ``CRS`` instances,
    WKT strings, or ``None``. Two CRS are considered equivalent when pyproj
    decides they describe the same coordinate reference, which is independent
    of the textual encoding. Falls back to plain string equality if pyproj
    fails (corrupt or non-CRS input) so callers never see a spurious raise.
    """
    try:
        from pyproj import CRS as _PyprojCRS

        ca = _PyprojCRS.from_user_input(a) if a is not None else None
        cb = _PyprojCRS.from_user_input(b) if b is not None else None
        if ca is None and cb is None:
            return True
        if ca is None or cb is None:
            return False
        return ca == cb  # pyproj CRS __eq__ compares by WKT semantics
    except Exception:
        return str(a) == str(b)


def default_resampling_for_dtype(dtype: Any) -> Resampling:
    """Pick a sensible default resampling for the given dtype.

    Integer and boolean rasters (classification maps, cloud masks, scene
    classification layers) must use nearest-neighbour to keep their class
    identifiers intact. Float rasters fall through to bilinear, which is the
    historical default for continuous reflectance / radiance data.
    """
    import numpy as np

    try:
        kind = np.dtype(dtype).kind
    except TypeError:
        return Resampling.bilinear
    return Resampling.nearest if kind in {"i", "u", "b"} else Resampling.bilinear
