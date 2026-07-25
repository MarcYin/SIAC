"""Shared DEM-derived terrain features for L2A/current-RT experiments.

The terrain quantities are intentionally geometric: elevation, slope and the
local solar-incidence cosine.  They do not apply a terrain correction by
themselves.  This keeps the paired-surface experiment separate from a future
RT sidecar that explicitly supports a per-pixel elevation axis.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

TERRAIN_SOURCE = "COPERNICUS/DEM/GLO30_2024_1"


@dataclass(frozen=True)
class TerrainFields:
    """DEM elevation and local surface-normal terms on one projected grid."""

    elevation_m: np.ndarray
    slope_deg: np.ndarray
    gradient_east: np.ndarray
    gradient_north: np.ndarray
    valid: np.ndarray
    source: str = TERRAIN_SOURCE


def fetch_glo30_terrain(ee: Any, grid: dict[str, Any]) -> TerrainFields:
    """Fetch GLO-30 and derive slope/normal gradients on ``grid``.

    Gradients are evaluated after replacing DEM nodata by the local finite
    median.  ``valid`` retains the original DEM support so downstream paired
    samples never silently use those filled pixels.
    """
    from bestpixel._gee import get_patch

    elevation_m = np.asarray(
        get_patch(
            ee,
            ee.ImageCollection(TERRAIN_SOURCE).select("DEM").mosaic(),
            ["DEM"],
            grid=grid,
            is_asset=False,
        )[0],
        dtype=np.float32,
    )
    valid = np.isfinite(elevation_m)
    if not np.any(valid):
        raise RuntimeError("Copernicus GLO-30 has no finite elevation on the target grid")
    filled = np.where(valid, elevation_m, np.nanmedian(elevation_m[valid]))
    resolution_m = float(grid["res"])
    if not np.isfinite(resolution_m) or resolution_m <= 0.0:
        raise ValueError(f"invalid terrain grid resolution: {resolution_m!r}")
    gradient_south, gradient_east = np.gradient(filled, resolution_m, resolution_m)
    gradient_north = -gradient_south
    slope_deg = np.degrees(
        np.arctan(np.hypot(gradient_east, gradient_north))
    ).astype(np.float32)
    return TerrainFields(
        elevation_m=elevation_m,
        slope_deg=slope_deg,
        gradient_east=np.asarray(gradient_east, dtype=np.float32),
        gradient_north=np.asarray(gradient_north, dtype=np.float32),
        valid=valid,
    )


def local_solar_incidence(
    terrain: TerrainFields,
    *,
    sza_deg: float,
    saa_deg: float,
) -> np.ndarray:
    """Return the cosine of local solar incidence from terrain and sun angles."""
    sza = np.radians(float(sza_deg))
    saa = np.radians(float(saa_deg))
    sun_east = np.sin(sza) * np.sin(saa)
    sun_north = np.sin(sza) * np.cos(saa)
    sun_up = np.cos(sza)
    normal_scale = np.sqrt(
        1.0 + terrain.gradient_east**2 + terrain.gradient_north**2
    )
    incidence = (
        -terrain.gradient_east * sun_east
        - terrain.gradient_north * sun_north
        + sun_up
    ) / normal_scale
    return np.clip(incidence, -1.0, 1.0).astype(np.float32)
