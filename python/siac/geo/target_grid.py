"""The canonical target grid contract used by the SIAC pipeline.

Wave 18 (opt 3) follow-up: every pipeline stage that resamples an
upstream raster (priors, cloud bands, geometry) ultimately needs to
agree on the same target grid for the solver to consume. That grid is
defined by the scene's bounds, CRS, and the configured aerosol
retrieval resolution. Previously each call site reconstructed those
parameters ad-hoc from whatever ``target`` DataArray the caller had at
hand, which made it hard to (a) cache reprojection outputs across runs
and (b) reason about grid mismatches.

``TargetGrid`` captures that contract in one place. The signature is
hashable so it can serve as a cache key, and ``as_template_da`` produces
a backwards-compatible ``xr.DataArray`` for the legacy
``reproject_match`` callers that still expect a template raster.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    import xarray as xr


@dataclass(frozen=True)
class TargetGrid:
    """The (bounds, CRS, resolution, shape) tuple every pipeline stage aligns to.

    Attributes:
        bounds: ``(x_min, y_min, x_max, y_max)`` in ``crs``.
        crs: CRS as a string identifier (e.g. ``"EPSG:32633"``).
        resolution_m: Pixel resolution in CRS units (metres for projected CRS).
        shape: ``(height, width)`` — equals ``ceil((y_max - y_min) / res),
            ceil((x_max - x_min) / res)`` for the canonical SIAC grid.
    """

    bounds: tuple[float, float, float, float]
    crs: str
    resolution_m: float
    shape: tuple[int, int]

    def signature(self) -> tuple[Any, ...]:
        """Hashable identity used as a cache key component.

        Rounded to micrometre precision in bounds and millimetre in
        resolution to absorb floating-point round-trip noise from
        config parsing / rasterio coordinate transforms. Anything finer
        than that is below the float32 quantisation in the underlying
        rasters.
        """
        bounds_rounded = tuple(round(float(v), 6) for v in self.bounds)
        return (
            "TargetGrid/v1",
            bounds_rounded,
            self.crs,
            round(float(self.resolution_m), 3),
            tuple(int(v) for v in self.shape),
        )

    @classmethod
    def from_template(cls, template: xr.DataArray) -> TargetGrid:
        """Build a TargetGrid from an existing template DataArray.

        Backwards-compatible entry point — lets existing code that has a
        template ``xr.DataArray`` in hand opt into the canonical contract
        without restructuring the surrounding pipeline.
        """
        import rioxarray  # noqa: F401  # registers .rio accessor

        bounds = tuple(float(v) for v in template.rio.bounds())
        crs_obj = template.rio.crs
        if crs_obj is None:
            raise ValueError(
                "Cannot derive TargetGrid: template has no CRS metadata. "
                "Either set it via .rio.write_crs(...) or build the "
                "TargetGrid explicitly."
            )
        crs = crs_obj.to_string()
        # Use the y-axis pixel size (positive value) as the canonical
        # resolution. For SIAC's solver grids these are square pixels.
        try:
            transform = template.rio.transform()
            resolution = float(abs(transform.a))
        except Exception:
            # Fallback to bounds + shape arithmetic if the affine
            # transform isn't recoverable. Identical for square pixels.
            height = int(template.shape[-2])
            width = int(template.shape[-1])
            x_min, y_min, x_max, y_max = bounds
            resolution = float((x_max - x_min) / max(width, 1))

        if template.ndim < 2:
            raise ValueError(
                f"Template DataArray must be at least 2-D; got shape {template.shape}"
            )
        shape = (int(template.shape[-2]), int(template.shape[-1]))

        return cls(
            bounds=(
                float(bounds[0]),
                float(bounds[1]),
                float(bounds[2]),
                float(bounds[3]),
            ),
            crs=crs,
            resolution_m=resolution,
            shape=shape,
        )

    def as_template_da(self) -> xr.DataArray:
        """Build a minimal DataArray matching the grid for rioxarray APIs.

        Returns an empty float32 DataArray with the right shape, CRS, and
        transform so ``source.rio.reproject_match(target)`` works against
        it. Allocates ``height * width * 4`` bytes — for a 10980×10980
        S2 finest grid that's ~480 MB. The caller is responsible for
        discarding the template promptly.
        """
        import rioxarray  # noqa: F401
        import xarray as xr
        from rasterio.transform import from_bounds

        height, width = self.shape
        x_min, y_min, x_max, y_max = self.bounds
        transform = from_bounds(x_min, y_min, x_max, y_max, width, height)
        # Build x/y coords (cell centres) so xarray dim alignment works.
        x_coords = x_min + (np.arange(width) + 0.5) * (x_max - x_min) / width
        y_coords = y_max - (np.arange(height) + 0.5) * (y_max - y_min) / height
        da = xr.DataArray(
            np.zeros(self.shape, dtype=np.float32),
            dims=("y", "x"),
            coords={"y": y_coords, "x": x_coords},
        )
        da = da.rio.write_crs(self.crs)
        da = da.rio.write_transform(transform)
        return da
