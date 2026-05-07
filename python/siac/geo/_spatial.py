"""Pure-xarray spatial-metadata helpers.

This module hosts utility functions that copy CRS, transform, and spatial
coordinate metadata between ``xarray.DataArray`` instances. It deliberately
has zero dependencies on ``siac.runtime`` so that ``siac.geo`` can be a
true upstream layer in the package's import graph.

REVIEW.md §1.4 flagged the previous home of ``copy_spatial_metadata_like``
in ``siac.runtime.models`` as a layering inversion: ``siac.runtime``
imports from ``siac.geo`` indirectly through several call paths, while
``siac.geo.resample`` and ``siac.algorithms.surface.*`` reach back into
``siac.runtime.models`` to grab this utility. Moving the function here
breaks the cycle. ``siac.runtime.models`` re-exports the symbol for
backward compatibility, so callers that still do ``from siac.runtime.models
import copy_spatial_metadata_like`` keep working.
"""

from __future__ import annotations

from contextlib import suppress

import xarray as xr


def copy_spatial_metadata_like(data: xr.DataArray, reference: xr.DataArray) -> xr.DataArray:
    """Copy spatial coords/CRS/transform from *reference* onto *data* when possible."""
    import rioxarray  # noqa: F401

    out = data
    x_dim: str | None = None
    y_dim: str | None = None
    try:
        x_dim = reference.rio.x_dim
        y_dim = reference.rio.y_dim
    except Exception:
        if "x" in reference.dims and "y" in reference.dims:
            x_dim, y_dim = "x", "y"
        elif "longitude" in reference.dims and "latitude" in reference.dims:
            x_dim, y_dim = "longitude", "latitude"

    coord_updates = {
        dim: reference.coords[dim]
        for dim in (x_dim, y_dim)
        if dim is not None
        and dim in out.dims
        and dim in reference.coords
        and out.sizes[dim] == reference.sizes[dim]
    }
    if coord_updates:
        out = out.assign_coords(coord_updates)

    spatial_sizes_match = (
        x_dim is not None
        and y_dim is not None
        and x_dim in out.dims
        and y_dim in out.dims
        and out.sizes[x_dim] == reference.sizes[x_dim]
        and out.sizes[y_dim] == reference.sizes[y_dim]
    )
    if spatial_sizes_match:
        with suppress(Exception):
            out = out.rio.set_spatial_dims(x_dim=x_dim, y_dim=y_dim)

    try:
        ref_crs = reference.rio.crs
    except Exception:
        ref_crs = None
    if ref_crs is not None and spatial_sizes_match:
        out = out.rio.write_crs(ref_crs)
        with suppress(Exception):
            out = out.rio.write_transform(reference.rio.transform(recalc=True))
    return out
