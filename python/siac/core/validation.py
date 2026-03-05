"""
Contract validation functions for SIAC pipeline module outputs.

Each function validates one module's output contract before it is passed
downstream. These are called automatically by the pipeline orchestrator.
"""

from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import xarray as xr

    from siac.core.types import (
        AtmosphericState,
        ObservationBundle,
        SolverInputBundle,
        SurfacePrior,
    )


def _spatial_shape(ds: xr.Dataset) -> tuple[int, ...]:
    """Extract (y, x) shape from an xr.Dataset."""
    if "y" in ds.dims and "x" in ds.dims:
        return (ds.sizes["y"], ds.sizes["x"])
    # Fallback: take last two dims
    dims = list(ds.dims)
    if len(dims) >= 2:
        return tuple(ds.sizes[d] for d in dims[-2:])
    raise ValueError("Dataset has fewer than 2 dimensions")


def _validate_observation_bundle(obs: ObservationBundle) -> None:
    """Validate M1 output before passing to M4.

    Raises:
        AssertionError: If any validation rule is violated.
    """
    # metadata must include observation_time as datetime
    assert "observation_time" in obs.metadata, (
        "metadata must include 'observation_time'"
    )
    assert isinstance(obs.metadata["observation_time"], datetime), (
        f"metadata['observation_time'] must be a datetime, "
        f"got {type(obs.metadata['observation_time']).__name__}"
    )

    # toa must have at least one band variable and spatial dimensions
    assert len(obs.toa.data_vars) > 0, "toa must have at least one band variable"
    assert obs.toa.sizes, "toa must have spatial dimensions"

    # cloud_mask must match toa spatial shape
    toa_shape = _spatial_shape(obs.toa)
    cloud_shape = obs.cloud_mask.shape
    assert cloud_shape == toa_shape, (
        f"cloud_mask shape {cloud_shape} must match TOA spatial shape {toa_shape}"
    )

    # bounds must be a 4-tuple of floats
    assert len(obs.bounds) == 4, f"bounds must have 4 elements, got {len(obs.bounds)}"

    # crs must be a non-empty string
    assert isinstance(obs.crs, str) and len(obs.crs) > 0, "crs must be a non-empty string"


def _validate_atmospheric_state(atmo: AtmosphericState) -> None:
    """Validate M2 output before passing to M4.

    Raises:
        AssertionError: If any validation rule is violated.
    """
    # Uncertainty arrays must be non-negative
    for name, field_val in [
        ("aot_unc", atmo.aot_unc),
        ("tcwv_unc", atmo.tcwv_unc),
        ("tco3_unc", atmo.tco3_unc),
    ]:
        vals = field_val.values
        assert (vals[np.isfinite(vals)] >= 0).all(), (
            f"{name} uncertainties must be non-negative"
        )


def _validate_surface_prior(prior: SurfacePrior) -> None:
    """Validate M3 output before passing to M4.

    Raises:
        AssertionError: If any validation rule is violated.
    """
    # boa and boa_unc must have matching shapes
    assert prior.boa.shape == prior.boa_unc.shape, (
        f"boa shape {prior.boa.shape} must match boa_unc shape {prior.boa_unc.shape}"
    )

    # mask must be broadcastable to boa spatial dimensions
    try:
        np.broadcast_shapes(prior.mask.shape, prior.boa.shape[-2:])
    except ValueError as e:
        raise AssertionError(
            f"mask shape {prior.mask.shape} must be broadcastable to "
            f"boa spatial shape {prior.boa.shape[-2:]}"
        ) from e


def _validate_solver_input_bundle(sib: SolverInputBundle) -> None:
    """Validate M4 output before passing to M5.

    Raises:
        AssertionError: If any validation rule is violated.
    """
    # bands must be a subset of sensor_config.bands
    config_bands = {b.name for b in sib.sensor_config.bands}
    solver_bands = {b.name for b in sib.bands}
    assert solver_bands <= config_bands, (
        f"Solver bands {solver_bands - config_bands} not in sensor_config"
    )

    # Resolution metadata must be positive
    assert sib.aux_resolution_m > 0, "aux_resolution_m must be positive"
    assert sib.aerosol_resolution_m > 0, "aerosol_resolution_m must be positive"
