"""Validation helpers for SIAC runtime payloads."""

from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING

import numpy as np

from siac.errors import ValidationError

if TYPE_CHECKING:
    import xarray as xr

    from siac.runtime.models import (
        AtmosphericState,
        CorrectionResult,
        ObservationBundle,
        SolvedAtmosphere,
        SolverInputBundle,
        SurfacePrior,
    )


def spatial_shape(ds: xr.Dataset) -> tuple[int, ...]:
    """Extract ``(y, x)`` shape from an ``xr.Dataset``."""
    if "y" in ds.dims and "x" in ds.dims:
        return (ds.sizes["y"], ds.sizes["x"])
    dims = list(ds.dims)
    if len(dims) >= 2:
        return tuple(ds.sizes[d] for d in dims[-2:])
    raise ValueError("Dataset has fewer than 2 dimensions")


def validate_observation_bundle(obs: ObservationBundle) -> None:
    if "observation_time" not in obs.metadata:
        raise ValidationError("metadata must include 'observation_time'")
    if not isinstance(obs.metadata["observation_time"], datetime):
        raise ValidationError(
            f"metadata['observation_time'] must be a datetime, "
            f"got {type(obs.metadata['observation_time']).__name__}"
        )

    if len(obs.toa.data_vars) == 0:
        raise ValidationError("toa must have at least one band variable")
    if not obs.toa.sizes:
        raise ValidationError("toa must have spatial dimensions")

    toa_shape = spatial_shape(obs.toa)
    cloud_shape = obs.cloud_mask.shape
    if cloud_shape != toa_shape:
        raise ValidationError(
            f"cloud_mask shape {cloud_shape} must match TOA spatial shape {toa_shape}"
        )

    if len(obs.bounds) != 4:
        raise ValidationError(f"bounds must have 4 elements, got {len(obs.bounds)}")

    if not isinstance(obs.crs, str) or len(obs.crs) == 0:
        raise ValidationError("crs must be a non-empty string")

    for angle_name in ("sza", "saa", "vza", "vaa"):
        angle = getattr(obs.geometry, angle_name)
        if angle.ndim != 2 or angle.size == 0:
            raise ValidationError(
                f"geometry.{angle_name} must be a non-empty 2-D array, "
                f"got shape {angle.shape}"
            )


def validate_atmospheric_state(atmo: AtmosphericState) -> None:
    for name, field_val in [
        ("aot_unc", atmo.aot_unc),
        ("tcwv_unc", atmo.tcwv_unc),
        ("tco3_unc", atmo.tco3_unc),
    ]:
        vals = field_val.values
        if not (vals[np.isfinite(vals)] >= 0).all():
            raise ValidationError(f"{name} uncertainties must be non-negative")


def validate_surface_prior(prior: SurfacePrior) -> None:
    if prior.boa.shape != prior.boa_unc.shape:
        raise ValidationError(
            f"boa shape {prior.boa.shape} must match boa_unc shape {prior.boa_unc.shape}"
        )

    try:
        np.broadcast_shapes(prior.mask.shape, prior.boa.shape[-2:])
    except ValueError as err:
        raise ValidationError(
            f"mask shape {prior.mask.shape} must be broadcastable to "
            f"boa spatial shape {prior.boa.shape[-2:]}"
        ) from err


def validate_solver_input_bundle(sib: SolverInputBundle) -> None:
    config_bands = {b.name for b in sib.sensor_config.bands}
    solver_bands = {b.name for b in sib.bands}
    if not solver_bands <= config_bands:
        raise ValidationError(
            f"Solver bands {solver_bands - config_bands} not in sensor_config"
        )

    if sib.aux_resolution_m <= 0:
        raise ValidationError("aux_resolution_m must be positive")
    if sib.aerosol_resolution_m <= 0:
        raise ValidationError("aerosol_resolution_m must be positive")
    if sib.sharp_transition_mask is not None:
        expected_shape = sib.cloud_mask.shape
        actual_shape = sib.sharp_transition_mask.shape
        if actual_shape != expected_shape:
            raise ValidationError(
                "sharp_transition_mask shape "
                f"{actual_shape} must match cloud_mask shape {expected_shape}"
            )
        if tuple(sib.sharp_transition_mask.dims) != tuple(sib.cloud_mask.dims):
            raise ValidationError(
                "sharp_transition_mask dims "
                f"{tuple(sib.sharp_transition_mask.dims)} must match "
                f"cloud_mask dims {tuple(sib.cloud_mask.dims)}"
            )


def _validate_mask_dataset_shape(name: str, dataset: xr.Dataset, template: xr.DataArray) -> None:
    for var_name, field in dataset.data_vars.items():
        if field.shape != template.shape or tuple(field.dims) != tuple(template.dims):
            raise ValidationError(
                f"{name}.{var_name} shape {field.shape} dims {tuple(field.dims)} "
                f"must match {tuple(template.dims)} shape {template.shape}"
            )


def validate_solved_atmosphere(solved: SolvedAtmosphere) -> None:
    if not isinstance(solved.converged, bool):
        raise ValidationError("converged must be a boolean")
    if not isinstance(solved.n_iterations, int):
        raise ValidationError("n_iterations must be an int")
    if not isinstance(solved.cost_final, (int, float)):
        raise ValidationError("cost_final must be numeric")
    if solved.n_iterations < 0:
        raise ValidationError("n_iterations must be non-negative")

    aot_vals = solved.aot.values
    if not (aot_vals[np.isfinite(aot_vals)] >= 0).all():
        raise ValidationError("solved AOT must be non-negative")

    tcwv_vals = solved.tcwv.values
    if not (tcwv_vals[np.isfinite(tcwv_vals)] >= 0).all():
        raise ValidationError("solved TCWV must be non-negative")

    if solved.qa is not None:
        _validate_mask_dataset_shape("qa", solved.qa, solved.aot)


def validate_correction_result(result: CorrectionResult) -> None:
    if len(result.boa.data_vars) == 0:
        raise ValidationError("BOA must have at least one band variable")
    if not isinstance(result.metadata, dict):
        raise ValidationError("metadata must be a dictionary")
    processing_time_s = result.diagnostics.processing_time_s
    if processing_time_s is not None and (
        not isinstance(processing_time_s, (int, float)) or not np.isfinite(processing_time_s)
    ):
        raise ValidationError("diagnostics.processing_time_s must be finite when provided")
    if processing_time_s is not None and processing_time_s < 0:
        raise ValidationError("diagnostics.processing_time_s must be non-negative")

    if result.solver_qa is not None:
        _validate_mask_dataset_shape("solver_qa", result.solver_qa, result.aot)


__all__ = [
    "spatial_shape",
    "validate_atmospheric_state",
    "validate_correction_result",
    "validate_observation_bundle",
    "validate_solved_atmosphere",
    "validate_solver_input_bundle",
    "validate_surface_prior",
]
