"""Diagnostic builders for SIAC pipeline outputs."""

from __future__ import annotations

from typing import TypeAlias, cast

import numpy as np
import xarray as xr
from numpy import typing as npt

from siac.algorithms.solver import build_solver_valid_mask
from siac.runtime import (
    AOTScatterBandDiagnostics,
    AtmosphericState,
    SolvedAtmosphere,
    SolverInputBundle,
)
from siac.workflows._pipeline_config import _MAX_SCATTER_POINTS_PER_BAND

Float32Array: TypeAlias = npt.NDArray[np.float32]


def _sample_scatter_values(values: np.ndarray, *, max_points: int) -> Float32Array:
    if values.size <= max_points:
        return cast("Float32Array", values.astype(np.float32, copy=False))
    indices = np.linspace(0, values.size - 1, max_points, dtype=np.int64)
    return cast("Float32Array", values[indices].astype(np.float32, copy=False))


def _select_band_slice(
    data: xr.DataArray,
    *,
    band_name: str,
    band_index: int,
) -> xr.DataArray | None:
    if "band" not in data.dims:
        return data
    band_coord = data.coords.get("band")
    if band_coord is not None:
        band_values = [str(value) for value in np.asarray(band_coord.values).tolist()]
        if band_name in band_values:
            return cast("xr.DataArray", data.sel(band=band_name, drop=True))
        if np.asarray(band_coord.values).dtype.kind in {"U", "S", "O"}:
            return None
    return cast("xr.DataArray", data.isel(band=band_index, drop=True))


def _finite_diagnostic_field(
    values: xr.DataArray,
    fallback: xr.DataArray,
) -> xr.DataArray:
    source = np.asarray(values.values, dtype=np.float32)
    if np.all(np.isfinite(source)):
        return values

    filled = source.copy()
    missing = ~np.isfinite(filled)
    if fallback.shape == values.shape:
        fallback_values = np.asarray(fallback.values, dtype=np.float32)
        fallback_finite = missing & np.isfinite(fallback_values)
        filled[fallback_finite] = fallback_values[fallback_finite]
        missing = ~np.isfinite(filled)

    if np.any(missing):
        finite_values = filled[np.isfinite(filled)]
        fill_value = float(np.mean(finite_values)) if finite_values.size else 0.0
        filled[missing] = np.float32(fill_value)

    return xr.DataArray(
        filled,
        dims=values.dims,
        coords=values.coords,
        attrs=values.attrs,
        name=values.name,
    )


def build_aot_scatter_diagnostics(
    solver_inputs: SolverInputBundle,
    solved: SolvedAtmosphere,
    *,
    max_points_per_band: int = _MAX_SCATTER_POINTS_PER_BAND,
) -> tuple[AOTScatterBandDiagnostics, ...]:
    valid_mask = build_solver_valid_mask(
        solver_inputs.cloud_mask,
        solver_inputs.toa,
        solver_inputs.surface_prior,
        sharp_transition_mask=solver_inputs.sharp_transition_mask,
        water_mask=solver_inputs.water_mask,
    ).values.astype(bool)
    atmo_state = solved.atmo_state
    atmo_finite_mask = (
        np.isfinite(atmo_state.aot.values)
        & np.isfinite(atmo_state.tcwv.values)
        & np.isfinite(atmo_state.tco3.values)
        & np.isfinite(atmo_state.elevation.values)
    )
    valid_mask = valid_mask & atmo_finite_mask
    diagnostic_atmo_state = AtmosphericState(
        aot=_finite_diagnostic_field(atmo_state.aot, solver_inputs.atmo_prior.aot),
        tcwv=_finite_diagnostic_field(atmo_state.tcwv, solver_inputs.atmo_prior.tcwv),
        tco3=_finite_diagnostic_field(atmo_state.tco3, solver_inputs.atmo_prior.tco3),
        elevation=_finite_diagnostic_field(
            atmo_state.elevation, solver_inputs.atmo_prior.elevation
        ),
        aot_unc=atmo_state.aot_unc,
        tcwv_unc=atmo_state.tcwv_unc,
        tco3_unc=atmo_state.tco3_unc,
    )
    diagnostics: list[AOTScatterBandDiagnostics] = []

    for band_index, band in enumerate(solver_inputs.bands):
        toa_band = _select_band_slice(solver_inputs.toa, band_name=band.name, band_index=band_index)
        surface_band = _select_band_slice(
            solver_inputs.surface_prior.boa,
            band_name=band.name,
            band_index=band_index,
        )
        if toa_band is None or surface_band is None:
            continue
        coeffs = solver_inputs.rt_model.compute_coefficients(
            solver_inputs.geometry,
            diagnostic_atmo_state,
            band,
            False,
        )
        simulated_toa = coeffs.simulate_toa(surface_band)

        band_valid = (
            valid_mask
            & np.isfinite(surface_band.values)
            & np.isfinite(toa_band.values)
            & np.isfinite(simulated_toa.values)
        )
        if not np.any(band_valid):
            continue

        surface_values = surface_band.values[band_valid].astype(np.float32, copy=False)
        observed_values = toa_band.values[band_valid].astype(np.float32, copy=False)
        simulated_values = simulated_toa.values[band_valid].astype(np.float32, copy=False)

        if surface_values.size > max_points_per_band:
            order = np.argsort(surface_values, kind="mergesort")
            surface_values = surface_values[order]
            observed_values = observed_values[order]
            simulated_values = simulated_values[order]
        diagnostics.append(
            AOTScatterBandDiagnostics(
                band_name=band.name,
                surface_reflectance=_sample_scatter_values(
                    surface_values, max_points=max_points_per_band
                ),
                observed_toa=_sample_scatter_values(
                    observed_values, max_points=max_points_per_band
                ),
                simulated_toa=_sample_scatter_values(
                    simulated_values, max_points=max_points_per_band
                ),
                total_valid_count=int(np.count_nonzero(band_valid)),
            )
        )

    return tuple(diagnostics)
