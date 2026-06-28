"""Surface-driven aerosol solver.

An alternative to the Bayesian multi-grid inversion: instead of an L-BFGS-B
descent, sweep a fixed AOT axis at the prior TCWV, atmospherically correct the
TOA to BOA at each node, and score that against the surface prior. The per-node
cost cube is built here (vectorised numpy), then the heavy spatial median-pool
+ Gaussian-backstop + per-pixel arg-min is delegated to the Rust kernel
``surface_driven_pool_argmin``. This is the cheap, robust scheme the
experimental harness validated; it is selected by
``algorithms.solver.method = "surface_driven"``.

TCWV is held at the prior for the sweep and carried through to the solved state
(the surface-driven retrieval constrains AOT, not column water vapour).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import xarray as xr

from siac._rust_compat import surface_driven_pool_argmin
from siac.algorithms.solver._grid_search import (
    build_candidate_coeff_provider,
    prepare_grid_search_observations,
)
from siac.algorithms.solver.multigrid import SolverResult, _log_aot_axis

if TYPE_CHECKING:
    from siac.config.algorithms import SolverAlgorithmConfig
    from siac.domain import SensorBand
    from siac.domain.protocols import RTModelBackend
    from siac.runtime import (
        AtmosphericState,
        GeometryAngles,
        SolverInputBundle,
        SurfacePrior,
    )

#: BOA-uncertainty floor for the per-band chi-square (matches the harness
#: ``UNC_FLOOR``; keeps a single clean band from dominating the cost).
_MIN_BOA_UNC = 0.006


class SurfaceDrivenSolver:
    """Surface-prior-driven AOT solver over a fixed AOT axis."""

    def __init__(self, config: SolverAlgorithmConfig):
        self.config = config

    # -- backstop ----------------------------------------------------------
    def _backstop_unc(self, aot_prior: np.ndarray, aot_prior_unc_raw: np.ndarray) -> np.ndarray:
        """Per-pixel AOT-prior (CAMS) backstop sigma.

        Calibrated mode (default) widens the prior when clean and at the
        high-AOD tail and tightens it through the moderate band where the
        surface signal is shallow; otherwise a flat 50 % relative prior.
        """
        if not self.config.surface_driven_backstop_calibrated:
            unc = np.maximum(0.5 * aot_prior, 0.02)
        else:
            m = np.asarray(aot_prior, dtype=np.float64)
            loose = np.maximum(0.5 * m, 0.02)
            mid = np.maximum(0.07, 0.5 * m / (1.0 + np.exp(-(m - 0.5) / 0.15)))
            unc = np.where(m < 0.15, loose, mid)
        unc = np.where(np.isfinite(unc), unc, np.maximum(aot_prior_unc_raw, 0.02))
        return np.asarray(unc, dtype=np.float32)

    # -- solve -------------------------------------------------------------
    def solve(
        self,
        toa: xr.DataArray,
        surface_prior: SurfacePrior,
        geometry: GeometryAngles,
        atmo_prior: AtmosphericState,
        rt_model: RTModelBackend,
        cloud_mask: xr.DataArray,
        bands: list[SensorBand],
        sharp_transition_mask: xr.DataArray | None = None,
        water_mask: xr.DataArray | None = None,
    ) -> SolverResult:
        # sharp_transition_mask is accepted for solver-interface parity but the
        # surface-driven sweep does not use it (no cross-edge smoothness term).
        del sharp_transition_mask
        ny = int(cloud_mask.sizes["y"])
        nx = int(cloud_mask.sizes["x"])
        shape = (ny, nx)
        n_bands = len(bands)

        # Valid = not cloudy/invalid, optionally also not water.
        valid_mask = ~cloud_mask.values.astype(bool)
        if water_mask is not None:
            valid_mask = valid_mask & ~water_mask.values.astype(bool)

        aot_prior = atmo_prior.aot.values.astype(np.float32)
        tcwv_prior = atmo_prior.tcwv.values.astype(np.float32)
        aot_prior_unc_raw = atmo_prior.aot_unc.values.astype(np.float32)
        tcwv_prior_unc = np.maximum(atmo_prior.tcwv_unc.values.astype(np.float32), 1e-3)

        # AOT axis: log-spaced over the configured bounds. TCWV is fixed at the
        # scene-median prior for the sweep.
        n_aot = max(3, int(self.config.grid_search_aot_points))
        aot_axis = _log_aot_axis(self.config.aot_bounds[0], self.config.aot_bounds[1], n_aot)
        finite_tcwv = tcwv_prior[np.isfinite(tcwv_prior)]
        tcwv_fixed = (
            float(np.median(finite_tcwv))
            if finite_tcwv.size
            else float(np.mean(self.config.tcwv_bounds))
        )

        (
            toa_values,
            _no_observation_mask,
            boa_prior,
            boa_unc,
            support_mask,
        ) = prepare_grid_search_observations(
            toa=toa,
            surface_prior=surface_prior,
            n_bands=n_bands,
            shape=shape,
            min_boa_unc=_MIN_BOA_UNC,
            aot_prior=aot_prior,
            tcwv_prior=tcwv_prior,
            aot_prior_unc=np.maximum(aot_prior_unc_raw, 1e-3),
            tcwv_prior_unc=tcwv_prior_unc,
        )
        solve_valid = valid_mask & support_mask

        # Coefficient provider with TCWV fixed at the prior; one shared
        # (band, y, x) stack per AOT node.
        coeff_provider = build_candidate_coeff_provider(
            rt_model=rt_model,
            bands=bands,
            geometry=geometry,
            atmo_prior=atmo_prior,
            aot_axis=aot_axis,
            tcwv_axis=np.asarray([tcwv_fixed], dtype=np.float32),
            solve_aot=True,
            solve_tcwv=False,
            aot_bounds=self.config.aot_bounds,
            tcwv_bounds=self.config.tcwv_bounds,
            rt_sample_step=1,
            shape=shape,
        )

        # Build the raw per-node surface-cost cube (vectorised chi-square). f64
        # and the 1.0 (no-0.5) chi-square convention so the kernel's backstop
        # ``d^2/unc^2`` is on the same scale as the obs cost — a bit-exact match
        # to the reference numpy resolve.
        inv_var = 1.0 / np.square(np.maximum(boa_unc, 1e-6))  # (band, y, x)
        cube = np.full((n_aot, ny, nx), np.nan, dtype=np.float64)
        for k in range(n_aot):
            xap, xbp, xcp = coeff_provider(float(aot_axis[k]), tcwv_fixed)
            y = np.asarray(xap, np.float64) * toa_values - np.asarray(xbp, np.float64)
            denom = 1.0 + np.asarray(xcp, np.float64) * y
            with np.errstate(invalid="ignore", divide="ignore"):
                boa = np.where(np.abs(denom) > 1e-12, y / denom, np.nan)
            diff = boa - boa_prior
            term = inv_var * diff * diff  # (band, y, x)
            finite = np.isfinite(term)
            used = np.count_nonzero(finite, axis=0)
            cost = np.where(finite, term, 0.0).sum(axis=0)
            cube[k] = np.where(used > 0, cost, np.nan)

        aot_prior_unc = self._backstop_unc(aot_prior, aot_prior_unc_raw)
        # Full rolling-window edge length (pandas/xarray center convention) from
        # the ~1.2 km radius; min_periods matches the reference resolve.
        pool_window = max(
            1,
            int(
                round(
                    2.0
                    * self.config.surface_driven_pool_radius_m
                    / max(self.config.aerosol_resolution, 1e-3)
                )
            ),
        )
        min_count = (
            max(int(self.config.surface_driven_pool_min_count), 20, pool_window * pool_window // 5)
            if pool_window > 1
            else 1
        )

        aot_arr, aot_unc_arr, cost_arr = surface_driven_pool_argmin(
            np.ascontiguousarray(cube, dtype=np.float64),
            aot_axis.astype(np.float64, copy=False),
            np.ascontiguousarray(aot_prior, dtype=np.float64),
            np.ascontiguousarray(aot_prior_unc, dtype=np.float64),
            np.ascontiguousarray(solve_valid, dtype=bool),
            int(pool_window),
            int(min_count),
        )

        aot_arr = np.asarray(aot_arr, dtype=np.float32)
        aot_unc_arr = np.asarray(aot_unc_arr, dtype=np.float32)
        cost_arr = np.asarray(cost_arr, dtype=np.float32)

        # Pixels with no usable observation fall back to the prior (backstop).
        solved = np.isfinite(aot_arr)
        aot_filled = np.where(solved, aot_arr, aot_prior).astype(np.float32)
        aot_unc_filled = np.where(solved, aot_unc_arr, aot_prior_unc).astype(np.float32)

        template = atmo_prior.aot
        aot_da = xr.DataArray(aot_filled, dims=template.dims, coords=template.coords, name="aot")
        aot_unc_da = xr.DataArray(
            aot_unc_filled, dims=template.dims, coords=template.coords, name="aot_unc"
        )

        finite_cost = cost_arr[np.isfinite(cost_arr)]
        final_cost = float(np.median(finite_cost)) if finite_cost.size else float("inf")
        n_solved = int(np.count_nonzero(solved))

        solved_atmo = atmo_prior.with_updated_aot_tcwv(
            aot=aot_da,
            tcwv=atmo_prior.tcwv,
            aot_unc=aot_unc_da,
            tcwv_unc=atmo_prior.tcwv_unc,
        )
        return SolverResult(
            aot=solved_atmo.aot,
            tcwv=solved_atmo.tcwv,
            aot_unc=solved_atmo.aot_unc,
            tcwv_unc=solved_atmo.tcwv_unc,
            n_iterations=1,
            final_cost=final_cost,
            success=n_solved > 0,
            message=f"surface_driven: {n_solved}/{ny * nx} pixels solved",
            atmo_state=solved_atmo,
        )

    def solve_bundle(self, bundle: SolverInputBundle) -> SolverResult:
        return self.solve(
            bundle.toa,
            bundle.surface_prior,
            bundle.geometry,
            bundle.atmo_prior,
            bundle.rt_model,
            bundle.cloud_mask,
            bundle.bands,
            sharp_transition_mask=bundle.sharp_transition_mask,
            water_mask=bundle.water_mask,
        )
