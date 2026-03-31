"""
Multi-grid L-BFGS-B solver for aerosol retrieval.

This module implements a multi-grid optimization approach that solves
for atmospheric parameters (AOT, TCWV) at progressively finer resolutions.

The multi-grid approach:
1. Start at coarse resolution (fast, captures large-scale structure)
2. Progressively refine to finer resolutions
3. Use solution from coarser level as initial guess for finer level
4. Efficient handling of the ill-posed inversion problem

The solver minimizes the total cost function:
    J = J_obs + J_prior + J_smooth

using the L-BFGS-B quasi-Newton method with box constraints.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr
from scipy import optimize

from siac._rust_compat import (
    evaluate_grid_search_cost_cube_with_provider_qa,
    interpolate_to_fine_grid,
    quadratic_refine_grid_search_qa,
    remap_to_coarse_grid,
)
from siac.algorithms.solver.cost import CostFunction, CostFunctionConfig
from siac.domain.protocols import RTModelBackend
from siac.runtime import AtmosphericState, GeometryAngles, SurfacePrior
from siac.runtime.models import copy_spatial_metadata_like

if TYPE_CHECKING:
    from siac.domain import SensorBand

logger = logging.getLogger(__name__)


def build_solver_valid_mask(
    cloud_mask: xr.DataArray,
    toa: xr.DataArray,
    surface_prior: SurfacePrior,
) -> xr.DataArray:
    """Build the same valid-pixel mask used by the aerosol solver."""
    valid = ~cloud_mask.values.astype(bool)

    if toa.ndim == 3:
        valid = valid & np.all((toa.values > 0) & (toa.values < 1), axis=0)
    else:
        valid = valid & (toa.values > 0) & (toa.values < 1)

    if surface_prior.mask is not None:
        surface_mask = surface_prior.mask.values
        if surface_mask.ndim == 3:
            surface_mask = np.all(surface_mask, axis=0)
        valid = valid & surface_mask.astype(bool)

    return xr.DataArray(valid, dims=cloud_mask.dims, coords=cloud_mask.coords)


@dataclass
class MultiGridConfig:
    """Configuration for multi-grid solver."""

    # Grid levels (from coarse to fine)
    n_levels: int = 6

    # Minimum grid size (at coarsest level)
    min_grid_size: int = 8

    # Maximum iterations per level
    max_iter_per_level: int = 300

    # L-BFGS-B parameters
    maxcor: int = 30  # Memory for L-BFGS
    gtol: float = 1e-2  # Gradient tolerance
    ftol: float = 1e-7 * np.finfo(float).eps  # Function tolerance

    # Parameter bounds
    aot_bounds: tuple[float, float] = (0.001, 2.5)
    tcwv_bounds: tuple[float, float] = (0.0, 8.0)

    # Cost function config
    aot_gamma: float = 10.0
    tcwv_gamma: float = 5.0

    # Convergence threshold for early stopping
    rel_tol: float = 1e-4

    # Alternative solver path when RT Jacobians are unavailable.
    use_grid_search_when_no_jacobian: bool = True
    grid_search_aot_points: int = 11
    grid_search_tcwv_points: int = 11


@dataclass
class SolverResult:
    """Result from multi-grid solver."""

    aot: xr.DataArray
    tcwv: xr.DataArray
    aot_unc: xr.DataArray
    tcwv_unc: xr.DataArray
    n_iterations: int
    final_cost: float
    success: bool
    message: str
    qa: xr.Dataset | None = None
    level_history: list[dict] = field(default_factory=list)


class MultiGridSolver:
    """
    Multi-grid L-BFGS-B solver for atmospheric parameter retrieval.

    Implements the AerosolSolver protocol.

    The solver progressively refines the solution from coarse to fine
    resolutions, using the coarse solution as initial guess for the
    finer levels. This approach is efficient for the ill-posed
    atmospheric inversion problem.
    """

    def __init__(self, config: MultiGridConfig | None = None):
        """
        Initialize solver.

        Args:
            config: Solver configuration
        """
        self.config = config or MultiGridConfig()

    def solve(
        self,
        toa: xr.DataArray,
        surface_prior: SurfacePrior,
        geometry: GeometryAngles,
        atmo_prior: AtmosphericState,
        rt_model: Any,  # RTModelBackend protocol
        cloud_mask: xr.DataArray,
        bands: list[SensorBand],
    ) -> SolverResult:
        """
        Solve for atmospheric parameters (AOT, TCWV).

        Args:
            toa: Top-of-atmosphere reflectance
            surface_prior: Surface reflectance prior from BRDF
            geometry: Viewing geometry
            atmo_prior: Prior atmospheric state from CAMS/MERRA-2
            rt_model: Radiative transfer model backend
            cloud_mask: Boolean mask (True = cloudy)
            bands: List of sensor bands to use

        Returns:
            SolverResult with solved AOT, TCWV and diagnostics.
        """
        if not isinstance(rt_model, RTModelBackend):
            raise TypeError(
                f"rt_model must implement RTModelBackend protocol, "
                f"got {type(rt_model).__name__}"
            )
        # Get image dimensions
        full_shape = self._get_shape(cloud_mask)
        logger.info(f"Starting multi-grid solver for {full_shape} image")

        # Create valid pixel mask (exclude clouds and invalid data)
        mask = self._create_mask(cloud_mask, toa, surface_prior)
        n_valid = int(np.count_nonzero(mask.values))
        logger.info(f"Valid pixels: {n_valid} ({100*n_valid/mask.size:.1f}%)")

        if n_valid == 0:
            logger.warning(
                "No valid pixels after cloud/quality masking; "
                "returning atmospheric prior as solver output."
            )
            template = cloud_mask
            return SolverResult(
                aot=xr.DataArray(atmo_prior.aot.values, dims=template.dims, coords=template.coords),
                tcwv=xr.DataArray(atmo_prior.tcwv.values, dims=template.dims, coords=template.coords),
                aot_unc=xr.DataArray(atmo_prior.aot_unc.values, dims=template.dims, coords=template.coords),
                tcwv_unc=xr.DataArray(atmo_prior.tcwv_unc.values, dims=template.dims, coords=template.coords),
                n_iterations=0,
                final_cost=float("nan"),
                success=False,
                message="No valid pixels after cloud/quality masking",
                qa=None,
                level_history=[],
            )

        has_rt_jacobian = self._rt_model_supports_jacobian(rt_model)
        use_grid_search = (
            self.config.use_grid_search_when_no_jacobian and not has_rt_jacobian
        )
        if use_grid_search:
            logger.info(
                "RT backend does not provide Jacobians; using single-level "
                "grid-search + quadratic-fit solver at full resolution."
            )
            grid_shapes = [full_shape]
        else:
            grid_shapes = self._compute_grid_levels(full_shape)

        logger.info(f"Grid levels: {grid_shapes}")

        # Initialize solution with prior
        aot = atmo_prior.aot.values.copy()
        tcwv = atmo_prior.tcwv.values.copy()
        aot_unc_final = np.maximum(atmo_prior.aot_unc.values.copy(), 0.02)
        tcwv_unc_final = np.maximum(atmo_prior.tcwv_unc.values.copy(), 0.1)

        # Track history
        level_history = []
        total_iterations = 0
        final_cost = float("nan")
        final_success = True
        final_message = "ok"
        cost_func_last: CostFunction | None = None
        final_invalid_mask: np.ndarray | None = None
        final_zero_obs_mask: np.ndarray | None = None

        # Multi-grid solve from coarse to fine
        for level, shape in enumerate(grid_shapes):
            logger.info(f"Level {level}: {shape}")

            # Create cost function at this resolution
            cost_config = CostFunctionConfig(
                aot_gamma=self.config.aot_gamma,
                tcwv_gamma=self.config.tcwv_gamma,
                aot_min=self.config.aot_bounds[0],
                aot_max=self.config.aot_bounds[1],
                tcwv_min=self.config.tcwv_bounds[0],
                tcwv_max=self.config.tcwv_bounds[1],
            )

            # Resample data to current grid
            toa_level = self._resample_to_grid(toa, shape)
            mask_level = self._resample_mask_to_grid(mask, shape)
            geometry_level = self._resample_geometry_to_grid(geometry, shape)
            atmo_prior_level = self._resample_atmo_to_grid(atmo_prior, shape)
            surface_prior_level = self._resample_surface_prior_to_grid(surface_prior, shape)

            # Resample current solution to this grid
            aot_level = self._resample_field(aot, shape)
            tcwv_level = self._resample_field(tcwv, shape)

            if use_grid_search:
                (
                    aot_solved,
                    tcwv_solved,
                    aot_unc_level,
                    tcwv_unc_level,
                    level_diag,
                ) = self._solve_level_grid_search(
                    toa=toa_level,
                    surface_prior=surface_prior_level,
                    geometry=geometry_level,
                    atmo_prior=atmo_prior_level,
                    rt_model=rt_model,
                    mask=mask_level,
                    bands=bands,
                    cost_config=cost_config,
                )
                level_iterations = int(level_diag["evaluations"])
                level_cost = float(level_diag["cost"])
                valid_pixels = int(level_diag.get("valid_pixels", 0))
                invalid_pixels = int(level_diag.get("qa_invalid_pixels", 0))
                final_invalid_mask = self._resample_boolean_mask(
                    np.asarray(level_diag["qa_invalid_mask"], dtype=bool),
                    full_shape,
                )
                final_zero_obs_mask = self._resample_boolean_mask(
                    np.asarray(level_diag["qa_zero_obs_mask"], dtype=bool),
                    full_shape,
                )
                level_success = valid_pixels > 0 and invalid_pixels < valid_pixels
                level_message = (
                    "grid-search "
                    f"invalid={invalid_pixels} "
                    f"zero_obs={int(level_diag.get('qa_zero_obs_pixels', 0))} "
                    f"boundary={int(level_diag.get('qa_boundary_pixels', 0))} "
                    f"lower_aot_boundary={int(level_diag.get('qa_lower_aot_boundary_pixels', 0))}"
                )
                level_diag = {
                    key: value
                    for key, value in level_diag.items()
                    if key not in {"qa_invalid_mask", "qa_zero_obs_mask"}
                }
            else:
                # Create cost function
                cost_func = CostFunction(
                    toa=toa_level,
                    surface_prior=surface_prior_level,
                    geometry=geometry_level,
                    atmo_prior=atmo_prior_level,
                    rt_model=rt_model,
                    bands=bands,
                    mask=mask_level,
                    config=cost_config,
                )
                cost_func_last = cost_func

                # Optimize at this level
                aot_solved, tcwv_solved, result = self._optimize_level(
                    cost_func, aot_level, tcwv_level, level
                )
                level_iterations = int(result.nit)
                level_cost = float(result.fun)
                level_success = bool(result.success)
                level_message = str(result.message)

            # Update solution
            aot = self._resample_field(aot_solved, full_shape)
            tcwv = self._resample_field(tcwv_solved, full_shape)
            if use_grid_search:
                aot_unc_final = self._resample_field(aot_unc_level, full_shape)
                tcwv_unc_final = self._resample_field(tcwv_unc_level, full_shape)

            # Record history
            history_entry = {
                "level": level,
                "shape": shape,
                "iterations": level_iterations,
                "cost": level_cost,
                "success": level_success,
                "method": "grid_search" if use_grid_search else "lbfgsb",
            }
            if use_grid_search:
                history_entry.update(level_diag)
            level_history.append(history_entry)

            total_iterations += level_iterations
            logger.info(
                f"Level {level} complete: iter={level_iterations}, "
                f"cost={level_cost:.2e}, success={level_success}"
            )
            final_cost = level_cost
            final_success = level_success
            final_message = level_message

        # Compute uncertainties
        if use_grid_search:
            aot_unc, tcwv_unc = aot_unc_final, tcwv_unc_final
        else:
            if cost_func_last is None:
                raise RuntimeError("Internal solver error: missing cost function for uncertainty estimation")
            aot_unc, tcwv_unc = self._estimate_uncertainties(
                aot, tcwv, atmo_prior, cost_func_last
            )

        # Create result
        template = cloud_mask
        qa = self._build_solver_qa_dataset(
            template=template,
            valid_mask=mask.values.astype(bool),
            aot=aot,
            tcwv=tcwv,
            invalid_mask=final_invalid_mask,
            zero_obs_mask=final_zero_obs_mask,
        )
        if level_history:
            level_history[-1].update(self._summarize_solver_qa(qa))
        result = SolverResult(
            aot=xr.DataArray(aot, dims=template.dims, coords=template.coords),
            tcwv=xr.DataArray(tcwv, dims=template.dims, coords=template.coords),
            aot_unc=xr.DataArray(aot_unc, dims=template.dims, coords=template.coords),
            tcwv_unc=xr.DataArray(tcwv_unc, dims=template.dims, coords=template.coords),
            n_iterations=total_iterations,
            final_cost=final_cost,
            success=final_success,
            message=final_message,
            qa=qa,
            level_history=level_history,
        )

        logger.info(
            f"Solver complete: AOT={aot.mean():.3f}, TCWV={tcwv.mean():.2f}"
        )

        return result

    @staticmethod
    def _rt_model_supports_jacobian(rt_model: Any) -> bool:
        """Return whether backend can provide per-pixel RT Jacobians."""
        fn = getattr(rt_model, "supports_jacobian", None)
        if callable(fn):
            try:
                return bool(fn())
            except Exception:
                return False
        return False

    @staticmethod
    def _bound_tolerance(bounds: tuple[float, float]) -> float:
        span = max(float(bounds[1]) - float(bounds[0]), 0.0)
        return max(1.0e-6, span * 1.0e-5)

    def _boundary_hit_masks(
        self,
        values: np.ndarray,
        bounds: tuple[float, float],
        valid_mask: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        atol = self._bound_tolerance(bounds)
        finite = np.isfinite(values)
        lower = valid_mask & finite & np.isclose(values, float(bounds[0]), rtol=0.0, atol=atol)
        upper = valid_mask & finite & np.isclose(values, float(bounds[1]), rtol=0.0, atol=atol)
        return lower, upper

    def _resample_boolean_mask(
        self,
        mask: np.ndarray,
        target_shape: tuple[int, int],
    ) -> np.ndarray:
        resampled = self._resample_mask_to_grid(
            xr.DataArray(np.asarray(mask, dtype=bool), dims=["y", "x"]),
            target_shape,
        )
        return np.asarray(resampled.values, dtype=bool)

    def _mask_to_data_array(
        self,
        mask: np.ndarray,
        template: xr.DataArray,
    ) -> xr.DataArray:
        coords = {
            dim: template.coords[dim]
            for dim in template.dims
            if dim in template.coords
        }
        out = xr.DataArray(np.asarray(mask, dtype=bool), dims=template.dims, coords=coords)
        return copy_spatial_metadata_like(out, template)

    def _build_solver_qa_dataset(
        self,
        *,
        template: xr.DataArray,
        valid_mask: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
        invalid_mask: np.ndarray | None,
        zero_obs_mask: np.ndarray | None,
    ) -> xr.Dataset:
        valid = np.asarray(valid_mask, dtype=bool)
        invalid = np.zeros_like(valid, dtype=bool) if invalid_mask is None else np.asarray(invalid_mask, dtype=bool) & valid
        zero_obs = np.zeros_like(valid, dtype=bool) if zero_obs_mask is None else np.asarray(zero_obs_mask, dtype=bool) & valid
        aot_lower, aot_upper = self._boundary_hit_masks(
            np.asarray(aot, dtype=np.float32),
            self.config.aot_bounds,
            valid,
        )
        tcwv_lower, tcwv_upper = self._boundary_hit_masks(
            np.asarray(tcwv, dtype=np.float32),
            self.config.tcwv_bounds,
            valid,
        )
        parameter_boundary = aot_lower | aot_upper | tcwv_lower | tcwv_upper
        low_quality = invalid | parameter_boundary

        return xr.Dataset(
            {
                "invalid_retrieval": self._mask_to_data_array(invalid, template),
                "zero_obs_support": self._mask_to_data_array(zero_obs, template),
                "aot_lower_boundary": self._mask_to_data_array(aot_lower, template),
                "aot_upper_boundary": self._mask_to_data_array(aot_upper, template),
                "tcwv_lower_boundary": self._mask_to_data_array(tcwv_lower, template),
                "tcwv_upper_boundary": self._mask_to_data_array(tcwv_upper, template),
                "parameter_boundary": self._mask_to_data_array(parameter_boundary, template),
                "low_quality": self._mask_to_data_array(low_quality, template),
            }
        )

    @staticmethod
    def _summarize_solver_qa(qa: xr.Dataset) -> dict[str, float]:
        return {
            "qa_final_invalid_pixels": float(np.count_nonzero(np.asarray(qa["invalid_retrieval"].values, dtype=bool))),
            "qa_final_zero_obs_pixels": float(np.count_nonzero(np.asarray(qa["zero_obs_support"].values, dtype=bool))),
            "qa_final_aot_lower_boundary_pixels": float(np.count_nonzero(np.asarray(qa["aot_lower_boundary"].values, dtype=bool))),
            "qa_final_aot_upper_boundary_pixels": float(np.count_nonzero(np.asarray(qa["aot_upper_boundary"].values, dtype=bool))),
            "qa_final_tcwv_lower_boundary_pixels": float(np.count_nonzero(np.asarray(qa["tcwv_lower_boundary"].values, dtype=bool))),
            "qa_final_tcwv_upper_boundary_pixels": float(np.count_nonzero(np.asarray(qa["tcwv_upper_boundary"].values, dtype=bool))),
            "qa_final_parameter_boundary_pixels": float(np.count_nonzero(np.asarray(qa["parameter_boundary"].values, dtype=bool))),
            "qa_final_low_quality_pixels": float(np.count_nonzero(np.asarray(qa["low_quality"].values, dtype=bool))),
        }

    @staticmethod
    def _compute_band_weights(bands: list[SensorBand], power: float) -> np.ndarray:
        """Compute normalized spectral weights used in observation cost."""
        wavelengths = np.array([b.center_wavelength for b in bands], dtype=np.float32)
        wl_um = np.maximum(wavelengths / 1000.0, 1e-6)
        weights = wl_um ** power
        total = float(np.sum(weights))
        if total <= 0:
            return np.full(len(bands), 1.0 / max(len(bands), 1), dtype=np.float32)
        return (weights / total).astype(np.float32)

    def _solve_level_grid_search(
        self,
        *,
        toa: xr.DataArray,
        surface_prior: SurfacePrior,
        geometry: GeometryAngles,
        atmo_prior: AtmosphericState,
        rt_model: Any,
        mask: xr.DataArray,
        bands: list[SensorBand],
        cost_config: CostFunctionConfig,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
        """
        Solve one grid level via exhaustive AOT/TCWV candidates + local quadratic fit.

        This path avoids Jacobian-based optimization and returns per-pixel
        uncertainty from the fitted local Hessian.
        """
        shape = self._get_shape(mask)
        valid_mask = mask.values.astype(bool)

        n_aot = max(3, int(self.config.grid_search_aot_points))
        n_tcwv = max(3, int(self.config.grid_search_tcwv_points))
        aot_axis = np.linspace(
            self.config.aot_bounds[0], self.config.aot_bounds[1], n_aot, dtype=np.float32
        )
        tcwv_axis = np.linspace(
            self.config.tcwv_bounds[0], self.config.tcwv_bounds[1], n_tcwv, dtype=np.float32
        )
        band_weights = self._compute_band_weights(
            bands, power=cost_config.band_weight_power
        )

        # Priors / uncertainties (with floors) for per-pixel prior term.
        aot_prior = atmo_prior.aot.values.astype(np.float32)
        tcwv_prior = atmo_prior.tcwv.values.astype(np.float32)
        aot_prior_unc = np.maximum(atmo_prior.aot_unc.values.astype(np.float32), cost_config.min_aot_unc)
        tcwv_prior_unc = np.maximum(atmo_prior.tcwv_unc.values.astype(np.float32), cost_config.min_tcwv_unc)

        n_bands = len(bands)
        toa_values = toa.values.astype(np.float32)
        if toa_values.ndim == 2:
            toa_values = toa_values[np.newaxis, ...]
        if toa_values.shape != (n_bands, *shape):
            raise ValueError(
                f"TOA shape {toa_values.shape} incompatible with {n_bands} bands and grid {shape}"
            )
        toa_values = np.ascontiguousarray(toa_values, dtype=np.float32)

        boa_prior = surface_prior.boa.values.astype(np.float32)
        if boa_prior.ndim == 2:
            boa_prior = np.broadcast_to(boa_prior, (n_bands, *shape))
        if boa_prior.ndim != 3 or boa_prior.shape[-2:] != shape:
            raise ValueError(
                f"BOA prior shape {boa_prior.shape} incompatible with {n_bands} bands and grid {shape}"
            )
        if boa_prior.shape[0] < n_bands:
            raise ValueError(
                f"BOA prior has {boa_prior.shape[0]} bands, needs at least {n_bands}"
            )
        if boa_prior.shape[0] > n_bands:
            boa_prior = boa_prior[:n_bands]
        boa_prior = np.ascontiguousarray(boa_prior, dtype=np.float32)

        boa_unc = np.maximum(surface_prior.boa_unc.values.astype(np.float32), cost_config.min_boa_unc)
        if boa_unc.ndim == 2:
            boa_unc = np.broadcast_to(boa_unc, (n_bands, *shape))
        if boa_unc.ndim != 3 or boa_unc.shape[-2:] != shape:
            raise ValueError(
                f"BOA uncertainty shape {boa_unc.shape} incompatible with {n_bands} bands and grid {shape}"
            )
        if boa_unc.shape[0] < n_bands:
            raise ValueError(
                f"BOA uncertainty has {boa_unc.shape[0]} bands, needs at least {n_bands}"
            )
        if boa_unc.shape[0] > n_bands:
            boa_unc = boa_unc[:n_bands]
        boa_unc = np.ascontiguousarray(boa_unc, dtype=np.float32)

        xap_stack = np.empty((n_bands, shape[0], shape[1]), dtype=np.float32)
        xbp_stack = np.empty((n_bands, shape[0], shape[1]), dtype=np.float32)
        xcp_stack = np.empty((n_bands, shape[0], shape[1]), dtype=np.float32)

        def _candidate_coeff_provider(
            aot_val: float, tcwv_val: float
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            aot_field = xr.full_like(atmo_prior.aot, fill_value=float(aot_val), dtype=np.float32)
            tcwv_field = xr.full_like(atmo_prior.tcwv, fill_value=float(tcwv_val), dtype=np.float32)
            atmo_state = atmo_prior.with_updated_aot_tcwv(aot=aot_field, tcwv=tcwv_field)

            for ib, band in enumerate(bands):
                coeffs = rt_model.compute_coefficients(
                    geometry, atmo_state, band, compute_jacobian=False
                )
                xap_stack[ib] = np.asarray(coeffs.xap.values, dtype=np.float32)
                xbp_stack[ib] = np.asarray(coeffs.xbp.values, dtype=np.float32)
                xcp_stack[ib] = np.asarray(coeffs.xcp.values, dtype=np.float32)
            return xap_stack, xbp_stack, xcp_stack

        costs_raw, obs_counts_raw = evaluate_grid_search_cost_cube_with_provider_qa(
            _candidate_coeff_provider,
            aot_axis.astype(np.float32, copy=False),
            tcwv_axis.astype(np.float32, copy=False),
            toa_values,
            boa_prior,
            boa_unc,
            band_weights,
            valid_mask.astype(bool, copy=False),
            aot_prior,
            tcwv_prior,
            aot_prior_unc,
            tcwv_prior_unc,
        )
        costs = np.asarray(costs_raw, dtype=np.float32)
        obs_counts = np.asarray(obs_counts_raw, dtype=np.uint16)

        (
            aot_best,
            tcwv_best,
            aot_unc,
            tcwv_unc,
            invalid_mask,
            boundary_mask,
            lower_aot_boundary_mask,
            zero_obs_mask,
        ) = quadratic_refine_grid_search_qa(
            costs.astype(np.float32, copy=False),
            obs_counts.astype(np.uint16, copy=False),
            aot_axis.astype(np.float32, copy=False),
            tcwv_axis.astype(np.float32, copy=False),
            valid_mask.astype(bool, copy=False),
        )
        aot_best = np.asarray(aot_best, dtype=np.float32)
        tcwv_best = np.asarray(tcwv_best, dtype=np.float32)
        aot_unc = np.asarray(aot_unc, dtype=np.float32)
        tcwv_unc = np.asarray(tcwv_unc, dtype=np.float32)
        invalid_mask = np.asarray(invalid_mask, dtype=bool)
        boundary_mask = np.asarray(boundary_mask, dtype=bool)
        lower_aot_boundary_mask = np.asarray(lower_aot_boundary_mask, dtype=bool)
        zero_obs_mask = np.asarray(zero_obs_mask, dtype=bool)

        invalid_floor_pixels = int(
            np.count_nonzero(
                invalid_mask
                & valid_mask
                & np.isclose(aot_best, float(aot_axis[0]), rtol=0.0, atol=1.0e-6)
            )
        )
        prior_floor_pixels = int(
            np.count_nonzero(
                invalid_mask
                & valid_mask
                & np.isclose(aot_prior, float(aot_axis[0]), rtol=0.0, atol=1.0e-6)
            )
        )
        if np.any(invalid_mask):
            aot_best = np.where(invalid_mask, aot_prior, aot_best).astype(np.float32, copy=False)
            tcwv_best = np.where(invalid_mask, tcwv_prior, tcwv_best).astype(np.float32, copy=False)
            aot_unc = np.where(invalid_mask, aot_prior_unc, aot_unc).astype(np.float32, copy=False)
            tcwv_unc = np.where(invalid_mask, tcwv_prior_unc, tcwv_unc).astype(np.float32, copy=False)

        flat = costs.reshape(n_aot * n_tcwv, -1)
        safe_flat = np.where(np.isfinite(flat), flat, np.inf)
        best_flat = np.argmin(safe_flat, axis=0)
        selected_costs = safe_flat[best_flat, np.arange(safe_flat.shape[1])]
        supported = (valid_mask & ~invalid_mask).reshape(-1)
        mean_cost = (
            float(np.mean(selected_costs[supported]))
            if np.any(supported)
            else float("nan")
        )
        valid_pixels = int(np.count_nonzero(valid_mask))
        invalid_pixels = int(np.count_nonzero(invalid_mask & valid_mask))
        zero_obs_pixels = int(np.count_nonzero(zero_obs_mask & valid_mask))
        boundary_pixels = int(np.count_nonzero(boundary_mask & valid_mask))
        lower_aot_boundary_pixels = int(
            np.count_nonzero(lower_aot_boundary_mask & valid_mask)
        )
        if invalid_pixels:
            logger.warning(
                "Grid-search QA flagged %d/%d valid pixels as invalid; %d had zero observation support and %d collapsed to the AOT floor.",
                invalid_pixels,
                valid_pixels,
                zero_obs_pixels,
                invalid_floor_pixels,
            )
        return (
            aot_best.astype(np.float32),
            tcwv_best.astype(np.float32),
            aot_unc.astype(np.float32),
            tcwv_unc.astype(np.float32),
            {
                "cost": mean_cost,
                "evaluations": float(n_aot * n_tcwv),
                "valid_pixels": float(valid_pixels),
                "qa_invalid_pixels": float(invalid_pixels),
                "qa_zero_obs_pixels": float(zero_obs_pixels),
                "qa_boundary_pixels": float(boundary_pixels),
                "qa_lower_aot_boundary_pixels": float(lower_aot_boundary_pixels),
                "qa_invalid_floor_pixels": float(invalid_floor_pixels),
                "qa_prior_floor_pixels": float(prior_floor_pixels),
                "qa_invalid_mask": invalid_mask.astype(bool, copy=False),
                "qa_zero_obs_mask": zero_obs_mask.astype(bool, copy=False),
            },
        )

    def _get_shape(self, arr: xr.DataArray) -> tuple[int, int]:
        """Get (ny, nx) shape from DataArray."""
        if "y" in arr.dims and "x" in arr.dims:
            return (arr.sizes["y"], arr.sizes["x"])
        return arr.shape

    def _create_mask(
        self,
        cloud_mask: xr.DataArray,
        toa: xr.DataArray,
        surface_prior: SurfacePrior,
    ) -> xr.DataArray:
        """Create combined valid pixel mask."""
        return build_solver_valid_mask(cloud_mask, toa, surface_prior)

    def _compute_grid_levels(
        self, full_shape: tuple[int, int]
    ) -> list[tuple[int, int]]:
        """Compute grid shapes for each level (coarse to fine)."""
        ny, nx = full_shape
        min_size = self.config.min_grid_size

        # Compute number of levels based on image size
        ratio = min(ny, nx) / min_size
        max_levels = min(
            self.config.n_levels,
            max(1, int(np.log2(max(ratio, 1.0))) + 1),
        )

        shapes = []
        for i in range(max_levels):
            # Exponential spacing from min to full
            scale = 2 ** (max_levels - 1 - i)
            shape_y = max(min_size, ny // scale)
            shape_x = max(min_size, nx // scale)
            shapes.append((shape_y, shape_x))

        return shapes

    def _resample_field(
        self, field: np.ndarray, target_shape: tuple[int, int]
    ) -> np.ndarray:
        """Resample 2D field to target shape."""
        if field.shape == target_shape:
            return field

        data = np.ascontiguousarray(field, dtype=np.float64)
        if target_shape[0] < field.shape[0]:
            return np.asarray(remap_to_coarse_grid(data, target_shape[0], target_shape[1]))
        return np.asarray(interpolate_to_fine_grid(data, target_shape[0], target_shape[1]))

    def _resample_to_grid(
        self, data: xr.DataArray, shape: tuple[int, int]
    ) -> xr.DataArray:
        """Resample DataArray to target grid."""
        if data.ndim == 3:
            # Multi-band
            result = np.stack([
                self._resample_field(data.values[i], shape)
                for i in range(data.shape[0])
            ])
            return xr.DataArray(result, dims=["band", "y", "x"])
        else:
            result = self._resample_field(data.values, shape)
            return xr.DataArray(result, dims=["y", "x"])

    def _resample_mask_to_grid(
        self, mask: xr.DataArray, shape: tuple[int, int]
    ) -> xr.DataArray:
        """Resample mask to target grid using max pooling."""
        if mask.shape == shape:
            return mask

        if shape[0] <= mask.shape[0] and shape[1] <= mask.shape[1]:
            # Match the center-based coarse-grid assignment used for numeric
            # field remapping so coarse validity covers the same source pixels.
            result = np.ones(shape, dtype=bool)
            src = np.asarray(mask.values, dtype=bool)
            for iy in range(src.shape[0]):
                dst_y = min(((2 * iy + 1) * shape[0]) // (2 * src.shape[0]), shape[0] - 1)
                for ix in range(src.shape[1]):
                    dst_x = min(((2 * ix + 1) * shape[1]) // (2 * src.shape[1]), shape[1] - 1)
                    result[dst_y, dst_x] &= bool(src[iy, ix])
            return xr.DataArray(result, dims=["y", "x"])

        upsampled = np.asarray(
            interpolate_to_fine_grid(
                np.asarray(mask.values, dtype=np.float64),
                shape[0],
                shape[1],
            )
        )
        return xr.DataArray(upsampled > 0.5, dims=["y", "x"])

    def _resample_geometry_to_grid(
        self, geometry: GeometryAngles, shape: tuple[int, int]
    ) -> GeometryAngles:
        """Resample geometry to target grid."""
        return GeometryAngles(
            sza=self._resample_to_grid(geometry.sza, shape),
            saa=self._resample_to_grid(geometry.saa, shape),
            vza=self._resample_to_grid(geometry.vza, shape),
            vaa=self._resample_to_grid(geometry.vaa, shape),
        )

    def _resample_atmo_to_grid(
        self, atmo: AtmosphericState, shape: tuple[int, int]
    ) -> AtmosphericState:
        """Resample atmospheric state to target grid."""
        return AtmosphericState(
            aot=self._resample_to_grid(atmo.aot, shape),
            tcwv=self._resample_to_grid(atmo.tcwv, shape),
            tco3=self._resample_to_grid(atmo.tco3, shape),
            aot_unc=self._resample_to_grid(atmo.aot_unc, shape),
            tcwv_unc=self._resample_to_grid(atmo.tcwv_unc, shape),
            tco3_unc=self._resample_to_grid(atmo.tco3_unc, shape),
            elevation=self._resample_to_grid(atmo.elevation, shape),
        )

    def _resample_surface_prior_to_grid(
        self, prior: SurfacePrior, shape: tuple[int, int]
    ) -> SurfacePrior:
        """Resample surface prior to target grid."""
        return SurfacePrior(
            boa=self._resample_to_grid(prior.boa, shape),
            boa_unc=self._resample_to_grid(prior.boa_unc, shape),
            kernels=prior.kernels,  # Keep original kernels
            mask=self._resample_mask_to_grid(prior.mask, shape),
        )

    def _optimize_level(
        self,
        cost_func: CostFunction,
        aot_init: np.ndarray,
        tcwv_init: np.ndarray,
        level: int,
    ) -> tuple[np.ndarray, np.ndarray, optimize.OptimizeResult]:
        """
        Run L-BFGS-B optimization at one grid level.

        Args:
            cost_func: Cost function for this level
            aot_init: Initial AOT guess
            tcwv_init: Initial TCWV guess
            level: Grid level index

        Returns:
            Tuple of (solved_aot, solved_tcwv, optimize_result)
        """
        # Pack initial guess
        p0 = np.concatenate([aot_init.ravel(), tcwv_init.ravel()])

        # Setup bounds
        n = aot_init.size
        bounds_aot = [self.config.aot_bounds] * n
        bounds_tcwv = [self.config.tcwv_bounds] * n
        bounds = bounds_aot + bounds_tcwv

        # Adjust parameters based on level
        maxiter = self.config.max_iter_per_level
        if level < 2:
            maxiter = min(maxiter, 100)  # Fewer iterations at coarse levels

        # Run optimization
        result = optimize.minimize(
            cost_func,
            p0,
            jac=True,
            method="L-BFGS-B",
            bounds=bounds,
            options={
                "disp": False,
                "maxcor": self.config.maxcor,
                "gtol": self.config.gtol,
                "ftol": self.config.ftol,
                "maxiter": maxiter,
                "maxls": 100,
            },
        )

        # Unpack solution
        aot_solved = result.x[:n].reshape(aot_init.shape)
        tcwv_solved = result.x[n:].reshape(tcwv_init.shape)

        return aot_solved, tcwv_solved, result

    def _estimate_uncertainties(
        self,
        aot: np.ndarray,
        tcwv: np.ndarray,
        _atmo_prior: AtmosphericState,
        cost_func: CostFunction,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Estimate posterior uncertainties using prior and smoothness info.

        A simplified uncertainty estimate based on the prior uncertainty
        scaled by the smoothness regularization.
        """
        # Prior uncertainties

        # Scale by smoothness gamma (higher gamma = more smoothing = higher unc)
        gamma_aot = cost_func.config.aot_gamma
        gamma_tcwv = cost_func.config.tcwv_gamma

        # Compute adjustment factor from DCT (Gamma filter mean)
        lambda_vals = cost_func.lambda_smooth
        gamma_filter_aot = 1.0 / (1.0 + gamma_aot ** 2 * lambda_vals)
        gamma_filter_tcwv = 1.0 / (1.0 + gamma_tcwv ** 2 * lambda_vals)

        adj_aot = 1.0 / np.sqrt(gamma_filter_aot.mean())
        adj_tcwv = 1.0 / np.sqrt(gamma_filter_tcwv.mean())

        # Posterior uncertainty
        aot_unc = np.maximum(aot * 0.5, 0.02) * adj_aot
        tcwv_unc = np.maximum(tcwv * 0.1, 0.1) * adj_tcwv

        # Ensure shape matches
        aot_unc = self._resample_field(aot_unc, aot.shape)
        tcwv_unc = self._resample_field(tcwv_unc, tcwv.shape)

        return aot_unc, tcwv_unc


def solve_atmospheric_parameters(
    toa: xr.DataArray,
    surface_prior: SurfacePrior,
    geometry: GeometryAngles,
    atmo_prior: AtmosphericState,
    rt_model: Any,
    cloud_mask: xr.DataArray,
    bands: list[SensorBand],
    config: MultiGridConfig | None = None,
) -> SolverResult:
    """
    Convenience function to solve for atmospheric parameters.

    Args:
        toa: Top-of-atmosphere reflectance
        surface_prior: Surface reflectance prior from BRDF
        geometry: Viewing geometry
        atmo_prior: Prior atmospheric state
        rt_model: Radiative transfer model backend
        cloud_mask: Boolean mask (True = cloudy)
        bands: List of sensor bands to use
        config: Solver configuration

    Returns:
        SolverResult with solved AOT, TCWV and diagnostics.
    """
    solver = MultiGridSolver(config)
    return solver.solve(
        toa=toa,
        surface_prior=surface_prior,
        geometry=geometry,
        atmo_prior=atmo_prior,
        rt_model=rt_model,
        cloud_mask=cloud_mask,
        bands=bands,
    )
