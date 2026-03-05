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
from typing import Any

import numpy as np
import xarray as xr
from scipy import ndimage, optimize

from siac.core.protocols import RTModelBackend
from siac.core.types import (
    AtmosphericState,
    GeometryAngles,
    SensorBand,
    SurfacePrior,
)
from siac.solver.cost import CostFunction, CostFunctionConfig

logger = logging.getLogger(__name__)

# Try to import Rust optimization functions
try:
    from siac._rust import interpolate_to_fine_grid, remap_to_coarse_grid
    _HAS_RUST_OPT = True
    logger.debug("Using Rust optimization functions")
except ImportError:
    _HAS_RUST_OPT = False
    logger.debug("Rust optimization not available, using scipy fallback")


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
        n_valid = mask.sum()
        logger.info(f"Valid pixels: {n_valid} ({100*n_valid/mask.size:.1f}%)")

        # Compute grid levels
        grid_shapes = self._compute_grid_levels(full_shape)
        logger.info(f"Grid levels: {grid_shapes}")

        # Initialize solution with prior
        aot = atmo_prior.aot.values.copy()
        tcwv = atmo_prior.tcwv.values.copy()

        # Track history
        level_history = []
        total_iterations = 0

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

            # Optimize at this level
            aot_solved, tcwv_solved, result = self._optimize_level(
                cost_func, aot_level, tcwv_level, level
            )

            # Update solution
            aot = self._resample_field(aot_solved, full_shape)
            tcwv = self._resample_field(tcwv_solved, full_shape)

            # Record history
            level_history.append({
                "level": level,
                "shape": shape,
                "iterations": result.nit,
                "cost": result.fun,
                "success": result.success,
            })

            total_iterations += result.nit
            logger.info(
                f"Level {level} complete: iter={result.nit}, "
                f"cost={result.fun:.2e}, success={result.success}"
            )

        # Compute uncertainties
        aot_unc, tcwv_unc = self._estimate_uncertainties(
            aot, tcwv, atmo_prior, cost_func
        )

        # Create result
        template = cloud_mask
        result = SolverResult(
            aot=xr.DataArray(aot, dims=template.dims, coords=template.coords),
            tcwv=xr.DataArray(tcwv, dims=template.dims, coords=template.coords),
            aot_unc=xr.DataArray(aot_unc, dims=template.dims, coords=template.coords),
            tcwv_unc=xr.DataArray(tcwv_unc, dims=template.dims, coords=template.coords),
            n_iterations=total_iterations,
            final_cost=level_history[-1]["cost"],
            success=level_history[-1]["success"],
            message=str(result.message),
            level_history=level_history,
        )

        logger.info(
            f"Solver complete: AOT={aot.mean():.3f}, TCWV={tcwv.mean():.2f}"
        )

        return result

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
        # Not cloudy
        valid = ~cloud_mask.values

        # Valid TOA
        if toa.ndim == 3:
            valid = valid & np.all((toa.values > 0) & (toa.values < 1), axis=0)
        else:
            valid = valid & (toa.values > 0) & (toa.values < 1)

        # Valid surface prior
        if surface_prior.mask is not None:
            valid = valid & surface_prior.mask.values

        return xr.DataArray(valid, dims=cloud_mask.dims, coords=cloud_mask.coords)

    def _compute_grid_levels(
        self, full_shape: tuple[int, int]
    ) -> list[tuple[int, int]]:
        """Compute grid shapes for each level (coarse to fine)."""
        ny, nx = full_shape
        min_size = self.config.min_grid_size

        # Compute number of levels based on image size
        max_levels = min(
            self.config.n_levels,
            int(np.log2(min(ny, nx) / min_size)) + 1
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

        if _HAS_RUST_OPT:
            data = np.ascontiguousarray(field, dtype=np.float64)
            if target_shape[0] < field.shape[0]:
                return np.asarray(remap_to_coarse_grid(data, target_shape[0], target_shape[1]))
            else:
                return np.asarray(interpolate_to_fine_grid(data, target_shape[0], target_shape[1]))

        # Fallback: scipy zoom
        zoom_factors = (
            target_shape[0] / field.shape[0],
            target_shape[1] / field.shape[1],
        )
        return ndimage.zoom(field, zoom_factors, order=1)

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

        # Block reduce with logical AND (valid only if all sub-pixels valid)
        block_y = mask.shape[0] // shape[0]
        block_x = mask.shape[1] // shape[1]

        # Reshape and reduce
        result = np.zeros(shape, dtype=bool)
        for i in range(shape[0]):
            for j in range(shape[1]):
                block = mask.values[
                    i * block_y : (i + 1) * block_y,
                    j * block_x : (j + 1) * block_x,
                ]
                result[i, j] = block.mean() > 0.5  # Majority valid

        return xr.DataArray(result, dims=["y", "x"])

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
