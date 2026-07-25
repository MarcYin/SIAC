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
import time as _time
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Any, Literal, TypeAlias, cast

import numpy as np
import xarray as xr
from numpy import typing as npt
from scipy import optimize

from siac._rust_compat import (
    evaluate_block_grid_search_cost_cube_with_provider_qa,
    evaluate_grid_search_cost_cube_with_provider_qa,
    interpolate_to_fine_grid,
    quadratic_refine_grid_search_qa,
    remap_to_coarse_grid,
)
from siac.algorithms.solver._grid_search import (
    aggregate_block_pixel_counts,
    aggregate_valid_counts,
    build_candidate_coeff_provider,
    evaluate_candidate_cost_cube,
    fixed_axis_from_prior,
    prepare_grid_search_observations,
)
from siac.algorithms.solver._qa import build_solver_qa_dataset
from siac.algorithms.solver.cost import CostFunction, CostFunctionConfig
from siac.domain.protocols import (
    RTModelBackend,
    rt_supports_jacobian,
)
from siac.runtime import AtmosphericState, BRDFKernelWeights, GeometryAngles, SurfacePrior

if TYPE_CHECKING:
    from siac.domain import SensorBand
    from siac.runtime import SolverInputBundle

logger = logging.getLogger(__name__)

BoolArray: TypeAlias = npt.NDArray[np.bool_]
Float32Array: TypeAlias = npt.NDArray[np.float32]
Float64Array: TypeAlias = npt.NDArray[np.float64]
AtmosphericParameterName: TypeAlias = Literal["aot", "tcwv", "tco3"]
FixedAtmosphericParameter: TypeAlias = Literal["none", "aot", "tcwv"]
StageInitialState: TypeAlias = Literal["prior", "previous"]


# Second-pass high-AOD refine. The coarse log AOT axis leaves a wide gap at its
# top (1.143<->2.506 at the default n=11 over [0.001, 2.5]), so a thick-aerosol
# cost minimum inside that gap is unsampled and the retrieval pins at the
# second-highest node (LAMTO 1.52->1.14, NGHIA 1.83->1.14) even though the
# observation cost cleanly prefers the higher value. When any pixel lands above
# ``_HIGH_AOT_REFINE_MIN`` the solver re-searches a fine axis spanning the gap
# and adopts it per-pixel ONLY where it yields a strictly lower grid cost --
# self-gating, so already-resolved low/mid pixels never switch and the axis-length
# invariant of the primary pass is untouched.
_HIGH_AOT_REFINE_MIN = 0.8
_HIGH_AOT_REFINE_LO = 0.5
_HIGH_AOT_REFINE_POINTS = 9


def _log_aot_axis(lo: float, hi: float, n: int) -> np.ndarray:
    """Log-spaced AOT grid-search axis.

    A linear axis over [0.001, 2.5] put its second node at 0.25, leaving the
    entire typical AOT range (~0.01-0.2) unsampled; the grid search then picked
    the 0.001 lower bound and the quadratic refine pinned retrievals at the
    floor. Log spacing places nodes near 0.025/0.056/0.125 so low-AOT scenes
    are resolvable. The lower bound is clamped above zero for geomspace.
    """
    lo_safe = max(float(lo), 1e-4)
    return np.geomspace(lo_safe, float(hi), int(n), dtype=np.float32)


def derive_grid_search_axes(
    config: MultiGridConfig,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Return the (aot_axis, tcwv_axis) that would be used by the block-grid
    search for the given config, or None when the axes depend on the per-pixel
    prior (fixed-parameter mode) or vary across stages.

    Used by the pipeline preload step to compute the joint LUT in parallel
    with prior-fetching, before the solver itself runs. When this returns
    None, the joint LUT can't be safely precomputed — it must wait for the
    solver to inspect the prior and decide.
    """
    # Staged solvers may use different fixed_parameter per stage; each stage
    # rebuilds its own grid axes. The conservative move is to skip preload
    # there — the savings exist only in the simple single-pass case anyway.
    if getattr(config, "stages", None):
        return None
    fixed = str(getattr(config, "fixed_atmospheric_parameter", "none"))
    if fixed in {"aot", "tcwv"}:
        # When one of the parameters is fixed, its axis derives from the
        # per-pixel prior — we can't know it until the prior is in hand.
        return None

    n_aot = max(3, int(config.grid_search_aot_points))
    n_tcwv = max(3, int(config.grid_search_tcwv_points))
    aot_axis = _log_aot_axis(config.aot_bounds[0], config.aot_bounds[1], n_aot)
    tcwv_axis = np.linspace(config.tcwv_bounds[0], config.tcwv_bounds[1], n_tcwv, dtype=np.float32)
    return aot_axis, tcwv_axis


def build_solver_valid_mask(
    cloud_mask: xr.DataArray,
    toa: xr.DataArray,
    surface_prior: SurfacePrior,
    *,
    sharp_transition_mask: xr.DataArray | None = None,
    water_mask: xr.DataArray | None = None,
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
    if sharp_transition_mask is not None:
        valid = valid & ~sharp_transition_mask.values.astype(bool)
    if water_mask is not None:
        valid = valid & ~water_mask.values.astype(bool)

    return xr.DataArray(valid, dims=cloud_mask.dims, coords=cloud_mask.coords)


@dataclass(frozen=True)
class SolverStageConfig:
    """One atmospheric retrieval pass in a staged solver chain."""

    name: str = "default"
    solve: tuple[AtmosphericParameterName, ...] = ("aot", "tcwv")
    fixed: tuple[AtmosphericParameterName, ...] = ("tco3",)
    bands: tuple[str, ...] | None = None
    initial_state: StageInitialState = "previous"


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
    # ``ftol`` is intentionally set to a value below machine epsilon so that
    # L-BFGS-B convergence is governed entirely by ``gtol`` (the gradient
    # norm). REVIEW.md §1.1 #6 flagged this as ``effectively zero''; the
    # current behaviour IS gradient-only convergence by design — switching to
    # a non-zero ftol would change retrievals across the regression suite and
    # should only be done with paired numerical verification on a real scene.
    ftol: float = 1e-7 * np.finfo(float).eps

    # Parameter bounds
    aot_bounds: tuple[float, float] = (0.001, 2.5)
    tcwv_bounds: tuple[float, float] = (0.0, 8.0)

    # Cost function config
    aot_gamma: float = 10.0
    tcwv_gamma: float = 5.0

    # Band weighting power (negative values weight shorter wavelengths more)
    band_weight_power: float = -1.6

    # Pseudo-Huber transition threshold for smoothness regularization.
    smoothness_delta: float = 0.02

    # Convergence threshold for early stopping
    rel_tol: float = 1e-4

    # Alternative solver path when RT Jacobians are unavailable.
    use_grid_search_when_no_jacobian: bool = True
    grid_search_aot_points: int = 11
    grid_search_tcwv_points: int = 11
    fixed_atmospheric_parameter: FixedAtmosphericParameter = "none"

    # Solve one shared AOT/TCWV pair per NxN block in the no-Jacobian
    # grid-search path, then broadcast that block solution back to the
    # full solver grid. The same block size also controls RT coefficient
    # sampling in that path.
    quadratic_block_size: int = 1

    # Minimum fraction of pixels in each quadratic block that must have valid
    # observation and surface-prior support before the block is solved.
    quadratic_block_min_valid_fraction: float = 0.5

    # Optional sequence of solver stages. Empty uses the single AOT/TCWV solve
    # configured by the fields above.
    stages: tuple[SolverStageConfig, ...] = field(default_factory=tuple)


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
    atmo_state: AtmosphericState | None = None
    diagnostics: dict[str, Any] = field(default_factory=dict)


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
        rt_model: RTModelBackend,
        cloud_mask: xr.DataArray,
        bands: list[SensorBand],
        sharp_transition_mask: xr.DataArray | None = None,
        water_mask: xr.DataArray | None = None,
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
                f"rt_model must implement RTModelBackend protocol, got {type(rt_model).__name__}"
            )
        # Get image dimensions
        full_shape = self._get_shape(cloud_mask)
        logger.info(f"Starting multi-grid solver for {full_shape} image")

        # Create valid pixel mask (exclude clouds and invalid data)
        mask = self._create_mask(
            cloud_mask,
            toa,
            surface_prior,
            sharp_transition_mask=sharp_transition_mask,
            water_mask=water_mask,
        )
        full_no_observation_mask = ~self._observation_presence_mask(toa)
        n_valid = int(np.count_nonzero(mask.values))
        logger.info(f"Valid pixels: {n_valid} ({100 * n_valid / mask.size:.1f}%)")

        if n_valid == 0:
            logger.warning(
                "No valid pixels after cloud/quality masking; "
                "returning atmospheric prior as solver output."
            )
            template = cloud_mask
            aot_out = np.where(full_no_observation_mask, np.nan, atmo_prior.aot.values)
            tcwv_out = np.where(full_no_observation_mask, np.nan, atmo_prior.tcwv.values)
            aot_unc_out = np.where(full_no_observation_mask, np.nan, atmo_prior.aot_unc.values)
            tcwv_unc_out = np.where(full_no_observation_mask, np.nan, atmo_prior.tcwv_unc.values)
            return SolverResult(
                aot=xr.DataArray(aot_out, dims=template.dims, coords=template.coords),
                tcwv=xr.DataArray(tcwv_out, dims=template.dims, coords=template.coords),
                aot_unc=xr.DataArray(aot_unc_out, dims=template.dims, coords=template.coords),
                tcwv_unc=xr.DataArray(tcwv_unc_out, dims=template.dims, coords=template.coords),
                n_iterations=0,
                final_cost=float("nan"),
                success=False,
                message="No valid pixels after cloud/quality masking",
                qa=None,
                level_history=[],
            )

        has_rt_jacobian = self._rt_model_supports_jacobian(rt_model)
        use_grid_search = self.config.use_grid_search_when_no_jacobian and not has_rt_jacobian
        fixed_parameter = self._fixed_parameter()
        if fixed_parameter != "none":
            logger.info(
                "Fixing %s to the atmospheric prior during solver optimization.",
                fixed_parameter.upper(),
            )
        if use_grid_search:
            block_size = max(1, self.config.quadratic_block_size)
            if block_size > 1:
                logger.info(
                    "RT backend does not provide Jacobians; using "
                    "single-level grid-search + quadratic-fit block solver "
                    "(block=%d).",
                    block_size,
                )
            else:
                logger.info(
                    "RT backend does not provide Jacobians; using "
                    "single-level grid-search + quadratic-fit solver."
                )
            grid_shapes = [full_shape]
        else:
            grid_shapes = self._compute_grid_levels(full_shape)

        logger.info(f"Grid levels: {grid_shapes}")

        # Initialize solution with views; copies are made at mutation points below.
        aot = np.array(atmo_prior.aot.values, dtype=np.float32, copy=True)
        tcwv = np.array(atmo_prior.tcwv.values, dtype=np.float32, copy=True)
        aot_unc_final = np.maximum(atmo_prior.aot_unc.values, 0.02)
        tcwv_unc_final = np.maximum(atmo_prior.tcwv_unc.values, 0.1)

        # Track history
        level_history = []
        total_iterations = 0
        final_cost = float("nan")
        final_success = True
        final_message = "ok"
        cost_func_last: CostFunction | None = None
        final_invalid_mask: np.ndarray | None = None
        final_zero_obs_mask: np.ndarray | None = None
        final_insufficient_support_mask: np.ndarray | None = None
        final_fitting_cost: np.ndarray | None = None

        # Multi-grid solve from coarse to fine. The non-Jacobian grid-search
        # branch does not consume coarse solutions as initial guesses, so it
        # only runs once on the finest grid.
        _solve_t0 = _time.monotonic()
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
                band_weight_power=self.config.band_weight_power,
                smoothness_delta=self.config.smoothness_delta,
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
                solve_invalid_pixels = int(
                    level_diag.get("qa_solve_invalid_pixels", invalid_pixels)
                )
                final_invalid_mask = self._resample_boolean_mask(
                    np.asarray(level_diag["qa_invalid_mask"], dtype=bool),
                    full_shape,
                )
                final_zero_obs_mask = self._resample_boolean_mask(
                    np.asarray(level_diag["qa_zero_obs_mask"], dtype=bool),
                    full_shape,
                )
                final_insufficient_support_mask = self._resample_boolean_mask(
                    np.asarray(level_diag["qa_insufficient_support_mask"], dtype=bool),
                    full_shape,
                )
                final_fitting_cost = self._resample_field(
                    np.asarray(level_diag["qa_fitting_cost_map"], dtype=np.float32),
                    full_shape,
                ).astype(np.float32)
                level_success = valid_pixels > 0 and solve_invalid_pixels < valid_pixels
                level_message = (
                    "grid-search "
                    f"invalid={invalid_pixels} "
                    f"solve_invalid={solve_invalid_pixels} "
                    f"zero_obs={int(level_diag.get('qa_zero_obs_pixels', 0))} "
                    f"no_observation={int(level_diag.get('qa_no_observation_pixels', 0))} "
                    f"insufficient_support={int(level_diag.get('qa_insufficient_support_pixels', 0))} "
                    f"boundary={int(level_diag.get('qa_boundary_pixels', 0))} "
                    f"lower_aot_boundary={int(level_diag.get('qa_lower_aot_boundary_pixels', 0))}"
                )
                level_diag = {
                    key: value
                    for key, value in level_diag.items()
                    if key
                    not in {
                        "qa_invalid_mask",
                        "qa_zero_obs_mask",
                        "qa_insufficient_support_mask",
                        "qa_fitting_cost_map",
                    }
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

        _solve_elapsed = _time.monotonic() - _solve_t0
        logger.info(
            "Multi-grid solver complete: %d levels, %d iterations, cost=%.4g (%.2fs)",
            len(grid_shapes),
            total_iterations,
            final_cost,
            _solve_elapsed,
        )

        # Compute uncertainties
        if use_grid_search:
            aot_unc, tcwv_unc = aot_unc_final, tcwv_unc_final
        else:
            if cost_func_last is None:
                raise RuntimeError(
                    "Internal solver error: missing cost function for uncertainty estimation"
                )
            aot_unc, tcwv_unc = self._estimate_uncertainties(aot, tcwv, atmo_prior, cost_func_last)

        if fixed_parameter == "aot":
            aot = atmo_prior.aot.values.astype(np.float32, copy=True)
            aot_unc = atmo_prior.aot_unc.values.astype(np.float32, copy=True)
        elif fixed_parameter == "tcwv":
            tcwv = atmo_prior.tcwv.values.astype(np.float32, copy=True)
            tcwv_unc = atmo_prior.tcwv_unc.values.astype(np.float32, copy=True)

        gap_fill_mask = ~mask.values.astype(bool) & ~full_no_observation_mask
        if np.any(gap_fill_mask):
            trusted_mask = mask.values.astype(bool)
            if fixed_parameter != "aot":
                trusted_aot_mask = trusted_mask & np.isfinite(aot)
                if np.any(trusted_aot_mask):
                    aot_gap_filled = self._smooth_grid_search_field(
                        aot,
                        gamma=self.config.aot_gamma,
                        delta=self.config.smoothness_delta,
                        n_iter=40,
                        trusted_mask=trusted_aot_mask,
                    )
                    aot = np.where(gap_fill_mask, aot_gap_filled, aot).astype(
                        np.float32, copy=False
                    )
                    aot_unc = np.where(
                        gap_fill_mask,
                        np.maximum(atmo_prior.aot_unc.values.astype(np.float32), np.float32(0.1)),
                        aot_unc,
                    ).astype(np.float32, copy=False)
            if fixed_parameter != "tcwv":
                trusted_tcwv_mask = trusted_mask & np.isfinite(tcwv)
                if np.any(trusted_tcwv_mask):
                    tcwv_gap_filled = self._smooth_grid_search_field(
                        tcwv,
                        gamma=self.config.tcwv_gamma,
                        delta=self.config.smoothness_delta,
                        n_iter=40,
                        trusted_mask=trusted_tcwv_mask,
                    )
                    tcwv = np.where(gap_fill_mask, tcwv_gap_filled, tcwv).astype(
                        np.float32, copy=False
                    )
                    tcwv_unc = np.where(
                        gap_fill_mask,
                        np.maximum(atmo_prior.tcwv_unc.values.astype(np.float32), np.float32(0.5)),
                        tcwv_unc,
                    ).astype(np.float32, copy=False)

        if np.any(full_no_observation_mask):
            aot = np.where(full_no_observation_mask, np.nan, aot).astype(np.float32, copy=False)
            tcwv = np.where(full_no_observation_mask, np.nan, tcwv).astype(np.float32, copy=False)
            aot_unc = np.where(full_no_observation_mask, np.nan, aot_unc).astype(
                np.float32, copy=False
            )
            tcwv_unc = np.where(full_no_observation_mask, np.nan, tcwv_unc).astype(
                np.float32, copy=False
            )
            final_zero_obs_mask = (
                full_no_observation_mask.copy()
                if final_zero_obs_mask is None
                else final_zero_obs_mask | full_no_observation_mask
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
            insufficient_support_mask=final_insufficient_support_mask,
            no_observation_mask=full_no_observation_mask,
            sharp_transition_mask=sharp_transition_mask.values.astype(bool)
            if sharp_transition_mask is not None
            else None,
            water_mask=water_mask.values.astype(bool) if water_mask is not None else None,
            fitting_cost=final_fitting_cost if use_grid_search else None,
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
            "Solver complete: AOT=%.3f, TCWV=%.2f",
            self._finite_mean(aot),
            self._finite_mean(tcwv),
        )

        return result

    def solve_bundle(self, bundle: SolverInputBundle) -> SolverResult:
        """Solve for atmospheric parameters from a SolverInputBundle.

        This is the preferred entry point. It unpacks the bundle fields and
        delegates to :meth:`solve`.
        """

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

    @staticmethod
    def _observation_presence_mask(toa: xr.DataArray) -> BoolArray:
        """Return pixels with at least one finite, physically valid TOA band."""
        values = np.asarray(toa.values, dtype=np.float32)
        if values.ndim == 2:
            values = values[np.newaxis, ...]
        band_valid = np.isfinite(values) & (values > 0.0) & (values < 1.0)
        return cast("BoolArray", np.any(band_valid, axis=0))

    @staticmethod
    def _finite_mean(values: np.ndarray) -> float:
        array = np.asarray(values, dtype=np.float32)
        finite = np.isfinite(array)
        if not np.any(finite):
            return float("nan")
        return float(np.mean(array[finite]))

    @staticmethod
    def _rt_model_supports_jacobian(rt_model: RTModelBackend) -> bool:
        """Return whether backend can provide per-pixel RT Jacobians."""
        return rt_supports_jacobian(rt_model)

    def _fixed_parameter(self) -> FixedAtmosphericParameter:
        value = str(self.config.fixed_atmospheric_parameter).strip().lower()
        if value in {"", "none", "false", "no"}:
            return "none"
        if value in {"aot", "tcwv"}:
            return cast("FixedAtmosphericParameter", value)
        raise ValueError("fixed_atmospheric_parameter must be one of 'none', 'aot', or 'tcwv'")

    def _resample_boolean_mask(
        self,
        mask: np.ndarray,
        target_shape: tuple[int, int],
    ) -> BoolArray:
        resampled = self._resample_mask_to_grid(
            xr.DataArray(np.asarray(mask, dtype=bool), dims=["y", "x"]),
            target_shape,
        )
        return cast("BoolArray", np.asarray(resampled.values, dtype=bool))

    def _build_solver_qa_dataset(
        self,
        *,
        template: xr.DataArray,
        valid_mask: np.ndarray,
        aot: np.ndarray,
        tcwv: np.ndarray,
        invalid_mask: np.ndarray | None,
        zero_obs_mask: np.ndarray | None,
        insufficient_support_mask: np.ndarray | None,
        no_observation_mask: np.ndarray | None,
        sharp_transition_mask: np.ndarray | None,
        water_mask: np.ndarray | None,
        fitting_cost: np.ndarray | None = None,
    ) -> xr.Dataset:
        """Assemble the per-pixel QA dataset (see :mod:`._qa`)."""
        return build_solver_qa_dataset(
            template=template,
            valid_mask=valid_mask,
            aot=aot,
            tcwv=tcwv,
            invalid_mask=invalid_mask,
            zero_obs_mask=zero_obs_mask,
            insufficient_support_mask=insufficient_support_mask,
            no_observation_mask=no_observation_mask,
            sharp_transition_mask=sharp_transition_mask,
            water_mask=water_mask,
            aot_bounds=self.config.aot_bounds,
            tcwv_bounds=self.config.tcwv_bounds,
            fitting_cost=fitting_cost,
        )

    @staticmethod
    def _summarize_solver_qa(qa: xr.Dataset) -> dict[str, float]:
        return {
            "qa_final_invalid_pixels": float(
                np.count_nonzero(np.asarray(qa["invalid_retrieval"].values, dtype=bool))
            ),
            "qa_final_zero_obs_pixels": float(
                np.count_nonzero(np.asarray(qa["zero_obs_support"].values, dtype=bool))
            ),
            "qa_final_insufficient_support_pixels": float(
                np.count_nonzero(
                    np.asarray(qa["insufficient_observation_support"].values, dtype=bool)
                )
            ),
            "qa_final_no_observation_pixels": float(
                np.count_nonzero(np.asarray(qa["no_observation"].values, dtype=bool))
            ),
            "qa_final_sharp_transition_pixels": float(
                np.count_nonzero(np.asarray(qa["sharp_transition_excluded"].values, dtype=bool))
            ),
            "qa_final_water_excluded_pixels": float(
                np.count_nonzero(np.asarray(qa["water_mask_excluded"].values, dtype=bool))
            ),
            "qa_final_aot_lower_boundary_pixels": float(
                np.count_nonzero(np.asarray(qa["aot_lower_boundary"].values, dtype=bool))
            ),
            "qa_final_aot_upper_boundary_pixels": float(
                np.count_nonzero(np.asarray(qa["aot_upper_boundary"].values, dtype=bool))
            ),
            "qa_final_tcwv_lower_boundary_pixels": float(
                np.count_nonzero(np.asarray(qa["tcwv_lower_boundary"].values, dtype=bool))
            ),
            "qa_final_tcwv_upper_boundary_pixels": float(
                np.count_nonzero(np.asarray(qa["tcwv_upper_boundary"].values, dtype=bool))
            ),
            "qa_final_parameter_boundary_pixels": float(
                np.count_nonzero(np.asarray(qa["parameter_boundary"].values, dtype=bool))
            ),
            "qa_final_low_quality_pixels": float(
                np.count_nonzero(np.asarray(qa["low_quality"].values, dtype=bool))
            ),
        }

    @staticmethod
    def _compute_band_weights(bands: list[SensorBand], power: float) -> Float32Array:
        """Compute normalized spectral weights used in observation cost."""
        wavelengths = np.array([b.center_wavelength for b in bands], dtype=np.float32)
        wl_um = np.maximum(wavelengths / 1000.0, 1e-6)
        weights = wl_um**power
        total = float(np.sum(weights))
        if total <= 0:
            return cast(
                "Float32Array",
                np.full(len(bands), 1.0 / max(len(bands), 1), dtype=np.float32),
            )
        return cast("Float32Array", (weights / total).astype(np.float32))

    @staticmethod
    def _adopt_lower_cost_high_aot(
        *,
        near_top: np.ndarray,
        block_size: int,
        shape: tuple[int, int],
        pass1_costs: np.ndarray,
        pass2_costs: np.ndarray,
        pass1: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        pass2: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Adopt the fine high-AOD solution where it strictly lowers the grid cost.

        Compares the per-pixel (or per-block) minimum of each candidate cost cube;
        a pixel switches to the high-AOD pass only if it was already retrieved
        high (``near_top``) *and* the fine axis found a strictly lower cost. This
        is self-gating: low/mid pixels, already resolved on the primary axis,
        keep a lower primary cost and never switch.
        """

        def _cube_min(costs: np.ndarray) -> np.ndarray:
            finite = np.where(np.isfinite(costs), costs, np.float32(np.inf))
            return np.min(finite, axis=(0, 1))

        improve_block = _cube_min(pass2_costs) < (_cube_min(pass1_costs) - np.float32(1e-6))
        if block_size > 1:
            improve = np.repeat(np.repeat(improve_block, block_size, axis=0), block_size, axis=1)
            improve = improve[: shape[0], : shape[1]]
        else:
            improve = improve_block
        # Only ever resolve *upward*: the fine pass exists to fill the sparse top
        # of the log axis (pixels pinned low at the second-highest node), so it
        # must not pull a pixel down into a biased-low minimum the dense lower
        # grid already resolves -- that is the primary pass's responsibility.
        moves_up = pass2[0] > pass1[0]
        switch = near_top & improve & moves_up
        out = []
        for field1, field2 in zip(pass1, pass2, strict=True):
            out.append(np.where(switch, field2, field1).astype(np.float32))
        return out[0], out[1], out[2], out[3]

    def _solve_level_grid_search(
        self,
        *,
        toa: xr.DataArray,
        surface_prior: SurfacePrior,
        geometry: GeometryAngles,
        atmo_prior: AtmosphericState,
        rt_model: RTModelBackend,
        mask: xr.DataArray,
        bands: list[SensorBand],
        cost_config: CostFunctionConfig,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
        """
        Solve one grid level via exhaustive AOT/TCWV candidates + local quadratic fit.

        This path avoids Jacobian-based optimization. Depending on
        ``quadratic_block_size``, it either solves each pixel independently
        or solves one shared AOT/TCWV pair per NxN block and broadcasts the
        block solution back to the full solver grid. RT coefficients are sampled
        on the same block grid. Uncertainty is taken from the fitted local
        Hessian on the active solve grid.

        The observation stacking, candidate coefficient provider, and
        cost-cube evaluation steps live in :mod:`._grid_search`; this method
        orchestrates them and post-processes the refined solution.
        """
        shape = self._get_shape(mask)
        valid_mask = mask.values.astype(bool)

        band_weights = self._compute_band_weights(bands, power=cost_config.band_weight_power)

        # Priors / uncertainties (with floors) for per-pixel prior term.
        aot_prior = atmo_prior.aot.values.astype(np.float32)
        tcwv_prior = atmo_prior.tcwv.values.astype(np.float32)
        aot_prior_unc = np.maximum(
            atmo_prior.aot_unc.values.astype(np.float32), cost_config.min_aot_unc
        )
        tcwv_prior_unc = np.maximum(
            atmo_prior.tcwv_unc.values.astype(np.float32), cost_config.min_tcwv_unc
        )
        fixed_parameter = self._fixed_parameter()
        solve_aot = fixed_parameter != "aot"
        solve_tcwv = fixed_parameter != "tcwv"
        n_aot = max(3, int(self.config.grid_search_aot_points)) if solve_aot else 1
        n_tcwv = max(3, int(self.config.grid_search_tcwv_points)) if solve_tcwv else 1
        aot_axis = (
            _log_aot_axis(self.config.aot_bounds[0], self.config.aot_bounds[1], n_aot)
            if solve_aot
            else fixed_axis_from_prior(aot_prior, self.config.aot_bounds)
        )
        tcwv_axis = (
            np.linspace(
                self.config.tcwv_bounds[0], self.config.tcwv_bounds[1], n_tcwv, dtype=np.float32
            )
            if solve_tcwv
            else fixed_axis_from_prior(tcwv_prior, self.config.tcwv_bounds)
        )

        n_bands = len(bands)
        (
            toa_values,
            no_observation_mask,
            boa_prior,
            boa_unc,
            support_mask,
        ) = prepare_grid_search_observations(
            toa=toa,
            surface_prior=surface_prior,
            n_bands=n_bands,
            shape=shape,
            min_boa_unc=cost_config.min_boa_unc,
            aot_prior=aot_prior,
            tcwv_prior=tcwv_prior,
            aot_prior_unc=aot_prior_unc,
            tcwv_prior_unc=tcwv_prior_unc,
        )
        solve_valid_mask = valid_mask & support_mask

        # RT coefficients are sampled on the block grid; the provider fills one
        # shared coefficient stack per (aot, tcwv) candidate, served from a
        # joint LUT when the backend offers one (see _grid_search for details).
        block_size = max(1, self.config.quadratic_block_size)
        candidate_coeff_provider = build_candidate_coeff_provider(
            rt_model=rt_model,
            bands=bands,
            geometry=geometry,
            atmo_prior=atmo_prior,
            aot_axis=aot_axis,
            tcwv_axis=tcwv_axis,
            solve_aot=solve_aot,
            solve_tcwv=solve_tcwv,
            aot_bounds=self.config.aot_bounds,
            tcwv_bounds=self.config.tcwv_bounds,
            rt_sample_step=block_size,
            shape=shape,
        )

        block_valid_counts = aggregate_valid_counts(solve_valid_mask, block_size)
        block_total_counts = aggregate_block_pixel_counts(shape, block_size)
        min_valid_fraction = min(
            max(float(self.config.quadratic_block_min_valid_fraction), 0.0),
            1.0,
        )
        block_required_counts = np.maximum(
            1,
            np.ceil(block_total_counts.astype(np.float32) * np.float32(min_valid_fraction)).astype(
                np.int32
            ),
        )
        block_support_mask = block_valid_counts >= block_required_counts

        # The Rust kernels are passed from this module's namespace so that
        # monkeypatching the ``multigrid`` attributes keeps intercepting them.
        costs, obs_counts, block_valid_mask, refine_valid_mask = evaluate_candidate_cost_cube(
            coeff_provider=candidate_coeff_provider,
            aot_axis=aot_axis,
            tcwv_axis=tcwv_axis,
            toa_values=toa_values,
            boa_prior=boa_prior,
            boa_unc=boa_unc,
            band_weights=band_weights,
            solve_valid_mask=solve_valid_mask,
            aot_prior=aot_prior,
            tcwv_prior=tcwv_prior,
            aot_prior_unc=aot_prior_unc,
            tcwv_prior_unc=tcwv_prior_unc,
            block_size=block_size,
            fixed_parameter=fixed_parameter,
            block_support_mask=block_support_mask,
            block_cost_cube_fn=evaluate_block_grid_search_cost_cube_with_provider_qa,
            pixel_cost_cube_fn=evaluate_grid_search_cost_cube_with_provider_qa,
        )

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
            refine_valid_mask,
            fixed_parameter,
        )
        aot_best = np.asarray(aot_best, dtype=np.float32)
        tcwv_best = np.asarray(tcwv_best, dtype=np.float32)
        aot_unc = np.asarray(aot_unc, dtype=np.float32)
        tcwv_unc = np.asarray(tcwv_unc, dtype=np.float32)

        # Second-pass high-AOD refine (see _HIGH_AOT_REFINE_* above): where the
        # primary solve pins a pixel high, re-search a fine axis over the sparse
        # top of the log grid and keep it only where the grid cost strictly drops.
        if solve_aot and bool(np.any(aot_best >= np.float32(_HIGH_AOT_REFINE_MIN))):
            near_top = solve_valid_mask & (aot_best >= np.float32(_HIGH_AOT_REFINE_MIN))
            fine_axis = np.geomspace(
                _HIGH_AOT_REFINE_LO,
                float(self.config.aot_bounds[1]),
                _HIGH_AOT_REFINE_POINTS,
                dtype=np.float32,
            )
            fine_provider = build_candidate_coeff_provider(
                rt_model=rt_model,
                bands=bands,
                geometry=geometry,
                atmo_prior=atmo_prior,
                aot_axis=fine_axis,
                tcwv_axis=tcwv_axis,
                solve_aot=solve_aot,
                solve_tcwv=solve_tcwv,
                aot_bounds=self.config.aot_bounds,
                tcwv_bounds=self.config.tcwv_bounds,
                rt_sample_step=block_size,
                shape=shape,
            )
            costs_high, obs_counts_high, _block_valid_high, refine_valid_high = (
                evaluate_candidate_cost_cube(
                    coeff_provider=fine_provider,
                    aot_axis=fine_axis,
                    tcwv_axis=tcwv_axis,
                    toa_values=toa_values,
                    boa_prior=boa_prior,
                    boa_unc=boa_unc,
                    band_weights=band_weights,
                    solve_valid_mask=solve_valid_mask,
                    aot_prior=aot_prior,
                    tcwv_prior=tcwv_prior,
                    aot_prior_unc=aot_prior_unc,
                    tcwv_prior_unc=tcwv_prior_unc,
                    block_size=block_size,
                    fixed_parameter=fixed_parameter,
                    block_support_mask=block_support_mask,
                    block_cost_cube_fn=evaluate_block_grid_search_cost_cube_with_provider_qa,
                    pixel_cost_cube_fn=evaluate_grid_search_cost_cube_with_provider_qa,
                )
            )
            refined_high = quadratic_refine_grid_search_qa(
                costs_high.astype(np.float32, copy=False),
                obs_counts_high.astype(np.uint16, copy=False),
                fine_axis.astype(np.float32, copy=False),
                tcwv_axis.astype(np.float32, copy=False),
                refine_valid_high,
                fixed_parameter,
            )
            aot_best, tcwv_best, aot_unc, tcwv_unc = self._adopt_lower_cost_high_aot(
                near_top=near_top,
                block_size=block_size,
                shape=shape,
                pass1_costs=costs,
                pass2_costs=np.asarray(costs_high, dtype=np.float32),
                pass1=(aot_best, tcwv_best, aot_unc, tcwv_unc),
                pass2=(
                    np.asarray(refined_high[0], dtype=np.float32),
                    np.asarray(refined_high[1], dtype=np.float32),
                    np.asarray(refined_high[2], dtype=np.float32),
                    np.asarray(refined_high[3], dtype=np.float32),
                ),
            )

        invalid_mask = np.asarray(invalid_mask, dtype=bool)
        boundary_mask = np.asarray(boundary_mask, dtype=bool)
        lower_aot_boundary_mask = np.asarray(lower_aot_boundary_mask, dtype=bool)
        zero_obs_mask = np.asarray(zero_obs_mask, dtype=bool)
        invalid_block_mask = invalid_mask.copy()
        insufficient_support_mask = ~support_mask

        if block_size > 1:
            block_solver_invalid_mask = ~block_valid_mask
            block_insufficient_support_mask = ~block_support_mask
            invalid_mask = invalid_mask | block_solver_invalid_mask
            zero_obs_mask = zero_obs_mask | block_insufficient_support_mask
            trusted_block_mask = (
                block_valid_mask
                & ~invalid_mask
                & ~zero_obs_mask
                & ~boundary_mask
                & ~lower_aot_boundary_mask
                & np.isfinite(aot_best)
                & np.isfinite(tcwv_best)
            )
            if solve_aot:
                aot_best = self._smooth_grid_search_field(
                    aot_best,
                    gamma=cost_config.aot_gamma,
                    delta=cost_config.smoothness_delta,
                    n_iter=40,
                    trusted_mask=trusted_block_mask,
                ).astype(np.float32)
            if solve_tcwv:
                tcwv_best = self._smooth_grid_search_field(
                    tcwv_best,
                    gamma=cost_config.tcwv_gamma,
                    delta=cost_config.smoothness_delta,
                    n_iter=40,
                    trusted_mask=trusted_block_mask,
                ).astype(np.float32)

            invalid_mask_full = self._broadcast_to_full(invalid_mask, shape, block_size) > 0.5
            zero_obs_mask_full = self._broadcast_to_full(zero_obs_mask, shape, block_size) > 0.5
            boundary_mask_full = self._broadcast_to_full(boundary_mask, shape, block_size) > 0.5
            lower_aot_boundary_mask_full = (
                self._broadcast_to_full(
                    lower_aot_boundary_mask,
                    shape,
                    block_size,
                )
                > 0.5
            )

            aot_best = self._broadcast_to_full(aot_best, shape, block_size).astype(np.float32)
            tcwv_best = self._broadcast_to_full(tcwv_best, shape, block_size).astype(np.float32)
            aot_unc = self._broadcast_to_full(aot_unc, shape, block_size).astype(np.float32)
            tcwv_unc = self._broadcast_to_full(tcwv_unc, shape, block_size).astype(np.float32)
            invalid_mask = invalid_mask_full
            zero_obs_mask = zero_obs_mask_full
            boundary_mask = boundary_mask_full
            lower_aot_boundary_mask = lower_aot_boundary_mask_full
            insufficient_support_mask = (
                self._broadcast_to_full(
                    block_insufficient_support_mask,
                    shape,
                    block_size,
                )
                > 0.5
            )
        else:
            zero_obs_mask = zero_obs_mask | insufficient_support_mask
            invalid_mask = invalid_mask | insufficient_support_mask

        invalid_floor_pixels = int(
            np.count_nonzero(
                solve_aot
                & invalid_mask
                & solve_valid_mask
                & np.isclose(aot_best, float(aot_axis[0]), rtol=0.0, atol=1.0e-6)
            )
        )
        prior_floor_pixels = int(
            np.count_nonzero(
                solve_aot
                & invalid_mask
                & solve_valid_mask
                & np.isclose(aot_prior, float(aot_axis[0]), rtol=0.0, atol=1.0e-6)
            )
        )

        if block_size <= 1:
            trusted_pixel_mask = (
                solve_valid_mask
                & ~invalid_mask
                & ~zero_obs_mask
                & ~boundary_mask
                & ~lower_aot_boundary_mask
                & np.isfinite(aot_best)
                & np.isfinite(tcwv_best)
            )
            # Apply edge-preserving spatial smoothing using only trusted
            # retrievals as seeds. Invalid, zero-support, and boundary-hit
            # retrievals are filled from neighbouring trusted seeds instead of
            # being allowed to influence the smoothed field.
            if solve_aot:
                aot_best = self._smooth_grid_search_field(
                    aot_best,
                    gamma=cost_config.aot_gamma,
                    delta=cost_config.smoothness_delta,
                    n_iter=40,
                    trusted_mask=trusted_pixel_mask,
                ).astype(np.float32)
            if solve_tcwv:
                tcwv_best = self._smooth_grid_search_field(
                    tcwv_best,
                    gamma=cost_config.tcwv_gamma,
                    delta=cost_config.smoothness_delta,
                    n_iter=40,
                    trusted_mask=trusted_pixel_mask,
                ).astype(np.float32)

        # Handle QA-invalid pixels.  If there are enough valid neighbours the
        # smoothing has already propagated sensible values into them.  When ALL
        # pixels are invalid (complete grid-search failure) the smoothed field
        # is meaningless, so fall back to the prior.
        if np.any(invalid_mask):
            all_invalid = np.all(invalid_mask | ~solve_valid_mask)
            if all_invalid:
                aot_best = np.where(invalid_mask, aot_prior, aot_best).astype(
                    np.float32, copy=False
                )
                tcwv_best = np.where(invalid_mask, tcwv_prior, tcwv_best).astype(
                    np.float32, copy=False
                )
            # Inflate uncertainty at invalid pixels — their values come from
            # spatial interpolation (or prior fallback), not direct retrieval.
            aot_unc = np.where(
                invalid_mask,
                np.maximum(aot_prior_unc, np.float32(0.1)),
                aot_unc,
            ).astype(np.float32, copy=False)
            tcwv_unc = np.where(
                invalid_mask,
                np.maximum(tcwv_prior_unc, np.float32(0.5)),
                tcwv_unc,
            ).astype(np.float32, copy=False)

        if fixed_parameter == "aot":
            aot_best = aot_prior.astype(np.float32, copy=True)
            aot_unc = aot_prior_unc.astype(np.float32, copy=True)
        elif fixed_parameter == "tcwv":
            tcwv_best = tcwv_prior.astype(np.float32, copy=True)
            tcwv_unc = tcwv_prior_unc.astype(np.float32, copy=True)

        if np.any(no_observation_mask):
            invalid_mask = invalid_mask | no_observation_mask
            zero_obs_mask = zero_obs_mask | no_observation_mask
            aot_best = np.where(no_observation_mask, np.nan, aot_best).astype(
                np.float32, copy=False
            )
            tcwv_best = np.where(no_observation_mask, np.nan, tcwv_best).astype(
                np.float32, copy=False
            )
            aot_unc = np.where(no_observation_mask, np.nan, aot_unc).astype(np.float32, copy=False)
            tcwv_unc = np.where(no_observation_mask, np.nan, tcwv_unc).astype(
                np.float32, copy=False
            )

        flat = costs.reshape(n_aot * n_tcwv, -1)
        safe_flat = np.where(np.isfinite(flat), flat, np.inf)
        best_flat = np.argmin(safe_flat, axis=0)
        selected_costs = safe_flat[best_flat, np.arange(safe_flat.shape[1])]
        if block_size > 1:
            supported = block_valid_mask & ~invalid_block_mask
            selected_costs = selected_costs.reshape(costs.shape[-2:])
            mean_cost = (
                float(np.sum(selected_costs[supported]))
                / float(np.sum(block_valid_counts[supported]))
                if np.any(supported)
                else float("nan")
            )
        else:
            supported = (solve_valid_mask & ~invalid_mask).reshape(-1)
            mean_cost = (
                float(np.mean(selected_costs[supported])) if np.any(supported) else float("nan")
            )
        valid_pixels = int(np.count_nonzero(solve_valid_mask))
        solve_invalid_pixels = int(np.count_nonzero(invalid_mask & solve_valid_mask))
        invalid_pixels = int(np.count_nonzero(invalid_mask))
        zero_obs_pixels = int(np.count_nonzero(zero_obs_mask))
        no_observation_pixels = int(np.count_nonzero(no_observation_mask))
        boundary_pixels = int(np.count_nonzero(boundary_mask & solve_valid_mask))
        lower_aot_boundary_pixels = int(
            np.count_nonzero(lower_aot_boundary_mask & solve_valid_mask)
        )
        if invalid_pixels:
            logger.warning(
                "Grid-search QA flagged %d pixels as invalid with %d solve-supported pixels; %d had zero observation support and %d collapsed to the AOT floor.",
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
                "qa_solve_invalid_pixels": float(solve_invalid_pixels),
                "qa_zero_obs_pixels": float(zero_obs_pixels),
                "qa_no_observation_pixels": float(no_observation_pixels),
                "qa_boundary_pixels": float(boundary_pixels),
                "qa_lower_aot_boundary_pixels": float(lower_aot_boundary_pixels),
                "qa_invalid_floor_pixels": float(invalid_floor_pixels),
                "qa_prior_floor_pixels": float(prior_floor_pixels),
                "qa_insufficient_support_pixels": float(
                    np.count_nonzero(insufficient_support_mask)
                ),
                "qa_invalid_mask": invalid_mask.astype(bool, copy=False),
                "qa_zero_obs_mask": zero_obs_mask.astype(bool, copy=False),
                "qa_insufficient_support_mask": insufficient_support_mask.astype(bool, copy=False),
                "qa_fitting_cost_map": (
                    selected_costs.reshape(costs.shape[-2:])
                    if block_size > 1
                    else selected_costs.reshape(shape)
                ).astype(np.float32, copy=False),
            },
        )

    def _get_shape(self, arr: xr.DataArray) -> tuple[int, int]:
        """Get (ny, nx) shape from DataArray."""
        if "y" in arr.dims and "x" in arr.dims:
            return (int(arr.sizes["y"]), int(arr.sizes["x"]))
        return (int(arr.shape[0]), int(arr.shape[1]))

    def _create_mask(
        self,
        cloud_mask: xr.DataArray,
        toa: xr.DataArray,
        surface_prior: SurfacePrior,
        *,
        sharp_transition_mask: xr.DataArray | None = None,
        water_mask: xr.DataArray | None = None,
    ) -> xr.DataArray:
        """Create combined valid pixel mask."""
        return build_solver_valid_mask(
            cloud_mask,
            toa,
            surface_prior,
            sharp_transition_mask=sharp_transition_mask,
            water_mask=water_mask,
        )

    def _compute_grid_levels(self, full_shape: tuple[int, int]) -> list[tuple[int, int]]:
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
            # Exponential spacing from min to full, capped at input size
            scale = 2 ** (max_levels - 1 - i)
            shape_y = min(ny, max(min_size, ny // scale))
            shape_x = min(nx, max(min_size, nx // scale))
            shapes.append((shape_y, shape_x))

        # Deduplicate while preserving order (small inputs may collapse levels)
        seen: set[tuple[int, int]] = set()
        unique: list[tuple[int, int]] = []
        for s in shapes:
            if s not in seen:
                seen.add(s)
                unique.append(s)

        return unique

    # -----------------------------------------------------------------
    # Block-grid helpers
    # -----------------------------------------------------------------

    @staticmethod
    def _subsample_da(da: xr.DataArray, step: int) -> xr.DataArray:
        """Subsample a 2-D DataArray by picking every *step*-th pixel."""
        return da[::step, ::step]

    @staticmethod
    def _broadcast_to_full(
        coarse: np.ndarray,
        full_shape: tuple[int, int],
        step: int,
    ) -> np.ndarray:
        """Broadcast a coarse-grid 2-D array back to *full_shape* via nearest-neighbour repeat."""
        return np.repeat(
            np.repeat(coarse, step, axis=0),
            step,
            axis=1,
        )[: full_shape[0], : full_shape[1]]

    @staticmethod
    def _smooth_grid_search_field(
        values: np.ndarray,
        *,
        gamma: float,
        delta: float,
        n_iter: int,
        trusted_mask: np.ndarray,
    ) -> Float32Array:
        """Smooth a retrieved field using only trusted retrievals as seeds."""
        from siac.algorithms.solver.aod_smoothing import nearest_seed_fill
        from siac.algorithms.solver.cost import apply_smoothness_filter

        source = np.asarray(values, dtype=np.float32)
        trusted = np.asarray(trusted_mask, dtype=bool) & np.isfinite(source)
        if source.shape != trusted.shape:
            raise ValueError(
                f"Value/trusted mask shape mismatch: {source.shape} vs {trusted.shape}"
            )
        if not np.any(trusted):
            return cast("Float32Array", source.astype(np.float32, copy=True))

        filled = nearest_seed_fill(source, trusted)
        if float(gamma) > 0.0:
            smoothed = apply_smoothness_filter(
                filled,
                gamma=float(gamma),
                delta=float(delta),
                n_iter=int(n_iter),
            )
        else:
            smoothed = filled
        return cast("Float32Array", np.asarray(smoothed, dtype=np.float32))

    def _resample_field(self, field: np.ndarray, target_shape: tuple[int, int]) -> Float64Array:
        """Resample 2D field to target shape.

        REVIEW.md §1.2 #9: previously every call promoted ``field`` to
        ``float64`` even when the shape was unchanged, producing a
        ``float32 → float64 → float32`` round-trip at every multigrid
        level when ``aot``/``tcwv`` are stored as float32. The same-shape
        fast path now returns the input view directly when it's already
        float64, and only allocates a copy when a dtype promotion is
        actually needed for the Rust kernels (which take f64 — promoting
        to f64 there is unavoidable).
        """
        if field.shape == target_shape:
            if field.dtype == np.float64:
                return cast("Float64Array", field)
            return cast("Float64Array", np.asarray(field, dtype=np.float64))

        data = np.ascontiguousarray(field, dtype=np.float64)
        if target_shape[0] < field.shape[0]:
            return cast(
                "Float64Array",
                np.asarray(remap_to_coarse_grid(data, target_shape[0], target_shape[1])),
            )
        return cast(
            "Float64Array",
            np.asarray(interpolate_to_fine_grid(data, target_shape[0], target_shape[1])),
        )

    def _resample_to_grid(self, data: xr.DataArray, shape: tuple[int, int]) -> xr.DataArray:
        """Resample DataArray to target grid."""
        if data.ndim == 3:
            # Multi-band
            result = np.stack(
                [self._resample_field(data.values[i], shape) for i in range(data.shape[0])]
            )
            return xr.DataArray(result, dims=["band", "y", "x"])
        else:
            result = self._resample_field(data.values, shape)
            return xr.DataArray(result, dims=["y", "x"])

    def _resample_mask_to_grid(self, mask: xr.DataArray, shape: tuple[int, int]) -> xr.DataArray:
        """Resample mask to target grid using max pooling."""
        if mask.shape == shape:
            return mask

        if shape[0] <= mask.shape[0] and shape[1] <= mask.shape[1]:
            # Match the center-based coarse-grid assignment used for numeric
            # field remapping so coarse validity covers the same source pixels.
            result: BoolArray = np.ones(shape, dtype=bool)
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
        n = aot_init.size
        bounds_aot = [self.config.aot_bounds] * n
        bounds_tcwv = [self.config.tcwv_bounds] * n
        fixed_parameter = self._fixed_parameter()
        fixed_aot = np.asarray(cost_func.aot_prior, dtype=np.float64).reshape(aot_init.shape)
        fixed_tcwv = np.asarray(cost_func.tcwv_prior, dtype=np.float64).reshape(tcwv_init.shape)

        if fixed_parameter == "aot":
            p0 = tcwv_init.ravel()
            bounds = bounds_tcwv

            def objective(p: np.ndarray) -> tuple[float, np.ndarray]:
                full_p = np.concatenate([fixed_aot.ravel(), p])
                cost, grad = cost_func(full_p)
                return cost, grad[n:]

        elif fixed_parameter == "tcwv":
            p0 = aot_init.ravel()
            bounds = bounds_aot

            def objective(p: np.ndarray) -> tuple[float, np.ndarray]:
                full_p = np.concatenate([p, fixed_tcwv.ravel()])
                cost, grad = cost_func(full_p)
                return cost, grad[:n]

        else:
            p0 = np.concatenate([aot_init.ravel(), tcwv_init.ravel()])
            bounds = bounds_aot + bounds_tcwv
            objective = cost_func

        # Adjust parameters based on level
        maxiter = self.config.max_iter_per_level
        if level < 2:
            maxiter = min(maxiter, 100)  # Fewer iterations at coarse levels

        # Run optimization
        result = optimize.minimize(
            objective,
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
        if fixed_parameter == "aot":
            aot_solved = fixed_aot.astype(np.float32, copy=False)
            tcwv_solved = result.x.reshape(tcwv_init.shape)
        elif fixed_parameter == "tcwv":
            aot_solved = result.x.reshape(aot_init.shape)
            tcwv_solved = fixed_tcwv.astype(np.float32, copy=False)
        else:
            aot_solved = result.x[:n].reshape(aot_init.shape)
            tcwv_solved = result.x[n:].reshape(tcwv_init.shape)

        return aot_solved, tcwv_solved, result

    @staticmethod
    def _local_huber_weight(field: np.ndarray, delta: float) -> np.ndarray:
        """Compute mean local Pseudo-Huber weight w = 1/√(1+(∇f/δ)²).

        Returns a per-pixel value in [0, 1].  Near 1 in smooth regions (full
        smoothing applied, lower uncertainty), near 0 at sharp edges (less
        smoothing, higher uncertainty).
        """
        delta2 = delta * delta
        dy = np.diff(field, axis=0)
        dx = np.diff(field, axis=1)
        wy = 1.0 / np.sqrt(1.0 + (dy * dy) / delta2)
        wx = 1.0 / np.sqrt(1.0 + (dx * dx) / delta2)
        # Average the four neighbour weights at each interior pixel.
        w = np.ones_like(field)
        count = np.ones_like(field)
        w[:-1, :] += wy
        count[:-1, :] += 1
        w[1:, :] += wy
        count[1:, :] += 1
        w[:, :-1] += wx
        count[:, :-1] += 1
        w[:, 1:] += wx
        count[:, 1:] += 1
        return w / count

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
        scaled by the effective local smoothing weight.  Smooth regions
        (low gradients) have higher effective smoothing and therefore lower
        uncertainty; sharp features (hotspots) have less smoothing and
        higher uncertainty.
        """
        delta = cost_func.config.smoothness_delta
        gamma_aot = cost_func.config.aot_gamma
        gamma_tcwv = cost_func.config.tcwv_gamma

        # Local Pseudo-Huber weight: 1 in smooth regions, →0 at edges.
        w_aot = self._local_huber_weight(aot, delta)
        w_tcwv = self._local_huber_weight(tcwv, delta)

        # Effective filter strength per pixel: 1/(1 + γ² * w).
        # Where w≈1 (smooth), filter is strong → low unc adjustment.
        # Where w≈0 (edge), filter is weak → high unc adjustment.
        eff_aot = 1.0 / (1.0 + gamma_aot**2 * w_aot)
        eff_tcwv = 1.0 / (1.0 + gamma_tcwv**2 * w_tcwv)

        adj_aot = 1.0 / np.sqrt(np.clip(eff_aot, 1e-6, None))
        adj_tcwv = 1.0 / np.sqrt(np.clip(eff_tcwv, 1e-6, None))

        # Normalize so that the mean adjustment in smooth regions is ~1.
        adj_aot = adj_aot / np.maximum(adj_aot.mean(), 1e-6)
        adj_tcwv = adj_tcwv / np.maximum(adj_tcwv.mean(), 1e-6)

        # Posterior uncertainty
        aot_unc = np.maximum(aot * 0.5, 0.02) * adj_aot
        tcwv_unc = np.maximum(tcwv * 0.1, 0.1) * adj_tcwv

        # Ensure shape matches
        aot_unc = self._resample_field(aot_unc, aot.shape)
        tcwv_unc = self._resample_field(tcwv_unc, tcwv.shape)

        return aot_unc, tcwv_unc


class StagedMultiGridSolver:
    """Chain one or more atmospheric solver stages.

    The staged wrapper is deliberately conservative: it can chain today's
    supported AOT/TCWV combinations and carries TCO3 through the atmospheric
    state, but it fails early for requested TCO3 solving until the RT Jacobian
    and Rust/grid-search contracts expose ozone as a solved parameter.
    """

    _VALID_PARAMETERS = {"aot", "tcwv", "tco3"}

    def __init__(self, config: MultiGridConfig | None = None):
        self.config = config or MultiGridConfig()

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
        stages = self._normalized_stages()
        if not stages:
            return MultiGridSolver(self.config).solve(
                toa=toa,
                surface_prior=surface_prior,
                geometry=geometry,
                atmo_prior=atmo_prior,
                rt_model=rt_model,
                cloud_mask=cloud_mask,
                bands=bands,
                sharp_transition_mask=sharp_transition_mask,
                water_mask=water_mask,
            )

        current_state = atmo_prior
        total_iterations = 0
        final_cost = float("nan")
        final_success = True
        messages: list[str] = []
        level_history: list[dict[str, Any]] = []
        final_qa: xr.Dataset | None = None
        final_result: SolverResult | None = None

        for index, stage in enumerate(stages):
            stage_solve = set(stage.solve)
            if not stage_solve:
                logger.info("Skipping solver stage %s: no solved parameters.", stage.name)
                continue

            fixed_parameter = self._fixed_parameter_for_stage(stage)
            stage_prior = atmo_prior if stage.initial_state == "prior" else current_state
            stage_toa, stage_surface_prior, stage_bands = self._select_stage_inputs(
                toa=toa,
                surface_prior=surface_prior,
                bands=bands,
                stage=stage,
            )
            stage_config = replace(
                self.config,
                fixed_atmospheric_parameter=fixed_parameter,
                stages=(),
            )

            logger.info(
                "Solver stage %s: solving %s with fixed %s using %d band(s).",
                stage.name,
                ", ".join(stage.solve),
                ", ".join(sorted(self._fixed_parameters_for_log(stage))),
                len(stage_bands),
            )
            result = MultiGridSolver(stage_config).solve(
                toa=stage_toa,
                surface_prior=stage_surface_prior,
                geometry=geometry,
                atmo_prior=stage_prior,
                rt_model=rt_model,
                cloud_mask=cloud_mask,
                bands=stage_bands,
                sharp_transition_mask=sharp_transition_mask,
                water_mask=water_mask,
            )

            updates: dict[str, xr.DataArray] = {}
            uncertainties: dict[str, xr.DataArray] = {}
            if "aot" in stage_solve:
                updates["aot"] = result.aot
                uncertainties["aot"] = result.aot_unc
            if "tcwv" in stage_solve:
                updates["tcwv"] = result.tcwv
                uncertainties["tcwv"] = result.tcwv_unc
            current_state = stage_prior.with_updated_parameters(updates, uncertainties)

            total_iterations += int(result.n_iterations)
            final_cost = float(result.final_cost)
            final_success = final_success and bool(result.success)
            messages.append(f"{stage.name}: {result.message}")
            final_qa = result.qa
            final_result = result
            for entry in result.level_history:
                tagged_entry = dict(entry)
                tagged_entry["stage"] = stage.name
                tagged_entry["stage_index"] = index
                tagged_entry["stage_solve"] = tuple(stage.solve)
                level_history.append(tagged_entry)

        if final_result is None:
            return SolverResult(
                aot=current_state.aot,
                tcwv=current_state.tcwv,
                aot_unc=current_state.aot_unc,
                tcwv_unc=current_state.tcwv_unc,
                n_iterations=0,
                final_cost=float("nan"),
                success=True,
                message="No solver stages configured with active parameters",
                qa=None,
                level_history=[],
                atmo_state=current_state,
            )

        return SolverResult(
            aot=current_state.aot,
            tcwv=current_state.tcwv,
            aot_unc=current_state.aot_unc,
            tcwv_unc=current_state.tcwv_unc,
            n_iterations=total_iterations,
            final_cost=final_cost,
            success=final_success,
            message="; ".join(messages),
            qa=final_qa,
            level_history=level_history,
            atmo_state=current_state,
        )

    def solve_bundle(self, bundle: SolverInputBundle) -> SolverResult:
        """Solve for atmospheric parameters from a SolverInputBundle.

        This is the preferred entry point. It unpacks the bundle fields and
        delegates to :meth:`solve`.
        """

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

    def _normalized_stages(self) -> tuple[SolverStageConfig, ...]:
        raw_stages = tuple(getattr(self.config, "stages", ()) or ())
        return tuple(self._coerce_stage(stage, index) for index, stage in enumerate(raw_stages))

    @classmethod
    def _coerce_stage(cls, stage: Any, index: int) -> SolverStageConfig:
        if isinstance(stage, SolverStageConfig):
            candidate = stage
        elif isinstance(stage, dict):
            candidate = SolverStageConfig(
                name=str(stage.get("name") or f"stage_{index + 1}"),
                solve=cls._parameter_tuple(stage.get("solve"), ("aot", "tcwv")),
                fixed=cls._parameter_tuple(stage.get("fixed"), ("tco3",)),
                bands=cls._band_tuple(stage.get("bands")),
                initial_state=cast(
                    "StageInitialState",
                    str(stage.get("initial_state") or "previous").lower(),
                ),
            )
        else:
            candidate = SolverStageConfig(
                name=str(getattr(stage, "name", None) or f"stage_{index + 1}"),
                solve=cls._parameter_tuple(getattr(stage, "solve", None), ("aot", "tcwv")),
                fixed=cls._parameter_tuple(getattr(stage, "fixed", None), ("tco3",)),
                bands=cls._band_tuple(getattr(stage, "bands", None)),
                initial_state=cast(
                    "StageInitialState",
                    str(getattr(stage, "initial_state", "previous")).lower(),
                ),
            )

        cls._validate_stage(candidate)
        return candidate

    @classmethod
    def _parameter_tuple(
        cls,
        value: Any,
        default: tuple[AtmosphericParameterName, ...],
    ) -> tuple[AtmosphericParameterName, ...]:
        if value is None:
            raw = default
        elif isinstance(value, str):
            raw = (value,)
        else:
            raw = tuple(value)
        normalized = tuple(str(item).strip().lower() for item in raw if str(item).strip())
        return cast("tuple[AtmosphericParameterName, ...]", normalized)

    @staticmethod
    def _band_tuple(value: Any) -> tuple[str, ...] | None:
        if value is None:
            return None
        if isinstance(value, str):
            return (value,)
        return tuple(str(item) for item in value)

    @classmethod
    def _validate_stage(cls, stage: SolverStageConfig) -> None:
        solve = set(stage.solve)
        fixed = set(stage.fixed)
        unknown = (solve | fixed) - cls._VALID_PARAMETERS
        if unknown:
            raise ValueError(
                f"Solver stage {stage.name!r} references unknown parameter(s): "
                f"{', '.join(sorted(unknown))}"
            )
        overlap = solve & fixed
        if overlap:
            raise ValueError(
                f"Solver stage {stage.name!r} cannot both solve and fix "
                f"{', '.join(sorted(overlap))}"
            )
        if stage.initial_state not in {"prior", "previous"}:
            raise ValueError(
                f"Solver stage {stage.name!r} initial_state must be 'prior' or 'previous'"
            )
        if "tco3" in solve:
            raise NotImplementedError(
                "Solver stage "
                f"{stage.name!r} requests TCO3 solving, but ozone is currently "
                "only a carried atmospheric state field. Add RT TCO3 Jacobian/grid-search "
                "support before enabling solve=['tco3']."
            )

    @classmethod
    def _fixed_parameter_for_stage(
        cls,
        stage: SolverStageConfig,
    ) -> FixedAtmosphericParameter:
        cls._validate_stage(stage)
        solve = set(stage.solve)
        solves_aot = "aot" in solve
        solves_tcwv = "tcwv" in solve
        if solves_aot and solves_tcwv:
            return "none"
        if solves_aot:
            return "tcwv"
        if solves_tcwv:
            return "aot"
        raise ValueError(
            f"Solver stage {stage.name!r} must solve at least one supported parameter: aot or tcwv"
        )

    @classmethod
    def _fixed_parameters_for_log(cls, stage: SolverStageConfig) -> set[str]:
        return cls._VALID_PARAMETERS - set(stage.solve)

    @classmethod
    def _select_stage_inputs(
        cls,
        *,
        toa: xr.DataArray,
        surface_prior: SurfacePrior,
        bands: list[SensorBand],
        stage: SolverStageConfig,
    ) -> tuple[xr.DataArray, SurfacePrior, list[SensorBand]]:
        if not stage.bands:
            return toa, surface_prior, bands

        requested = tuple(stage.bands)
        index_by_name = {band.name: index for index, band in enumerate(bands)}
        missing = [name for name in requested if name not in index_by_name]
        if missing:
            raise ValueError(
                f"Solver stage {stage.name!r} requested unknown band(s): {', '.join(missing)}"
            )

        indices = [index_by_name[name] for name in requested]
        selected_bands = [bands[index] for index in indices]
        selected_toa = cls._select_band_axis(toa, indices, len(bands))
        selected_surface_prior = cls._select_surface_prior_bands(
            surface_prior,
            indices,
            len(bands),
        )
        return selected_toa, selected_surface_prior, selected_bands

    @staticmethod
    def _select_band_axis(
        data: xr.DataArray,
        indices: list[int],
        source_band_count: int,
    ) -> xr.DataArray:
        if data.ndim >= 3 and data.shape[0] == source_band_count:
            return data.isel({data.dims[0]: indices})
        return data

    @classmethod
    def _select_surface_prior_bands(
        cls,
        surface_prior: SurfacePrior,
        indices: list[int],
        source_band_count: int,
    ) -> SurfacePrior:
        kernels = surface_prior.kernels
        if kernels is not None:
            kernels = BRDFKernelWeights(
                f0=cls._select_band_axis(kernels.f0, indices, source_band_count),
                f1=cls._select_band_axis(kernels.f1, indices, source_band_count),
                f2=cls._select_band_axis(kernels.f2, indices, source_band_count),
                f0_unc=cls._select_band_axis(kernels.f0_unc, indices, source_band_count),
                f1_unc=cls._select_band_axis(kernels.f1_unc, indices, source_band_count),
                f2_unc=cls._select_band_axis(kernels.f2_unc, indices, source_band_count),
                reflectance_unc=cls._select_band_axis(
                    kernels.reflectance_unc,
                    indices,
                    source_band_count,
                )
                if kernels.reflectance_unc is not None
                else None,
            )

        return SurfacePrior(
            boa=cls._select_band_axis(surface_prior.boa, indices, source_band_count),
            boa_unc=cls._select_band_axis(surface_prior.boa_unc, indices, source_band_count),
            kernels=kernels,
            mask=cls._select_band_axis(surface_prior.mask, indices, source_band_count),
            monthly_composites=surface_prior.monthly_composites,
        )


def solve_atmospheric_parameters(
    toa: xr.DataArray,
    surface_prior: SurfacePrior,
    geometry: GeometryAngles,
    atmo_prior: AtmosphericState,
    rt_model: RTModelBackend,
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
    solver_cls = (
        StagedMultiGridSolver
        if config is not None and getattr(config, "stages", ())
        else MultiGridSolver
    )
    solver = solver_cls(config)
    return solver.solve(
        toa=toa,
        surface_prior=surface_prior,
        geometry=geometry,
        atmo_prior=atmo_prior,
        rt_model=rt_model,
        cloud_mask=cloud_mask,
        bands=bands,
    )
