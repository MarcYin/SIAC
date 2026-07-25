"""Solver assembly for the atmospheric inversion stage.

The concrete solver is chosen by ``config.algorithms.solver.method`` and looked
up in :data:`SOLVER_METHOD_REGISTRY`. ``multigrid`` (default) is the Bayesian
multi-grid inversion; ``surface_driven`` sweeps a fixed AOT axis at the prior
TCWV and picks the per-pixel arg-min of the surface-prior mismatch.
"""

from __future__ import annotations

import inspect
from typing import TYPE_CHECKING, Any

from siac.algorithms.solver import MultiGridConfig, MultiGridSolver, StagedMultiGridSolver
from siac.app.registry import SOLVER_METHOD_REGISTRY
from siac.runtime import SolvedAtmosphere, SolverInputBundle

if TYPE_CHECKING:
    from siac.workflows.pipeline import SolverFn


def _wrap_solved(result: Any, inputs: SolverInputBundle) -> SolvedAtmosphere:
    """Wrap an internal ``SolverResult`` as the M5 ``SolvedAtmosphere`` payload."""
    solved_atmo = getattr(result, "atmo_state", None)
    if solved_atmo is None:
        solved_atmo = inputs.atmo_prior.with_updated_aot_tcwv(
            aot=result.aot,
            tcwv=result.tcwv,
            aot_unc=result.aot_unc,
            tcwv_unc=result.tcwv_unc,
        )
    return SolvedAtmosphere(
        atmo_state=solved_atmo,
        aot=solved_atmo.aot,
        tcwv=solved_atmo.tcwv,
        aot_unc=solved_atmo.aot_unc,
        tcwv_unc=solved_atmo.tcwv_unc,
        cost_final=float(result.final_cost),
        n_iterations=result.n_iterations,
        converged=result.success,
        qa=getattr(result, "qa", None),
        level_history=tuple(getattr(result, "level_history", ())),
        diagnostics=dict(getattr(result, "diagnostics", {}) or {}),
    )


@SOLVER_METHOD_REGISTRY.register("multigrid")
def _build_multigrid_solver(config: Any) -> SolverFn:
    def _solver(inputs: SolverInputBundle, _config: Any) -> SolvedAtmosphere:
        sc = config.algorithms.solver
        solver_config = MultiGridConfig(
            aot_gamma=sc.aot_gamma,
            tcwv_gamma=sc.tcwv_gamma,
            aot_bounds=tuple(sc.aot_bounds),
            tcwv_bounds=tuple(sc.tcwv_bounds),
            band_weight_power=getattr(sc, "alpha", -1.6),
            smoothness_delta=getattr(sc, "smoothness_delta", 0.02),
            grid_search_aot_points=getattr(sc, "grid_search_aot_points", 11),
            grid_search_tcwv_points=getattr(sc, "grid_search_tcwv_points", 11),
            fixed_atmospheric_parameter=getattr(sc, "fixed_atmospheric_parameter", "none"),
            stages=tuple(getattr(sc, "stages", ()) or ()),
            quadratic_block_size=getattr(sc, "quadratic_block_size", 1),
            quadratic_block_min_valid_fraction=getattr(
                sc, "quadratic_block_min_valid_fraction", 0.5
            ),
        )
        # ``solver_config`` is a MultiGridConfig instance (constructed above);
        # the previous ``isinstance(solver_config, dict)`` branch was unreachable
        # (REVIEW.md §3.6 _assembly_solver.py).
        solver_stages = getattr(solver_config, "stages", ())
        solver_cls = StagedMultiGridSolver if solver_stages else MultiGridSolver
        mg_solver = solver_cls(solver_config)
        solve_kwargs: dict[str, Any] = {}
        try:
            signature = inspect.signature(mg_solver.solve)
        except (TypeError, ValueError):
            # Wave 15: inspect.signature can fail on wrapped / C-implemented
            # callables. If we couldn't introspect we just won't forward the
            # optional kwargs. Note in the log so a debugger isn't confused
            # by sharp_transition_mask / water_mask silently dropping.
            import logging as _logging

            _logging.getLogger(__name__).warning(
                "Could not inspect solver signature; sharp_transition_mask "
                "and water_mask kwargs will not be forwarded. The solver "
                "may treat those regions less defensively than configured."
            )
        else:
            if "sharp_transition_mask" in signature.parameters:
                solve_kwargs["sharp_transition_mask"] = inputs.sharp_transition_mask
            if "water_mask" in signature.parameters:
                solve_kwargs["water_mask"] = inputs.water_mask
        result = mg_solver.solve(
            inputs.toa,
            inputs.surface_prior,
            inputs.geometry,
            inputs.atmo_prior,
            inputs.rt_model,
            inputs.cloud_mask,
            inputs.bands,
            **solve_kwargs,
        )
        return _wrap_solved(result, inputs)

    return _solver


@SOLVER_METHOD_REGISTRY.register("surface_driven")
def _build_surface_driven_solver(config: Any) -> SolverFn:
    from siac.algorithms.solver.surface_driven import SurfaceDrivenSolver

    solver = SurfaceDrivenSolver(config.algorithms.solver)

    def _solver(inputs: SolverInputBundle, _config: Any) -> SolvedAtmosphere:
        return _wrap_solved(solver.solve_bundle(inputs), inputs)

    return _solver


def resolve_solver(config: Any) -> SolverFn:
    method = getattr(config.algorithms.solver, "method", "multigrid")
    key = str(getattr(method, "value", method))
    builder = SOLVER_METHOD_REGISTRY.get(key)
    return builder(config)
