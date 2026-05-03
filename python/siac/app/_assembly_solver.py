"""Solver assembly for the atmospheric inversion stage."""

from __future__ import annotations

import inspect
from contextlib import suppress
from typing import TYPE_CHECKING, Any

from siac.algorithms.solver import MultiGridConfig, MultiGridSolver, StagedMultiGridSolver
from siac.runtime import SolvedAtmosphere, SolverInputBundle

if TYPE_CHECKING:
    from siac.workflows.pipeline import SolverFn


def resolve_solver(config: Any) -> SolverFn:
    def _default_solver(inputs: SolverInputBundle, _config: Any) -> SolvedAtmosphere:
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
        solver_stages = (
            solver_config.get("stages", ())
            if isinstance(solver_config, dict)
            else getattr(solver_config, "stages", ())
        )
        solver_cls = StagedMultiGridSolver if solver_stages else MultiGridSolver
        mg_solver = solver_cls(solver_config)
        solve_kwargs: dict[str, Any] = {}
        with suppress(TypeError, ValueError):
            signature = inspect.signature(mg_solver.solve)
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
        )

    return _default_solver
