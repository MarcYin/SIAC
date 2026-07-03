"""Scene-processing workflows."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

from siac.app.planning import build_execution_plan
from siac.observability import derive_execution_report_path
from siac.workflows.pipeline import (
    _OUTPUTS_WRITTEN_METADATA_KEY,
    _resolve_execution_settings,
    run_pipeline,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from siac.app.planning import ExecutionPlan
    from siac.app.requests import SceneProcessRequest
    from siac.domain.protocols import OutputWriter
    from siac.runtime import CorrectionResult
    from siac.workflows.pipeline import SurfacePriorFn


AODRailStatistic = Literal["max", "p95", "mean", "median"]


def _finite_metadata_float(value: Any, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def write_output(
    result: CorrectionResult,
    output_path: Path | str,
    *,
    output_writer: OutputWriter,
) -> None:
    """Persist SIAC outputs through the configured output adapter."""
    output_writer.write(result, Path(output_path))


def _finalize_scene_result(
    result: CorrectionResult,
    *,
    plan: ExecutionPlan,
) -> CorrectionResult:
    if (
        plan.output_path is not None
        and plan.output_writer is not None
        and not getattr(result, "metadata", {}).get(_OUTPUTS_WRITTEN_METADATA_KEY, False)
    ):
        write_output(result, plan.output_path, output_writer=plan.output_writer)
    return result


def execute_plan(
    plan: ExecutionPlan,
    *,
    execution: Any | None = None,
) -> CorrectionResult:
    """Execute a fully resolved :class:`ExecutionPlan`."""
    return run_pipeline(
        plan.input_path,
        plan.runtime_aoi,
        plan.config,
        preprocessor=plan.preprocessor,
        atmo_provider=plan.atmo_provider,
        surface_prior_provider=plan.surface_prior_provider,
        grid_assembler=plan.grid_assembler,
        solver=plan.solver,
        corrector=plan.corrector,
        rt_model=plan.rt_model,
        execution=execution,
        output_path=getattr(plan, "output_path", None),
        output_writer=getattr(plan, "output_writer", None),
    )


def process_scene(
    request: SceneProcessRequest,
    *,
    preprocessor: Any | None = None,
    atmo_provider: Any | None = None,
    surface_prior_provider: Any | None = None,
    grid_assembler: Any | None = None,
    solver: Any | None = None,
    corrector: Any | None = None,
    rt_model: Any | None = None,
    aoi_resolver: Any | None = None,
    build_preprocessor_runtime_fn: Any | None = None,
    resolve_atmo_provider_fn: Any | None = None,
    resolve_surface_prior_provider_fn: Any | None = None,
    resolve_grid_assembler_fn: Any | None = None,
    resolve_solver_fn: Any | None = None,
    resolve_corrector_fn: Any | None = None,
    resolve_rt_model_fn: Any | None = None,
) -> CorrectionResult:
    """Build a runtime plan and execute it for one scene."""
    plan = build_execution_plan(
        request,
        preprocessor=preprocessor,
        atmo_provider=atmo_provider,
        surface_prior_provider=surface_prior_provider,
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
        aoi_resolver=aoi_resolver,
        build_preprocessor_runtime_fn=build_preprocessor_runtime_fn,
        resolve_atmo_provider_fn=resolve_atmo_provider_fn,
        resolve_surface_prior_provider_fn=resolve_surface_prior_provider_fn,
        resolve_grid_assembler_fn=resolve_grid_assembler_fn,
        resolve_solver_fn=resolve_solver_fn,
        resolve_corrector_fn=resolve_corrector_fn,
        resolve_rt_model_fn=resolve_rt_model_fn,
    )
    execution_override: dict[str, Any] | None = None
    settings_config = getattr(plan, "config", request.config)
    settings = _resolve_execution_settings(settings_config, execution=None, max_workers=None)
    if settings["show_progress"] and settings["performance_report_path"] is None:
        default_report_path = derive_execution_report_path(plan.output_path)
        if default_report_path is not None:
            execution_override = {"performance_report_path": default_report_path}

    result = execute_plan(plan, execution=execution_override)
    if plan.output_path is not None and plan.output_writer is not None:
        result = _finalize_scene_result(result, plan=plan)
    return result


def aod_rail_value(
    result: CorrectionResult,
    *,
    statistic: AODRailStatistic = "max",
) -> float:
    """Return a finite AOD summary value used by rail-fallback orchestration."""
    values = np.asarray(result.aot.values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    if statistic == "max":
        return float(np.max(finite))
    if statistic == "p95":
        return float(np.nanpercentile(finite, 95.0))
    if statistic == "mean":
        return float(np.mean(finite))
    if statistic == "median":
        return float(np.median(finite))
    raise ValueError(f"Unsupported AOD rail statistic {statistic!r}")


def aod_scene_mean(result: CorrectionResult) -> float:
    """Return the finite scene/AOI mean AOD from a processed result."""
    return aod_rail_value(result, statistic="mean")


def aod_quality_score(
    result: CorrectionResult,
    *,
    cost_weight: float = 1.0,
    uncertainty_weight: float = 1.0,
    spatial_weight: float = 0.25,
    non_converged_penalty: float = 10.0,
) -> float:
    """Return a generic lower-is-better quality score for an AOD candidate."""
    solver = result.metadata.get("solver", {})
    if not isinstance(solver, dict):
        solver = {}
    cost = max(
        _finite_metadata_float(
            solver.get("cost_final_per_band", solver.get("cost_final")),
            default=np.inf,
        ),
        0.0,
    )
    uncertainty = max(_finite_metadata_float(solver.get("aot_unc_median")), 0.0)
    spatial = max(_finite_metadata_float(solver.get("aot_std")), 0.0)
    finite_fraction = _finite_metadata_float(
        solver.get("aot_finite_fraction"),
        default=1.0,
    )
    score = (
        cost_weight * float(np.log1p(cost))
        + uncertainty_weight * uncertainty
        + spatial_weight * spatial
        + max(0.0, 1.0 - finite_fraction) * non_converged_penalty
    )
    if not bool(solver.get("converged", True)):
        score += non_converged_penalty
    return float(score)


def _metadata_with_aod_rail_fallback(
    result: CorrectionResult,
    *,
    enabled: bool,
    triggered: bool,
    threshold: float,
    statistic_name: str,
    statistic_value: float,
    source: str,
) -> CorrectionResult:
    return replace(
        result,
        metadata={
            **result.metadata,
            "aod_rail_fallback": {
                "enabled": enabled,
                "triggered": triggered,
                "threshold": threshold,
                "statistic": statistic_name,
                "statistic_value": statistic_value,
                "source": source,
            },
        },
    )


def _metadata_with_aod_selector(
    result: CorrectionResult,
    *,
    selected: str,
    selector_name: str,
    decision_metadata: dict[str, Any],
) -> CorrectionResult:
    return replace(
        result,
        metadata={
            **result.metadata,
            "aod_selector": {
                "enabled": True,
                "selected": selected,
                "selector": selector_name,
                **decision_metadata,
            },
        },
    )


def _candidate_aod_statistics(result: CorrectionResult) -> dict[str, float | None]:
    quality_score = aod_quality_score(result)
    return {
        "median": aod_rail_value(result, statistic="median"),
        "mean": aod_rail_value(result, statistic="mean"),
        "p95": aod_rail_value(result, statistic="p95"),
        "max": aod_rail_value(result, statistic="max"),
        "quality_score": quality_score if np.isfinite(quality_score) else None,
    }


def process_scene_with_aod_ensemble(
    request: SceneProcessRequest,
    *,
    candidate_surface_prior_providers: dict[str, SurfacePriorFn | None],
    candidate_requests: dict[str, SceneProcessRequest] | None = None,
    candidate_weights: dict[str, float] | None = None,
    ensemble_method: Literal["mean", "power_mean", "median"] | None = None,
    ensemble_power: float = 1.0,
    preprocessor: Any | None = None,
    atmo_provider: Any | None = None,
    grid_assembler: Any | None = None,
    solver: Any | None = None,
    corrector: Any | None = None,
    rt_model: Any | None = None,
    aoi_resolver: Any | None = None,
    build_preprocessor_runtime_fn: Any | None = None,
    resolve_atmo_provider_fn: Any | None = None,
    resolve_surface_prior_provider_fn: Any | None = None,
    resolve_grid_assembler_fn: Any | None = None,
    resolve_solver_fn: Any | None = None,
    resolve_corrector_fn: Any | None = None,
    resolve_rt_model_fn: Any | None = None,
) -> CorrectionResult:
    """Run named AOD candidates and return a fixed AOD ensemble.

    This is a generic fusion primitive: it does not inspect AERONET truth and it
    does not switch by AOD thresholds. Callers provide the candidate surface
    prior providers, optionally with per-candidate requests/configs, and this
    workflow combines the candidate AOD rasters with a fixed arithmetic mean,
    power mean, or per-pixel median.
    """
    if not candidate_surface_prior_providers:
        raise ValueError("candidate_surface_prior_providers must contain at least one candidate.")
    if ensemble_method is None:
        method = "mean" if abs(float(ensemble_power) - 1.0) < 1.0e-9 else "power_mean"
    else:
        method = ensemble_method
    if method not in {"mean", "power_mean", "median"}:
        raise ValueError(f"Unsupported ensemble_method {method!r}.")
    if float(ensemble_power) <= 0.0:
        raise ValueError("ensemble_power must be positive.")
    names = list(candidate_surface_prior_providers)
    if method == "median" and candidate_weights is not None:
        raise ValueError("candidate_weights are not supported with ensemble_method='median'.")
    if method == "median":
        weights: dict[str, float] | None = None
    elif candidate_weights is None:
        weights = {name: 1.0 / float(len(names)) for name in names}
    else:
        missing_weights = [name for name in names if name not in candidate_weights]
        if missing_weights:
            raise ValueError(f"candidate_weights missing entries for {missing_weights!r}.")
        weights = {name: float(candidate_weights[name]) for name in names}
        if any(weight < 0.0 for weight in weights.values()):
            raise ValueError("candidate_weights must be non-negative.")
        weight_sum = float(sum(weights.values()))
        if weight_sum <= 0.0:
            raise ValueError("candidate_weights must contain a positive total weight.")
        weights = {name: weight / weight_sum for name, weight in weights.items()}

    common_kwargs: dict[str, Any] = {
        "preprocessor": preprocessor,
        "atmo_provider": atmo_provider,
        "grid_assembler": grid_assembler,
        "solver": solver,
        "corrector": corrector,
        "rt_model": rt_model,
        "aoi_resolver": aoi_resolver,
        "build_preprocessor_runtime_fn": build_preprocessor_runtime_fn,
        "resolve_atmo_provider_fn": resolve_atmo_provider_fn,
        "resolve_surface_prior_provider_fn": resolve_surface_prior_provider_fn,
        "resolve_grid_assembler_fn": resolve_grid_assembler_fn,
        "resolve_solver_fn": resolve_solver_fn,
        "resolve_corrector_fn": resolve_corrector_fn,
        "resolve_rt_model_fn": resolve_rt_model_fn,
    }

    candidate_requests = candidate_requests or {}
    candidate_results: dict[str, CorrectionResult] = {}
    for name, provider in candidate_surface_prior_providers.items():
        candidate_request = candidate_requests.get(name, request)
        candidate_results[name] = process_scene(
            replace(candidate_request, output_path=None),
            surface_prior_provider=provider,
            **common_kwargs,
        )

    first_name = names[0]
    first_result = candidate_results[first_name]
    if method == "median":
        aot_stack = np.stack(
            [np.asarray(candidate_results[name].aot.values, dtype=np.float64) for name in names],
            axis=0,
        )
        ensemble_aot = first_result.aot.copy(data=np.nanmedian(aot_stack, axis=0))
    elif method == "mean":
        assert weights is not None
        ensemble_aot = first_result.aot * weights[first_name]
        for name in names[1:]:
            ensemble_aot = ensemble_aot + candidate_results[name].aot * weights[name]
    else:
        assert weights is not None
        powered = (first_result.aot ** float(ensemble_power)) * weights[first_name]
        for name in names[1:]:
            powered = (
                powered + (candidate_results[name].aot ** float(ensemble_power)) * weights[name]
            )
        ensemble_aot = powered ** (1.0 / float(ensemble_power))

    template = first_result
    result = replace(
        template,
        aot=ensemble_aot.astype(template.aot.dtype),
        metadata={
            **template.metadata,
            "aod_ensemble": {
                "enabled": True,
                "method": method,
                "power": float(ensemble_power),
                "candidates": names,
                "weights": weights,
                "candidate_statistics": {
                    name: _candidate_aod_statistics(candidate)
                    for name, candidate in candidate_results.items()
                },
            },
        },
    )
    if request.output_path is None:
        return result

    output_plan = build_execution_plan(
        request,
        preprocessor=preprocessor,
        atmo_provider=atmo_provider,
        surface_prior_provider=next(iter(candidate_surface_prior_providers.values())),
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
        aoi_resolver=aoi_resolver,
        build_preprocessor_runtime_fn=build_preprocessor_runtime_fn,
        resolve_atmo_provider_fn=resolve_atmo_provider_fn,
        resolve_surface_prior_provider_fn=resolve_surface_prior_provider_fn,
        resolve_grid_assembler_fn=resolve_grid_assembler_fn,
        resolve_solver_fn=resolve_solver_fn,
        resolve_corrector_fn=resolve_corrector_fn,
        resolve_rt_model_fn=resolve_rt_model_fn,
    )
    return _finalize_scene_result(result, plan=output_plan)


def process_scene_with_aod_selector(
    request: SceneProcessRequest,
    *,
    fallback_surface_prior_provider: SurfacePriorFn,
    selector: Callable[[CorrectionResult, CorrectionResult], bool | tuple[bool, dict[str, Any]]],
    fallback_request: SceneProcessRequest | None = None,
    preprocessor: Any | None = None,
    atmo_provider: Any | None = None,
    surface_prior_provider: Any | None = None,
    grid_assembler: Any | None = None,
    solver: Any | None = None,
    corrector: Any | None = None,
    rt_model: Any | None = None,
    aoi_resolver: Any | None = None,
    build_preprocessor_runtime_fn: Any | None = None,
    resolve_atmo_provider_fn: Any | None = None,
    resolve_surface_prior_provider_fn: Any | None = None,
    resolve_grid_assembler_fn: Any | None = None,
    resolve_solver_fn: Any | None = None,
    resolve_corrector_fn: Any | None = None,
    resolve_rt_model_fn: Any | None = None,
) -> CorrectionResult:
    """Run base and fallback passes, then return the selector-chosen result."""
    common_kwargs: dict[str, Any] = {
        "preprocessor": preprocessor,
        "atmo_provider": atmo_provider,
        "grid_assembler": grid_assembler,
        "solver": solver,
        "corrector": corrector,
        "rt_model": rt_model,
        "aoi_resolver": aoi_resolver,
        "build_preprocessor_runtime_fn": build_preprocessor_runtime_fn,
        "resolve_atmo_provider_fn": resolve_atmo_provider_fn,
        "resolve_surface_prior_provider_fn": resolve_surface_prior_provider_fn,
        "resolve_grid_assembler_fn": resolve_grid_assembler_fn,
        "resolve_solver_fn": resolve_solver_fn,
        "resolve_corrector_fn": resolve_corrector_fn,
        "resolve_rt_model_fn": resolve_rt_model_fn,
    }
    fallback_request = request if fallback_request is None else fallback_request

    base_result = process_scene(
        replace(request, output_path=None),
        surface_prior_provider=surface_prior_provider,
        **common_kwargs,
    )
    fallback_result = process_scene(
        replace(fallback_request, output_path=None),
        surface_prior_provider=fallback_surface_prior_provider,
        **common_kwargs,
    )

    decision = selector(base_result, fallback_result)
    if isinstance(decision, tuple):
        use_fallback, decision_metadata = decision
    else:
        use_fallback = decision
        decision_metadata = {}
    use_fallback = bool(use_fallback)
    selected = "fallback" if use_fallback else "base"
    selected_result = fallback_result if use_fallback else base_result
    selected_request = fallback_request if use_fallback else request
    selected_provider = fallback_surface_prior_provider if use_fallback else surface_prior_provider
    selector_name = getattr(selector, "__name__", "callable")

    result = _metadata_with_aod_selector(
        selected_result,
        selected=selected,
        selector_name=selector_name,
        decision_metadata=decision_metadata,
    )
    if request.output_path is None:
        return result

    output_plan = build_execution_plan(
        replace(selected_request, output_path=request.output_path),
        preprocessor=preprocessor,
        atmo_provider=atmo_provider,
        surface_prior_provider=selected_provider,
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
        aoi_resolver=aoi_resolver,
        build_preprocessor_runtime_fn=build_preprocessor_runtime_fn,
        resolve_atmo_provider_fn=resolve_atmo_provider_fn,
        resolve_surface_prior_provider_fn=resolve_surface_prior_provider_fn,
        resolve_grid_assembler_fn=resolve_grid_assembler_fn,
        resolve_solver_fn=resolve_solver_fn,
        resolve_corrector_fn=resolve_corrector_fn,
        resolve_rt_model_fn=resolve_rt_model_fn,
    )
    return _finalize_scene_result(result, plan=output_plan)


def process_scene_with_aod_rail_fallback(
    request: SceneProcessRequest,
    *,
    fallback_surface_prior_provider: SurfacePriorFn,
    aod_rail_threshold: float = 1.7,
    aod_rail_statistic: AODRailStatistic | Callable[[CorrectionResult], float] = "max",
    preprocessor: Any | None = None,
    atmo_provider: Any | None = None,
    surface_prior_provider: Any | None = None,
    grid_assembler: Any | None = None,
    solver: Any | None = None,
    corrector: Any | None = None,
    rt_model: Any | None = None,
    aoi_resolver: Any | None = None,
    build_preprocessor_runtime_fn: Any | None = None,
    resolve_atmo_provider_fn: Any | None = None,
    resolve_surface_prior_provider_fn: Any | None = None,
    resolve_grid_assembler_fn: Any | None = None,
    resolve_solver_fn: Any | None = None,
    resolve_corrector_fn: Any | None = None,
    resolve_rt_model_fn: Any | None = None,
) -> CorrectionResult:
    """Run a production two-pass AOD rail fallback.

    The first pass uses the normal surface prior. If the retrieved AOD summary
    reaches ``aod_rail_threshold``, the scene is rerun with
    ``fallback_surface_prior_provider`` and the fallback result is returned.
    This is an orchestration primitive; it does not construct the fallback
    prior itself.
    """
    if aod_rail_threshold <= 0.0:
        raise ValueError("aod_rail_threshold must be positive.")

    common_kwargs: dict[str, Any] = {
        "preprocessor": preprocessor,
        "atmo_provider": atmo_provider,
        "grid_assembler": grid_assembler,
        "solver": solver,
        "corrector": corrector,
        "rt_model": rt_model,
        "aoi_resolver": aoi_resolver,
        "build_preprocessor_runtime_fn": build_preprocessor_runtime_fn,
        "resolve_atmo_provider_fn": resolve_atmo_provider_fn,
        "resolve_surface_prior_provider_fn": resolve_surface_prior_provider_fn,
        "resolve_grid_assembler_fn": resolve_grid_assembler_fn,
        "resolve_solver_fn": resolve_solver_fn,
        "resolve_corrector_fn": resolve_corrector_fn,
        "resolve_rt_model_fn": resolve_rt_model_fn,
    }

    # Keep the exploratory first pass from writing user-visible outputs. Only
    # the selected result is persisted below.
    base_request = replace(request, output_path=None)
    base_result = process_scene(
        base_request,
        surface_prior_provider=surface_prior_provider,
        **common_kwargs,
    )
    if callable(aod_rail_statistic):
        rail_value = float(aod_rail_statistic(base_result))
        statistic_name = getattr(aod_rail_statistic, "__name__", "callable")
    else:
        rail_value = aod_rail_value(base_result, statistic=aod_rail_statistic)
        statistic_name = aod_rail_statistic
    triggered = bool(np.isfinite(rail_value) and rail_value >= aod_rail_threshold)

    if triggered:
        fallback_result = process_scene(
            replace(request, output_path=None),
            surface_prior_provider=fallback_surface_prior_provider,
            **common_kwargs,
        )
        result = _metadata_with_aod_rail_fallback(
            fallback_result,
            enabled=True,
            triggered=True,
            threshold=aod_rail_threshold,
            statistic_name=statistic_name,
            statistic_value=rail_value,
            source="fallback",
        )
    else:
        result = _metadata_with_aod_rail_fallback(
            base_result,
            enabled=True,
            triggered=False,
            threshold=aod_rail_threshold,
            statistic_name=statistic_name,
            statistic_value=rail_value,
            source="base",
        )

    if request.output_path is None:
        return result

    output_plan = build_execution_plan(
        request,
        preprocessor=preprocessor,
        atmo_provider=atmo_provider,
        surface_prior_provider=(
            fallback_surface_prior_provider if triggered else surface_prior_provider
        ),
        grid_assembler=grid_assembler,
        solver=solver,
        corrector=corrector,
        rt_model=rt_model,
        aoi_resolver=aoi_resolver,
        build_preprocessor_runtime_fn=build_preprocessor_runtime_fn,
        resolve_atmo_provider_fn=resolve_atmo_provider_fn,
        resolve_surface_prior_provider_fn=resolve_surface_prior_provider_fn,
        resolve_grid_assembler_fn=resolve_grid_assembler_fn,
        resolve_solver_fn=resolve_solver_fn,
        resolve_corrector_fn=resolve_corrector_fn,
        resolve_rt_model_fn=resolve_rt_model_fn,
    )
    return _finalize_scene_result(result, plan=output_plan)
