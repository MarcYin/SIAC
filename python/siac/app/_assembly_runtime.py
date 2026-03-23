"""Runtime-oriented assembly helpers for preprocessors, solver, and correction."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from siac.adapters.output import ConfiguredOutputWriter
from siac.adapters.rt import build_rt_model
from siac.adapters.satellite import detect_sensor, get_preprocessor
from siac.algorithms.correction import AtmosphericCorrector, CorrectionResult
from siac.algorithms.solver import MultiGridConfig, MultiGridSolver
from siac.domain.aoi import AOI
from siac.runtime import AtmosphericState, ObservationBundle, SolvedAtmosphere, SolverInputBundle

if TYPE_CHECKING:
    import xarray as xr

    from siac.domain.sensors import SensorConfig
    from siac.workflows.pipeline import CorrectorFn, GridAssemblerFn, PreprocessorFn, SolverFn


@dataclass(frozen=True)
class PreprocessorRuntime:
    """Preprocessor callable plus its sensor configuration."""

    preprocessor: PreprocessorFn
    sensor_config: SensorConfig


def build_preprocessor_runtime(
    config: Any,
    *,
    input_path: Path | None = None,
    sensor: str | None = None,
    default_aoi_resolver: Any | None = None,
    detect_sensor_fn: Any | None = None,
    get_preprocessor_fn: Any | None = None,
) -> PreprocessorRuntime:
    """Build the M1 preprocessor runtime."""
    if detect_sensor_fn is None:
        detect_sensor_fn = detect_sensor
    if get_preprocessor_fn is None:
        get_preprocessor_fn = get_preprocessor

    sensor_name = sensor or config.sensor
    if sensor_name == "auto":
        if input_path is None:
            raise ValueError("Cannot resolve preprocessor for sensor='auto' without an input path.")
        sensor_name = detect_sensor_fn(input_path)

    cloud_mask_config = config.cloud_mask.model_dump(exclude={"user_callable"})
    paths = getattr(config, "paths", None)
    rsrf_root = getattr(paths, "rsrf_root", None)
    preprocessor_config = {"cloud_mask": cloud_mask_config}
    if rsrf_root is not None:
        preprocessor_config["rsrf_root"] = rsrf_root
    try:
        preprocessor_obj = get_preprocessor_fn(sensor_name, config=preprocessor_config)
    except KeyError as exc:
        raise ValueError(f"Unknown sensor: {sensor_name!r}") from exc
    except TypeError:
        try:
            preprocessor_obj = get_preprocessor_fn(sensor_name)
        except KeyError as exc:
            raise ValueError(f"Unknown sensor: {sensor_name!r}") from exc
        if hasattr(preprocessor_obj, "config") and isinstance(preprocessor_obj.config, dict):
            preprocessor_obj.config.setdefault("cloud_mask", cloud_mask_config)
            if rsrf_root is not None:
                preprocessor_obj.config.setdefault("rsrf_root", rsrf_root)

    if input_path is not None and hasattr(preprocessor_obj, "get_metadata"):
        preprocessor_obj.get_metadata(Path(input_path))

    sensor_config = preprocessor_obj.sensor_config

    def _resolve_default_aoi(toa: xr.Dataset) -> AOI:
        if callable(default_aoi_resolver):
            return cast("AOI", default_aoi_resolver(toa))
        return AOI.from_raster(toa[list(toa.data_vars)[0]])

    def _preprocessor(path: Path, aoi: AOI | None = None) -> ObservationBundle:
        raw = preprocessor_obj.preprocess(path)
        toa = raw["toa"]
        resolved_aoi = aoi or _resolve_default_aoi(toa)
        active_sensor_config = getattr(preprocessor_obj, "sensor_config", sensor_config)
        return ObservationBundle(
            toa=toa,
            geometry=raw["geometry"],
            cloud_mask=raw["cloud_mask"],
            sensor_config=active_sensor_config,
            metadata=raw["metadata"],
            crs=str(resolved_aoi.crs),
            bounds=resolved_aoi.get_bounds(),
        )

    return PreprocessorRuntime(preprocessor=_preprocessor, sensor_config=sensor_config)


def resolve_preprocessor(
    config: Any,
    *,
    input_path: Path | None = None,
    sensor: str | None = None,
    default_aoi_resolver: Any | None = None,
) -> PreprocessorFn:
    return build_preprocessor_runtime(
        config,
        input_path=input_path,
        sensor=sensor,
        default_aoi_resolver=default_aoi_resolver,
    ).preprocessor


def resolve_output_writer(config: Any) -> ConfiguredOutputWriter:
    return ConfiguredOutputWriter(config.output.defaults)


def resolve_grid_assembler() -> GridAssemblerFn:
    from siac.algorithms.grid.assembler import assemble_grids

    return assemble_grids


def resolve_solver(config: Any) -> SolverFn:
    def _default_solver(inputs: SolverInputBundle, _config: Any) -> SolvedAtmosphere:
        solver_config = MultiGridConfig(
            aot_gamma=config.solver.aot_gamma,
            tcwv_gamma=config.solver.tcwv_gamma,
            aot_bounds=tuple(config.solver.aot_bounds),
            tcwv_bounds=tuple(config.solver.tcwv_bounds),
        )
        mg_solver = MultiGridSolver(solver_config)
        result = mg_solver.solve(
            inputs.toa,
            inputs.surface_prior,
            inputs.geometry,
            inputs.atmo_prior,
            inputs.rt_model,
            inputs.cloud_mask,
            inputs.bands,
        )
        solved_atmo = inputs.atmo_prior.with_updated_aot_tcwv(
            aot=result.aot,
            tcwv=result.tcwv,
            aot_unc=result.aot_unc,
            tcwv_unc=result.tcwv_unc,
        )
        return SolvedAtmosphere(
            atmo_state=solved_atmo,
            aot=result.aot,
            tcwv=result.tcwv,
            aot_unc=result.aot_unc,
            tcwv_unc=result.tcwv_unc,
            cost_final=float(result.final_cost),
            n_iterations=result.n_iterations,
            converged=result.success,
        )

    return _default_solver


def _shares_template_grid(field: Any, template: Any) -> bool:
    """Return True when a field already matches the template's spatial grid."""
    if tuple(getattr(field, "shape", ())) != tuple(getattr(template, "shape", ())):
        return False

    field_dims = tuple(getattr(field, "dims", ()))
    template_dims = tuple(getattr(template, "dims", ()))
    if field_dims != template_dims:
        return False

    field_coords = getattr(field, "coords", {})
    template_coords = getattr(template, "coords", {})
    for axis in template_dims:
        field_has_coord = axis in field_coords
        template_has_coord = axis in template_coords
        if field_has_coord != template_has_coord:
            return False
        if not field_has_coord:
            continue
        field_values = np.asarray(field.coords[axis].values)
        template_values = np.asarray(template.coords[axis].values)
        if not np.array_equal(field_values, template_values, equal_nan=True):
            return False
    return True


def _resample_field_to_template(field: Any, template: Any) -> Any:
    if _shares_template_grid(field, template):
        return field
    if (
        len(field.dims) == 2
        and field.dims == template.dims
        and all(dim in field.coords for dim in field.dims)
        and all(dim in template.coords for dim in template.dims)
    ):
        try:
            return field.interp(
                coords={dim: template.coords[dim] for dim in template.dims},
                method="linear",
            )
        except Exception:
            pass

    src = np.asarray(field.values, dtype=np.float32)
    if src.ndim != 2 or len(template.dims) != 2:
        return field

    from scipy import ndimage

    h_out = int(template.sizes[template.dims[0]])
    w_out = int(template.sizes[template.dims[1]])
    if src.shape[0] == 0 or src.shape[1] == 0:
        out: np.ndarray[Any, Any] = np.full((h_out, w_out), np.nan, dtype=np.float32)
    else:
        out = ndimage.zoom(src, (h_out / src.shape[0], w_out / src.shape[1]), order=1)
        out = out[:h_out, :w_out]
        if out.shape != (h_out, w_out):
            padded: np.ndarray[Any, Any] = np.full((h_out, w_out), np.nan, dtype=np.float32)
            padded[: out.shape[0], : out.shape[1]] = out
            out = padded

    return field.__class__(
        out,
        dims=template.dims,
        coords={d: template.coords[d] for d in template.dims if d in template.coords},
    )


def resolve_corrector(_config: Any) -> CorrectorFn:
    def _default_corrector(
        obs: ObservationBundle,
        solved: SolvedAtmosphere,
        rt_model: Any,
    ) -> CorrectionResult:
        corrector_obj = AtmosphericCorrector(rt_model, obs.sensor_config)
        first_band = obs.toa[next(iter(obs.toa.data_vars))]
        atmo = solved.atmo_state
        matched_atmo = AtmosphericState(
            aot=_resample_field_to_template(atmo.aot, first_band),
            tcwv=_resample_field_to_template(atmo.tcwv, first_band),
            tco3=_resample_field_to_template(atmo.tco3, first_band),
            aot_unc=_resample_field_to_template(atmo.aot_unc, first_band),
            tcwv_unc=_resample_field_to_template(atmo.tcwv_unc, first_band),
            tco3_unc=_resample_field_to_template(atmo.tco3_unc, first_band),
            elevation=_resample_field_to_template(atmo.elevation, first_band),
        )
        return corrector_obj.correct(obs.toa, obs.geometry, matched_atmo, obs.cloud_mask)

    return _default_corrector


def resolve_rt_model_for_pipeline(
    config: Any,
    auth: Any = None,
    *,
    sensor_config: SensorConfig | None = None,
) -> Any:
    return build_rt_model(config, auth=auth, sensor_config=sensor_config)


__all__ = [
    "PreprocessorRuntime",
    "_resample_field_to_template",
    "_shares_template_grid",
    "build_preprocessor_runtime",
    "resolve_corrector",
    "resolve_grid_assembler",
    "resolve_output_writer",
    "resolve_preprocessor",
    "resolve_rt_model_for_pipeline",
    "resolve_solver",
]
