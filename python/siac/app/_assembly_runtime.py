"""Runtime-oriented assembly helpers for preprocessors, solver, and correction."""

from __future__ import annotations

import inspect
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from siac.adapters.output import ConfiguredOutputWriter
from siac.adapters.rt import build_rt_model
from siac.adapters.satellite import detect_sensor, get_preprocessor
from siac.algorithms.correction import AtmosphericCorrector, CorrectionResult
from siac.algorithms.solver import MultiGridConfig, MultiGridSolver, StagedMultiGridSolver
from siac.domain.aoi import AOI
from siac.geo.resample import (
    resample_field_for_correction as _resample_field_for_correction,
)
from siac.geo.resample import (
    resample_field_to_template as _resample_field_to_template,  # noqa: F401
)
from siac.geo.resample import (
    shares_template_grid as _shares_template_grid,  # noqa: F401
)
from siac.runtime import (
    AtmosphericState,
    GeometryAngles,
    ObservationBundle,
    SolvedAtmosphere,
    SolverInputBundle,
)

if TYPE_CHECKING:
    import xarray as xr

    from siac.domain.protocols import RTModelBackend
    from siac.domain.sensors import SensorConfig
    from siac.workflows.pipeline import CorrectorFn, GridAssemblerFn, PreprocessorFn, SolverFn


@dataclass(frozen=True)
class PreprocessorRuntime:
    """Preprocessor callable plus its sensor configuration."""

    preprocessor: PreprocessorFn
    sensor_config: SensorConfig


def _ordered_unique_band_names(names: list[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(name for name in names if name))


def _toa_preload_band_names(config: Any, sensor_config: SensorConfig) -> tuple[str, ...] | None:
    if not hasattr(sensor_config, "select_nearest_band") or not hasattr(sensor_config, "bands"):
        return None
    cloud_mask = getattr(config, "cloud_mask", None)
    cloud_mode = str(getattr(cloud_mask, "mode", "auto")).lower()
    if cloud_mode == "user_callable":
        return None

    names: list[str] = []
    ref_band = sensor_config.select_nearest_band(665.0, tolerance_nm=80.0)
    if ref_band is None:
        ref_band = sensor_config.bands[0]
    names.append(ref_band.name)

    names.extend(band.name for band in sensor_config.default_aerosol_solver_bands())

    if cloud_mode == "auto":
        for wl_min, wl_max in ((530.0, 590.0), (630.0, 690.0), (760.0, 900.0)):
            names.extend(
                band.name
                for band in sensor_config.select_bands_in_range(wl_min, wl_max)
            )

    surface_prior = getattr(config, "surface_prior", None)
    if str(getattr(surface_prior, "method", "")).lower() == "monthly_database":
        from siac.app._assembly_surface import _select_route_b_query_bands

        names.extend(
            band.name for band in _select_route_b_query_bands(sensor_config)
        )

    resolved = _ordered_unique_band_names(names)
    return resolved or None


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
    preload_toa_bands = _toa_preload_band_names(config, sensor_config)
    if hasattr(preprocessor_obj, "config") and isinstance(preprocessor_obj.config, dict):
        if rsrf_root is not None:
            preprocessor_obj.config.setdefault("rsrf_root", rsrf_root)
        preprocessor_obj.config.setdefault("cloud_mask", cloud_mask_config)
        if preload_toa_bands is not None:
            preprocessor_obj.config.setdefault("preload_toa_bands", preload_toa_bands)

    def _resolve_default_aoi(toa: xr.Dataset) -> AOI:
        if callable(default_aoi_resolver):
            return cast("AOI", default_aoi_resolver(toa))
        return AOI.from_raster(toa[list(toa.data_vars)[0]])

    def _clip_dataarray_to_aoi(field: object, *, bounds: tuple[float, float, float, float], crs: str) -> object:
        rio = getattr(field, "rio", None)
        if rio is None:
            return field
        try:
            return rio.clip_box(*bounds, crs=crs)
        except Exception:
            return field

    def _clip_raw_output(raw: dict[str, Any], *, aoi_obj: AOI) -> dict[str, Any]:
        import xarray as xr

        bounds = aoi_obj.get_bounds()
        crs = str(aoi_obj.crs)
        clipped = dict(raw)

        toa = raw.get("toa")
        if isinstance(toa, xr.Dataset):
            clipped_vars = {
                name: cast("xr.DataArray", _clip_dataarray_to_aoi(data, bounds=bounds, crs=crs))
                for name, data in toa.data_vars.items()
            }
            clipped["toa"] = xr.Dataset(clipped_vars, attrs=toa.attrs)

        geometry = raw.get("geometry")
        if geometry is not None:
            field_names = ("sza", "saa", "vza", "vaa")
            if all(hasattr(geometry, name) for name in field_names):
                clipped["geometry"] = geometry.__class__(
                    **{
                        name: _clip_dataarray_to_aoi(getattr(geometry, name), bounds=bounds, crs=crs)
                        for name in field_names
                    }
                )

        cloud_mask = raw.get("cloud_mask")
        clipped["cloud_mask"] = _clip_dataarray_to_aoi(cloud_mask, bounds=bounds, crs=crs)
        cloud_classes = raw.get("cloud_classes")
        if cloud_classes is not None:
            clipped["cloud_classes"] = _clip_dataarray_to_aoi(cloud_classes, bounds=bounds, crs=crs)
        return clipped

    def _preprocessor(path: Path, aoi: AOI | None = None) -> ObservationBundle:
        if aoi is not None:
            set_subset = getattr(preprocessor_obj, "set_spatial_subset", None)
            if callable(set_subset):
                set_subset(aoi.get_bounds(), crs=str(aoi.crs))
        try:
            raw = preprocessor_obj.preprocess(path)
        finally:
            clear_subset = getattr(preprocessor_obj, "clear_spatial_subset", None)
            if callable(clear_subset):
                clear_subset()
        toa = raw["toa"]
        resolved_aoi = aoi or _resolve_default_aoi(toa)
        if aoi is not None:
            raw = _clip_raw_output(raw, aoi_obj=resolved_aoi)
            toa = raw["toa"]
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
            band_weight_power=getattr(config.solver, "alpha", -1.6),
            smoothness_delta=getattr(config.solver, "smoothness_delta", 0.02),
            grid_search_aot_points=getattr(config.solver, "grid_search_aot_points", 11),
            grid_search_tcwv_points=getattr(config.solver, "grid_search_tcwv_points", 11),
            fixed_atmospheric_parameter=getattr(config.solver, "fixed_atmospheric_parameter", "none"),
            stages=tuple(getattr(config.solver, "stages", ()) or ()),
            quadratic_block_size=getattr(config.solver, "quadratic_block_size", 1),
            quadratic_block_min_valid_fraction=getattr(
                config.solver,
                "quadratic_block_min_valid_fraction",
                0.5,
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


def resolve_corrector(_config: Any) -> CorrectorFn:
    def _default_corrector(
        obs: ObservationBundle,
        solved: SolvedAtmosphere,
        rt_model: Any,
        output_stream: Any | None = None,
    ) -> CorrectionResult:
        corrector_obj = AtmosphericCorrector(rt_model, obs.sensor_config)
        atmo = solved.atmo_state
        matched_atmo = AtmosphericState(
            aot=_resample_field_for_correction(atmo.aot, atmo.aot, field_name="aot"),
            tcwv=_resample_field_for_correction(atmo.tcwv, atmo.aot, field_name="tcwv"),
            tco3=_resample_field_for_correction(atmo.tco3, atmo.aot, field_name="tco3"),
            aot_unc=_resample_field_for_correction(atmo.aot_unc, atmo.aot, field_name="aot_unc"),
            tcwv_unc=_resample_field_for_correction(atmo.tcwv_unc, atmo.aot, field_name="tcwv_unc"),
            tco3_unc=_resample_field_for_correction(atmo.tco3_unc, atmo.aot, field_name="tco3_unc"),
            elevation=_resample_field_for_correction(atmo.elevation, atmo.aot, field_name="elevation"),
        )
        coeff_geometry = GeometryAngles(
            sza=_resample_field_for_correction(obs.geometry.sza, atmo.aot, field_name="sza"),
            saa=_resample_field_for_correction(obs.geometry.saa, atmo.aot, field_name="saa"),
            vza=_resample_field_for_correction(obs.geometry.vza, atmo.aot, field_name="vza"),
            vaa=_resample_field_for_correction(obs.geometry.vaa, atmo.aot, field_name="vaa"),
        )
        correct_kwargs: dict[str, Any] = {}
        if output_stream is not None:
            with suppress(TypeError, ValueError):
                signature = inspect.signature(corrector_obj.correct)
                if "boa_band_writer" in signature.parameters or any(
                    param.kind == inspect.Parameter.VAR_KEYWORD
                    for param in signature.parameters.values()
                ):
                    correct_kwargs["boa_band_writer"] = output_stream.write_boa_band
        corrected = corrector_obj.correct(
            obs.toa,
            coeff_geometry,
            matched_atmo,
            obs.cloud_mask,
            **correct_kwargs,
        )
        if not isinstance(corrected, CorrectionResult):
            return corrected
        return CorrectionResult(
            boa=corrected.boa,
            boa_unc=corrected.boa_unc,
            aot=atmo.aot,
            tcwv=atmo.tcwv,
            cloud_mask=corrected.cloud_mask,
            surface_prior=corrected.surface_prior,
            surface_prior_unc=corrected.surface_prior_unc,
            solver_qa=corrected.solver_qa if corrected.solver_qa is not None else solved.qa,
            monthly_composites=corrected.monthly_composites,
            metadata=corrected.metadata,
            diagnostics=corrected.diagnostics,
        )

    return _default_corrector


def resolve_rt_model_for_pipeline(
    config: Any,
    auth: Any = None,
    *,
    sensor_config: SensorConfig | None = None,
) -> RTModelBackend:
    return build_rt_model(config, auth=auth, sensor_config=sensor_config)  # type: ignore[no-any-return]


__all__ = [
    "PreprocessorRuntime",
    "build_preprocessor_runtime",
    "resolve_corrector",
    "resolve_grid_assembler",
    "resolve_output_writer",
    "resolve_preprocessor",
    "resolve_rt_model_for_pipeline",
    "resolve_solver",
]
