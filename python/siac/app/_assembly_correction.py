"""Correction assembly for the M6 atmospheric-correction stage."""

from __future__ import annotations

import inspect
from contextlib import suppress
from typing import TYPE_CHECKING, Any

from siac.algorithms.correction import AtmosphericCorrector, CorrectionResult
from siac.geo.resample import resample_field_for_correction as _resample_field_for_correction
from siac.runtime import AtmosphericState, GeometryAngles, ObservationBundle, SolvedAtmosphere

if TYPE_CHECKING:
    from siac.workflows.pipeline import CorrectorFn


def _callable_supports_kwarg(target: Any, keyword: str) -> bool:
    with suppress(TypeError, ValueError):
        signature = inspect.signature(target)
        return keyword in signature.parameters or any(
            param.kind == inspect.Parameter.VAR_KEYWORD for param in signature.parameters.values()
        )
    return False


def resolve_corrector(config: Any) -> CorrectorFn:
    execution = getattr(getattr(config, "runtime", None), "execution", None)
    raw_correction_workers = getattr(execution, "correction_max_workers", None)
    if raw_correction_workers is None:
        raw_correction_workers = (
            getattr(execution, "max_workers", 1) if execution is not None else 1
        )
    correction_workers = max(1, int(raw_correction_workers or 1))

    def _default_corrector(
        obs: ObservationBundle,
        solved: SolvedAtmosphere,
        rt_model: Any,
        output_stream: Any | None = None,
    ) -> CorrectionResult:
        init_kwargs: dict[str, Any] = {}
        if _callable_supports_kwarg(AtmosphericCorrector, "correction_workers"):
            init_kwargs["correction_workers"] = correction_workers
        corrector_obj = AtmosphericCorrector(rt_model, obs.sensor_config, **init_kwargs)
        atmo = solved.atmo_state
        matched_atmo = AtmosphericState(
            aot=_resample_field_for_correction(atmo.aot, atmo.aot, field_name="aot"),
            tcwv=_resample_field_for_correction(atmo.tcwv, atmo.aot, field_name="tcwv"),
            tco3=_resample_field_for_correction(atmo.tco3, atmo.aot, field_name="tco3"),
            aot_unc=_resample_field_for_correction(atmo.aot_unc, atmo.aot, field_name="aot_unc"),
            tcwv_unc=_resample_field_for_correction(atmo.tcwv_unc, atmo.aot, field_name="tcwv_unc"),
            tco3_unc=_resample_field_for_correction(atmo.tco3_unc, atmo.aot, field_name="tco3_unc"),
            elevation=_resample_field_for_correction(
                atmo.elevation, atmo.aot, field_name="elevation"
            ),
        )
        coeff_geometry = GeometryAngles(
            sza=_resample_field_for_correction(obs.geometry.sza, atmo.aot, field_name="sza"),
            saa=_resample_field_for_correction(obs.geometry.saa, atmo.aot, field_name="saa"),
            vza=_resample_field_for_correction(obs.geometry.vza, atmo.aot, field_name="vza"),
            vaa=_resample_field_for_correction(obs.geometry.vaa, atmo.aot, field_name="vaa"),
        )
        correct_kwargs: dict[str, Any] = {}
        if output_stream is not None and _callable_supports_kwarg(
            corrector_obj.correct, "boa_band_writer"
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
