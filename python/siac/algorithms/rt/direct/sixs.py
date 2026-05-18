"""Native 6SV2.1 RT backend."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

from siac.algorithms.rt.direct.sixs_native import (
    JointGridSearchLUT,
    PreparedCorrectionScene,
    SixSNativeRunner,
)
from siac.rt_setup import resolve_backend_rt_setup
from siac.runtime import RTCoefficients

if TYPE_CHECKING:
    from datetime import datetime

    import numpy as np

    from siac.domain.sensors import SensorBand, SensorConfig
    from siac.runtime import AtmosphericState, GeometryAngles

logger = logging.getLogger(__name__)


class SixSBackend:
    """Direct native 6SV2.1 backend with batched array execution."""

    def __init__(
        self,
        *,
        sixs_config: Any,
        sensor_config: SensorConfig | None = None,
        rt_setup: Any | None = None,
        runner: SixSNativeRunner | None = None,
    ) -> None:
        self._config = sixs_config
        self._sensor_config = sensor_config
        self._rt_setup = resolve_backend_rt_setup("sixs", rt_setup)
        self._runner = runner or SixSNativeRunner(
            sixs_config=sixs_config,
            sensor_config=sensor_config,
            rt_setup=self._rt_setup,
        )

    @property
    def backend_name(self) -> str:
        return "sixs"

    @property
    def requested_output_variables(self) -> tuple[str, ...]:
        return tuple(getattr(self._config, "output_variables", ("xap", "xbp", "xcp")))

    @property
    def rt_setup(self) -> Any | None:
        return self._rt_setup

    def supports_jacobian(self) -> bool:
        return False

    def is_available_for_sensor(self, sensor_id: str, satellite_id: str) -> bool:
        _ = (sensor_id, satellite_id)
        return True

    def set_observation_time(self, observation_time: datetime | None) -> None:
        self._runner.set_observation_time(observation_time)

    def compute_coefficients(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        band: SensorBand,
        compute_jacobian: bool = False,
    ) -> RTCoefficients:
        if compute_jacobian:
            raise NotImplementedError("The native 6S backend does not expose Jacobians yet.")
        outputs = self._runner.compute_coefficients(
            geometry=geometry,
            atmo_state=atmo_state,
            band=band,
            output_variables=self.requested_output_variables,
        )
        return self._coerce_result(outputs)

    def compute_coefficients_multi(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand],
        compute_jacobian: bool = False,
    ) -> list[RTCoefficients]:
        outputs = self._runner.compute_coefficients_multi(
            geometry=geometry,
            atmo_state=atmo_state,
            bands=bands,
            compute_jacobian=compute_jacobian,
            output_variables=self.requested_output_variables,
        )
        return [self._coerce_result(item) for item in outputs]

    def preload_scene_subset(self, *args: Any, **kwargs: Any) -> Any:
        return self._runner.preload_scene_subset(*args, **kwargs)

    def build_joint_grid_search_lut(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        aot_axis: np.ndarray,
        tcwv_axis: np.ndarray,
        bands: list[SensorBand],
    ) -> JointGridSearchLUT | None:
        """Build a joint (aot × tcwv × geometry) LUT for grid-search reuse.

        The block-grid-search calls the RT model with the same scene-level
        geometry/atmosphere but many different ``(aot, tcwv)`` candidate
        pairs. Computing a single LUT spanning the whole grid amortises the
        6S work across all candidates. Returns ``None`` when the
        optimization is disabled or unsupported by the underlying runner —
        the caller should then fall back to per-candidate
        :meth:`compute_coefficients` invocations.
        """
        return self._runner.build_joint_grid_search_lut(
            geometry=geometry,
            atmo_state=atmo_state,
            aot_axis=aot_axis,
            tcwv_axis=tcwv_axis,
            bands=bands,
            output_variables=self.requested_output_variables,
        )

    def prepare_correction_scene(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
    ) -> PreparedCorrectionScene:
        """Share the scene prep + LUT plan across per-band correction calls.

        Wave 18e: M6 correction calls compute_coefficients per band, which
        previously redid the same per-scene prep work 13 times. With this
        helper the prep runs once and is reused by
        :meth:`compute_coefficients_with_prepared` for each band.
        """
        return self._runner.prepare_correction_scene(
            geometry=geometry, atmo_state=atmo_state
        )

    def compute_coefficients_with_prepared(
        self,
        prepared_scene: PreparedCorrectionScene,
        band: SensorBand,
    ) -> RTCoefficients:
        """Run RT for a single band using shared scene prep (wave 18e)."""
        outputs = self._runner.compute_coefficients_with_prepared(
            prepared_scene=prepared_scene,
            band=band,
            output_variables=self.requested_output_variables,
        )
        return self._coerce_result(outputs)

    def preload_joint_grid_search_lut(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        aot_axis: np.ndarray,
        tcwv_axis: np.ndarray,
        bands: list[SensorBand],
    ) -> JointGridSearchLUT | None:
        """Build the joint LUT eagerly and cache it for the next solver call.

        Wave 18 (opt 4): called during prior-fetch so the joint LUT's 6S
        work overlaps with the surface-prior reprojection — saves the
        single-pass cost of the LUT build off the critical path. The next
        ``build_joint_grid_search_lut`` call with matching inputs picks
        up the cached LUT (single-shot — staged solvers that rebuild with
        different inputs miss and recompute fresh).
        """
        return self._runner.preload_joint_grid_search_lut(
            geometry=geometry,
            atmo_state=atmo_state,
            aot_axis=aot_axis,
            tcwv_axis=tcwv_axis,
            bands=bands,
            output_variables=self.requested_output_variables,
        )

    @staticmethod
    def _coerce_result(outputs: Any) -> RTCoefficients:
        if isinstance(outputs, RTCoefficients):
            return outputs
        if not isinstance(outputs, dict):
            raise TypeError(
                "SixS runner must return RTCoefficients or dict[str, xr.DataArray], "
                f"got {type(outputs).__name__}"
            )

        missing = {"xap", "xbp", "xcp"} - set(outputs)
        if missing:
            raise KeyError(
                "SixS runner did not return required coefficients: " + ", ".join(sorted(missing))
            )

        extras = {
            name: value for name, value in outputs.items() if name not in {"xap", "xbp", "xcp"}
        }
        return RTCoefficients(
            xap=outputs["xap"],
            xbp=outputs["xbp"],
            xcp=outputs["xcp"],
            extras=extras,
        )
