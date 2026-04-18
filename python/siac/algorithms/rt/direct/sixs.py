"""Native 6SV2.1 RT backend."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

from siac.algorithms.rt.direct.sixs_native import SixSNativeRunner
from siac.runtime import RTCoefficients

if TYPE_CHECKING:
    from datetime import datetime

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
        runner: SixSNativeRunner | None = None,
    ) -> None:
        self._config = sixs_config
        self._sensor_config = sensor_config
        self._runner = runner or SixSNativeRunner(
            sixs_config=sixs_config,
            sensor_config=sensor_config,
        )

    @property
    def backend_name(self) -> str:
        return "sixs"

    @property
    def requested_output_variables(self) -> tuple[str, ...]:
        return tuple(getattr(self._config, "output_variables", ("xap", "xbp", "xcp")))

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
                "SixS runner did not return required coefficients: "
                + ", ".join(sorted(missing))
            )

        extras = {name: value for name, value in outputs.items() if name not in {"xap", "xbp", "xcp"}}
        return RTCoefficients(
            xap=outputs["xap"],
            xbp=outputs["xbp"],
            xcp=outputs["xcp"],
            extras=extras,
        )


__all__ = ["SixSBackend"]
