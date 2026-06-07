"""libRadtran (``uvspec``) RT backend.

A runtime :class:`~siac.domain.protocols.RTModelBackend` that drives the real
libRadtran engine. It mirrors :class:`~siac.algorithms.rt.direct.sixs.SixSBackend`:
a thin protocol-conforming adapter delegating to an injectable
:class:`~siac.algorithms.rt.direct.libradtran_runner.LibRadtranRunner`, which
builds a scene-scoped spectral grid via ``uvspec`` and derives the standard
``RTCoefficients`` through the reused LUT spectral path.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

from siac.algorithms.rt.direct.libradtran_runner import LibRadtranRunner
from siac.rt_setup import resolve_backend_rt_setup
from siac.runtime import RTCoefficients

if TYPE_CHECKING:
    from datetime import datetime

    from siac.domain.sensors import SensorBand, SensorConfig
    from siac.runtime import AtmosphericState, GeometryAngles

logger = logging.getLogger(__name__)


class LibRadtranBackend:
    """Direct libRadtran backend building scene-scoped grids via ``uvspec``."""

    def __init__(
        self,
        *,
        libradtran_config: Any,
        sensor_config: SensorConfig | None = None,
        rt_setup: Any | None = None,
        runner: LibRadtranRunner | None = None,
    ) -> None:
        self._config = libradtran_config
        self._sensor_config = sensor_config
        self._rt_setup = resolve_backend_rt_setup("libradtran", rt_setup)
        self._runner = runner or LibRadtranRunner(
            libradtran_config=libradtran_config,
            sensor_config=sensor_config,
            rt_setup=self._rt_setup,
        )

    @property
    def backend_name(self) -> str:
        return "libradtran"

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
            raise NotImplementedError("The libRadtran backend does not expose Jacobians yet.")
        outputs = self._runner.compute_coefficients(
            geometry=geometry,
            atmo_state=atmo_state,
            band=band,
        )
        return self._coerce_result(outputs)

    def compute_coefficients_multi(
        self,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        bands: list[SensorBand],
        compute_jacobian: bool = False,
    ) -> list[RTCoefficients]:
        # The runner caches the assembled scene LUT keyed on the scene
        # signature, so each band reuses one set of uvspec runs.
        if compute_jacobian:
            raise NotImplementedError("The libRadtran backend does not expose Jacobians yet.")
        return [self.compute_coefficients(geometry, atmo_state, band) for band in bands]

    @staticmethod
    def _coerce_result(outputs: Any) -> RTCoefficients:
        if isinstance(outputs, RTCoefficients):
            return outputs
        if not isinstance(outputs, dict):
            raise TypeError(
                "LibRadtran runner must return RTCoefficients or dict[str, xr.DataArray], "
                f"got {type(outputs).__name__}"
            )
        missing = {"xap", "xbp", "xcp"} - set(outputs)
        if missing:
            raise KeyError(
                "LibRadtran runner did not return required coefficients: "
                + ", ".join(sorted(missing))
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
