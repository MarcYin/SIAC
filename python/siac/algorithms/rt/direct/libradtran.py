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

from siac.algorithms.rt.direct._common import coerce_rt_coefficients
from siac.algorithms.rt.direct.libradtran_runner import LibRadtranRunner
from siac.rt_setup import resolve_backend_rt_setup

if TYPE_CHECKING:
    from datetime import datetime

    from siac.domain.sensors import SensorBand, SensorConfig
    from siac.runtime import AtmosphericState, GeometryAngles, RTCoefficients

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

    def build_joint_grid_search_lut(
        self,
        *,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        aot_axis: Any,
        tcwv_axis: Any,
        bands: list[SensorBand],
    ) -> Any:
        """Amortise the block-grid-search RT cost across (aot, tcwv) candidates.

        The solver probes the RT model with many candidate ``(aot, tcwv)`` pairs.
        Without this hook each call would rebuild the (multi-minute) ``uvspec``
        scene LUT. Here the runner builds one scene LUT spanning the grid-search
        range and returns a joint LUT that serves each candidate by
        interpolation - the same amortisation 6S uses. Mirrors ``SixSBackend``.
        """
        return self._runner.build_joint_grid_search_lut(
            geometry=geometry,
            atmo_state=atmo_state,
            aot_axis=aot_axis,
            tcwv_axis=tcwv_axis,
            bands=bands,
        )

    @staticmethod
    def _coerce_result(outputs: Any) -> RTCoefficients:
        return coerce_rt_coefficients(outputs, runner_name="LibRadtran")
