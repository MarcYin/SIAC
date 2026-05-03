"""Radiative-transfer backend adapters."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

from siac.rt_setup import resolve_effective_rt_setup

if TYPE_CHECKING:
    from siac.adapters.auth import CredentialManager
    from siac.domain.sensors import SensorConfig

logger = logging.getLogger(__name__)
TwoLayerNNEmulator: Any | None = None
ZarrLUTBackend: Any | None = None
SixSBackend: Any | None = None

_SENSOR_DEFAULTS: dict[str, tuple[str, str]] = {
    "s2": ("MSI", "S2A"),
    "s2a": ("MSI", "S2A"),
    "s2b": ("MSI", "S2B"),
    "s2c": ("MSI", "S2C"),
    "sentinel2": ("MSI", "S2A"),
    "l8": ("OLI", "L8"),
    "l9": ("OLI", "L9"),
    "auto": ("MSI", "S2A"),
}


def build_rt_model(
    config: Any,
    auth: CredentialManager | None = None,
    *,
    sensor_config: SensorConfig | None = None,
) -> Any:
    """Build the RT backend for the resolved config."""
    rt_config = config.algorithms.rt
    if sensor_config is not None:
        sensor_id, satellite_id = sensor_config.sensor_id, sensor_config.satellite_id
    else:
        sensor_name = str(getattr(config, "sensor", "auto") or "auto").lower()
        sensor_id, satellite_id = _SENSOR_DEFAULTS.get(sensor_name, ("MSI", "S2A"))

    paths = getattr(config, "paths", None)
    emulator_dir = getattr(rt_config, "emulator_dir", None) or getattr(paths, "emulator_dir", None)
    lut_path = getattr(rt_config, "lut_path", None) or getattr(paths, "lut_path", None)
    backend = rt_config.backend

    if backend == "emulator":
        try:
            emulator_cls = TwoLayerNNEmulator
            if emulator_cls is None:
                from siac.algorithms.rt.emulator import TwoLayerNNEmulator as emulator_cls

            return emulator_cls(
                emulator_dir=emulator_dir,
                sensor_id=sensor_id,
                satellite_id=satellite_id,
            )
        except Exception as exc:
            if not rt_config.fallback_to_lut or not lut_path:
                raise ValueError(
                    "Cannot resolve emulator RT model and LUT fallback is unavailable."
                ) from exc
            logger.warning("Emulator RT model unavailable; falling back to LUT backend.")
            backend = "lut"

    if backend == "lut" and lut_path:
        lut_backend_cls = ZarrLUTBackend
        if lut_backend_cls is None:
            from siac.algorithms.rt.lut import ZarrLUTBackend as lut_backend_cls

        storage_options = dict(rt_config.lut_storage_options)
        if (
            auth is not None
            and auth.aws().has_credentials()
            and "key" not in storage_options
            and str(lut_path).startswith("s3://")
        ):
            storage_options.update(auth.aws().storage_options())
        return lut_backend_cls(
            lut_path,
            interpolation_method=rt_config.lut_interpolation,
            storage_options=storage_options,
            rt_setup=resolve_effective_rt_setup(rt_config, "lut"),
        )

    if backend == "sixs":
        sixs_backend_cls = SixSBackend
        if sixs_backend_cls is None:
            from siac.algorithms.rt.direct.sixs import SixSBackend as sixs_backend_cls

        return sixs_backend_cls(
            sixs_config=rt_config.sixs,
            sensor_config=sensor_config,
            rt_setup=resolve_effective_rt_setup(rt_config, "sixs"),
        )

    raise ValueError(f"Cannot resolve RT model from config: backend={rt_config.backend!r}")
