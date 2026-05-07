"""Built-in sensor registry.

Provides the ``SENSOR_CONFIGS`` mapping plus a small extension API so external
packages can register additional sensors without monkey-patching the dict
directly. See REVIEW.md §3.7.

Usage from a downstream package::

    from siac.catalog.sensors import register, SensorConfig
    register("MYSENSOR", "SAT-1", my_config)

Subsequent calls to :func:`get_sensor_config` resolve the new entry.
"""

from __future__ import annotations

from siac.catalog.sensors.landsat import LANDSAT8_OLI_CONFIG, LANDSAT9_OLI2_CONFIG
from siac.catalog.sensors.sentinel2 import (
    SENTINEL2A_CONFIG,
    SENTINEL2B_CONFIG,
    SENTINEL2C_CONFIG,
)
from siac.domain.sensors import SensorConfig

SENSOR_CONFIGS: dict[tuple[str, str], SensorConfig] = {
    ("MSI", "S2A"): SENTINEL2A_CONFIG,
    ("MSI", "S2B"): SENTINEL2B_CONFIG,
    ("MSI", "S2C"): SENTINEL2C_CONFIG,
    ("OLI", "L8"): LANDSAT8_OLI_CONFIG,
    ("OLI", "L9"): LANDSAT9_OLI2_CONFIG,
}


def register(
    sensor_id: str,
    satellite_id: str,
    config: SensorConfig,
    *,
    overwrite: bool = False,
) -> None:
    """Register an external sensor configuration.

    Args:
        sensor_id: Instrument identifier (e.g. ``"MSI"``, ``"OLI"``).
        satellite_id: Platform identifier (e.g. ``"S2A"``, ``"L8"``).
        config: The :class:`SensorConfig` to register. Its ``sensor_id`` /
            ``satellite_id`` fields must match the keys.
        overwrite: When ``False`` (default), raises ``KeyError`` if the
            ``(sensor_id, satellite_id)`` pair is already registered. Set to
            ``True`` to replace an existing entry.

    Raises:
        ValueError: If ``config.sensor_id`` / ``config.satellite_id`` does not
            match the supplied keys.
        KeyError: If the entry already exists and ``overwrite`` is ``False``.
    """
    if config.sensor_id != sensor_id or config.satellite_id != satellite_id:
        raise ValueError(
            "register(): config.sensor_id/satellite_id "
            f"({config.sensor_id!r}/{config.satellite_id!r}) must match "
            f"keys ({sensor_id!r}/{satellite_id!r})"
        )
    key = (sensor_id, satellite_id)
    if key in SENSOR_CONFIGS and not overwrite:
        raise KeyError(
            f"Sensor {sensor_id}/{satellite_id} is already registered; "
            "pass overwrite=True to replace."
        )
    SENSOR_CONFIGS[key] = config


def get_sensor_config(sensor_id: str, satellite_id: str) -> SensorConfig:
    """Get sensor configuration by ID."""
    key = (sensor_id, satellite_id)
    if key not in SENSOR_CONFIGS:
        raise KeyError(
            f"Unknown sensor {sensor_id}/{satellite_id}. Available: {list(SENSOR_CONFIGS.keys())}"
        )
    return SENSOR_CONFIGS[key]
