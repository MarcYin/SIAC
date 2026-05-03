"""Radiative-transfer model assembly helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from siac.adapters.rt import build_rt_model

if TYPE_CHECKING:
    from siac.domain.protocols import RTModelBackend
    from siac.domain.sensors import SensorConfig


def resolve_rt_model_for_pipeline(
    config: Any,
    auth: Any = None,
    *,
    sensor_config: SensorConfig | None = None,
) -> RTModelBackend:
    return build_rt_model(config, auth=auth, sensor_config=sensor_config)
