"""Radiative-transfer backend adapters.

Backend selection is registry-driven: each ``[algorithms.rt] backend`` value
maps to a builder in :data:`_RT_BACKEND_BUILDERS`. Adding a backend means
registering one builder (see :func:`register_rt_backend`) instead of editing
an if/elif chain — the builder receives a :class:`RTBuildContext` with the
resolved config, sensor identity, and credential manager.

The module-level ``TwoLayerNNEmulator`` / ``ZarrLUTBackend`` / ``SixSBackend``
/ ``LibRadtranBackend`` globals are deliberate test seams: they default to
``None`` (builders lazy-import the real class) and tests monkeypatch them to
inject fakes without touching the heavy backend modules.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from siac.rt_setup import resolve_effective_rt_setup

if TYPE_CHECKING:
    from collections.abc import Callable

    from siac.adapters.auth import CredentialManager
    from siac.domain.sensors import SensorConfig

logger = logging.getLogger(__name__)
TwoLayerNNEmulator: Any | None = None
ZarrLUTBackend: Any | None = None
SixSBackend: Any | None = None
LibRadtranBackend: Any | None = None

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


@dataclass(frozen=True)
class RTBuildContext:
    """Resolved inputs a backend builder needs (config + identity + auth)."""

    rt_config: Any
    sensor_config: SensorConfig | None
    sensor_id: str
    satellite_id: str
    emulator_dir: Any
    lut_path: Any
    auth: CredentialManager | None

    def unresolvable(self) -> ValueError:
        return ValueError(
            f"Cannot resolve RT model from config: backend={self.rt_config.backend!r}"
        )


def _build_emulator(ctx: RTBuildContext) -> Any:
    try:
        emulator_cls = TwoLayerNNEmulator
        if emulator_cls is None:
            from siac.algorithms.rt.emulator import TwoLayerNNEmulator as emulator_cls

        return emulator_cls(
            emulator_dir=ctx.emulator_dir,
            sensor_id=ctx.sensor_id,
            satellite_id=ctx.satellite_id,
        )
    except Exception as exc:
        if not ctx.rt_config.fallback_to_lut or not ctx.lut_path:
            raise ValueError(
                "Cannot resolve emulator RT model and LUT fallback is unavailable."
            ) from exc
        # exc_info=True so the operator can see why the emulator failed,
        # not just that we silently degraded (REVIEW.md §2.1, §3.3 rt.py).
        logger.warning(
            "Emulator RT model unavailable; falling back to LUT backend.",
            exc_info=True,
        )
        return _build_lut(ctx)


def _build_lut(ctx: RTBuildContext) -> Any:
    if not ctx.lut_path:
        raise ctx.unresolvable()
    lut_backend_cls = ZarrLUTBackend
    if lut_backend_cls is None:
        from siac.algorithms.rt.lut import ZarrLUTBackend as lut_backend_cls

    rt_config = ctx.rt_config
    storage_options = dict(rt_config.lut_storage_options)
    if (
        ctx.auth is not None
        and ctx.auth.aws().has_credentials()
        and "key" not in storage_options
        and str(ctx.lut_path).startswith("s3://")
    ):
        storage_options.update(ctx.auth.aws().storage_options())
    return lut_backend_cls(
        ctx.lut_path,
        interpolation_method=rt_config.lut_interpolation,
        storage_options=storage_options,
        rt_setup=resolve_effective_rt_setup(rt_config, "lut"),
    )


def _build_sixs(ctx: RTBuildContext) -> Any:
    sixs_backend_cls = SixSBackend
    if sixs_backend_cls is None:
        from siac.algorithms.rt.direct.sixs import SixSBackend as sixs_backend_cls

    return sixs_backend_cls(
        sixs_config=ctx.rt_config.sixs,
        sensor_config=ctx.sensor_config,
        rt_setup=resolve_effective_rt_setup(ctx.rt_config, "sixs"),
    )


def _build_libradtran(ctx: RTBuildContext) -> Any:
    libradtran_backend_cls = LibRadtranBackend
    if libradtran_backend_cls is None:
        from siac.algorithms.rt.direct.libradtran import (
            LibRadtranBackend as libradtran_backend_cls,
        )

    # Construction is cheap and engine-free (libRadtran is fetched/compiled
    # lazily on the first compute_coefficients call), so a missing uvspec
    # surfaces as a clear runtime error rather than here.
    return libradtran_backend_cls(
        libradtran_config=ctx.rt_config.libradtran,
        sensor_config=ctx.sensor_config,
        rt_setup=resolve_effective_rt_setup(ctx.rt_config, "libradtran"),
    )


_RT_BACKEND_BUILDERS: dict[str, Callable[[RTBuildContext], Any]] = {
    "emulator": _build_emulator,
    "lut": _build_lut,
    "sixs": _build_sixs,
    "libradtran": _build_libradtran,
}


def register_rt_backend(name: str, builder: Callable[[RTBuildContext], Any]) -> None:
    """Register (or replace) a backend builder for ``[algorithms.rt] backend``."""
    _RT_BACKEND_BUILDERS[str(name)] = builder


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
    ctx = RTBuildContext(
        rt_config=rt_config,
        sensor_config=sensor_config,
        sensor_id=sensor_id,
        satellite_id=satellite_id,
        emulator_dir=(
            getattr(rt_config, "emulator_dir", None) or getattr(paths, "emulator_dir", None)
        ),
        lut_path=getattr(rt_config, "lut_path", None) or getattr(paths, "lut_path", None),
        auth=auth,
    )
    builder = _RT_BACKEND_BUILDERS.get(str(rt_config.backend))
    if builder is None:
        raise ctx.unresolvable()
    return builder(ctx)
