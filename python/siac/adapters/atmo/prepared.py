"""Aerosol prior read from a prepared per-scene store.

Deriving the aerosol prior live from satellite granules is fragile: over a small
AOI the granule may be missing, cloud-screened away, or QA-rejected entirely, and
the provider then substitutes a default. For a surface-driven retrieval that is
the worst possible failure, because the retrieval is prior-limited wherever the
visible bands constrain AOT weakly — it will simply return a value near the
fabricated prior, and the result looks plausible while being wrong.

Preparing the aerosol prior offline removes that failure mode: the extraction
runs once, per scene, where it can aggregate over a generous window and be
checked, and the retrieval then reads a known-good scalar. It is the aerosol
counterpart of :class:`siac.adapters.surface_library.PreparedSurfaceLibrary`.

Only aerosol optical depth is prepared. Water vapour, ozone and terrain come
from a base provider, since those are not the prior-limiting term.
"""

from __future__ import annotations

import json
import logging
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import xarray as xr

if TYPE_CHECKING:
    from datetime import datetime

    from siac.domain.protocols import AtmosphericPriorProvider
    from siac.runtime.models import AtmosphericState

logger = logging.getLogger(__name__)

__all__ = ["PreparedScalarAODProvider"]


class PreparedScalarAODProvider:
    """Aerosol prior whose AOD comes from a prepared per-scene scalar.

    The store holds one ``<scene_key>.json`` per scene with an ``aot`` value and
    optionally ``aot_unc``. ``base`` supplies the rest of the atmospheric state
    and the grid the scalar is broadcast onto.
    """

    def __init__(
        self,
        base: AtmosphericPriorProvider,
        root: str | Path,
        *,
        scene_key: str | None = None,
        required: bool = True,
        tcwv_cm: float | None = None,
    ) -> None:
        self._base = base
        self._root = Path(root).expanduser()
        self._scene_key = scene_key
        self._tcwv_cm = tcwv_cm
        self._required = required

    @property
    def source_name(self) -> str:
        return f"prepared_aod[{self._root.name}]"

    def _read(self) -> tuple[float, float | None] | None:
        if not self._scene_key:
            if self._required:
                raise ValueError(
                    "Prepared aerosol prior needs a scene key; set "
                    "providers.atmo.prepared_scalar_scene_key."
                )
            return None
        path = self._root / f"{self._scene_key}.json"
        if not path.is_file():
            if self._required:
                raise FileNotFoundError(
                    f"Prepared aerosol prior has no entry for scene "
                    f"{self._scene_key!r} under {self._root}."
                )
            logger.warning("No prepared aerosol prior for %s; using the base provider.", path)
            return None
        payload = json.loads(path.read_text())
        aot = payload.get("aot")
        if aot is None or not np.isfinite(float(aot)):
            if self._required:
                raise ValueError(f"Prepared aerosol prior {path} has no usable 'aot'.")
            return None
        unc = payload.get("aot_unc")
        return float(aot), (float(unc) if unc is not None and np.isfinite(float(unc)) else None)

    def get_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState:
        state = self._base.get_prior(bounds, crs, obs_time, resolution)
        prepared = self._read()
        if prepared is None:
            return state

        aot_value, unc_value = prepared
        aot = xr.full_like(state.aot, np.float32(aot_value))
        logger.info(
            "Aerosol prior for %s: prepared AOD %.3f (base %s reported %.3f).",
            self._scene_key,
            aot_value,
            getattr(self._base, "source_name", "base"),
            float(np.nanmedian(np.asarray(state.aot.values, dtype=float))),
        )
        updates: dict[str, Any] = {"aot": aot}
        if unc_value is not None:
            updates["aot_unc"] = xr.full_like(state.aot_unc, np.float32(unc_value))
        if self._tcwv_cm is not None:
            updates["tcwv"] = xr.full_like(state.tcwv, np.float32(self._tcwv_cm))
        return replace(state, **updates)
