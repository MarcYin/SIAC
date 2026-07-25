"""Fused aerosol-prior provider combining several AOD sources.

The surface-driven retrieval defers to the aerosol prior wherever the surface
signal is weak, which is most of the clean and moderate AOD range: the visible
bands constrain AOT only weakly there, so the retrieval inherits the prior's
regional bias. That makes the *centre* of the aerosol prior the dominant error
term, and it is measurable — replacing it with the truth lifts within-EE
retrieval from 79% to 87% on the 152-matchup campaign.

Both routinely available AOD sources under-estimate, and they do so
independently: a satellite retrieval (MAIAC/MCD19) reads low by about 0.05 on
average and 0.17 above AOD 0.5, and an assimilating model (CAMS) reads low by
about 0.08 and 0.27. Because the two errors share a direction but not a
mechanism, taking their maximum removes most of the common bias by
construction, with no fitting against reference data: as a point predictor
``max`` scores 73% within EE at a bias of −0.004, against 70% / −0.047 for
MAIAC alone. Used as the retrieval's aerosol prior it lifts within-EE from
79.2% to 84.6% (RMSE 0.120 → 0.110, bias −0.021 → +0.001), recovering roughly
70% of the headroom to a perfect prior, with the gain concentrated exactly
where the retrieval is prior-limited (the moderate band, 69% → 84%).

The fusion applies to the aerosol optical depth only. Water vapour, ozone and
terrain come from the primary provider, so the sources being fused need supply
nothing but AOD.
"""

from __future__ import annotations

import logging
from dataclasses import replace
from typing import TYPE_CHECKING, Literal

import numpy as np
import xarray as xr

if TYPE_CHECKING:
    from collections.abc import Sequence
    from datetime import datetime

    from siac.domain.protocols import AtmosphericPriorProvider
    from siac.runtime.models import AtmosphericState

logger = logging.getLogger(__name__)

__all__ = ["AODFusionOp", "FusedAODProvider"]

AODFusionOp = Literal["max", "min", "mean"]


def _combine(values: Sequence[xr.DataArray], op: AODFusionOp) -> xr.DataArray:
    """Combine aligned AOD fields, ignoring sources that are missing a pixel."""

    stacked = xr.concat(values, dim="_aod_source")
    if op == "max":
        combined = stacked.max(dim="_aod_source", skipna=True)
    elif op == "min":
        combined = stacked.min(dim="_aod_source", skipna=True)
    else:
        combined = stacked.mean(dim="_aod_source", skipna=True)
    # Where every source is missing, ``skipna`` yields NaN; keep that rather
    # than inventing a value, so downstream masking still sees a gap.
    return combined


def _as_reference_grid(reference: xr.DataArray, values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(
        np.asarray(values, dtype=np.float32).reshape(reference.shape),
        dims=reference.dims,
        coords=reference.coords,
    )


def _align_to(reference: xr.DataArray, other: xr.DataArray) -> xr.DataArray | None:
    """Put ``other`` on ``reference``'s grid, or return ``None`` if it cannot be.

    Aerosol providers return their own native grids — a satellite retrieval on a
    projected grid, a model on a geographic one — so sources generally differ in
    resolution, CRS and dimension names. Concatenating mismatched grids directly
    would broadcast them into a spurious extra dimension rather than align them.

    Returning ``None`` rather than raising keeps an unalignable source from
    failing the scene: the retrieval falls back to the primary prior, which is
    the behaviour without fusion at all.
    """

    if other.dims == reference.dims and other.shape == reference.shape:
        return _as_reference_grid(reference, other.values)

    # Geographic reprojection is the correct alignment when both sides carry a
    # CRS, which is the normal case for real providers.
    try:
        from siac.geo.reprojection import reproject_match

        aligned = reproject_match(other.astype(np.float32), reference)
        if aligned.size == reference.size:
            return _as_reference_grid(reference, aligned.values)
    except Exception as exc:  # noqa: BLE001 - try coordinate interpolation next
        logger.debug("AOD fusion source reprojection failed (%s); interpolating.", exc)

    # Without a CRS, interpolate by coordinate value, matching axes positionally
    # so differing dimension names (latitude/longitude vs y/x) still align.
    try:
        renamed = other.astype(np.float32).rename(
            dict(zip(other.dims, reference.dims, strict=True))
        )
        interpolated = renamed.interp(
            coords={dim: reference.coords[dim] for dim in reference.dims},
            method="linear",
            kwargs={"fill_value": None},
        )
        if interpolated.size == reference.size:
            return _as_reference_grid(reference, interpolated.values)
    except Exception as exc:  # noqa: BLE001 - source is unusable; skip it
        logger.debug("AOD fusion source interpolation failed (%s).", exc)

    return None


class FusedAODProvider:
    """Atmospheric prior whose AOD is the combination of several sources.

    ``primary`` supplies the full atmospheric state (and the grid every other
    source is aligned to); ``sources`` supply additional AOD fields that are
    combined with it under ``op``.
    """

    def __init__(
        self,
        primary: AtmosphericPriorProvider,
        sources: Sequence[AtmosphericPriorProvider],
        op: AODFusionOp = "max",
    ) -> None:
        if not sources:
            raise ValueError(
                "FusedAODProvider needs at least one additional AOD source to fuse "
                "with the primary provider."
            )
        self._primary = primary
        self._sources = tuple(sources)
        self._op: AODFusionOp = op

    @property
    def source_name(self) -> str:
        names = [getattr(self._primary, "source_name", "unknown")]
        names += [getattr(source, "source_name", "unknown") for source in self._sources]
        return f"{self._op}({', '.join(names)})"

    def get_prior(
        self,
        bounds: tuple[float, float, float, float],
        crs: str,
        obs_time: datetime,
        resolution: float,
    ) -> AtmosphericState:
        state = self._primary.get_prior(bounds, crs, obs_time, resolution)

        fields = [state.aot]
        for source in self._sources:
            other = source.get_prior(bounds, crs, obs_time, resolution)
            aligned = _align_to(state.aot, other.aot)
            if aligned is None:
                logger.warning(
                    "Could not align AOD source %s onto the %s grid; "
                    "continuing without it.",
                    getattr(source, "source_name", "unknown"),
                    getattr(self._primary, "source_name", "primary"),
                )
                continue
            fields.append(aligned)

        if len(fields) == 1:
            return state

        fused = _combine(fields, self._op)
        fused = fused.rename(state.aot.name) if state.aot.name else fused
        return replace(state, aot=fused)
