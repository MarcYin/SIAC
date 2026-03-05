"""Atmospheric correction: TOA to BOA conversion."""
from __future__ import annotations

import time
from typing import Any

import numpy as np
import xarray as xr

from siac.core.protocols import RTModelBackend
from siac.core.types import (
    AtmosphericState,
    CorrectionResult,
    GeometryAngles,
    SensorConfig,
)


class AtmosphericCorrector:
    def __init__(self, rt_model: Any, sensor_config: SensorConfig):
        if not isinstance(rt_model, RTModelBackend):
            raise TypeError(
                f"rt_model must implement RTModelBackend protocol, "
                f"got {type(rt_model).__name__}"
            )
        self.rt_model = rt_model
        self.sensor_config = sensor_config

    def correct(self, toa: xr.Dataset, geometry: GeometryAngles, atmo_state: AtmosphericState,
                cloud_mask: xr.DataArray | None = None) -> CorrectionResult:
        t0 = time.monotonic()
        boa_vars = {}
        for band_name in toa.data_vars:
            try:
                band_spec = self.sensor_config.get_band(band_name)
            except KeyError:
                continue
            coeffs = self.rt_model.compute_coefficients(geometry, atmo_state, band_spec, False)
            boa = coeffs.apply_correction(toa[band_name])
            boa_vars[band_name] = boa.where((boa > 0) & (boa < 1.5))

        boa_ds = xr.Dataset(boa_vars)
        first = list(boa_vars.keys())[0]
        mask = np.isfinite(boa_ds[first]) & (boa_ds[first] > 0)
        if cloud_mask is not None:
            mask = mask & ~cloud_mask
        elapsed = time.monotonic() - t0
        return CorrectionResult(
            boa=boa_ds,
            boa_unc=None,
            aot=atmo_state.aot,
            tcwv=atmo_state.tcwv,
            cloud_mask=mask,
            metadata={"processing_time_s": elapsed},
        )
