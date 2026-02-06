"""Atmospheric correction: TOA to BOA conversion."""
from __future__ import annotations
import numpy as np
import xarray as xr
from dataclasses import dataclass
from typing import Any
from siac.core.types import AtmosphericState, GeometryAngles, SensorConfig

@dataclass
class CorrectionResult:
    boa: xr.Dataset
    boa_unc: xr.Dataset | None
    aot: xr.DataArray
    tcwv: xr.DataArray
    mask: xr.DataArray

class AtmosphericCorrector:
    def __init__(self, rt_model: Any, sensor_config: SensorConfig):
        self.rt_model = rt_model
        self.sensor_config = sensor_config

    def correct(self, toa: xr.Dataset, geometry: GeometryAngles, atmo_state: AtmosphericState,
                cloud_mask: xr.DataArray | None = None) -> CorrectionResult:
        boa_vars = {}
        for band_name in toa.data_vars:
            try:
                band_spec = self.sensor_config.get_band(band_name)
            except KeyError:
                continue
            coeffs = self.rt_model.compute_coefficients(geometry, atmo_state, band_spec, False)
            y = coeffs.xap * toa[band_name] - coeffs.xbp
            boa = y / (1.0 + coeffs.xcp * y)
            boa_vars[band_name] = boa.where((boa > 0) & (boa < 1.5))

        boa_ds = xr.Dataset(boa_vars)
        first = list(boa_vars.keys())[0]
        mask = np.isfinite(boa_ds[first]) & (boa_ds[first] > 0)
        if cloud_mask is not None:
            mask = mask & ~cloud_mask
        return CorrectionResult(boa_ds, None, atmo_state.aot, atmo_state.tcwv, mask)
