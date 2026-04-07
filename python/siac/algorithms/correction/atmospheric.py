"""Atmospheric correction: TOA to BOA conversion."""
from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable

import numpy as np
import xarray as xr

from siac.domain.protocols import RTModelBackend
from siac.geo.resample import (
    resample_coefficients_to_template,
    resample_mask_to_template,
)
from siac.runtime import (
    AtmosphericState,
    CorrectionDiagnostics,
    CorrectionResult,
    GeometryAngles,
)

if TYPE_CHECKING:
    from siac.domain import SensorConfig

logger = logging.getLogger(__name__)

_TOA_BAND_LOADER_ATTR = "_siac_toa_band_loader"


class AtmosphericCorrector:
    def __init__(self, rt_model: RTModelBackend, sensor_config: SensorConfig):
        if not isinstance(rt_model, RTModelBackend):
            raise TypeError(
                f"rt_model must implement RTModelBackend protocol, "
                f"got {type(rt_model).__name__}"
            )
        self.rt_model = rt_model
        self.sensor_config = sensor_config

    def correct(self, toa: xr.Dataset, geometry: GeometryAngles, atmo_state: AtmosphericState,
                cloud_mask: xr.DataArray | None = None,
                boa_band_writer: Callable[[str, xr.DataArray], xr.DataArray] | None = None) -> CorrectionResult:
        t0 = time.monotonic()
        boa_vars = {}
        invalid_boa_mask: xr.DataArray | None = None
        band_loader = toa.attrs.get(_TOA_BAND_LOADER_ATTR)
        for band_name in self.sensor_config.band_names:
            try:
                band_spec = self.sensor_config.get_band(band_name)
            except KeyError:
                continue
            band_data = toa.data_vars.get(band_name)
            if band_data is None and callable(band_loader):
                try:
                    band_data = band_loader(band_name)
                except KeyError:
                    band_data = None
            if band_data is None:
                continue
            coeffs = self.rt_model.compute_coefficients(geometry, atmo_state, band_spec, False)
            coeffs = resample_coefficients_to_template(coeffs, band_data)
            boa = coeffs.apply_correction(band_data)
            band_valid = np.isfinite(boa) & (boa > -0.05) & (boa < 1.5)
            masked_boa = boa.where(band_valid)
            if boa_band_writer is not None:
                masked_boa = boa_band_writer(band_name, masked_boa)
            boa_vars[band_name] = masked_boa
            band_invalid = ~band_valid
            invalid_boa_mask = (
                band_invalid
                if invalid_boa_mask is None
                else (invalid_boa_mask | band_invalid)
            )

        if not boa_vars:
            raise ValueError(
                "No bands in TOA dataset matched sensor_config bands. "
                f"TOA vars: {list(toa.data_vars)}, "
                f"sensor bands: {self.sensor_config.band_names}"
            )

        boa_ds = xr.Dataset(boa_vars)
        first = list(boa_vars.keys())[0]
        if invalid_boa_mask is None:
            invalid_boa_mask = xr.zeros_like(boa_ds[first], dtype=bool)

        # cloud_mask contract: True = cloudy/invalid, preserved from M1
        if cloud_mask is not None:
            final_cloud_mask = resample_mask_to_template(cloud_mask, invalid_boa_mask) | invalid_boa_mask.astype(bool)
        else:
            final_cloud_mask = invalid_boa_mask.astype(bool)
        elapsed = time.monotonic() - t0
        return CorrectionResult(
            boa=boa_ds,
            boa_unc=None,
            aot=atmo_state.aot,
            tcwv=atmo_state.tcwv,
            cloud_mask=final_cloud_mask,
            diagnostics=CorrectionDiagnostics(processing_time_s=elapsed),
        )
