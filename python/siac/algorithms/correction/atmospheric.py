"""Atmospheric correction: TOA to BOA conversion."""

from __future__ import annotations

import logging
import time
from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING, Any

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
    def __init__(
        self,
        rt_model: RTModelBackend,
        sensor_config: SensorConfig,
        correction_workers: int = 1,
    ):
        if not isinstance(rt_model, RTModelBackend):
            raise TypeError(
                f"rt_model must implement RTModelBackend protocol, got {type(rt_model).__name__}"
            )
        if correction_workers < 1:
            raise ValueError("correction_workers must be >= 1")
        self.rt_model = rt_model
        self.sensor_config = sensor_config
        self.correction_workers = int(correction_workers)

    def correct(
        self,
        toa: xr.Dataset,
        geometry: GeometryAngles,
        atmo_state: AtmosphericState,
        cloud_mask: xr.DataArray | None = None,
        boa_band_writer: Callable[[str, xr.DataArray], xr.DataArray] | None = None,
    ) -> CorrectionResult:
        t0 = time.monotonic()
        logger.info(
            "M6 correction starting: %d sensor bands, %d workers",
            len(self.sensor_config.band_names),
            self.correction_workers,
        )
        boa_vars = {}
        invalid_boa_mask: xr.DataArray | None = None
        band_loader = toa.attrs.get(_TOA_BAND_LOADER_ATTR)
        compute_time_by_band: dict[str, float] = {}
        write_time_by_band: dict[str, float] = {}

        work_items: list[tuple[str, Any, xr.DataArray]] = []
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

            work_items.append((band_name, band_spec, band_data))

        def _correct_single_band(
            band_spec: Any,
            band_data: xr.DataArray,
        ) -> tuple[xr.DataArray, xr.DataArray, float]:
            t_band = time.perf_counter()
            coeffs = self.rt_model.compute_coefficients(geometry, atmo_state, band_spec, False)
            coeffs = resample_coefficients_to_template(coeffs, band_data)
            boa = coeffs.apply_correction(band_data)
            band_valid = np.isfinite(boa) & (boa > -0.05) & (boa < 1.5)
            masked_boa = boa.where(band_valid)
            return masked_boa, (~band_valid), (time.perf_counter() - t_band)

        per_band_results: dict[str, tuple[xr.DataArray, xr.DataArray, float]] = {}
        t_compute_phase = time.perf_counter()
        if self.correction_workers > 1 and len(work_items) > 1:
            max_workers = min(self.correction_workers, len(work_items))
            # RT coefficient computation is CPU-bound; use threads (ProcessPoolExecutor
            # would require pickling xarray/numpy objects, which adds overhead).
            # ThreadPoolExecutor still helps when the RT model releases the GIL
            # (e.g. Rust-backed emulator) or when I/O is mixed in.
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    band_name: executor.submit(_correct_single_band, band_spec, band_data)
                    for band_name, band_spec, band_data in work_items
                }
                for band_name, future in futures.items():
                    per_band_results[band_name] = future.result()
        else:
            for band_name, band_spec, band_data in work_items:
                per_band_results[band_name] = _correct_single_band(band_spec, band_data)
        compute_phase_s = time.perf_counter() - t_compute_phase

        t_write_phase = time.perf_counter()
        for band_name, _band_spec, _band_data in work_items:
            masked_boa, band_invalid, compute_s = per_band_results[band_name]
            compute_time_by_band[band_name] = compute_s
            t_write_band = time.perf_counter()
            if boa_band_writer is not None:
                masked_boa = boa_band_writer(band_name, masked_boa)
            write_s = time.perf_counter() - t_write_band
            write_time_by_band[band_name] = write_s
            boa_vars[band_name] = masked_boa
            invalid_boa_mask = (
                band_invalid if invalid_boa_mask is None else (invalid_boa_mask | band_invalid)
            )

            logger.info(
                "M6 band %s timings: compute=%.3fs write=%.3fs total=%.3fs",
                band_name,
                compute_s,
                write_s,
                compute_s + write_s,
            )
        write_phase_s = time.perf_counter() - t_write_phase

        if work_items:
            logger.info(
                (
                    "M6 correction timing summary: bands=%d workers=%d "
                    "compute_phase=%.3fs write_phase=%.3fs "
                    "compute_band_sum=%.3fs write_band_sum=%.3fs"
                ),
                len(work_items),
                self.correction_workers,
                compute_phase_s,
                write_phase_s,
                float(sum(compute_time_by_band.values())),
                float(sum(write_time_by_band.values())),
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
            final_cloud_mask = resample_mask_to_template(
                cloud_mask, invalid_boa_mask
            ) | invalid_boa_mask.astype(bool)
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
