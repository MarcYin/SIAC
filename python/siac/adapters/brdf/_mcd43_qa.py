"""QA-to-uncertainty decoding for the Earthaccess BRDF providers.

Pure array logic moved out of ``siac.adapters.brdf.mcd43_earthaccess``:
mapping MCD43/VNP43/MCD19 per-band QA values to reflectance uncertainties.
The provider classes bind these functions as staticmethods so the existing
method seams (``MCD43EarthAccessProvider._qa_to_uncertainty`` etc.) keep
resolving. This module must not import ``mcd43_earthaccess`` (one-way).
"""

from __future__ import annotations

import numpy as np
import xarray as xr

_BEST_QA_REFLECTANCE_UNCERTAINTY = 0.015
_QA_UNCERTAINTY_POWER = 1.6


def qa_to_uncertainty(qa: xr.DataArray) -> xr.DataArray:
    qa_values = qa_values_to_uncertainty(qa.values)
    return xr.DataArray(qa_values, dims=qa.dims, coords=qa.coords)


def qa_values_to_uncertainty(qa_values: np.ndarray) -> np.ndarray:
    qa_values = np.asarray(qa_values, dtype=np.float32)
    unc = np.full(qa_values.shape, np.nan, dtype=np.float32)
    valid = np.isfinite(qa_values) & (qa_values >= 0.0)
    unc = np.where(
        valid,
        _BEST_QA_REFLECTANCE_UNCERTAINTY * np.power(qa_values + 1.0, _QA_UNCERTAINTY_POWER),
        unc,
    )
    unc_array: np.ndarray = np.asarray(unc, dtype=np.float32)
    return unc_array
