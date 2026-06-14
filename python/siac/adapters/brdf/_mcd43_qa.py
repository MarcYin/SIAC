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

# Best-QA (QA==0) reflectance uncertainty when no reflectance is available to
# scale against — an *absolute* fallback used by the legacy earthaccess path.
_BEST_QA_REFLECTANCE_UNCERTAINTY = 0.015
_QA_UNCERTAINTY_POWER = 1.6

# Reflectance-relative best-QA uncertainty, used when a per-pixel reflectance is
# supplied: a fixed absolute floor (dark targets) combined in quadrature with a
# relative fraction (bright targets). The flat 0.015 absolute value above is
# ~50% relative over dark land (~0.03 reflectance), which makes the surface
# prior far too loose to constrain the aerosol solve, so the solver defaults to
# the (sometimes biased) CAMS AOT prior. Scaling with reflectance brings
# dark-target uncertainty down to ~0.003 — matching the Route-B
# monthly_database neighbour spread — where dark-target AOT retrieval needs it,
# while keeping bright surfaces near the legacy ~0.015.
_QA_REFLECTANCE_UNCERTAINTY_FLOOR = 0.003
_QA_REFLECTANCE_UNCERTAINTY_FRACTION = 0.05


def qa_to_uncertainty(
    qa: xr.DataArray,
    reflectance: xr.DataArray | np.ndarray | None = None,
) -> xr.DataArray:
    refl_values = reflectance.values if isinstance(reflectance, xr.DataArray) else reflectance
    qa_values = qa_values_to_uncertainty(qa.values, reflectance=refl_values)
    return xr.DataArray(qa_values, dims=qa.dims, coords=qa.coords)


def qa_values_to_uncertainty(
    qa_values: np.ndarray,
    reflectance: np.ndarray | None = None,
) -> np.ndarray:
    qa_values = np.asarray(qa_values, dtype=np.float32)
    unc = np.full(qa_values.shape, np.nan, dtype=np.float32)
    valid = np.isfinite(qa_values) & (qa_values >= 0.0)
    if reflectance is None:
        best_qa_uncertainty: np.ndarray | float = _BEST_QA_REFLECTANCE_UNCERTAINTY
    else:
        refl = np.asarray(reflectance, dtype=np.float32)
        relative = np.sqrt(
            _QA_REFLECTANCE_UNCERTAINTY_FLOOR**2
            + (_QA_REFLECTANCE_UNCERTAINTY_FRACTION * np.abs(refl)) ** 2
        )
        # Where reflectance is missing, fall back to the absolute best-QA value.
        best_qa_uncertainty = np.where(
            np.isfinite(refl), relative, np.float32(_BEST_QA_REFLECTANCE_UNCERTAINTY)
        )
    unc = np.where(
        valid,
        best_qa_uncertainty * np.power(qa_values + 1.0, _QA_UNCERTAINTY_POWER),
        unc,
    )
    unc_array: np.ndarray = np.asarray(unc, dtype=np.float32)
    return unc_array
