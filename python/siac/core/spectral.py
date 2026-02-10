"""
Sensor-agnostic spectral model for SIAC.

Provides:
- ``SpectralBandDescriptor``: band definition with Gaussian or full SRF
- ``sensor_to_reference()``: convolve sensor reflectance to a reference basis
- ``reference_to_sensor()``: project reference-basis reflectance back to sensor bands
- ``load_reference_rsr()``: load tabulated spectral response functions

See PLANS.md §9 for design rationale.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import xarray as xr
from numpy.typing import NDArray


def _trapezoid(y: NDArray, x: NDArray) -> float:
    """Compatibility wrapper for NumPy 1.x/2.x integration API."""
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


# ── Reference sensor definitions ──────────────────────────────────────

# MODIS land bands: (band_name, center_nm, fwhm_nm)
_MODIS_LAND_BANDS = (
    ("Band1", 645.0, 50.0),   # Red
    ("Band2", 858.5, 35.0),   # NIR
    ("Band3", 469.0, 20.0),   # Blue
    ("Band4", 555.0, 20.0),   # Green
    ("Band5", 1240.0, 20.0),  # SWIR-0
    ("Band6", 1640.0, 24.0),  # SWIR-1
    ("Band7", 2130.0, 50.0),  # SWIR-2
)


# ── SpectralBandDescriptor ────────────────────────────────────────────

@dataclass(frozen=True)
class SpectralBandDescriptor:
    """Description of a single spectral band.

    Supports two modes:
    1. **Gaussian approximation** (multispectral): center wavelength + FWHM
    2. **Full SRF** (hyperspectral / precision): tabulated (wavelength, response)
    """

    name: str
    center_wavelength_nm: float
    fwhm_nm: float
    resolution_m: float
    srf_wavelengths_nm: NDArray[np.floating] | None = None
    srf_response: NDArray[np.floating] | None = None

    @property
    def has_srf(self) -> bool:
        """Whether a tabulated SRF is available."""
        return self.srf_wavelengths_nm is not None and self.srf_response is not None

    @property
    def wavelength_um(self) -> float:
        """Center wavelength in micrometers."""
        return self.center_wavelength_nm / 1000.0

    def gaussian_response(self, wavelengths_nm: NDArray) -> NDArray:
        """Compute Gaussian spectral response at given wavelengths.

        Uses the stored center wavelength and FWHM.
        """
        sigma = self.fwhm_nm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        return np.exp(-0.5 * ((wavelengths_nm - self.center_wavelength_nm) / sigma) ** 2)

    def effective_response(self, wavelengths_nm: NDArray) -> NDArray:
        """Return the best available response function at given wavelengths.

        If a tabulated SRF exists, interpolate it; otherwise use Gaussian.
        """
        if self.has_srf:
            return np.interp(
                wavelengths_nm,
                self.srf_wavelengths_nm,
                self.srf_response,
                left=0.0,
                right=0.0,
            )
        return self.gaussian_response(wavelengths_nm)


# ── Reference RSR loading ─────────────────────────────────────────────

def load_reference_rsr(
    reference_sensor: str = "MODIS",
) -> dict[str, tuple[NDArray, NDArray]]:
    """Load tabulated spectral response functions for a reference sensor.

    For MODIS the RSR is approximated as a Gaussian from the known center
    wavelength and FWHM (sufficient for the band-mapping use-case).

    Returns:
        dict mapping band name -> (wavelengths_nm, response) arrays.

    Raises:
        ValueError: If the reference sensor is unknown.
    """
    if reference_sensor.upper() == "MODIS":
        result: dict[str, tuple[NDArray, NDArray]] = {}
        for band_name, center_nm, fwhm_nm in _MODIS_LAND_BANDS:
            sigma = fwhm_nm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
            wl = np.linspace(
                center_nm - 3.0 * fwhm_nm,
                center_nm + 3.0 * fwhm_nm,
                101,
            )
            resp = np.exp(-0.5 * ((wl - center_nm) / sigma) ** 2)
            result[band_name] = (wl, resp)
        return result

    raise ValueError(
        f"Unknown reference sensor: {reference_sensor!r}. Available: 'MODIS'"
    )


# ── Spectral convolution functions ────────────────────────────────────

def _build_conversion_matrix(
    sensor_bands: Sequence[SpectralBandDescriptor],
    ref_bands: dict[str, tuple[NDArray, NDArray]],
    common_wl: NDArray | None = None,
) -> tuple[NDArray, list[str]]:
    """Build a (n_ref, n_sensor) spectral conversion matrix.

    Each entry ``M[i, j]`` is the overlap integral between reference band *i*
    and sensor band *j*, normalised so that columns sum to 1 for a flat
    spectrum.  When the overlap is negligible (< 1e-6) the entry is zero.

    Returns:
        (matrix, ref_band_names)
    """
    if common_wl is None:
        common_wl = np.linspace(350.0, 2500.0, 2151)

    ref_names = list(ref_bands.keys())
    n_ref = len(ref_names)
    n_sensor = len(sensor_bands)

    # Pre-compute normalised reference responses on the common grid
    ref_resp = np.zeros((n_ref, len(common_wl)))
    for i, name in enumerate(ref_names):
        wl_r, resp_r = ref_bands[name]
        interp_r = np.interp(common_wl, wl_r, resp_r, left=0.0, right=0.0)
        norm = _trapezoid(interp_r, common_wl)
        ref_resp[i] = interp_r / norm if norm > 0 else interp_r

    # Pre-compute normalised sensor responses
    sensor_resp = np.zeros((n_sensor, len(common_wl)))
    for j, sb in enumerate(sensor_bands):
        sr = sb.effective_response(common_wl)
        norm = _trapezoid(sr, common_wl)
        sensor_resp[j] = sr / norm if norm > 0 else sr

    # Compute overlap matrix: integral of ref_i * sensor_j over wavelength
    matrix = np.zeros((n_ref, n_sensor))
    for i in range(n_ref):
        for j in range(n_sensor):
            overlap = _trapezoid(ref_resp[i] * sensor_resp[j], common_wl)
            matrix[i, j] = overlap if overlap > 1e-6 else 0.0

    return matrix, ref_names


def sensor_to_reference(
    sensor_reflectance: xr.Dataset,
    sensor_bands: Sequence[SpectralBandDescriptor],
    reference_sensor: str = "MODIS",
) -> NDArray:
    """Convolve sensor-band reflectance to a reference-sensor basis.

    Args:
        sensor_reflectance: Dataset with one variable per sensor band.
        sensor_bands: SpectralBandDescriptors matching the Dataset variables.
        reference_sensor: Name of the reference sensor (default "MODIS").

    Returns:
        Array of shape ``(n_ref_bands, *spatial)`` with reference-basis
        reflectance values.
    """
    ref_rsr = load_reference_rsr(reference_sensor)
    matrix, ref_names = _build_conversion_matrix(sensor_bands, ref_rsr)

    # Stack sensor reflectance into (n_sensor, *spatial)
    band_names = [b.name for b in sensor_bands]
    sensor_vals = np.stack(
        [sensor_reflectance[bn].values for bn in band_names if bn in sensor_reflectance],
        axis=0,
    )

    # matrix: (n_ref, n_sensor)  ×  sensor_vals: (n_sensor, *spatial) -> (n_ref, *spatial)
    spatial_shape = sensor_vals.shape[1:]
    flat = sensor_vals.reshape(sensor_vals.shape[0], -1)  # (n_sensor, N)
    ref_flat = matrix @ flat  # (n_ref, N)
    return ref_flat.reshape(len(ref_names), *spatial_shape)


def reference_to_sensor(
    ref_reflectance: NDArray,
    target_bands: Sequence[SpectralBandDescriptor],
    reference_sensor: str = "MODIS",
) -> NDArray:
    """Project reference-basis reflectance back to sensor bands.

    Uses the pseudo-inverse of the conversion matrix for a least-squares
    mapping.

    Args:
        ref_reflectance: Array of shape ``(n_ref_bands, *spatial)``.
        target_bands: SpectralBandDescriptors of the target sensor.
        reference_sensor: Name of the reference sensor.

    Returns:
        Array of shape ``(n_target_bands, *spatial)``.
    """
    ref_rsr = load_reference_rsr(reference_sensor)
    matrix, _ = _build_conversion_matrix(target_bands, ref_rsr)

    # matrix: (n_ref, n_sensor) -> pseudo-inverse: (n_sensor, n_ref)
    pinv = np.linalg.pinv(matrix)

    spatial_shape = ref_reflectance.shape[1:]
    flat = ref_reflectance.reshape(ref_reflectance.shape[0], -1)
    sensor_flat = pinv @ flat
    return sensor_flat.reshape(len(target_bands), *spatial_shape)
