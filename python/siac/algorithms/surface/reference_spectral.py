"""Reference-basis spectral projection helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

    import xarray as xr
    from numpy.typing import NDArray

    from siac.domain import SensorBand


def _trapezoid(y: NDArray, x: NDArray) -> float:
    """Compatibility wrapper for NumPy 1.x/2.x integration API."""
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.add.reduce((y[1:] + y[:-1]) * np.diff(x) * 0.5))


_MODIS_LAND_BANDS = (
    ("Band1", 645.0, 50.0),
    ("Band2", 858.5, 35.0),
    ("Band3", 469.0, 20.0),
    ("Band4", 555.0, 20.0),
    ("Band5", 1240.0, 20.0),
    ("Band6", 1640.0, 24.0),
    ("Band7", 2130.0, 50.0),
)


def load_reference_rsr(
    reference_sensor: str = "MODIS",
) -> dict[str, tuple[NDArray, NDArray]]:
    """Load tabulated spectral response functions for a reference sensor."""
    if reference_sensor.upper() == "MODIS":
        result: dict[str, tuple[NDArray, NDArray]] = {}
        for band_name, center_nm, fwhm_nm in _MODIS_LAND_BANDS:
            sigma = fwhm_nm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
            wl = np.linspace(center_nm - 3.0 * fwhm_nm, center_nm + 3.0 * fwhm_nm, 101)
            resp = np.exp(-0.5 * ((wl - center_nm) / sigma) ** 2)
            result[band_name] = (wl, resp)
        return result

    raise ValueError(f"Unknown reference sensor: {reference_sensor!r}. Available: 'MODIS'")


def _build_conversion_matrix(
    sensor_bands: Sequence[SensorBand],
    ref_bands: dict[str, tuple[NDArray, NDArray]],
    common_wl: NDArray | None = None,
) -> tuple[NDArray, list[str]]:
    """Build a spectral conversion matrix."""
    if common_wl is None:
        common_wl = np.linspace(350.0, 2500.0, 2151)

    ref_names = list(ref_bands.keys())
    n_ref = len(ref_names)
    n_sensor = len(sensor_bands)

    ref_resp = np.zeros((n_ref, len(common_wl)))
    for i, name in enumerate(ref_names):
        wl_r, resp_r = ref_bands[name]
        interp_r = np.interp(common_wl, wl_r, resp_r, left=0.0, right=0.0)
        norm = _trapezoid(interp_r, common_wl)
        ref_resp[i] = interp_r / norm if norm > 0 else interp_r

    sensor_resp = np.zeros((n_sensor, len(common_wl)))
    for j, band in enumerate(sensor_bands):
        sr = band.effective_response(common_wl)
        norm = _trapezoid(sr, common_wl)
        sensor_resp[j] = sr / norm if norm > 0 else sr

    matrix = np.zeros((n_ref, n_sensor))
    for i in range(n_ref):
        for j in range(n_sensor):
            overlap = _trapezoid(ref_resp[i] * sensor_resp[j], common_wl)
            matrix[i, j] = overlap if overlap > 1e-6 else 0.0

    return matrix, ref_names


def sensor_to_reference(
    sensor_reflectance: xr.Dataset,
    sensor_bands: Sequence[SensorBand],
    reference_sensor: str = "MODIS",
) -> NDArray:
    """Convolve sensor-band reflectance to a reference basis."""
    ref_rsr = load_reference_rsr(reference_sensor)
    matrix, ref_names = _build_conversion_matrix(sensor_bands, ref_rsr)

    band_names = [band.name for band in sensor_bands]
    sensor_vals = np.stack(
        [sensor_reflectance[band_name].values for band_name in band_names if band_name in sensor_reflectance],
        axis=0,
    )

    spatial_shape = sensor_vals.shape[1:]
    flat = sensor_vals.reshape(sensor_vals.shape[0], -1)
    ref_flat = matrix @ flat
    return ref_flat.reshape(len(ref_names), *spatial_shape)


def reference_to_sensor(
    ref_reflectance: NDArray,
    target_bands: Sequence[SensorBand],
    reference_sensor: str = "MODIS",
) -> NDArray:
    """Project reference-basis reflectance back to sensor bands."""
    ref_rsr = load_reference_rsr(reference_sensor)
    matrix, _ = _build_conversion_matrix(target_bands, ref_rsr)
    pinv = np.linalg.pinv(matrix)

    spatial_shape = ref_reflectance.shape[1:]
    flat = ref_reflectance.reshape(ref_reflectance.shape[0], -1)
    sensor_flat = pinv @ flat
    return sensor_flat.reshape(len(target_bands), *spatial_shape)


__all__ = [
    "load_reference_rsr",
    "sensor_to_reference",
    "reference_to_sensor",
]
