"""
BRDF kernel calculations.

Implements Ross-Thick and Li-Sparse BRDF kernels used in the MODIS
MCD43 BRDF model via the Rust extension shipped with SIAC.

The BRDF model is:
    ρ(θv, θs, φ) = f0 + f1 * K_vol + f2 * K_geo

where:
    f0, f1, f2 are the isotropic, volumetric, and geometric kernel weights
    K_vol is the Ross-Thick volumetric scattering kernel
    K_geo is the Li-Sparse geometric scattering kernel
    θv, θs are view and solar zenith angles
    φ is the relative azimuth angle
"""

from __future__ import annotations

import logging
from typing import Any, Literal, cast

import numpy as np
import xarray as xr

from siac._rust_compat import RossThickLiSparse as _RustKernels

logger = logging.getLogger(__name__)
FloatArray = np.ndarray[Any, np.dtype[np.float64]]


def _data_array_template(*arrays: object) -> xr.DataArray | None:
    for array in arrays:
        if isinstance(array, xr.DataArray):
            return array
    return None


def _to_numpy(array: np.ndarray | xr.DataArray) -> FloatArray:
    if isinstance(array, xr.DataArray):
        return cast("FloatArray", np.asarray(array.values, dtype=np.float64))
    return cast("FloatArray", np.asarray(array, dtype=np.float64))


class BRDFKernels:
    """
    Ross-Thick Li-Sparse BRDF kernel calculator.

    Computes the volumetric (Ross-Thick) and geometric (Li-Sparse) kernel
    values for given viewing geometry.

    Args:
        hb: Height-to-breadth ratio for Li kernel (default: 2.0 for MODIS)
        br: Crown relative height for Li kernel (default: 1.0 for sparse)
    """

    def __init__(
        self,
        hb: float = 2.0,
        br: float = 1.0,
    ) -> None:
        self.hb = hb
        self.br = br
        self._rust_kernels = _RustKernels(hb, br)

    def compute(
        self,
        vza: np.ndarray | xr.DataArray,
        sza: np.ndarray | xr.DataArray,
        raa: np.ndarray | xr.DataArray,
    ) -> tuple[np.ndarray | xr.DataArray, np.ndarray | xr.DataArray]:
        """
        Compute Ross-Thick and Li-Sparse kernel values.

        Args:
            vza: View zenith angle in radians
            sza: Solar zenith angle in radians
            raa: Relative azimuth angle in radians

        Returns:
            Tuple of (Ross kernel, Li kernel) arrays
        """
        template = _data_array_template(vza, sza, raa)
        vza_np = _to_numpy(vza)
        sza_np = _to_numpy(sza)
        raa_np = _to_numpy(raa)

        if vza_np.shape != sza_np.shape or vza_np.shape != raa_np.shape:
            raise ValueError("vza, sza, and raa must have the same shape")

        original_shape = vza_np.shape
        vza_in = np.ascontiguousarray(vza_np.reshape(1, -1), dtype=np.float64)
        sza_in = np.ascontiguousarray(sza_np.reshape(1, -1), dtype=np.float64)
        raa_in = np.ascontiguousarray(raa_np.reshape(1, -1), dtype=np.float64)
        k_vol, k_geo = self._rust_kernels.compute(vza_in, sza_in, raa_in)
        k_vol = np.asarray(k_vol, dtype=np.float64).reshape(original_shape)
        k_geo = np.asarray(k_geo, dtype=np.float64).reshape(original_shape)

        # Convert back to xarray if needed
        if template is not None:
            k_vol = xr.DataArray(k_vol, dims=template.dims, coords=template.coords)
            k_geo = xr.DataArray(k_geo, dims=template.dims, coords=template.coords)

        return k_vol, k_geo


def compute_kernels(
    vza: np.ndarray | xr.DataArray,
    sza: np.ndarray | xr.DataArray,
    raa: np.ndarray | xr.DataArray,
    kernel_type: Literal["modis", "roujean"] = "modis",
) -> tuple[np.ndarray | xr.DataArray, np.ndarray | xr.DataArray]:
    """
    Convenience function to compute BRDF kernels.

    Args:
        vza: View zenith angle in radians
        sza: Solar zenith angle in radians
        raa: Relative azimuth angle in radians
        kernel_type: Kernel type ("modis" for Ross-Thick Li-Sparse)

    Returns:
        Tuple of (volumetric kernel, geometric kernel)
    """
    if kernel_type == "modis":
        kernels = BRDFKernels(hb=2.0, br=1.0)
    else:
        raise ValueError(f"Unknown kernel type: {kernel_type}")

    return kernels.compute(vza, sza, raa)


def compute_reflectance(
    f0: np.ndarray | xr.DataArray,
    f1: np.ndarray | xr.DataArray,
    f2: np.ndarray | xr.DataArray,
    k_vol: np.ndarray | xr.DataArray,
    k_geo: np.ndarray | xr.DataArray,
) -> np.ndarray | xr.DataArray:
    """
    Compute surface reflectance from BRDF parameters and kernels.

    Args:
        f0: Isotropic kernel weight
        f1: Volumetric kernel weight
        f2: Geometric kernel weight
        k_vol: Volumetric kernel values
        k_geo: Geometric kernel values

    Returns:
        Surface reflectance
    """
    reflectance: np.ndarray | xr.DataArray = f0 + f1 * k_vol + f2 * k_geo
    return reflectance


def compute_white_sky_albedo(
    f0: np.ndarray | xr.DataArray,
    f1: np.ndarray | xr.DataArray,
    f2: np.ndarray | xr.DataArray,
) -> np.ndarray | xr.DataArray:
    """
    Compute white-sky (bi-hemispherical) albedo.

    This is the albedo under completely diffuse illumination.

    Args:
        f0, f1, f2: BRDF kernel weights

    Returns:
        White-sky albedo
    """
    # Pre-computed analytical integrals of Ross-Thick and Li-Sparse kernels
    # over the full hemisphere (Lucht et al. 2000, Table 1).
    g0 = 1.0  # Integral of isotropic kernel (≡ 1)
    g1 = 0.189184  # Integral of Ross-Thick volumetric kernel
    g2 = -1.377622  # Integral of Li-Sparse geometric kernel (negative is expected)

    albedo: np.ndarray | xr.DataArray = f0 * g0 + f1 * g1 + f2 * g2

    # Flag pixels where the albedo falls outside the physically plausible
    # range [0, 1].  The negative g2 weight on a large geometric kernel can
    # drive the result below zero for certain BRDF parameter combinations.
    albedo_arr = np.asarray(albedo)
    n_negative = int(np.sum(albedo_arr[np.isfinite(albedo_arr)] < 0.0))
    n_over_one = int(np.sum(albedo_arr[np.isfinite(albedo_arr)] > 1.0))
    if n_negative > 0 or n_over_one > 0:
        logger.warning(
            "White-sky albedo out of [0, 1]: %d negative, %d > 1 "
            "(total finite pixels: %d). Clamping to [0, 1].",
            n_negative,
            n_over_one,
            int(np.sum(np.isfinite(albedo_arr))),
        )
        albedo = np.clip(albedo, 0.0, 1.0)

    return albedo


def compute_black_sky_albedo(
    f0: np.ndarray | xr.DataArray,
    f1: np.ndarray | xr.DataArray,
    f2: np.ndarray | xr.DataArray,
    sza: np.ndarray | xr.DataArray,
) -> np.ndarray | xr.DataArray:
    """
    Compute black-sky (directional-hemispherical) albedo.

    This is the albedo under direct illumination at solar zenith angle sza.

    Args:
        f0, f1, f2: BRDF kernel weights
        sza: Solar zenith angle in radians

    Returns:
        Black-sky albedo
    """
    # Polynomial coefficients for black-sky albedo approximation
    # (Lucht et al. 2000, Table 2).
    # a_k(sza) = g0_k + g1_k * sza + g2_k * sza^2 + g3_k * sza^3

    # Ross-Thick volumetric polynomial coefficients
    g_ross = np.asarray([-0.007574, -0.070987, 0.307588, 0.0], dtype=np.float64)

    # Li-Sparse geometric polynomial coefficients
    g_li = np.asarray([-1.284909, -0.166314, 0.041840, 0.0], dtype=np.float64)

    sza2: np.ndarray | xr.DataArray = sza * sza
    sza3: np.ndarray | xr.DataArray = sza2 * sza

    a_ross = g_ross[0] + g_ross[1] * sza + g_ross[2] * sza2 + g_ross[3] * sza3
    a_li = g_li[0] + g_li[1] * sza + g_li[2] * sza2 + g_li[3] * sza3

    albedo: np.ndarray | xr.DataArray = f0 + f1 * a_ross + f2 * a_li
    return albedo
