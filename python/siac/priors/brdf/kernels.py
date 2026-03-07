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

from typing import Any, Literal

import numpy as np
import xarray as xr


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
    ):
        self.hb = hb
        self.br = br
        self._rust_kernels: Any | None = None  # lazily initialized on first compute()

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
        # Convert xarray to numpy if needed
        is_xarray = isinstance(vza, xr.DataArray)
        if is_xarray:
            vza_np = vza.values
            sza_np = sza.values
            raa_np = raa.values
            template = vza
        else:
            vza_np = np.asarray(vza)
            sza_np = np.asarray(sza)
            raa_np = np.asarray(raa)
            template = None

        if vza_np.shape != sza_np.shape or vza_np.shape != raa_np.shape:
            raise ValueError("vza, sza, and raa must have the same shape")

        original_shape = vza_np.shape
        vza_in = np.ascontiguousarray(vza_np.reshape(1, -1), dtype=np.float64)
        sza_in = np.ascontiguousarray(sza_np.reshape(1, -1), dtype=np.float64)
        raa_in = np.ascontiguousarray(raa_np.reshape(1, -1), dtype=np.float64)
        if self._rust_kernels is None:
            from siac._rust import RossThickLiSparse as _RustKernels  # noqa: PLC0415 - lazy; siac._rust is optional at import time

            self._rust_kernels = _RustKernels(self.hb, self.br)
        k_vol, k_geo = self._rust_kernels.compute(vza_in, sza_in, raa_in)
        k_vol = np.asarray(k_vol, dtype=np.float64).reshape(original_shape)
        k_geo = np.asarray(k_geo, dtype=np.float64).reshape(original_shape)

        # Convert back to xarray if needed
        if is_xarray:
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
    return f0 + f1 * k_vol + f2 * k_geo


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
    # Integration weights for white-sky albedo
    # These are pre-computed integrals of the kernels
    g0 = 1.0
    g1 = 0.189184  # Integral of Ross-Thick kernel
    g2 = -1.377622  # Integral of Li-Sparse kernel

    return f0 * g0 + f1 * g1 + f2 * g2


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
    # Polynomial coefficients for black-sky albedo
    # a_k(sza) = g0_k + g1_k * sza + g2_k * sza^2 + g3_k * sza^3

    # Ross-Thick polynomial coefficients
    g_ross = np.array([-0.007574, -0.070987, 0.307588, 0.0])

    # Li-Sparse polynomial coefficients
    g_li = np.array([-1.284909, -0.166314, 0.041840, 0.0])

    sza2 = sza**2
    sza3 = sza**3

    a_ross = g_ross[0] + g_ross[1] * sza + g_ross[2] * sza2 + g_ross[3] * sza3
    a_li = g_li[0] + g_li[1] * sza + g_li[2] * sza2 + g_li[3] * sza3

    return f0 + f1 * a_ross + f2 * a_li
