"""
BRDF kernel calculations.

Implements Ross-Thick and Li-Sparse BRDF kernels used in the MODIS
MCD43 BRDF model. Provides both pure Python and Rust-accelerated
implementations.

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
from typing import Literal

import numpy as np
import xarray as xr

logger = logging.getLogger(__name__)

# Try to import Rust implementation
try:
    from siac._rust import RossThickLiSparse as _RustKernels

    _HAS_RUST = True
    logger.debug("Using Rust BRDF kernel implementation")
except ImportError:
    _HAS_RUST = False
    logger.debug("Rust BRDF kernels not available, using Python implementation")


class BRDFKernels:
    """
    Ross-Thick Li-Sparse BRDF kernel calculator.

    Computes the volumetric (Ross-Thick) and geometric (Li-Sparse) kernel
    values for given viewing geometry.

    Args:
        hb: Height-to-breadth ratio for Li kernel (default: 2.0 for MODIS)
        br: Crown relative height for Li kernel (default: 1.0 for sparse)
        use_rust: Whether to use Rust acceleration if available
    """

    def __init__(
        self,
        hb: float = 2.0,
        br: float = 1.0,
        use_rust: bool = True,
    ):
        self.hb = hb
        self.br = br
        self._use_rust = use_rust and _HAS_RUST

        if self._use_rust:
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

        # Compute kernels
        if self._use_rust and vza_np.ndim == 2:
            k_vol, k_geo = self._rust_kernels.compute(vza_np, sza_np, raa_np)
        else:
            k_vol, k_geo = self._compute_python(vza_np, sza_np, raa_np)

        # Convert back to xarray if needed
        if is_xarray:
            k_vol = xr.DataArray(k_vol, dims=template.dims, coords=template.coords)
            k_geo = xr.DataArray(k_geo, dims=template.dims, coords=template.coords)

        return k_vol, k_geo

    def _compute_python(
        self,
        vza: np.ndarray,
        sza: np.ndarray,
        raa: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Pure Python kernel computation."""
        cos_sza = np.cos(sza)
        cos_vza = np.cos(vza)
        sin_sza = np.sin(sza)
        sin_vza = np.sin(vza)
        cos_raa = np.cos(np.abs(raa))

        # Phase angle
        cos_phase = np.clip(
            cos_sza * cos_vza + sin_sza * sin_vza * cos_raa,
            -1.0,
            1.0,
        )
        phase = np.arccos(cos_phase)

        # Ross-Thick kernel
        k_vol = self._ross_thick(cos_sza, cos_vza, cos_phase, phase)

        # Li-Sparse kernel
        k_geo = self._li_sparse(
            cos_sza, cos_vza, sin_sza, sin_vza, cos_raa, cos_phase
        )

        return k_vol, k_geo

    def _ross_thick(
        self,
        cos_sza: np.ndarray,
        cos_vza: np.ndarray,
        cos_phase: np.ndarray,
        phase: np.ndarray,
    ) -> np.ndarray:
        """Ross-Thick volumetric scattering kernel."""
        denom = cos_sza + cos_vza

        # Avoid division by zero
        with np.errstate(divide="ignore", invalid="ignore"):
            k_vol = ((np.pi / 2 - phase) * cos_phase + np.sin(phase)) / denom - np.pi / 2

        # Handle edge cases
        k_vol = np.where(np.abs(denom) < 1e-10, 0.0, k_vol)

        return k_vol

    def _li_sparse(
        self,
        cos_sza: np.ndarray,
        cos_vza: np.ndarray,
        sin_sza: np.ndarray,
        sin_vza: np.ndarray,
        cos_raa: np.ndarray,
        _cos_phase: np.ndarray,
    ) -> np.ndarray:
        """Li-Sparse geometric scattering kernel."""
        # Prime angles (scaled by br for sparse vegetation)
        with np.errstate(divide="ignore", invalid="ignore"):
            tan_sza = sin_sza / np.maximum(cos_sza, 1e-10)
            tan_vza = sin_vza / np.maximum(cos_vza, 1e-10)

        tan_sza_prime = self.br * tan_sza
        tan_vza_prime = self.br * tan_vza

        sza_prime = np.arctan(tan_sza_prime)
        vza_prime = np.arctan(tan_vza_prime)

        cos_sza_prime = np.cos(sza_prime)
        cos_vza_prime = np.cos(vza_prime)
        sin_sza_prime = np.sin(sza_prime)
        sin_vza_prime = np.sin(vza_prime)

        # Prime phase angle
        cos_phase_prime = np.clip(
            cos_sza_prime * cos_vza_prime + sin_sza_prime * sin_vza_prime * cos_raa,
            -1.0,
            1.0,
        )

        # Distance term
        d2 = (
            tan_sza_prime**2
            + tan_vza_prime**2
            - 2 * tan_sza_prime * tan_vza_prime * cos_raa
        )
        d2 = np.maximum(d2, 0.0)

        # Secant values
        sec_sza_prime = 1.0 / np.maximum(cos_sza_prime, 1e-10)
        sec_vza_prime = 1.0 / np.maximum(cos_vza_prime, 1e-10)

        # Overlap function
        sin_raa = np.sin(np.abs(np.arccos(cos_raa)))
        cost_arg = self.hb * np.sqrt(
            d2 + (tan_sza_prime * tan_vza_prime * sin_raa) ** 2
        ) / (sec_sza_prime + sec_vza_prime)
        cost = np.clip(cost_arg, -1.0, 1.0)
        t = np.arccos(cost)
        overlap = (1.0 / np.pi) * (t - np.sin(t) * np.cos(t)) * (
            sec_sza_prime + sec_vza_prime
        )

        # Final Li-Sparse kernel
        k_geo = (
            overlap
            - sec_sza_prime
            - sec_vza_prime
            + 0.5 * (1 + cos_phase_prime) * sec_sza_prime * sec_vza_prime
        )

        return k_geo


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
