"""BRDF product providers (MCD43, VNP43, MCD19)."""

from siac.priors.brdf.kernels import (
    BRDFKernels,
    compute_kernels,
    compute_reflectance,
    compute_white_sky_albedo,
    compute_black_sky_albedo,
)
from siac.priors.brdf.mcd43_earthaccess import (
    MCD43EarthAccessProvider,
    VNP43EarthAccessProvider,
)

__all__ = [
    "BRDFKernels",
    "compute_kernels",
    "compute_reflectance",
    "compute_white_sky_albedo",
    "compute_black_sky_albedo",
    "MCD43EarthAccessProvider",
    "VNP43EarthAccessProvider",
]
