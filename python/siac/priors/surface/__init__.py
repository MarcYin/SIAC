"""Surface prior derivation from BRDF parameters."""

from siac.priors.surface.brdf_monthly_composite import (
    MonthlyBestPixelComposite,
    build_monthly_best_pixel_composite,
)
from siac.priors.surface.brdf_monthly_database import (
    MonthlyCompositeDatabase,
    build_monthly_composite_database,
)
from siac.priors.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.priors.surface.kernel_model import KernelModelDeriver, PSFConvolver
from siac.priors.surface.spectral_mapping import (
    HyperspectralLibrary,
    SpectralMapper,
    build_default_spectral_library,
    convolve_hyperspectral_reflectance,
    map_multispectral_reflectance,
    needs_spectral_mapping,
)
from siac.priors.surface.swir_refine import (
    build_monthly_surface_prior_database,
    query_surface_prior_from_monthly_database,
    resample_geometry_for_surface_prior,
)

__all__ = [
    "KernelModelDeriver",
    "BRDFWhittakerDeriver",
    "PSFConvolver",
    "MonthlyBestPixelComposite",
    "build_monthly_best_pixel_composite",
    "MonthlyCompositeDatabase",
    "build_monthly_composite_database",
    "build_monthly_surface_prior_database",
    "HyperspectralLibrary",
    "SpectralMapper",
    "build_default_spectral_library",
    "convolve_hyperspectral_reflectance",
    "map_multispectral_reflectance",
    "needs_spectral_mapping",
    "query_surface_prior_from_monthly_database",
    "resample_geometry_for_surface_prior",
]
