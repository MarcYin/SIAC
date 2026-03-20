"""Surface prior derivation from BRDF parameters."""

from siac.algorithms.surface.brdf_monthly_composite import (
    MonthlyBestPixelComposite,
    build_monthly_best_pixel_composite,
)
from siac.algorithms.surface.brdf_monthly_database import (
    MonthlyCompositeDatabase,
    build_monthly_composite_database,
)
from siac.algorithms.surface.brdf_whittaker import BRDFWhittakerDeriver
from siac.algorithms.surface.kernel_model import KernelModelDeriver, PSFConvolver
from siac.algorithms.surface.reference_spectral import (
    load_reference_rsr,
    load_reference_rsrf,
    reference_to_sensor,
    sensor_to_reference,
)
from siac.algorithms.surface.spectral_mapping import (
    HyperspectralLibrary,
    SpectralMapper,
    SpectralMappingConfig,
    convolve_hyperspectral_reflectance,
    map_multispectral_reflectance,
    needs_spectral_mapping,
)
from siac.algorithms.surface.swir_refine import (
    build_monthly_surface_prior_database,
    query_surface_prior_from_monthly_database,
    resample_geometry_for_surface_prior,
)

__all__ = [
    "KernelModelDeriver",
    "BRDFWhittakerDeriver",
    "PSFConvolver",
    "load_reference_rsr",
    "load_reference_rsrf",
    "sensor_to_reference",
    "reference_to_sensor",
    "MonthlyBestPixelComposite",
    "build_monthly_best_pixel_composite",
    "MonthlyCompositeDatabase",
    "build_monthly_composite_database",
    "build_monthly_surface_prior_database",
    "HyperspectralLibrary",
    "SpectralMappingConfig",
    "SpectralMapper",
    "convolve_hyperspectral_reflectance",
    "map_multispectral_reflectance",
    "needs_spectral_mapping",
    "query_surface_prior_from_monthly_database",
    "resample_geometry_for_surface_prior",
]
