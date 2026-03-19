"""Remote and local data-access adapters."""

from siac.adapters.data.copernicus_dataspace import (
    CopernicusDataspaceBackend,
    download_cdse,
    search_cdse,
)
from siac.adapters.data.earthaccess_catalog import (
    EarthAccessCatalog,
    EarthaccessProduct,
    ProductValidationResult,
    default_products,
)
from siac.adapters.data.earthaccess_source import EarthAccessSource
from siac.adapters.data.gcs_sentinel2 import GCSSentinel2Backend, download_gcs, search_gcs
from siac.adapters.data.s2_data_source import (
    S2DataAccess,
    S2Product,
    S2Query,
    deduplicate_products,
    fetch_s2,
    search_s2,
    select_best_product,
)

__all__ = [
    "CopernicusDataspaceBackend",
    "EarthAccessCatalog",
    "EarthAccessSource",
    "EarthaccessProduct",
    "GCSSentinel2Backend",
    "ProductValidationResult",
    "S2DataAccess",
    "S2Product",
    "S2Query",
    "deduplicate_products",
    "default_products",
    "download_cdse",
    "download_gcs",
    "fetch_s2",
    "search_cdse",
    "search_gcs",
    "search_s2",
    "select_best_product",
]
