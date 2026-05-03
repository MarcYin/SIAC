"""Remote and local data-access adapters."""

from siac.adapters.data.copernicus_dataspace import (
    CopernicusDataspaceBackend as CopernicusDataspaceBackend,
)
from siac.adapters.data.copernicus_dataspace import (
    download_cdse as download_cdse,
)
from siac.adapters.data.copernicus_dataspace import (
    search_cdse as search_cdse,
)
from siac.adapters.data.earthaccess_catalog import (
    EarthAccessCatalog as EarthAccessCatalog,
)
from siac.adapters.data.earthaccess_catalog import (
    EarthaccessProduct as EarthaccessProduct,
)
from siac.adapters.data.earthaccess_catalog import (
    ProductValidationResult as ProductValidationResult,
)
from siac.adapters.data.earthaccess_catalog import (
    default_products as default_products,
)
from siac.adapters.data.earthaccess_source import EarthAccessSource as EarthAccessSource
from siac.adapters.data.gcs_sentinel2 import GCSSentinel2Backend as GCSSentinel2Backend
from siac.adapters.data.gcs_sentinel2 import download_gcs as download_gcs
from siac.adapters.data.gcs_sentinel2 import search_gcs as search_gcs
from siac.adapters.data.s2_data_source import (
    S2DataAccess as S2DataAccess,
)
from siac.adapters.data.s2_data_source import (
    S2Product as S2Product,
)
from siac.adapters.data.s2_data_source import (
    S2Query as S2Query,
)
from siac.adapters.data.s2_data_source import (
    deduplicate_products as deduplicate_products,
)
from siac.adapters.data.s2_data_source import (
    fetch_s2 as fetch_s2,
)
from siac.adapters.data.s2_data_source import (
    search_s2 as search_s2,
)
from siac.adapters.data.s2_data_source import (
    select_best_product as select_best_product,
)
from siac.adapters.data.water_mask import (
    DEFAULT_WATER_MASK_CACHE_DIR as DEFAULT_WATER_MASK_CACHE_DIR,
)
from siac.adapters.data.water_mask import (
    DEFAULT_WATER_MASK_VRT_URL as DEFAULT_WATER_MASK_VRT_URL,
)
from siac.adapters.data.water_mask import (
    default_water_mask_cache_dir as default_water_mask_cache_dir,
)
from siac.adapters.data.water_mask import (
    ensure_local_water_mask_source as ensure_local_water_mask_source,
)
from siac.adapters.data.water_mask import (
    load_water_mask_subset as load_water_mask_subset,
)
from siac.adapters.data.water_mask import (
    required_water_mask_tiles as required_water_mask_tiles,
)
