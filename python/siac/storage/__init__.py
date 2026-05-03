"""Public storage entrypoints for raster, product, and STAC I/O."""

from siac.storage.product_writers import (
    write_aot_scatter_plot as write_aot_scatter_plot,
)
from siac.storage.product_writers import (
    write_auxiliary_products as write_auxiliary_products,
)
from siac.storage.product_writers import (
    write_boa_products as write_boa_products,
)
from siac.storage.product_writers import (
    write_cloud_mask_preview as write_cloud_mask_preview,
)
from siac.storage.product_writers import (
    write_dataset as write_dataset,
)
from siac.storage.product_writers import (
    write_false_colour_preview as write_false_colour_preview,
)
from siac.storage.product_writers import (
    write_field_preview as write_field_preview,
)
from siac.storage.product_writers import (
    write_rgb_quicklook as write_rgb_quicklook,
)
from siac.storage.raster_writers import write_cog as write_cog
from siac.storage.raster_writers import write_netcdf as write_netcdf
from siac.storage.raster_writers import write_raster as write_raster
from siac.storage.raster_writers import write_zarr as write_zarr
from siac.storage.readers import (
    check_rasters_aligned as check_rasters_aligned,
)
from siac.storage.readers import (
    get_raster_info as get_raster_info,
)
from siac.storage.readers import (
    read_hdf_subdataset as read_hdf_subdataset,
)
from siac.storage.readers import (
    read_jp2 as read_jp2,
)
from siac.storage.readers import (
    read_multiband as read_multiband,
)
from siac.storage.readers import (
    read_multiband_stack as read_multiband_stack,
)
from siac.storage.readers import (
    read_netcdf_variable as read_netcdf_variable,
)
from siac.storage.readers import (
    read_raster as read_raster,
)
from siac.storage.readers import (
    read_raster_at_resolution as read_raster_at_resolution,
)
from siac.storage.readers import (
    read_raster_window as read_raster_window,
)
from siac.storage.readers import (
    read_zarr_array as read_zarr_array,
)
from siac.storage.stac import build_stac_item_from_result as build_stac_item_from_result
