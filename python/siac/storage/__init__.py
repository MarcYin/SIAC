"""Public storage entrypoints for raster, product, and STAC I/O."""

from siac.storage.product_writers import (
    write_aot_scatter_plot,
    write_auxiliary_products,
    write_boa_products,
    write_dataset,
    write_rgb_quicklook,
)
from siac.storage.raster_writers import write_cog, write_netcdf, write_raster, write_zarr
from siac.storage.readers import (
    check_rasters_aligned,
    get_raster_info,
    read_hdf_subdataset,
    read_jp2,
    read_multiband,
    read_multiband_stack,
    read_netcdf_variable,
    read_raster,
    read_raster_at_resolution,
    read_raster_window,
    read_zarr_array,
)
from siac.storage.stac import build_stac_item, write_stac_item

__all__ = [
    "build_stac_item",
    "check_rasters_aligned",
    "get_raster_info",
    "read_hdf_subdataset",
    "read_jp2",
    "read_multiband",
    "read_multiband_stack",
    "read_netcdf_variable",
    "read_raster",
    "read_raster_at_resolution",
    "read_raster_window",
    "read_zarr_array",
    "write_aot_scatter_plot",
    "write_auxiliary_products",
    "write_boa_products",
    "write_cog",
    "write_dataset",
    "write_netcdf",
    "write_raster",
    "write_rgb_quicklook",
    "write_stac_item",
    "write_zarr",
]
