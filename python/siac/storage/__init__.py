"""Raster and product storage helpers."""

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
from siac.storage.writers import (
    write_auxiliary_products,
    write_boa_products,
    write_cog,
    write_dataset,
    write_netcdf,
    write_raster,
    write_rgb_quicklook,
    write_zarr,
)

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
