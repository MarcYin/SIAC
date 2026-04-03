"""Public storage entrypoints for raster, product, and STAC I/O."""

from siac.storage.product_writers import (
    write_aot_scatter_plot,
    write_auxiliary_products,
    write_boa_products,
    write_cloud_mask_preview,
    write_dataset,
    write_false_colour_preview,
    write_field_preview,
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
from siac.storage.stac import build_stac_item, build_stac_item_from_result, write_stac_item

__all__ = [
    "build_stac_item",
    "build_stac_item_from_result",
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
    "write_cloud_mask_preview",
    "write_cog",
    "write_dataset",
    "write_false_colour_preview",
    "write_field_preview",
    "write_netcdf",
    "write_raster",
    "write_rgb_quicklook",
    "write_stac_item",
    "write_zarr",
]
