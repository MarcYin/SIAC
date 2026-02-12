"""
I/O module for reading and writing geospatial data.

This module provides rioxarray-based utilities for:
- Reading rasters in various formats (GeoTIFF, COG, VRT, JP2, HDF, Zarr)
- Writing rasters with compression (GeoTIFF, COG, Zarr, NetCDF)
- Reprojection and resampling
- Geometry and vector operations

Example:
    >>> from siac.io import read_raster, write_cog, reproject_match
    >>>
    >>> # Read and process
    >>> da = read_raster("/path/to/input.tif")
    >>> aligned = reproject_match(da, target)
    >>> write_cog(aligned, "/path/to/output.tif")
"""

from siac.io.readers import (
    read_raster,
    read_raster_window,
    read_raster_at_resolution,
    read_multiband,
    read_multiband_stack,
    read_jp2,
    read_hdf_subdataset,
    read_netcdf_variable,
    read_zarr_array,
    get_raster_info,
    check_rasters_aligned,
)

from siac.io.writers import (
    write_raster,
    write_cog,
    write_dataset,
    write_zarr,
    write_netcdf,
    write_boa_products,
    write_auxiliary_products,
    write_rgb_quicklook,
)

from siac.io.reprojection import (
    reproject_to_crs,
    reproject_match,
    reproject_dataset_match,
    resample,
    resample_to_shape,
    clip_to_bounds,
    clip_to_geometry,
    transform_bounds,
    transform_points,
    get_resolution,
    get_bounds,
    get_transform,
    get_crs,
    align_grids,
    compute_common_bounds,
)

from siac.io.geometry import (
    load_aoi,
    bounds_to_polygon,
    polygon_to_bounds,
    get_raster_bounds,
    get_raster_footprint,
    get_raster_center,
    raster_to_geojson_feature,
    bounds_intersect,
    bounds_contains,
    bounds_intersection,
    bounds_union,
    bounds_area,
    create_grid_from_bounds,
    pixel_coords_to_geo,
    geo_to_pixel_coords,
)
from siac.io.s2_data_source import (
    S2DataAccess,
    S2Product,
    S2Query,
    deduplicate_products,
    fetch_s2,
    search_s2,
    select_best_product,
)
from siac.io.copernicus_dataspace import (
    CopernicusDataspaceBackend,
    download_cdse,
    search_cdse,
)
from siac.io.gcs_sentinel2 import GCSSentinel2Backend, download_gcs, search_gcs

__all__ = [
    # Readers
    "read_raster",
    "read_raster_window",
    "read_raster_at_resolution",
    "read_multiband",
    "read_multiband_stack",
    "read_jp2",
    "read_hdf_subdataset",
    "read_netcdf_variable",
    "read_zarr_array",
    "get_raster_info",
    "check_rasters_aligned",
    # Writers
    "write_raster",
    "write_cog",
    "write_dataset",
    "write_zarr",
    "write_netcdf",
    "write_boa_products",
    "write_auxiliary_products",
    "write_rgb_quicklook",
    # Reprojection
    "reproject_to_crs",
    "reproject_match",
    "reproject_dataset_match",
    "resample",
    "resample_to_shape",
    "clip_to_bounds",
    "clip_to_geometry",
    "transform_bounds",
    "transform_points",
    "get_resolution",
    "get_bounds",
    "get_transform",
    "get_crs",
    "align_grids",
    "compute_common_bounds",
    # Geometry
    "load_aoi",
    "bounds_to_polygon",
    "polygon_to_bounds",
    "get_raster_bounds",
    "get_raster_footprint",
    "get_raster_center",
    "raster_to_geojson_feature",
    "bounds_intersect",
    "bounds_contains",
    "bounds_intersection",
    "bounds_union",
    "bounds_area",
    "create_grid_from_bounds",
    "pixel_coords_to_geo",
    "geo_to_pixel_coords",
    # Sentinel-2 data access
    "S2Query",
    "S2Product",
    "S2DataAccess",
    "search_s2",
    "fetch_s2",
    "deduplicate_products",
    "select_best_product",
    # CDSE backend
    "CopernicusDataspaceBackend",
    "search_cdse",
    "download_cdse",
    # GCS backend
    "GCSSentinel2Backend",
    "search_gcs",
    "download_gcs",
]
