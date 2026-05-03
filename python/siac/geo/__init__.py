"""Geospatial utilities with no remote I/O concerns."""

from siac.geo.geometry import (
    bounds_area as bounds_area,
)
from siac.geo.geometry import (
    bounds_contains as bounds_contains,
)
from siac.geo.geometry import (
    bounds_intersect as bounds_intersect,
)
from siac.geo.geometry import (
    bounds_intersection as bounds_intersection,
)
from siac.geo.geometry import (
    bounds_to_polygon as bounds_to_polygon,
)
from siac.geo.geometry import (
    bounds_union as bounds_union,
)
from siac.geo.geometry import (
    create_grid_from_bounds as create_grid_from_bounds,
)
from siac.geo.geometry import (
    geo_to_pixel_coords as geo_to_pixel_coords,
)
from siac.geo.geometry import (
    get_raster_bounds as get_raster_bounds,
)
from siac.geo.geometry import (
    get_raster_center as get_raster_center,
)
from siac.geo.geometry import (
    get_raster_footprint as get_raster_footprint,
)
from siac.geo.geometry import (
    load_aoi as load_aoi,
)
from siac.geo.geometry import (
    pixel_coords_to_geo as pixel_coords_to_geo,
)
from siac.geo.geometry import (
    polygon_to_bounds as polygon_to_bounds,
)
from siac.geo.geometry import (
    raster_to_geojson_feature as raster_to_geojson_feature,
)
from siac.geo.reprojection import (
    align_grids as align_grids,
)
from siac.geo.reprojection import (
    clip_to_bounds as clip_to_bounds,
)
from siac.geo.reprojection import (
    clip_to_geometry as clip_to_geometry,
)
from siac.geo.reprojection import (
    compute_common_bounds as compute_common_bounds,
)
from siac.geo.reprojection import (
    get_bounds as get_bounds,
)
from siac.geo.reprojection import (
    get_crs as get_crs,
)
from siac.geo.reprojection import (
    get_resolution as get_resolution,
)
from siac.geo.reprojection import (
    get_transform as get_transform,
)
from siac.geo.reprojection import (
    reproject_dataset_match as reproject_dataset_match,
)
from siac.geo.reprojection import (
    reproject_match as reproject_match,
)
from siac.geo.reprojection import (
    reproject_to_crs as reproject_to_crs,
)
from siac.geo.reprojection import (
    resample as resample,
)
from siac.geo.reprojection import (
    resample_to_shape as resample_to_shape,
)
from siac.geo.reprojection import (
    transform_bounds as transform_bounds,
)
from siac.geo.reprojection import (
    transform_points as transform_points,
)
from siac.geo.resample import (
    axis_resolution as axis_resolution,
)
from siac.geo.resample import (
    fill_nonfinite_like_template as fill_nonfinite_like_template,
)
from siac.geo.resample import (
    resample_coefficients_to_template as resample_coefficients_to_template,
)
from siac.geo.resample import (
    resample_field_for_correction as resample_field_for_correction,
)
from siac.geo.resample import (
    resample_field_to_template as resample_field_to_template,
)
from siac.geo.resample import (
    resample_mask_to_template as resample_mask_to_template,
)
from siac.geo.resample import (
    shares_template_grid as shares_template_grid,
)
