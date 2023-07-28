# copyied and adapted from https://gist.github.com/PawaritL/ec7136c0b718ca65db6df1c33fd1bb11

from enum import Enum
from typing import List, Union
import shapely
from shapely import affinity
from shapely.geometry import GeometryCollection, Polygon, mapping

import math
import copy
import json
from typing import Union, List
from shapely.geometry import Polygon, LineString, GeometryCollection
from shapely.ops import split



class OutputFormat(Enum):
    Geojson = 'geojson'
    Polygons = 'polygons'
    GeometryCollection = 'geometrycollection'


def check_crossing(lon1: float, lon2: float, validate: bool = True, dlon_threshold: float = 180.0):
    """
    Assuming a minimum travel distance between two provided longitude coordinates,
    checks if the 180th meridian (antimeridian) is crossed.
    """
    if validate and (any(abs(x) > dlon_threshold) for x in [lon1, lon2]):
        raise ValueError("longitudes must be in degrees [-180.0, 180.0]")   
    return abs(lon2 - lon1) > dlon_threshold


def translate_polygons(
    geometry_collection: GeometryCollection,
) -> List[Polygon]:

    geo_polygons = []
    for polygon in geometry_collection.geoms:
        (minx, _, maxx, _) = polygon.bounds
        if minx < -180:
            geo_polygon = affinity.translate(polygon, xoff = 360)
        elif maxx > 180:
            geo_polygon = affinity.translate(polygon, xoff = -360)
        else:
            geo_polygon = polygon
        geo_polygon = shapely.geometry.polygon.orient(geo_polygon, sign=1.0)
        geo_polygons.append(geo_polygon)

    return geo_polygons


def split_polygon(coords, validate: bool = False):
    """
    Given a GeoJSON representation of a Polygon, returns a collection of
    'antimeridian-safe' constituent polygons split at the 180th meridian, 
    ensuring compliance with GeoJSON standards (https://tools.ietf.org/html/rfc7946#section-3.1.9)
    Assumptions:
      - Any two consecutive points with over 180 degrees difference in
        longitude are assumed to cross the antimeridian
      - The polygon spans less than 360 degrees in longitude (i.e. does not wrap around the globe)
      - However, the polygon may cross the antimeridian on multiple occasions
    Parameters:
        geojson (dict): GeoJSON of input polygon to be split. For example:
                        {
                        "type": "Polygon",
                        "coordinates": [
                          [
                            [179.0, 0.0], [-179.0, 0.0], [-179.0, 1.0],
                            [179.0, 1.0], [179.0, 0.0]
                          ]
                        ]
                        }
        output_format (str): Available options: "geojson", "polygons", "geometrycollection"
                             If "geometrycollection" returns a Shapely GeometryCollection.
                             Otherwise, returns a list of either GeoJSONs or Shapely Polygons
        validate (bool): Checks if all longitudes are within [-180.0, 180.0]
      
    Returns:
        List[dict]/List[Polygon]/GeometryCollection: antimeridian-safe polygon(s)
    """

    coords_shift = copy.deepcopy(coords)
    shell_minx = shell_maxx = None
    split_meridians = set()
    splitter = None

    for ring_index, ring in enumerate(coords_shift):
        if len(ring) < 1: 
            continue
        else:
            ring_minx = ring_maxx = ring[0][0]
            crossings = 0

        for coord_index, (lon, _) in enumerate(ring[1:], start=1):
            lon_prev = ring[coord_index - 1][0] # [0] corresponds to longitude coordinate
            if check_crossing(lon, lon_prev, validate=validate):
                direction = math.copysign(1, lon - lon_prev)
                coords_shift[ring_index][coord_index][0] = lon - (direction * 360.0)
                crossings += 1

            x_shift = coords_shift[ring_index][coord_index][0]
            if x_shift < ring_minx: ring_minx = x_shift
            if x_shift > ring_maxx: ring_maxx = x_shift

        # Ensure that any holes remain contained within the (translated) outer shell
        if (ring_index == 0): # by GeoJSON definition, first ring is the outer shell
            shell_minx, shell_maxx = (ring_minx, ring_maxx)
        elif (ring_minx < shell_minx):
            ring_shift = [[x + 360, y] for (x, y) in coords_shift[ring_index]]
            coords_shift[ring_index] = ring_shift
            ring_minx, ring_maxx = (x + 360 for x in (ring_minx, ring_maxx))
        elif (ring_maxx > shell_maxx):
            ring_shift = [[x - 360, y] for (x, y) in coords_shift[ring_index]]
            coords_shift[ring_index] = ring_shift
            ring_minx, ring_maxx = (x - 360 for x in (ring_minx, ring_maxx))

        if crossings: # keep track of meridians to split on
            if ring_minx < -180: split_meridians.add(-180)
            if ring_maxx > 180: split_meridians.add(180)

    n_splits = len(split_meridians)
    if n_splits > 1:
        raise NotImplementedError(
            """Splitting a Polygon by multiple meridians (MultiLineString) 
               not supported by Shapely"""
        )
    elif n_splits == 1:
        split_lon = next(iter(split_meridians))
        meridian = [[split_lon, -90.0], [split_lon, 90.0]]
        splitter = LineString(meridian)

    shell, *holes = coords_shift if splitter else coords
    if splitter: split_polygons = split(Polygon(shell, holes), splitter)
    else: split_polygons = GeometryCollection([Polygon(shell, holes)])
        
    geo_polygons = list(translate_polygons(split_polygons))  
    list_of_coords = [list(polygon.exterior.coords) for polygon in geo_polygons]
    return list_of_coords


if __name__ == "__main__":
        
    coords = [[[ 178.92616652,   63.88152371],
    [-178.83953018,   63.92717598],
    [-178.90672471,   64.91194028],
    [ 178.77777827,   64.8642401 ]]]
    
    # Polygon(tuple(coords[0]))

    list_of_coords = split_polygon(coords)