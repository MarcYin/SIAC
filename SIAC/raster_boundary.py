from osgeo import osr, ogr, gdal
import json
import numpy  as np

def get_boundary(raster, to_wgs84 = True):
    '''
    Take a raster and find boundary of the raster 
    if to_wgs84 = True, transform them into wgs84 coordinates
    with at least meters resoluton.
    Return geojson(s)
    '''
    g    = gdal.Open(raster)
    geom = g.GetGeoTransform()
    x_coords = geom[0] + np.arange(g.RasterXSize + 1) * geom[1]
    y_coords = geom[3] + np.arange(g.RasterYSize + 1) * geom[5]
    north_boundary = np.dstack([x_coords,       np.ones(g.RasterXSize + 1) * y_coords[0]])
    south_boundary = np.dstack([x_coords[::-1], np.ones(g.RasterXSize + 1) * y_coords[-1]])
    east_boundary  = np.dstack([x_coords[-1] * np.ones(g.RasterYSize + 1), y_coords])
    west_boundary  = np.dstack([x_coords[0]  * np.ones(g.RasterYSize + 1), y_coords[::-1]])
    boundary = np.hstack([north_boundary, east_boundary, south_boundary, west_boundary]).tolist()   
    polygon  = {'type': 'Polygon', 'coordinates': boundary}
    geometry = ogr.CreateGeometryFromJson(json.dumps(polygon))
    sgeometry = geometry.SimplifyPreserveTopology(abs(geom[1])/1000)
    spolygon = json.loads(sgeometry.ExportToJson())
    sourceprj = osr.SpatialReference()
    sourceprj.ImportFromWkt(g.GetProjection())
    try:
        urn = sourceprj.ExportToXML().split('srsID')[1].split('=')[1].split('</gml:name>')[0].replace('"', '').replace('>', '')
        geojson = {                       
          "type": "FeatureCollection",                               
          "name": "AOI",          
          "crs": { "type": "name", "properties": { "name": urn} },
          "features": [           
          { "type": "Feature", "properties": {}, "geometry": spolygon}
                      ]                       
         }
    except:
        proj4 = sourceprj.ExportToProj4()
        geojson = {                       
          "type": "FeatureCollection",                               
          "name": "AOI",          
          "crs": { "type": "name", "properties": { "name": proj4, "type": "proj4"} },
          "features": [           
          { "type": "Feature", "properties": {}, "geometry": spolygon}
                      ]                       
         }

    geojson = json.dumps(geojson)

    if to_wgs84:
        sourceprj = osr.SpatialReference()
        sourceprj.ImportFromWkt(g.GetProjection())

        targetprj = osr.SpatialReference()
        targetprj.ImportFromProj4('+proj=longlat +datum=WGS84 +no_defs ')
        if int(gdal.VersionInfo()) > 2500000:
          targetprj.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        transform = osr.CoordinateTransformation(sourceprj, targetprj)
        sgeometry.Transform(transform)

        # stgeometry = geometry.SimplifyPreserveTopology(abs(geom[1]) / 10. * 0.000001) 
        stpolygon = json.loads(sgeometry.ExportToJson())
        
        tgeojson = {                       
          "type": "FeatureCollection",                               
          "name": "AOI",          
          "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::4326"} },
          "features": [           
          { "type": "Feature", "properties": {}, "geometry": stpolygon}
                      ]                       
         }
        tgeojson = json.dumps(tgeojson)

        return tgeojson, geojson
    else:
        return geojson

