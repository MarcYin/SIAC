import os
import json
import numpy as np
from osgeo import gdal, osr, ogr
from SIAC.split_polygon import split_polygon

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


def subset_raster(raster_file, subset_file, out_raster_file, xRes, yRes):
    
    g = gdal.Open(raster_file)
    geo_t = g.GetGeoTransform()
    x_size, y_size = g.RasterXSize, g.RasterYSize     

    xmin, xmax = min(geo_t[0], geo_t[0] + x_size * geo_t[1]), \
                    max(geo_t[0], geo_t[0] + x_size * geo_t[1])  
    ymin, ymax = min(geo_t[3], geo_t[3] + y_size * geo_t[5]), \
                max(geo_t[3], geo_t[3] + y_size * geo_t[5])
    raster_bounds = [xmin, ymin, xmax, ymax]

    dstSRS     = osr.SpatialReference( )
    raster_wkt = g.GetProjection()
    dstSRS.ImportFromWkt(raster_wkt)

    geojson = get_boundary(raster_file)[0]    
    coords = np.array(json.loads(geojson)['features'][0]['geometry']['coordinates'])
    list_of_coords = split_polygon(coords)

    g = gdal.Open(subset_file)
    no_data_value = g.GetRasterBand(1).GetNoDataValue()
    if no_data_value is None:
        no_data_value = -9999

    dats = []
    for coords in list_of_coords:
        xmax, ymax = np.max(coords, axis=0)
        xmin, ymin = np.min(coords, axis=0)
        bounds = [xmin, ymin, xmax, ymax]
        g  = gdal.Warp('', subset_file, format='MEM', outputBounds=bounds,  warpOptions = [ 'CUTLINE_ALL_TOUCHED=TRUE' ], dstNodata = no_data_value)
        gg = gdal.Warp('', g, format='MEM', dstSRS = dstSRS, outputBounds = raster_bounds, xRes = xRes, yRes = yRes)
        data = gg.ReadAsArray()
        subset_array_geotrans = gg.GetGeoTransform()
        data = np.ma.array(data, mask = data == no_data_value)
        dats.append(data)

    data = np.ma.mean(dats, axis=0)
    data = np.array(data)

    g = gdal.Open(subset_file)
    if data.ndim == 2:
        data = data[np.newaxis, :, :]
    nbands = data.shape[0]
    y_size, x_size = data.shape[1], data.shape[2]

    # save the data to a raster file
    driver = gdal.GetDriverByName('GTiff')
    dst_ds = driver.Create(out_raster_file, x_size, y_size, nbands, gdal.GDT_Float32, options = [ 'COMPRESS=LZW', 'TILED=YES'])
    dst_ds.SetGeoTransform(subset_array_geotrans)
    dst_ds.SetProjection(raster_wkt)
    for i in range(nbands):
        metadata = g.GetRasterBand(i+1).GetMetadata()
        dst_ds.GetRasterBand(i+1).WriteArray(data[i,:,:])
        dst_ds.GetRasterBand(i+1).SetNoDataValue(no_data_value)
        dst_ds.GetRasterBand(i+1).SetMetadata(metadata)
        if 'scale_factor' in metadata:
            scale = metadata['scale_factor']
            dst_ds.GetRasterBand(i+1).SetScale(float(scale))
        if 'add_offset' in metadata:
            offset = metadata['add_offset']
            dst_ds.GetRasterBand(i+1).SetOffset(float(offset))
    dst_ds.FlushCache()
    dst_ds = None


def subset_dem(raster_file, global_dem):
    auxdata_folder = os.path.dirname(raster_file).replace('IMG_DATA', 'AUX_DATA')
    out_raster_file = os.path.join(auxdata_folder, 'dem.tif')
    if not os.path.exists(os.path.dirname(out_raster_file)):
        os.makedirs(os.path.dirname(out_raster_file))
    subset_raster(raster_file, global_dem, out_raster_file, xRes=30, yRes=30)
    return out_raster_file

def subset_atmos(raster_file, atmo_file):
    auxdata_folder = os.path.dirname(raster_file).replace('IMG_DATA', 'AUX_DATA')
    out_raster_file = os.path.join(auxdata_folder, atmo_file.split('/')[-2], atmo_file.split('/')[-1])
    if not os.path.exists(os.path.dirname(out_raster_file)):
        os.makedirs(os.path.dirname(out_raster_file))
    subset_raster(raster_file, atmo_file, out_raster_file, xRes=27450, yRes=27450)
    return out_raster_file

# # out_raster_file = '/Users/fengyin/Downloads/2021_12_22_aod550.tif'
# global_dem = '/vsicurl/https://raw.githubusercontent.com/MarcYin/Copernicus_GLO_30_DEM_VRT/main/copernicus_GLO_30_dem.vrt'
# raster_file = '/Users/fengyin/S2B_MSIL1C_20220801T233659_N0400_R030_T01WCM_20220801T235506.SAFE/GRANULE/L1C_T01WCM_A028227_20220801T233653/IMG_DATA/T01WCM_20220801T233659_B02.jp2' 
# atmo_file = '/vsicurl/https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/2021_12_22/2021_12_22_aod550.tif'
# subset_atmos(raster_file, atmo_file)
# subset_dem(raster_file, global_dem)