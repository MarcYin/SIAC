#!/usr/bin/env python  
import os
import osr
import ogr
import gdal
import numpy as np
from pyproj import Proj, transform

'''
This is a function used for the calculation of MODIS
tile names from lat and lon coordinates. get_raster_hv 
is used for the calculation of MODIS hv from a raster 
file and get_vector_hv form a raster file

'''

x_step = -463.31271653   
y_step = 463.31271653    
#m_y0, m_x0 = -20015109.354, 10007554.677

tile_width = 1111950.5196666666
# m_x0, m_y0 = -20015109.35579742, -10007554.677898709
m_x0, m_y0 = -20015109.35579742, -10007554.677898709

def get_raster_hv(example_file):
    try:
        g = gdal.Open(example_file)
    except:
        try:
            g = example_file
            g.GetGeoTransform()[0]
        except:
            raise IOError('aoi has to be raster file or a gdal object.')
    geo_t = g.GetGeoTransform()
    x_size, y_size = g.RasterYSize, g.RasterXSize
 
    H_res_geo = osr.SpatialReference()
    H_res_geo.ImportFromWkt(g.GetProjection())
    modisProj= Proj('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs')
    wgs84Proj = Proj(init='epsg:4326')
    raster_wkt = Proj(H_res_geo.ExportToProj4())
    xs = [geo_t[0], geo_t[0] + geo_t[1]*x_size, geo_t[0], geo_t[0] + geo_t[1]*x_size]
    ys = [geo_t[3], geo_t[3] + geo_t[5]*y_size, geo_t[3] + geo_t[5]*y_size, geo_t[3]]

    lon, lat = transform(raster_wkt, wgs84Proj, np.array(xs).ravel(), np.array(ys).ravel())
    hh, vv = mtile_cal(lat, lon)
    tiles  = ['h%02dv%02d'%(hh[i], vv[i]) for i in range(len(hh))]
    tiles  = np.unique(tiles)
    return tiles 


def get_vector_hv(aoi):
    try:
        og = ogr.Open(aoi)
    except:
        try:
            og = aoi
            l = og.GetLayer(0)
        except:
            raise IOError('aoi has to be vector file or a ogr object')
    feature = og.GetLayer(0).GetFeature(0)
    coordinates = feature.geometry().GetGeometryRef(-0).GetPoints()
    gg = feature.GetGeometryRef()  
    sp = gg.GetSpatialReference()

    wgs84Proj = Proj(init='epsg:4326')
    vector_wkt = Proj(sp.ExportToProj4())
    xs, ys  = np.array(coordinates).T
    lon, lat = transform(vector_wkt, wgs84Proj, np.array(xs).ravel(), np.array(ys).ravel())
    hh, vv = mtile_cal(lat, lon)
    tiles  = ['h%02dv%02d'%(hh[i], vv[i]) for i in range(len(hh))]
    tiles  = np.unique(tiles)
    return tiles    

def mtile_cal(lat, lon):
    outProj= Proj('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs')
    inProj = Proj(init='epsg:4326')
    ho, vo = transform(inProj,outProj, np.array(lon).ravel(), np.array(lat).ravel())
    h = ((ho-m_x0)/tile_width).astype(int)
    v = 17 - ((vo-m_y0)/tile_width).astype(int)
    return h, v

if __name__ == '__main__':
    aoi          = '~/DATA/S2_MODIS/l_data/LC08_L1TP_014034_20170831_20170915_01_T1/AOI.json'
    example_file = '~/DATA/S2_MODIS/l_data/LC08_L1TP_014034_20170831_20170915_01_T1/aot.tif'
    print(get_vector_hv(aoi))
    print(get_vector_hv(ogr.Open(aoi)))
    print(get_raster_hv(example_file))
    print(get_raster_hv(gdal.Open(example_file)))

