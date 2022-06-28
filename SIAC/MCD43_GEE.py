import os
import requests
import datetime
from io import BytesIO
from zipfile import ZipFile
from functools import partial
from concurrent.futures import ThreadPoolExecutor
from osgeo import gdal, osr
import numpy as np
from retry import retry

# Stop GDAL printing both warnings and errors to STDERR
gdal.PushErrorHandler('CPLQuietErrorHandler')

# Make GDAL raise python exceptions for errors (warnings won't raise an exception)
gdal.UseExceptions()


import ee
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')

@retry(tries=10, delay=1, backoff=2)
def downloader(band, geometry, date_start, date_end, folder = './'):
    '''
    band: 
        BRDF_Albedo_Parameters_Band[1-7]_[iso, geo, vol]
        BRDF_Albedo_Band_Mandatory_Quality_Band[1-7]
    geometry:
        ee.Geometry
    date_start:
        'YYYY-MM-DD', e.g. '2018-01-01'
    date_end:
        'YYYY-MM-DD', e.g. '2018-10-10'
    '''
    paras = {
             'name': 'test', 
             'crs': 'SR-ORG:6974',
             'filePerBand': False,
             'scale': 463.312716528,
             'region': geometry.buffer(1000)
             }
    mcd43a1 = ee.ImageCollection("MODIS/006/MCD43A1").filterDate(date_start, date_end)
    band_time_series = mcd43a1.select([band])
    img = band_time_series.toBands()
    header = 'MCD43_'
    
    paras['name'] =  header + band
    url = img.getDownloadURL(paras)
    resp = requests.get(url)
    if resp.status_code != 200:
        resp.raise_for_status()
    zipfile = ZipFile(BytesIO(resp.content))
    zip_names = zipfile.namelist()
    # print(zip_names)
    zipfile.extractall(folder)
    return folder + '/' + header + band + '.tif'


def get_MCD43_GEE(obs_time, coords, temporal_window, folder):
    geometry = ee.Geometry.Polygon(coords[0])

    first_half_folder = folder + '/first_half/'
    if not os.path.exists(first_half_folder):
        os.mkdir(first_half_folder)

    date_start  = (obs_time - datetime.timedelta(days = int(temporal_window))).strftime('%Y-%m-%d') 
    date_end    = obs_time.strftime('%Y-%m-%d')
    par = partial(downloader, geometry=geometry, date_start=date_start, date_end=date_end, folder = first_half_folder)
    
    bands = ['BRDF_Albedo_Parameters_Band%d_%s'%(band, ker) for band in [1, 2, 3, 4, 6, 7] for ker in ['iso', 'vol', 'geo']]
    qa_bands = ['BRDF_Albedo_Band_Mandatory_Quality_Band%d'%(band) for band in [1, 2, 3, 4, 6, 7]]
    all_bands = bands + qa_bands

    with ThreadPoolExecutor() as executor:
        result = executor.map(par, all_bands)


    second_half_folder = folder + '/second_half/'
    if not os.path.exists(second_half_folder):
        os.mkdir(second_half_folder)

    date_start  = obs_time.strftime('%Y-%m-%d') 
    date_end    = (obs_time + datetime.timedelta(days = int(temporal_window) + 1)).strftime('%Y-%m-%d')
    par = partial(downloader, geometry=geometry, date_start=date_start, date_end=date_end, folder = second_half_folder)
    
    bands = ['BRDF_Albedo_Parameters_Band%d_%s'%(band, ker) for band in [1, 2, 3, 4, 6, 7] for ker in ['iso', 'vol', 'geo']]
    qa_bands = ['BRDF_Albedo_Band_Mandatory_Quality_Band%d'%(band) for band in [1, 2, 3, 4, 6, 7]]
    all_bands = bands + qa_bands

    with ThreadPoolExecutor() as executor:
        result = executor.map(par, all_bands)
    
    for band in all_bands:
        first_half_name = os.path.join(first_half_folder, 'MCD43_' + band + '.tif')
        second_half_name = os.path.join(first_half_folder, 'MCD43_' + band + '.tif')
        fname = os.path.join(folder, 'MCD43_' + band + '.tif')
    
        g = gdal.Open(first_half_name)
        gg = g.GetRasterBand(1)
        datatype = gg.DataType
        nodatavalue = gg.GetNoDataValue() 

        geotransform = g.GetGeoTransform()
        srs = g.GetProjectionRef()

        first_half_data  = g.ReadAsArray()
        second_half_data = gdal.Open(second_half_name).ReadAsArray()
        rows, cols = second_half_data.shape[1:]
        nbands = len(first_half_data) + len(second_half_data)
        data = np.concatenate([first_half_data, second_half_data], axis=0)

        driver = gdal.GetDriverByName('GTiff')
        raster = driver.Create(fname, cols, rows, nbands, datatype)
        raster.SetGeoTransform(geotransform)
        for i in range(nbands):
            raster_band = raster.GetRasterBand(i+1)
            raster_band.WriteArray(data[i])
            raster_band.SetNoDataValue(nodatavalue)
        raster.SetProjection(srs)
        raster.FlushCache()
        


def read_MCD3_GEE(folder, boa_bands, outputBounds, xRes, yRes, dstSRS):

    modis_sinu = osr.SpatialReference() 
    sinu = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    modis_sinu.ImportFromProj4 (sinu)
    
    reader = lambda fname: gdal.Warp('', fname, 
                          format='MEM', 
                          xRes=xRes, 
                          yRes=yRes, 
                          dstSRS = dstSRS, 
                          resampleAlg = 0, 
                          outputType = gdal.GDT_Float32,
                          outputBounds = outputBounds, 
                          srcSRS=sinu,
                          warpOptions = ['NUM_THREADS=ALL_CPUS'])

    data_bands = []
    qa_bands = []
    for band in boa_bands:
        band_name = 'BRDF_Albedo_Band_Mandatory_Quality_Band%d'%(band)
        fname = folder + '/' + 'MCD43_' + band_name + '.tif'
        qa_bands.append(reader(fname).ReadAsArray())

        for ker in ['iso', 'vol', 'geo']:
            band_name = 'BRDF_Albedo_Parameters_Band%d_%s'%(band, ker)
            fname = folder + '/' + 'MCD43_' + band_name + '.tif'
            data_bands.append(reader(fname).ReadAsArray())
    
    mg = reader(fname)
    data_bands = np.array(data_bands)
    qa_bands = np.array(qa_bands)

    new_shape = (len(boa_bands), 3) + data_bands.shape[1:]
    data_bands = data_bands.reshape(new_shape)

    data_bands = data_bands.transpose(2, 0, 1, 3, 4)
    qa_bands   = qa_bands.transpose(1, 0, 2, 3)
    ws  = 0.618034**qa_bands       
    ws[qa_bands==255]    = 0 
    data_bands[data_bands==32767] = 0

    return data_bands, ws, mg

def get_bounds(aoi, toa, pix_res):
    g = gdal.Warp('', toa, format = 'MEM', cutlineDSName=aoi, cropToCutline=True, xRes = pix_res, yRes = pix_res)
    dstSRS = g.GetProjectionRef()
    geo_t = g.GetGeoTransform()
    x_size, y_size = g.RasterXSize, g.RasterYSize     
    xmin, xmax = min(geo_t[0], geo_t[0] + x_size * geo_t[1]), max(geo_t[0], geo_t[0] + x_size * geo_t[1])  
    ymin, ymax = min(geo_t[3], geo_t[3] + y_size * geo_t[5]), max(geo_t[3], geo_t[3] + y_size * geo_t[5])
    outputBounds = [xmin, ymin, xmax, ymax]
    return dstSRS, outputBounds


# from osgeo import gdal

# fname = '/Users/fengyin/Downloads/LC09_L1TP_025030_20220514_20220514_02_T1/MCD43/first_half/MCD43_BRDF_Albedo_Parameters_Band1_iso.tif'
# fname2 = fname.replace('first_half', 'second_half')


# g1 = gdal.Open(fname)
# x_size, y_size = g1.RasterXSize, g1.RasterYSize
# g2 = gdal.Open(fname2)

# nbands1 = g1.RasterCount
# nbands2 = g2.RasterCount

# nbands = nbands1 + nbands2


# geotransform = g.GetGeoTransform()
# srs = g.GetProjectionRef()
    
# drv = gdal.GetDriverByName("VRT")
# vrt = drv.Create("test.vrt", x_size, y_size, nbands, g1.GetRasterBand(1).DataType)
# vrt.SetGeoTransform(geotransform)
# vrt.SetProjection(srs)

# source_files = [fname] * nbands1 + [fname2] * nbands2

# x_offset, y_offset, x_source_size, y_source_size = 0, 0, x_size, y_size
# dest_x_offset, dest_y_offset, x_dest_size, y_dest_size = 0, 0, x_size, y_size


# for band_ind in range(nbands1):
#     band = g1.GetRasterBand(band_ind+1)
#     x_block, y_block = band.GetBlockSize()
    
    
#     new_band = vrt.GetRasterBand(band_ind+1)

#     source_path = fname
#     source_band = band_ind + 1
    
#     simple_source = '<SourceFilename relativeToVRT="1">%s</SourceFilename>' % source_path + \
#         '<SourceBand>%i</SourceBand>' % source_band + \
#         '<SourceProperties RasterXSize="%i" RasterYSize="%i" DataType="Real" BlockXSize="%i" BlockYSize="%i"/>' % (x_size, y_size, x_block, y_block) + \
#         '<SrcRect xOff="%i" yOff="%i" xSize="%i" ySize="%i"/>' % (x_offset, y_offset, x_source_size, y_source_size) + \
#         '<DstRect xOff="%i" yOff="%i" xSize="%i" ySize="%i"/>' % (dest_x_offset, dest_y_offset, x_dest_size, y_dest_size)
#     new_band.SetMetadataItem("ComplexSource", simple_source)
#     new_band.SetNoDataValue(band.GetNoDataValue())
    
# for band_ind in range(nbands2):
#     band = g2.GetRasterBand(band_ind+1)
#     x_block, y_block = band.GetBlockSize()
    
#     new_band = vrt.GetRasterBand(band_ind+1 + nbands1)

#     source_path = fname
#     source_band = band_ind + 1
    
#     simple_source = '<SourceFilename relativeToVRT="1">%s</SourceFilename>' % source_path + \
#         '<SourceBand>%i</SourceBand>' % source_band + \
#         '<SourceProperties RasterXSize="%i" RasterYSize="%i" DataType="Real" BlockXSize="%i" BlockYSize="%i"/>' % (x_size, y_size, x_block, y_block) + \
#         '<SrcRect xOff="%i" yOff="%i" xSize="%i" ySize="%i"/>' % (x_offset, y_offset, x_source_size, y_source_size) + \
#         '<DstRect xOff="%i" yOff="%i" xSize="%i" ySize="%i"/>' % (dest_x_offset, dest_y_offset, x_dest_size, y_dest_size)
#     new_band.SetMetadataItem("ComplexSource", simple_source)
#     new_band.SetNoDataValue(band.GetNoDataValue())
    
# vrt.FlushCache()







# # vrt_template =  '''<VRTDataset rasterXSize="1048" rasterYSize="521">
# #   <SRS dataAxisToSRSAxisMapping="1,2">PROJCS["MODIS Sinusoidal",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH]]</SRS>
# #   <GeoTransform> -7.6205675614525443e+06,  4.6331271652800001e+02,  0.0000000000000000e+00,  4.9213076749604158e+06,  0.0000000000000000e+00, -4.6331271652800001e+02</GeoTransform>
# # </VRTDataset>'''

# vrt_fname = fname.replace('/first_half', '').replace('.tif', '.vrt')
# g = gdal.BuildVRT(vrt_fname, fname)
# g.FlushCache()


# with open(vrt_fname, 'r') as f:
#     txt = f.read()
 
    