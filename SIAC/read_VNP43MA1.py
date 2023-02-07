import h5py
import numpy as np
from osgeo import gdal, osr
from functools import partial
from concurrent.futures import ThreadPoolExecutor

def GetGeoTransform(fname):
    
    ds = gdal.OpenEx(fname, gdal.OF_MULTIDIM_RASTER)
    rootGroup = ds.GetRootGroup()
    info = rootGroup.OpenGroup('HDFEOS INFORMATION')
    StructMetadata = info.OpenMDArray('StructMetadata.0')
    StructMetadata = StructMetadata.Read()[0].split()

    
    ulc = [i for i in StructMetadata if 'UpperLeftPointMtrs' in i]
    ulcX, ulcY = eval(ulc[0].split('=')[1])
    
    lrc = [i for i in StructMetadata if 'LowerRightMtrs' in i]
    lrcX, lrcY = eval(lrc[0].split('=')[1])
    
    xRes = 926.6254330555555
    yRes = -926.6254330555555   
    # minX, minY, maxX, maxY = ulcX, ulcY + yRes * 1200, ulcX + xRes * 1200, ulcY
    # print(minX, minY, maxX, maxY)
    # 
    geotransform = (ulcX, xRes, 0, ulcY, 0, yRes)
    outputBounds = [ulcX, ulcY, lrcX, lrcY]
    bounds = (ulcX, lrcY, lrcX, ulcY)
    
    return geotransform


def get_bounds(aoi, toa, pix_res):
    g = gdal.Warp('', toa, format = 'MEM', cutlineDSName=aoi, cropToCutline=True, xRes = pix_res, yRes = pix_res)
    dstSRS = g.GetProjectionRef()
    geo_t = g.GetGeoTransform()
    x_size, y_size = g.RasterXSize, g.RasterYSize     
    xmin, xmax = min(geo_t[0], geo_t[0] + x_size * geo_t[1]), max(geo_t[0], geo_t[0] + x_size * geo_t[1])  
    ymin, ymax = min(geo_t[3], geo_t[3] + y_size * geo_t[5]), max(geo_t[3], geo_t[3] + y_size * geo_t[5])
    outputBounds = [xmin, ymin, xmax, ymax]
    return dstSRS, outputBounds

def clip_VNP43(fnames, bands, outputBounds, xRes, yRes, dstSRS, folder = None):
    
    # brdf_qa_template = 'HDF5:"%s"://HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data_Fields/BRDF_Albedo_Band_Mandatory_Quality_%s'
    # brdf_para_template = 'HDF5:"%s"://HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data_Fields/BRDF_Albedo_Parameters_%s'
    
    
    brdf_paras = []
    brdf_qas = []
    
    for fname in fnames:

        geotransform = GetGeoTransform(fname)

        input_file = h5py.File(fname, 'r')
        
        
        # BRDF in memory file
        brdf_paras_band_num = len(bands) * 3
        driver = gdal.GetDriverByName("MEM")
        dst_ds = driver.Create('', 1200, 1200, brdf_paras_band_num, gdal.GDT_Int16)
        
        dst_ds.SetGeoTransform(geotransform)
        
        spatialRef = osr.SpatialReference()
        spatialRef.ImportFromProj4('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs')
        dst_ds.SetProjection(spatialRef.ExportToWkt())
    
        # temp = brdf_para_template
        i = 1
        for band in bands:
            
            # band_fname = temp%(fname, band)
            # g = gdal.Open(band_fname)
            # data = g.ReadAsArray()
            
            raster_bands_fname = 'BRDF_Albedo_Parameters_%s'%band
            data = input_file['HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/%s'%raster_bands_fname][:]
            
            raster_band_name_prefix = raster_bands_fname
            brdf_parameters = ['fiso', 'fvol', 'fgeo']
            for j in range(3):
                raster_band = dst_ds.GetRasterBand(i)
                raster_band.SetNoDataValue(32767)
    
                raster_band_name = raster_band_name_prefix + '_' + brdf_parameters[j]
                raster_band.SetDescription(raster_band_name)
                raster_band.WriteArray(data[:, :, j])
                i += 1

        
        brdf_paras.append(dst_ds)
        
        brdf_qa_band_num = len(bands)
        driver = gdal.GetDriverByName("MEM")
        dst_ds = driver.Create('', 1200, 1200, brdf_qa_band_num, gdal.GDT_Byte)
        
        dst_ds.SetGeoTransform(geotransform)
        spatialRef = osr.SpatialReference()
        spatialRef.ImportFromProj4('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs')
        dst_ds.SetProjection(spatialRef.ExportToWkt())
        
        # temp = brdf_qa_template
        i = 1
        for band in bands:
            # band_fname = temp%(fname, band)
            # g = gdal.Open(band_fname)
            # data = g.ReadAsArray()
            
            raster_band_name = 'BRDF_Albedo_Band_Mandatory_Quality_%s'%band
            data = input_file['HDFEOS/GRIDS/VIIRS_Grid_BRDF/Data Fields/%s'%raster_band_name][:]
            
            # raster_band_name = band_fname.split('/')[-1]
            raster_band = dst_ds.GetRasterBand(i)
            raster_band.SetNoDataValue(255)
            # raster_band_name = raster_band_name + '_' + brdf_parameters[j]
            raster_band.SetDescription(raster_band_name)
            raster_band.WriteArray(data)
            i += 1
        brdf_qas.append(dst_ds)
        
        
    brdf_paras_vrt = gdal.BuildVRT('', brdf_paras)
    brdf_qas_vrt = gdal.BuildVRT('', brdf_qas)
    
    modis_sinu = osr.SpatialReference() 
    sinu = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    modis_sinu.ImportFromProj4 (sinu)

    year_doy = fname.split('/')[-1].split('.')[1][1:]
    
    para_filename = '_'.join(['VNP43MA1', year_doy] + bands + ['paras.tif'])
    qa_filename = '_'.join(['VNP43MA1', year_doy] + bands + ['qas.tif'])
    
    # reader = lambda inps: gdal.Warp(inps[0], inps[1], 
    #                   format='GTiff', 
    #                   xRes=xRes, 
    #                   yRes=yRes, 
    #                   dstSRS = dstSRS, 
    #                   resampleAlg = 0, 
    #                   outputBounds = outputBounds, 
    #                   srcSRS=sinu,
    #                   )
    # brdf_ds = reader([para_filename, brdf_paras_vrt])
    # qa_ds = reader([qa_filename, brdf_qas_vrt])
    # brdf_ds = None
    # qa_ds = None

    reader = lambda vrt: gdal.Warp('', vrt, 
                      format='MEM', 
                      xRes=xRes, 
                      yRes=yRes, 
                      dstSRS = dstSRS, 
                      resampleAlg = 0, 
                      outputBounds = outputBounds, 
                      srcSRS=sinu,
                      )
    brdf_ds = reader(brdf_paras_vrt)
    qa_ds = reader(brdf_qas_vrt)
    
    geoTransoformation = qa_ds.GetGeoTransform()
    
    dst_ds = None
    brdf_paras_vrt = None
    brdf_qas_vrt = None

    brdf_data = brdf_ds.ReadAsArray()
    qa_data = qa_ds.ReadAsArray()

    # print(geoTransoformation)
    
    return brdf_data, qa_data, geoTransoformation


def read_VNP43MA1(fnames_dates, bands, outputBounds, xRes,  yRes, dstSRS):
    
    par = partial(clip_VNP43, bands = bands, outputBounds = outputBounds, xRes = xRes,  yRes = yRes, dstSRS = dstSRS)
    
    all_fnames = [fname_date[0] for fname_date in fnames_dates]

    with ThreadPoolExecutor(17) as executor:
        ret = executor.map(par, all_fnames)
    ret = list(ret)
    
    brdf_dats = np.array([i[0] for i in ret])
    qa_dats   = np.array([i[1] for i in ret])
    geotransform = ret[0][2]
    
    driver = gdal.GetDriverByName("MEM")
    mg = driver.Create('', brdf_dats.shape[-2], brdf_dats.shape[-1], 1, gdal.GDT_Byte)
    mg.SetGeoTransform(geotransform)
    mg.SetProjection(dstSRS)
    
    new_shape = (brdf_dats.shape[0], len(bands), 3, brdf_dats.shape[-2], brdf_dats.shape[-1])
    brdf_dats = brdf_dats.reshape(new_shape)
    ws  = 0.618034**qa_dats       
    ws[qa_dats==255]    = 0 
    brdf_dats[brdf_dats==32767] = 0
    return brdf_dats, ws, mg
