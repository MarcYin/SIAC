import numpy as np
from osgeo import gdal, osr


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