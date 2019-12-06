import gdal
from osgeo import osr 
import numpy as np
import numpy.ma as ma
#gdalwarp -of VRT -t_srs "+proj=sinu" -te xmin ymin xmax ymax -ts 2400 2400 global_dem.vrt modis_cutoff.vrt

class reproject_data(object):
    '''
    A function uses a source and a target images and 
    and clip the source image to match the extend, 
    projection and resolution as the target image.

    '''
    def __init__(self, source_img,
                 target_img   = None,
                 dstSRS       = None,
                 srcNodata    = np.nan,
                 dstNodata    = np.nan,
                 outputType   = None,
                 verbose      = False,
                 xmin         = None,
                 xmax         = None,
                 ymin         = None, 
                 ymax         = None,
                 xRes         = None,
                 yRes         = None,
                 xSize        = None,
                 ySize        = None,
                 resample     = 1
                 ):

        self.source_img = source_img
        self.target_img = target_img
        self.verbose    = verbose
        self.dstSRS     = dstSRS
        self.srcNodata  = srcNodata
        self.dstNodata  = dstNodata
        self.outputType = gdal.GDT_Unknown if outputType is None else outputType
        self.xmin       = xmin
        self.xmax       = xmax
        self.ymin       = ymin
        self.ymax       = ymax
        self.xRes       = xRes
        self.yRes       = yRes
        self.xSize      = xSize
        self.ySize      = ySize
        self.resample   = resample
        if self.srcNodata is None:
            try:                           
                self.srcNodata = ' '.join([i.split("=")[1] for i in gdal.Info(self.source_img).split('\n') if' NoData' in i])
            except:                        
                self.srcNodata = None
        if (self.target_img is None) & (self.dstSRS is None):
            raise IOError('Projection should be specified ether from a file or a projection code.')
        elif self.target_img is not None:
            try:
                g     = gdal.Open(self.target_img)
            except:
                g     = target_img
            geo_t = g.GetGeoTransform()
            x_size, y_size = g.RasterXSize, g.RasterYSize     

            if self.xRes is None:
                self.xRes = abs(geo_t[1])
            if self.yRes is None:
                self.yRes = abs(geo_t[5])

            if self.xSize is not None: 
                x_size = 1. * self.xSize * self.xRes / abs(geo_t[1])
            if self.ySize is not None: 
                y_size = 1. * self.ySize * self.yRes / abs(geo_t[5])

            xmin, xmax = min(geo_t[0], geo_t[0] + x_size * geo_t[1]), \
                         max(geo_t[0], geo_t[0] + x_size * geo_t[1])  
            ymin, ymax = min(geo_t[3], geo_t[3] + y_size * geo_t[5]), \
                         max(geo_t[3], geo_t[3] + y_size * geo_t[5])
            dstSRS     = osr.SpatialReference( )
            raster_wkt = g.GetProjection()
            dstSRS.ImportFromWkt(raster_wkt)
            self.g = gdal.Warp('', self.source_img, format = 'MEM', outputBounds = [xmin, ymin, xmax, ymax], dstNodata=self.dstNodata, warpOptions = ['NUM_THREADS=ALL_CPUS'],\
                                xRes = self.xRes, yRes = self.yRes, dstSRS = dstSRS, outputType = self.outputType, srcNodata = self.srcNodata, resampleAlg = self.resample)
            
        else:
            self.g = gdal.Warp('', self.source_img, format = 'MEM', outputBounds = [self.xmin, self.ymin, \
                               self.xmax, self.ymax], xRes = self.xRes, yRes = self.yRes, dstSRS = self.dstSRS, warpOptions = ['NUM_THREADS=ALL_CPUS'],\
                               copyMetadata=True, outputType = self.outputType, dstNodata=self.dstNodata, srcNodata = self.srcNodata, resampleAlg = self.resample)
        if self.g.RasterCount <= 3:
            self.data = self.g.ReadAsArray()
            #return self.data
        elif self.verbose:
            print('There are %d bands in this file, use g.GetRasterBand(<band>) to avoid reading the whole file.'%self.g.RasterCount)

def array_to_raster(array, example_file):                                                                                                               
    if array.ndim == 2:                  
        bands = 1                        
    elif array.ndim ==3:                 
        bands = min(array.shape)     
        t = np.argsort(array.shape)
        array = array.transpose(t)
    else:                                
        raise IOError('Only 2 or 3 D array is supported.')
    try:                                 
        g = gdal.Open(example_file)      
    except:                              
        g = example_file                 
    driver = gdal.GetDriverByName('MEM') 
    ds = driver.Create('', array.shape[-1], array.shape[-2], bands, gdal.GDT_Float32)
    ds.SetProjection(g.GetProjection())  
    geotransform    = list(g.GetGeoTransform())
    geotransform[1] = geotransform[1] * g.RasterXSize / (1. * array.shape[-1])
    geotransform[5] = geotransform[5] * g.RasterYSize / (1. * array.shape[-2])
    ds.SetGeoTransform(geotransform)     
    if array.ndim == 3:                  
        for i in range(bands):           
            ds.GetRasterBand(i+1).WriteArray(array[i])
    else:                                
         ds.GetRasterBand(1).WriteArray(array)                                                                                                           
    return ds


if __name__=='__main__':
    ele = reproject_data('~/DATA/Multiply/eles/global_dem.vrt','/data/nemesis/S2_data/32/U/PU/2017/12/15/0/B02.jp2', outputType= gdal.GDT_Float32, ) 
    #ele.get_it()
    mask = (ele.data == -32768) | (~np.isfinite(ele.data))
    ele.data = ma.array(ele.data, mask = mask)
