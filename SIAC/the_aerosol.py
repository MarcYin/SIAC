from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import gc
import sys
import SIAC.kernels as kernels 
import logging
import platform
import warnings
warnings.filterwarnings("ignore") 
import subprocess
import numpy as np
from glob import glob
try:     
    import cPickle as pkl
except:  
    import pickle as pkl
from functools import partial
from osgeo import ogr, osr, gdal
from SIAC.smoothn import smoothn
from scipy.fftpack import dct, idct
from SIAC.multi_process import parmap
from scipy.interpolate import griddata
from datetime import datetime, timedelta
from SIAC.psf_optimize import psf_optimize
from SIAC.create_logger import create_logger
from SIAC.raster_boundary import get_boundary
from SIAC.atmo_solver import solving_atmo_paras
from sklearn.linear_model import HuberRegressor 
from SIAC.reproject import reproject_data, array_to_raster
from scipy.ndimage import binary_dilation, binary_erosion

class solve_aerosol(object):
    '''
    Prepareing modis data to be able to pass into 
    atmo_cor for the retrieval of atmospheric parameters.
    '''
    def __init__(self,
                 sensor_sat,
                 toa_bands, 
                 band_wv,
                 band_index,
                 view_angles, 
                 sun_angles,
                 obs_time,
                 cloud_mask,
                 gamma       = 10.,
                 a_z_order   = 1,
                 pixel_res   = None,
                 aoi         = None,
                 aot_prior   = None,
                 wv_prior    = None,
                 o3_prior    = None,
                 aot_unc     = None,
                 wv_unc      = None,
                 o3_unc      = None, 
                 log_file    = None,
                 ref_scale   = 0.0001,
                 ref_off     = 0,
                 ang_scale   = 0.01,
                 ele_scale   = 0.001,
                 prior_scale = [1., 0.1, 46.698, 1., 1., 1.],
                 emus_dir    = 'SIAC/emus/',
                 mcd43_dir   = '~/DATA/Multiply/MCD43/',
                 global_dem  = '/vsicurl/http://www2.geog.ucl.ac.uk/~ucfafyi/eles/global_dem.vrt',
                 cams_dir    = '/vsicurl/http://www2.geog.ucl.ac.uk/~ucfafyi/cams/',
                 spec_m_dir  = 'SIAC/spectral_mapping',
                 aero_res    = 1000):                                 
        self.sensor      = sensor_sat[0]
        self.satellite   = sensor_sat[1]
        self.toa_bands   = toa_bands
        self.band_wv     = band_wv
        self.band_index  = band_index
        self.view_angles = view_angles
        self.sun_angles  = sun_angles 
        self.obs_time    = obs_time
        self.cloud_mask  = cloud_mask
        self.gamma       = gamma
        self.a_z_order   = a_z_order 
        self.pixel_res   = pixel_res
        self.aoi         = aoi
        self.aot_prior   = aot_prior
        self.wv_prior    = wv_prior
        self.o3_prior    = o3_prior
        self.aot_unc     = aot_unc
        self.wv_unc      = wv_unc
        self.o3_unc      = o3_unc
        self.log_file    = log_file
        self.ref_scale   = ref_scale
        self.ref_off     = ref_off
        self.ang_scale   = ang_scale
        self.ele_scale   = ele_scale
        self.prior_scale = prior_scale
        self.mcd43_dir   = mcd43_dir
        self.emus_dir    = emus_dir
        self.global_dem  = global_dem
        self.cams_dir    = cams_dir
        self.boa_wv      = [645, 859, 469, 555, 1640, 2130]
        self.aero_res    = aero_res
        self.mcd43_tmp   = '%s/MCD43A1.A%d%03d.%s.006.*.hdf'
        self.toa_dir     =  os.path.abspath('/'.join(toa_bands[0].split('/')[:-1]))
        try:
            spec_map     = np.loadtxt(spec_m_dir + '/TERRA_%s_spectral_mapping.txt'%self.satellite).T
        except:
            spec_map     = np.loadtxt(spec_m_dir + '/TERRA_%s_spectral_mapping.txt'%self.sensor).T
        self.spec_slope  = spec_map[0]
        self.spec_off    = spec_map[1]
        self.logger      = create_logger(self.log_file)

    def _create_base_map(self,):
        '''
        Deal with different types way to define the AOI, if none is specified, then the image bound is used.
        '''
        gdal.UseExceptions()
        ogr.UseExceptions()
        if self.aoi is not None:
            if os.path.exists(self.aoi):
                try:
                    g = gdal.Open(self.aoi)
                    #subprocess.call(['gdaltindex', '-f', 'GeoJSON', '-t_srs', 'EPSG:4326', self.toa_dir + '/AOI.json', self.aoi])
                    geojson = get_boundary(self.aoi)[0]
                    with open(self.toa_dir + '/AOI.json', 'wb') as f:
                        f.write(geojson.encode())
                except:
                    try:
                        gr = ogr.Open(str(self.aoi))
                        l = gr.GetLayer(0)
                        f = l.GetFeature(0)
                        g = f.GetGeometryRef()
                    except:
                        raise IOError('AOI file cannot be opened by gdal, please check it or transform into format can be opened by gdal')
            else:
                try:
                    g = ogr.CreateGeometryFromJson(self.aoi)
                except:
                    try:
                        g = ogr.CreateGeometryFromGML(self.aoi)
                    except:
                        try:
                            g = ogr.CreateGeometryFromWkt(self.aoi)
                        except:
                            try:
                                g = ogr.CreateGeometryFromWkb(self.aoi)
                            except:
                                raise IOError('The AOI has to be one of GeoJSON, GML, Wkt or Wkb.')
            gjson_str = '''{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},"geometry":%s}]}'''% g.ExportToJson()
            with open(self.toa_dir + '/AOI.json', 'wb') as f:
                f.write(gjson_str.encode())

        ogr.DontUseExceptions()
        gdal.DontUseExceptions()
        if not os.path.exists(self.toa_dir + '/AOI.json'):
            #g = gdal.Open(self.toa_bands[0])    
            #proj = g.GetProjection()            
            #if 'WGS 84' in proj:                
            #    subprocess.call(['gdaltindex', '-f', 'GeoJSON', self.toa_dir +'/AOI.json', self.toa_bands[0]])
            #else:                               
            #    subprocess.call(['gdaltindex', '-f', 'GeoJSON', '-t_srs', 'EPSG:4326', self.toa_dir +'/AOI.json', self.toa_bands[0]])
            geojson = get_boundary(self.toa_bands[0])[0]
            with open(self.toa_dir + '/AOI.json', 'wb') as f:                   
                f.write(geojson.encode())

            self.logger.warning('AOI is not created and full band extend is used')
            self.aoi = self.toa_dir + '/AOI.json'
        else:
            self.aoi = self.toa_dir + '/AOI.json'

        if self.pixel_res is None:
            self.pixel_res = abs(gdal.Open(self.toa_bands[0]).GetGeoTransform()[1])

        self.psf_xstd = 260 / self.pixel_res
        self.psf_ystd = 340 / self.pixel_res

    def _get_psf(self,):
        '''
        Get the PSF parameters
        '''
        toa     = self._toa_bands[-1].ReadAsArray()*self.ref_scale+self.ref_off
        index   = [self.hx, self.hy]
        boa     = np.ma.array(self.boa[-1])
        boa_unc = self.boa_unc[-1]
        mask    = self.bad_pix
        thresh   = 0.1
        ang     = 0
        psf     = psf_optimize(toa, index, boa, boa_unc, mask,thresh, xstd=self.psf_xstd, ystd= self.psf_ystd)
        xs, ys  = psf.fire_shift_optimize()
        self.logger.info('Solved PSF: %.02f, %.02f, %d, %d, %d, and R value is: %.03f.' \
                              %(self.psf_xstd, self.psf_ystd, 0, xs, ys, 1-psf.costs.min()))  
        shifted_mask = np.logical_and.reduce(((self.hx+int(xs)>=0),
                                              (self.hx+int(xs)<self.full_res[0]),
                                              (self.hy+int(ys)>=0),
                                              (self.hy+int(ys)<self.full_res[1])))
                               
        self.hx, self.hy = self.hx[shifted_mask]+int(xs), self.hy[shifted_mask]+int(ys)
        self.boa         = self.boa    [:, shifted_mask]
        self.boa_unc     = self.boa_unc[:, shifted_mask]

    def _mask_bad_pix(self):
        #snow_mask = blue > 0.6
        if os.path.exists(str(self.cloud_mask)):  
            cloud_g = gdal.Open(str(self.cloud_mask)) 
        elif type(self.cloud_mask).__module__== 'numpy':
            cloud_g = array_to_raster(self.cloud_mask, self.example_file)
        cloud = reproject_data(cloud_g, \
                                self.example_file, \
                                xRes=self.pixel_res, \
                                yRes=self.pixel_res, \
                                xSize=self.full_res[1], \
                                ySize=self.full_res[0], \
                                srcNodata = np.nan,\
                                outputType= gdal.GDT_Float32,\
                                resample = 0).data
        cloud  = cloud.astype(bool)
        RED    = None
        BLUE   = None
        SWIR_1 = None
        NIR    = None
        GREEN  = None

        if 3 in self.boa_bands:
            BLUE   = self._toa_bands[np.where(self.boa_bands==3)[0][0]].ReadAsArray()*self.ref_scale+self.ref_off
        if BLUE is not None:
            water_mask = BLUE < 0.05
            snow_mask  = BLUE > 0.6
            #del BLUE
        else:
            self.logger.error('Blue band is needed for the retirval of aerosol.')
        
        if 2 in self.boa_bands:
            NIR    = self._toa_bands[np.where(self.boa_bands==2)[0][0]].ReadAsArray()*self.ref_scale+self.ref_off
        if 1 in self.boa_bands:
            RED    = self._toa_bands[np.where(self.boa_bands==1)[0][0]].ReadAsArray()*self.ref_scale+self.ref_off
        if (RED is not None) & (NIR is not None):
            NDVI = (NIR - RED) / (NIR + RED)
            water_mask = ((NDVI < 0.01) & (NIR < 0.11)) | ((NDVI < 0.1) & (NIR < 0.05)) \
                        | (NIR <= 0.0001) | (RED <= 0.0001) | np.isnan(NIR) | np.isnan(RED)
            del NIR; del RED; del NDVI
        elif NIR is not None:
            water_mask = NIR < 0.05
            del NIR
        elif RED is not None:
            water_mask = RED < 0.05
            del RED

        if 6 in self.boa_bands:
            SWIR_1 = self._toa_bands[np.where(self.boa_bands==6)[0][0]].ReadAsArray()*self.ref_scale+self.ref_off
        if 1 in self.boa_bands:
            GREEN  = self._toa_bands[np.where(self.boa_bands==4)[0][0]].ReadAsArray()*self.ref_scale+self.ref_off
        if (SWIR_1 is not None) & (GREEN is not None):
            NDSI      = (GREEN - SWIR_1) / (SWIR_1 + GREEN)
            snow_mask = (NDSI > 0.42) | (SWIR_1 <= 0.0001) | (GREEN <= 0.0001) | np.isnan(SWIR_1) | np.isnan(GREEN) 
            del SWIR_1; del GREEN; del NDSI
        mask     = water_mask | snow_mask | cloud | (BLUE > 1)
        ker_size = int(2 * 1.96 * self.psf_ystd)
        mask     = binary_erosion (mask, structure = np.ones((3,3)).astype(bool), iterations=5).astype(bool)
        #mask     = binary_dilation(mask, structure = np.ones((3,3)).astype(bool), iterations=30).astype(bool)
        mask = binary_dilation(mask, structure = np.ones((3,3)).astype(bool), iterations=30+ker_size).astype(bool)
        mask[:30,:]  = mask[:,:30] = mask[:,-30:] = mask[-30:,:] =  True
        self.bad_pix = mask

    def _create_band_gs(self,):
        '''
        Create a lost of boa gs and cut them with AOI.
        '''
        self._toa_bands = []
        for band in self.toa_bands:
            g = gdal.Warp('', str(band), format = 'MEM', xRes = self.pixel_res, yRes = self.pixel_res, warpOptions = ['NUM_THREADS=ALL_CPUS'],\
                          srcNodata = 0, dstNodata=0, cutlineDSName= self.aoi, cropToCutline=True, resampleAlg = 0)
            self._toa_bands.append(g)
        self.example_file = self._toa_bands[0]
        x_max_pix         = self.example_file.RasterXSize * self.pixel_res
        y_max_pix         = self.example_file.RasterYSize * self.pixel_res
        self.xSize = int(np.ceil(x_max_pix/ (1. * self.aero_res))) 
        self.ySize = int(np.ceil(y_max_pix/ (1. * self.aero_res)))
        self.full_res   = self.example_file.RasterYSize, self.example_file.RasterXSize


    def _resamplers(self,):
        '''
        Define the sresamplers used for resampling different gdal files or objects to required resolution and size.
        '''
        self.nearest_resampler = lambda fname: reproject_data(fname, \
                                                                self.example_file, \
                                                                xRes=self.aero_res*0.5, \
                                                                yRes=self.aero_res*0.5, \
                                                                srcNodata = None,\
                                                                outputType= gdal.GDT_Float32, \
                                                                resample = 0 ).g # GRIORA_NearestNeighbour due to version changes...

        self.bilinear_resampler = lambda fname: reproject_data(fname, \
                                                                self.example_file, \
                                                                xRes=self.aero_res*0.5, \
                                                                yRes=self.aero_res*0.5, \
                                                                srcNodata = 0,\
                                                                outputType= gdal.GDT_Float32,\
                                                                resample = 1 ).g # GRIORA_Bilinear

    def _var_parser(self, var):
        if os.path.exists(str(var)):    
            var_g = gdal.Open(str(var)) 
        elif type(var).__module__== 'numpy':
            var_g = array_to_raster(var, self.example_file) 
        elif ('/vsicurl/' in str(var)) or ('/vsizip/') in str(var):                                                    
            var_g = gdal.Open(str(var))
        else:              
            var = float(var) 
            var_array = np.zeros((10,10))      
            var_array[:] = var
            var_g = array_to_raster(var_array, self.example_file)
        return var_g

    def _parse_angles(self,):
        '''
        Parsing angles
        '''
        self._view_angles = []
        if len(self.view_angles)==1:
            ang_g = self._var_parser(self.view_angles[0])
            self.view_angles = [ang_g,]
        for i in self.view_angles:
            ang_g = self._var_parser(i)
            self._view_angles.append(ang_g)

        self._sun_angles = []            
        if os.path.exists(str(self.sun_angles)):
            ang_g = self._var_parser(self.sun_angles)
            self._sun_angles = [ang_g,]#self.nearest_resampler(ang_g).ReadAsArray()
        else:
            for i in self.sun_angles:        
                ang_g = self._var_parser(i)
                self._sun_angles.append(ang_g)#self.nearest_resampler(ang_g).ReadAsArray())

      
    def _read_aux(self,):
        '''
        Create a list of AUX gdal objects, like priors and DEM, and cutted with AOI.
        '''
        self._ele   = self.bilinear_resampler(self.global_dem).ReadAsArray() * self.ele_scale
        
        time_ind    = np.abs((self.obs_time.hour  + self.obs_time.minute/60. + \
                              self.obs_time.second/3600.) - np.arange(0,25,3)).argmin()

        use_cams    = [False, False, False, False, False, False]
        priors      = [self.aot_prior, self.wv_prior, self.o3_prior, self.aot_unc, self.wv_unc, self.o3_unc]
        cams_names  = ['aod550', 'tcwv', 'gtco3'] 
        defalt_uncs = [0.4, 0.1, 0.05]
        for i in range(3):
            if priors[i] is None:
                priors[i] = self.cams_dir + '/'.join([datetime.strftime(self.obs_time, '%Y_%m_%d'),\
                                                       datetime.strftime(self.obs_time, '%Y_%m_%d')+'_%s.tif'%cams_names[i]])
                use_cams[i] = True
                if priors[i+3] is None:
                    priors[i+3] = defalt_uncs[i]
        temp = []
        for _, i in enumerate(priors):
            var_g = self._var_parser(i) 
            prior_g = self.bilinear_resampler(var_g)
            if use_cams[_]:
                g      = var_g.GetRasterBand(int(time_ind+1))
                offset = g.GetOffset()            
                scale  = g.GetScale()             
                data   = prior_g.GetRasterBand(int(time_ind+1)).ReadAsArray() * scale + offset
            else:
                data   = prior_g.ReadAsArray()
            temp.append(data * self.prior_scale[_])

        self._annoying_angles(prior_g)
        self._aot, self._tcwv, self._tco3, self._aot_unc, self._tcwv_unc, self._tco3_unc = temp
  
    def _find_boa_bands(self,):
        '''
        Find the closest MODIS bands to the TOA bands based on the Central wavelength of each band, 
        also reject bands are far from (more than 150nm) MODIS bands.
        '''
        self.band_wv     = np.array(self.band_wv)
        self.boa_wv      = np.array(self.boa_wv)
        self.toa_bands   = np.array(self.toa_bands)
        self.view_angles = np.array(self.view_angles)        

        if len(self.view_angles) == len(self.toa_bands):
            sa_va_seperate = False
        elif len(self.view_angles) == 2*len(self.toa_bands):
            sa_va_seperate = True

        mask = np.any(abs(self.band_wv[...,None] - self.boa_wv) <150, axis=1)
        band_index       = np.argmin(abs(self.band_wv[...,None] - self.boa_wv), axis=1)
        band_index       = band_index[mask]
        self.boa_bands   = np.array([1, 2, 3, 4, 6, 7])[band_index,]
        self.boa_wv      = self.boa_wv[band_index,]
        self.toa_bands   = self.toa_bands[mask,]

        order            = np.argsort(self.band_wv[mask,])
        self.boa_bands   = self.boa_bands[order,]
        self.boa_wv      = self.boa_wv[order,]
        self.toa_bands   = self.toa_bands[order,]
        if sa_va_seperate is True:
            self.view_angles = self.view_angles.reshape(2, -1)[:,mask]
            self.view_angles = self.view_angles[:,order].ravel()
        elif sa_va_seperate is False:
            self.view_angles = self.view_angles[mask,]
            self.view_angles = self.view_angles[order,]
        
        
    def _get_bounds(self,):
        geo_t = self.example_file.GetGeoTransform()
        raster_wkt =  self.example_file.GetProjection()
        x_size, y_size = self.example_file.RasterXSize, self.example_file.RasterYSize
        xmin, xmax = min(geo_t[0], geo_t[0] + x_size * geo_t[1]), \
                     max(geo_t[0], geo_t[0] + x_size * geo_t[1])  
        ymin, ymax = min(geo_t[3], geo_t[3] + y_size * geo_t[5]), \
                     max(geo_t[3], geo_t[3] + y_size * geo_t[5])  
        bounds = [xmin, ymin, xmax, ymax]
        return bounds, raster_wkt
    
    def _read_MCD43(self,fnames):
        def warp_data(fname, aoi,  xRes, yRes):
            g = gdal.Warp('',fname, format = 'MEM', srcNodata = 32767, dstNodata=0, \
                          cutlineDSName=aoi, xRes = xRes, yRes = yRes, cropToCutline=True, resampleAlg = 0) # weird adaptation for gdal 2.3, this should be a bug in gdal 2.3
            return g.ReadAsArray()                                                                          # no reason you have to specify the srcNodata to use dstNodata
        par = partial(warp_data, aoi = self.aoi, xRes = self.aero_res*0.5, yRes = self.aero_res*0.5) 
        n_files = int(len(fnames)/2)
        ret = parmap(par, fnames) 
        das = np.array(ret[:n_files])     
        qas = np.array(ret[n_files:]) 
        ws  = 0.618034**qas       
        ws[qas==255]     = 0 
        das[das==32767] = 0
        return das, ws

    def _get_boa(self, temporal_filling = 16):
        
        qa_temp = 'MCD43_%s_BRDF_Albedo_Band_Mandatory_Quality_Band%d.vrt'
        da_temp = 'MCD43_%s_BRDF_Albedo_Parameters_Band%d.vrt'
        doy = self.obs_time.strftime('%Y%j')
        if temporal_filling == True:
            temporal_filling = 16
        if temporal_filling:
            days   = [(self.obs_time - timedelta(days = int(i))) for i in np.arange(temporal_filling, 0, -1)] + \
                     [(self.obs_time + timedelta(days = int(i))) for i in np.arange(0, temporal_filling+1,  1)]
            fnames = []
            for temp in [da_temp, qa_temp]:                                                                                    
                for day in days:
                    for band in self.boa_bands:
                        fname = self.mcd43_dir + '/'.join([day.strftime('%Y-%m-%d'), temp%(day.strftime('%Y%j'), band)])
                        fnames.append(fname) 
        else:
            fnames = []
            for temp in [da_temp, qa_temp]:
                for band in self.boa_bands:
                    fname = self.MCD43_dir + '/'.join([datetime.strftime(self.obs_time, '%Y-%m-%d'), temp%(doy, band)])
                    fnames.append(fname)
        das, ws = self._read_MCD43(fnames)   
        mg = gdal.Warp('',fnames[0], format = 'MEM', dstNodata= None, xRes = self.aero_res*0.5, yRes = \
                       self.aero_res*0.5, cutlineDSName=self.aoi, cropToCutline=True, resampleAlg = 0)
        hg = self.example_file
        self.hx, self.hy, hmask, rmask = self._get_index(mg, hg)

        No_band = len(self.toa_bands)
        mask = ~(mg.ReadAsArray()[0]==0)
        self._annoying_angles(mg)
        sza = np.repeat(self._sza[None, ...], len(self._vza), axis = 0)[:,mask][:, hmask][:, rmask]
        saa = np.repeat(self._saa[None, ...], len(self._vza), axis = 0)[:,mask][:, hmask][:, rmask]
        angles = self._vza[:,mask][:, hmask][:, rmask], sza ,self._vaa[:,mask][:, hmask][:, rmask] - saa
        kk     = get_kk(angles)   
        k_vol  = kk.Ross          
        k_geo  = kk.Li  
        kers = np.array([np.ones(k_vol.shape), k_vol, k_geo])
        surs = []
        for i in range(No_band):
            surs.append((das[i::No_band][:,:,mask][:,:,hmask][:,:,rmask] * kers[:,i] * 0.001).sum(axis=1))
        if temporal_filling:
            
            def smooth(da_w):
                da, w = da_w
                mid = int(da.shape[0]/2)
                if (da.shape[-1]==0) | (w.shape[-1]==0):
                    return da[mid], w[mid]
                data  = np.array(smoothn(da, s=10., smoothOrder=1., axis=0, TolZ=0.001, verbose=False, isrobust=True, W = w))[[0, 3],]
                return data[0][mid], data[1][mid]
            boa = []
            w = []
            for i in range(No_band):
                das = surs[i]
                Ws  = ws[i::No_band][:,mask][:, hmask][:, rmask]
                chunks = zip(np.array_split(das, 18, axis=1), np.array_split(Ws, 18,  axis=1))
                ret = parmap(smooth, chunks)
                _b  = np.hstack([i[0] for i in ret])                   
                _w  = np.hstack([i[1] for i in ret]) 
                boa.append(_b)
                w.append(_w)    
            boa = np.array(boa)
            w = np.array(w)
            unc = 0.015 / w
        else:
            boa = np.array(surs)
            unc = 0.015 / ws
        self.boa     = boa
        self.boa_unc = np.minimum(unc, 0.5)

    def _annoying_angles(self, destination):
        mg = destination
        _sun_angles = [] 
        _view_angles = []
        for fname in self._sun_angles:     
            try:                           
                nodatas = ' '.join([i.split("=")[1] for i in gdal.Info(fname).split('\n') if' NoData' in i])
            except:                        
                nodatas = None 
            ang = reproject_data(fname, mg, srcNodata = None, resample = \
                                 0, dstNodata=np.nan, outputType= gdal.GDT_Float32).data
            _sun_angles.append(ang)        
        for fname in self._view_angles:    
            try:                           
                nodatas = ' '.join([i.split("=")[1] for i in gdal.Info(fname).split('\n') if' NoData' in i])
            except:                        
                nodatas = None 
            ang = reproject_data(fname, mg, srcNodata = None, resample = \
                                 0, dstNodata=np.nan, outputType= gdal.GDT_Float32).data
            _view_angles.append(ang)       
        _view_angles = np.squeeze(np.array(_view_angles))
        _sun_angles = np.squeeze(np.array(_sun_angles))
        if len(self.view_angles) == len(self.toa_bands):
            if self.a_z_order == 1:        
                self._vaa = _view_angles[:,0] 
                self._vza = _view_angles[:,1]
            else:                          
                self._vaa = _view_angles[:,1] 
                self._vza = _view_angles[:,0]
        elif len(self.view_angles) == 2*len(self.toa_bands):
            self._vaa = _view_angles[:len(self.toa_bands)]           
            self._vza = _view_angles[len(self.toa_bands):]           
         
        if os.path.exists(str(self.sun_angles)):
            if self.a_z_order == 1:        
                self._saa = _sun_angles[0] 
                self._sza = _sun_angles[1] 
            else:                          
                self._saa = _sun_angles[1] 
                self._sza = _sun_angles[0] 
        elif len(self.sun_angles) == 2:    
            self._saa = _sun_angles[0]     
            self._sza = _sun_angles[1]
        self._sza = self._sza * self.ang_scale
        self._saa = self._saa * self.ang_scale
        self._vza = self._vza * self.ang_scale
        self._vaa = self._vaa * self.ang_scale
        

    def _get_index(self, mg, hg):
        '''
        Pixel indexes between high resolution image and MCD43.
        '''
        temp_data = ~(mg.ReadAsArray()[0]==0)
        geotransform = mg.GetGeoTransform() 
        xgeo = geotransform[0] + np.arange(0.5, mg.RasterXSize+0.5, 1) * geotransform[1]
        ygeo = geotransform[3] + np.arange(0.5, mg.RasterYSize+0.5, 1) * geotransform[5]
        xgeo = np.repeat(xgeo[None,...], mg.RasterYSize, axis=0)
        ygeo = np.repeat(ygeo[...,None], mg.RasterXSize, axis=1)
        m_proj = modis_sinu = osr.SpatialReference()
        m_proj.ImportFromProj4("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
        h_proj = osr.SpatialReference() 
        h_proj.ImportFromWkt(hg.GetProjection())
        h_xgeo, h_ygeo, _ = np.array(osr.CoordinateTransformation(m_proj, \
                            h_proj).TransformPoints(list(zip(xgeo[temp_data], ygeo[temp_data])))).T
        geotransform = hg.GetGeoTransform()
        hy = ((h_xgeo - geotransform[0])/geotransform[1]).astype(int)
        hx = ((h_ygeo - geotransform[3])/geotransform[5]).astype(int)
        hmask = (hx>=0) & (hx<hg.RasterYSize) & (hy>=0) & (hy<hg.RasterXSize)
        hy = hy[hmask]
        hx = hx[hmask]
        rmask = ~self.bad_pix[hx, hy]
        hy = hy[rmask]
        hx = hx[rmask]
        return hx, hy, hmask, rmask

    def _load_xa_xb_xc_emus(self,):
        xap_emu = glob(self.emus_dir + '/isotropic_%s_emulators_correction_xap_%s.pkl'%(self.sensor, self.satellite))[0]
        xbp_emu = glob(self.emus_dir + '/isotropic_%s_emulators_correction_xbp_%s.pkl'%(self.sensor, self.satellite))[0]
        xcp_emu = glob(self.emus_dir + '/isotropic_%s_emulators_correction_xcp_%s.pkl'%(self.sensor, self.satellite))[0]
        if sys.version_info >= (3,0):
            f = lambda em: pkl.load(open(em, 'rb'), encoding = 'latin1')
        else:     
            f = lambda em: pkl.load(open(str(em), 'rb'))
        self.emus = parmap(f, [xap_emu, xbp_emu, xcp_emu])

    def _get_convolved_toa(self,):       
                                         
        imgs = [band_g.ReadAsArray() for band_g in self._toa_bands]                       
        self.bad_pixs = self.bad_pix[self.hx, self.hy]
        xgaus  = np.exp(-2.*(np.pi**2)*(self.psf_xstd**2)*((0.5 * np.arange(self.full_res[0]) /self.full_res[0])**2))
        ygaus  = np.exp(-2.*(np.pi**2)*(self.psf_ystd**2)*((0.5 * np.arange(self.full_res[1]) /self.full_res[1])**2))
        gaus_2d = np.outer(xgaus, ygaus) 
        def convolve(img, gaus_2d, hx, hy):
            dat = idct(idct(dct(dct(img, axis=0, norm = 'ortho'), axis=1, \
                  norm='ortho') * gaus_2d, axis=1, norm='ortho'), axis=0, norm='ortho')[hx, hy]
            return dat
        par = partial(convolve, gaus_2d = gaus_2d, hx = self.hx, hy = self.hy)
        if np.array(self.ref_scale).ndim ==2:
            self.ref_scale = self.ref_scale[self.hx, self.hy]
        if np.array(self.ref_off).ndim == 2:
            self.ref_off = self.ref_off[self.hx, self.hy]
        self.toa  = np.array(parmap(par,imgs)) * self.ref_scale+self.ref_off


    def _re_mask(self,):
        boa_mask = np.all(self.boa >= 0.001, axis = 0) &\
                   np.all(self.boa < 1,      axis = 0)
        toa_mask = ~self.bad_pix[self.hx, self.hy]
        _mask    = boa_mask & toa_mask
        self.hx      = self.hx       [_mask]
        self.hy      = self.hy       [_mask]
        self.toa     = self.toa    [:, _mask] 
        self.boa     = self.boa    [:, _mask] * self.spec_slope[...,None] + self.spec_off[...,None]
        self.boa_unc = self.boa_unc[:, _mask]
        eps=1.35
        mask = True                                 
        if self.boa.shape[1] > 3: 
            for i in range(len(self.toa)):           
                x,y = self.boa[i][...,None], self.toa[i][...,None]
                huber = HuberRegressor(fit_intercept=True, alpha=0.0, max_iter=100,epsilon=eps)
                huber.fit(x,y)                          
                mask *= ~huber.outliers_ 
            self.mask = mask #& boa_mask & toa_mask
        else:
            self.mask = False
    
    def _fill_nan(self,):
        def fill_nan(array):                        
            x_shp, y_shp = array.shape                     
            mask  = ~np.isnan(array)                       
            valid = np.array(np.where(mask)).T             
            value = array[mask]                            
            mesh  = np.repeat(range(x_shp), y_shp).reshape(x_shp, y_shp), \
                    np.tile  (range(y_shp), x_shp).reshape(x_shp, y_shp)
            array = griddata(valid, value, mesh, method='nearest')
            return array
        self._vza = np.array(parmap(fill_nan, list(self._vza)))
        self._vaa = np.array(parmap(fill_nan, list(self._vaa)))
        self._saa, self._sza, self._ele, self._aot, self._tcwv, self._tco3 = \
        parmap(fill_nan, [self._saa, self._sza, self._ele, self._aot, self._tcwv, self._tco3])
        self._aot = self._aot * 1.3 - 0.08
        self._aot = np.maximum(self._aot, 0)

    def _solving(self,):
        self.logger.propagate = False
        self.logger.info('Set AOI.')
        self._create_base_map()
        self.logger.info('Get corresponding bands.')
        self._find_boa_bands()
        self.logger.info('Slice TOA bands based on AOI.')
        self._create_band_gs()
        self._resamplers()
        self.logger.info('Parsing angles.')
        self._parse_angles()
        self.logger.info('Mask bad pixeles.')
        self._mask_bad_pix()
        if not np.all(self.bad_pix):
            self.logger.info('Get simulated BOA.')
            self._get_boa()
            self.logger.info('Get PSF.')
            self._get_psf()
            self.logger.info('Get simulated TOA reflectance.')
            self._get_convolved_toa()
            self.logger.info('Filtering data.')
            self._re_mask()
            if self.mask is not False:
                self.logger.info('Loading emulators.')
                self._load_xa_xb_xc_emus()
                self.logger.info('Reading priors and elevation.')
                self._read_aux()
                self._fill_nan()
                self.logger.info('Mean values for prior AOT: %.02f and TCWV: %.02f'%(self._aot.mean(), self._tcwv.mean()))
                if self.mask.sum() ==0:
                    self.logger.info('No valid value is found for retrieval of atmospheric parameters and priors are stored.')
                    ret = np.array([[self._aot, self._tcwv, self._tco3], [self._aot_unc, self._tcwv_unc, self._tco3_unc]])
                    self.aero_res /=2
                    #self.ySize *=2
                    #self.xSize *=2
                    self.ySize, self.xSize = self._aot.shape
                else:
                    self.aero = solving_atmo_paras(self.boa, 
                                                   self.toa,
                                                   self._sza, 
                                                   self._vza,
                                                   self._saa, 
                                                   self._vaa,
                                                   self._aot, 
                                                   self._tcwv,
                                                   self._tco3, 
                                                   self._ele,
                                                   self._aot_unc,
                                                   self._tcwv_unc,
                                                   self._tco3_unc,
                                                   self.boa_unc,
                                                   self.hx, self.hy,
                                                   self.mask,
                                                   self.full_res,
                                                   self.aero_res,
                                                   self.emus,
                                                   self.band_index,
                                                   self.boa_wv,
                                                   pix_res = self.pixel_res,
                                                   gamma = self.gamma,
                                                   log_file = self.log_file
                                                   )
                    ret = self.aero._multi_grid_solver()
            else:
                self.logger.info('No valid value is found for retrieval of atmospheric parameters and priors are stored.')
                self._read_aux()       
                self._fill_nan()
                self.logger.info('Mean values for prior AOT: %.02f and TCWV: %.02f'%(self._aot.mean(), self._tcwv.mean()))
                ret = np.array([[self._aot, self._tcwv, self._tco3], [self._aot_unc, self._tcwv_unc, self._tco3_unc]])
                self.aero_res /=2
                #self.ySize *=2
                #self.xSize *=2
                self.ySize, self.xSize = self._aot.shape
        else:
            self.logger.info('No valid value is found for retrieval of atmospheric parameters and priors are stored.')
            self._read_aux()       
            self._fill_nan()
            self.logger.info('Mean values for prior AOT: %.02f and TCWV: %.02f'%(self._aot.mean(), self._tcwv.mean()))
            ret = np.array([[self._aot, self._tcwv, self._tco3], [self._aot_unc, self._tcwv_unc, self._tco3_unc]])
            self.aero_res /=2
            self.ySize, self.xSize = self._aot.shape
            
        solved     = ret[0].reshape(3, self.ySize, self.xSize)
        unc        = ret[1].reshape(3, self.ySize, self.xSize)
        self.logger.info('Finished retrieval and saving them into local files.')
        para_names = 'aot', 'tcwv', 'tco3', 'aot_unc', 'tcwv_unc', 'tco3_unc'
        toa_dir    = self.toa_dir + '/' +'B'.join(self.toa_bands[0].split('/')[-1].split('B')[:-1])
        name_arrays     = zip(para_names, list(solved ) + list(unc))
        par = partial(save_posterior, g = self.example_file, aero_res = self.aero_res, toa_dir = toa_dir)
        parmap(par, name_arrays)
        self.post_aot,     self.post_tcwv,     self.post_tco3,    = solved
        self.post_aot_unc, self.post_tcwv_unc, self.post_tco3_unc = unc
        handlers = self.logger.handlers[:]
        for handler in handlers:
            handler.close()
            self.logger.removeHandler(handler)
    
def get_kk(angles):
    vza ,sza,raa = angles
    kk = kernels.Kernels(vza ,sza,raa,\
                         RossHS=False,MODISSPARSE=True,\
                         RecipFlag=True,normalise=1,\
                         doIntegrals=False,LiType='Sparse',RossType='Thick')
    return kk


def save_posterior(name_array, g, aero_res, toa_dir):
    name, array  = name_array
    xmin, ymax   =  g.GetGeoTransform()[0], g.GetGeoTransform()[3]
    projection   = g.GetProjection()
    xres, yres   = aero_res, aero_res
    geotransform = (xmin, xres, 0, ymax, 0, -yres)
    nx, ny       = array.shape
    outputFileName = toa_dir + '%s.tif'%name
    if os.path.exists(outputFileName):
        os.remove(outputFileName)
    dst_ds = gdal.GetDriverByName('GTiff').Create(outputFileName, ny, nx, 1, gdal.GDT_Float32, options=["TILED=YES", "COMPRESS=DEFLATE"])
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(projection)
    dst_ds.GetRasterBand(1).WriteArray(array)
    dst_ds.FlushCache()
    dst_ds = None

def test_s2():
    from datetime import datetime
    sensor_sat = 'MSI', 'S2A'
    base = '/data/nemesis/tmp/S2A_MSIL1C_20171201T094341_N0206_R036_T34TDT_20171201T114354.SAFE/'\
           +'GRANULE/L1C_T34TDT_A012760_20171201T094519/IMG_DATA/T34TDT_20171201T094341_'
    toa_bands  = [base + i + '.jp2' for i in ['B02', 'B03', 'B04', 'B08', 'B11', 'B12']]
    band_wv    = [469, 555, 645, 859, 1640, 2130]
    band_index = [1,2,3,7,11,12]
    ang_base = '/data/nemesis/tmp/S2A_MSIL1C_20171201T094341_N0206_R036_T34TDT_20171201T114354.SAFE/'+\
               'GRANULE/L1C_T34TDT_A012760_20171201T094519/ANG_DATA/VAA_VZA_'
    view_angles = [ang_base + i + '.tif' for i in ['B02', 'B03', 'B04', 'B08', 'B11', 'B12']]
    sun_angles = ang_base.replace('VAA_VZA_', 'SAA_SZA.tif')
    obs_time = datetime(2017, 12, 1, 9, 45, 19)
    cloud_mask = gdal.Open(base.replace('IMG_DATA/T34TDT_20171201T094341_', 'cloud.tif')).ReadAsArray() / 100.
    cloud_mask = cloud_mask > 0.6
    aero = solve_aerosol(sensor_sat,toa_bands,band_wv, band_index,view_angles,sun_angles,obs_time,cloud_mask, gamma=10.)
    aero._solving()
    return aero

def test_l8():
    sensor_sat = 'OLI', 'L8'
    base = '~/DATA/S2_MODIS/l_data/LC08_L1TP_014034_20170831_20170915_01_T1/LC08_L1TP_014034_20170831_20170915_01_T1_'
    toa_bands  = [base + i + '.TIF' for i in ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']]
    view_angles = [base + 'sensor_%s'%i +'.tif'  for i in ['B02', 'B03', 'B04', 'B05', 'B06', 'B07']]
    sun_angles = base + 'solar_B01.tif'
    sza = gdal.Open(sun_angles).ReadAsArray()[1] * 0.01
    scale = 2.0000E-05 / np.cos(np.deg2rad(sza))
    off   = -0.100000 / np.cos(np.deg2rad(sza))
    cloud_mask = gdal.Open(base + 'BQA.TIF').ReadAsArray()   
    cloud_mask = ~((cloud_mask  >= 2720) & ( cloud_mask <= 2732))
    band_wv    = [469, 555, 645, 859, 1640, 2130]
    obs_time = datetime(2017, 8, 31, 15, 40, 46)  
    band_index = [1,2,3,4,5,6]
    aero = solve_aerosol(sensor_sat,toa_bands,band_wv, band_index,view_angles,sun_angles,obs_time,cloud_mask, gamma=10., ref_scale = scale, ref_off = off)
    aero._solving()
    return aero

def test_modis():
    import sys
    import os
    sys.path.insert(0, '/data/store01/data_dirs/students/ucfafyi/Multiply/Atmo_cor')
    from modis_l1b_reader import MODIS_L1b_reader
    from datetime import datetime
    modis_l1b = MODIS_L1b_reader("/data/selene/ucfajlg/Ujia/MODIS_L1b/GRIDDED/", "h17v05",2017)
    obs_time = datetime(2017, 9, 4, 11, 25)
    tile = "h17v05"
    a = modis_l1b.granules[obs_time]
    scales = a.scale
    bands = [a.b1, a.b2, a.b3, a.b4, a.b5, a.b6, a.b7]
    scales = a.scale
    sun_angles = [a.saa, a.sza]
    view_angles = [a.vaa] * 6 + [a.vza] * 6
    sza = gdal.Open(a.sza).ReadAsArray()/100.
    cos_sza = np.cos(np.deg2rad(sza))
    toa_dir = './MODIS/'+ obs_time.strftime('%Y_%m_%d')
    if not os.path.exists(toa_dir):
        os.mkdir(toa_dir)
    for i in range(7):
        g = gdal.Open(bands[i])
        array = g.ReadAsArray().astype(float) * scales[i] / cos_sza
        driver = gdal.GetDriverByName('GTiff')
        band_name = '/MODIS_' + tile + obs_time.strftime('_%Y%m%d_%H%M_') + 'band%d.tif'%(i+1)
        if os.path.exists(toa_dir + band_name):                                                        
            os.remove(toa_dir + band_name) 
        ds = driver.Create(toa_dir + band_name, 2400, 2400, 1, gdal.GDT_Float32, options=["TILED=YES", "COMPRESS=DEFLATE"])
        geotransform = g.GetGeoTransform()
        projection  = g.GetProjection()
        ds.SetProjection(projection)
        ds.SetGeoTransform(geotransform)
        ds.GetRasterBand(1).WriteArray(array)
        ds.FlushCache()
        ds = None
    toa_bands  =  [toa_dir + '/MODIS_' + tile + obs_time.strftime('_%Y%m%d_%H%M_') + 'band%d.tif'%i for i in [3,4,1,2,6,7]]
    cloud_mask = np.zeros((2400, 2400)).astype(bool)
    scale = 1
    off   = 0
    band_wv    = [469, 555, 645, 859, 1640, 2130]
    gamma = 10.
    sensor_sat = 'TERRA', 'MODIS'
    band_index = [2, 3,0,1,5,6]
    emus_dir   = '/data/store01/data_dirs/students/ucfafyi/Multiply/emus/old_emus/'
    aero = solve_aerosol(sensor_sat,toa_bands,band_wv, band_index,view_angles,sun_angles,obs_time,cloud_mask,gamma=gamma, emus_dir =emus_dir, ref_scale = scale, ref_off = off, aero_res = 3000)
    aero._solving() 
    return aero
if __name__ == '__main__':
    #l8_aero = test_l8()
    s2_aero = test_s2()
    #m_aero = test_modis()
