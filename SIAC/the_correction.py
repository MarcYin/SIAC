#!/usr/bin/env python 
import os
import sys
import psutil
import logging
import warnings
import subprocess
import numpy as np
from copy import copy
from glob import glob
try:
    import cPickle as pkl
except:
    import pickle as pkl
from osgeo import gdal, ogr
from numpy import clip, uint8
from SIAC.multi_process import parmap
from scipy.interpolate import griddata
from scipy.ndimage import binary_dilation
from SIAC.create_logger import create_logger
from SIAC.raster_boundary import get_boundary
from SIAC.reproject import reproject_data, array_to_raster

warnings.filterwarnings("ignore")

class atmospheric_correction(object):
    ''' 
    A class doing the atmospheric coprrection with the input of TOA reflectance
    angles, elevation and emulators of 6S from TOA to surface reflectance.
    ''' 
    def __init__(self,
                 sensor_sat,
                 toa_bands, 
                 band_index,
                 view_angles,
                 sun_angles,
                 aoi         = None,
                 aot         = None, 
                 tcwv        = None,
                 tco3        = None,
                 aot_unc     = None,
                 tcwv_unc    = None,
                 tco3_unc    = None,
                 cloud_mask  = None,
                 obs_time    = None,
                 log_file    = None,
                 a_z_order   = 1,
                 ref_scale   = 0.0001,
                 ref_off     = 0,
                 ang_scale   = 0.01,
                 ele_scale   = 0.001,
                 atmo_scale  = [1., 1., 1., 1., 1., 1.],
                 global_dem  = '/vsicurl/http://www2.geog.ucl.ac.uk/~ucfafyi/eles/global_dem.vrt',
                 cams_dir    = '/vsicurl/http://www2.geog.ucl.ac.uk/~ucfafyi/cams/',
                 emus_dir    = 'SIAC/emus/',
                 cams_scale  = [1., 0.1, 46.698, 1., 1., 1.],
                 block_size  = 600,
                 rgb         = [None, None, None]
                
                 ):      
    

        self.sensor      = sensor_sat[0]
        self.satellite   = sensor_sat[1]
        self.toa_bands   = toa_bands
        self.band_index  = band_index
        self.view_angles = view_angles
        self.sun_angles  = sun_angles 
        self.aoi         = aoi
        self.aot         = aot
        self.tcwv        = tcwv
        self.tco3        = tco3
        self.aot_unc     = aot_unc
        self.tcwv_unc    = tcwv_unc
        self.tco3_unc    = tco3_unc
        self.cloud_mask  = cloud_mask
        self.obs_time    =  obs_time
        self.a_z_order   = a_z_order
        self.ref_scale   = ref_scale
        self.ref_off     = ref_off
        self.ang_scale   = ang_scale
        self.ele_scale   = ele_scale
        self.atmo_scale  = atmo_scale
        self.global_dem  = global_dem
        self.emus_dir    = emus_dir
        self.cams_dir    = cams_dir
        self.cams_scale  = cams_scale
        self.block_size  = block_size
        self.toa_dir     =  os.path.abspath('/'.join(toa_bands[0].split('/')[:-1]))
        self.rgb         = rgb
        r, g, b          = self.rgb
        self._do_rgb = False
        self.ri, self.gi, self.bi = None, None, None
        if (r is not None) & (g is not None) & (b is not None):
            self.ri, self.gi, self.bi = self.toa_bands.index(r), self.toa_bands.index(g), self.toa_bands.index(b)
            self._do_rgb = True

        self.logger = create_logger(log_file)
         
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
                    #subprocess.call(['gdaltindex', '-f', 'GeoJSON',  '-t_srs', 'EPSG:4326', self.toa_dir + '/AOI.json', self.aoi])
                    geojson = get_boundary(self.aoi)[0]
                    with open(self.toa_dir + '/AOI.json', 'wb') as f:
                        f.write(geojson.encode())
                except:  
                    try: 
                        gr = ogr.Open(self.aoi)
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

    def _create_band_gs(self,):
        '''
        Create a lost of boa gs and cut them with AOI.
        '''
        self._toa_bands = []
        for band in self.toa_bands:
            g = gdal.Warp('', band, format = 'MEM', srcNodata = 0, dstNodata=0, warpOptions = \
                          ['NUM_THREADS=ALL_CPUS'], cutlineDSName= self.aoi, cropToCutline=True, resampleAlg = 0)
            self._toa_bands.append(g)
        self.example_file = self._toa_bands[0]
        if self.cloud_mask is None:
            self.cloud_mask = False
        self.mask = self.example_file.ReadAsArray() >0# & (~self.cloud_mask)
        self.mask_g = array_to_raster(self.mask.astype(float), self.example_file)
        self.mask = self.mask.astype(bool)
        self.mg    = gdal.Warp('', band, format = 'MEM', srcNodata = None, xRes = self.block_size, yRes = \
                               self.block_size, cutlineDSName= self.aoi, cropToCutline=True, resampleAlg = 0)        
  
    def _load_xa_xb_xc_emus(self,):
        xap_emu = glob(self.emus_dir + '/isotropic_%s_emulators_correction_xap_%s.pkl'%(self.sensor, self.satellite))[0]
        xbp_emu = glob(self.emus_dir + '/isotropic_%s_emulators_correction_xbp_%s.pkl'%(self.sensor, self.satellite))[0]
        xcp_emu = glob(self.emus_dir + '/isotropic_%s_emulators_correction_xcp_%s.pkl'%(self.sensor, self.satellite))[0]
        if sys.version_info >= (3,0):
            f = lambda em: pkl.load(open(em, 'rb'), encoding = 'latin1')
        else:     
            f = lambda em: pkl.load(open(str(em), 'rb'))
        #print([xap_emu, xbp_emu, xcp_emu])
        self.xap_emus, self.xbp_emus, self.xcp_emus = parmap(f, [xap_emu, xbp_emu, xcp_emu])
    
    def _var_parser(self, var):
        if os.path.exists(str(var)):    
            var_g = gdal.Open(str(var)) 
        elif type(var).__module__== 'numpy':
            var_g = array_to_raster(var, self.example_file) 
        elif ('/vsicurl/' in str(var)) or ('/vsizip/') in str(var):
            var_g = gdal.Open(str(var))
        else:              
            var          = float(var) 
            var_array    = np.zeros((10,10))      
            var_array[:] = var
            var_g        = array_to_raster(var_array, self.example_file)
        return var_g
         
    def _parse_angles(self,):
        '''
        Parsing angles
        '''
        self._view_angles = []
        if os.path.exists(str(self.view_angles)):
            ang_g = self._var_parser(str(self.view_angles))
            self._view_angles = [ang_g,]
        else:
            for i in self.view_angles:
                ang_g = self._var_parser(i)
                self._view_angles.append(ang_g)
         
        self._sun_angles = []            
        if os.path.exists(str(self.sun_angles)):
            ang_g = self._var_parser(str(self.sun_angles))
            self._sun_angles = [ang_g,]
        else:
            for i in self.sun_angles:        
                ang_g = self._var_parser(i)
                self._sun_angles.append(ang_g)

    def _annoying_angles(self,):
        _sun_angles = [] 
        _view_angles = []
        for fname in self._sun_angles:     
            #nodatas = [float(i.split("=")[1]) for i in gdal.Info(fname).split('\n') if' NoData' in i] 
            try:
                nodatas = ' '.join([i.split("=")[1] for i in gdal.Info(fname).split('\n') if' NoData' in i])
            except:
                nodatas = None
            ang = reproject_data(fname, self.mg, srcNodata = nodatas, resample = \
                                 0, dstNodata=np.nan, outputType= gdal.GDT_Float32).data
            _sun_angles.append(ang)        
        for fname in self._view_angles:    
            try:                           
                nodatas = ' '.join([i.split("=")[1] for i in gdal.Info(fname).split('\n') if' NoData' in i])
            except:                        
                nodatas = None 
            ang = reproject_data(fname, self.mg, srcNodata = nodatas, resample = \
                                 0, dstNodata=np.nan, outputType= gdal.GDT_Float32).data
            _view_angles.append(ang)       
        _view_angles = np.array(_view_angles)
        _sun_angles = np.squeeze(np.array(_sun_angles))
        if len(self._view_angles) == len(self.toa_bands):
            if self.a_z_order == 1:        
                self._vaa = _view_angles[:,0] 
                self._vza = _view_angles[:,1]
            else:                          
                self._vaa = _view_angles[:,1] 
                self._vza = _view_angles[:,0]
        elif len(self._view_angles) == 2*len(self.toa_bands):
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
        self._vaa, self._vza, self._saa, self._sza = map(np.deg2rad, [self._vaa, self._vza, self._saa, self._sza])
    
    def _parse_atmo(self,):
        atmos       = [self.aot, self.tcwv, self.tco3]
        atmo_uncs   = [self.aot_unc, self.tcwv_unc, self.tco3_unc]
        atmo_scale  = self.atmo_scale
        cams_names  = ['aod550', 'tcwv', 'gtco3']
        defalt_uncs = [0.4, 0.1, 0.05]
        if self.obs_time is not None:
            time_ind    = np.abs((self.obs_time.hour  + self.obs_time.minute/60. + \
                                  self.obs_time.second/3600.) - np.arange(0,25,3)).argmin()
            for i in range(3):
                if atmos[i] is None:
                    atmos[i] = self.cams_dir + '/'.join([datetime.strftime(self.obs_time, '%Y_%m_%d'),\
                                                         datetime.strftime(self.obs_time, '%Y_%m_%d')+'_%s.tif'%cams_names[i]])
                    var_g  = self._var_parser(atmos[i])
                    _g     = reproject_data(var_g, self.mg, srcNodata = None, resample = \
                                            0, dstNodata=np.nan, outputType= gdal.GDT_Float32).g
                    offset   = var_g.GetOffset()            
                    scale    = var_g.GetScale()             
                    data     = _g.GetRasterBand(int(time_ind+1)).ReadAsArray() * scale + offset
                    atmos[i] = data * self.cams_scale[i]
                    if atmo_uncs[i] is None:
                        atmo_uncs[i] = defalt_uncs[i]
        else:
            for i in range(3):
                var_g    = self._var_parser(atmos[i])
                data     = reproject_data(var_g, self.mg, srcNodata = np.nan, resample = \
                                     0, dstNodata=np.nan, outputType= gdal.GDT_Float32).data
                atmos[i] = data * self.atmo_scale[i]
                unc_g    = self._var_parser(atmo_uncs[i])
                ug     = reproject_data(unc_g, self.mg, srcNodata = np.nan, resample = \
                                          0, dstNodata=np.nan, outputType= gdal.GDT_Float32).data
                atmo_uncs[i] = ug

        self._aot, self._tcwv, self._tco3             = atmos
        self._aot_unc, self._tcwv_unc, self._tco3_unc = atmo_uncs

    def _parse_aux(self,):
        auxs =  [self.global_dem, self.mask.astype(float)]
        scales = [self.ele_scale, 1]
        for i  in range(len(auxs)):
            var_g = self._var_parser(auxs[i])
            dat = reproject_data(var_g, self.mg, srcNodata = 0, resample = \
                                 0, dstNodata=np.nan, outputType= gdal.GDT_Float32).data
            auxs[i] = dat * scales[i]
        self._ele, self._cmask = auxs
        self._ele[np.isnan(self._ele)] = 0
        self._cmask[np.isnan(self._cmask)] = 0
        self._cmask = binary_dilation(self._cmask, structure = np.ones((3,3)).astype(bool)).astype(bool)
        #self._cmask  = self._cmask.astype(bool)

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
        self._saa, self._sza, self._ele, self._aot, self._tcwv, self._tco3, self._aot_unc, self._tcwv_unc, self._tco3_unc = \
        parmap(fill_nan, [self._saa, self._sza, self._ele, self._aot, self._tcwv, self._tco3, self._aot_unc, self._tcwv_unc, self._tco3_unc])
        self._aot_unc = array_to_raster(self._aot_unc, self.example_file)
        self._tcwv_unc = array_to_raster(self._tcwv_unc, self.example_file)
        self._tco3_unc = array_to_raster(self._tco3_unc, self.example_file)

    def _get_xps(self,):
        p   = [np.cos(self._sza), None, None, self._aot, self._tcwv, self._tco3, self._ele]   
        raa = np.cos(self._saa - self._vaa)
        vza = np.cos(self._vza)
        xap, xap_dH = [], []               
        xbp, xbp_dH = [], []               
        xcp, xcp_dH = [], []               
        xps = [xap, xbp, xcp]              
        xhs = [xap_dH, xbp_dH, xcp_dH]
        for i in range(len(self.toa_bands)):
            emus_index = self.band_index[i]
            X = copy(p)
            X[1:3] = vza[i], raa[i]
            X = np.array(X).reshape(7, -1)[:, self._cmask.ravel()]
            for ei, emu in enumerate([self.xap_emus, self.xbp_emus, self.xcp_emus]):
                H, _, dH = emu[emus_index].predict(X.T, do_unc=True)

                temp = np.zeros_like(self._sza)
                temp[self._cmask] = H
                xps[ei].append(array_to_raster(temp, self.example_file))

                temp = np.zeros(self._sza.shape + (3,))
                temp[self._cmask] = dH[:, 3:6]
                xhs[ei].append(array_to_raster(temp, self.example_file))
            
        self._saa; del self._sza; del self._ele; del self._aot; del self._tcwv; del self._tco3#; del self._aot_unc; del self._tcwv_unc; del self._tco3_unc
        self.xap_H, self.xbp_H, self.xcp_H    = np.array(xps)
        self.xap_dH, self.xbp_dH, self.xcp_dH = np.array(xhs)


    def _get_bounds(self, g):                                                                                                                                 
        geo_t = g.GetGeoTransform()
        #raster_wkt =  self.example_file.GetProjection()
        x_size, y_size = g.RasterXSize, g.RasterYSize
        xmin, xmax = min(geo_t[0], geo_t[0] + x_size * geo_t[1]), \
                     max(geo_t[0], geo_t[0] + x_size * geo_t[1])  
        ymin, ymax = min(geo_t[3], geo_t[3] + y_size * geo_t[5]), \
                     max(geo_t[3], geo_t[3] + y_size * geo_t[5])  
        bounds = [xmin, ymin, xmax, ymax]
        return bounds, x_size, y_size #, raster_wkt
    
    def _get_boa(self,):

        pix_mem = 180.
        av_ram = psutil.virtual_memory().available
        needed = np.array([i.RasterXSize * i.RasterYSize * pix_mem for i in self._toa_bands])
        u_need = np.unique(needed)
        procs  = av_ram / u_need
        if av_ram > sum(needed):
            #ret = parmap(self._do_band, range(len(self.toa_bands)))
            self._chunks = 1
            ret = parmap(self._do_chunk, range(len(self.toa_bands)))
        else:
            ret = []
            index = []
            for i, proc in enumerate(procs):
                bands_to_do = np.arange(len(self.toa_bands))[needed==u_need[i]]
                if int(proc) >= 1:
                    self._chunks = 1
                    re = parmap(self._do_chunk, bands_to_do, min(int(proc), len(bands_to_do)))
                else:
                    self._chunks = int(np.ceil(1. / proc))
                    re = list(map(self._do_chunk, bands_to_do))
                ret   += re
                index += (np.where(needed==u_need[i])[0]).tolist()
            ret = list(zip(*sorted(zip(index, ret))))[1]
        self.boa_rgb = ret[self.ri], ret[self.gi],ret[self.bi]
        #return ret


    def _do_chunk(self, ind): 
        toa_g = self._toa_bands[ind]       
        geo = toa_g.GetGeoTransform()
        x_size, y_size = toa_g.RasterXSize, toa_g.RasterYSize
        boas = []
        uncs = []
        for i in range(self._chunks):
            r = i * 1. /  self._chunks
            x_start = int(r * x_size)
            xmin = geo[0] + r * x_size * geo[1]

            r = (i+1.)/ self._chunks
            x_end = int(r * x_size)
            xmax = geo[0] + r * x_size * geo[1]

            x_off = x_end - x_start
            toa = toa_g.ReadAsArray(x_start, 0, x_off, int(y_size))
            _g = gdal.Warp('', toa_g, format = 'MEM',  outputBounds = [xmin, geo[3] + y_size * geo[5], xmax, geo[3]], resampleAlg = 0)
            xRes = yRes = int(geo[1])
            xps  = [self.xap_H[ind], self.xbp_H[ind], self.xcp_H[ind]] 
            for j in range(3): 
                xps[j]  = reproject_data(xps[j], _g,xRes=xRes, yRes = yRes,xSize=x_off, \
                                         ySize=y_size, resample = 0, outputType= gdal.GDT_Float32).data
             
            xap_H, xbp_H, xcp_H = xps
            toa = toa * self.ref_scale + self.ref_off                             
            y   = xap_H * toa - xbp_H
            boa = y / (1 + xcp_H * y)
             
            xhs = [self.xap_dH[ind], self.xbp_dH[ind], self.xcp_dH[ind]]
            for j in range(3):         
                #var_g = xhs[j]
                xhs[j]  = reproject_data(xhs[j], _g, xRes=xRes, yRes = yRes, xSize=x_off, \
                                         ySize=y_size, resample = 0, \
                                         outputType= gdal.GDT_Float32).data.transpose(1,2,0)
             
            xap_dH, xbp_dH, xcp_dH = xhs
            del xps; del xhs   
             
            dH           = -1 * (-toa[...,None] * xap_dH - \
                            2 * toa[...,None] * xap_H[...,None] * xbp_H[...,None] * xcp_dH + \
                            toa[...,None]**2 * xap_H[...,None]**2 * xcp_dH + \
                            xbp_dH + \
                            xbp_H[...,None]**2 * xcp_dH) / \
                            (toa[...,None] * xap_H[...,None] * xcp_H[...,None] - \
                                xbp_H[...,None] * xcp_H[...,None] + 1)**2
             
            toa_dH = xap_H / (xcp_H*(toa * xap_H - xbp_H) + 1)**2
            del xap_H; del xbp_H; del xcp_H; del xap_dH; del xbp_dH; del xcp_dH
            aot_dH, tcwv_dH, tco3_dH = [dH[:, :,i] for i in range(3)]
            del dH
             
            tmp = []
            for unc_g in [self._aot_unc, self._tcwv_unc, self._tco3_unc]:
                unc = reproject_data(unc_g, _g, xSize=x_off, ySize=y_size, resample = 0, outputType= gdal.GDT_Float32).data
                tmp.append(unc)
            aot_unc, tcwv_unc, tco3_unc = tmp; del tmp
            unc  = np.sqrt(aot_dH ** 2 * aot_unc**2 + tcwv_dH ** 2 * tcwv_unc**2 + tco3_dH ** 2 * tco3_unc**2 + toa_dH**2 * 0.015**2)
            del aot_unc; del tcwv_unc; del tco3_unc
            boas.append(boa)
            uncs.append(unc)
        boas, uncs = np.hstack(boas), np.hstack(uncs)
        boa_name = self.toa_dir + '/' + '.'.join(self.toa_bands[ind].split('/')[-1].split('.')[0:-1]) + '_sur.tif'
        unc_name  = boa_name.replace('_sur.tif', '_sur_unc.tif')
        projectionRef = toa_g.GetProjectionRef()
        geotransform  = toa_g.GetGeoTransform()
        self._save_band(boas, boa_name, projectionRef, geotransform)
        self._save_band(uncs, unc_name, projectionRef, geotransform)
        if ind in [self.ri, self.gi, self.bi]:
            return boas
        else:
            return None
            
    def _save_band(self, array, outputFileName, projectionRef, geotransform):
        nx, ny = array.shape
        if os.path.exists(outputFileName):
            os.remove(outputFileName)
        dst_ds = gdal.GetDriverByName('GTiff').Create(outputFileName, ny, nx, 1, gdal.GDT_Int16, options=["TILED=YES", "COMPRESS=DEFLATE"])
        dst_ds.SetGeoTransform(geotransform)
        dst_ds.SetProjection(projectionRef)
        array = array * 10000
        array[~(array>0)] = -9999
        sur = array.astype(np.int16)
        dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
        dst_ds.GetRasterBand(1).WriteArray(array)
        dst_ds.FlushCache()                  
        dst_ds = None

    def _save_rgb(self, rgba_array, outputFileName, projection, geotransform):
        nx, ny       = rgba_array.shape[1:]
        #outputFileName = self.s2_file_dir+'/%s'%name
        if os.path.exists(outputFileName):
            os.remove(outputFileName)
        dst_ds = gdal.GetDriverByName('GTiff').Create(outputFileName, ny, nx, 4, gdal.GDT_Byte, options=["TILED=YES", "COMPRESS=JPEG"])
        dst_ds.SetGeoTransform(geotransform)
        dst_ds.SetProjection(projection)
        dst_ds.GetRasterBand(1).WriteArray(rgba_array[0])
        dst_ds.GetRasterBand(1).SetScale(self.rgb_scale)
        dst_ds.GetRasterBand(2).WriteArray(rgba_array[1])
        dst_ds.GetRasterBand(2).SetScale(self.rgb_scale)
        dst_ds.GetRasterBand(3).WriteArray(rgba_array[2])
        dst_ds.GetRasterBand(4).SetScale(self.rgb_scale)
        dst_ds.GetRasterBand(4).WriteArray(rgba_array[3])
        dst_ds.FlushCache()
        dst_ds = None 
     
    def _compose_rgb(self,):
        self.rgb_scale = 4
        r, g, b = self._toa_bands[self.ri].ReadAsArray() * self.ref_scale + self.ref_off, \
                  self._toa_bands[self.gi].ReadAsArray() * self.ref_scale + self.ref_off, \
                  self._toa_bands[self.bi].ReadAsArray() * self.ref_scale + self.ref_off
        alpha   = (r>0) & (g>0) & (b>0)
        rgba_array = np.clip([r * self.rgb_scale * 255, g * self.rgb_scale * 255, \
                              b * self.rgb_scale * 255, alpha * self.rgb_scale * 255], 0, 255).astype(np.uint8)
        name = self.toa_dir + '/TOA_RGB.tif'
        projection   = self._toa_bands[self.ri].GetProjectionRef()
        geotransform = self._toa_bands[self.ri].GetGeoTransform() 
        self._save_rgb(rgba_array, name, projection, geotransform)
        gdal.Translate(self.toa_dir +'/TOA_overview.png', self.toa_dir+'/TOA_RGB.tif', \
                       format = 'PNG', widthPct=25, heightPct=25, resampleAlg=gdal.GRA_Bilinear ).FlushCache()
        gdal.Translate(self.toa_dir +'/TOA_ovr.png', self.toa_dir+'/TOA_RGB.tif', \
                       format = 'PNG', widthPct=10, heightPct=10, resampleAlg=gdal.GRA_Bilinear ).FlushCache()

        rgba_array = np.clip([self.boa_rgb[0] * self.rgb_scale * 255, self.boa_rgb[1] * self.rgb_scale * 255, \
                              self.boa_rgb[2] * self.rgb_scale * 255, alpha * self.rgb_scale * 255], 0, 255).astype(np.uint8)
        name = self.toa_dir + '/BOA_RGB.tif'
        self._save_rgb(rgba_array, name, projection, geotransform)
        gdal.Translate(self.toa_dir+'/BOA_overview.png', self.toa_dir+'/BOA_RGB.tif', \
                       format = 'PNG', widthPct=25, heightPct=25, resampleAlg=gdal.GRA_Bilinear ).FlushCache()
        gdal.Translate(self.toa_dir+'/BOA_ovr.png', self.toa_dir+'/BOA_RGB.tif', \
                       format = 'PNG', widthPct=10, heightPct=10, resampleAlg=gdal.GRA_Bilinear ).FlushCache()

    def _doing_correction(self,):
        self.logger.propagate = False
        self.logger.info('Set AOI.')
        self._create_base_map()
        self.logger.info('Slice TOA bands based on AOI.')
        self._create_band_gs()
        self.logger.info('Parsing angles.')
        self._parse_angles()
        self._annoying_angles()
        self.logger.info('Parsing auxs.')
        self._parse_aux()
        self.logger.info('Parsing atmo parameters.')
        self._parse_atmo()
        self._fill_nan()
        self.logger.info('Loading emus.')
        self._load_xa_xb_xc_emus()
        self.logger.info('Get correction coefficients.')
        self._get_xps()
        self.logger.info('Doing corrections.')
        ret = self._get_boa()
        if self._do_rgb:
            self.logger.info('Composing RGB.')
            self._compose_rgb()
        self.logger.info('Done.')
        handlers = self.logger.handlers[:]
        for handler in handlers:
            handler.close()
            self.logger.removeHandler(handler)
        return ret
         
def test_s2():
         
    from datetime import datetime
    sensor_sat = 'MSI', 'S2A'
    base = '/data/nemesis/tmp/S2A_MSIL1C_20171201T094341_N0206_R036_T34TDT_20171201T114354.SAFE/'\
           +'GRANULE/L1C_T34TDT_A012760_20171201T094519/IMG_DATA/T34TDT_20171201T094341_'
    bs = ['B01','B02', 'B03', 'B04','B05', 'B06', 'B07', 'B08','B8A','B09', 'B10','B11', 'B12']#[:1]
    toa_bands  = [base + i + '.jp2' for i in bs]
    #band_wv    = [469, 555, 645, 859, 1640, 2130]
    band_index = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    ang_base = '/data/nemesis/tmp/S2A_MSIL1C_20171201T094341_N0206_R036_T34TDT_20171201T114354.SAFE/'+\
               'GRANULE/L1C_T34TDT_A012760_20171201T094519/ANG_DATA/VAA_VZA_'
    view_angles = [ang_base + i + '.tif' for i in bs]
    sun_angles = ang_base.replace('VAA_VZA_', 'SAA_SZA.tif')
    obs_time = datetime(2017, 12, 1, 9, 45, 19)
    cloud_mask = gdal.Open(base.replace('IMG_DATA/T34TDT_20171201T094341_', 'cloud.tif')).ReadAsArray() / 100.
    aot = base.replace('T34TDT_20171201T094341_', 'aot.tif')
    tcwv = base.replace('T34TDT_20171201T094341_', 'tcwv.tif')
    tco3 = base.replace('T34TDT_20171201T094341_', 'tco3.tif')
    aot_unc = base.replace('T34TDT_20171201T094341_', 'aot_unc.tif')  
    tcwv_unc = base.replace('T34TDT_20171201T094341_', 'tcwv_unc.tif')
    tco3_unc = base.replace('T34TDT_20171201T094341_', 'tco3_unc.tif')
    cloud_mask = cloud_mask > 0.6
    rgb = [toa_bands[3], toa_bands[2], toa_bands[1]]
    atmo = atmospheric_correction(sensor_sat,toa_bands, band_index,view_angles,sun_angles, aot = aot, cloud_mask = cloud_mask,\
                                  tcwv = tcwv, tco3 = tco3, aot_unc = aot_unc, tcwv_unc = tcwv_unc, tco3_unc = tco3_unc, rgb = rgb)
    ret = atmo._doing_correction()
    return ret, atmo
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
    band_index = [1,2,3,4,5,6]
    aot = base + 'aot.tif'
    tcwv = base + 'tcwv.tif'
    tco3 = base + 'tco3.tif'
    aot_unc = base + 'aot_unc.tif'
    tcwv_unc = base + 'tcwv_unc.tif'
    tco3_unc = base + 'tco3_unc.tif'
    #cloud_mask = cloud_mask > 0.6
    cloud_mask = gdal.Open(base + 'BQA.TIF').ReadAsArray()
    cloud_mask = ~((cloud_mask  >= 2720) & ( cloud_mask <= 2732))

    rgb = [toa_bands[2], toa_bands[1], toa_bands[0]]
    atmo = atmospheric_correction(sensor_sat, toa_bands, band_index,view_angles,sun_angles, aot = aot, cloud_mask = cloud_mask,\
                                  tcwv = tcwv, tco3 = tco3, aot_unc = aot_unc, tcwv_unc = tcwv_unc, tco3_unc = tco3_unc, rgb = rgb, ref_scale = scale, ref_off = off)
    ret = atmo._doing_correction()
    return ret, atmo    

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
    sensor_sat = 'TERRA', 'MODIS'   
    band_index = [2, 3,0,1,5,6]
    aot = toa_dir + '/aot.tif'
    tcwv = toa_dir + '/tcwv.tif'
    tco3 = toa_dir + '/tco3.tif'
    aot_unc = toa_dir + '/aot_unc.tif'
    tcwv_unc = toa_dir + '/tcwv_unc.tif'
    tco3_unc = toa_dir + '/tco3_unc.tif'
    emus_dir   = '~/DATA/Multiply/emus/old_emus/'
    rgb = [toa_bands[2], toa_bands[1], toa_bands[0]]
    atmo = atmospheric_correction(sensor_sat, toa_bands, band_index,view_angles,sun_angles, aot = aot, cloud_mask = cloud_mask, ref_scale = scale, ref_off = off,\
                                  tcwv = tcwv, tco3 = tco3, aot_unc = aot_unc, tcwv_unc = tcwv_unc, tco3_unc = tco3_unc, rgb = rgb, emus_dir =emus_dir,block_size=3000)
    ret = atmo._doing_correction()
    return ret, atmo

if __name__ == '__main__':
    s_ret, s_atmo = test_s2()
    #l_ret, l_atmo = test_l8()
    #m_ret, m_atmo = test_modis()

