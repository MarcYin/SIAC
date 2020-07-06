#!/usr/bin/env python
import os
import gc
import sys
import gdal
import logging
import lightgbm
import datetime
import numpy as np
from glob import glob
import multiprocessing
from SIAC.Two_NN import Two_NN
from SIAC.smoothn import smoothn
#from sklearn.externals import joblib 
from SIAC.reproject import reproject_data
from SIAC.s2_angle import resample_s2_angles
from SIAC.create_logger import create_logger
from scipy.interpolate import NearestNDInterpolator
from skimage.morphology import disk, binary_dilation, binary_erosion


import warnings
warnings.filterwarnings("ignore")

try:
    import cPickle as pkl
except:
    import pickle as pkl
file_path = os.path.dirname(os.path.realpath(__file__))
gc.disable()
cl = pkl.load(open(file_path + '/data/sen2cloud_detector.pkl', 'rb'))
gc.enable()
cl.n_jobs = multiprocessing.cpu_count()


def do_cloud(cloud_bands, cloud_name = None):
    toas = [reproject_data(str(band), cloud_bands[0], dstNodata=0, resample=5).data for band in cloud_bands]
    toas = np.array(toas)/10000.
    mask = np.all(toas >= 0.0001, axis=0)
    valid_pixel = mask.sum()
    cloud_proba = cl.predict_proba(toas[:, mask].T)
    cloud_mask = np.zeros_like(toas[0])
    cloud_mask[mask] = cloud_proba[:,1]
    cloud_mask[~mask] = -256
    if cloud_name is None:
        return cloud_mask
    else:
        g =  gdal.Open(cloud_bands[0])
        dst = gdal.GetDriverByName('GTiff').Create(cloud_name, g.RasterXSize, g.RasterYSize, 1, gdal.GDT_Byte,  options=["TILED=YES", "COMPRESS=DEFLATE"])
        dst.SetGeoTransform(g.GetGeoTransform())
        dst.SetProjection  (g.GetProjection())
        dst.GetRasterBand(1).WriteArray((cloud_mask * 100).astype(int))
        dst.GetRasterBand(1).SetNoDataValue(-256)
        dst=None; g=None
        return cloud_mask


def do_aot_tcwv(aot_bands, tcwv_bands, cams_dir, obs_time, sun_ang_name, view_ang_name, dem, tcwv_name = None, aot_name = None):

    toas = [reproject_data(str(band), tcwv_bands[0], dstNodata=0, resample=5, xRes=120, yRes=120).data for band in aot_bands]
    toas = np.array(toas)/10000.

    mask = np.all(toas >= 0.0001, axis=0)


    time_ind    = np.abs((obs_time.hour  + obs_time.minute/60. + obs_time.second/3600.) - np.arange(0,25,3)).argmin()
    prior_uncs  = 0.1
    prior_scale = 46.698
    prior_f = cams_dir + '/'.join([datetime.datetime.strftime(obs_time, '%Y_%m_%d'),
                                       datetime.datetime.strftime(obs_time, '%Y_%m_%d')+'_gtco3.tif'])
    var_g   = gdal.Open(prior_f)
    prior_g = reproject_data(prior_f , tcwv_bands[0], dstNodata=0, resample=1, xRes=120, yRes=120).g
    g       = var_g.GetRasterBand(int(time_ind+1))
    offset  = g.GetOffset()
    scale   = g.GetScale()
    tco3    = prior_g.GetRasterBand(int(time_ind+1)).ReadAsArray() * scale + offset
    tco3[:] = np.nanmean(tco3) * prior_scale

    saa, sza = reproject_data(str(sun_ang_name) , tcwv_bands[0], dstNodata=0, resample=1, xRes=120, yRes=120).data / 100.
    vaa, vza = reproject_data(str(view_ang_name), tcwv_bands[0], dstNodata=0, resample=1, xRes=120, yRes=120).data / 100.
    raa      = vaa - saa

    sza      = np.cos(np.deg2rad(sza))
    vza      = np.cos(np.deg2rad(vza))
    raa      = np.cos(np.deg2rad(raa))

    ele      = reproject_data(dem, tcwv_bands[0], dstNodata=0, resample=1, xRes=120,  yRes=120).data / 10000.

    X        = np.vstack([toas[[9,8]], np.array([sza, vza, raa, tco3, ele])]).reshape(7, -1).T
    iso_tcwv = Two_NN(np_model_file = file_path + '/emus/S2_TCWV.npz')
    tcwv     = iso_tcwv.predict(X)[0].reshape(toas[0].shape).astype(float)

    bad = (tcwv<0) |  (tcwv>8) | (~mask ) | np.isnan(tcwv)
    tcwv[bad] = np.nanmedian(tcwv[~bad])

    bands_min = [0.11572497, 0.08986528, 0.07280412, 0.05007033, 0.06228712, 0.06849915, 0.07015634, 0.06554198, 0.06809415, 0.02089167, 0.00080242, 0.04392302, 0.03012728]
    bands_max = [0.32072931, 0.30801709, 0.32420026, 0.38214899, 0.3951169 , 0.41631785, 0.44391895, 0.42444658, 0.46161492, 0.19927658, 0.00870671, 0.39533241, 0.32505772]
    
    for _, toa in enumerate(toas):
            mask = mask & (toa >= bands_min[_]) & (toa <= bands_max[_])
            
    X        = np.vstack([toas, np.array([sza, vza, raa, tco3, ele])]).reshape(18, -1).T
    gbm      = joblib.load(file_path + '/emus/lgb.pkl')
    aot      = gbm.predict(X).reshape(toas[0].shape).astype(float)
    aot      = np.exp(-1*aot)
    shape    = toas[0].shape
    if mask.sum()>3:
        aot_min, aot_median, aot_max = np.nanpercentile(aot[mask], [5, 50, 95])
        good_aot =  (aot >= aot_min) & (aot <= aot_max) 
        indx, indy = np.where(good_aot.reshape(shape))
        myInterpolator = NearestNDInterpolator((indx, indy), aot[good_aot]) 
        grid_x, grid_y = np.mgrid[0:shape[0]:1, 0:shape[1]:1,]
        aot = myInterpolator(grid_x, grid_y)
        aot = smoothn(aot, isrobust=True, TolZ=1e-6, W = 100*((aot >= aot_min) & (aot <=aot_max)), s=10, MaxIter=1000)[0]
    else:
        aot = np.nanmedian(aot)*np.ones(shape)    

    if tcwv_name is not None:
        g =  gdal.Open(tcwv_bands[0])
        ySize, xSize = tcwv.shape
        dst = gdal.GetDriverByName('GTiff').Create(tcwv_name, xSize, ySize, 1, gdal.GDT_UInt16,  options=["TILED=YES", "COMPRESS=DEFLATE"])
        dst.SetGeoTransform(g.GetGeoTransform())
        dst.SetProjection  (g.GetProjection())
        dst.GetRasterBand(1).WriteArray((tcwv*1000).astype(int))
        dst.GetRasterBand(1).SetNoDataValue(65535)
        dst=None; g=None

    if aot_name is not None:
        g =  gdal.Open(aot_bands[0])
        ySize, xSize = aot.shape
        dst = gdal.GetDriverByName('GTiff').Create(aot_name, xSize, ySize, 1, gdal.GDT_UInt16,  options=["TILED=YES", "COMPRESS=DEFLATE"])
        dst.SetGeoTransform(g.GetGeoTransform())
        dst.SetProjection  (g.GetProjection())
        dst.GetRasterBand(1).WriteArray((aot * 1000).astype(int))
        dst.GetRasterBand(1).SetNoDataValue(65535)
        dst=None; g=None

    return aot, tcwv


def s2_pre_processing(s2_dir, cams_dir, dem):
    s2_dir = os.path.abspath(s2_dir)
    scihub = []
    aws    = []
    for (dirpath, dirnames, filenames)  in os.walk(s2_dir):
        if len(filenames)>0:
            temp = [dirpath + '/' + i for i in filenames]
            for j in temp:
                if ('MTD' in j) & ('TL' in j) & ('xml' in j):
                    scihub.append(j)
                if 'metadata.xml' in j:
                    aws.append(j)
    s2_tiles = []
    for metafile in scihub + aws:
        with open(metafile) as f:
            for i in f.readlines():
                if 'TILE_ID' in i:
                    tile  = i.split('</')[0].split('>')[-1]
                if 'SENSING_TIME' in i:
                    sensing_time = i.split('</')[0].split('>')[-1]
                    obs_time = datetime.datetime.strptime(sensing_time, u'%Y-%m-%dT%H:%M:%S.%fZ')

        log_file = os.path.dirname(metafile) + '/SIAC_S2.log'
        if os.path.exists(log_file):
            os.remove(log_file)
        logger = create_logger(log_file)
        logger.propagate = False
        logger.info('Preprocessing for %s'%tile)
        logger.info('Doing per pixel angle resampling.')
        sun_ang_name, view_ang_names, toa_refs, cloud_name = resample_s2_angles(metafile)
        cloud_bands = np.array(toa_refs)[[0,1,3,4,7,8,9,10,11,12]]
        logger.info('Getting cloud mask.')
        cloud = do_cloud(cloud_bands, cloud_name)
        cloud_mask = cloud > 0.6# 0.4 is the sentinel hub default
        clean_pixel = ((cloud>=0).sum() - cloud_mask.sum()) / 3348900. * 100.
        valid_pixel = (cloud>=0).sum() / 3348900. * 100.
        logger.info('Valid pixel percentage: %.02f'%(valid_pixel))
        logger.info('Clean pixel percentage: %.02f'%(clean_pixel))
        cloud_mask = binary_dilation(binary_erosion (cloud_mask, disk(2)), disk(3)) | (cloud < 0)
        
        #logger.info('Getting first guess AOT and TCWV.')
        #tcwv_bands = np.array(toa_refs)[[9,8]]
        #aot_bands = np.array(toa_refs)
        #view_ang_name = view_ang_names[8]
        #aot, tcwv     = do_aot_tcwv(aot_bands, tcwv_bands, cams_dir, obs_time, sun_ang_name, view_ang_name, dem, tcwv_name = None, aot_name = None)
        aot = None
        tcwv = None
        s2_tiles.append([sun_ang_name, view_ang_names, toa_refs, cloud_name, cloud_mask, aot, tcwv, metafile])
        handlers = logger.handlers[:]
        for handler in handlers:
            handler.close()
            logger.removeHandler(handler)
    return s2_tiles
