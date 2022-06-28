import os
import sys
import psutil
import numpy as np
from numba import jit
import multiprocessing
from functools import partial
import SIAC.kernels as kernels
from osgeo import ogr, osr, gdal
from multiprocessing import Pool
from SIAC.smoothn import smoothn
from scipy import ndimage, signal
from datetime import datetime, timedelta
from scipy.linalg import get_lapack_funcs
if (sys.version_info[0] == 3) & (sys.version_info[1] >= 4):
    multiprocessing.set_start_method('spawn', force=True)

procs =  psutil.cpu_count()
# def warp_data(fname, aoi,  xRes, yRes):
#     g = gdal.Warp('',fname, format = 'MEM', srcNodata = 32767, dstNodata=0, outputType = gdal.GDT_Float32,\
# 		  cutlineDSName=aoi, xRes = xRes, yRes = yRes, cropToCutline=True, resampleAlg = 0, warpOptions = ['NUM_THREADS=ALL_CPUS']) # weird adaptation for gdal 2.3, this should be a bug in gdal 2.3
#     return g.ReadAsArray()                                                                          # no reason you have to specify the srcNodata to use dstNodata

def get_kk(angles):
    vza ,sza,raa = angles
    kk = kernels.Kernels(vza ,sza,raa,\
                         RossHS=False,MODISSPARSE=True,\
                         RecipFlag=True,normalise=1,\
                         doIntegrals=False,LiType='Sparse',RossType='Thick')
    return kk

def warp_data(fname, outputBounds,  xRes, yRes,  dstSRS):
    modis_sinu = osr.SpatialReference() 
    sinu = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    modis_sinu.ImportFromProj4 (sinu)
    g = gdal.Warp('',fname, format = 'MEM',  
                    dstSRS = dstSRS, 
                    resampleAlg = 0, 
                    # srcNodata = 32767, 
                    # dstNodata=0, 
                    outputType = gdal.GDT_Float32,
                    outputBounds = outputBounds,
                    xRes = xRes, 
                    yRes = yRes, 
                    srcSRS=modis_sinu, 
                    warpOptions = ['NUM_THREADS=ALL_CPUS'])  # 
    return g.ReadAsArray() 

def read_MCD43(fnames, dstSRS, outputBounds, xRes, yRes, logger, parallel = False):
    
    par = partial(warp_data, outputBounds = outputBounds, xRes = xRes, yRes = yRes, dstSRS = dstSRS) 
    if not parallel:
        ret = []
        for _, fname in enumerate(fnames):
            if _ % 20 == 0:
                logger.info('Reading %.01f'%(_/len(fnames) * 100) + ' %...')
            ret.append(par(fname))
    else:
        procs = 6
        p = Pool(procs)
        ret = p.map(par,  fnames)
        p.close()
        p.join()
    
    n_files = int(len(fnames)/2)
    #ret = parmap(par, fnames) 
    das = np.array(ret[:n_files])     
    qas = np.array(ret[n_files:]) 
    ws  = 0.618034**qas       
    ws[qas==255]    = 0 
    das[das==32767] = 0
    return das, ws
     
def get_kernel_weights(mcd43_dir, boa_bands, obs_time, outputBounds, xRes, yRes, dstSRS, logger, temporal_filling = 16, cache_mcd43 = True, cache_name = '/mcd43_cache.npz'):
     
    qa_temp = 'MCD43_%s_BRDF_Albedo_Band_Mandatory_Quality_Band%d.vrt'
    da_temp = 'MCD43_%s_BRDF_Albedo_Parameters_Band%d.vrt'
    doy = obs_time.strftime('%Y%j')
    #print(temporal_filling)
    #if temporal_filling == True:
    #    temporal_filling = 16
    #print(temporal_filling)
    if temporal_filling:
        days   = [(obs_time - timedelta(days = int(i))) for i in np.arange(temporal_filling, 0, -1)] + \
                 [(obs_time + timedelta(days = int(i))) for i in np.arange(0, temporal_filling+1,  1)]
        #print(days)
        fnames = []
        for temp in [da_temp, qa_temp]:                                                                                    
            for day in days:
                for band in boa_bands:
                    fname = mcd43_dir + '/'.join([day.strftime('%Y-%m-%d'), temp%(day.strftime('%Y%j'), band)])
                    fnames.append(fname) 
    else:
        fnames = []
        for temp in [da_temp, qa_temp]:
            for band in boa_bands:
                    fname = mcd43_dir + '/'.join([datetime.strftime(obs_time, '%Y-%m-%d'), temp%(doy, band)])
                    fnames.append(fname)
    if not os.path.exists(cache_name):
        das, ws = read_MCD43(fnames, dstSRS, outputBounds, xRes, yRes, logger) 
    else:
        f = np.load(cache_name, allow_pickle=True)
        das = f.f.das
        ws = f.f.ws
    modis_sinu = osr.SpatialReference() 
    sinu = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    modis_sinu.ImportFromProj4 (sinu)
    mg = gdal.Warp('',fnames[0], 
                      format = 'MEM',  
                      dstSRS = dstSRS, 
                      resampleAlg = 1, 
                      xRes = xRes, 
                      yRes = yRes,
                      outputType = gdal.GDT_Float32,
                      outputBounds = outputBounds,  
                      srcSRS=sinu) 
    if cache_mcd43:
        np.savez(cache_name, 
                das = das, ws = ws, 
                fname = fnames[0], 
                xRes = xRes, 
                yRes = yRes,
                dstSRS = dstSRS,
                srcSRS=sinu,
                outputBounds = outputBounds
                )
    return das, ws, mg

def get_bounds(aoi, toa, pix_res):
    g = gdal.Warp('', toa, format = 'MEM', cutlineDSName=aoi, cropToCutline=True, xRes = pix_res, yRes = pix_res)
    dstSRS = g.GetProjectionRef()
    geo_t = g.GetGeoTransform()
    x_size, y_size = g.RasterXSize, g.RasterYSize     
    xmin, xmax = min(geo_t[0], geo_t[0] + x_size * geo_t[1]), max(geo_t[0], geo_t[0] + x_size * geo_t[1])  
    ymin, ymax = min(geo_t[3], geo_t[3] + y_size * geo_t[5]), max(geo_t[3], geo_t[3] + y_size * geo_t[5])
    outputBounds = [xmin, ymin, xmax, ymax]
    return dstSRS, outputBounds

def gaussian(xstd, ystd, norm = True):
    winx = 2*int(np.ceil(1.96*xstd))
    winy = 2*int(np.ceil(1.96*ystd))
    xgaus = signal.gaussian(winx, xstd)
    ygaus = signal.gaussian(winy, ystd)
    gaus  = np.outer(xgaus, ygaus)
    if norm:
        return gaus/gaus.sum()
    else:
        return gaus 

def smooth(da_w):
	da, w = da_w
	mid = int(da.shape[0]/2)
	if (da.shape[-1]==0) | (w.shape[-1]==0):
		return da[mid], w[mid]
	data  = np.array(smoothn(da, s=10., smoothOrder=1., axis=0, TolZ=0.001, verbose=False, isrobust=True, W = w))[[0, 3],]
	return data[0][mid], data[1][mid]


def _getKernelWeights(das, ws, mask):
    # reading MCD43 time series
    ws = ws.reshape(33, 6, -1).copy()
    f0 = das[:, 0].reshape(33, 6, -1).copy()
    f1 = das[:, 1].reshape(33, 6, -1).copy()
    f2 = das[:, 2].reshape(33, 6, -1).copy()

    w0 = np.zeros_like(f0)
    w1 = np.zeros_like(f1)
    w2 = np.zeros_like(f2)

    for i in range(6):
        ret0 = smoothn(f0[:, i][:, mask.ravel()], W = ws[:, i][:, mask.ravel()], axis=0, isrobust=True, s=0.5)
        f0[:, i][:, mask.ravel()] = ret0[0].copy()
        w0[:, i][:, mask.ravel()] = ret0[3].copy()
        
        ret1 = smoothn(f1[:, i][:, mask.ravel()], W = ws[:, i][:, mask.ravel()], axis=0, isrobust=True, s=0.5)
        f1[:, i][:, mask.ravel()] = ret1[0].copy()
        w1[:, i][:, mask.ravel()] = ret1[3].copy()
        
        ret2 = smoothn(f2[:, i][:, mask.ravel()], W = ws[:, i][:, mask.ravel()], axis=0, isrobust=True, s=0.5)
        f2[:, i][:, mask.ravel()] = ret2[0].copy()
        w2[:, i][:, mask.ravel()] = ret2[3].copy()
    
    # ind = 167
    # plt.plot(f0[:, 1, ind])
    # plt.plot(das[:, 0].reshape(33, 6, -1)[:, 1, ind])
    # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
    # plt.bar(range(33), w0[:, 1, ind]*30)

    # ind = 167
    # plt.plot(f1[:, 1, ind])
    # plt.plot(das[:, 1].reshape(33, 6, -1)[:, 1, ind])
    # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
    # plt.bar(range(33), w1[:, 1, ind]*30)

    # ind = 167
    # plt.plot(f2[:, 1, ind])
    # plt.plot(das[:, 2].reshape(33, 6, -1)[:, 1, ind])
    # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
    # plt.bar(range(33), w2[:, 1, ind]*30)

    mid = int(f0.shape[0]/2)

    f0 = f0[mid].reshape((6, ) + das.shape[-2:]) * 0.001
    f1 = f1[mid].reshape((6, ) + das.shape[-2:]) * 0.001
    f2 = f2[mid].reshape((6, ) + das.shape[-2:]) * 0.001

    w0 = w0[mid].reshape((6, ) + das.shape[-2:])
    w1 = w1[mid].reshape((6, ) + das.shape[-2:])
    w2 = w2[mid].reshape((6, ) + das.shape[-2:])

    for i in range(6):
        sret0 = smoothn(f0[i], W = w0[i], isrobust=True, s=0.1)
        w0[i] = sret0[3].copy()
        f0[i] = sret0[0].copy()
        
        sret1 = smoothn(f1[i], W = w1[i], isrobust=True, s=0.1)
        w1[i] = sret1[3].copy()
        f1[i] = sret1[0].copy()

        sret2 = smoothn(f2[i], W = w2[i], isrobust=True, s=0.1)
        w2[i] = sret2[3].copy()
        f2[i] = sret2[0].copy()

    return f0, f1, f2, w0, w1, w2

def getKernelWeights(das, ws, mask, sur):
    # reading MCD43 time series
    # shape = (33, 6,) + das.shape[-2:]
    # ws = ws.reshape(shape)
    f0 = das[:, :, 0]
    f1 = das[:, :, 1]
    f2 = das[:, :, 2]

    bad_ones = np.all(ws==0, axis=0)

    part = 10

    f0 = np.array_split(f0, part, axis=2)
    w0 = np.array_split(ws.copy(), part, axis=2)
    for i in range(10):
        ret = robust_smooth(f0[i], w0[i], 2, 0.5, axis=0)
        f0[i] = ret[0]
        w0[i] = ret[1]
    f0 = np.concatenate(f0, axis=2)
    w0 = np.concatenate(w0, axis=2)

    f1 = np.array_split(f1, part, axis=2)
    w1 = np.array_split(ws.copy(), part, axis=2)
    for i in range(10):
        ret = robust_smooth(f1[i], w1[i], 2, 0.5, axis=0)
        f1[i] = ret[0]
        w1[i] = ret[1]
    f1 = np.concatenate(f1, axis=2)
    w1 = np.concatenate(w1, axis=2)

    f2 = np.array_split(f2, part, axis=2)
    w2 = np.array_split(ws.copy(), part, axis=2)
    for i in range(10):
        ret = robust_smooth(f2[i], w2[i], 2, 0.5, axis=0)
        f2[i] = ret[0]
        w2[i] = ret[1]
    f2 = np.concatenate(f2, axis=2)
    w2 = np.concatenate(w2, axis=2)

    # import pdb;pdb.set_trace()
    sur = np.array_split(sur, part, axis=2)
    wsur = np.array_split(ws.copy(), part, axis=2)
    for i in range(10):
        ret = robust_smooth(sur[i], wsur[i], 2, 0.5, axis=0)
        sur[i] = ret[0]
        wsur[i] = ret[1]
    sur = np.concatenate(sur, axis=2)
    wsur = np.concatenate(wsur, axis=2)
    # less than 5% change is classified as stable starget
    mid = int(f0.shape[0]/2)
    stable_target  = np.std(sur[mid - 10: mid + 10], axis=0) / np.mean(sur[mid - 10: mid + 10], axis=0)
    stable_target_back_up = np.std(sur[mid - 4: mid + 4], axis=0) / np.mean(sur[mid - 4: mid + 4], axis=0)
    # cv  = np.std(sur[mid - 10: mid + 10], axis=0) / np.mean(sur[mid - 10: mid + 10], axis=0)
    # stable_target = np.all((cv < 0.05) & (cv>0), axis=0)

    # import pylab as plt
    # plt.plot(f0[:, 0, 400, 50])
    # plt.plot(f1[:, 0, 400, 50])
    # plt.plot(f2[:, 0, 400, 50])
    # f0, w0 = robust_smooth(f0, ws, 0.5, d=1)
    # f1, w1 = robust_smooth(f1, ws, 0.5, d=1)
    # f2, w2 = robust_smooth(f2, ws, 0.5, d=1)
    
    # ind = 167
    # plt.plot(f0[:, 1, ind])
    # plt.plot(das[:, 0].reshape(33, 6, -1)[:, 1, ind])
    # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
    # plt.bar(range(33), w0[:, 1, ind]*30)

    # ind = 167
    # plt.plot(f1[:, 1, ind])
    # plt.plot(das[:, 1].reshape(33, 6, -1)[:, 1, ind])
    # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
    # plt.bar(range(33), w1[:, 1, ind]*30)

    # ind = 167
    # plt.plot(f2[:, 1, ind])
    # plt.plot(das[:, 2].reshape(33, 6, -1)[:, 1, ind])
    # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
    # plt.bar(range(33), w2[:, 1, ind]*30)

    

    f0 = f0[mid]
    f1 = f1[mid]
    f2 = f2[mid]
    sur = sur[mid]

    w0 = w0[mid]
    w1 = w1[mid]
    w2 = w2[mid]
    wsur = wsur[mid]

    
    # for i in range(6):
    #     sret0 = smoothn(f0[i], W = w0[i], isrobust=True, s=0.1)
    #     w0[i] = sret0[3].copy()
    #     f0[i] = sret0[0].copy()
        
    #     sret1 = smoothn(f1[i], W = w1[i], isrobust=True, s=0.1)
    #     w1[i] = sret1[3].copy()
    #     f1[i] = sret1[0].copy()

    #     sret2 = smoothn(f2[i], W = w2[i], isrobust=True, s=0.1)
    #     w2[i] = sret2[3].copy()
    #     f2[i] = sret2[0].copy()

    #     sret3 = smoothn(sur[i], W = wsur[i], isrobust=True, s=0.1)
    #     wsur[i] = sret3[3].copy()
    #     sur[i] = sret3[0].copy()

    w0[bad_ones] = 0.00001
    w1[bad_ones] = 0.00001
    w2[bad_ones] = 0.00001
    wsur[bad_ones] = 0.00001

    return f0, f1, f2, w0, w1, w2, sur, wsur, stable_target, stable_target_back_up



# toa_dir = '/data/store01/data_dirs/students/ucfafyi/devSIAC/SIAC/SIAC/tests/S2B_MSIL1C_20171123T111349_N0206_R137_T31UCU_20171123T130706.SAFE/GRANULE/L1C_T31UCU_A003738_20171123T111347/IMG_DATA/'

# @jit(nopython=True)
# def diffMatrix(arr, d):
#     # arr = np.array(arr)
#     diff = arr[1:] - arr[:-1]
#     if d > 1:
#         diff = diffMatrix(diff, d-1)
#     return diff
# @jit(nopython=True)
# def whitsmddw(y, w, s, DTD, d=2, iteration=3):
#     W = np.diag(w)
#     wnew = w
#     for i in range(iteration):
#         left  = W + DTD
#         right = wnew * y
#         smoothed = np.linalg.solve(left, right)        
#         residuals = np.abs(smoothed - y)
#         tmp = np.sqrt(1+16*s)
#         h   = np.sqrt(1+tmp)/np.sqrt(2)/tmp
#         MAD = np.median(residuals[w>0])
#         if MAD < 0.00001:
#             MAD = 0.00001
#         u = np.abs(residuals / (1.4826*MAD) / np.sqrt(1-h))
#         wnew = (1 - (u / 4.685) ** 2)**2 * ((u/4.685)<1)
        
#         W = np.diag(wnew)
         
#     return wnew, smoothed

# @jit(nopython=True)
# def robust_smooth(array, Warray, s, d=2, iteration=3, axis=0):
#     shape = array.shape
#     sarray = np.zeros(shape)
#     warray = np.zeros(shape)

#     D = diffMatrix(np.eye(len(array)), d)
#     DTD = D.T.dot(D) * s

#     for i in range(shape[1]):
#         for j in range(shape[2]):
#             for k in range(shape[3]):
#                 y = array[:, i, j, k]
#                 w = Warray[:, i, j, k]
#                 if (y.sum() < 0.00001) | (w.sum() < 0.00001):
#                     ret = (w, y)
#                 else:
#                     ret = whitsmddw(y, w, s, DTD, d=d, iteration=iteration)
#                 warray[:, i, j, k] = ret[0]
#                 sarray[:, i, j, k] = ret[1]
#     return sarray, warray

def diffMatrix(arr, d):
    # arr = np.array(arr)
    diff = arr[1:] - arr[:-1]
    if d > 1:
        diff = diffMatrix(diff, d-1)
    return diff
def cholesky_banded_solve(left, right):
    upper  = np.diagonal(left, offset=-1, axis1=1, axis2=2)
    center = np.diagonal(left, offset= 0, axis1=1, axis2=2)

    ptsv, = get_lapack_funcs(('ptsv',), (center[0], right[0]))
    xs = [ptsv(center[i], upper[i], right[i], False, False,False)[2] for i in range(len(left))]
    xs = np.vstack(xs)
    
    return xs
def robust_smooth(array, Warray, s, d, iterations=3, axis=0):
    
    shape = list(array.shape)
    t = list(range(array.ndim))
    t.remove(axis)
    t = [axis, ] + t

    mask = np.sum(Warray, axis = 1) > 0.000001 #, axis=1) #np.all(array > 0.00001, axis=1) | 
    if mask.sum() > 0:
        array  =  array.transpose(t).reshape(shape[axis], -1).T
        Warray = Warray.transpose(t).reshape(shape[axis], -1).T
        mask = np.sum(Warray, axis = 1) > 0.000001 #, axis=1) #np.all(array > 0.00001, axis=1) | 
        left = np.zeros((mask.sum(), shape[axis], shape[axis]))
        diag_x_ind, diag_y_ind = np.arange(shape[axis]), np.arange(shape[axis])
        left[:, diag_x_ind, diag_y_ind] = Warray[mask]

        D = diffMatrix(np.eye(shape[axis]), d)
        DTD = D.T.dot(D) * s
        DTD = DTD[None, ...]

        tmp = np.sqrt(1+16*s)
        h   = np.sqrt(1+tmp)/np.sqrt(2)/tmp
        bottom = np.sqrt(1-h)

        wnew = Warray[mask]
        # iii = 86
        for i in range(iterations):
            right = array[mask] *  wnew 
            # smoothed = np.linalg.solve(left + DTD, right)
            smoothed =  cholesky_banded_solve(left + DTD, right)    
            residuals = np.abs(smoothed - array[mask])
            
            residuals = np.ma.array(residuals, mask = Warray[mask]<=0)
            MAD = np.ma.median(residuals, axis=1).data
            MAD = np.maximum(MAD, 0.00001)[...,None]
            u = residuals.data / (1.4826*MAD) / bottom
            wnew = (1 - (u / 4.685) ** 2)**2 * ((u/4.685)<1)
            left[:, diag_x_ind, diag_y_ind] = wnew
            # plt.plot(smoothed[iii], label='%d'%(i+1))

        # plt.plot(array[mask][iii], 'o-', label='orig')
        # plt.bar(range(33), Warray[mask][iii]*20)
        # plt.bar(range(33), wnew[iii]*20)

        new_shape = shape.copy()
        del shape[axis]
        new_shape = tuple([new_shape[axis], ] + shape)

        t = list(range(1, len(new_shape)))
        t.insert(axis, 0)

        array[mask] = smoothed
        Warray[mask] = wnew

        array = array.T.reshape(new_shape).transpose(t)
        Warray = Warray.T.reshape(new_shape).transpose(t)

    return array, Warray

# def getKernelWeights(das, ws, mask):
#     # reading MCD43 time series
#     shape = (33, 6,) + das.shape[-2:]
#     ws = ws.reshape(shape)
#     f0 = das[:, 0].reshape(shape)
#     f1 = das[:, 1].reshape(shape)
#     f2 = das[:, 2].reshape(shape)

#     f0, w0 = robust_smooth(f0.astype(np.float64), ws.astype(np.float64), 0.5, d=1)
#     f1, w1 = robust_smooth(f1.astype(np.float64), ws.astype(np.float64), 0.5, d=1)
#     f2, w2 = robust_smooth(f2.astype(np.float64), ws.astype(np.float64), 0.5, d=1)
    
#     # ind = 167
#     # plt.plot(f0[:, 1, ind])
#     # plt.plot(das[:, 0].reshape(33, 6, -1)[:, 1, ind])
#     # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
#     # plt.bar(range(33), w0[:, 1, ind]*30)

#     # ind = 167
#     # plt.plot(f1[:, 1, ind])
#     # plt.plot(das[:, 1].reshape(33, 6, -1)[:, 1, ind])
#     # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
#     # plt.bar(range(33), w1[:, 1, ind]*30)

#     # ind = 167
#     # plt.plot(f2[:, 1, ind])
#     # plt.plot(das[:, 2].reshape(33, 6, -1)[:, 1, ind])
#     # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
#     # plt.bar(range(33), w2[:, 1, ind]*30)

#     mid = int(f0.shape[0]/2)

#     f0 = f0[mid] * 0.001
#     f1 = f1[mid] * 0.001
#     f2 = f2[mid] * 0.001

#     w0 = w0[mid]
#     w1 = w1[mid]
#     w2 = w2[mid]

#     for i in range(6):
#         sret0 = smoothn(f0[i], W = w0[i], isrobust=True, s=0.1)
#         w0[i] = sret0[3].copy()
#         f0[i] = sret0[0].copy()
        
#         sret1 = smoothn(f1[i], W = w1[i], isrobust=True, s=0.1)
#         w1[i] = sret1[3].copy()
#         f1[i] = sret1[0].copy()

#         sret2 = smoothn(f2[i], W = w2[i], isrobust=True, s=0.1)
#         w2[i] = sret2[3].copy()
#         f2[i] = sret2[0].copy()

#     return f0, f1, f2, w0, w1, w2

# toa_dir = '/data/store01/data_dirs/students/ucfafyi/devSIAC/SIAC/SIAC/tests/S2B_MSIL1C_20171123T111349_N0206_R137_T31UCU_20171123T130706.SAFE/GRANULE/L1C_T31UCU_A003738_20171123T111347/IMG_DATA/'

# @jit(nopython=True)
# def diffMatrix(arr, d):
#     # arr = np.array(arr)
#     diff = arr[1:] - arr[:-1]
#     if d > 1:
#         diff = diffMatrix(diff, d-1)
#     return diff
# @jit(nopython=True, nogil=True)
# def whitsmddw(y, w, s, DTD, d=2, iteration=3):
#     W = np.diag(w)
#     wnew = w
#     for i in range(iteration):
#         left  = W + DTD
#         right = wnew * y
#         smoothed = np.linalg.solve(left, right)        
#         residuals = np.abs(smoothed - y)
#         tmp = np.sqrt(1+16*s)
#         h   = np.sqrt(1+tmp)/np.sqrt(2)/tmp
#         MAD = np.median(residuals[w>0])
#         if MAD < 0.00001:
#             MAD = 0.00001
#         u = np.abs(residuals / (1.4826*MAD) / np.sqrt(1-h))
#         wnew = (1 - (u / 4.685) ** 2)**2 * ((u/4.685)<1)
        
#         W = np.diag(wnew)
         
#     return wnew, smoothed

# @jit(nopython=True, nogil=True)
# def robust_smooth(array, Warray, s, d=2, iteration=3, axis=0):
#     shape = array.shape
#     sarray = np.zeros(shape)
#     warray = np.zeros(shape)

#     D = diffMatrix(np.eye(len(array)), d)
#     DTD = D.T.dot(D) * s

#     for i in range(shape[1]):
#         for j in range(shape[2]):
#             for k in range(shape[3]):
#                 y = array[:, i, j, k]
#                 w = Warray[:, i, j, k]
#                 if (y.sum() < 0.00001) | (w.sum() < 0.00001):
#                     ret = (w, y)
#                 else:
#                     ret = whitsmddw(y, w, s, DTD, d=d, iteration=iteration)
#                 warray[:, i, j, k] = ret[0]
#                 sarray[:, i, j, k] = ret[1]
#     return sarray, warray

# import numpy as np
# f = np.load('/home/users/marcyin/devSIAC/devSIAC/SIAC/SIAC/tests/LC08_L1TP_200023_20180505_20180517_01_T1/mcd43_cache.npz')
# das = f.f.das
# ws = f.f.ws

# ret = getKernelWeights(das, ws, 0)

# shape = (33, 6,) + das.shape[-2:]
# ws = ws.reshape(shape)
# f0 = das[:, 0].reshape(shape)
# f1 = das[:, 1].reshape(shape)
# f2 = das[:, 2].reshape(shape)

# f0, w0 = robust_smooth(f0.astype(np.float64), ws.astype(np.float64), 0.5, d=1)
# f1, w1 = robust_smooth(f1.astype(np.float64), ws.astype(np.float64), 0.5, d=1)
# f2, w2 = robust_smooth(f2.astype(np.float64), ws.astype(np.float64), 0.5, d=1)

# # ind = 167
# # plt.plot(f0[:, 1, ind])
# # plt.plot(das[:, 0].reshape(33, 6, -1)[:, 1, ind])
# # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
# # plt.bar(range(33), w0[:, 1, ind]*30)

# # ind = 167
# # plt.plot(f1[:, 1, ind])
# # plt.plot(das[:, 1].reshape(33, 6, -1)[:, 1, ind])
# # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
# # plt.bar(range(33), w1[:, 1, ind]*30)

# # ind = 167
# # plt.plot(f2[:, 1, ind])
# # plt.plot(das[:, 2].reshape(33, 6, -1)[:, 1, ind])
# # plt.bar(range(33), ws.reshape(33, 6, -1)[:, 1, ind]*30)
# # plt.bar(range(33), w2[:, 1, ind]*30)

# mid = int(f0.shape[0]/2)

# f0 = f0[mid] * 0.001
# f1 = f1[mid] * 0.001
# f2 = f2[mid] * 0.001

# w0 = w0[mid]
# w1 = w1[mid]
# w2 = w2[mid]

# das = das.reshape((33, 6, 3, 484, 478))
# ws = ws.reshape((33, 6,  484, 478))
# y = das[:, 0, 0, 478, 457] / 20
# w = ws[:, 0, 478, 457]

# y = '0. 0. 0. 0. 0. 0. 0. 2. 2. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'.split()
# w = '''0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 1. 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401 0.61803401'''.split()

# y = np.array(y).astype(float)
# w = np.array(w).astype(float)
# # ret = whitsmddw(y.astype(np.float64), w.astype(np.float64), 10., DTD, d=1, iteration=3)

# i = 463
# j = 460
# ret = robust_smooth(das[:, :, 0, i:470, 450:j].astype(np.float64), mask[:, :,i:470, 450:j].astype(np.float64), 1)

# %time ret = robust_smooth(das[:, :, 0, ].astype(np.float64), mask.astype(np.float64), 1)



# plt.bar(x, ret[0])
# plt.plot(x, ret[1])
# plt.plot(x, y, label = 'original')
# plt.legend()


# getKernelWeights(das, ws, mask)


#     def _get_boa(self, temporal_filling = 16):
#         self.logger.info('Reading MCD43 files.')
#         dstSRS = self.example_file.GetProjectionRef()
#         dstSRS, outputBounds = get_bounds(self.aoi, self.example_file, self.pixel_res)
#         das, ws, mg = get_kernel_weights(self.mcd43_dir, 
#                                         self.boa_bands, 
#                                         self.obs_time, 
#                                         outputBounds, 
#                                         xRes = self.aero_res*0.5, 
#                                         yRes = self.aero_res*0.5, 
#                                         dstSRS = dstSRS, 
#                                         logger = self.logger,
#                                         temporal_filling = 16,
#                                         cache_mcd43 = True,
#                                         cache_name = self.toa_dir + '/mcd43_cache.npz')
#         hg = self.example_file
        
#         geotransform = hg.GetGeoTransform()
#         h_res = geotransform[1]
#         geotransform = mg.GetGeoTransform() 
#         m_res = geotransform[1]

#         toa_mask = array_to_raster(self.bad_pix * 1., self.example_file)
#         toa_mask = reproject_data(toa_mask, mg, resample=0).data.astype(bool)

#         hx, hy = cal_psf_points(h_res, m_res, mg.RasterYSize, mg.RasterXSize)
#         hmask = (hx>=0) & (hx<hg.RasterYSize) & (hy>=0) & (hy<hg.RasterXSize)
#         hmask = hmask.reshape(toa_mask.shape)

#         mask = ~np.all(mg.ReadAsArray() == 0, axis=0) & (~toa_mask) & hmask
        
#         self.logger.info('Filling MCD43 gaps by temporal interpolation.')
#         self._annoying_angles(mg)
#         sza = np.repeat(self._sza[None, ...], len(self._vza), axis = 0)[:,mask]
#         saa = np.repeat(self._saa[None, ...], len(self._vza), axis = 0)[:,mask]
#         angles = self._vza[:,mask], sza ,self._vaa[:,mask] - saa
#         kk     = get_kk(angles)   
#         k_vol  = kk.Ross          
#         k_geo  = kk.Li  
#         kers = np.array([np.ones(k_vol.shape), k_vol, k_geo])
#         surs = []
#         No_band = len(self.toa_bands)
#         for i in range(No_band):
#             surs.append((das[i::No_band][:,:,mask] * kers[:,i] * 0.001).sum(axis=1))
#         if temporal_filling:
            
#             boa = []
#             w = []
#             self.logger.info('Filling MCD43 gaps by temporal interpolation.')
#             for i in range(No_band):
#                 das = surs[i]
#                 Ws  = ws[i::No_band][:,mask]
#                 chunks = zip(np.array_split(das, 18, axis=1), np.array_split(Ws, 18,  axis=1))
# #                 ret = parmap(smooth, chunks)
#                 ret = list(map(smooth, chunks))
#                 _b  = np.hstack([i[0] for i in ret])                   
#                 _w  = np.hstack([i[1] for i in ret]) 
#                 boa.append(_b)
#                 w.append(_w)    
#             boa = np.array(boa)
#             w = np.array(w)
#             w = np.maximum(w, 0.0001)
#             unc = 0.015 / w
#         else:
#             boa = np.array(surs)
#             ws = np.maximum(ws, 0.0001)
#             unc = 0.015 / ws
#         self.boa     = boa
#         self.boa_unc = np.minimum(unc, 0.5)
  
#         self.hx = hx[mask.ravel()]
#         self.hy = hy[mask.ravel()]



#     def _get_boa(self, temporal_filling = 16):
    
#         self.logger.info('Reading MCD43 files.')
#         dstSRS = self.example_file.GetProjectionRef()
#         dstSRS, outputBounds = get_bounds(self.aoi, self.example_file, self.pixel_res)
#         das, ws, mg = get_kernel_weights(self.mcd43_dir, 
#                                         self.boa_bands, 
#                                         self.obs_time, 
#                                         outputBounds, 
#                                         xRes = self.aero_res*0.5, 
#                                         yRes = self.aero_res*0.5, 
#                                         dstSRS = dstSRS, 
#                                         logger = self.logger,
#                                         temporal_filling = 16,
#                                         cache_mcd43 = True,
#                                         cache_name = self.toa_dir + '/mcd43_cache.npz')
#         hg = self.example_file
#         self.logger.info('Getting indexes.')
#         # self.hx, self.hy, hmask, rmask = get_index(mg, hg, self.bad_pix) #self._get_index(mg, hg)
#         self.hx, self.hy, hmask, rmask = get_index(mg, hg, self.bad_pix) #self._get_index(mg, hg)
#         No_band = len(self.toa_bands)
#         mask = ~(mg.ReadAsArray()[0]==0)
#         self._annoying_angles(mg)
#         sza = np.repeat(self._sza[None, ...], len(self._vza), axis = 0)[:,mask][:, hmask][:, rmask]
#         saa = np.repeat(self._saa[None, ...], len(self._vza), axis = 0)[:,mask][:, hmask][:, rmask]
#         angles = self._vza[:,mask][:, hmask][:, rmask], sza ,self._vaa[:,mask][:, hmask][:, rmask] - saa
#         kk     = get_kk(angles)   
#         k_vol  = kk.Ross          
#         k_geo  = kk.Li  
#         kers = np.array([np.ones(k_vol.shape), k_vol, k_geo])
#         surs = []
#         for i in range(No_band):
#             surs.append((das[i::No_band][:,:,mask][:,:,hmask][:,:,rmask] * kers[:,i] * 0.001).sum(axis=1))
#         if temporal_filling:
            
#             boa = []
#             w = []
#             self.logger.info('Filling MCD43 gaps by temporal interpolation.')
#             for i in range(No_band):
#                 das = surs[i]
#                 Ws  = ws[i::No_band][:,mask][:, hmask][:, rmask]
#                 chunks = zip(np.array_split(das, 18, axis=1), np.array_split(Ws, 18,  axis=1))
# #                 ret = parmap(smooth, chunks)
#                 ret = list(map(smooth, chunks))
#                 _b  = np.hstack([i[0] for i in ret])                   
#                 _w  = np.hstack([i[1] for i in ret]) 
#                 boa.append(_b)
#                 w.append(_w)    
#             boa = np.array(boa)
#             w = np.array(w)
#             unc = 0.015 / w
#         else:
#             boa = np.array(surs)
#             unc = 0.015 / ws
#         self.boa     = boa
#         self.boa_unc = np.minimum(unc, 0.5)


