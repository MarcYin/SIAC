import os
import kernels
import datetime
import numpy as np
# import pynndescent
from glob import glob
from osgeo import gdal, osr
from pyproj import Proj, transform, CRS
from SIAC.read_MCD43 import robust_smooth
from SIAC.reproject import array_to_raster, reproject_data

from SIAC.smoothn import smoothn

import os
import logging

from copy import deepcopy

# import numba

# @numba.njit(
#     [
#         "f4(f4[::1],f4[::1])",
#         numba.types.float32(
#             numba.types.Array(numba.types.float32, 1, "C", readonly=True),
#             numba.types.Array(numba.types.float32, 1, "C", readonly=True),
#         ),
#     ],
#     fastmath=True,
#     locals={
#         "result": numba.types.float32,
#         "diff": numba.types.float32,
#         "dim": numba.types.uint32,
#         "i": numba.types.uint16,
#     },
# )
# def squared_euclidean(x, y):
#     r"""Squared euclidean distance.

#     .. math::
#         D(x, y) = \sum_i (x_i - y_i)^2
#     """
#     result = 0.0
#     dim = x.shape[0]
#     for i in range(dim):
#         diff = x[i] - y[i]
#         result += diff * diff

#     return result #/ 30.

def create_logger(fname = None):
    logger = logging.getLogger('SIAC NBAR')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    if not logger.handlers:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    if fname is not None:
        fh = logging.FileHandler(fname)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)    
        logger.addHandler(fh)
    return logger



def solar_geometry(year, month, day, hourGMT, lat, lon):
    '''
    Solar geometry unsing subsolar point and atan2: https://doi.org/10.1016/j.renene.2021.03.047
    
    Parameters:
    -----------
    year: scalar
        4 digit year (1950 to 2050), e.g. 2021
    month: scalar
        month of year, e.g. 1-12
    day: scalar
        day of month, e.g. 1-31
    hourGMT: scalar
        hour of the day, 15.51
    lat: scalar or nD array
        latitude in degrees, positive for northern hemisphere
    lon: scalar or nD array
        longitude in degrees, positve for East hemisphere
    
    Return:
    -------
    sza: same shape as lat and lon
        solar zenith angle in degrees
    vza: same shape as lat and lon
        solar azimuth angle in degrees
    '''
    
    
    if (year % 4 == 0) | (year % 100 == 0) | (((year % 100) & ((year % 400)))==0):
        nday = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
        nday = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    julday = np.cumsum([0, ] + nday)
    
    dyear = year - 2000
    
    dayofyr = julday[month - 1] + day

    if dyear <= 0:
        xleap = int(dyear / 4)
    elif dyear > 0:
        if dyear % 4 ==0:
            xleap = int(dyear / 4)
        else:
            xleap = int(dyear / 4) + 1
    
    n = 1.5 + dyear * 365 + xleap * 1 + dayofyr + hourGMT / 24
    L = abs((280.466 + 0.9856474 * n) % 360)
    g = abs((357.528 + 0.9856003 * n) % 360)
    lambda_ = abs((L + 1.915 * np.sin(np.deg2rad(g)) + 0.020 * np.sin(2 * np.deg2rad(g))) % 360)
    epsilon = (23.439 - 0.0000004 * n)
    
    alpha = abs(np.arctan2(np.cos(np.deg2rad(epsilon)) * np.sin(np.deg2rad(lambda_)),  np.cos(np.deg2rad(lambda_))) * 180 / np.pi % 360)
    delta = np.arcsin(np.sin(np.deg2rad(epsilon)) * np.sin(np.deg2rad(lambda_))) * 180 / np.pi
    
    R = 1.00014 - 0.01671 * np.cos(np.deg2rad(g)) - 0.00014 * np.cos(2 * np.deg2rad(g))
    EoT = abs(((L - alpha) + 180) % 360) - 180
    
    sunlat = delta
    sunlon = -15.0 * (hourGMT - 12 + EoT * 4 / 60)
    PHIo = np.deg2rad(lat)
    PHIs = np.deg2rad(sunlat)
    LAMo = np.deg2rad(lon)
    LAMs = np.deg2rad(sunlon)
    
    Sx = np.cos(PHIs) * np.sin(LAMs - LAMo)
    Sy = np.cos(PHIo) * np.sin(PHIs) - np.sin(PHIo) * np.cos(PHIs) * np.cos(LAMs - LAMo)
    Sz = np.sin(PHIo) * np.sin(PHIs) + np.cos(PHIo) * np.cos(PHIs) * np.cos(LAMs - LAMo)
    sza = np.rad2deg(np.arccos(Sz))
    saa = np.rad2deg(np.arctan2(-Sx, -Sy))
    return sza, saa

def get_kk(angles):

    vza ,sza,raa = angles
    kk = kernels.Kernels(vza ,sza,raa,\
                         RossHS=False,MODISSPARSE=True,\
                         RecipFlag=True,normalise=1,\
                         doIntegrals=False,LiType='Sparse',RossType='Thick')
    return kk

def getKernelWeights(toa_dir):
    # reading MCD43 time series
    
    f = np.load(toa_dir + '/kernel_weights.npz')
    fs = f.f.fs
    uncs = f.f.uncs
    uncs = np.minimum(uncs, 1)
    # f = np.load(toa_dir + '/mcd43_cache.npz')
    # modis_sinu = osr.SpatialReference() 
    # sinu = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
    # modis_sinu.ImportFromProj4 (sinu)

    # dstSRS = osr.SpatialReference() 
    # dstSRS.ImportFromProj4 (str(f.f.dstSRS))

        # mg = gdal.Warp('',fname, format = 'MEM',  dstSRS = dstSRS, resampleAlg = 0, outputType = gdal.GDT_Float32,
        #               cutlineDSName=aoi, cropToCutline=True, xRes = xRes, yRes = yRes, srcSRS=modis_sinu, ) 

    # mg = gdal.Warp('', str(f.f.fname), format = 'MEM', dstSRS = dstSRS, xRes = float(f.f.xRes), yRes = \
    #                 float(f.f.yRes), outputBounds=f.f.outputBounds.tolist(), resampleAlg = 0, srcSRS=modis_sinu, )    
    f0, f1, f2 = fs
    f0_unc, f1_unc, f2_unc = uncs
    # f0_unc, f1_unc, f2_unc = 0.3 * fs
    # for i in range(6):
    #     # plt.figure()
    #     # plt.imshow(f0[i]);plt.colorbar()
    #     sret0 = smoothn(f0[i], isrobust=True, s=0.5)
    #     w = sret0[3]
    #     mask = w < 0.2
    #     # plt.figure()
    #     # plt.imshow(mask);plt.colorbar()
    #     f0[i][mask] = sret0[0][mask]
    #     f0[i][:] = sret0[0][:]
    #     # plt.figure()
    #     # plt.imshow(f0[i]);plt.colorbar()
    #     # plt.figure()
    #     # plt.imshow(sret[0]);plt.colorbar()
        
    #     sret1 = smoothn(f1[i], isrobust=True, s=0.5)
    #     w = sret1[3]
    #     mask = w < 0.2
    #     f1[i][mask] = sret1[0][mask]
    #     f1[i][:] = sret1[0][:]

    #     sret2 = smoothn(f2[i], isrobust=True, s=0.5)
    #     w = sret2[3]
    #     mask = w < 0.2
    #     f2[i][mask] = sret2[0][mask]
    #     f2[i][:] = sret2[0][:]
        
    return f0, f1, f2, f0_unc, f1_unc, f2_unc

def interp_fs(obs_cwv, modis_cwv, f0, f1, f2, f0_unc, f1_unc, f2_unc):
    lr_wvs = []
    for wv in obs_cwv:
        if wv < 469:
            lr_wvs.append([modis_cwv[0], modis_cwv[0]])
        elif (wv > 858) & (wv <= 1000):
            lr_wvs.append([modis_cwv[3], modis_cwv[3]])
        elif (wv > 645) & (wv<700):
            lr_wvs.append([modis_cwv[2], modis_cwv[2]])
        elif (wv > 645) & (wv <= 858):
            lr_wvs.append([modis_cwv[2], modis_cwv[3]])
        elif (wv > 1580) & (wv <= 1740):
            lr_wvs.append([modis_cwv[4], modis_cwv[4]])
        elif wv > 2130:
            lr_wvs.append([modis_cwv[5], modis_cwv[5]])
        else:
            wv_diff = modis_cwv - wv
            ind = np.argsort(abs(wv_diff)).astype(int)
            lr_wv = modis_cwv[ind[:2]]
            lr_wvs.append(lr_wv)

    interped_f0 = []
    interped_f1 = []
    interped_f2 = []
    interped_f0_unc = []
    interped_f1_unc = []
    interped_f2_unc = []

    for i in range(len(obs_cwv)):
        lr_wv = lr_wvs[i]
        wv = obs_cwv[i]
        l_modis_ind = modis_cwv.tolist().index(lr_wv[0])
        r_modis_ind = modis_cwv.tolist().index(lr_wv[1])
        
        if l_modis_ind == r_modis_ind :
            f0b = f0[l_modis_ind]
            f1b = f1[l_modis_ind]
            f2b = f2[l_modis_ind]
            # increase uncertainty by 20% for MODIS bands that are not interpolated
            f0b_unc = f0_unc[l_modis_ind] * 1.2
            f1b_unc = f1_unc[l_modis_ind] * 1.2
            f2b_unc = f2_unc[l_modis_ind] * 1.2

        else:
            slope = (wv - lr_wv[0]) / (lr_wv[1] - lr_wv[0])
            
            f0b = (f0[r_modis_ind] - f0[l_modis_ind]) * slope + f0[l_modis_ind]
            f1b = (f1[r_modis_ind] - f1[l_modis_ind]) * slope + f1[l_modis_ind]
            f2b = (f2[r_modis_ind] - f2[l_modis_ind]) * slope + f2[l_modis_ind]

            f0b_unc =  np.sqrt(slope**2 * f0_unc[l_modis_ind]**2 + (1-slope)**2 * f0_unc[r_modis_ind]**2)
            f1b_unc =  np.sqrt(slope**2 * f1_unc[l_modis_ind]**2 + (1-slope)**2 * f1_unc[r_modis_ind]**2)
            f2b_unc =  np.sqrt(slope**2 * f2_unc[l_modis_ind]**2 + (1-slope)**2 * f2_unc[r_modis_ind]**2)



            # slope = (f0[r_modis_ind] - f0[l_modis_ind]) / (lr_wv[1] - lr_wv[0])
            # f0b = (wv - lr_wv[0]) * slope + f0[l_modis_ind]

            # slope = (f1[r_modis_ind] - f1[l_modis_ind]) / (lr_wv[1] - lr_wv[0])
            # f1b = (wv - lr_wv[0]) * slope + f1[l_modis_ind]
            
            # slope = (f2[r_modis_ind] - f2[l_modis_ind]) / (lr_wv[1] - lr_wv[0])
            # f2b = (wv - lr_wv[0]) * slope + f2[l_modis_ind]


            # slope = (wv - modis_cwv[l_modis_ind]) / (modis_cwv[r_modis_ind] - modis_cwv[l_modis_ind] )
            # c_factor = (c_factor_bands[r_modis_ind] - c_factor_bands[l_modis_ind]) * slope + c_factor_bands[l_modis_ind]
            # c_factor_std = np.sqrt(slope**2 * c_factor_bands_std[l_modis_ind]**2 + (1-slope)**2 * c_factor_bands_std[r_modis_ind]**2)


        interped_f0.append(f0b)
        interped_f1.append(f1b)
        interped_f2.append(f2b)
        interped_f0_unc.append(f0b_unc)
        interped_f1_unc.append(f1b_unc)
        interped_f2_unc.append(f2b_unc)
    
    interped_f0 = np.array(interped_f0)
    interped_f1 = np.array(interped_f1)
    interped_f2 = np.array(interped_f2)
    interped_f0_unc = np.array(interped_f0_unc)
    interped_f1_unc = np.array(interped_f1_unc)
    interped_f2_unc = np.array(interped_f2_unc)


    return interped_f0, interped_f1, interped_f2, interped_f0_unc, interped_f1_unc, interped_f2_unc

def save_band(array, outputFileName, projectionRef, geotransform):
    '''
    This function is used to save the array as a geotiff file.
    Parameters:
        array: 2D array
        outputFileName: output file name
        projectionRef: projection reference
        geotransform: geotransform
    '''
    nx, ny = array.shape
    if os.path.exists(outputFileName):
        os.remove(outputFileName)
    dst_ds = gdal.GetDriverByName('GTiff').Create(outputFileName, ny, nx, 1, gdal.GDT_Int16, options=["TILED=YES", "COMPRESS=DEFLATE"])
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(projectionRef)
    array = array
    array[~(array>=0)] = -9999
    array = array.astype(np.int16)
    dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
    dst_ds.GetRasterBand(1).WriteArray(array)
    dst_ds.FlushCache()                  
    dst_ds = None


def gap_filling_kernel_weights(das, ws):
    """
    This function is used to calculate the weights for the gap filling kernel.
    The weights are calculated based on the distance between the pixel and its
    neighbors. The weights are normalized so that the sum of the weights is 1.
    The weights are calculated for each band separately.
    Parameters
    ----------
    das : numpy array
        The array of the data to be smoothed. The shape of the array is
        (nrow, ncol, nband).
    ws : numpy array
        The array of the weights
    Returns 
    -------
    fs : numpy array
        The array of the weights for the middle of the time series of das
    ws : numpy array
        The array of the weights for the middle of the time series of ws
    """
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

    mid = int(f0.shape[0]/2)

    # # import pdb;pdb.set_trace()
    # sur = np.array_split(sur, part, axis=2)
    # wsur = np.array_split(ws.copy(), part, axis=2)
    # for i in range(10):
    #     ret = robust_smooth(sur[i], wsur[i], 2, 0.5, axis=0)
    #     sur[i] = ret[0]
    #     wsur[i] = ret[1]
    # sur = np.concatenate(sur, axis=2)
    # wsur = np.concatenate(wsur, axis=2)
    # # less than 5% change is classified as stable starget
    # mid = int(f0.shape[0]/2)
    # stable_target  = np.std(sur[mid - 10: mid + 10], axis=0) / np.mean(sur[mid - 10: mid + 10], axis=0)
    # stable_target_back_up = np.std(sur[mid - 4: mid + 4], axis=0) / np.mean(sur[mid - 4: mid + 4], axis=0)

    f0 = f0[mid]
    f1 = f1[mid]
    f2 = f2[mid]
    # sur = sur[mid]

    w0 = w0[mid]
    w1 = w1[mid]
    w2 = w2[mid]
    # wsur = wsur[mid]


    w0[bad_ones] = 0.00001
    w1[bad_ones] = 0.00001
    w2[bad_ones] = 0.00001
    # wsur[bad_ones] = 0.00001
    
    w0 = np.maximum(w0, 0.0001)
    w1 = np.maximum(w1, 0.0001)
    w2 = np.maximum(w2, 0.0001)
    
    f0_unc = 0.015 / w0
    f1_unc = 0.015 / w1
    f2_unc = 0.015 / w2

    uncs = np.array([f0_unc, f1_unc, f2_unc])
    fs = np.array([f0, f1, f2])

    return fs, uncs #, sur, wsur, stable_target, stable_target_back_up

def get_BRDF_product(s2_file_dir, temporal_window = 16, Gee = True, use_VIIRS = False, vnp43_folder=None, logger=None):
    
    """
    Compute the BRDF product for a Sentinel-2 L1C image for a given temporal window (default 16 days)
    You can use either GEE or VIIRS to download the MCD43 product
    The output is a 3D array with the following bands:
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
    The output is saved in the same folder as the Sentinel-2 L1C data as a .npz file
    agrs:
        s2_file_dir: str, path to the folder of Sentinel-2 L1C data
        temporal_window: int, number of days to use for the temporal window
        Gee: bool, if True, use GEE to download the MCD43 product
        use_VIIRS: bool, if True, use VIIRS to download the VNP43 product
        vnp43_folder: str, path to the folder of VNP43 product
        logger: logging object, if None, create a new logger
    return:
        None

    """    
    if logger is None:
        logger = create_logger()
        print = logger.info
    else:
        print = logger.info

    obs_time = datetime.datetime.strptime(s2_file_dir.split('_MSIL1C_')[1][:15], '%Y%m%dT%H%M%S')
    mcd43_date = datetime.datetime(obs_time.year, obs_time.month, obs_time.day)

    B2 = glob(s2_file_dir + '/GRANULE/*/IMG_DATA/*B02_sur.tif')[0]
    example_file = gdal.Open(B2)

    dstSRS = example_file.GetProjectionRef()
    
    geo_t = example_file.GetGeoTransform()
    x_size, y_size = example_file.RasterXSize, example_file.RasterYSize     
    xmin, xmax = min(geo_t[0], geo_t[0] + x_size * geo_t[1]), max(geo_t[0], geo_t[0] + x_size * geo_t[1])  
    ymin, ymax = min(geo_t[3], geo_t[3] + y_size * geo_t[5]), max(geo_t[3], geo_t[3] + y_size * geo_t[5])
    outputBounds = [xmin, ymin, xmax, ymax]


    if Gee:
        import json
        from SIAC.MCD43_GEE import get_MCD43_GEE
        from SIAC.read_mcd43GEE import read_MCD3_GEE
        from SIAC.raster_boundary import get_boundary
        
        print('Getting MCD43A1 from GEE')
        geojson = get_boundary(B2, to_wgs84 = True)[0]
        coords = json.loads(geojson)['features'][0]['geometry']['coordinates']
        mcd43_gee_folder = os.path.dirname(os.path.dirname(B2)) + '/MCD43/'
        if not os.path.exists(mcd43_gee_folder):
            os.mkdir(mcd43_gee_folder)
        get_MCD43_GEE(mcd43_date, coords, temporal_window, mcd43_gee_folder)

        print('Reading MCD43A1 files.')
        boa_bands = [3, 4, 1, 2, 6, 7]
        das, ws, mg = read_MCD3_GEE(mcd43_gee_folder, boa_bands, outputBounds,  500, 500, dstSRS, logger)
        
    elif use_VIIRS:
        from SIAC.get_VNP43MA1 import download_VNP43MA1
        from SIAC.read_VNP43MA1 import read_VNP43MA1
        
        if vnp43_folder is None:
            vnp43_folder = os.path.dirname(os.path.dirname(B2)) + '/VNP43/'
            if not os.path.exists(vnp43_folder):
                os.mkdir(vnp43_folder)

        print('Getting VNP43MA1 from NASA')
        filenames = download_VNP43MA1(B2, mcd43_date, vnp43_folder, temporal_window = temporal_window)
        filenames = np.array(filenames)
        all_dates = np.array([i.split('/')[-1] .split('.')[1][1:9] for i in filenames])          
        udates = np.unique(all_dates)  
        VNP43_fnames_dates =  [[filenames[all_dates==date].tolist(),date] for date in udates]

        viirs_bands = ['M3', 'M4', 'M5', 'M7', 'M10', 'M11']
        modis_bands = [3, 4, 1, 2, 6, 7]
        modis_to_viirs = dict(zip(modis_bands, viirs_bands))
        bands = [modis_to_viirs[i] for i in modis_bands]

        print('Reading VNP43MA1 files.')
        das, ws, mg = read_VNP43MA1(VNP43_fnames_dates, bands, outputBounds, 500, 500, dstSRS, logger)
        

    else:
        # vrt_dir = get_mcd43(B2, mcd43_date, mcd43_dir = mcd43, vrt_dir = vrt_dir, logger = logger, jasmin = False)
        # mcd43_gee_folder = None
        print('Only GEE and VIIRS are supported for now.')
        return

    toa_dir = os.path.dirname(B2)
    np.savez(toa_dir + '/kernel_weights_time_series.npz', das = das.astype(int), ws = ws)
    print('Gap filling kernel weights.')
    das = das * 0.001
    fs, uncs = gap_filling_kernel_weights(das, ws)
    print('Saving kernel weights.')
    np.savez(toa_dir + '/kernel_weights.npz', fs = fs, uncs = uncs)


def create_nbar(s2_file_dir, nbar_sza='atan2', logger=None, mosaic_start_date=None, mosaic_end_date=None, mosaic_hour = None, Gee = True, use_VIIRS = False, vnp43_folder=None, temporal_window = 16):
    """
    Create NBAR from Sentinel-2 data.
    
    Solar zenith angle `nbar_sza`:
    
    It can be either the mean of the SZA from S2 ('use_s2') 
    Or the SZA at the subsolar point ('atan2', default) from https://doi.org/10.1016/j.renene.2021.03.047.
    If mosaic over a period of time, use the mean SZA for the whole period calculated
    from the SZA with the subsolar point ('temporal_average_sza').
    Or to any user defined sza (float number between 0-60 is suggested).
    
    Args:
        s2_file_dir (str): path to the Sentinel-2 data.
        nbar_sza (str, optional): nbar sza. Defaults to 'atan2'.
        logger (logging, optional): logger for function. Defaults to None.
        mosaic_start_date (datetime, optional): mosaic starting date. Defaults to None.
        mosaic_end_date (datetime, optional): mosaic ending date. Defaults to None.
        mosaic_hour (float, optional): float hour between 0-24. Defaults to None.
        Gee (bool, optional): Whether to use GEE or not. Defaults to True.
        use_VIIRS (bool, optional): whether to use VIIRS. Defaults to False.
        vnp43_folder (str, optional): path for saving VIIRS data. Defaults to None.
        temporal_window (int, optional): days before and after the obervation date. Defaults to 16.
    return:
        None
    """    

    if logger is None:
        logger = create_logger()
        print = logger.info
    else:
        print = logger.info

    B2 = glob(s2_file_dir + '/GRANULE/*/IMG_DATA/*B02_sur.tif')
    if len(B2) == 0:
        print('No S2 surface reflectance found.')
        return
    else:
        B2 = B2[0]
    
    toa_dir = os.path.dirname(B2)
    ang_dir = os.path.dirname(toa_dir) + '/ANG_DATA/'
    
    kernel_weights_fname = toa_dir + '/kernel_weights.npz'
    if not os.path.exists(kernel_weights_fname):
        print('No kernel weights file found, creating one.')
        get_BRDF_product(s2_file_dir, Gee = Gee, use_VIIRS = use_VIIRS, vnp43_folder=vnp43_folder, logger=logger, temporal_window = temporal_window)

    s2_cwv = [442.7, 492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 945.1, 1373.5, 1613.7, 2202.4]
    modis_cwv = [469, 555, 645, 858, 1640, 2130]
    modis_cwv = np.array(modis_cwv)
    f0, f1, f2, f0_unc, f1_unc, f2_unc = getKernelWeights(toa_dir)
    interped_f0, interped_f1, interped_f2, interped_f0_unc, interped_f1_unc, interped_f2_unc = interp_fs(s2_cwv, modis_cwv, f0, f1, f2, f0_unc, f1_unc, f2_unc)
    
    kernel_g = array_to_raster(interped_f0[0], B2)

    saa_sza = ang_dir + 'SAA_SZA.tif'
    saa, sza = reproject_data(saa_sza, kernel_g, resample=1).data / 100

    # determine the solar zenith angle
    # it can be either the mean of the SZA from S2 or the SZA at the subsolar point (atan2)
    # if mosaic over a period of time, use the mean SZA for the whole period calculated
    # from the SAA and SZA at the center of the mosaic with the subsolar point
    if nbar_sza == 'atan2':
        # print(sza, vza, saa, vaa)
        print('Using solar zenith angle based on subsolar point and atan2 function')
        g = gdal.Open(B2)
        geom = g.GetGeoTransform()
        x = np.repeat(np.arange(saa.shape[0]), saa.shape[1])
        y = np.tile(np.arange(saa.shape[0]), saa.shape[1])

        geom = kernel_g.GetGeoTransform()
        x, y = geom[0] + geom[1] * x, geom[3] + geom[5] * y
        
        proj = osr.SpatialReference(wkt=g.GetProjection())

        inProj = CRS.from_proj4(proj.ExportToProj4())

        # epsg = proj.GetAttrValue('AUTHORITY',1)
        # inProj = Proj(init='epsg:%s'%epsg)
        
        outProj = Proj(init='epsg:4326')
        lon, lat = transform(inProj,outProj, x, y)

        # print(lon, lat)
        
        date = datetime.datetime.strptime(s2_file_dir.split('_MSIL1C_')[1][:15], '%Y%m%dT%H%M%S')
        # print(date.year, date.month, date.day, date.hour + date.minute/60 + date.second/3600 , lat.mean(), lon.mean())
        sza_avg, saa_avg = solar_geometry(date.year, date.month, date.day, date.hour + date.minute/60 + date.second/3600 , lat, lon)

        sza_avg = sza_avg.reshape(sza.shape)
        # print(sza_avg.mean(), sza.mean())
        fname_postfix = 'sza_avg'
    if nbar_sza == 'temporal_average_sza':
        g = gdal.Open(B2)
        geom = g.GetGeoTransform()
        x = np.repeat(np.arange(saa.shape[0]), saa.shape[1])
        y = np.tile(np.arange(saa.shape[0]), saa.shape[1])

        geom = kernel_g.GetGeoTransform()
        x, y = geom[0] + geom[1] * x, geom[3] + geom[5] * y
        
        proj = osr.SpatialReference(wkt=g.GetProjection())

        inProj = CRS.from_proj4(proj.ExportToProj4())

        # epsg = proj.GetAttrValue('AUTHORITY',1)
        # inProj = Proj(init='epsg:%s'%epsg)
        
        outProj = Proj(init='epsg:4326')
        lon, lat = transform(inProj,outProj, x, y)
        if (mosaic_start_date is None) | (mosaic_end_date is None):
            print('mosaic_start_date and mosaic_end_date should be provided for temporal_average_sza')
            return
        else:
            print('Doing NBAR with temporal average of SZA between %s and %s'%(mosaic_start_date, mosaic_end_date))
        date = datetime.datetime.strptime(s2_file_dir.split('_MSIL1C_')[1][:15], '%Y%m%dT%H%M%S')
        if mosaic_hour is None:
            # Sentinel-2 orbits are sun-synchronous, with a period of 98 minutes.
            # The satellite Local Time at Descending Node (LTAN) is 10:30 UTC.
            # convert to GMT time based on lat lon of the image
            # Mean Local Solar Time (MLST)
            S2_local_time = 10.5
            GMT_hour = S2_local_time - lon.mean() / 180 * 12

            image_hour = date.hour + date.minute/60 + date.second/3600 

            print('GMT hour from mean solar time: %.02f'%GMT_hour)
            print('GMT hour from image: %.02f'%image_hour)
            
            mosaic_hour = GMT_hour
            print('mosaic_hour is not provided, use GMT hour calculated from mean soloar time: GMT hour %.01f'%mosaic_hour)
            
        
        # print(date.year, date.month, date.day, date.hour + date.minute/60 + date.second/3600 , lat.mean(), lon.mean())
        doy_gap = (mosaic_end_date - mosaic_start_date).days

        # sza_avg = 0
        # # saa_avg = 0
        # for doy in range(doy_gap):
        #     mosaic_date = mosaic_start_date + datetime.timedelta(days=doy)
        #     mosaic_sza, mosaic_saa = solar_geometry(mosaic_date.year, mosaic_date.month, mosaic_date.day, mosaic_hour, lat, lon)
        #     mosaic_sza = mosaic_sza.reshape(sza.shape)
        #     sza_avg += mosaic_sza
        #         # mosaic_saa_avg += mosaic_saa
        # sza_avg /= doy_gap
        # # saa_avg /= doy_gap
        half_doy_gap = int(doy_gap/2)
        mosaic_date = mosaic_start_date + datetime.timedelta(days=half_doy_gap)
        mosaic_sza, mosaic_saa = solar_geometry(mosaic_date.year, mosaic_date.month, mosaic_date.day, mosaic_hour, lat, lon)
        mosaic_sza = mosaic_sza.reshape(sza.shape)
        sza_avg = mosaic_sza


        fname_postfix = 'sza_temporal_avg'
    elif nbar_sza == 'use_s2':
        sza_avg = sza
        fname_postfix = 'sza_s2'
    elif isinstance(nbar_sza, float):
        sza_avg = np.ones_like(sza) * nbar_sza
        fname_postfix = 'sza_%02d'%int(nbar_sza)
    else:
        print('nbar_sza should be one of the following: use_s2, temporal_average_sza, or a float number between 0 and 60')
        return

    s2_bands = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']
    band_inds = [1, 2, 3, 4, 5, 6, 7, 8, 11, 12]

    for ind in band_inds:
        band = s2_bands[ind]
        print('-----------------')
        print(band)
        band_file = B2.replace('B02_sur.tif', band + '_sur.tif')

        vaa_vza = ang_dir + 'VAA_VZA_%s.tif'%band

        vaa, vza = reproject_data(vaa_vza, kernel_g, resample=1).data / 100

        raa = (vaa - saa) % 360
        
        kk = get_kk([vza, sza, raa])
        k_vol  = kk.Ross          
        k_geo  = kk.Li  
        kers = np.array([np.ones(k_vol.shape), k_vol, k_geo])
        rho_mod_s2_ang = kers * np.array([interped_f0[ind], interped_f1[ind], interped_f2[ind]])

        rho_mod_s2_ang = rho_mod_s2_ang.sum(axis=0)

        rho_mod_s2_ang_unc = np.sqrt(np.sum(kers**2 * np.array([interped_f0_unc[ind], interped_f1_unc[ind], interped_f2_unc[ind]])**2, axis=0))

        kk = get_kk([vza*0, np.ones_like(sza) * sza_avg, raa*0])
        k_vol  = kk.Ross          
        k_geo  = kk.Li  
        kers = np.array([np.ones(k_vol.shape), k_vol, k_geo])
        rho_mod_s2_nadir = kers * np.array([interped_f0[ind], interped_f1[ind], interped_f2[ind]])
        rho_mod_s2_nadir = rho_mod_s2_nadir.sum(axis=0)
        
        rho_mod_s2_nadir_unc = np.sqrt(np.sum(kers**2 * np.array([interped_f0_unc[ind], interped_f1_unc[ind], interped_f2_unc[ind]])**2, axis=0))

        rho_mod_s2_nadir[rho_mod_s2_nadir<=0] = np.nan
        rho_mod_s2_nadir_unc[rho_mod_s2_nadir<=0] = np.nan

        rho_mod_s2_ang[rho_mod_s2_ang<=0] = np.nan
        rho_mod_s2_ang_unc[rho_mod_s2_ang<=0] = np.nan

        c_factor = rho_mod_s2_nadir / rho_mod_s2_ang
        
        # https://www.phenix.bnl.gov/WWW/publish/elke/EIC/Files-for-Wiki/lara.02-008.errors.pdf
        # rho_mod_s2_nadir and rho_mod_s2_ang are highly correlated (r>0.99), 
        # so we use the formula for correlated uncertainties from the above link
        # sigma_ab = np.cov(rho_mod_s2_nadir.ravel(), rho_mod_s2_ang.ravel())[0,1]
        # c_factor_unc = np.abs(c_factor) * np.sqrt((rho_mod_s2_nadir_unc / rho_mod_s2_nadir)**2 + (rho_mod_s2_ang_unc / rho_mod_s2_ang)**2  - 2 * sigma_ab / (rho_mod_s2_nadir * rho_mod_s2_ang))
        valid_mask = np.isfinite(rho_mod_s2_nadir) & np.isfinite(rho_mod_s2_ang)
        corrcoef = np.corrcoef(rho_mod_s2_nadir[valid_mask], rho_mod_s2_ang[valid_mask])[0,1]
        # covariance not available, so using correlation coefficient of between rho_mod_s2_nadir and rho_mod_s2_ang
        # sigma_rho_mod_s2_nadir_rho_mod_s2_ang = corrcoef * rho_mod_s2_nadir_unc * rho_mod_s2_ang_unc
        c_factor_unc = np.abs(c_factor) * np.sqrt((rho_mod_s2_nadir_unc / rho_mod_s2_nadir)**2 + (rho_mod_s2_ang_unc / rho_mod_s2_ang)**2  - 2 * rho_mod_s2_nadir_unc * rho_mod_s2_ang_unc * corrcoef / (rho_mod_s2_nadir * rho_mod_s2_ang))
        
        max_c_factor = 2
        max_c_factor_unc = 2
        good_pixels = (c_factor > 0) & (c_factor <= max_c_factor) & np.isfinite(c_factor_unc) & (c_factor_unc < max_c_factor_unc)
        
        c_factor[~good_pixels] = 1
        c_factor_unc[~good_pixels] = 2

        # import pdb; pdb.set_trace()
        w = 1 / c_factor_unc**2
        w[~good_pixels] = 0
        
        
        arr = deepcopy(c_factor)
        ret = smoothn(arr, isrobust=True, s=0.1, W = w)
        c_factor = ret[0]
        c_factor_unc = 1 / np.sqrt(ret[3]) * c_factor_unc
        c_factor_unc = c_factor_unc.clip(0, 2)

        print('Mean c factor, vza, sza, raa, nbar_sza: %.03f, %.03f, %.03f, %.03f, %.03f'%(np.nanmean(c_factor), vza.mean(), sza.mean(), raa.mean(), sza_avg.mean()))


        # data = gdal.Open(band_file)
        geom = gdal.Open(band_file).GetGeoTransform()
        projectionRef = gdal.Open(band_file).GetProjectionRef()

        xRes = abs(geom[1])
        yRes = abs(geom[5])

        c_factor = array_to_raster(c_factor, B2)
        c_factor_unc = array_to_raster(c_factor_unc, B2)
        
        c_factor_fname = band_file.replace('_sur', '_cfactor_%s'%(fname_postfix))
        c_factor_unc_fname = band_file.replace('_sur', '_cfactor_unc_%s'%(fname_postfix))
        print('save band')
        # print(c_factor_fname)
        save_band(c_factor.ReadAsArray() * 10000, c_factor_fname, c_factor.GetProjectionRef(), c_factor.GetGeoTransform())
        save_band(c_factor_unc.ReadAsArray() * 10000, c_factor_unc_fname, c_factor.GetProjectionRef(), c_factor.GetGeoTransform())

        good_pixels_g = array_to_raster((good_pixels).astype(int), B2)
        c_factor_mask = reproject_data(good_pixels_g, B2, xRes = xRes, yRes = yRes, resample= gdal.GRIORA_NearestNeighbour).data
        mask_fname = band_file.replace('_sur', '_cfactor_mask_%s'%(fname_postfix))
        # print(mask_fname)
        save_band((c_factor_mask).astype(int), mask_fname, projectionRef, geom)


        c_factor = reproject_data(c_factor, B2, xRes = xRes, yRes = yRes).data
        # mask = ~np.isfinite(c_factor)
        # c_factor[~np.isfinite(c_factor)] = 1
        
        data = gdal.Open(band_file).ReadAsArray() / 10000
        mask = data<0
        ### NBAR
        img = (gdal.Open(band_file).ReadAsArray() * c_factor).astype(int)
        nbar_fname = band_file.replace('_sur', '_nbar_%s'%(fname_postfix))
        # print(nbar_fname)
        img[mask] = -9999
        save_band(img, nbar_fname, projectionRef, geom)
        
        # c_factor_unc = reproject_data(c_factor_unc, B2, xRes = xRes, yRes = yRes).data
        
        band_unc_file = B2.replace('B02_sur.tif', band + '_sur_unc.tif')
        # band_unc = gdal.Open(band_unc_file).ReadAsArray() / 10000.

        ### NBAR_UNC
        img = np.sqrt(c_factor**2 * (gdal.Open(band_unc_file).ReadAsArray() / 10000.)**2
                      + reproject_data(c_factor_unc, B2, xRes = xRes, yRes = yRes).data**2 * (gdal.Open(band_file).ReadAsArray() / 10000)**2)
        
        rho_mod_s2_ang = array_to_raster(rho_mod_s2_ang, B2)
        # rho_mod_s2_ang = reproject_data(rho_mod_s2_ang, B2, xRes = xRes, yRes = yRes).data
        # distance = (reproject_data(rho_mod_s2_ang, B2, xRes = xRes, yRes = yRes).data - (gdal.Open(band_file).ReadAsArray() / 10000))**2
       
        img = np.sqrt((reproject_data(rho_mod_s2_ang, B2, xRes = xRes, yRes = yRes).data - (gdal.Open(band_file).ReadAsArray() / 10000))**2 + img**2)
        img[~np.isfinite(img)] = 1
        img = np.minimum(img, 1)
        img = (img * 10000).astype(int)
        nbar_unc_fname = band_file.replace('_sur', '_nbar_unc_%s'%(fname_postfix))
        # print(nbar_unc_fname)
        img[mask] = -9999
        print('Saving NBAR Uncertainty')
        save_band(img, nbar_unc_fname, projectionRef, geom)

    # savig nbar sza
    sza_avg = array_to_raster(sza_avg, B2)
    sza_fname = B2.replace('B02_sur.tif', 'nbar_sza_%s.tif'%(fname_postfix))
    # print(sza_fname)
    geom = sza_avg.GetGeoTransform()
    projectionRef = sza_avg.GetProjectionRef()
    sza_avg = (sza_avg.ReadAsArray() * 100).astype(int)
    save_band(sza_avg, sza_fname, projectionRef, geom)

# s2_file_dir = '/home/users/marcyin/nceo_ard/public/S2/20/M/LE/S2A_MSIL1C_20190728T143751_N0208_R096_T20MLE_20190728T180109.SAFE'
# logger=None
# nbar_sza = 'average'
# mosaic_start_date = datetime.datetime(2021, 1, 1)
# mosaic_end_date = datetime.datetime(2021, 4, 1)

if __name__ == '__main__':
    s2_dir = "/mnt/d/Test_Data/JRC_testing_data/luxcarta_testing_data/SIAC/S2A_MSIL1C_20220825T155541_N0400_R011_T17QRG_20220825T205824.SAFE"
    s2_dir = "S2A_MSIL1C_20220808T142721_N0400_R053_T20MRD_20220808T192835.SAFE"
    viirs_dir = ""
    create_nbar(
        s2_dir,
        nbar_sza='temporal_average_sza',
        mosaic_start_date=datetime.datetime.strptime('2022-08-25', '%Y-%m-%d'),
        mosaic_end_date=datetime.datetime.strptime('2022-08-25', '%Y-%m-%d'),
        Gee = False,
        use_VIIRS = True
        )