import os
import kernels
import datetime
import numpy as np
import pynndescent
from glob import glob
from osgeo import gdal, osr
from pyproj import Proj, transform, CRS
from SIAC.read_MCD43 import robust_smooth
from reproject import array_to_raster, reproject_data

import os
import logging


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
    Solar geometry unsing subsolar point and atan2: https://www.sciencedirect.com/science/article/pii/S0960148121004031?via%3Dihub#bib1
    
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
        
    return f0, f1, f2

def interp_fs(obs_cwv, modis_cwv, f0, f1, f2):
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
    for i in range(len(obs_cwv)):
        lr_wv = lr_wvs[i]
        wv = obs_cwv[i]
        l_modis_ind = modis_cwv.tolist().index(lr_wv[0])
        r_modis_ind = modis_cwv.tolist().index(lr_wv[1])
        
        if l_modis_ind == r_modis_ind :
            f0b = f0[l_modis_ind]
            f1b = f1[l_modis_ind]
            f2b = f2[l_modis_ind]
        else:
            slope = (f0[r_modis_ind] - f0[l_modis_ind]) / (lr_wv[1] - lr_wv[0])
            f0b = (wv - lr_wv[0]) * slope + f0[l_modis_ind]

            slope = (f1[r_modis_ind] - f1[l_modis_ind]) / (lr_wv[1] - lr_wv[0])
            f1b = (wv - lr_wv[0]) * slope + f1[l_modis_ind]
            
            slope = (f2[r_modis_ind] - f2[l_modis_ind]) / (lr_wv[1] - lr_wv[0])
            f2b = (wv - lr_wv[0]) * slope + f2[l_modis_ind]

        
        interped_f0.append(f0b)
        interped_f1.append(f1b)
        interped_f2.append(f2b)
    
    interped_f0 = np.array(interped_f0)
    interped_f1 = np.array(interped_f1)
    interped_f2 = np.array(interped_f2)
    
    return interped_f0, interped_f1, interped_f2

def save_band(array, outputFileName, projectionRef, geotransform):
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
        das, ws, mg = read_MCD3_GEE(mcd43_gee_folder, boa_bands, outputBounds,  500, 500, dstSRS)
        
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
        das, ws, mg = read_VNP43MA1(VNP43_fnames_dates, bands, outputBounds, 500, 500, dstSRS)
        

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


def get_time_series_kernel_weights(s2_file_dir, logger=None, nbar_sza = 'average'):
    '''
    s2_file_dir: directory of Sentinel-2 file
    logger: logger object
    '''
    if logger is None:
        logger = create_logger()
        print = logger.info
    else:
        print = logger.info

    B2 = glob(s2_file_dir + '/GRANULE/*/IMG_DATA/*B02_sur.tif')[0]
    toa_dir = os.path.dirname(B2)
    ang_dir = os.path.dirname(toa_dir) + '/ANG_DATA/'

    kernel_weights_fname = toa_dir + '/kernel_weights_time_series.npz'
    if not os.path.exists(kernel_weights_fname):
        print('No kernel weights file found.')
        # return
    print('Reading kernel weights.')
    f = np.load(kernel_weights_fname)
    das = f.f.das * 0.001
    mask = f.f.ws > 0.5
    bands_combined_mask = np.all(mask, axis=1)
    das = das.transpose(1, 2, 0, 3, 4)[:, :, bands_combined_mask]
    
    s2_bands = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']
    # band_inds = [1, 2, 3, 4, 5, 6, 7, 8, 11, 12]
    
    saa_sza = ang_dir + 'SAA_SZA.tif'
    g = gdal.Warp('', saa_sza, format = 'MEM', xRes = 1000, yRes = 1000, resampleAlg = gdal.GRA_Bilinear)
    saa, sza = g.ReadAsArray() / 100.

    if nbar_sza == 'average':
        # print(sza, vza, saa, vaa)
        g = gdal.Open(B2)
        geom = g.GetGeoTransform()
        x = np.repeat(np.arange(saa.shape[0]), saa.shape[1])
        y = np.tile(np.arange(saa.shape[0]), saa.shape[1])

        geom = g.GetGeoTransform()
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

        sza_avg = np.nanmean(sza_avg)
        saa_avg = np.nanmean(saa_avg)
        # print(sza_avg.mean(), sza.mean())
        fname_postfix = 'sza_avg'
    
    elif nbar_sza == 'use_s2':
        sza_avg = np.nanmean(sza)
        saa_avg = np.nanmean(saa)
        fname_postfix = 'sza_s2'
    elif (nbar_sza > 0) & (nbar_sza < 60):
        sza_avg = nbar_sza
        saa_avg = np.nanmean(saa)
        fname_postfix = 'sza_%02d'%int(nbar_sza)
    else:
        print('nbar_sza should be "average" or float value between 0 and 60')
        # return
    print('Compute kernel weights ensemble.')
    c_factors = []
    rho_mod_s2_angs = []
    for ii, ind in enumerate([1, 2, 3, 7, 11, 12]):
        band = s2_bands[ind]

        vaa_vza = ang_dir + 'VAA_VZA_%s.tif'%band

        vaa, vza = gdal.Open(vaa_vza).ReadAsArray() / 100
        vaa_avg = np.nanmean(vaa)
        vza_avg = np.nanmean(vza)

        raa_avg = (vaa_avg - np.nanmean(saa)) % 360
        
        angs = np.atleast_2d([vza_avg, np.nanmean(sza), raa_avg]).T

        kk = get_kk(angs)
        k_vol  = kk.Ross          
        k_geo  = kk.Li  
        kers = np.array([np.ones(k_vol.shape), k_vol, k_geo])
        rho_mod_s2_ang = kers * das[ii]
        rho_mod_s2_ang = rho_mod_s2_ang.sum(axis=0)
        
        angs = np.atleast_2d([0., sza_avg, 0.]).T
        kk = get_kk(angs)
        k_vol  = kk.Ross          
        k_geo  = kk.Li  
        kers = np.array([np.ones(k_vol.shape), k_vol, k_geo])
        rho_mod_s2_nadir = kers * das[ii]
        rho_mod_s2_nadir = rho_mod_s2_nadir.sum(axis=0)

        c_factor = rho_mod_s2_nadir / rho_mod_s2_ang
        c_factors.append(c_factor)
        rho_mod_s2_angs.append(rho_mod_s2_ang)
        # print('Mean c factor, vza, sza, raa, nbar_sza: %.03f, %.03f, %.03f, %.03f, %.03f'%(np.nanmean(c_factor), vza.mean(), sza.mean(), raa.mean(), sza_avg.mean()))

    c_factors = np.array(c_factors)
    rho_mod_s2_angs = np.array(rho_mod_s2_angs)
    good_ones = np.all(np.isfinite(c_factors), axis=0)
    c_factors = c_factors[:, good_ones]
    rho_mod_s2_angs = rho_mod_s2_angs[:, good_ones] 

    closest_s2_to_modis = [B2.replace('B02_sur.tif', s2_bands[ind] + '_sur.tif') for ind in [1, 2, 3, 7, 11, 12]]

    print('Create index for nearest neighbour search.')
    index = pynndescent.NNDescent(rho_mod_s2_angs.T)    

    g = gdal.BuildVRT('', closest_s2_to_modis, xRes = 60, yRes = 60, resampleAlg = gdal.GRA_Average, separate = True)
    s2_data = g.ReadAsArray()
    s2_data_mask = np.all(s2_data > 0, axis=0)
    s2_data = s2_data[:, s2_data_mask] / 10000.

    s2_data = np.array_split(s2_data.T, 100)
    c_factor_slice_weighted_means = []
    c_factor_slice_weighted_stds = []
    c_factor_slice_mean_dists = []

    ref_range = np.max(rho_mod_s2_angs, axis=1) - np.min(rho_mod_s2_angs, axis=1)
    
    c_factor_range = np.max(c_factors, axis=1) - np.min(c_factors, axis=1)
    
    diff_cfactor_to_ref = c_factor_range / ref_range

    print('Compute weighted mean and std.')
    for ii, s2_data_slice in enumerate(s2_data):
        if ii % 10 == 0:
            # print out progress in percent
            print('%.02f%%'%(ii/len(s2_data)*100))
        inds, dists = index.query(s2_data_slice, k=50)
        
        s2_data[ii] = None
        c_factors_slice = c_factors[:, inds]

        weights = 1./dists**2
        weights = weights / weights.sum(axis=1)[:, None]
        c_factor_slice_weighted_mean = np.sum(c_factors_slice * weights[None], axis=2)
        c_factor_slice_weighted_std = np.sqrt(np.sum(weights[None] * (c_factors_slice - c_factor_slice_weighted_mean[:, :, None])**2, axis=2))

        # ensemble_y_include_observation = np.concatenate([c_factors_slice, c_factor_slice_weighted_mean[:, :, None]], axis=2)
        # ensemble_y_include_pred_mean_std = np.std(ensemble_y_include_observation, axis=2)
        
        mean_dists = np.mean(dists, axis=1)
        # min_dist = np.min(dists, axis=1)
        # min_dist = mean_dists

        c_factor_slice_mean_dists.append(mean_dists)
        
        # c_factor_slice_weighted_std = c_factor_slice_weighted_std + min_dist[None] / 6 * diff_cfactor_to_ref[:, None] * ensemble_x_include_pred_mean_std
        # https://link.springer.com/chapter/10.1007/978-3-319-91479-4_40
        c_factor_slice_weighted_std = c_factor_slice_weighted_std  \
                                    + mean_dists[None] / 6 * diff_cfactor_to_ref[:, None] #\
                                    # + min_dist[None] / 6 * ref_range[:, None] * ensemble_y_include_pred_mean_std

        c_factor_slice_weighted_means.append(c_factor_slice_weighted_mean)
        c_factor_slice_weighted_stds.append(c_factor_slice_weighted_std)

    c_factor_slice_weighted_means = np.concatenate(c_factor_slice_weighted_means, axis=1)
    c_factor_bands = np.zeros((6, s2_data_mask.shape[0], s2_data_mask.shape[1]))
    c_factor_bands[:, s2_data_mask] = c_factor_slice_weighted_means
    del c_factor_slice_weighted_means

    c_factor_slice_weighted_stds = np.concatenate(c_factor_slice_weighted_stds, axis=1)
    c_factor_bands_std = np.zeros((6, s2_data_mask.shape[0], s2_data_mask.shape[1]))
    c_factor_bands_std[:, s2_data_mask] = c_factor_slice_weighted_stds
    del c_factor_slice_weighted_stds
    
    c_factor_slice_mean_dists = np.concatenate(c_factor_slice_mean_dists)
    c_factor_band_errs = np.zeros((s2_data_mask.shape[0], s2_data_mask.shape[1]))
    c_factor_band_errs[s2_data_mask] = c_factor_slice_mean_dists
    del c_factor_slice_mean_dists

    return c_factor_bands, c_factor_bands_std, c_factor_band_errs

def create_nbar_from_time_series(s2_file_dir, nbar_sza='average', logger=None):

    if logger is None:
        logger = create_logger()
        print = logger.info
    else:
        print = logger.info

    if nbar_sza == 'average':
        fname_postfix = 'sza_avg'    
    elif nbar_sza == 'use_s2':
        fname_postfix = 'sza_s2'
    elif (nbar_sza > 0) & (nbar_sza < 60):
        fname_postfix = 'sza_%02d'%int(nbar_sza)
    else:
        print('nbar_sza should be "average" or float value between 0 and 60')
    
    B2 = glob(s2_file_dir + '/GRANULE/*/IMG_DATA/*B02_sur.tif')[0]
    # toa_dir = os.path.dirname(B2)
    
    c_factor_bands, c_factor_bands_std, c_factor_band_errs = get_time_series_kernel_weights(s2_file_dir, logger=logger, nbar_sza = nbar_sza)
    s2_bands  = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']
    band_inds = [1, 2, 3, 4, 5, 6, 7, 8, 11, 12]
    s2_cwvs   = [442.7, 492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 945.1, 1373.5, 1613.7, 2202.4]
    modis_cwv = [469, 555, 645, 858, 1640, 2130]
    modis_cwv = np.array(modis_cwv)
    

    for ind in band_inds:
        band = s2_bands[ind]
        print('------------------------')
        print(band)
        wv = s2_cwvs[ind]

        if wv < 469:
            l_modis_ind = 0
            r_modis_ind = 0
        elif (wv > 858) & (wv <= 1000):
            l_modis_ind = 3
            r_modis_ind = 3
        elif (wv > 645) & (wv<700):
            l_modis_ind = 2
            r_modis_ind = 2
        elif (wv > 645) & (wv <= 858):
            l_modis_ind = 2
            r_modis_ind = 3
        elif (wv > 1580) & (wv <= 1740):
            l_modis_ind = 4
            r_modis_ind = 4
        elif wv > 2130:
            l_modis_ind = 5
            r_modis_ind = 5
        else:
            wv_diff = modis_cwv - wv
            ind = np.argsort(abs(wv_diff)).astype(int)
            l_modis_ind, r_modis_ind = ind[:2]

        if l_modis_ind == r_modis_ind :
            c_factor = c_factor_bands[l_modis_ind]
            c_factor_std = c_factor_bands_std[l_modis_ind]
        else:
            # slope = (c_factor_bands[r_modis_ind] - c_factor_bands[l_modis_ind]) / (modis_cwv[r_modis_ind] - modis_cwv[l_modis_ind] )
            # c_factor = (wv - modis_cwv[l_modis_ind]) * slope + c_factor_bands[l_modis_ind]
            slope = (wv - modis_cwv[l_modis_ind]) / (modis_cwv[r_modis_ind] - modis_cwv[l_modis_ind] )
            c_factor = (c_factor_bands[r_modis_ind] - c_factor_bands[l_modis_ind]) * slope + c_factor_bands[l_modis_ind]
            c_factor_std = np.sqrt(slope**2 * c_factor_bands_std[l_modis_ind]**2 + (1-slope)**2 * c_factor_bands_std[r_modis_ind]**2)

        band_file = B2.replace('B02_sur.tif', band + '_sur.tif')
        band_g = gdal.Open(band_file)
        
        geom = band_g.GetGeoTransform()
        projectionRef = band_g.GetProjectionRef()

        xRes = abs(geom[1])
        yRes = abs(geom[5])

        c_factor = array_to_raster(c_factor, B2)
        c_factor = reproject_data(c_factor, B2, xRes = xRes, yRes = yRes).data

        c_factor_std = array_to_raster(c_factor_std, B2)
        c_factor_std = reproject_data(c_factor_std, B2, xRes = xRes, yRes = yRes).data

        mask = ~np.isfinite(c_factor)
        c_factor[mask] = 1
        mask_fname = band_file.replace('_sur', '_cfactor_mask_%s'%(fname_postfix))
        # print(mask_fname)
        print('Saving mask')
        save_band(mask.astype(int), mask_fname, projectionRef, geom)
        
        data = band_g.ReadAsArray()
        nbar = (data * c_factor).astype(int)
        nbar_fname = band_file.replace('_sur', '_nbar_%s'%(fname_postfix))
        # print(nbar_fname)
        print('Saving NBAR')
        save_band(nbar, nbar_fname, projectionRef, geom)
        
        band_unc_file = B2.replace('B02_sur.tif', band + '_sur_unc.tif')
        band_unc = gdal.Open(band_unc_file).ReadAsArray() / 10000.

        nbar_unc = np.sqrt(c_factor**2 * band_unc**2 + c_factor_std**2 * (data / 10000.)**2)
        nbar_unc = (nbar_unc * 10000).astype(int)
        nbar_unc_fname = band_file.replace('_sur', '_nbar_unc_%s'%(fname_postfix))
        # print(nbar_unc_fname)
        print('Saving NBAR Uncertainty')
        save_band(nbar_unc, nbar_unc_fname, projectionRef, geom)

    geom = list(band_g.GetGeoTransform())
    geom[1] = 60
    geom[5] = -60
    geom = tuple(geom)
    c_factor_band_errs = (c_factor_band_errs * 10000).astype(int) 
    nbar_ann_mean_dist_fname = B2.replace('_B02_sur', '_nbar_ann_mean_dist_%s'%(fname_postfix))
    # print(nbar_ann_mean_dist_fname)
    save_band(c_factor_band_errs, nbar_ann_mean_dist_fname, projectionRef, geom)
    
def create_nbar(s2_file_dir, nbar_sza='average', logger=None):
    '''
    nbar_sza: 'average' or float value between 0 and 60
    
    '''
    if logger is None:
        logger = create_logger()
        print = logger.info
    else:
        print = logger.info

    B2 = glob(s2_file_dir + '/GRANULE/*/IMG_DATA/*B02_sur.tif')[0]
    toa_dir = os.path.dirname(B2)
    ang_dir = os.path.dirname(toa_dir) + '/ANG_DATA/'

    kernel_weights_fname = toa_dir + '/kernel_weights.npz'
    if not os.path.exists(kernel_weights_fname):
        print('No kernel weights file found.')
        return

    s2_cwv = [442.7, 492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 945.1, 1373.5, 1613.7, 2202.4]
    modis_cwv = [469, 555, 645, 858, 1640, 2130]
    modis_cwv = np.array(modis_cwv)
    f0, f1, f2 = getKernelWeights(toa_dir)
    interped_f0, interped_f1, interped_f2 = interp_fs(s2_cwv, modis_cwv, f0, f1, f2)
    kernel_g = array_to_raster(interped_f0[0], B2)

    saa_sza = ang_dir + 'SAA_SZA.tif'
    saa, sza = reproject_data(saa_sza, kernel_g, resample=1).data / 100


    if nbar_sza == 'average':
        # print(sza, vza, saa, vaa)
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
    elif nbar_sza == 'use_s2':
        sza_avg = sza
        fname_postfix = 'sza_s2'
    elif (nbar_sza > 0) & (nbar_sza < 60):
        sza_avg = np.ones_like(sza) * nbar_sza
        fname_postfix = 'sza_%02d'%int(nbar_sza)
    else:
        print('nbar_sza should be "average" or float value between 0 and 60')
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
        
        kk = get_kk([vza*0, np.ones_like(sza) * sza_avg, raa*0])
        k_vol  = kk.Ross          
        k_geo  = kk.Li  
        kers = np.array([np.ones(k_vol.shape), k_vol, k_geo])
        rho_mod_s2_nadir = kers * np.array([interped_f0[ind], interped_f1[ind], interped_f2[ind]])
        rho_mod_s2_nadir = rho_mod_s2_nadir.sum(axis=0)
        
        c_factor = rho_mod_s2_nadir / rho_mod_s2_ang
        print('Mean c factor, vza, sza, raa, nbar_sza: %.03f, %.03f, %.03f, %.03f, %.03f'%(np.nanmean(c_factor), vza.mean(), sza.mean(), raa.mean(), sza_avg.mean()))


        band_g = gdal.Open(band_file)
        geom = band_g.GetGeoTransform()
        projectionRef = band_g.GetProjectionRef()

        xRes = abs(geom[1])
        yRes = abs(geom[5])

        c_factor = array_to_raster(c_factor, B2)
        
        c_factor_fname = band_file.replace('_sur', '_cfactor_%s'%(fname_postfix))
        print('save band')
        print(c_factor_fname)
        save_band(c_factor.ReadAsArray() * 10000, c_factor_fname, c_factor.GetProjectionRef(), c_factor.GetGeoTransform())
        
        c_factor = reproject_data(c_factor, B2,xRes = xRes, yRes = yRes).data
        
        mask = ~np.isfinite(c_factor)
        c_factor[mask] = 1
        mask_fname = band_file.replace('_sur', '_cfactor_mask_%s'%(fname_postfix))
        print(mask_fname)
        save_band(mask.astype(int), mask_fname, projectionRef, geom)
        
        data = band_g.ReadAsArray()
        nbar = (data * c_factor).astype(int)
        nbar_fname = band_file.replace('_sur', '_nbar_%s'%(fname_postfix))
        print(nbar_fname)
        save_band(nbar, nbar_fname, projectionRef, geom)