import os
import time
import gdal
import getpass
import logging
import requests
import tempfile 
import numpy as np
from osgeo import osr
from glob import glob
from shutil import copy2
from six.moves import input
from functools import partial
from os.path import expanduser
from multiprocessing import Pool
from datetime import datetime, timedelta
from SIAC.create_logger import create_logger
from SIAC.modis_tile_cal import get_vector_hv, get_raster_hv
#from create_logger import create_logger
#from modis_tile_cal import get_vector_hv, get_raster_hv

home = expanduser("~")
file_path = os.path.dirname(os.path.realpath(__file__))
test_url = 'https://e4ftl01.cr.usgs.gov/MOTA/MCD43A1.006/2000.02.24/MCD43A1.A2000055.h32v08.006.2016101152216.hdf.xml'

def get_auth(logger):
    try:
        auth = tuple([os.environ['Earthdata_user'], os.environ['Earthdata_pass']])
        with requests.Session() as s:                                            
            s.auth = auth                                                        
            r1     = s.get(test_url)                                             
            r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
        if r.status_code == 401:
            logger.error('Wrong username and password are set for Earthdata_user and Earthdata_pass Environment variables.') 
        else:
            with open(file_path + '/data/.earthdata_auth', 'wb') as f:
                for i in auth:               
                    f.write((i + '\n').encode())
    except:
        #logger.error('Environment variables Earthdata_user and Earthdata_pass for accessing NASA Earthdata are not set, and trying to read from file or input.')
        pass
    if os.path.exists(file_path + '/data/.earthdata_auth'):
        try:
            username, password = np.loadtxt(file_path + '/data/.earthdata_auth', dtype=str)
            auth = tuple([username, password])
            with requests.Session() as s:
                s.auth = auth
                r1     = s.get(test_url)
                r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
            while r.status_code == 401:                 
                logger.error('Wrong username and password stored, please enter again')
                username = input('Username for NASA Earthdata: ')      
                password = getpass.getpass('Password for NASA Earthdata: ')
                auth = tuple([username, password])
                with requests.Session() as s:
                    s.auth = auth
                    r1     = s.get(test_url)
                    r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
                
            os.remove(file_path + '/data/.earthdata_auth')    
            with open(file_path + '/data/.earthdata_auth', 'wb') as f:
                for i in auth:                                 
                    f.write((i + '\n').encode())              
            auth = tuple([username, password])

        except:
            logger.error('Please provide NASA Earthdata username and password for downloading MCD43 data, which can be applied here: https://urs.earthdata.nasa.gov.')
            username = input('Username for NASA Earthdata: ')
            password = getpass.getpass('Password for NASA Earthdata: ')
            auth = username, password                              
            with requests.Session() as s:
                s.auth = auth
                r1     = s.get(test_url)
                r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
            while r.status_code == 401:                 
                logger.error('Wrong username and password typed, please enter again')
                username = input('Username for NASA Earthdata: ')      
                password = getpass.getpass('Password for NASA Earthdata: ')
                auth = tuple([username, password])
                with requests.Session() as s:
                    s.auth = auth
                    r1     = s.get(test_url)
                    r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
            os.remove(file_path + '/data/.earthdata_auth')
            with open(file_path + '/data/.earthdata_auth', 'wb') as f:     
                for i in auth:                                     
                    f.write((i + '\n').encode())
            auth = tuple([username, password])
    else:
        username = input('Username for NASA Earthdata: ')
        password = getpass.getpass('Password for NASA Earthdata: ')
        auth = username, password
        with requests.Session() as s:
            s.auth = auth
            r1     = s.get(test_url)
            r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
        while r.status_code == 401:                 
            logger.error('Wrong username and password typed, please enter again')
            username = input('Username for NASA Earthdata: ')      
            password = getpass.getpass('Password for NASA Earthdata: ')
            auth = tuple([username, password])
            with requests.Session() as s:
                s.auth = auth
                r1     = s.get(test_url)
                r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'}) 
        with open(file_path + '/data/.earthdata_auth', 'wb') as f:
            for i in auth: 
                f.write((i + '\n').encode())
        auth = tuple([username, password])
    return auth


def find_files(aoi, obs_time, mcd43_dir, temporal_window = 16,jasmin = False):
    days   = [(obs_time - timedelta(days = int(i))).strftime('%Y.%m.%d') for i in np.arange(temporal_window, 0, -1)] + \
             [(obs_time + timedelta(days = int(i))).strftime('%Y.%m.%d') for i in np.arange(0, temporal_window+1,  1)]
    try:
        tiles = get_vector_hv(aoi)
    except:
        try:
            tiles = get_raster_hv(aoi)
        except:
            raise IOError('AOI has to be raster or vector object/files.')
    fls = zip(np.repeat(tiles, len(days)), np.tile(days, len(tiles)))
    ret = []
    need_grab = []
    
    for (tile, the_date) in fls:
        the_jday = datetime.strptime(the_date, '%Y.%m.%d').strftime("%Y%j")
        potential_fname = "MCD43A1.A{:s}.{:s}.*.hdf".format(the_jday, tile)
        the_files = glob(mcd43_dir + "/" + potential_fname)
        if jasmin:
            if len(the_files) == 0 :
                jasmin_mcd43_dir = '/neodc/modis/data/MCD43A1/collection6/'
                jasmin_date = datetime.strptime(the_date, '%Y.%m.%d').strftime('/%Y/%m/%d/')
                potential_fname = "MCD43A1.A{:s}.{:s}.*.hdf".format(the_jday, tile)
                jasminfs = glob(jasmin_mcd43_dir + "/" +  jasmin_date + potential_fname)
                for ii in jasminfs:
                    copy2(ii, mcd43_dir)
                the_files += glob(mcd43_dir + "/" + potential_fname)
        if len(the_files) > 0:
            ret.append(the_files[0])
        else:
            need_grab.append((tile, the_date))
    if len(need_grab) > 1:
        #p = Pool(min(len(need_grab), 8))
        ret_get = list(map(get_one_tile, need_grab))
        #p.close()
        #p.join()
        ret.extend(ret_get)
    elif len(need_grab) == 1:
        ret_get = get_one_tile(need_grab[0])
        ret.extend([ret_get])
        
    return ret

def get_one_tile(tile_date):
    base = 'https://e4ftl01.cr.usgs.gov/MOTA/MCD43A1.006/'
    tile, date = tile_date
    for j in range(100):
        r = requests.get(base + date)
        fname = [i.split('>')[-1] for i in r.content.decode().split('</') if (i[-3:]=='hdf') & (tile in i)]
        if len(fname) == 1:
            break
        else:
            #print(base + date)
            time.sleep(1)
    return base + date + '/' + fname[0]

def downloader(url_fname, auth):
    url, fname = url_fname
    fname = os.path.abspath(fname)
    with requests.Session() as s:
        s.max_redirects = 100000
        s.auth = auth
        r1     = s.get(url)
        r      = s.get(r1.url, stream=True, headers={'user-agent': 'My app'})
        if r.ok:
            remote_size = int(r.headers['Content-Length'])
            if os.path.exists(fname):
                local_size = os.path.getsize(fname)
                if local_size != remote_size:
                    os.remove(fname)
                    data = r.content
                    if len(data) == remote_size:
                        with open(fname, 'wb') as f:
                            f.write(data)
                    else:
                        raise IOError('Failed to download the whole file.')
            else:
                data = r.content
                if len(data) == remote_size:
                    with open(fname, 'wb') as f:
                        f.write(data)
                else:           
                    raise IOError('Failed to download the whole file.')
        else:
            print(r.content)
        
def daily_vrt_jasmin(fnames_date, vrt_dir = None):
    temp1 = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Band_Mandatory_Quality_%s'
    temp2 = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Parameters_%s'

    spatialRef = osr.SpatialReference()
    spatialRef.ImportFromProj4('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs')

    fnames, date = fnames_date
    fnames = list(map(os.path.abspath, fnames))
    all_files = fnames
    ''' this maybe problematic when reading and writing the same time
    all_files = []
    for fname in fnames:
        all_files += glob(os.path.dirname(fname) + '/MCD43A1.A%s.h??v??.006.*.hdf'%date)
    '''
    if vrt_dir is None:
        vrt_dir = './MCD43_VRT/' 
    if not os.path.exists(vrt_dir):
        os.mkdir(vrt_dir)
    d = datetime.strptime(date, '%Y%j').strftime('%Y-%m-%d')
    date_dir = vrt_dir + '/' + '%s/'%d
    if not os.path.exists(date_dir):                               
        os.mkdir(date_dir) 
    for temp in [temp1, temp2]:                                                      
        for band in ['Band1','Band2','Band3','Band4','Band5','Band6','Band7', 'vis', 'nir', 'shortwave']:
            bs = []                                                                  
            for fname in all_files:                                                     
                bs.append(temp%(fname, band))                                        
            gdal.BuildVRT(date_dir + '_'.join(['MCD43', date, bs[0].split(':')[-1]])+'.vrt', bs, outputSRS = spatialRef).FlushCache()

def daily_vrt(fnames_date, vrt_dir = None):
    temp1 = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Band_Mandatory_Quality_%s'
    temp2 = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Parameters_%s'

    spatialRef = osr.SpatialReference()
    spatialRef.ImportFromProj4('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs')

    fnames, date = fnames_date
    fnames = list(map(os.path.abspath, fnames))
    all_files = []
    for fname in fnames:
        all_files += glob(os.path.dirname(fname) + '/MCD43A1.A%s.h??v??.006.*.hdf'%date)
    if vrt_dir is None:
        vrt_dir = './MCD43_VRT/'
    if not os.path.exists(vrt_dir):
        os.mkdir(vrt_dir)
    d = datetime.strptime(date, '%Y%j').strftime('%Y-%m-%d')
    date_dir = vrt_dir + '/' + '%s/'%d
    if not os.path.exists(date_dir):
        os.mkdir(date_dir)
    for temp in [temp1, temp2]:
        for band in ['Band1','Band2','Band3','Band4','Band5','Band6','Band7', 'vis', 'nir', 'shortwave']:
            bs = []
            for fname in all_files:
                bs.append(temp%(fname, band))
            gdal.BuildVRT(date_dir + '_'.join(['MCD43', date, bs[0].split(':')[-1]])+'.vrt', bs, outputSRS = spatialRef).FlushCache()

def get_mcd43(aoi, obs_time, mcd43_dir = './MCD43/', vrt_dir = './MCD43_VRT/', logger = None, jasmin = False):
    mcd43_dir = os.path.expanduser(mcd43_dir) # incase ~ used
    vrt_dir   = os.path.expanduser(vrt_dir)
    if logger is None:
        logger = create_logger()
    logger.info('Querying MCD43 files...')
    ret = find_files(aoi, obs_time, mcd43_dir, temporal_window = 16, jasmin=jasmin)
    urls = [granule for granule in ret if granule.find("http") >= 0]
    url_fnames= [granule for granule in ret if granule.find("http") < 0]
    flist = url_fnames
    logger.info("Will need to download {:d} files, ".format(len(urls)) +
                "{:d} are already present".format(len(url_fnames)))
    if len(urls) > 0:
        auth = get_auth(logger)
        par = partial(downloader, auth = auth)
        logger.info('Start downloading MCD43 for the AOI, this may take some time.')
        url_fnames_to_get = [[i, mcd43_dir + '/' + i.split('/')[-1]] for i in urls]
        if jasmin:
            ret = list(map(par, url_fnames_to_get))
        else:
            p = Pool(5)
            logger.info('Start downloading...')
            ret = p.map(par, url_fnames_to_get)
            p.close()
            p.join()
        flist.extend([x[1] for x in url_fnames_to_get])
    flist = np.array(flist)
    all_dates = np.array([i.split('/')[-1] .split('.')[1][1:9] for i in flist])          
    udates = np.unique(all_dates)  
    fnames_dates =  [[flist[all_dates==date].tolist(),date] for date in udates]
    
    if jasmin:
#         vrt_dir = tempfile.TemporaryDirectory(dir  =  vrt_dir).name + '/'
#         while os.path.exists(vrt_dir):
#             vrt_dir = tempfile.TemporaryDirectory(dir  =  vrt_dir).name + '/'
#         os.mkdir(vrt_dir)
        vrt_dir = tempfile.mkdtemp(suffix="", prefix="tmp", dir  =  vrt_dir) + '/'
        logger.info('Creating daily VRT...')
        par = partial(daily_vrt_jasmin, vrt_dir = vrt_dir)
        njobs = min(len(fnames_dates), 4)
        #p = Pool(njobs) 
        #p.map(par, fnames_dates)
        #p.close()                           
        #p.join()
        list(map(par, fnames_dates))
        logger.info('Finished creating vrt...')
    else:
#         vrt_dir = tempfile.TemporaryDirectory(dir  =  vrt_dir).name + '/'
#         while os.path.exists(vrt_dir):
#             vrt_dir = tempfile.TemporaryDirectory(dir  =  vrt_dir).name + '/'
#         os.mkdir(vrt_dir)
        logger.info('Creating daily VRT...')
        vrt_dir = tempfile.mkdtemp(suffix="", prefix="tmp", dir  =  vrt_dir) + '/'
        par = partial(daily_vrt, vrt_dir = vrt_dir)
        njobs = min(len(fnames_dates), 4)
        p = Pool(njobs)
        p.map(par, fnames_dates)
        p.close()                           
        p.join()
        logger.info('Finished creating vrt...')
    #par(fnames_dates[0])
    handlers = logger.handlers[:]
    for handler in handlers:
        handler.close()
        logger.removeHandler(handler)
    return vrt_dir

if __name__ == '__main__':
    aoi = '~/DATA/S2_MODIS/l_data/LC08_L1TP_014034_20170831_20170915_01_T1/AOI.json'
    aoi = '~/DATA/S2_MODIS/l_data/LC08_L1TP_014034_20170831_20170915_01_T1/aot.tif'
    obs_time = datetime(2017, 7, 8, 10, 8, 20)
    ret = get_mcd43(aoi, obs_time, mcd43_dir = '~/hep/MCD43/', vrt_dir = '~/DATA/Multiply/MCD43/')




