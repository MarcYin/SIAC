import os
import time
import gdal
import getpass
import logging
import requests
import numpy as np
from glob import glob
from six.moves import input
from functools import partial
from os.path import expanduser
from multiprocessing import Pool
from datetime import datetime, timedelta
from SIAC.create_logger import create_logger
from SIAC.modis_tile_cal import get_vector_hv, get_raster_hv

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


def find_files(aoi, obs_time, temporal_window = 16):
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
    p = Pool(8)
    ret = p.map(get_one_tile, fls)
    p.close()
    p.join()
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
        
def daily_vrt(fnames_date, vrt_dir = None):
    temp1 = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Band_Mandatory_Quality_%s'
    temp2 = 'HDF4_EOS:EOS_GRID:"%s":MOD_Grid_BRDF:BRDF_Albedo_Parameters_%s'

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
            gdal.BuildVRT(date_dir + '_'.join(['MCD43', date, bs[0].split(':')[-1]])+'.vrt', bs).FlushCache()

def get_mcd43(aoi, obs_time, mcd43_dir = './MCD43/', vrt_dir = './MCD43_VRT/', log_file = None):
    mcd43_dir = os.path.expanduser(mcd43_dir) # incase ~ used
    vrt_dir   = os.path.expanduser(vrt_dir)
    logger = create_logger(log_file)
    logger.propagate = False
    logger.info('Start downloading MCD43 for the AOI, this may take some time.')
    logger.info('Query files...')
    ret = find_files(aoi, obs_time, temporal_window = 16)
    url_fnames = [[i, mcd43_dir + '/' + i.split('/')[-1]] for i in ret]
    p = Pool(5)
    logger.info('Start downloading...')
    auth = get_auth(logger)
    par = partial(downloader, auth = auth)
    ret = p.map(par, url_fnames)
    p.close()
    p.join()
    flist = np.array(url_fnames)[:,1]
    all_dates = np.array([i.split('/')[-1] .split('.')[1][1:9] for i in flist])          
    udates = np.unique(all_dates)  
    fnames_dates =  [[flist[all_dates==date].tolist(),date] for date in udates]
    logger.info('Creating daily VRT...')
    par = partial(daily_vrt, vrt_dir = vrt_dir)
    p = Pool(len(fnames_dates)) 
    p.map(par, fnames_dates)
    p.close()                           
    p.join()
    handlers = logger.handlers[:]
    for handler in handlers:
        handler.close()
        logger.removeHandler(handler)

if __name__ == '__main__':
    aoi = '~/DATA/S2_MODIS/l_data/LC08_L1TP_014034_20170831_20170915_01_T1/AOI.json'
    aoi = '~/DATA/S2_MODIS/l_data/LC08_L1TP_014034_20170831_20170915_01_T1/aot.tif'
    obs_time = datetime(2017, 7, 8, 10, 8, 20)
    ret = get_mcd43(aoi, obs_time, mcd43_dir = '~/hep/MCD43/', vrt_dir = '~/DATA/Multiply/MCD43/')




