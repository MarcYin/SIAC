import os
import sys
import argparse
import requests
import warnings
warnings.filterwarnings("ignore") 
import numpy as np
from glob import glob
from SIAC.get_MCD43 import get_mcd43
from datetime import datetime
from SIAC.the_aerosol import solve_aerosol
from SIAC.create_logger import create_logger
from SIAC.the_correction import atmospheric_correction
from SIAC.s2_preprocessing import s2_pre_processing
from SIAC.downloaders import downloader
from SIAC.multi_process import parmap
from os.path import expanduser
home = expanduser("~")
file_path = os.path.dirname(os.path.realpath(__file__))

def SIAC_S2(s2_t, send_back = False, mcd43 = home + '/MCD43/', vrt_dir = home + '/MCD43_VRT/', aoi = None):
    if not os.path.exists(file_path + '/emus/'):
        os.mkdir(file_path + '/emus/')
    if len(glob(file_path + '/emus/' + 'isotropic_MSI_emulators_*_x?p_S2?.pkl')) < 12:
        url = 'http://www2.geog.ucl.ac.uk/~ucfafyi/emus/'
        req = requests.get(url)
        to_down = []
        for line in req.text.split():
            if 'MSI' in line:
                fname   = line.split('"')[1].split('<')[0]
                if 'MSI' in fname:
                    to_down.append([fname, url])
        f = lambda fname_url: downloader(fname_url[0], fname_url[1], file_path + '/emus/')
        parmap(f, to_down)
    rets = s2_pre_processing(s2_t)
    aero_atmos = []
    for ret in rets:
        ret += (mcd43, vrt_dir, aoi)
        aero_atmo = do_correction(*ret)
        if send_back:
            aero_atmos.append(aero_atmo)
    if send_back:
        return aero_atmos

def do_correction(sun_ang_name, view_ang_names, toa_refs, cloud_name, \
                  cloud_mask, metafile, mcd43 = home + '/MCD43/', vrt_dir = home + '/MCD43_VRT/', aoi=None):

    if os.path.realpath(mcd43) in os.path.realpath(home + '/MCD43/'):
        if not os.path.exists(home + '/MCD43/'):
            os.mkdir(home + '/MCD43/')

    if os.path.realpath(vrt_dir) in os.path.realpath(home + '/MCD43_VRT/'):
        if not os.path.exists(home + '/MCD43_VRT/'):
            os.mkdir(home + '/MCD43_VRT/')

    base = os.path.dirname(toa_refs[0])
    base = toa_refs[0].replace('B01.jp2', '')
    with open(metafile) as f:
        for i in f.readlines():
            if 'SENSING_TIME' in i:
                sensing_time = i.split('</')[0].split('>')[-1]
                obs_time = datetime.strptime(sensing_time, u'%Y-%m-%dT%H:%M:%S.%fZ')
            if 'TILE_ID' in i:
                sat  = i.split('</')[0].split('>')[-1].split('_')[0]
                tile = i.split('</')[0].split('>')[-1]
    log_file = os.path.dirname(metafile) + '/SIAC_S2.log'
    logger = create_logger(log_file)
    logger.info('Starting atmospheric corretion for %s'%tile)
    if not np.all(cloud_mask):
        handlers = logger.handlers[:]
        for handler in handlers:
            handler.close()
            logger.removeHandler(handler)
        get_mcd43(toa_refs[0], obs_time, mcd43_dir = mcd43, vrt_dir = vrt_dir, log_file = log_file)
        #logger = create_logger(log_file)
    else:
        logger.info('No clean pixel in this scene and no MCD43 is downloaded.')
    sensor_sat = 'MSI', sat
    band_index  = [1,2,3,7,11,12]
    band_wv    = [469, 555, 645, 859, 1640, 2130]
    toa_bands   = (np.array(toa_refs)[band_index,]).tolist()
    view_angles = (np.array(view_ang_names)[band_index,]).tolist()
    sun_angles  = sun_ang_name
    #logger.info('Running SIAC for tile: %s on %s'%(tile, obs_time.strftime('%Y-%M-%d')))
    aero = solve_aerosol(sensor_sat,toa_bands,band_wv, band_index,view_angles,\
                         sun_angles,obs_time,cloud_mask, gamma=10., spec_m_dir= \
                         file_path+'/spectral_mapping/', emus_dir=file_path+'/emus/', \
                         mcd43_dir=vrt_dir, aoi=aoi, log_file = log_file)
    aero._solving()
    toa_bands  = toa_refs
    view_angles = view_ang_names
    aot = base + 'aot.tif'
    tcwv = base + 'tcwv.tif'
    tco3 = base + 'tco3.tif'
    aot_unc = base + 'aot_unc.tif'
    tcwv_unc = base + 'tcwv_unc.tif'
    tco3_unc = base + 'tco3_unc.tif'
    rgb = [toa_bands[3], toa_bands[2], toa_bands[1]]
    band_index = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    atmo = atmospheric_correction(sensor_sat,toa_bands, band_index,view_angles,\
                                  sun_angles, aot = aot, cloud_mask = cloud_mask,\
                                  tcwv = tcwv, tco3 = tco3, aot_unc = aot_unc, \
                                  tcwv_unc = tcwv_unc, tco3_unc = tco3_unc, rgb = \
                                  rgb, emus_dir=file_path+'/emus/', log_file = log_file)
    atmo._doing_correction()
    return aero, atmo

def exe():
    parser = argparse.ArgumentParser(description='Sentinel 2 Atmospheric Correction Excutable')
    parser.add_argument('-f', "--file_path",      help='Sentinel 2 file path', required=True)
    parser.add_argument("-m", "--MCD43_file_dir", help="Directory where you store MCD43A1.006 data", default = home + '/MCD43/')
    parser.add_argument("-v", "--vrt_dir",        help="Where MCD43 vrt stored.",                    default = home + '/MCD43_VRT/')
    parser.add_argument("-a", "--aoi",            help="Area of Interest.",                          default = None)
    args = parser.parse_args()
    SIAC_S2(s2_t=args.file_path, mcd43=args.MCD43_file_dir, vrt_dir=args.vrt_dir, aoi=args.aoi)

if __name__ == '__main__':
    exe()

