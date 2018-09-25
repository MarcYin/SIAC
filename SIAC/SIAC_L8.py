import os
import sys
import gdal
import requests
import numpy as np
from glob import glob
from datetime import datetime
from SIAC.the_aerosol import solve_aerosol
from SIAC.the_correction import atmospheric_correction
from SIAC.l8_preprocessing import l8_pre_processing
from SIAC.get_MCD43 import get_mcd43
from SIAC.downloaders import downloader
from SIAC.multi_process import parmap
from os.path import expanduser
home = expanduser("~")
file_path = os.path.dirname(os.path.realpath(__file__))

def SIAC_L8(l8_dir, send_back = False, mcd43 = home + '/MCD43/', vrt_dir = home + '/MCD43_VRT/', aoi = None):
    file_path = os.path.dirname(os.path.realpath(__file__))
    if not os.path.exists(file_path + '/emus/'):
        os.mkdir(file_path + '/emus/')
    #print(file_path)
    if len(glob(file_path + '/emus/' + 'isotropic_OLI_emulators_*_x?p_L8.pkl')) < 6:
        url = 'http://www2.geog.ucl.ac.uk/~ucfafyi/emus/'
        req = requests.get(url)
        to_down = []
        for line in req.text.split():
            if 'OLI' in line:
                fname   = line.split('"')[1].split('<')[0]
                if ('OLI' in fname) & ('L8' in fname):
                    to_down.append([fname, url])
        f = lambda fname_url: downloader(fname_url[0], fname_url[1], file_path + '/emus/')
        parmap(f, to_down)

    rets = l8_pre_processing(l8_dir)
    aero_atmos = []
    for ret in rets:
        ret += (mcd43, vrt_dir, aoi)
        #sun_ang_name, view_ang_names, toa_refs, cloud_name, cloud_mask, metafile = ret
        aero_atmo = do_correction(*ret)
        if send_back:
            aero_atmos.append(aero_atmo)
    if send_back:
        return aero_atmos

def do_correction(sun_ang_name, view_ang_names, toa_refs, qa_name, cloud_mask, \
                  metafile, mcd43 = home + '/MCD43/', vrt_dir = home + '/MCD43_VRT/', aoi = None):

    if os.path.realpath(mcd43) in os.path.realpath(home + '/MCD43/'):
        if not os.path.exists(home + '/MCD43/'):
            os.mkdir(home + '/MCD43/')
                 
    if os.path.realpath(vrt_dir) in os.path.realpath(home + '/MCD43_VRT/'):
        if not os.path.exists(home + '/MCD43_VRT/'):
            os.mkdir(home + '/MCD43_VRT/')

    base = os.path.dirname(toa_refs[0])
    with open(metafile) as f:
        for line in f:
            if 'REFLECTANCE_MULT_BAND' in line:
                scale = float(line.split()[-1])
            elif 'REFLECTANCE_ADD_BAND' in line:
                off = float(line.split()[-1])
            elif 'DATE_ACQUIRED' in line:
                date = line.split()[-1]
            elif 'SCENE_CENTER_TIME' in line:
                time = line.split()[-1]

    datetime_str= date + time
    obs_time    = datetime.strptime(datetime_str.split('.')[0], '%Y-%m-%d"%H:%M:%S')
    get_mcd43(toa_refs[0], obs_time, mcd43_dir = mcd43, vrt_dir = vrt_dir)
    sensor_sat  = 'OLI', 'L8'
    band_index  = [1,2,3,4,5,6]
    band_wv     = [469, 555, 645, 859, 1640, 2130]
    toa_bands   = (np.array(toa_refs)[band_index,]).tolist()
    view_angles = (np.array(view_ang_names)[band_index,]).tolist()
    sun_angles  = sun_ang_name
    sza         = gdal.Open(sun_angles).ReadAsArray()[1] * 0.01
    scale       = scale / np.cos(np.deg2rad(sza))
    off         = off / np.cos(np.deg2rad(sza))
    aero = solve_aerosol(sensor_sat,toa_bands,band_wv, band_index,view_angles,sun_angles,\
                         obs_time,cloud_mask, gamma=10., ref_scale = scale, ref_off = off, \
                         spec_m_dir=file_path+'/spectral_mapping/', emus_dir=file_path+'/emus/', mcd43_dir=vrt_dir)
    aero._solving()
    toa_bands   = toa_refs
    view_angles = view_ang_names
    base        = base + '/' +'B'.join(toa_bands[0].split('/')[-1].split('B')[:-1])
    aot         = base + 'aot.tif'
    tcwv        = base + 'tcwv.tif'
    tco3        = base + 'tco3.tif'
    aot_unc     = base + 'aot_unc.tif'
    tcwv_unc    = base + 'tcwv_unc.tif'
    tco3_unc    = base + 'tco3_unc.tif'
    rgb = [toa_bands[3], toa_bands[2], toa_bands[1]]
    band_index = [0,1,2,3,4,5,6]
    atmo = atmospheric_correction(sensor_sat,toa_bands, band_index,view_angles,sun_angles, \
                                  aot = aot, cloud_mask = cloud_mask,tcwv = tcwv, tco3 = tco3, \
                                  aot_unc = aot_unc, tcwv_unc = tcwv_unc, tco3_unc = tco3_unc, \
                                  rgb = rgb, ref_scale = scale, ref_off = off, emus_dir=file_path+'/emus/', aoi=aoi)
    atmo._doing_correction()
    return aero, atmo

