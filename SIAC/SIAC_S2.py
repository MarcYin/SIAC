import os
import sys
import numpy as np
from get_MCD43 import get_mcd43
from datetime import datetime
from the_aerosol import solve_aerosol
from the_correction import atmospheric_correction
#sys.path.insert(0, '/data/store01/data_dirs/students/ucfafyi/Atmospheric_correction/atmospheric_correction')
from s2_preprocessing import s2_pre_processing

def SIAC_S2(s2_t, send_back = False):
    rets = s2_pre_processing(s2_t)
    aero_atmos = []
    for ret in rets:
        #sun_ang_name, view_ang_names, toa_refs, cloud_name, cloud_mask, metafile = ret
        aero_atmo = do_correction(*ret)
        if send_back:
            aero_atmos.append(aero_atmo)
    if send_back:
        return aero_atmos

def do_correction(sun_ang_name, view_ang_names, toa_refs, cloud_name, cloud_mask, metafile):
    base = os.path.dirname(toa_refs[0])
    base = toa_refs[0].replace('B01.jp2', '')
    with open(metafile) as f:
        for i in f.readlines():
            if 'SENSING_TIME' in i:
                sensing_time = i.split('</')[0].split('>')[-1]
                obs_time = datetime.strptime(sensing_time, u'%Y-%m-%dT%H:%M:%S.%fZ')
            if 'TILE_ID' in i:
                sat = i.split('</')[0].split('>')[-1].split('_')[0]
    get_mcd43(toa_refs[0], obs_time, mcd43_dir = '/data/nemesis/MCD43/', vrt_dir = '/home/ucfafyi/DATA/Multiply/MCD43/')
    #get_mcd43(toa_refs[0], obs_time, mcd43_dir = '/home/ucfafyi/hep/MCD43/', vrt_dir = '/home/ucfafyi/DATA/Multiply/MCD43/')
    sensor_sat = 'MSI', sat
    band_index  = [1,2,3,7,11,12]
    band_wv    = [469, 555, 645, 859, 1640, 2130]
    toa_bands   = (np.array(toa_refs)[band_index,]).tolist()
    view_angles = (np.array(view_ang_names)[band_index,]).tolist()
    sun_angles  = sun_ang_name
    aero = solve_aerosol(sensor_sat,toa_bands,band_wv, band_index,view_angles,sun_angles,obs_time,cloud_mask, gamma=10.)
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
    atmo = atmospheric_correction(sensor_sat,toa_bands, band_index,view_angles,sun_angles, aot = aot, cloud_mask = cloud_mask,\
                                  tcwv = tcwv, tco3 = tco3, aot_unc = aot_unc, tcwv_unc = tcwv_unc, tco3_unc = tco3_unc, rgb = rgb)
    atmo._doing_correction()
    return aero, atmo

