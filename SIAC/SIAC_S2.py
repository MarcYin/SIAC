import os
import sys
import json
import shutil
import argparse
import requests
import warnings
warnings.filterwarnings("ignore") 
import numpy as np
from glob import glob
from osgeo import gdal
from SIAC.get_MCD43 import get_mcd43
from datetime import datetime
from SIAC.create_logger import create_logger
from SIAC.s2_preprocessing import s2_pre_processing
from SIAC.the_aerosol import solve_aerosol
from SIAC.the_correction import atmospheric_correction
from SIAC.downloaders import downloader
from SIAC.multi_process import parmap
from os.path import expanduser
from SIAC.raster_boundary import get_boundary
from SIAC.get_VNP43MA1 import download_VNP43MA1

home = expanduser("~")
file_path = os.path.dirname(os.path.realpath(__file__))

def SIAC_S2(s2_t, send_back = False, mcd43 = home + '/MCD43/', vrt_dir = home + '/MCD43_VRT/', aoi = None, 
             global_dem  = None, cams_dir = None, jasmin = False, Gee = True, use_VIIRS = False, do_rgb = True):
    '''
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
    '''
    rets = s2_pre_processing(s2_t, cams_dir, global_dem)
    aero_atmos = []
    for ret in rets:
        ret += (mcd43, vrt_dir, aoi, global_dem, cams_dir, jasmin, Gee, use_VIIRS, do_rgb)
        aero_atmo = do_correction(*ret)
        if send_back:
            aero_atmos.append(aero_atmo)
    if send_back:
        return aero_atmos

def do_correction(sun_ang_name, view_ang_names, toa_refs, cloud_name, \
                  cloud_mask, aot, tcwv, metafile, mcd43 = home + '/MCD43/', \
                  vrt_dir = home + '/MCD43_VRT/', aoi=None, \
                  global_dem  = None, cams_dir = None, jasmin = False, Gee = True, use_VIIRS = False, do_rgb = True
                  ):
    
    log_file = os.path.dirname(metafile) + '/SIAC_S2.log'
    logger = create_logger(log_file)

    jasmin_cams_dir = '/vsicurl/https://gws-access.jasmin.ac.uk/public/nceo_ard/cams/'
    if jasmin:
        if global_dem is None:
            global_dem  = '/work/scratch-pw/marcyin/DEM/global_dem.vrt'
        if cams_dir is None:
            cams_dir    = '/work/scratch-pw/marcyin/CAMS/'     
        os.environ['jasmin_memory_limit'] = '6.4e+10'
    else:
        if global_dem is None:
            # global_dem  = '/vsicurl/http://www2.geog.ucl.ac.uk/~ucfafyi/eles/global_dem.vrt'
            global_dem = '/vsicurl/https://gws-access.jasmin.ac.uk/public/nceo_ard/DEM_V3/global_dem.vrt'
        if cams_dir is None:
            # cams_dir    = '/vsicurl/http://www2.geog.ucl.ac.uk/~ucfafyi/cams/'
            cams_dir = jasmin_cams_dir
            
    if os.path.realpath(mcd43) in os.path.realpath(home + '/MCD43/'):
        if not os.path.exists(home + '/MCD43/'):
            os.mkdir(home + '/MCD43/')

    if os.path.realpath(vrt_dir) in os.path.realpath(home + '/MCD43_VRT/'):
        if not os.path.exists(home + '/MCD43_VRT/'):
            os.mkdir(home + '/MCD43_VRT/')

    logger.debug(f"Path to global dem: {global_dem}")
    logger.debug(f"Path to CAMS directory: {cams_dir}")
    logger.debug(f"Path to MCD43 directory: {mcd43}")
    logger.debug(f"Path to MCD43_VRT directory: {vrt_dir}")

    #base = os.path.dirname(toa_refs[0])
    base = toa_refs[0].replace('B01.jp2', '')
    with open(metafile) as f:
        for i in f.readlines():
            if 'SENSING_TIME' in i:
                sensing_time = i.split('</')[0].split('>')[-1]
                obs_time = datetime.strptime(sensing_time, u'%Y-%m-%dT%H:%M:%S.%fZ')
            if 'TILE_ID' in i:
                sat  = i.split('</')[0].split('>')[-1].split('_')[0]
                tile = i.split('</')[0].split('>')[-1]

    # check if prior dir is readable
    if cams_dir != jasmin_cams_dir:
        prior_dir = os.path.realpath(cams_dir + '/' + datetime.strftime(obs_time, '%Y_%m_%d'))
        dir = gdal.ReadDir(prior_dir)
        if dir is None:
            logger.debug(f"{prior_dir} does not exist.")
            logger.debug(f"Set CAMS directory to {jasmin_cams_dir}")
            cams_dir = jasmin_cams_dir
        dir = None

    
    s2_file_dir = os.path.dirname(metafile)
    PRODUCT_ID = s2_file_dir.split('/')[-3]
    processing_baseline = PRODUCT_ID.split('_')[3][1:]
    

    # getting offset for new processing baseline
    if int(processing_baseline) >= 400:
        
        tile_dir = '/'.join(s2_file_dir.split('/')[:-2])
        
        product_meta = tile_dir + '/MTD_MSIL1C.xml'
        offsets = []
        with open(product_meta) as f:
            for i in f.readlines():
                if 'RADIO_ADD_OFFSET' in i:
                    offset = i.replace('<RADIO_ADD_OFFSET band_id=', '').replace('</RADIO_ADD_OFFSET>', '').replace('"', '').split('>')
                    offset = i.replace('<RADIO_ADD_OFFSET band_id=', '').replace('</RADIO_ADD_OFFSET>', '').replace(' ', '').replace('\n', '').replace('"', '').split('>')
                    offsets.append(offset)
                if 'QUANTIFICATION_VALUE' in i:
                    QUANTIFICATION_VALUE = i.replace('<QUANTIFICATION_VALUE unit="none">', '').replace('</QUANTIFICATION_VALUE>', '').replace(' ', '').replace('\n', '')
                    QUANTIFICATION_VALUE = float(QUANTIFICATION_VALUE)
                    # print(QUANTIFICATION_VALUE)
        offsets = np.array(offsets).astype(int)
        inds = np.argsort(offsets[:,0])
        offsets = offsets[inds]
        offsets = offsets[:, 1]
        QUANTIFICATION_VALUE = np.array([QUANTIFICATION_VALUE] * len(toa_refs))
    else:
        offsets = np.array([0] * len(toa_refs))
        QUANTIFICATION_VALUE = np.array([10000.] * len(toa_refs))
    # L1C_TOAi = (L1C_DNi + RADIO_ADD_OFFSETi) / QUANTIFICATION_VALUEi
    #          = L1C_DNi / QUANTIFICATION_VALUEi  + RADIO_ADD_OFFSETi / QUANTIFICATION_VALUEi
    scale       = 1 / QUANTIFICATION_VALUE
    off         = offsets / QUANTIFICATION_VALUE

    logger.info('Starting atmospheric corretion for %s'%PRODUCT_ID)

    VNP43_fnames_dates = None
    mcd43_gee_folder = None
    mcd43_vrt_dir = None

    # if not np.all(cloud_mask):
        # handlers = logger.handlers[:]
        # for handler in handlers:
        #     handler.close()
        #     logger.removeHandler(handler)
    mcd43_date = datetime(obs_time.year, obs_time.month, obs_time.day)
    
    if Gee:
        from SIAC.MCD43_GEE import get_MCD43_GEE
        logger.info('MODIS BRDF product is chosen')
        logger.info('Getting MCD43 from GEE')
        geojson = get_boundary(toa_refs[0], to_wgs84 = True)[0]
        coords = json.loads(geojson)['features'][0]['geometry']['coordinates']
        mcd43_gee_folder = os.path.dirname(os.path.dirname(toa_refs[0])) + '/MCD43/'
        if not os.path.exists(mcd43_gee_folder):
            os.mkdir(mcd43_gee_folder)
        temporal_window = 16
        get_MCD43_GEE(mcd43_date, coords, temporal_window, mcd43_gee_folder)
    elif use_VIIRS:
        logger.info('VIIRS BRDF product is chosen')
        logger.info('Getting VNP43MA1 from NASA server')
        filenames = download_VNP43MA1(toa_refs[0], mcd43_date, mcd43, temporal_window = 16)
        filenames = np.array(filenames)
        all_dates = np.array([i.split('/')[-1] .split('.')[1][1:9] for i in filenames])          
        udates = np.unique(all_dates)  
        VNP43_fnames_dates =  [[filenames[all_dates==date].tolist(),date] for date in udates]
        
    else:
        logger.info('MODIS BRDF product is chosen')
        logger.info('Getting MCD43 from NASA server')
        mcd43_vrt_dir = get_mcd43(toa_refs[0], mcd43_date, mcd43_dir = mcd43, vrt_dir = vrt_dir, logger = logger, jasmin = jasmin)
            
        #logger = create_logger(log_file)
    # else:
    #     logger.info('No clean pixel in this scene and no MCD43 is downloaded.')
        

    sensor_sat = 'MSI', sat
    band_index  = [1,2,3,7,11,12]
    band_wv    = [469, 555, 645, 859, 1640, 2130]
    toa_bands   = (np.array(toa_refs)[band_index,]).tolist()
    view_angles = (np.array(view_ang_names)[band_index,]).tolist()
    sun_angles  = sun_ang_name
    ref_scale = scale[band_index, None, None]
    ref_off = off[band_index, None, None]
    #logger.info('First pass AOT and TCWV: %.02f, %.02f'%(aot.mean(), tcwv.mean()))
    #logger.info('Running SIAC for tile: %s on %s'%(tile, obs_time.strftime('%Y-%M-%d')))
    aero = solve_aerosol(sensor_sat,toa_bands,band_wv, band_index,view_angles,\
                         sun_angles,obs_time, gamma=10., spec_m_dir= file_path+'/spectral_mapping/',
                         ref_scale = ref_scale, ref_off = ref_off, emus_dir=file_path+'/emus/',\
                         mcd43_dir=mcd43_vrt_dir, aoi=aoi, log_file = log_file, global_dem  = global_dem, cams_dir = cams_dir, \
                         prior_scale = [1., 0.1, 46.698, 1., 1., 1.], mcd43_gee_folder = mcd43_gee_folder,
                         VNP43_fnames_dates = VNP43_fnames_dates
                         )
    example_file = aero._solving(cloud_mask)

    toa_bands  = toa_refs
    aot = base + 'aot.tif'
    tcwv = base + 'tcwv.tif'
    tco3 = base + 'tco3.tif'
    aot_unc = base + 'aot_unc.tif'
    tcwv_unc = base + 'tcwv_unc.tif'
    tco3_unc = base + 'tco3_unc.tif'
    rgb = [toa_bands[3], toa_bands[2], toa_bands[1]]
    band_index = [0,1,2,3,4,5,6,7,8,9,10,11,12]
    view_angles = (np.array(view_ang_names)[band_index,]).tolist()
    ref_scale = scale[band_index, None, None]
    ref_off = off[band_index, None, None]
    

    atmospheric_correction(sensor_sat, toa_bands, band_index, view_angles,\
                                  sun_angles, example_file, aot = aot, \
                                  tcwv = tcwv, tco3 = tco3, aot_unc = aot_unc, \
                                  tcwv_unc = tcwv_unc, tco3_unc = tco3_unc, rgb = rgb,  ref_scale = ref_scale, ref_off = ref_off,\
                                  emus_dir=file_path+'/emus/', log_file = log_file, global_dem  = global_dem, cams_dir = cams_dir,
                                  do_rgb = do_rgb)
    #atmo._doing_correction()

    logger.info('Generating summery Json.')
    print(os.path.dirname(base))
    summeryJson(os.path.dirname(base))

    if not np.all(cloud_mask):
#         if jasmin:
        try:
            shutil.rmtree(mcd43_vrt_dir)
        except:
            pass
    #return aero, atmo

def summeryJson(dest):

    B02 = glob(os.path.join(dest, '*B02_sur.tif'))[0]

    bNames = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']

    header = '/'.join(B02.split('/')[:-5]) + '/'

    imgs = []
    imgUncs = []
    for band in bNames:
        ff = glob(B02.replace('B02_sur.tif', band + '_sur.tif'))
        imgs.append(ff[0].replace(header, ''))

        ff = glob(B02.replace('B02_sur.tif', band + '_sur_unc.tif'))
        imgUncs.append(ff[0].replace(header, ''))

    viewAngles = []
    for band in bNames:
        ff = glob(dest.replace('IMG_DATA', 'ANG_DATA/VAA_VZA_' + band + '.tif'))
        viewAngles.append(ff[0].replace(header, ''))

    ff = glob(dest.replace('IMG_DATA', 'cloud.tif'))
    cloud = ff[0].replace(header, '')

    ff = glob(dest.replace('IMG_DATA', 'ANG_DATA/SAA_SZA.tif'))
    sunAngles = ff[0].replace(header, '')

    boaFull = (dest + '/BOA_RGB.tif').replace(header, '')
    toaFull = (dest + '/TOA_RGB.tif').replace(header, '')

    boaOvrs = [dest.replace(header, '') + '/BOA_ovr_large.png', 
               dest.replace(header, '') + '/BOA_ovr_medium.png', 
               dest.replace(header, '') + '/BOA_ovr_small.png']
    toaOvrs = [dest.replace(header, '') + '/TOA_ovr_large.png', 
               dest.replace(header, '') + '/TOA_ovr_medium.png', 
               dest.replace(header, '') + '/TOA_ovr_small.png']


    atmoBands = ['aot', 'tcwv', 'tco3']
    atmoParas = []
    atmoParasUncs = []

    for band in atmoBands:
        ff = glob(B02.replace('B02_sur.tif', band + '.tif'))
        atmoParas.append(ff[0].replace(header, ''))

        ff = glob(B02.replace('B02_sur.tif', band + '_unc.tif'))
        atmoParasUncs.append(ff[0].replace(header, ''))

    ff = glob(dest.replace('IMG_DATA', 'SIAC_S2.log'))
    siacLog =  ff[0]

    ff = glob(dest + '/AOI.json')
    aoi =  ff[0]

    with open(aoi, 'r') as f:
        txt = json.load(f)
    txt['name'] = 'SIAC outputs'

    with open(siacLog, 'r') as f:
        logstr = f.read().split('\n')
        version = logstr[0].split(' - ')[1]
        for i in logstr:
            if 'Clean pixel percentage' in i:
                CleanPixelPercentage = float(i.split('Clean pixel percentage: ')[1])
            if 'Valid pixel percentage' in i:
                ValidPixelPercentage = float(i.split('Valid pixel percentage: ')[1])

    txt.update({'Version': version, 
                'CleanPixelPercentage': CleanPixelPercentage,
                'ValidPixelPercentage': ValidPixelPercentage
            })
    txt['features'][0].update({'aoi': aoi.replace(header, ''), 
                                'siacLog': siacLog.replace(header, ''),
                                'toaOvrs': toaOvrs,
                                'boaOvrs': boaOvrs,
                                'toaOvrFull': toaFull,
                                'boaOvrFull': boaFull,
                                'viewAngles': viewAngles, 
                                'sunAngles': sunAngles,
                                'SurfaceReflectance': imgs,
                                'SurfaceReflectanceUncertainty': imgUncs,
                                'atmoParas': atmoParas,
                                'atmoParasUncs': atmoParasUncs,
                                'cloud': cloud
                            })
    
    s2_tile_dir = '/'.join(dest.split('/')[:-3])
    with open(s2_tile_dir + '/siac_output.json', 'w') as f:
        json.dump(txt, f, ensure_ascii=False, indent=4)
        
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

