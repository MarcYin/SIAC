#/usr/bin/env python
import os
import sys
import inspect
from osgeo import gdal
import numpy as np
from glob import glob
from tensorflow import keras
from scipy.ndimage import binary_erosion, binary_dilation

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath) 

model = keras.models.load_model(myPath + '/data/cloud_segmentation_3.h5')

def get_cloud_prob(toa_data):
    
    xsize, ysize = toa_data.shape[1:]
    nxs = np.ceil(xsize/200).astype(int)
    nys = np.ceil(ysize/200).astype(int)
    inds = np.indices((nxs, nys)).reshape(2, -1).T

    cloud_probs = np.zeros((nxs*200, nys*200))
    shadow_probs = np.zeros((nxs*200, nys*200))

    for xx, yy in inds:
        xoff = int(xx * 200)
        yoff = int(yy * 200)
        
        sub = np.zeros((6, 256, 256))

        xstart = np.max([0,     xoff - 28])
        xend   = np.min([xsize, xoff + 200 + 28])
        
        ystart = np.max([0,     yoff - 28])
        yend   = np.min([ysize, yoff + 200 + 28]) 
        
        sub_data = toa_data[:, xstart: xend, ystart: yend]
        sub[:, :sub_data.shape[1], :sub_data.shape[2]] = sub_data
        
        if sub.sum() == 0:
            clear_prob, cloud_prob, shadow_prob = np.zeros((256, 256)), np.zeros((256, 256)), np.zeros((256, 256))
        else:
            clear_prob, cloud_prob, shadow_prob = model.predict(sub.transpose(1,2,0)[None,])[0].transpose(2, 0, 1)
        
        if xoff - 28 < 0:
            xcut = 0
        else:
            xcut = 28
        if yoff - 28 < 0:
            ycut = 0
        else:
            ycut = 28    
        cloud_probs[xoff: xoff + 200, yoff: yoff + 200] = cloud_prob[xcut:xcut+200, ycut:ycut+200]
        shadow_probs[xoff: xoff + 200, yoff: yoff + 200] = shadow_prob[xcut:xcut+200, ycut:ycut+200]

    cloud_probs = cloud_probs[:xsize, :ysize]
    shadow_probs = shadow_probs[:xsize, :ysize]
    return cloud_probs, shadow_probs

def read_l8_cloud_bands(toa_dir, res = 30):
    l8_bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
    toa_bands = [glob(toa_dir + '/*%s.TIF'%band)[0] for band in l8_bands]
    header = toa_bands[0].replace('B2.TIF', '')
    with open( header + 'MTL.txt', 'r') as f:
        txt = f.readlines()

    toa_data = []
    for band in l8_bands:
        for line in txt:
            if 'REFLECTANCE_MULT_BAND_%s'%band[1] in line:
                scale = float(line.split(' = ')[1].strip())
            if 'REFLECTANCE_ADD_BAND_%s'%band[1] in line:
                offset = float(line.split(' = ')[1].strip())
            if 'SUN_ELEVATION' in line:
                sza = 90 - float(line.split(' = ')[1].strip())
        toa_band = header + band + '.TIF'
        data = (gdal.Warp('', toa_band, xRes= res, yRes = res, format='MEM').ReadAsArray() * scale + offset) /  np.cos(np.deg2rad(sza))
        toa_data.append(data)

    toa_data = np.array(toa_data) 
    toa_data[toa_data<=0] = 0
    return toa_data


def read_s2_cloud_bands_old(toa_dir, res = 20):
    s2_bands = ['B02', 'B03', 'B04', 'B8A', 'B11', 'B12']
    toa_bands = [glob(toa_dir + '/*%s.jp2'%band)[0] for band in s2_bands]
    header = toa_bands[0].replace('B02.jp2', '')
    toa_data = []
    for toa_band in toa_bands:
        data = gdal.Warp('', toa_band, xRes= res, yRes = res, format='MEM', resampleAlg=0).ReadAsArray() / 10000.
        toa_data.append(data)
    toa_data = np.array(toa_data) 
    toa_data[toa_data<=0] = 0
    return toa_data

def read_s2_cloud_bands(toa_dir, res = 20):
    s2_bands = ['B02', 'B03', 'B04', 'B8A', 'B11', 'B12']
    band_inds = [1, 2, 3, 8, 11, 12]

    toa_dir = toa_dir.replace('/IMG_DATA/', '/IMG_DATA')
    s2_file_dir = os.path.dirname(toa_dir)
    PRODUCT_ID = s2_file_dir.split('/')[-3]
    
    # PRODUCT_ID = [i for i in toa_dir.split('/') if ('_MSIL1C_' in i) & (i.endswith('.SAFE'))][0]
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
        offsets = offsets[band_inds, 1]
        QUANTIFICATION_VALUE = np.array([QUANTIFICATION_VALUE] * len(s2_bands))
    else:
        offsets = np.array([0] * len(s2_bands))
        QUANTIFICATION_VALUE = np.array([10000.] * len(s2_bands))
    # L1C_TOAi = (L1C_DNi + RADIO_ADD_OFFSETi) / QUANTIFICATION_VALUEi
    #          = L1C_DNi / QUANTIFICATION_VALUEi  + RADIO_ADD_OFFSETi / QUANTIFICATION_VALUEi
    scale       = 1 / QUANTIFICATION_VALUE
    off         = offsets / QUANTIFICATION_VALUE

    toa_bands = [glob(toa_dir + '/*%s.jp2'%band)[0] for band in s2_bands]
    header = toa_bands[0].replace('B02.jp2', '')
    toa_data = []
    for i, toa_band in enumerate(toa_bands):
        data = gdal.Warp('', toa_band, xRes= res, yRes = res, format='MEM', resampleAlg=0).ReadAsArray() * scale[i] + off[i]
        toa_data.append(data)
    toa_data = np.array(toa_data) 
    toa_data[toa_data<=0] = 0
    return toa_data
    


def get_cloud_and_shadow_mask(cloud_probs, shadow_probs):
    cloud = cloud_probs > 0.3 
    shadow = shadow_probs > 0.5

    cloud = binary_dilation(cloud, iterations = 10)
    shadow = binary_dilation(shadow, iterations = 10)
    return cloud, shadow

def save_mask(toa_dir, cloud_and_shadow, example_file, res):
    
    outputFileName = toa_dir + '/CNN_cloud_and_shadow_mask.tif'
    g = gdal.Open(example_file)
    projectionRef = g.GetProjectionRef()
    geotransform  = g.GetGeoTransform()
    
    print(geotransform)
    geotransform = list(geotransform)
    geotransform[1] =  res
    geotransform[5] = -res
    geotransform = tuple(geotransform)

    n_band, nx, ny = cloud_and_shadow.shape
    if os.path.exists(outputFileName):
        os.remove(outputFileName)
    dst_ds = gdal.GetDriverByName('GTiff').Create(outputFileName, ny, nx, n_band, gdal.GDT_Byte, options=["TILED=YES", "COMPRESS=DEFLATE"])
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(projectionRef)
    
    cloud_and_shadow = cloud_and_shadow.astype(np.byte)
    dst_ds.GetRasterBand(1).WriteArray(cloud_and_shadow[0])
    dst_ds.GetRasterBand(2).WriteArray(cloud_and_shadow[1])
    dst_ds.FlushCache()                  
    dst_ds = None
    return outputFileName

def save_prob(toa_dir, cloud_and_shadow_probs, example_file, res):
    
    outputFileName = toa_dir + '/CNN_cloud_and_shadow_prob.tif'
    g = gdal.Open(example_file)
    projectionRef = g.GetProjectionRef()
    geotransform  = g.GetGeoTransform()
    
    print(geotransform)
    geotransform = list(geotransform)
    geotransform[1] =  res
    geotransform[5] = -res
    geotransform = tuple(geotransform)

    scale = 100
    n_band, nx, ny = cloud_and_shadow_probs.shape
    if os.path.exists(outputFileName):
        os.remove(outputFileName)
    dst_ds = gdal.GetDriverByName('GTiff').Create(outputFileName, ny, nx, n_band, gdal.GDT_Byte, options=["TILED=YES", "COMPRESS=DEFLATE"])
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(projectionRef)
    
    cloud_and_shadow_probs = np.maximum(cloud_and_shadow_probs, 0)
    cloud_and_shadow_probs = np.minimum(cloud_and_shadow_probs, 1)

    cloud_and_shadow_probs = (cloud_and_shadow_probs * scale).astype(np.byte)

    dst_ds.GetRasterBand(1).WriteArray(cloud_and_shadow_probs[0])
    dst_ds.GetRasterBand(2).WriteArray(cloud_and_shadow_probs[1])
    dst_ds.FlushCache()                  
    dst_ds = None
    return outputFileName


def get_s2_cloud_old(toa_dir, do_plot = True):
    
    res = 20
    
    toa_data = read_s2_cloud_bands(toa_dir, res = res)
    cloud_probs, shadow_probs = get_cloud_prob(toa_data)
    cloud, shadow = get_cloud_and_shadow_mask(cloud_probs, shadow_probs)
    

    example_file = glob(toa_dir + '/*%s.jp2'%'B8A')[0]

    cloud_and_shadow_probs = np.array([cloud_probs, shadow_probs])
    save_prob(toa_dir, cloud_and_shadow_probs, example_file, res)
    
    cloud_and_shadow = np.array([cloud, shadow])
    save_mask(toa_dir, cloud_and_shadow, example_file, res)

    if do_plot:
        cloud = cloud * 1.
        cloud[cloud < 1] = np.nan
        shadow = shadow * 1.
        shadow[shadow < 1] = np.nan
        
        import matplotlib
        matplotlib.use('Agg')
        import pylab as plt
        plt.figure(figsize=(40,40))
        plt.imshow(toa_data[[2,1,0]].transpose(1,2,0)*4)
        plt.imshow(cloud, alpha=0.8, cmap = plt.cm.cool)
        plt.imshow(shadow, alpha=0.8, cmap = plt.cm.winter)
        plt.savefig(toa_dir + '/cloud_and_shadow_mask.png', bbox_inches='tight', pad_inches=0)


def get_s2_cloud(toa_dir, do_plot = True, res = 20):
    
    toa_dir = glob(f'{toa_dir}/GRANULE/*/IMG_DATA')[0]
    
    toa_data = read_s2_cloud_bands(toa_dir, res = res)
    cloud_probs, shadow_probs = get_cloud_prob(toa_data)
    cloud, shadow = get_cloud_and_shadow_mask(cloud_probs, shadow_probs)
    

    example_file = glob(toa_dir + '/*%s.jp2'%'B8A')[0]

    cloud_and_shadow_probs = np.array([cloud_probs, shadow_probs])
    save_prob(toa_dir, cloud_and_shadow_probs, example_file, res)
    
    cloud_and_shadow = np.array([cloud, shadow])
    save_mask(toa_dir, cloud_and_shadow, example_file, res)

    if do_plot:
        cloud = cloud * 1.
        cloud[cloud < 1] = np.nan
        shadow = shadow * 1.
        shadow[shadow < 1] = np.nan
        
        import matplotlib
        matplotlib.use('Agg')
        import pylab as plt
        plt.figure(figsize=(40,40))
        plt.imshow(toa_data[[2,1,0]].transpose(1,2,0)*4)
        plt.imshow(cloud, alpha=0.8, cmap = plt.cm.cool)
        plt.imshow(shadow, alpha=0.8, cmap = plt.cm.winter)
        plt.savefig(toa_dir + '/cloud_and_shadow_mask.png', bbox_inches='tight', pad_inches=0)
        
def get_l8_cloud(toa_dir, do_plot = True):
    
    res = 30
    toa_data = read_l8_cloud_bands(toa_dir, res = res)
    cloud_probs, shadow_probs = get_cloud_prob(toa_data)
    cloud, shadow = get_cloud_and_shadow_mask(cloud_probs, shadow_probs)
    
    example_file = glob(toa_dir + '/*%s.TIF'%'B2')[0]

    cloud_and_shadow_probs = np.array([cloud_probs, shadow_probs])
    save_prob(toa_dir, cloud_and_shadow_probs, example_file, res)
    
    cloud_and_shadow = np.array([cloud, shadow])
    save_mask(toa_dir, cloud_and_shadow, example_file, res)
    
    if do_plot:
        import matplotlib
        matplotlib.use('Agg')
        import pylab as plt
        cloud = cloud * 1.
        cloud[cloud < 1] = np.nan
        shadow = shadow * 1.
        shadow[shadow < 1] = np.nan
        plt.figure(figsize=(40,40))
        plt.imshow(toa_data[[2,1,0]].transpose(1,2,0)*4)
        plt.imshow(cloud, alpha=0.8, cmap = plt.cm.cool)
        plt.imshow(shadow, alpha=0.8, cmap = plt.cm.winter)
        plt.savefig(toa_dir + '/cloud_and_shadow_mask.png', bbox_inches='tight', pad_inches=0)



# toa_dir = '/gws/nopw/j04/nceo_ard/public/S2/13/S/CR/S2B_MSIL1C_20190114T174659_N0207_R098_T13SCR_20190114T194448.SAFE' 
# get_s2_cloud(toa_dir, do_plot = True)




#t_dir = '/gws/nopw/j04/nceo_ard/public/L8/176/036/LC08_L1TP_176036_20190414_20190422_01_T1'
#t_dir = '/gws/nopw/j04/nceo_ard/public/L8/006/068/LC08_L1TP_006068_20190829_20190903_01_T1'
#get_l8_cloud(t_dir, do_plot = True)

