#!/usr/bin/env python 
import os
import sys
from osgeo import ogr, gdal
import psutil
import errno 
import numpy as np
from glob import glob
from shutil import copyfile
from functools import partial
from multiprocessing import Pool
import xml.etree.ElementTree as ET
from scipy.ndimage import binary_dilation
#from scipy.interpolate import griddata

def parse_xml(meta_file, example_file, sun_ang_name):
    tree = ET.parse(meta_file)
    root = tree.getroot()
    saa =[]
    sza =[]
    msz = []
    msa = []
    vza = {}
    vaa = {}
    mvz = {}
    mva = {}
    for child in root:
        for j in child:
            for k in j.findall('Sun_Angles_Grid'):
                for l in k.findall('Zenith'):
                    for m in l.findall('Values_List'):
                        for x in m.findall('VALUES'):
                            sza.append(x.text.split())
     
                for n in k.findall('Azimuth'):
                    for o in n.findall('Values_List'):
                        for p in o.findall('VALUES'):
                            saa.append(p.text.split())
            for ms in j.findall('Mean_Sun_Angle'):
                msz = float(ms.find('ZENITH_ANGLE').text)
                msa = float(ms.find('AZIMUTH_ANGLE').text)
     
            for k in j.findall('Viewing_Incidence_Angles_Grids'):
                for l in k.findall('Zenith'):
                    for m in l.findall('Values_List'):
                        vza_sub = []
                        for x in m.findall('VALUES'):
                            vza_sub.append(x.text.split())
                        bi, di, angles = k.attrib['bandId'], \
                                         k.attrib['detectorId'], np.array(vza_sub).astype(float)
                        vza[(int(bi),int(di))] = angles
     
                for n in k.findall('Azimuth'):
                    for o in n.findall('Values_List'):
                        vaa_sub = []
                        for p in o.findall('VALUES'):
                            vaa_sub.append(p.text.split())
                        bi, di, angles = k.attrib['bandId'],\
                                         k.attrib['detectorId'], np.array(vaa_sub).astype(float)
                        vaa[(int(bi),int(di))] = angles
     
            for mvia in j.findall('Mean_Viewing_Incidence_Angle_List'):
                for i in mvia.findall('Mean_Viewing_Incidence_Angle'):
                    mvz[int(i.attrib['bandId'])] = float(i.find('ZENITH_ANGLE').text)
                    mva[int(i.attrib['bandId'])] = float(i.find('AZIMUTH_ANGLE').text)

    sza  = np.array(sza).astype(float)
    saa  = np.array(saa).astype(float)
    saa[saa>180] = saa[saa>180] - 360                 
    g                = gdal.Open(example_file)
    geo              = g.GetGeoTransform()
    projection       = g.GetProjection()
    geotransform     = (geo[0], 5000, geo[2], geo[3], geo[4], -5000)
    if os.path.exists(sun_ang_name):
        os.remove(sun_ang_name)
    dst_ds = gdal.GetDriverByName('GTiff').Create(sun_ang_name, 23, 23, 2, gdal.GDT_Int16, options=["TILED=YES", "COMPRESS=DEFLATE"])
    dst_ds.SetGeoTransform(geotransform)   
    dst_ds.SetProjection(projection) 
    dst_ds.GetRasterBand(1).WriteArray((saa * 100).astype(int))
    dst_ds.GetRasterBand(2).WriteArray((sza * 100).astype(int))
    dst_ds, g = None, None
    return vaa, vza, mvz, mva

def get_mean_angle(meta_file):
    tree = ET.parse(meta_file)
    root = tree.getroot()
    mvz = {} 
    mva = {} 
    msz = []
    msa = []
    for child in root: 
        for j in child: 
            for mvia in j.findall('Mean_Viewing_Incidence_Angle_List'): 
                for i in mvia.findall('Mean_Viewing_Incidence_Angle'): 
                    mvz[int(i.attrib['bandId'])] = float(i.find('ZENITH_ANGLE').text) 
                    mva[int(i.attrib['bandId'])] = float(i.find('AZIMUTH_ANGLE').text) 
            for ms in j.findall('Mean_Sun_Angle'):
                    msz = float(ms.find('ZENITH_ANGLE').text)
                    msa = float(ms.find('AZIMUTH_ANGLE').text)
    return msz, msa, mvz, mva

def _get_angle(view_ang_name_gml, vaa, vza, band_dict):
    band_name, view_ang_name, gml = view_ang_name_gml
    if gml[-7:-4] != 'L1C':
        band_ind = band_dict[gml[-7:-4]]
    else:
        band_ind = band_dict[gml[-14:-11]]
    _va = np.nanmean([vaa[i] for i in vaa.keys() if i[0] == band_ind])
    _vz = np.nanmean([vza[i] for i in vza.keys() if i[0] == band_ind])
    if np.isnan(_va) or np.isnan(_vz):
        return False   
    g = ogr.Open(gml)
    xRes = 10; yRes=10
    g1     = gdal.Open(band_name)
    geo_t = g1.GetGeoTransform()
    x_size, y_size = g1.RasterXSize, g1.RasterYSize
    vas = np.zeros((y_size, x_size), dtype=np.float32)
    vas[:] = -327.67
    vzs = np.zeros((y_size, x_size), dtype=np.float32)
    vzs[:] = -327.67
    x_min, x_max  = min(geo_t[0], geo_t[0] + x_size * geo_t[1]), \
      max(geo_t[0], geo_t[0] + x_size * geo_t[1])
    y_min, y_max  = min(geo_t[3], geo_t[3] + y_size * geo_t[5]), \
      max(geo_t[3], geo_t[3] + y_size * geo_t[5])
    xRes, yRes  = abs(geo_t[1]), abs(geo_t[5])
    x_scale = 5000. / xRes
    y_scale = 5000. / yRes
    
    foot1 = None                            
    foot2 = None                                   
    va1   = None
    vz1   = None         
    va2   = None                                                     
    vz2   = None     
    try:
        layer = g.GetLayer()
        dets = []
        for i in range(layer.GetFeatureCount()): 
            dets.append(layer.GetFeature(i).items()['gml_id'])
        dets = sorted(dets, key = lambda i: int(i.split('-')[2]))
        for i in range(len(dets)):
            det = dets[i]
            foot1 = gdal.Rasterize("", gml, format="MEM", xRes=xRes, yRes=yRes, \
                                   where="gml_id='%s'"%det, outputBounds=[ x_min, y_min, x_max, y_max], \
                                   noData=0, burnValues=1, outputType=gdal.GDT_Byte).ReadAsArray()
            foot1 = binary_dilation(foot1)
            key =  band_dict[det.split('-')[-3]], int(det.split('-')[-2])
            va1 = vaa[key]       
            vz1 = vza[key]                                                
            if i>0:
                overlap = foot1 * foot2
                if overlap.sum() < 10:
                    foot1 = foot2
                else:
                    x,y = np.where(overlap)
                    xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(),y.max()
                    ll = x[x==xmax][-1], y[x==xmax][-1]
                    lr = x[y==ymax][-1], y[y==ymax][-1]
                    ul = x[y==ymin][0], y[y==ymin][0]
                    ur = x[x==xmin][0], y[x==xmin][0]
                    p1 = np.mean([lr, ur], axis=0)
                    p2 = np.mean([ll, ul], axis=0)
                    x1,y1 = np.where(foot2)
                    vamax, vamin  = np.nanmax(va2), np.nanmin(va2)
                    vzmax, vzmin  = np.nanmax(vz2), np.nanmin(vz2)
                    if not (p1==p2).all():
                        p = np.poly1d(np.polyfit([p1[1], p2[1]],[p1[0], p2[0]],1))
                        foot2[x[x > p(y)], y[x > p(y)]] = False
                        minx, miny = np.where(va2 == vamin)
                        maxx, maxy = np.where(va2 == vamax)
                        min_dst = abs(p(miny*y_scale)-minx*x_scale)/(np.sqrt(1+p.c[0]**2))
                        max_dst = abs(p(maxy*y_scale)-maxx*x_scale)/(np.sqrt(1+p.c[0]**2))
                        if  (max_dst < min_dst).any():
                            tmp1 = vamin.copy()
                            vamin = vamax 
                            vamax = tmp1
                        minx, miny = np.where(vz2 == vzmin)
                        maxx, maxy = np.where(vz2 == vzmax)
                        min_dst = abs(p(miny*y_scale)-minx*x_scale)/(np.sqrt(1+p.c[0]**2))
                        max_dst = abs(p(maxy*y_scale)-maxx*x_scale)/(np.sqrt(1+p.c[0]**2))
                        if (max_dst < min_dst).any():
                            tmp2 = vzmin.copy()
                            vzmin = vzmax
                            vzmax = tmp2 
                        dist = abs(p(y1)-x1)/(np.sqrt(1+p.c[0]**2))
                        vas[x1,y1] = vamin + dist/(dist.max()-dist.min()) * (vamax-vamin)
                        vzs[x1,y1] = vzmin + dist/(dist.max()-dist.min()) * (vzmax-vzmin)
                    else:
                        vas[x1,y1] = vamin
                        vzs[x1,y1] = vzmin
                    x1,y1 = np.where(foot1)
                    if i == layer.GetFeatureCount()-1:
                        vamax, vamin  = np.nanmax(va1), np.nanmin(va1)
                        vzmax, vzmin  = np.nanmax(vz1), np.nanmin(vz1) 
                        if not (p1==p2).all():
                            foot1[x[x <= p(y)], y[x <= p(y)]] = False
                            minx, miny = np.where(va1 == vamin)         
                            maxx, maxy = np.where(va1 == vamax)
                            min_dst = abs(p(miny*y_scale)-minx*x_scale)/(np.sqrt(1+p.c[0]**2))
                            max_dst = abs(p(maxy*y_scale)-maxx*x_scale)/(np.sqrt(1+p.c[0]**2))
                            if  (max_dst < min_dst).any():                      
                                tmp1 = vamin.copy()                     
                                vamin = vamax                           
                                vamax = tmp1                            
                            minx, miny = np.where(vz1 == vzmin)         
                            maxx, maxy = np.where(vz1 == vzmax)         
                            min_dst = abs(p(miny*y_scale)-minx*x_scale)/(np.sqrt(1+p.c[0]**2))
                            max_dst = abs(p(maxy*y_scale)-maxx*x_scale)/(np.sqrt(1+p.c[0]**2))
                            if (max_dst < min_dst).any():                       
                                tmp2 = vzmin.copy()                     
                                vzmin = vzmax                           
                                vzmax = tmp2
                            dist = abs(p(y1)-x1)/(np.sqrt(1+p.c[0]**2))
                            vas[x1,y1] = vamin + dist/(dist.max()-dist.min()) * (vamax-vamin)
                            vzs[x1,y1] = vzmin + dist/(dist.max()-dist.min()) * (vzmax-vzmin)
                        else:
                            vas[x1,y1] = vamin 
                            vas[x1,y1] = vamin 
            foot2 = foot1
            va2   = va1
            vz2   = vz1
        #    vas[:] = np.nanmean(vaa.values())
        #    vzs[:] = np.nanmean(vza.values())
        # mask      = vas < -180
        # if (~mask).sum()<1:
        #     vas[:] = np.nanmean(va1)
        #     #vas = fill_bad(vas, ~mask)
        # mask      = vzs < 0                              
        # if (~mask).sum()<1:
        #     vzs[:] = np.nanmean(vz1)
            #vzs = fill_bad(vzs, ~mask)
        #vas[(vas>180) & (vas<=360)] = vas[(vas>180) & (vas<=360)].mean() - 360
        #vas[(vas>=0)  & (vas<=180)] = vas[(vas>=0)  & (vas<=180)].mean()
        # if os.path.exists(view_ang_name):                   
        #     os.remove(view_ang_name)                        
        # dst_ds = gdal.GetDriverByName('GTiff').Create(view_ang_name, x_size, y_size, 2, gdal.GDT_Int16, options=["TILED=YES", "COMPRESS=DEFLATE"])
        # dst_ds.SetGeoTransform(g1.GetGeoTransform())         
        # dst_ds.SetProjection(g1.GetProjection())             
        # dst_ds.GetRasterBand(1).WriteArray((vas * 100).astype(int))            
        # dst_ds.GetRasterBand(2).WriteArray((vzs * 100).astype(int))            
        # dst_ds.GetRasterBand(1).SetNoDataValue(-32767)       
        # dst_ds.GetRasterBand(2).SetNoDataValue(-32767)       
        # dst_ds.FlushCache()                                  
        # dst_ds = None  
        # g1 = None   
        # return True  
    except:
        band_ind = band_dict[gml[-7:-4]]
        vas[:] = np.nanmean([vaa[i] for i in vaa.keys() if i[0] == band_ind])
        vzs[:] = np.nanmean([vza[i] for i in vza.keys() if i[0] == band_ind])
    
    if os.path.exists(view_ang_name):                   
        os.remove(view_ang_name)                        
    dst_ds = gdal.GetDriverByName('GTiff').Create(view_ang_name, x_size, y_size, 2, gdal.GDT_UInt16, options=["TILED=YES", "COMPRESS=DEFLATE"])
    dst_ds.SetGeoTransform(g1.GetGeoTransform())         
    dst_ds.SetProjection(g1.GetProjection())             
    
    mask      = vas < -180
    if (~mask).sum()<1:
        vas[:] = np.nanmean(va1)
        #vas = fill_bad(vas, ~mask)
    mask      = vzs < 0                              
    if (~mask).sum()<1:
        vzs[:] = np.nanmean(vz1)
        #vzs = fill_bad(vzs, ~mask)
    #vas[(vas>180) & (vas<=360)] = vas[(vas>180) & (vas<=360)].mean() - 360
    #vas[(vas>=0)  & (vas<=180)] = vas[(vas>=0)  & (vas<=180)].mean()
    dst_ds.GetRasterBand(1).WriteArray((vas * 100).astype(int))            
    dst_ds.GetRasterBand(2).WriteArray((vzs * 100).astype(int))            
    # dst_ds.GetRasterBand(1).SetNoDataValue(-32767)       
    # dst_ds.GetRasterBand(2).SetNoDataValue(-32767)       
    dst_ds.FlushCache()                                  
    dst_ds = None  
    g1 = None  
    return True


def get_angle(toa_view_ang_name, mva, mvz):
    band_name, view_ang_name = toa_view_ang_name
    band = band_name[-7:-4]
    try:
        _va = mva[band]
        _vz = mvz[band]
    except:
        return False
    if np.isnan(_va) or np.isnan(_vz):
        return False   
    g1     = gdal.Open(band_name)
    x_size, y_size = g1.RasterXSize, g1.RasterYSize
    vaa = np.zeros((y_size, x_size), dtype=np.float32)
    vaa[:] = _va
    vza = np.zeros((y_size, x_size), dtype=np.float32)
    vza[:] = _vz
    
    vaa[vaa < 0] = vaa[vaa < 0] + 360
    vza[vza < 0] = vza[vza < 0] + 360
    
    if os.path.exists(view_ang_name):                   
        os.remove(view_ang_name)                        
    dst_ds = gdal.GetDriverByName('GTiff').Create(view_ang_name, x_size, y_size, 2, gdal.GDT_UInt16, options=["TILED=YES", "COMPRESS=DEFLATE"])
    dst_ds.SetGeoTransform(g1.GetGeoTransform())         
    dst_ds.SetProjection(g1.GetProjection())             
    
    #vzs = fill_bad(vzs, ~mask)
    #vas[(vas>180) & (vas<=360)] = vas[(vas>180) & (vas<=360)].mean() - 360
    #vas[(vas>=0)  & (vas<=180)] = vas[(vas>=0)  & (vas<=180)].mean()
    dst_ds.GetRasterBand(1).WriteArray((vaa * 100).astype(int))            
    dst_ds.GetRasterBand(2).WriteArray((vza * 100).astype(int))            
    # dst_ds.GetRasterBand(1).SetNoDataValue(65535)
    # dst_ds.GetRasterBand(2).SetNoDataValue(65535)
    dst_ds.FlushCache()                                  
    dst_ds = None  
    g1 = None  
    return True

def save_mean_view_angs(mean_vza, mean_vaa, example_file, mean_view_angle_name):
    
    g1     = gdal.Open(example_file)
    x_size, y_size = g1.RasterXSize, g1.RasterYSize
    vaa = np.zeros((y_size, x_size), dtype=np.float32)
    vaa[:] = mean_vaa
    vza = np.zeros((y_size, x_size), dtype=np.float32)
    vza[:] = mean_vza
    
    vaa[vaa < 0] = vaa[vaa < 0] + 360
    vza[vza < 0] = vza[vza < 0] + 360

    if os.path.exists(mean_view_angle_name):                   
        os.remove(mean_view_angle_name)                        
    dst_ds = gdal.GetDriverByName('GTiff').Create(mean_view_angle_name, x_size, y_size, 2, gdal.GDT_UInt16, options=["TILED=YES", "COMPRESS=DEFLATE"])
    dst_ds.SetGeoTransform(g1.GetGeoTransform())         
    dst_ds.SetProjection(g1.GetProjection())             
    dst_ds.GetRasterBand(1).WriteArray((vaa * 100).astype(int))            
    dst_ds.GetRasterBand(2).WriteArray((vza * 100).astype(int))            
    # dst_ds.GetRasterBand(1).SetNoDataValue(-32767)       
    # dst_ds.GetRasterBand(2).SetNoDataValue(-32767)       
    dst_ds.FlushCache()                                  
    dst_ds = None  
    g1 = None  

'''
def fill_bad(array, mask):                        
    x_shp, y_shp = array.shape                     
    valid = np.array(np.where(mask)).T             
    value = array[mask]                            
    mesh  = np.repeat(range(x_shp), y_shp).reshape(x_shp, y_shp), \
            np.tile  (range(y_shp), x_shp).reshape(x_shp, y_shp)
    array = griddata(valid, value, mesh, method='nearest')
    return array
'''
def resample_s2_angles_N0400(metafile):
    #check the available rams and decide cores can be used
    #start multiprocessing
    bands    = ['B01', 'B02', 'B03','B04','B05' ,'B06', 'B07', 'B08','B8A', 'B09', 'B10', 'B11', 'B12'] #all bands
    band_ram = 5e9
    av_ram = psutil.virtual_memory().available
    procs = np.min([int(av_ram / band_ram), psutil.cpu_count(), 4])
    if procs < 1:
        raise MemoryError('At least 500MB ram is needed.')
    s2_file_dir = os.path.dirname(metafile)
    PRODUCT_ID = s2_file_dir.split('/')[-3]
    processing_baseline = PRODUCT_ID.split('_')[3][1:]

    B2 = glob(s2_file_dir + '/IMG_DATA/*B02.jp2')[0]

    if int(processing_baseline) >= 400:
        if not os.path.exists(s2_file_dir + '/ANG_DATA/'):
            os.mkdir(s2_file_dir + '/ANG_DATA/')
        sun_ang_name = s2_file_dir + '/ANG_DATA/' + 'SAA_SZA.tif'
        view_ang_names = [s2_file_dir + '/ANG_DATA/' + 'VAA_VZA_%s.tif'%band for band in bands]
        toa_refs = [B2.replace('B02.jp2', band + '.jp2') for band in bands]
        mean_view_angle_name = s2_file_dir + '/ANG_DATA/' + 'Mean_VAA_VZA.tif'
    else:
        raise IOError('Only processing base line higher than 4.0 can use this function to generate angles')

    cloud_name = s2_file_dir+'/cloud.tif'
    example_file = toa_refs[1]
    
    vaa, vza, mvz, mva = parse_xml(metafile, example_file, sun_ang_name)
    
    mvz = [mvz[i] for i in range(13)]
    mvz = dict(zip(bands, mvz))

    mva = [mva[i] for i in range(13)]
    mva = dict(zip(bands, mva))


    mean_vza = np.mean(list(mvz.values()))
    mean_vaa = np.mean(list(mva.values()))
    
    save_mean_view_angs(mean_vza, mean_vaa, example_file, mean_view_angle_name)
    
    # some gmls may lost....
    toa_view_ang_names = list(zip(np.array(toa_refs), np.array(view_ang_names)))
    par = partial(get_angle, mva=mva, mvz=mvz)
    # p = Pool(procs)
    # #print(view_ang_name_gmls)
    # ret = p.map(par,  view_ang_name_gmls)
    ret  =list(map(par,  toa_view_ang_names))
    # #print(ret)
    # p.close()
    # p.join()
    ret = np.array(ret)
    view_ang_names = np.array(view_ang_names)
 
    if (ret.sum()<13) & (ret.sum()>0):
        inds = np.array(inds)[ret] # this is the good proccessed
        ret = np.zeros(13)
        ret[inds] = 1 # recover them to the original bands

    if ret.sum()>0:
        ret = ret.astype(bool)
        bad_angs       = view_ang_names[~ret]
        # src_files      = view_ang_names[ret][abs(np.arange(13)[~ret][...,None] - np.arange(13)[ret]).argmin(axis=1)]
        for i in range(len(bad_angs)):
            copyfile(mean_view_angle_name, bad_angs[i])
    else:
        raise LookupError('failed to reconstract angles...')
    return sun_ang_name, view_ang_names, mean_view_angle_name, toa_refs, cloud_name

if __name__ == "__main__":
    metafile = '/Users/fengyin/Downloads/S2B_MSIL1C_20220508T105619_N0400_R094_T30UYD_20220508T114205.SAFE/GRANULE/L1C_T30UYD_A027004_20220508T105909/MTD_TL.xml'
    resample_s2_angles_N0400(metafile)
    