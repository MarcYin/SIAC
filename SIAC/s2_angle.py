#!/usr/bin/env python 
import os
import sys
import ogr
import gdal
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
    return vaa, vza

def get_angle(view_ang_name_gml, vaa, vza, band_dict):
    band_name, view_ang_name, gml = view_ang_name_gml
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
    layer = g.GetLayer()
    foot1 = None                            
    foot2 = None                                   
    va1   = None
    vz1   = None         
    va2   = None                                                     
    vz2   = None     
    try:
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
        if os.path.exists(view_ang_name):                   
            os.remove(view_ang_name)                        
        dst_ds = gdal.GetDriverByName('GTiff').Create(view_ang_name, x_size, y_size, 2, gdal.GDT_Int16, options=["TILED=YES", "COMPRESS=DEFLATE"])
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
        vas[vas>180] = vas[vas>180] - 360         
        dst_ds.GetRasterBand(1).WriteArray((vas * 100).astype(int))            
        dst_ds.GetRasterBand(2).WriteArray((vzs * 100).astype(int))            
        dst_ds.GetRasterBand(1).SetNoDataValue(-32767)       
        dst_ds.GetRasterBand(2).SetNoDataValue(-32767)       
        dst_ds.FlushCache()                                  
        dst_ds = None  
        g1 = None   
        return True  
    except:
        return False

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

def resample_s2_angles(metafile):
    #check the available rams and decide cores can be used
    #start multiprocessing
    bands    = ['B01', 'B02', 'B03','B04','B05' ,'B06', 'B07', 'B08','B8A', 'B09', 'B10', 'B11', 'B12'] #all bands
    band_ram = 5e9
    av_ram = psutil.virtual_memory().available 
    procs = np.min([int(av_ram / band_ram), psutil.cpu_count(), len(bands)])
    if procs < 1:
        raise MemoryError('At least 500MB ram is needed.')
    s2_file_dir = os.path.dirname(metafile)
    if ('MTD' in metafile) & ('TL' in metafile) & ('xml' in metafile):
        if not os.path.exists(s2_file_dir + '/ANG_DATA/'):
            os.mkdir(s2_file_dir + '/ANG_DATA/')
        gmls = glob(s2_file_dir + '/QI_DATA/*MSK_DETFOO*.gml')
        sun_ang_name = s2_file_dir + '/ANG_DATA/' + 'SAA_SZA.tif'
        view_ang_names = [s2_file_dir + '/ANG_DATA/' + 'VAA_VZA_%s.tif'%band for band in bands]
        toa_refs = glob(s2_file_dir + '/IMG_DATA/*B??.jp2')
    elif 'metadata.xml' in metafile:
        if not os.path.exists(s2_file_dir + '/angles/'):
            os.mkdir(s2_file_dir + '/angles/')
        gmls = glob(s2_file_dir + '/qi/*MSK_DETFOO*.gml')
        sun_ang_name = s2_file_dir + '/angles/' + 'SAA_SZA.tif'
        view_ang_names = [s2_file_dir + '/angles/' + 'VAA_VZA_%s.tif'%band for band in bands]
        toa_refs = glob(s2_file_dir + '/*B??.jp2')
    else:                                               
        raise IOError('Invalid metafile please use the default AWS or SCIHUB format.')
    gmls = sorted(gmls, key = lambda gml: bands.index('B' + gml.split('B')[-1][:2]))
    toa_refs = sorted(toa_refs, key = lambda toa_ref: bands.index('B' + toa_ref.split('B')[-1][:2]))
    cloud_name = s2_file_dir+'/cloud.tif'
    example_file = toa_refs[1]
    vaa, vza = parse_xml(metafile, example_file, sun_ang_name)
    view_ang_name_gmls = zip(toa_refs, view_ang_names, gmls)
    band_dict = dict(zip(bands, range(13)))
    par = partial(get_angle, vaa=vaa, vza=vza, band_dict=band_dict)
    p = Pool(procs)
    ret = p.map(par,  view_ang_name_gmls)
    #ret  =list( map(par,  view_ang_name_gmls))
    p.close()
    p.join()
    ret = np.array(ret)
    view_ang_names = np.array(view_ang_names)
    if ret.sum()>0:
        bad_angs       = view_ang_names[~ret]
        src_files      = view_ang_names[ret][abs(np.arange(13)[~ret][...,None] - np.arange(13)[ret]).argmin(axis=1)]
        for i in range(len(bad_angs)):
            copyfile(src_files[i], bad_angs[i])
    else:
        raise LookupError('failed to reconstract angles...')
    return sun_ang_name, view_ang_names, toa_refs, cloud_name
