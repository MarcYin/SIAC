#/usr/bin/env python 
import os
import gdal
import warnings
import subprocess
import numpy as np
from glob import glob
from functools import partial
from multiprocessing import Pool

warnings.filterwarnings("ignore")

file_path = os.path.dirname(os.path.realpath(__file__))
angle_exe = file_path + '/util/l8_angles'

def post_view_angle(view_angle, out_name, delete_origin = False):
    g = gdal.Open(view_angle)
    vaa, vza = g.ReadAsArray()/100.
    vaa[vaa==-327.67] = np.nan
    vza[vza==-327.67] = np.nan
    x, y = np.where(vza==np.nanmin(vza))
    p = np.poly1d(np.polyfit(y,x,1))
    x1, y1 = np.where(~np.isnan(vza))
    dist = abs(p(y1)-x1)/(np.sqrt(1+p.c[0]**2))
    vza[x1, y1] = np.nanmin(vza) + dist/(dist.max()-dist.min()) * (np.nanmax(vza)-np.nanmin(vza))
    hist, bin_edges = np.histogram(vaa[vaa<0], bins=1000)
    rvaa = bin_edges[np.argmax(hist)] + 180
    hist, bin_edges = np.histogram(vaa[vaa>0], bins=1000)
    lvaa = bin_edges[np.argmax(hist)]
    mvaa = np.mean([rvaa, lvaa])
    vaa[vaa<0] = mvaa - 360
    vaa[vaa>0] = mvaa
    vaa[np.isnan(vaa)] = -327.67
    vza[np.isnan(vza)] = -327.67
    if os.path.exists(out_name):                   
        os.remove(out_name)                        
    x_size, y_size = g.RasterXSize, g.RasterYSize
    dst_ds = gdal.GetDriverByName('GTiff').Create(out_name, x_size, y_size, 2, gdal.GDT_Int16, options=["TILED=YES", "COMPRESS=DEFLATE"])
    dst_ds.SetGeoTransform(g.GetGeoTransform())         
    dst_ds.SetProjection(g.GetProjection())
    dst_ds.GetRasterBand(1).WriteArray((vaa * 100).astype(int))            
    dst_ds.GetRasterBand(2).WriteArray((vza * 100).astype(int))            
    dst_ds.GetRasterBand(1).SetNoDataValue(-32767)       
    dst_ds.GetRasterBand(2).SetNoDataValue(-32767)       
    dst_ds.FlushCache()                                  
    dst_ds = None  
    g = None 
    if delete_origin:
        os.remove(view_angle)    
        os.remove(view_angle + '.hdr')

def usgs_angle(band_angType, ang_file):
    band, ang_type = band_angType
    FNULL = open(os.devnull, 'w')
    subprocess.call([angle_exe, ang_file, ang_type, '1', '-f', '-32767', '-b', str(band)], stdout=FNULL, stderr=subprocess.STDOUT)
    header     = '_'.join(ang_file.split('/')[-1].split('_')[:-1])
    view_angle = header  + '_sensor_B%02d.img'%band
    out_name   = header  + '_VAA_VZA_B%02d.TIF'%band
    post_view_angle(view_angle, out_name, delete_origin = True)
    return out_name


def do_l8_angle(metafile):
    l8_file_dir = os.path.dirname(metafile)
    if not os.path.exists(l8_file_dir + '/angles'):
        os.mkdir(l8_file_dir + '/angles')
    header    = '_'.join(metafile.split('/')[-1].split('_')[:-1])
    bs        = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'BQA']
    bands = [l8_file_dir + '/' + header + '_%s.TIF'%i for i in bs]
    toa_bands = bands[:-1]
    qa_band  = bands[-1]
    ang_file = l8_file_dir + '/' + header + '_ANG.txt'
    bands    = np.arange(1, 8)
    ang_types = ['BOTH', ] + ['SATELLITE',] * 7
    band_angTypes = zip(bands, ang_types)
    par = partial(usgs_angle, ang_file = ang_file)
    cwd = os.getcwd()
    os.chdir(l8_file_dir + '/angles')
    p = Pool() 
    rets = p.map(par, band_angTypes) 
    p.close()
    p.join()  
    sun_angle  = header  + '_solar_B01.img'
    nsun_angle = header  + '_SAA_SZA.TIF'
    gdal.Translate(nsun_angle, sun_angle, creationOptions = ['COMPRESS=DEFLATE', 'TILED=YES'], noData='-32767').FlushCache()
    rets.append(nsun_angle)
    os.remove(sun_angle)
    os.remove(sun_angle + '.hdr')
    os.chdir(cwd)
    
    ang_names = [l8_file_dir + '/angles/' + i for i in rets]
    view_ang_names = ang_names[:-1]
    sun_ang_name = ang_names[-1]
    cloud_mask = gdal.Open(qa_band).ReadAsArray()
    cloud_mask = ~((cloud_mask  >= 2720) & ( cloud_mask <= 2732))

    return sun_ang_name, view_ang_names, toa_bands, qa_band, cloud_mask, metafile
