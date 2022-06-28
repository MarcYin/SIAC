#/usr/bin/env python
import os
from glob import glob
from SIAC.l8_angle import do_l8_angle, do_l8_angle_collection2

def l8_pre_processing(l8_dir):
    metafiles = []
    for (dirpath, dirnames, filenames)  in os.walk(l8_dir):
        if len(filenames)>0:
            temp = [dirpath + '/' + i for i in filenames]
            for j in temp:
                if 'mtl.' in j.lower():
                    metafiles.append(os.path.realpath(j))
    l8_tiles = []
    for metafile in metafiles:
        ret = do_l8_angle(metafile)
        l8_tiles.append(ret)
    return l8_tiles

def l8_pre_processing(l8_dir):
    metafiles = []
    for (dirpath, dirnames, filenames)  in os.walk(l8_dir):
        if len(filenames)>0:
            temp = [dirpath + '/' + i for i in filenames]
            for j in temp:
                if 'mtl.txt' in j.lower():
                    metafiles.append(os.path.realpath(j))
    l8_tiles = []
    for metafile in metafiles:
        base = os.path.dirname(metafile)
        collection = int(base.split('_')[-2])
        if collection == 2:
            ret = do_l8_angle_collection2(metafile)
        elif collection == 1:
            ret = do_l8_angle(metafile)
        else:
            raise IOError('Only collection <=2 is supported')
        l8_tiles.append(ret)
    return l8_tiles

