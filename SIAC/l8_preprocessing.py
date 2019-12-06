#/usr/bin/env python
import os
from glob import glob
from SIAC.l8_angle import do_l8_angle

def l8_pre_processing(l8_dir):
    metafiles = []
    for (dirpath, dirnames, filenames)  in os.walk(l8_dir):
        if len(filenames)>0:
            temp = [dirpath + '/' + i for i in filenames]
            for j in temp:
                if 'mtl.' in j.lower():
                    metafiles.append(j)
    l8_tiles = []
    for metafile in metafiles:
        ret = do_l8_angle(metafile)
        l8_tiles.append(ret)
    return l8_tiles
