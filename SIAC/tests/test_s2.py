import os, sys
import requests
import zipfile
import numpy as np
from functools import partial
from multiprocessing import Pool

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../../')
user_pass = os.environ['Earthdata_user'], os.environ['Earthdata_pass']
 
with open(myPath.replace('tests', '') + 'data/.earthdata_auth', 'wb') as f:
    for i in user_pass:
        f.write((i+'\n').encode())

def test_s2():
    from SIAC import SIAC_S2
    from SIAC.downloaders import downloader
    with open(myPath + '/s2_flists.txt', 'rb') as f:
        urls =  [i.decode().split('\n')[0] for i in f.readlines()]

    for url in urls:
        filename = '/'.join(url.split('/')[8:])
        if not os.path.exists(filename):
            if not os.path.exists(os.path.dirname(filename)):           
                try:               
                    os.makedirs(os.path.dirname(filename))
                except OSError as exc: # Guard against race condition               
                    if exc.errno != errno.EEXIST:          
                        raise              
            downloader(filename, '/'.join(url.split('/')[:8]) + '/', './')
        else:
            pass

    with open(myPath + '/MCD43.txt', 'rb') as f:             
        MCD43 =  [i.decode().split('\n')[0] for i in f.readlines()]
    if not os.path.exists(os.path.expanduser("~") + '/MCD43/'):
        os.makedirs(os.path.expanduser("~") + '/MCD43/')
    #downloader('MCD43.zip', url_root = 'http://www2.geog.ucl.ac.uk/~ucfafyi/mcd43/', file_dir = './')
   
    #par = partial(downloader, url_root = 'http://www2.geog.ucl.ac.uk/~ucfafyi/mcd43/MCD43/', file_dir = os.path.expanduser("~") + '/MCD43/')
    #p = Pool(4)
    #p.map(par, MCD43)
    #p.close()
    #p.join()
    #with zipfile.ZipFile("MCD43.zip","r") as zip_ref:
    #    zip_ref.extractall(os.path.expanduser("~"))
    #os.remove("MCD43.zip")

    s2_file_dir = filename.split('/')[0]
    SIAC_S2(s2_file_dir, aoi = str(myPath + '/aoi.geojson'))

    with open(myPath.replace('tests', '') + 'data/.earthdata_auth', 'wb') as f:
        f.write(('').encode())  
    assert True

if __name__ == '__main__':
    test_s2()
