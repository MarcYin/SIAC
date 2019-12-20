import os, sys
import time
import requests
import zipfile
import numpy as np
import subprocess
from functools import partial
from multiprocessing import Pool

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../../')
user_pass = os.environ['Earthdata_user'], os.environ['Earthdata_pass']
 
with open(myPath.replace('tests', '') + 'data/.earthdata_auth', 'wb') as f:
    for i in user_pass:
        f.write((i+'\n').encode())

def downloader(url):
    req = requests.get(url)
    filename = os.path.expanduser("~") + '/MCD43/' + url.split('/')[-1].split('?')[0]
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)

def test_s2():
    from SIAC import SIAC_S2
    #from SIAC.downloaders import downloader
    url = 'gs://gcp-public-data-sentinel-2/tiles/31/U/CU/S2B_MSIL1C_20171123T111349_N0206_R137_T31UCU_20171123T130706.SAFE'
    file_dir = './'
    subprocess.call(['gsutil', '-m', 'cp', '-n', '-r', url, file_dir],stderr=subprocess.STDOUT)
    #with open(myPath + '/s2_flists.txt', 'rb') as f:
    #    urls =  [i.decode().split('\n')[0] for i in f.readlines()]
    #headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.77 Safari/537.36'}
    #for url in urls:
    #    filename = '/'.join(url.split('/')[8:])
    #    if not os.path.exists(filename):
    #        req = requests.get(url)#, headers=headers)
    #        if req.ok:
    #            print('downloading %s' % filename)
    #            if not os.path.exists(os.path.dirname(filename)):
    #                try:
    #                    os.makedirs(os.path.dirname(filename))
    #                except OSError as exc: # Guard against race condition
    #                    if exc.errno != errno.EEXIST:
    #                        raise
    #           with open(filename, "wb") as f:
    #               f.write(req.content)
    #        else:
    #            print(req.content)
    #        time.sleep(1)
    #    else:
    #        pass
    with open(myPath + '/MCD43.txt', 'rb') as f:             
        MCD43 =  [i.decode().split('\n')[0] for i in f.readlines()]
    urls = ['https://zenodo.org/record/3588473/files/%s?download=1'%i for i in MCD43]

    if not os.path.exists(os.path.expanduser("~") + '/MCD43/'):
        os.makedirs(os.path.expanduser("~") + '/MCD43/')
    #downloader('MCD43.zip', url_root = 'http://www2.geog.ucl.ac.uk/~ucfafyi/mcd43/', file_dir = './')
    
    #par = partial(downloader, url_root = 'http://www2.geog.ucl.ac.uk/~ucfafyi/mcd43/MCD43/', file_dir = os.path.expanduser("~") + '/MCD43/')
    p = Pool(4)
    p.map(downloader, urls)
    p.close()
    p.join()
    #with zipfile.ZipFile("MCD43.zip","r") as zip_ref:
    #    zip_ref.extractall(os.path.expanduser("~"))
    #os.remove("MCD43.zip")

    s2_file_dir = './S2B_MSIL1C_20171123T111349_N0206_R137_T31UCU_20171123T130706.SAFE/' #filename.split('/')[0]
    SIAC_S2(s2_file_dir, aoi = str(myPath + '/aoi.geojson'))

    with open(myPath.replace('tests', '') + 'data/.earthdata_auth', 'wb') as f:
        f.write(('').encode())  
    assert True

if __name__ == '__main__':
    test_s2()
