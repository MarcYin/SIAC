import os, sys
import requests
import numpy as np

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../../')
user_pass = os.environ['Earthdata_user'], os.environ['Earthdata_pass']
 
with open(myPath.replace('tests', '') + 'data/.earthdata_auth', 'wb') as f:
    for i in user_pass:
        f.write((i+'\n').encode())

def test_s2():
    from SIAC import SIAC_S2
    with open(myPath + '/s2_flists.txt', 'rb') as f:
        urls =  [i.decode().split('\n')[0] for i in f.readlines()]

    import hashlib
    def file_as_bytes(file):
        with file:
            return file.read()

    for url in urls:
        url = url
        filename = '/'.join(url.split('/')[8:])
        if not os.path.exists(filename):
            req = requests.get(url)
            if req.ok:
                print('downloading %s' % filename)
                if not os.path.exists(os.path.dirname(filename)):
                    try:
                        os.makedirs(os.path.dirname(filename))
                    except OSError as exc: # Guard against race condition
                        if exc.errno != errno.EEXIST:
                            raise
                with open(filename, "wb") as f:
                    f.write(req.content)
        else:
            pass

    s2_file_dir = filename.split('/')[0]
    SIAC_S2(s2_file_dir, aoi = str(myPath + '/aoi.geojson'))

    with open(myPath.replace('tests', '') + 'data/.earthdata_auth', 'wb') as f:
        f.write(('').encode())  
    assert True
