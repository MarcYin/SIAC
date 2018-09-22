import os, sys
import requests
import numpy as np

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../../')
user_pass = os.environ['Earthdata_user'], os.environ['Earthdata_pass']
 
with open(myPath.replace('tests', '') + 'data/.earthdata_auth', 'wb') as f:
    for i in user_pass:
        f.write((i+'\n').encode())

from SIAC import SIAC_S2

with open(myPath + '/s2_flists.txt', 'rb') as f:
    urls =  [(i.split('\n')[0]).decode() for i in f.readlines()]

import hashlib
def file_as_bytes(file):
    with file:
        return file.read()

for url in urls:
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
            with open(filename, "w") as f:
                f.write(req.content)
    else:
        pass

flists_md5 = np.loadtxt(myPath + '/md5Checksum', dtype='str')

for f_md5 in flists_md5:
    filename, md5 = f_md5
    if os.path.exists(filename):
        assert (hashlib.md5(file_as_bytes(open(filename, 'rb'))).hexdigest() == md5, 'file size is not right')
    else:
        raise IOError
s2_file_dir = filename.split('/')[0]
SIAC_S2(s2_file_dir)


with open(myPath.replace('tests', '') + 'data/.earthdata_auth', 'wb') as f:
    f.write(('').encode())  
