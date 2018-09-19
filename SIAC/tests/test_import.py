import sys, os
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../../')
user_pass = 'user', 'pass'

with open(myPath.replace('tests', '') + 'data/.earthdata_auth', 'wb') as f:
    for i in user_pass:
        f.write((i+'\n').encode())

from SIAC import SIAC_S2
from SIAC import SIAC_L8


with open(myPath.replace('tests', '') + 'data/.earthdata_auth', 'wb') as f:
    f.write(('').encode())
