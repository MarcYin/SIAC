import os
import requests
import logging
from SIAC.create_logger import create_logger
'''
logger = logging.getLogger('SIAC')
logger.setLevel(logging.INFO)
if not logger.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
'''
logger = create_logger()

def downloader(fname, url_root, file_dir):
    logger.propagate = False
    new_url = url_root + fname
    new_req = requests.get(new_url)
    if new_req.ok:   
        logger.info('downloading %s and save it at %s' % (fname, os.path.join(file_dir, fname)))
        with open(os.path.join(file_dir, fname), 'wb') as fp:
            #for chunk in new_req.iter_content(chunk_size=1024):
            #    if chunk:
            fp.write(new_req.content)
    else:            
        logger.error(new_req.content)
