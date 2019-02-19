import os
import logging
file_path = os.path.dirname(os.path.realpath(__file__))                                                                                                                                                                 
version_file = file_path + '/VERSION'
with open(version_file, 'rb') as f:
    version = f.read().decode().strip()
def create_logger(fname = None):
    logger = logging.getLogger('SIAC-V%s'%version)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    if not logger.handlers:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    if fname is not None:
        fh = logging.FileHandler(fname)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)    
        logger.addHandler(fh)
    return logger
