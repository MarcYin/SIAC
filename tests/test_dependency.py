import os
from os.path import expanduser
home = expanduser("~")

def test_sen2cloud():
    sen2cloud_path = file_path = os.path.dirname(os.path.realpath(__file__))
    assert os.path.exists(file_path.replace('tests', 'SIAC/data/sen2cloud_detector.pkl'))

def test_earthdata_auth():
    sen2cloud_path = file_path = os.path.dirname(os.path.realpath(__file__))
    assert os.path.exists(file_path.replace('tests', 'SIAC/data/.earthdata_auth'))

def test_spectral_mapping_1():                                                  
    sen2cloud_path = file_path = os.path.dirname(os.path.realpath(__file__))
    assert os.path.exists(file_path.replace('tests', 'SIAC/spectral_mapping/TERRA_S2A_spectral_mapping.txt'))

def test_spectral_mapping_2():                                                  
    sen2cloud_path = file_path = os.path.dirname(os.path.realpath(__file__))
    assert os.path.exists(file_path.replace('tests', 'SIAC/spectral_mapping/TERRA_S2B_spectral_mapping.txt'))

def test_spectral_mapping_3():                                                  
    sen2cloud_path = file_path = os.path.dirname(os.path.realpath(__file__))
    assert os.path.exists(file_path.replace('tests', 'SIAC/spectral_mapping/TERRA_L8_spectral_mapping.txt'))

def test_spectral_mapping_4():                                                  
    sen2cloud_path = file_path = os.path.dirname(os.path.realpath(__file__))
    assert os.path.exists(file_path.replace('tests', 'SIAC/spectral_mapping/TERRA_TERRA_spectral_mapping.txt'))
def test_l8_angle():
    sen2cloud_path = file_path = os.path.dirname(os.path.realpath(__file__))
    assert os.path.exists(file_path.replace('tests', 'SIAC/util/l8_angles'))


if __name__ == '__main__':
    test_sen2cloud()
    test_earthdata_auth()
    test_spectral_mapping_1()
    test_spectral_mapping_2()
    test_spectral_mapping_3()
    test_spectral_mapping_4()
    test_l8_angle()
