import os
from os.path import expanduser
home = expanduser("~")
file_path = os.path.dirname(os.path.realpath(__file__))

def test_sen2cloud():
    assert os.path.exists(file_path.replace('tests', 'data/sen2cloud_detector.pkl'))

def test_earthdata_auth():
    assert os.path.exists(file_path.replace('tests', 'data/.earthdata_auth'))

def test_spectral_mapping_1():                                                  
    assert os.path.exists(file_path.replace('tests', 'spectral_mapping/TERRA_S2A_spectral_mapping.txt'))

def test_spectral_mapping_2():                                                  
    assert os.path.exists(file_path.replace('tests', 'spectral_mapping/TERRA_S2B_spectral_mapping.txt'))

def test_spectral_mapping_3():                                                  
    assert os.path.exists(file_path.replace('tests', 'spectral_mapping/TERRA_L8_spectral_mapping.txt'))

def test_spectral_mapping_4():                                                  
    assert os.path.exists(file_path.replace('tests', 'spectral_mapping/TERRA_TERRA_spectral_mapping.txt'))

def test_l8_angle():
    assert os.path.exists(file_path.replace('tests', 'util/l8_angles'))

if __name__ == '__main__':
    test_sen2cloud()
    test_earthdata_auth()
    test_spectral_mapping_1()
    test_spectral_mapping_2()
    test_spectral_mapping_3()
    test_spectral_mapping_4()
    test_l8_angle()
