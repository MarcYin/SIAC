def test_import_SIAC_S2():
    try:
        from SIAC import SIAC_S2
        ret = True
    except:
        ret = False
    assert ret

def test_import_SIAC_L8():      
    try:      
        from SIAC import SIAC_L8
        ret = True 
    except:   
        ret = False
    assert ret
if __name__ == '__main__':
    test_import_SIAC_S2()
    test_import_SIAC_L8()
