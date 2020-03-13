from pyscf import dft

def is_dft(xc_code):
    try:
        dft.libxc.xc_type(xc_code)
    except:
        return False
    else:
        return True
