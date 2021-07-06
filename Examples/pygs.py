from numpy.ctypeslib import ndpointer
import ctypes
import numpy as np

def single_thread_wrapper(mwfunc, Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, verbose=False):
    _intp = ndpointer(dtype=ctypes.c_int32, flags='C')
    _doublep = ndpointer(dtype=ctypes.c_double, flags='C')
    mwfunc.argtypes = [_intp, _doublep, _doublep, _doublep, _doublep, _doublep, _doublep]
    mwfunc.restype = ctypes.c_int
    # force the data types to what they need to be
    Lparms = np.array(Lparms, dtype='int32')
    Rparms = np.array(Rparms, dtype='double')
    Parms = np.array(Parms, dtype='double')
    E_arr = np.array(E_arr, dtype='double')
    mu_arr = np.array(mu_arr, dtype='double')
    f_arr = np.array(f_arr, dtype='double')

    # define output
    Nf = Lparms[1]
    RL = np.zeros((Nf, 7), dtype='double')

    void = mwfunc(Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, RL)
    return RL

def get_mw_main(libname, Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, verbose=False):
    # load the library
    try:
        libc_mw = ctypes.CDLL(libname)
        if verbose:
            print("Successfully loaded ", libc_mw)
        mwfunc = libc_mw.pyGET_MW_main
        RL = single_thread_wrapper(mwfunc, Lparms, Rparms, Parms, E_arr, mu_arr, f_arr)
        return RL
    except Exception as e:
        print(e)

def get_mw(libname, Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, verbose=False):
    # load the library
    try:
        libc_mw = ctypes.CDLL(libname)
        if verbose:
            print("Successfully loaded ", libc_mw)
        mwfunc = libc_mw.pyGET_MW
        RL = single_thread_wrapper(mwfunc, Lparms, Rparms, Parms, E_arr, mu_arr, f_arr)
        return RL
    except Exception as e:
        print(e)
