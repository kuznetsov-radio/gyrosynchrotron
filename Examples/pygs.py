from numpy.ctypeslib import ndpointer
import ctypes
import numpy as np

def single_thread_wrapper(mwfunc, Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, verbose=False):
    _intp = ndpointer(dtype=ctypes.c_int32, flags='C')
    _doublep = ndpointer(dtype=ctypes.c_double, flags='C')
    mwfunc.argtypes = [_intp, _doublep, _doublep, _doublep, _doublep, _doublep, _doublep]
    mwfunc.restype = ctypes.c_int
    # force the data types
    Lparms = Lparms.astype('int32')
    Rparms = Rparms.astype('double')
    Parms = Parms.astype('double')
    if np.isscalar(E_arr):
        E_arr = float(E_arr)
    else:
        E_arr = E_arr.astype('double')
    if np.isscalar(mu_arr):
        mu_arr = float(mu_arr)
    else:
        mu_arr = mu_arr.astype('double')
    if np.isscalar(f_arr):
        mu_arr = float(mu_arr)
    else:
        f_arr = f_arr.astype('double')

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
