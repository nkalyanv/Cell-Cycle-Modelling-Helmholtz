import scipy.io as sio
import matlab.engine
import numpy as np

def extend(arr):
    for i in range(len(arr)):
        c = arr[i]
        arr[i] = c + c[:, -1]

def _check_keys( dict):
#checks if entries in dictionary are mat-objects. If yes
#todict is called to change them to nested dictionaries
    for key in dict:
        if isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


def loadmat(filename):
    """
    this function should be called instead of direct scipy.io .loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def KS_bindata_mean_20140916(a, b, val, nargs_out = 3):
    eng = matlab.engine.start_matlab()
    a = matlab.double(a.tolist())
    b = matlab.double(b.tolist())
    [mean, error, binsmean] = eng.KS_bindata_mean_20140916(a, b, val, nargout = nargs_out)
    [mean] = mean
    [error] = error
    [binsmean] = binsmean
    error = np.array(error)
    binsmean = np.array(binsmean)
    mean = np.array(mean)
    return mean, error, binsmean


def KS_bindata_mean_20140506(a, b, val, nargs_out = 3):
    eng = matlab.engine.start_matlab()
    a = matlab.double(a.tolist())
    b = matlab.double(b.tolist())
    [mean, error, binsmean] = eng.KS_bindata_mean_20140506(a, b, val, nargout = nargs_out)
    [mean] = mean
    [error] = error
    [binsmean] = binsmean
    error = np.array(error)
    binsmean = np.array(binsmean)
    mean = np.array(mean)
    return mean, error, binsmean

def fillForward(arr):
    for i in range(arr.shape[0]):
        val = 0
        for j in range(arr.shape[1]):
            if(arr[i][j] != 0):
                val = arr[i][j]
            else:
                arr[i][j] = val