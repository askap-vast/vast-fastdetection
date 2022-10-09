import numpy as np
import time

from astropy.convolution import convolve, Gaussian1DKernel

import dask.array as da
import dask

import numpy as np

FILE_0 = "tests/fixtures/sigcube.npy"
# FILE_0 = "tests/fixtures/small_sigcube.npy"


sigcube = np.load(FILE_0)

da_sigcube = da.from_array(sigcube)

def get_smocube(sigcube):
    kernel = Gaussian1DKernel(stddev=4)
    kernel = da.from_array(kernel.array)
    smocube =  da.apply_along_axis(lambda m: convolve(m, kernel), 
                                   axis=0, arr=sigcube)
    # np.save("smocube", smocube)
    # np.save("small_smocube", smocube)
    return smocube

def get_gmap(smocube):
    res = np.nanmax(smocube, axis=0) - np.nanmean(smocube, axis=0)
    
    # np.save("gaussianmap", res)
    # np.save("small_gmap", res)

@dask.delayed
def conv_1d(arr, stddev):
    kernel = Gaussian1DKernel(stddev=stddev)
    arr_res = convolve(arr, kernel)
    return arr_res

def get_smocube_2(sigcube):
    # kernel = Gaussian1DKernel(stddev=4)
    smocube =  np.apply_along_axis(conv_1d, 
                                   axis=0, arr=sigcube, stddev=4)
    # np.save("smocube", smocube)
    # np.save("small_smocube", smocube)
    return smocube

def run_block(arr1, block_info=None):
    print(block_info)
    arr2 = Gaussian1DKernel(stddev=4)
    return np.apply_along_axis(lambda m: convolve(m, arr2), axis=0, arr=arr1)

def get_smocube_3(sigcube):
    tt = da.map_blocks(run_block, sigcube)
    return tt

if __name__ == "__main__":
    # start = time.time()
    # smocube = get_smocube(sigcube)
    # mid = time.time()
    # get_gmap(smocube)
    # end = time.time()
    # print("generate smocube: ", mid - start)
    # print("generate gmap: ", end - mid)
    # print("total: ", end - start)
    start = time.time()
    da_smocube = get_smocube_3(da_sigcube)
    print("da_smocube shape: ", da_smocube.shape)
    print("chunksize: ", da_smocube.chunksize)
    smocube = da_smocube.compute()
    np.save("smocube", smocube)
    # np.save("small_smocube", smocube)
    mid = time.time()
    print("generate smocube: ", mid - start)

    
    
