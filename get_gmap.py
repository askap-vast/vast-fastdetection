import numpy as np
import time

from astropy.convolution import convolve, Gaussian1DKernel

import dask.array as da
import dask

import numpy as np
import os

from memory_profiler import profile
import gc

FILE_0 = "tests/fixtures/sigcube.npy"
# FILE_0 = "tests/fixtures/small_sigcube.npy"


sigcube = np.load(FILE_0)
sigcube_t = sigcube.transpose(1,2,0).copy(order="C")
del sigcube
gc.collect()

# da_sigcube = da.from_array(sigcube)
da_sigcube_t = da.from_array(sigcube_t, chunks=(100,100,40))

# blocks = da_sigcube_t.to_delayed().ravel()
# print("blocks: ", blocks)

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


def conv_1d(arr, kernel):
    # kernel = Gaussian1DKernel(stddev=stddev)
    arr_res = convolve(arr, kernel)
    # res = np.nanmax(arr_res) - np.nanmean(arr_res)
    # print("called: ", os.getpid())
    return arr_res

def get_smocube_2(sigcube, stddev):
    kernel = Gaussian1DKernel(stddev=stddev)
    smocube =  np.apply_along_axis(conv_1d, 
                                   axis=2, arr=sigcube, kernel=kernel)
    # np.save("smocube", smocube)
    # np.save("small_smocube", smocube)
    return smocube


# @profile
def run_block(arr1, block_info=None):
    print(block_info)
    kernel = Gaussian1DKernel(stddev=4)
    # res = np.zeros(arr1.shape)
    # res[:] = np.nan
    # print("res shape: ", res.shape)
    # for i in range(arr1.shape[0]):
    #     for j in range(arr1.shape[1]):
    #         res[i,j,:] = convolve(arr1[i,j,:],arr2)
    res = np.apply_along_axis(conv_1d, 
                                   axis=2, arr=arr1, kernel=kernel)
    print("process: ", os.getegid())
    # res = da.nanmax(res_0, axis=2) - np.nanmean(res_0, axis=2)
    # time.sleep(20)
    return res

# @profile
def get_smocube_3(sigcube):
    tt = da.map_blocks(run_block, sigcube, chunks=(100,100,40))
    res = da.nanmax(tt, axis=2) - da.nanmean(tt, axis=2)
    print("res chunksize: ", res.chunksize)
    return res

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
    da_smocube = get_smocube_3(da_sigcube_t)
    print("da_smocube shape: ", da_smocube.shape)
    print("chunksize: ", da_smocube.chunksize)
    smocube = da_smocube.compute(scheduler='processes', num_workers=4)
    # smocube = da_smocube.compute(num_workers=1)
    # res = smocube.transpose(2,0,1).copy(order="C")
    np.save("smocube", smocube)
    # np.save("small_smocube", smocube)
    mid = time.time()
    print("generate smocube: ", mid - start)

    
    
