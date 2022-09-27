import numpy as np
# from memory_profiler import profile
import dask.array as da
from dask.distributed import Client

FILE_0 = "tests/fixtures/sigcube.npy"
FILE_1 = "tests/fixtures/rmscube.npy"

sigcube = np.load(FILE_0)
rmscube = np.load(FILE_1)

def spliterate(buf, chunk):
    for start in range(0, buf.shape[0], chunk):
        yield buf[start:start + chunk,:,:]
        print("sp: ", start)


# @profile
def _chimap(sig, rms):
    """Chi-square map
    """
    # local_rms = np.std(self.sigcube, axis=(1, 2))
    # return np.apply_along_axis(lambda m: self._chisquare(m, local_rms), axis=0, arr=self.sigcube)
    
    # freedom
    nu = sig.shape[0] - 1
    # mean, rms
    # mean = np.nanmean(sigcube, axis=0)
    # rms = np.nanstd(self.sigcube, axis=(1, 2))
    # rms = self.cube_local_rms()
    
    # for each data point
    # data = (self.sigcube - mean) / rms
    # data = (sigcube - mean) / rmscube
    # res = np.nanmean(sigcube, axis=0)
    # sigcube = np.load(FILE_0)
    # sigcube = sigcube - np.nanmean(sigcube, axis=0)
    # sigcube =  (sigcube - np.nanmean(sigcube, axis=0)) / rmscube
    # print("res: ", res.shape, res.dtype)
    tmp = np.zeros((3000,3000), dtype=np.float32)
    tmp[:] = np.nan
    
    res = np.sum(np.square((sig - np.nanmean(sig, axis=0)) / rms),axis=0)/nu
    
    # np.save("chimap", res)
    
    return res
    


if __name__ == "__main__":
    # _chimap()
    client = Client(n_workers=5)
    da_sig = da.from_array(sigcube, chunks=(40,600,600))
    da_rms = da.from_array(rmscube, chunks=(40,600,600))
    res = _chimap(da_sig, da_rms)
    res.compute()


        