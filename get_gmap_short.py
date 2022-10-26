import numpy as np

from vastfast.cube import Filter, _gmap

import multiprocessing as mp

# ctx = mp.get_context("fork")

# from memory_profiler import profile

FILE_0 = "tests/fixtures/sigcube.npy"
# FILE_0 = "small_sigcube.npy"
# FILE_1 = "tests/fixtures/rmscube.npy"

sigcube = np.load(FILE_0)
# rmscube = np.load(FILE_1)



   

if __name__ == "__main__":
    # get the matched filter in time axis
    # f = Filter(sigcube)
    # ##f = Filter(cube.oricube)
    # print("before gaussian map")
    # ktype = "gaussian"
    # f.fmap(ktype, width=4)
    # gmap = f.map
    # np.save("gaussianmap", gmap)

    _gmap(sigcube)
    


        