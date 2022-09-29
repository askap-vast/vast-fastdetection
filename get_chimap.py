import numpy as np

from vastfast.cube import Filter

# from memory_profiler import profile

FILE_0 = "tests/fixtures/sigcube.npy"
# FILE_1 = "tests/fixtures/rmscube.npy"

sigcube = np.load(FILE_0)
# rmscube = np.load(FILE_1)



   


 # get the matched filter in time axis
f = Filter(sigcube)
##f = Filter(cube.oricube)

ktype = "chisquare"
f.fmap(ktype, width=1)
chimap = f.map
np.save("chimap", chimap)




        