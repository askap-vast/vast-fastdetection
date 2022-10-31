from glob import glob
import os
import re

path = os.getcwd()

file_path = path + "/../test_data/SB9596_beam14_casa/beam*"

files = glob(file_path)

for ff in files:
    path_dir = ff.split("/")
    fname = path_dir[-1]
    prefix = "/".join(path_dir[:-1])
    print(fname)
    print(prefix)

    d = fname.split(".")
    nn = d[0].split("_")
    print(nn[-1])

    out ="{}_{:04}.{}".format(nn[0],int(nn[-1]),d[-1])
    new = "/".join([prefix, out])
    print(new)
    # os.rename(ff, new)



