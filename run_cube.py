#!/usr/bin/env python3

import sys
import glob
import matplotlib.pyplot as plt

from vastfast.cube import Cube, Filter



import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
sh = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
sh.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(sh)


# the folder location with processed short images
folder = sys.argv[-3]
name = sys.argv[-1]
beam = sys.argv[-2]
# save folder 
outdir = "/o9000/ASKAP/VAST/fast_test/detection_results"


# get the imagelist with correct order
imagelist = []
for size in ["?", "??", "???", "????"]:
    tmp = glob.glob(folder + f'{beam}_{size}.fits')
    tmp.sort()
    imagelist += tmp

#imagelist = glob.glob(folder + f'{beam}_?.fits') + glob.glob(folder + f'{beam}_??.fits') + glob.glob(folder + f'{beam}_???.fits')
#imagelist = glob.glob(folder + 'image_?.fits') + glob.glob(folder + 'image_??.fits') + glob.glob(folder + 'image_???.fits') + glob.glob(folder + 'image_????.fits')
#imagelist.sort()
logger.info("Loading foler {}".format(folder))
logger.info("Processing {} images...".format(len(imagelist)))
logger.info(imagelist)

# get the psf list with correct order
#psflist = glob.glob(folder + 'beam??_?.psf.fits') + glob.glob(folder + 'beam??_??.psf.fits')
#logger.info("Processing {} psf...".format(len(psflist)))
#logger.info(psflist)

# get the significance cube (do smooth with each 2d image)
ktype = 'gaussian'
logger.info("============")
logger.info("Starting to build the cube...")

cube = Cube(imagelist)
cube.icube(ktype, 19, 19)

logger.info(cube.sigcube.shape)
logger.info("Finish to create the cube.")
logger.info("============")

# get the matched filter in time axis
f = Filter(cube.sigcube)

logger.info("===== Matched Filter =====")
ktype = "chisquare"
logger.info("Kernel match filter '{}'...".format(ktype))
f.fmap(ktype, width=1)
logger.info("Kernel match Done")

f.tofits(fitsname="{}/{}_{}.fits".format(outdir, name, ktype))
logger.info("Save the results to {}_{}.fits".format(name, ktype))


logger.info("Finding local maximum...")
f.local_max(imagename=imagelist[0], sigma=5, min_distance=120)
logger.info("Finding local maximum: Done. ")


#logger.info("Save the smooth cube...")
#import numpy as np
#np.save("../test_products/{}_{}_smocube.npy".format(name, ktype), f.smocube)



# get the position
for i in range(len(f.coord)):
    print(i, f.coord[i].ra.hms, f.coord[i].dec.dms)

logger.info("====== Finished. =====")













