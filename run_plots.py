#!/usr/bin/env python

import os
import sys
import glob

from vastfast.cube import Cube, Filter
from vastfast import plot
from vastfast.plot import Candidates, Products


import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
sh = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
sh.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(sh)


## ========= input ================
## deep image
deepimage = sys.argv[-6]
## catalgoue csv file contain interested source info 
catalogue = sys.argv[-5]
## the folder location with processed short images
folder = sys.argv[-4]
## beam number - format of "beam00"
beam = sys.argv[-3]
## output folder 
outdir = sys.argv[-2]
## output file name 
name = sys.argv[-1]


## ====================================
## get the imagelist with correct order
imagelist = []
for size in ["?", "??", "???", "????"]:
    tmp = glob.glob(os.path.join(folder, f'*{beam}_{size}.fits'))
    # tmp = glob.glob(folder + f'image_{size}.fits') # for FRB field 
    tmp.sort()
    imagelist += tmp

imagelist = imagelist[:-1]

logger.info("Loading foler {}".format(folder))
logger.info("Processing {} images...".format(len(imagelist)))
logger.info(imagelist)


# =====================
# plot final candidates 
logger.info("========= Plotting =============")

# final_csv = "{}/{}_final.csv".format(outdir, name)
final_csv = catalogue

p = Products(final_csv)
p.generate_slices(imagelist=imagelist, 
                  savename='{}/{}_slices'.format(outdir, name))
p.generate_cutout(fitsname=deepimage, 
                  savename='{}/{}_deepcutout'.format(outdir, name))
p.generate_lightcurve(imagelist=imagelist, 
                      deepname=deepimage, 
                      savename='{}/{}_lightcurve'.format(outdir, name))





logger.info("====== Finished. =====")













