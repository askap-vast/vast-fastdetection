#!/usr/bin/env python

import sys
import glob
import matplotlib.pyplot as plt

from vastfast.cube import Cube, Filter
from vastfast import plot
from vastfast.plot import Candidates



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


logger.info("======== Select candidates ==========")

# read fits
chisq_map = "{}/{}_{}.fits".format(outdir, name, 'chisquare')
peak_map = "{}/{}_{}.fits".format(outdir, name, 'peak')
std_map = "{}/{}_{}.fits".format(outdir, name, 'std')
gaussian_map = "{}/{}_{}.fits".format(outdir, name, 'gaussian')

c = Candidates(chisq_map, peak_map, std_map)
#c = Candidates(chisq_map, peak_map, std_map, gaussian_map=gaussian_map)


# find local maximum
logger.info("Find local maximum....")
#c.local_max(min_distance=30, sigma=5, data='peak')
c.local_max(min_distance=30, sigma=5)
logger.info("Find local maximum done. ")

# plot a map
#c.plot_fits(fitsname=chisq_map, imagename="{}/{}_map1".format(outdir, name))

# read deep information to select candidates
#catalogue = glob.glob(folder + "*_{}*comp.vot".format(beam))
catalogue = sys.argv[-4]

logger.info("Deep image catalogue {}".format(catalogue))
c.select_candidates(deepcatalogue=catalogue)

# save the table
c.save_csvtable(tablename="{}/{}_cand".format(outdir, name), savevot=True)

# plot a final map
#c.plot_fits(fitsname=chisq_map, imagename="{}/{}_map2".format(outdir, name))

# plot!
for i, candname in enumerate(c.cand_name):
    logger.info("Plot slices {}/{}: {}".format(i, len(c.cand_name), candname))
    plot.plot_slices(src_name=candname, 
                     imagelist=imagelist, 
                     name="{}/{}_{}".format(outdir, name, candname))



logger.info("====== Finished. =====")













