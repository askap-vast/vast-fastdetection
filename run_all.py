#!/usr/bin/env python

import sys
import glob

from vastfast.cube import Cube, Filter
# from vastfast import plot
from vastfast.plot import Candidates


import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
sh = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
sh.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(sh)


## ========= input ================
## deep catalogue
catalogue = sys.argv[-5]
## the folder location with processed short images
folder = sys.argv[-4]
## beam number - format of "beam00"
beam = sys.argv[-3]
## output folder 
outdir = sys.argv[-2]
## output file name 
name = sys.argv[-1]

## generating how many different types of map
ktypelist = ['chisquare', 'peak', 'std'] 
# ktypelist = ['chisquare', 'peak', 'std', 'gaussian']

## how many different types of map to select candidates 
# maplist = ['chisquare']
maplist = ['chisquare', 'peak']
# maplist = ['chisquare', 'peak', 'gaussian']


## ====================================
## get the imagelist with correct order
imagelist = []
for size in ["?", "??", "???", "????"]:
    tmp = glob.glob(folder + f'{beam}_{size}.fits')
    # tmp = glob.glob(folder + f'image_{size}.fits') # for FRB field 
    tmp.sort()
    imagelist += tmp


logger.info("Loading foler {}".format(folder))
logger.info("Processing {} images...".format(len(imagelist)))
logger.info(imagelist)


## ===================================
## get the psf list with correct order
# psflist = glob.glob(folder + 'beam??_?.psf.fits') + glob.glob(folder + 'beam??_??.psf.fits')
# logger.info("Processing {} psf...".format(len(psflist)))
# logger.info(psflist)
## ===================================


logger.info("============")
logger.info("Starting to build the cube...")
cube = Cube(imagelist)

## get the significance cube (do smooth with each 2d image)
# ktype = 'gaussian' # 'psf'
# cube.icube(ktype, 19, 19)

cube.icube()
logger.info(cube.sigcube.shape)
logger.info("Finish to create the cube.")
logger.info("============")

logger.info("Remove bad images...")
cube.remove_bad_images()
logger.info(cube.sigcube.shape)


## ====================================
## get the matched filter in time axis
f = Filter(cube.sigcube)



for ktype in ktypelist:

    logger.info("===== Matched Filter =====")
    logger.info("Kernel match filter '{}'...".format(ktype))
    f.fmap(ktype, width=4)
    logger.info("Kernel match Done")
    
    f.tofits(fitsname="{}/{}_{}.fits".format(outdir, name, ktype), imagename=imagelist[0])
    logger.info("Save the results to {}_{}.fits".format(name, ktype))





logger.info("======== Select candidates ==========")

# read fits
chisquare_map = "{}/{}_{}.fits".format(outdir, name, 'chisquare')
peak_map = "{}/{}_{}.fits".format(outdir, name, 'peak')
std_map = "{}/{}_{}.fits".format(outdir, name, 'std')


# if 'gaussian' in maplist:
#     ## include Gaussian map during candidates selection 
#     gaussian_map = "{}/{}_{}.fits".format(outdir, name, 'gaussian')
    
#     c = Candidates(chisquare_map, peak_map, std_map, gaussian_map=gaussian_map)
# else:
    
#     c = Candidates(chisquare_map, peak_map, std_map)




## =============== select candidates =================
for maptype in maplist:
    
    if maptype != 'gaussian':
        c = Candidates(chisquare_map, peak_map, std_map)
        
    else:
        ## include Gaussian map during candidates selection 
        gaussian_map = "{}/{}_{}.fits".format(outdir, name, 'gaussian')
        c = Candidates(chisquare_map, peak_map, std_map, gaussian_map=gaussian_map)
        

    # find local maximum
    logger.info("Find local maximum....")
    c.local_max(min_distance=30, sigma=5, data=maptype)
    logger.info("Find local maximum done. ")
    
    ## plot a map with all of candidates above the threshold 
    c.plot_fits(fitsname=vars()[maptype+'_map'], 
                imagename="{}/{}_{}_map1".format(outdir, name, maptype))
        
    
    logger.info("Deep image catalogue {}".format(catalogue))
    c.select_candidates(deepcatalogue=catalogue)
    
    # save the table
    c.save_csvtable(tablename="{}/{}_{}_cand".format(outdir, name, maptype), savevot=True)
    
    ## plot a final map with promising candidates 
    c.plot_fits(fitsname=vars()[maptype+'_map'], 
                imagename="{}/{}_{}_map2".format(outdir, name, maptype))



# # plot!
# for i, candname in enumerate(c.cand_name):
#     logger.info("Plot slices {}/{}: {}".format(i, len(c.cand_name), candname))
#     plot.plot_slices(src_name=candname, 
#                      imagelist=imagelist, 
#                      name="{}/{}_{}".format(outdir, name, candname))



logger.info("====== Finished. =====")













