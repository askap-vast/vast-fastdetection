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


### exclude Gaussian mao 
## generating how many different types of map
ktypelist = ['chisquare', 'peak', 'std'] 

## how many different types of map to select candidates 
# maplist = ['chisquare']
maplist = ['chisquare', 'peak']

### include Gaussian map
#ktypelist = ['chisquare', 'peak', 'std', 'gaussian']
#maplist = ['chisquare', 'peak', 'gaussian']



## ====================================
## get the imagelist with correct order
imagelist = []
for size in ["?", "??", "???", "????"]:
    tmp = glob.glob(os.path.join(folder, f'*{beam}_{size}.fits'))
    # tmp = glob.glob(folder + f'image_{size}.fits') # for FRB field 
    tmp.sort()
    imagelist += tmp

# imagelist = glob.glob(os.path.join(folder, f"{beam}*cutout.fits"))
# imagelist.sort()

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
num = cube.sigcube.shape[0]


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
chisquare_map = os.path.join(outdir, name+'_chisquare.fits')
peak_map = os.path.join(outdir, name+'_peak.fits')
std_map = os.path.join(outdir, name+'_std.fits')



## =============== select candidates =================
for maptype in maplist:
    
    if 'gaussian' not in maplist:
        c = Candidates(chisquare_map, peak_map, std_map, num=num)
        
    else:
        ## include Gaussian map during candidates selection 
        gaussian_map = os.path.join(outdir, name+'_gaussian.fits')
        c = Candidates(chisquare_map, peak_map, std_map, gaussian_map=gaussian_map, 
                       num=num)
        
    if maptype == "chisquare":
        sigma = 5
    else:
        sigma = 6


    # find local maximum
    logger.info("Find local maximum....")
    c.local_max(min_distance=30, sigma=sigma, data=maptype)
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



    
# =====================
# combine those three cand list to one
logger.info("=========Combine catalogue==========")


namelist = ['{}/{}_{}_cand.csv'.format(outdir, name, maptype) for maptype in maplist ]

plot.combine_csv(namelist, tablename="{}/{}_final".format(outdir, name), 
             savevot=True)



# =====================
# plot final candidates 
logger.info("========= Plotting =============")

final_csv = "{}/{}_final.csv".format(outdir, name)

if os.path.exists(final_csv):
    p = Products(final_csv, limit=50)
    p.generate_slices(imagelist=imagelist, 
                      savename='{}/{}_slices'.format(outdir, name))
    p.generate_cutout(fitsname=deepimage, 
                      savename='{}/{}_deepcutout'.format(outdir, name))
    p.generate_lightcurve(imagelist=imagelist, 
                          deepname=deepimage, 
                          savename='{}/{}_lightcurve'.format(outdir, name))
else:
    logger.info("Final csv {} not exists. ".format(final_csv))


logger.info("====== Finished. =====")













