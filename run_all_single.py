#!/usr/bin/env python

import sys
import os
import glob

from vastfast.cube import Cube, Filter
from vastfast import plot
from vastfast.plot import Candidates, Products
from vastfast.input import *
from vastfast.procedures import *
from vastfast.exceptions import *

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
sh = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
sh.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(sh)


## ========= input ================
INPUT_PATH = "../test_data/SB9596_beam14_casa/"
# INPUT_PATH = "../test_data/SB12704/"

# output folder
outdir = "./output"
if not os.path.exists(outdir):
    os.mkdir(outdir)

# output file prefix
out_prefix = "output"

### include Gaussian map
# ktypelist = ['chisquare', 'peak', 'std', 'gaussian']
ktypelist = ['chisquare', 'peak', 'std']
# maplist = ['chisquare', 'peak', 'gaussian']
maplist = ['chisquare', 'peak']



def process_beam(beam):
    input_path = INPUT_PATH

    # input data: short images
    try:
        imagelist = read_shortimage(input_path, beam)
    except NoInputError:
        logger.error("Images are not available; skip the beam{:02}...".format(beam))
        return -1
    
    # input data: catalogue
    try:
        catalogue = read_deepcat(input_path, beam)
    except NoInputError:
        logger.error("Deep catalogue is not available; skip the beam{:02}...".format(beam))
        return -1


    # input data: deep image
    deepimage_file = input_path + "/" + "*beam{:02}*.tt0.fits".format(beam)
    deepimage = glob.glob(deepimage_file)[0]
    
    # output beam dir
    outdir_beam = outdir + "/output_beam{:02}".format(beam)
    # output file prefix
    name = out_prefix + "_beam{:02}".format(beam)
    # output tar file
    tar_name = outdir + "/output_beam{:02}.tar.gz".format(beam)
    
    # logger.info("============")
    # logger.info("Starting to build the cube...")
    # cube = Cube(imagelist)


    # cube.icube()
    # logger.info(cube.sigcube.shape)
    # logger.info("Finish to create the cube.")
    # logger.info("============")

    # logger.info("Remove bad images...")
    # cube.remove_bad_images()
    # logger.info(cube.sigcube.shape)
    try:
        sigcube = get_sigcube(imagelist)
    except Exception:
        logger.error("Fail to generate sugcube; skip the beam{:02}".format(beam))
        return -1

    
    # ====================================
    # get the matched filter in time axis
    f = Filter(sigcube)



    for ktype in ktypelist:

        logger.info("===== Matched Filter =====")
        logger.info("Kernel match filter '{}'...".format(ktype))
        f.fmap(ktype, width=4)
        logger.info("Kernel match Done")
        
        f.tofits(fitsname="{}/{}_{}.fits".format(outdir_beam, name, ktype), imagename=imagelist[0])
        logger.info("Save the results to {}_{}.fits".format(name, ktype))


    # get_map(sigcube, ktypelist, imagelist)


    logger.info("======== Select candidates ==========")

    # read fits
    chisquare_map = "{}/{}_{}.fits".format(outdir_beam, name, 'chisquare')
    peak_map = "{}/{}_{}.fits".format(outdir_beam, name, 'peak')
    std_map = "{}/{}_{}.fits".format(outdir_beam, name, 'std')



    ## =============== select candidates =================
    for maptype in maplist:
        
        if 'gaussian' not in maplist:
            c = Candidates(chisquare_map, peak_map, std_map)
            
        else:
            ## include Gaussian map during candidates selection 
            gaussian_map = "{}/{}_{}.fits".format(outdir_beam, name, 'gaussian')
            c = Candidates(chisquare_map, peak_map, std_map, gaussian_map=gaussian_map)
            

        # find local maximum
        logger.info("Find local maximum....")
        c.local_max(min_distance=30, sigma=5, data=maptype)
        logger.info("Find local maximum done. ")
        
        ## plot a map with all of candidates above the threshold 
        c.plot_fits(fitsname=vars()[maptype+'_map'], 
                imagename="{}/{}_{}_map1".format(outdir_beam, name, maptype))
            
        
        logger.info("Deep image catalogue {}".format(catalogue))
        c.select_candidates(deepcatalogue=catalogue)
        
        # save the table
        c.save_csvtable(tablename="{}/{}_{}_cand".format(outdir_beam, name, maptype), savevot=True)
        
        ## plot a final map with promising candidates 
        c.plot_fits(fitsname=vars()[maptype+'_map'], 
                imagename="{}/{}_{}_map2".format(outdir_beam, name, maptype))
        
        # # plot!
        # for i, candname in enumerate(c.cand_name):
        #     logger.info("Plot slices {}/{}: {}".format(i, len(c.cand_name), candname))
        #     plot.plot_slices(src_name=candname, 
        #                      imagelist=imagelist, 
        #                      name="{}/{}_{}".format(outdir_beam, name, candname))
        



        
    # =====================
    # combine those three cand list to one
    logger.info("=========Combine catalogue==========")


    namelist = ['{}/{}_{}_cand.csv'.format(outdir_beam, name, maptype) for maptype in maplist ]

    plot.combine_csv(namelist, tablename="{}/{}_final".format(outdir_beam, name), 
                savevot=True)



    # =====================
    # plot final candidates 
    logger.info("========= Plotting =============")

    final_csv = "{}/{}_final.csv".format(outdir_beam, name)

    if os.path.exists(final_csv):
        p = Products(final_csv)
        p.generate_slices(imagelist=imagelist, 
                        savename='{}/{}_slices'.format(outdir_beam, name))
        p.generate_cutout(fitsname=deepimage, 
                        savename='{}/{}_deepcutout'.format(outdir_beam, name))
        p.generate_lightcurve(imagelist=imagelist, 
                            deepname=deepimage, 
                            savename='{}/{}_lightcurve'.format(outdir_beam, name))


    
    
    finalize_output(tar_name, outdir_beam)

    logger.info("====== Finished. =====")

if __name__ == "__main__":
    try:
        process_beam(14)
    except Exception:
        logger.exception("run into error {}: {}".format(sys.exc_info()[0], sys.exc_info()[1]))
        pass
   
        
    # for i in range(36):
        # process_beam(i)











