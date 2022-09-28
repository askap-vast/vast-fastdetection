#!/usr/bin/env python

import sys
import os
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

# import dask
# dask.config.set({"multiprocessing.context": "fork"})
# from dask.distributed import Client

import multiprocessing






# input data
# INPUT_PATH = "../test_data/SB9596_beam14_casa/"
INPUT_PATH = "../test_data/SB12704/"
# output folder
outdir = "./output"
if not os.path.exists(outdir):
    os.mkdir(outdir)

# output file prefix
out_prefix = "output"

# @profile
def process_beam(beam):
    input_path = INPUT_PATH
    # input data: short images
    image_files = input_path +"/" + "beam{:02}*".format(beam)
    print("image name: ", image_files)
    imagelist = sorted(glob.glob(image_files))

    logger.info("Loading foler {}".format(input_path))
    logger.info("Processing {} images...".format(len(imagelist)))
    logger.info(imagelist)
    
    # input data: catalogue
    catalogue_file = input_path + "/" +"*beam{:02}*.vot".format(beam)
    catalogue = glob.glob(catalogue_file)[0]

    # TODO: RAISE WARNING IF CATALOGUE HAVE MORE THAN ONE FILE

    # output file prefix
    name = out_prefix +"_beam{:02}".format(beam)
    
    # get the psf list with correct order
    #psflist = glob.glob(folder + 'beam??_?.psf.fits') + glob.glob(folder + 'beam??_??.psf.fits')
    #logger.info("Processing {} psf...".format(len(psflist)))
    #logger.info(psflist)

    # get the significance cube (do smooth with each 2d image)
    #ktype = 'gaussian'
    logger.info("============")
    logger.info("Starting to build the cube...")

    cube = Cube(imagelist)
    #cube.icube(ktype, 19, 19)
    #cube.save_oricube()
    cube.icube()

    logger.info(cube.sigcube.shape)
    #logger.info(cube.oricube.shape)
    logger.info("Finish to create the cube.")
    logger.info("============")

    logger.info("Remove bad images...")
    cube.remove_bad_images()
    logger.info(cube.sigcube.shape)

    # get the matched filter in time axis
    f = Filter(cube.sigcube)
    ##f = Filter(cube.oricube)

    logger.info("===== Matched Filter =====")
    ktype = "chisquare"
    logger.info("Kernel match filter '{}'...".format(ktype))
    f.fmap(ktype, width=1)
    logger.info("Kernel match Done")

    fitsname = "{}/{}_{}.fits".format(outdir, name, ktype)
    if os.path.exists(fitsname):
        os.remove(fitsname)

    f.tofits(fitsname=fitsname, imagename=imagelist[0])
    logger.info("Save the results to {}_{}.fits".format(name, ktype))



    logger.info("======== Matched Filter ========")
    ktype = "peak"
    logger.info("Kernel match filter '{}'...".format(ktype))
    f.fmap(ktype, width=1)
    logger.info("Kernel match Done")

    fitsname = "{}/{}_{}.fits".format(outdir, name, ktype)
    if os.path.exists(fitsname):
        os.remove(fitsname)

    f.tofits(fitsname=fitsname, imagename=imagelist[0])
    logger.info("Save the results to {}_{}.fits".format(name, ktype))



    logger.info("======== Matched Filter ========")
    ktype = "std"
    logger.info("Kernel match filter '{}'...".format(ktype))
    f.fmap(ktype, width=1)
    logger.info("Kernel match Done")

    fitsname = "{}/{}_{}.fits".format(outdir, name, ktype)
    if os.path.exists(fitsname):
        os.remove(fitsname)

    f.tofits(fitsname=fitsname, imagename=imagelist[0])
    logger.info("Save the results to {}_{}.fits".format(name, ktype))



    logger.info("======== Matched Filter ========")
    ktype = "gaussian"
    logger.info("Kernel match filter '{}'...".format(ktype))
    #f.fmap(ktype, width=4)
    logger.info("Kernel match Done")

    fitsname = "{}/{}_{}.fits".format(outdir, name, ktype)
    if os.path.exists(fitsname):
        os.remove(fitsname)

    #f.tofits(fitsname=fitsname, imagename=imagelist[0])
    logger.info("Save the results to {}_{}.fits".format(name, ktype))




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
    c.local_max(min_distance=30, sigma=5, data='peak')
    logger.info("Find local maximum done. ")

    # plot a map
    #c.plot_fits(fitsname=chisq_map, imagename="{}/{}_map1".format(outdir, name))

    # read deep information to select candidates
    #catalogue = glob.glob(folder + "*_{}*comp.vot".format(beam))
    

    logger.info("Deep image catalogue {}".format(catalogue))
    c.select_candidates(deepcatalogue=catalogue)

    # save the table

    tablename = "{}/{}_cand".format(outdir, name)
    if os.path.exists(tablename):
        os.remove(tablename)
    try:
        c.save_csvtable(tablename="{}/{}_cand".format(outdir, name), savevot=True)
        for i, candname in enumerate(c.cand_name):
            logger.info("Plot slices {}/{}: {}".format(i, len(c.cand_name), candname))
            plot.plot_slices(src_name=candname, 
                        imagelist=imagelist, 
                        name="{}/{}_{}".format(outdir, name, candname))

    except TypeError:
        logger.info("no candidates found")

    # plot a final map
    #c.plot_fits(fitsname=chisq_map, imagename="{}/{}_map2".format(outdir, name))

    # plot!
    


    logger.info("====== Finished. =====")
    return 1


if __name__ == "__main__":
    # process_beam(INPUT_PATH,14)

    # results = []
    # beams = [8,16]
    # for beam in beams:
    #     results.append(process_beam(INPUT_PATH, beam))

    # dask.compute(results)
    pool = multiprocessing.Pool(processes=4)
    pool.map(process_beam, list(range(36)), chunksize=9)
    # print(type(ll[0]))
   




