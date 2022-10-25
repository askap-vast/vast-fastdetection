import logging

import os
import glob
import tarfile
import  shutil

from .cube import Cube, Filter
from .utils import *

logger = logging.getLogger(__name__)

def get_sigcube(imagelist):
    """get sigcube"""
    logger.info("============")
    logger.info("Starting to build the cube...")
    cube = Cube(imagelist)
    cube.icube()
    logger.info("cube shape: {}".format(cube.sigcube.shape))
    logger.info("Finish to create the cube.")
    logger.info("============")
    logger.info("Remove bad images...")
    cube.remove_bad_images()
    logger.info("cube shape: {}".format(cube.sigcube.shape))
    return cube.sigcube
    

def get_map(sigcube, ktypelist, imagelist):
    """get maps of given ktype"""
    f = Filter(sigcube)

    for ktype in ktypelist:
        logger.info("===== Matched Filter =====")
        logger.info("Kernel match filter '{}'...".format(ktype))
        f.fmap(ktype, width=4)
        logger.info("Kernel match Done")
        
        f.tofits(fitsname="{}/{}_{}.fits".format(outdir, name, ktype), imagename=imagelist[0])
        logger.info("Save the results to {}_{}.fits".format(name, ktype))


def finalize_output(tar_name, dir):
    """remove temporary files and generate a tarfile"""
    logger.info("Removing fits files...")
    remove_fits(dir)

    logger.info("Creating tarfile...")
    tar_file(tar_name, dir)

    logger.info("Cleaning up output directory...")
    shutil.rmtree(dir)

