import logging

import os
import glob
import tarfile
import  shutil

from .cube import Cube, Filter
from .plot import Candidates
from .utils import *
from .exceptions import *

logger = logging.getLogger(__name__)

def get_sigcube(imagelist):
    """get sigcube"""
    logger.info("============")
    logger.info("Building the cube...")
    cube = Cube(imagelist)
    cube.icube()
    logger.info("Finish buidling the cube.")
    logger.info("============")
    logger.info("Removing bad images...")
    cube.remove_bad_images()
    logger.info("cube shape: {}".format(cube.sigcube.shape))
    return cube.sigcube
    
def get_filter(sigcube):
    f = Filter(sigcube)
    return f

def get_map_all(f, ktypelist, imagelist, outdir, name):
    """get maps of given ktype"""
    for ktype in ktypelist:
        try:
            get_map(f, ktype, imagelist, outdir, name)
        except Exception:
            logger.exception("Fail to generate {} map".format(ktype))
        
def get_map(f, ktype, imagelist, outdir, name):
    logger.info("===== Matched Filter =====")
    logger.info("Kernel match filter '{}'...".format(ktype))
    f.fmap(ktype, width=4)
    logger.info("Kernel match Done")
    f.tofits(fitsname="{}/{}_{}.fits".format(outdir, name, ktype), imagename=imagelist[0])
    logger.info("Save the results to {}_{}.fits".format(name, ktype))

def read_map_list(ktypelist, outdir_beam, name):
    maps = {}
    for ktype in ktypelist:
        map_name = ktype+ "_map"
        maps[map_name] = "{}/{}_{}.fits".format(outdir_beam, name, ktype)

    return maps

def get_candidates(map_dict, outdir_beam, name, catalogue):
    valid_map = Candidates(map_dict).valid_map
    for maptype in valid_map:
        _get_candidates_single_map(maptype, map_dict, outdir_beam, name, catalogue)
        
def _get_candidates_single_map(maptype, map_dict, outdir_beam, name, catalogue):
    c = Candidates(map_dict)

    # find local maximum
    logger.info("Finding local maximum for {}...".format(maptype))
    c.local_max(maptype=maptype, min_distance=30, sigma=5)
    logger.info("Finish finding local maximum. ")
    
    
    # plot a map with all of candidates above the threshold 
    c.plot_fits(fitsname=map_dict[maptype], 
            imagename="{}/{}_{}1".format(outdir_beam, name, maptype))
    
    logger.info("Deep image catalogue {}".format(catalogue))
    c.select_candidates(deepcatalogue=catalogue)
    
    # check 
    # save the table
    c.save_csvtable(tablename="{}/{}_{}_cand".format(outdir_beam, name, maptype), savevot=True)
    logger.info("Save candidates in table.")

    # plot a final map with promising candidates 
    c.plot_fits(fitsname=map_dict[maptype], 
            imagename="{}/{}_{}2".format(outdir_beam, name, maptype))



def finalize_output(tar_name, dir):
    """remove temporary files and generate a tarfile"""
    logger.info("Removing fits files...")
    remove_fits(dir)

    logger.info("Creating tarfile...")
    tar_file(tar_name, dir)

    logger.info("Cleaning up output directory...")
    shutil.rmtree(dir)

