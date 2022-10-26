import glob
import logging
import os

from .exceptions import *

logger = logging.getLogger(__name__)


def read_shortimage(path, beam):
    """Read in short images"""
    simage_files = path +"/" + "beam{:02}*".format(beam)
    simagelist = sorted(glob.glob(simage_files))
    if len(simagelist) > 1:
        logger.info("Load folder: {}".format(path))
        logger.info("Process {} short images".format(len(simagelist)))
        # logger.info(imagelist)
    else:
        raise NoInputError()
        
    return simagelist
    

def read_deepcat(path, beam):
    """Read in deep catalogue"""
    cat_file = path + "/" + "*beam{:02}*.vot".format(beam)
    cats = sorted(glob.glob(cat_file))
    if len(cats) == 1:
        cat = cats[0]
        logger.info("Load deep catalogue: {}".format(cat))
        return cat

    if len(cats) < 1:
        raise NoInputError()

    if len(cats) > 1:
        cat = cats[0]
        logger.warning("Multiple deep catalogues found; the frist one is used by default")
        logger.info("Load deep catalogue: {}".format(cat))
        return cat


def read_deepimage(path, beam):
    """Read in deep catalogue"""
    dimage_file = path + "/" + "*beam{:02}*.tt0.fits".format(beam)
    dimages = sorted(glob.glob(dimage_file))
    if len(dimages) == 1:
        dimage = dimages[0]
        logger.info("Load deep image: {}".format(dimage))
        return dimage

    if len(dimages) == 0:
        logger.warning("Deep image is not available")
        return None

    if len(dimages) > 1:
        dimage = dimages[0]
        logger.warning("Multiple deep images found; the frist one is used by default")
        logger.info("Load deep dimage: {}".format(dimage))
        return dimage


