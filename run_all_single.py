#!/usr/bin/env python

import sys
import os
import glob

from vastfast.cube import Filter
from vastfast import plot
from vastfast.plot import Candidates, Products
from vastfast.input import read_shortimage, read_deepcat, read_deepimage
from vastfast.procedures import *
from vastfast.exceptions import *
from vastfast.setting import INPUT_PATH, OUT_DIR, OUT_PREFIX, KTYPELIST

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
sh = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
sh.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(sh)



## ========= arguments ================
beam = sys.argv[-1]
beam = int(beam)

def main():
   
    logger.info("============")
    logger.info("Reading input files...")
    # input data: short images
    try:
        imagelist = read_shortimage(INPUT_PATH, beam)
    except NoInputError:
        logger.error("At least 2 images are needed; skip the beam{:02}...".format(beam))
        return -1
    
    # input data: catalogue
    try:
        catalogue = read_deepcat(INPUT_PATH, beam)
    except NoInputError:
        logger.error("Deep catalogue is not available; skip the beam{:02}...".format(beam))
        return -1

    # input data: deep image
    deepimage = read_deepimage(INPUT_PATH, beam)
    
    # process beam
    p = Procedures(imagelist, beam, catalogue, OUT_DIR, OUT_PREFIX, KTYPELIST, deepimage)
    p.run_all_steps()

      
if __name__ == "__main__":
    main()
    
   
        
   










