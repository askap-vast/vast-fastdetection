#!/usr/bin/env python

import sys
import logging


from vastfast.input import read_shortimage, read_deepcat, read_deepimage
from vastfast.procedures import create_out_beam, Procedures
from vastfast.exceptions import *
from vastfast.setting import INPUT_PATH, OUT_DIR, OUT_PREFIX, KTYPELIST


## ========= arguments ========= 
beam = sys.argv[-1]
beam = int(beam)

## ========= beam output directory ========= 
outdir_beam = create_out_beam(OUT_DIR, beam)

## =========  set up logger ========= 
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
file_handler = logging.FileHandler(outdir_beam + '/log.log', 'w')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(file_handler)
logger.addHandler(stream_handler)



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
    p = Procedures(imagelist, beam, catalogue, OUT_DIR, OUT_PREFIX, outdir_beam, KTYPELIST, deepimage)
    p.run_all_steps()

      
if __name__ == "__main__":
    main()
    
   
        
   










