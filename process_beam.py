#!/usr/bin/env python
import os
import sys
import logging
import configparser


from vastfast.input import read_shortimage, read_deepcat, read_deepimage
from vastfast.procedures import create_out_beam, Procedures
from vastfast.exceptions import *

## ========= arguments ========= 
beam = sys.argv[-1]
beam = int(beam)

## ========= read in config file =========
config = configparser.ConfigParser()		
config.read("/code/run_vastfast.ini")

input_path = config["PATH"]["input_dir"]

out_dir = config["PATH"]["out_dir"]

nprocess = int(config["RESOURCE"]["n_cpu"])


## ========= override defaults from docker config ===========

try:
    config = configparser.ConfigParser()
    config.read("/input/app_settings.ini")
except Exception as e:
    print(e)
else:
    try:
        nprocess = int(config["RESOURCE"]["n_cpu"])
    except:
        pass

    try:
        beam = int(config["PROCESS"]["beam"])
    except:
        pass

    try:
        input_directory_name = config["PROCESS"]["input_directory_name"]
        input_path = os.path.join(input_path, input_directory_name)
    except:
        pass


## ========= beam output directory ========= 
outdir_beam = create_out_beam(out_dir, beam)

## =========  set up logger ========= 
logger = logging.getLogger()
logger.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
file_handler = logging.FileHandler(outdir_beam + '/beam{:02}.log'.format(beam))
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
        imagelist = read_shortimage(input_path, beam)
    except NoInputError:
        logger.error("At least 2 images are needed; skip the beam{:02}...".format(beam))
        return -1
    
    # input data: catalogue
    try:
        catalogue = read_deepcat(input_path, beam)
    except NoInputError:
        logger.error("Deep catalogue is not available; skip the beam{:02}...".format(beam))
        return -1

    # input data: deep image
    deepimage = read_deepimage(input_path, beam)
    
    # process beam
    p = Procedures(imagelist, beam, catalogue, out_dir, outdir_beam, deepimage, nprocess=nprocess)
    p.run_all_steps()

      
if __name__ == "__main__":
    main()
    
   
        
   










