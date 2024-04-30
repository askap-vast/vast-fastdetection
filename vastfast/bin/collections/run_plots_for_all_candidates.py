#!/usr/bin/env python

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


# ## ========= input ================
# ## deep image
# deepimage = sys.argv[-6]
# ## catalgoue csv file contain interested source info 
# catalogue = sys.argv[-5]
# ## the folder location with processed short images
# folder = sys.argv[-4]
# ## beam number - format of "beam00"
# beam = sys.argv[-3]
# ## output folder 
# outdir = sys.argv[-2]
# ## output file name 
# name = sys.argv[-1]

# a csv file in format SOURCE, SBID, FIELD, BEAM, ...
filename = sys.argv[-2] 
outdir = sys.argv[-1]
# an example fake csv file to match the format 
final_csv = '/o9000/ASKAP/VAST/fast_test/detection_results/possum_transient.csv'


source_info = []

# read file
with open(filename, 'r') as f:
    for i, line in enumerate(f):
        # dont need the header line
        if i == 0:
            continue

        line = line.rstrip('\n')
        # SOURCE, SBID, FIELD, BEAM = line.split(',')
        source_info.append(line.split(','))


print('Reading Finished... ')




for i, (SOURCE, SBID, FIELD, beam, SEP) in enumerate(source_info):

    if SBID != "SB11816":
        continue 

    folder = '/o9000/ASKAP/VAST/fast_test/{}/images/{}/'.format(FIELD, SBID)
    name = '{}_{}_{}_{:.3f}d'.format(SOURCE, SBID, beam, float(SEP))
    deepimage = '/o9000/ASKAP/VAST/fast_test/{}/cutout_models/{}_{}.image.tt0.fits'.format(FIELD, SBID, beam)

    ## ====================================
    ## get the imagelist with correct order
    imagelist = []
    for size in ["?", "??", "???", "????"]:
        tmp = glob.glob(folder + f'{beam}_{size}.fits')
        # tmp = glob.glob(folder + f'image_{size}.fits') # for FRB field 
        tmp.sort()
        imagelist += tmp
    
    # imagelist = imagelist[:-1]
    
    logger.info("Loading foler {}".format(folder))
    logger.info("Processing {} images...".format(len(imagelist)))
    logger.info(imagelist)
    
    
    # =====================
    # plot final candidates 
    logger.info("========= Plotting =============")
    
    # final_csv = "{}/{}_final.csv".format(outdir, name)
    # final_csv = catalogue
    
    try: 

        p = Products(final_csv)
    
        p.cand_name = [SOURCE]
    
        p.generate_slices(imagelist=imagelist, 
                      savename='{}/{}_slices'.format(outdir, name))
        p.generate_cutout(fitsname=deepimage, 
                      savename='{}/{}_deepcutout'.format(outdir, name))
        p.generate_lightcurve(imagelist=imagelist, 
                          deepname=deepimage, 
                          savename='{}/{}_lightcurve'.format(outdir, name))
    except:
        logger.info("{} {} skipped. ".format(SBID, beam))
        continue
    
    
    
    
    logger.info("====== Finished. =====")













