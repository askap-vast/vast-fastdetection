#!/usr/bin/env python

import os
import sys
import glob
import math

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u

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
coord = sys.argv[-2]
sbid = sys.argv[-1]
num = 36 # number of beams

# just assume it follows standard naming convention 
#path = os.path.join('/o9000/ASKAP/VAST/fast_survey', "SB"+sbid)
path = os.path.join(os.getcwd(), "SB"+sbid)
path_models = os.path.join(path, 'models') 
path_images = os.path.join(path, 'images')
path_cand = os.path.join(path, 'candidates')

logging.info(path)
logging.info(path_models)
logging.info(path_images)
logging.info(path_cand)


def get_beam_position(imagename):
    '''
    Get beam coordinate and fwhm size for a certain fits (beam) file
    
    imagename: str 
        fits image of certain beam
    
    Return:
        ra: str in hourangle format, e.g. 20:50:34.450
        dec: str in degree format, e.g. -30:20:45.32
        fwhm: full width at half maximum, unit of degree (diameter)
    '''
    
    hdu = fits.open(imagename)[0].header
    beam_center = SkyCoord(hdu['CRVAL1'], 
                           hdu['CRVAL2'], unit=u.degree)
    fwhm = 3e8 / hdu['CRVAL3']/12 * 180/math.pi # unit of degree
    
    return beam_center, fwhm/2
    

# first check if this object is in the field, and if it is - which beam 
src = SkyCoord(coord, unit=(u.hourangle, u.degree))

beamlist = []
seplist = []
# run it in all beams to check if source is in there! 
for i in range(num):
    imagename = glob.glob(os.path.join(path_models, "*beam{:02d}*.fits".format(i)))
    
    if len(imagename) == 0:
        logging.info("ERROR: No fits file for beam{:02d}".format(i))
        continue
        
    elif len(imagename) != 1:
        logging.info('')
        logging.info('Number of matched images for beam{:02d} is'.format(i), len(imagename))
        logging.info(imagename[:3])
        
    imagename = imagename[0]
    beam_position, radius = get_beam_position(imagename)

    sep = src.separation(beam_position).degree
    if sep > 1.2*radius:
        continue
    
    logging.info(coord+' within beam{:02d} with separation {:.2f} degree'.format(i, sep))
    beamlist.append("beam{:02d}".format(i))
    seplist.append(sep)

logging.info('beamlist:')
logging.info(beamlist)
logging.info('seplist:')
logging.info(seplist)

if len(beamlist) == 0:
    sys.exit("Source is not in the field!!")


for i, beam in enumerate(beamlist):
    logging.info('Ploting the object in {} (Progress {}/{}...)'.format(beam, i+1, len(beamlist)))
    # just assume it follows standard naming convention 
    ## deep image
    deepimage = glob.glob(os.path.join(path_models, "*"+beam+"*.fits"))[0]
    ## catalgoue csv file contain interested source info 
    final_csv = glob.glob(os.path.join(path_cand, "*"+beam+"*_peak_cand.csv"))[0]
    ## the folder location with processed short images
    folder = path_images
    ## output folder 
    outdir = path_cand
    ## output file name 
    name = "SB" + sbid + '_' + beam + '_{:.2f}d'.format(seplist[i])

    ## ====================================
    ## get the imagelist with correct order
    imagelist = []
    for size in ["?", "??", "???", "????"]:
        tmp = glob.glob(os.path.join(folder, f'*{beam}_{size}.fits'))
        # tmp = glob.glob(folder + f'image_{size}.fits') # for FRB field 
        tmp.sort()
        imagelist += tmp

    imagelist = imagelist[:-1]

    if len(imagelist) == 0:
        logging.info("WARNING: No short images, skip {}".format(beam))
        continue

    logger.info("Loading foler {}".format(folder))
    logger.info("Processing {} images...".format(len(imagelist)))
    logger.info(imagelist)


    # =====================
    # plot final candidates 
    logger.info("========= Plotting =============")

    p = Products(final_csv, cand_name=[coord])
    p.generate_slices(imagelist=imagelist, 
                    savename='{}/{}_slices'.format(outdir, name))
    p.generate_cutout(fitsname=deepimage, 
                    savename='{}/{}_deepcutout'.format(outdir, name))
    p.generate_lightcurve(imagelist=imagelist, 
                        deepname=deepimage, 
                        savename='{}/{}_lightcurve'.format(outdir, name))


logger.info("====== Finished. =====")













