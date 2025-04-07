#!/usr/bin/env python
"""
Copyright (C) Swinburne 2024
"""
import time
import argparse
import os
import sys
import glob
import math

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u

from vaster.structure import DataBasic
from vaster.vtools import measure_running_time
from vaster.vastfast.plot import Products

import logging
logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        prog='VOevent', 
        description='VO Event trigger', 
        epilog='Example usage: python ~/scripts/notebooks/notes/template.py -h', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('coord', type=str, help='source coordinates')
    parser.add_argument('-s', '--sbids', type=int, nargs='+', help='input sbid for prepare scripts, number only')
    parser.add_argument('--dir', type=str, default='.', help='where those SBIDs folders are stored')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run')
    parser.add_argument('-v', '--verbose', action='store_true',help='make it verbose')
    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)

    for i, sbid in enumerate(args.sbids):
        args.sbid = sbid
        logger.info("Processing observation SB%s (%s/%s)", sbid, i+1, len(args.sbids))
        databasic= DataBasic(sbid, args.dir)
        args.paths = databasic.paths
        logger.debug(args.paths)
        args.nbeam = databasic.nbeam

        beamlist, seplist = check_src_in_which_beam(args, args.coord)
        plot_candidates(args, beamlist, seplist)


    end_time = time.time()
    measure_running_time(start_time, end_time)



def make_verbose(args):
    if args.verbose:
        logging.basicConfig(
            format='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
            level=logging.DEBUG,
            datefmt='%Y-%m-%d %H:%M:%S')
        logger.warning("verbose mode")
    else:
        logging.basicConfig(
            format='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')


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
    

def check_src_in_which_beam(args, coord):
    # first check if this object is in the field, and if it is - which beam 
    src = SkyCoord(coord, unit=(u.hourangle, u.degree))

    beamlist = []
    seplist = []
    # run it in all beams to check if source is in there! 
    for i in range(args.nbeam):
        imagename = []
        imagename += glob.glob(os.path.join(args.paths['path_models'], "*beam{:02d}*MFS-image.fits".format(i)))
        imagename += glob.glob(os.path.join(args.paths['path_models'], "*beam{:02d}*image.tt0.fits".format(i)))
        
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
            logging.info(src.to_string('hmsdms') + ' outside of beam{:02d} with separation {:.2f} degree'.format(i, sep))
            continue
        
        logging.info(src.to_string('hmsdms') +' within beam{:02d} with separation {:.2f} degree'.format(i, sep))
        beamlist.append("beam{:02d}".format(i))
        seplist.append(sep)

    logging.info('beamlist:')
    logging.info(beamlist)
    logging.info('seplist:')
    logging.info(seplist)

    if len(beamlist) == 0:
        sys.exit("Source is not in the field!!")

    return beamlist, seplist


def plot_candidates(args, beamlist, seplist):

    for i, beam in enumerate(beamlist):
        logging.info('Ploting the object in {} (Progress {}/{}...)'.format(beam, i+1, len(beamlist)))
        # just assume it follows standard naming convention 
        ## deep image
        imagename = []
        imagename += glob.glob(os.path.join(args.paths['path_models'], "*"+beam+"*image.tt0.fits"))
        imagename += glob.glob(os.path.join(args.paths['path_models'], "*"+beam+"*MFS-image.fits"))
        deepimage = imagename[0]
        logger.info('deepimage %s', deepimage)
        ## catalgoue csv file contain interested source info 
        final_csv = glob.glob(os.path.join(args.paths['path_cand'], "*"+beam+"*_peak_cand.csv"))[0]
        ## the folder location with processed short images
        folder = args.paths['path_images']
        ## output folder 
        outdir = args.paths['path_cand']
        ## output file name 
        name = "SB" + str(args.sbid) + '_' + beam + '_{:.2f}d'.format(seplist[i])

        ## ====================================
        ## get the imagelist with correct order
        imagelist = glob.glob(os.path.join(folder, f'*{beam}*image.fits'))
        imagelist.sort()

        if len(imagelist) == 0:
            logging.info("WARNING: No short images, skip {}".format(beam))
            continue

        logger.info("Loading foler {}".format(folder))
        logger.info("Processing {} images...".format(len(imagelist)))
        logger.info(imagelist)

        if args.dry_run:
            logger.info('Dry run - skip plotting SB%s %s', args.sbid, beam)
            continue


        # =====================
        # plot final candidates 
        logger.info("========= Plotting =============")

        p = Products(final_csv, cand_name=[args.coord])
        p.generate_slices(imagelist=imagelist, 
                        savename='{}/{}_slices'.format(outdir, name))
        p.generate_fits_cube(imagelist=imagelist, 
                        savename='{}/{}_slices'.format(outdir, name))

        p.generate_cutout(fitsname=deepimage, 
                        savename='{}/{}_deepcutout'.format(outdir, name))
        p.generate_fits_cutout(fitsname=deepimage, 
                        savename='{}/{}_deepcutout'.format(outdir, name))
                        
        p.generate_lightcurve(imagelist=imagelist, 
                            deepname=deepimage, 
                            savename='{}/{}_lightcurve'.format(outdir, name))


    logger.info("====== Finished. =====")




if __name__ == '__main__':
    _main()







