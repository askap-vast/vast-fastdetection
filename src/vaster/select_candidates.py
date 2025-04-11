#!/usr/bin/env python
"""
Copyright (C) VAST 2024
"""
import os
import sys
import glob
import argparse
import logging
import yaml

import warnings 
from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning

from vaster.vastfast.cube import Cube, Filter
from vaster.vastfast import plot
from vaster.vastfast.plot import Candidates, Products

logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():
    parser = argparse.ArgumentParser(
        prog='SELCAND', 
        description='select candidates from VASTER short images products', 
        epilog='Example usage: select_candidates.py -h', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('--outdir', type=str, default='.', help='output directory')
    parser.add_argument('--deepimage', type=str, help='deep image file')
    parser.add_argument('--catalogue', type=str, help='deep catalogue file')
    parser.add_argument('--folder', type=str, help='folder location with processed short images')
    parser.add_argument('--beam', type=str, help='beam number, format of "beam00"')
    parser.add_argument('--name', type=str, help='output file name')
    parser.add_argument('--config', type=str, default='config.yml', help='configuration file')
    parser.add_argument('--ignore-warning', action='store_true', help='suppress various astropy warning')
    parser.add_argument('-v', '--verbose', action='store_true', help='make it verbose')
    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)

    if args.ignore_warning:
        suppress_warnings()

    config = read_config(args.config)
    imagelist = get_imagelist(args)

    if isinstance(config['CANDIDATES']['CHUNK'], int):
        chunk_size = int(config['CANDIDATES']['CHUNK'])
        logger.info('Running in chunk mode, chunk size is %s', chunk_size)

        for i in range(0, len(imagelist), chunk_size):
            chunk = imagelist[i:i + chunk_size]
            oname = f'{args.name}_chunk{i}_{i+chunk_size}'
            logger.info('Processing chunk %s to %s for total %s images', i, i+chunk_size, len(imagelist))
            logger.info('    with oname %s', oname)
            run(args, config, chunk, oname)

    else:
        oname = args.name
        run(args, config, imagelist, oname)
    

    combine_cand_csv(args, config)
    plot_final_candidates(args, imagelist)

    logger.info("====== Finished. =====")


def run(args, config, imagelist, oname):
    f, num = get_cube(imagelist)
    generate_statistical_fits(args, config, f, imagename=imagelist[0], oname=oname)

    chisquare_map, peak_map, std_map, gaussian_map = read_statistical_fits(args, oname=oname)
    select_local_maximum(args, config, num, chisquare_map, peak_map, std_map, gaussian_map, oname=oname)



def make_verbose(args):
    if args.verbose:
        logger.warning("verbose mode")
        logging.basicConfig(
            format='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
            level=logging.DEBUG,
            datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(
            format='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')
        

def read_config(fname):
    # read configuration file 
    if os.path.isfile(fname):
        with open(fname, 'r') as yaml_file:
            config = yaml.safe_load(yaml_file)
        return config
    
    else:
        logger.error('Cannot find configuration file')
        logger.error('%s does not exists', fname)
        sys.exit()


def suppress_warnings():
    logger.warning('** NOTICE: IGNORE ASTROPY WARNING **')
    warnings.filterwarnings('ignore', category=AstropyWarning, append=True)
    warnings.filterwarnings('ignore',
                        category=AstropyDeprecationWarning, append=True)



def get_imagelist(args, ):
    imagelist = glob.glob(os.path.join(args.folder, f'*{args.beam}*image.fits'))
    imagelist.sort()

    logger.info("Loading foler {}".format(args.folder))
    logger.info("Processing {} images...".format(len(imagelist)))
    logger.info(imagelist)

    return imagelist


def get_psflist(args, ):
    # ===================================
    # get the psf list with correct order
    psflist = glob.glob(args.folder + 'beam??_?.psf.fits') + glob.glob(args.folder + 'beam??_??.psf.fits')
    logger.info("Processing {} psf...".format(len(psflist)))
    logger.info(psflist)
    # ===================================
    return psflist


def get_cube(imagelist):
    logger.info("============")
    logger.info("Starting to build the cube...")
    cube = Cube(imagelist)

    ## get the significance cube (do smooth with each 2d image)
    # ktype = 'gaussian' # 'psf'
    # cube.icube(ktype, 19, 19)

    cube.icube()
    logger.info(cube.sigcube.shape)
    logger.info("Finish to create the cube.")
    logger.info("============")

    logger.info("Remove bad images...")
    cube.remove_bad_images()
    logger.info(cube.sigcube.shape)
    num = cube.sigcube.shape[0]

    ## ====================================
    ## get the matched filter in time axis
    f = Filter(cube.sigcube)
    return f, num



def generate_statistical_fits(args, config, f, imagename, oname):
    ktypelist = config['CANDIDATES']['KTYPE']
    for ktype in ktypelist:
        logger.info("===== Matched Filter =====")
        logger.info("Kernel match filter '{}'...".format(ktype))

        f.fmap(ktype, width=config['CANDIDATES']['GAUSSIAN_WIDTH'])
        logger.info("Kernel match Done")
        
        fitsname = os.path.join(args.outdir, f'{oname}_{ktype}.fits')
        f.tofits(fitsname=fitsname, imagename=imagename)
        logger.info("Save the results to %s", fitsname)


def read_statistical_fits(args, oname):
    logger.info("======== Select candidates ==========")
    # read fits
    chisquare_map = os.path.join(args.outdir, oname+'_chisquare.fits')
    peak_map = os.path.join(args.outdir, oname+'_peak.fits')
    std_map = os.path.join(args.outdir, oname+'_std.fits')
    gaussian_map = os.path.join(args.outdir, oname+'_gaussian.fits')
    return chisquare_map, peak_map, std_map, gaussian_map


def select_local_maximum(args, config, num, chisquare_map, peak_map, std_map, gaussian_map, oname):
    maplist = config['CANDIDATES']['MAP']
    ## =============== select candidates =================
    for maptype in maplist:
        if 'gaussian' not in maplist:
            c = Candidates(chisquare_map, peak_map, std_map, num=num)
        else:
            ## include Gaussian map during candidates selection 
            c = Candidates(chisquare_map, peak_map, std_map, gaussian_map=gaussian_map, 
                        num=num)
            
        sigma = config['CANDIDATES']['THRESHOLD_SIGMA'][maptype.upper()]
        log_sigma = config['CANDIDATES']['LOG_THRESHOLD_SIGMA'][maptype.upper()]
        min_distance = config['CANDIDATES']['MIN_DISTANCE']
        tabletype = config['SOURCE_FINDER']

        sep_deep = config['CANDIDATES']['SEP_DEEP']
        mdlim = config['CANDIDATES']['MD_LIMIT']
        extlim = config['CANDIDATES']['EXT_LIMIT']
        beamlim = config['CANDIDATES']['SEP_BEAM']
        bright = config['CANDIDATES']['BRIGHT_LIMIT']

        # get theortical threshold
        if maptype == 'chisquare':
            threshold = c.get_threshold_chisquare(sigma=sigma, num=num)
        elif maptype == 'peak':
            threshold = c.get_threshold_peak(sigma=sigma, num=num)
        else:
            threshold = sigma

        # find local maximum
        logger.info("Find local maximum....")
        c.local_max(min_distance=min_distance, sigma=log_sigma, threshold=threshold, data=maptype)
        logger.info("Find local maximum done. ")
        
        ## plot a map with all of candidates above the threshold 
        imagename = os.path.join(args.outdir, f'{oname}_{maptype}_map1')
        c.plot_fits(fitsname=vars()[maptype+'_map'], 
                    imagename=imagename)
            
        logger.info("Deep image catalogue {}".format(args.catalogue))
        c.select_candidates(deepcatalogue=args.catalogue, tabletype=tabletype, 
                            sep=sep_deep, mdlim=mdlim, extlim=extlim, beamlim=beamlim, 
                            bright=bright)
        
        # save the table
        tablename = os.path.join(args.outdir, f'{oname}_{maptype}_cand')
        c.save_csvtable(tablename=tablename, savevot=config['CANDIDATES']['SAVEVOT'])
        
        ## plot a final map with promising candidates 
        imagename = os.path.join(args.outdir, f'{oname}_{maptype}_map2')
        c.plot_fits(fitsname=vars()[maptype+'_map'], 
                    imagename=imagename)
        

def combine_cand_csv(args, config):
    maplist = config['CANDIDATES']['MAP']
    # =====================
    # combine those three cand list to one
    logger.info("=========Combine catalogue==========")
    # namelist = [os.path.join(args.outdir, f'{args.name}_{maptype}_cand.csv') for maptype in maplist ]
    namelist = []
    for maptype in maplist:
        pattern = os.path.join(args.outdir, f'{args.name}*{maptype}_cand.csv')
        namelist += glob.glob(pattern)

    tablename = os.path.join(args.outdir, f'{args.name}_final')
    logger.info('Combine tables:')
    logger.info(namelist)
    plot.combine_csv(namelist, tablename=tablename, savevot=config['CANDIDATES']['SAVEVOT'])



def plot_final_candidates(args, imagelist):
    # =====================
    # plot final candidates 
    logger.info("========= Plotting =============")
    final_csv = os.path.join(args.outdir, f'{args.name}_final.csv')

    if os.path.exists(final_csv):
        p = Products(final_csv, limit=50)
        savename = os.path.join(args.outdir, f'{args.name}_slices')
        p.generate_slices(imagelist=imagelist, savename=savename)
        p.generate_fits_cube(imagelist=imagelist, savename=savename)
        
        savename = os.path.join(args.outdir, f'{args.name}_deepcutout')
        p.generate_cutout(fitsname=args.deepimage, savename=savename)
        p.generate_fits_cutout(fitsname=args.deepimage, savename=savename)

        savename = os.path.join(args.outdir, f'{args.name}_lightcurve')
        p.generate_lightcurve(imagelist=imagelist, 
                            deepname=args.deepimage, 
                            savename=savename)
    else:
        logger.info("Final csv {} not exists. ".format(final_csv))



if __name__ == '__main__':
    _main()




