#!/usr/bin/env python
"""
Copyright (C) Swinburne 2024
"""
from vaster.vtools import measure_running_time 
from vaster.structure import DataBasic

import glob
import time
import os
import sys
import argparse

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
    parser.add_argument('-s', '--sbids', type=int, nargs='+', default=None, 
                        help='input sbid for prepare scripts, number only, emplty for all sbids in this folder')
    parser.add_argument('--dir', type=str, default='.', help='output directory')
    # parser.add_argument('--dry-run', action='store_true', help='Perform a dry run')
    parser.add_argument('-v', '--verbose', action='store_true',help='make it verbose')
    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)

    if args.sbids is None:
        folderlist = glob.glob(os.path.join(args.dir, "SB*"))
        logger.info(folderlist)
        sbidlist = [ int(folder.split('SB')[-1]) for folder in folderlist]
        sbidlist.sort()
        logger.info('Found %s SBIDs: %s', len(folderlist), sbidlist)
        args.sbids = sbidlist


    for i, sbid in enumerate(args.sbids):
        logger.debug("Processing observation SB%s (%s/%s)", sbid, i+1, len(args.sbids))
        databasic = DataBasic(sbid, args.dir)
        paths = databasic.paths
        nbeam = databasic.nbeam

        sbid_complete = check_sbid_compelte(args, sbid, paths, nbeam)
        if sbid_complete:
            num_cand = measure_final_candidates(args, sbid, paths)
            logger.info('SB%s completed: final_cand=%s', sbid, num_cand)


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
        


def check_sbid_compelte(args, sbid, paths, nbeam):
    peak_cand = glob.glob( os.path.join(paths['path_cand'], "*peak*cand.csv" ))
    logger.debug(peak_cand)
    chisq_cand = glob.glob( os.path.join(paths['path_cand'], "*chisquare*cand.csv" ))
    logger.debug(chisq_cand)

    if len(peak_cand) != nbeam:
        logger.warning('** SB%s does not complete: peak_cand=%s **', sbid, len(peak_cand))
        return False
    elif len(chisq_cand) != nbeam:
        logger.warning('** SB%s does not complete: chisq_cand=%s **', sbid, len(chisq_cand))
        return False
    else:
        return True
    

def measure_final_candidates(args, sbid, paths):
    num_cand = len(glob.glob( os.path.join(paths['path_cand'], "*lightcurve*.png") ))
    return num_cand


if __name__ == '__main__':
    _main()
