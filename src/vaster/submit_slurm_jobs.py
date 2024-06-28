#!/usr/bin/env python
"""
Copyright (C) VAST 2024
"""
from vaster.structure import DataDir

import sys
import os
import argparse
import logging

logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():
    parser = argparse.ArgumentParser(
        prog='VOevent', 
        description='VO Event trigger', 
        epilog='Example usage: python ~/scripts/notebooks/notes/template.py -h', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('sbids', type=int, nargs='+', help='input sbid for processing, number only')
    parser.add_argument('--beams', type=int, nargs='+', default=None, 
                        help='input beams for processing, number only')
    parser.add_argument('--dir', type=str, default='.', help='where those SBIDs folders are stored')
    parser.add_argument('--mode', type=str, choices=['median', 'mean', 'max'], default='median', help='select coadd mode')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='make it verbose')
    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)
    logger.info('Total of %s SBIDs to prepare', len(args.sbids))

    for i, sbid in enumerate(args.sbids):
        logger.info("Processing observation SB%s (%s/%s)", sbid, i+1, len(args.sbids))
        datadir= DataDir(sbid, args.dir)
        args.paths = datadir.paths
        logger.debug(args.paths)



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




if __name__ == '__main__':
    _main()
