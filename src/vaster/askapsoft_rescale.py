#!/usr/bin/env python
"""
Copyright (C) VAST 2024
"""
from casacore.tables import *
import numpy as np

# import sys
import os
import argparse
import logging

logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"

def main():
    parser = argparse.ArgumentParser(
        prog='askapsoft_rescale', 
        description='rescale data processed by ASKAPsoft', 
        epilog='Example usage: askapsoft_rescale ms_in ms_out', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('ms_in', type=str, help='input sbid for prepare scripts, number only')
    parser.add_argument('ms_out', type=str, help='input sbid for prepare scripts, number only')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='make it verbose')
    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)

    # ms_file = sys.argv[-2]
    # ms_file_out = sys.argv[-1]
    ms_file = args.ms_in
    ms_file_out = args.ms_out

    text = "cp -R %s %s" % (ms_file, ms_file_out)
    logger.info('copy data: %s', text)
    os.system(text)
    logger.info('copy data finished')

    t = table(ms_file_out, readonly=False, ack=False)

    logger.info('Rescaling data...')
    nrows = t.nrows()
    for row in range(nrows):
        logger.debug("%d/%d" % (row, nrows))

        if(row % 1000 == 0):
            logger.info("%d/%d" % (row, nrows))

        cdata = t.getcol("DATA", startrow=row, nrow=1)
        cdata *= 2.0
        t.putcol("DATA", cdata, startrow=row, nrow=1)
    t.close()

    logger.info('Rescaling data finished. ')


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



if __name__ == '__main__':
    main()




