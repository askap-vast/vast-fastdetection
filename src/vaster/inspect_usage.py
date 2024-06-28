#!/usr/bin/env python
"""
Copyright (C) VAST 2024
"""
import pandas as pd

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
    parser.add_argument('tb', type=str, nargs='+', help='usage statistics generate from slurm output')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='make it verbose')
    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)

    colnames = ['MaxRSS', 'MaxVMSize']

    for usage in args.tb:
        tb = read_tb(usage)
        tb = convert_mem_unit(tb, colnames)
        logger.info(tb)




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

def read_tb(fname):
    return pd.read_csv(fname, delimiter='|')


def convert_mem_unit(tb, colnames):
    for colname in colnames:
        newcol = []
        
        for mem in tb[colname]:
            logger.debug('mem: %s', mem)
            if pd.isna(mem) or mem == '0' or isinstance(mem, float):
                newcol.append(mem)
                continue
        
            value, unit = float(mem[:-1]), mem[-1]
            logger.debug('mem: %s; value: %s; unit: %s', mem, value, unit)
            
            if unit == 'K' and value >= 1024**2:
                new_value = round(value / 1024**2, 1)
                new_unit = 'G'
                newcol.append(str(new_value) + new_unit)
            elif unit == 'K' and value >= 1024:
                new_value = round(value / 1024, 1)
                new_unit = 'M'
                newcol.append(str(new_value) + new_unit)
            else:
                newcol.append(mem)

        tb[colname] = newcol

    return tb



if __name__ == '__main__':
    _main()
