#!/usr/bin/env python
"""
Copyright (C) VAST 2024
"""
from vaster.structure import DataBasic

import pandas as pd
import glob
import sys
import os
import argparse
import logging

logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():
    parser = argparse.ArgumentParser(
        prog='InspectHPCUsage', 
        description='Inspect HPC usage', 
        epilog='Example usage: inspect_usage slurm_example.usage', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('-f', '--fname', type=str, nargs='+', default=None, 
                        help='usage statistics generate from slurm sacct output')
    parser.add_argument('-s', '--sbids', type=int, nargs='+', help='SBID, numbers')
    parser.add_argument('-b', '--beams', type=int, nargs='+', default=None, 
                        help='input beams for checking, number only')
    parser.add_argument('--steps', type=str, nargs='+', default=['FIXDATA', 'MODELING', 'IMGFAST', 'SELCAND', 'CLNDATA'], 
                        help='tasks to process, following the order')
    parser.add_argument('--dir', type=str, default='.', help='where those SBIDs folders are stored')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='make it verbose')
    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)

    colnames = ['MaxRSS', 'MaxVMSize']
    beamlist = get_beamlist(args, num=36)

    if args.fname is None:
        logger.info('Total of %s SBIDs to inspect', len(args.sbids))
        fname_list = []
        for i, sbid in enumerate(args.sbids):
            logger.info("Processing observation SB%s (%s/%s)", sbid, i+1, len(args.sbids))
            fname_sbid = get_logs(args, sbid, beamlist)
            fname_list += fname_sbid
            logger.info('SB%s: Found %s usage logs', sbid, len(fname_sbid))
        logger.info('Found total of %s usage logs', len(fname_list))
        logger.debug(fname_list)
    else:
        fname_list = args.fname

    # pd.options.display.max_columns = None
    for fname in fname_list:
        tb = read_tb(fname)
        tb = convert_mem_unit(tb, colnames)
        tb = clean_tb(tb)
        logger.info("=======")
        logger.info(fname)
        print(tb)


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


def get_beamlist(args, num=36):
    if args.beams is None:
        beamlist = [f'beam{idx:02d}' for idx in range(num)]
    else:
        beamlist = [f'beam{idx:02d}' for idx in args.beams]
    return beamlist


def get_logs(args, sbid, beamlist):
    databasic= DataBasic(sbid, args.dir)
    fname_list = []
    for beam in beamlist:
        for step in args.steps:
            fname_list += glob.glob(os.path.join(databasic.paths['path_logs'], f'*{step}*{beam}*.usage'))
    return fname_list


def clean_tb(tb, skip_jobs=['batch', 'extern']):
    # get rid of .batch and .extern jobs 
    pattern = '|'.join(skip_jobs)
    # Filter the DataFrame
    df_filtered = tb[~tb["JobID"].str.contains(pattern, na=False)]
    return df_filtered


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
