#!/usr/bin/env python
"""
Copyright (C) VAST 2024
"""
from vaster.structure import DataDir

import threading
import time
import glob
import sys
import os
import argparse
import logging

logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():
    parser = argparse.ArgumentParser(
        prog='DownloadASKAP', 
        description='Downloading ASKAP data, default will download visibilities and selavy catalogue', 
        epilog='Example usage: prepare_data 62268', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('sbids', type=int, nargs='+', help='input sbid for processing, number only')
    parser.add_argument('--beams', type=int, nargs='+', default=None, 
                        help='input beams for processing, number only')
    parser.add_argument('--dir', type=str, default='.', help='where those SBIDs folders are stored')
    parser.add_argument('--untar', action='store_true', help='untar visibilities')
    parser.add_argument('-p', '--parallel', type=int, default=1, help='parallel processing num')
    parser.add_argument('--no-vis', action='store_true', help='do not download visibilities')
    parser.add_argument('--no-selavy', action='store_true', help='do not download selavy')
    parser.add_argument('--mosaic', action='store_true', help='download mosaic images')
    parser.add_argument('--clean', action='store_true', help='Delete any relevant files before downloading')
    parser.add_argument('-v', '--verbose', action='store_true', help='make it verbose')
    parser.add_argument('--dry-run', action='store_true', help='perform a dry run, nothing will be downloaded')
    args = parser.parse_args()

    start_time = time.time()

    make_verbose(args)
    logger.info(args)
    logger.info('Total of %s SBIDs to prepare', len(args.sbids))

    beamlist = get_beamlist(args, )
    num_total_threads = len(args.sbids)*len(beamlist)
    num_threads = min(num_total_threads, args.parallel) 
    logger.info('Total of %s tasks', num_total_threads)
    logger.info('Number of parallel processing: %s', num_threads)

    # Initialize a list to keep track of active threads
    active_threads = []

    for i, sbid in enumerate(args.sbids):
        logger.info("Processing observation SB%s (%s/%s)", sbid, i+1, len(args.sbids))
        datadir= DataDir(sbid, args.dir)
        paths = datadir.paths
        logger.debug(paths)

        if args.no_selavy:
            logger.warning('SB%s: skip download selavy', sbid)
        else:
            if args.clean:
                clean_data(args, paths, sbid, affix="*.xml", command="rm")

            download_selavy(args, paths)

        for beam in beamlist:
            while len(active_threads) >= num_threads:
                for thread in active_threads:
                    if not thread.is_alive():
                        logger.info('%s thread finished', thread.name)
                        active_threads.remove(thread)
                time.sleep(1)  # Short sleep to avoid busy waiting

            name = f'SB{sbid}_{beam}'
            thread = threading.Thread(target=run, name=name, args=(args, paths, sbid, beam))
            active_threads.append(thread)
            thread.start()

    # Wait for all threads to complete
    for thread in active_threads:
        thread.join()

    end_time = time.time()
    measure_running_time(start_time, end_time)
    


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


def run(args, paths, sbid, beam):
    logger.info('SB%s: %s running', sbid, beam)
    if args.no_vis:
        logger.warning('SB%s: skip download visibilities', sbid)
    else:
        if args.clean:
            clean_data(args, paths, sbid, affix=f"*{beam}*.tar", command="rm")
            clean_data(args, paths, sbid, affix=f"*{beam}*.ms", command="rm -r")
            clean_data(args, paths, sbid, affix=f"*{beam}*.ms.corrected", command="rm -r")
            clean_data(args, paths, sbid, affix=f"*{beam}*.ms.corrected.flagversions", command="rm -r")

        download_visbility(args, paths, beam)
        if args.untar:
            untar_visibility(args, paths, beam)


def get_beamlist(args, num=36):
    if args.beams is None:
        beamlist = [f'beam{idx:02d}' for idx in range(num)]
    else:
        beamlist = [f'beam{idx:02d}' for idx in args.beams]
    return beamlist


def measure_running_time(start_time, end_time, nround=2):
    total_time = end_time - start_time
    if total_time <= 60:
        logger.info('Total running time %s seconds', round(total_time, nround))
    elif total_time <= 60*60:
        logger.info('Total running time %s minutes', round(total_time/60, nround))
    elif total_time <= 60*60*24:
        logger.info('Total running time %s hours', round(total_time/60/60, nround))
    else:
        logger.info('Total running time %s hours', round(total_time/60/60/24, nround))
        

def download_selavy(args, paths):
    fname = os.path.join(paths['path_scripts'], "download_selavy.sh")
    txt = 'bash ' + fname
    process_txt(args, txt)


def download_visbility(args, paths, beam):
    fname = os.path.join(paths['path_scripts'], f"bash_GETDATA_{beam}.sh")
    txt = 'bash ' + fname
    process_txt(args, txt)


def untar_visibility(args, paths, beam):
    fname = os.path.join(paths['path_scripts'], f"bash_UNTAR_{beam}.sh")
    txt = 'bash ' + fname
    process_txt(args, txt)


def clean_data(args, paths, sbid, affix, command):
    logger.debug('SB%s: clean %s', sbid, affix)
    fnamelist = glob.glob(os.path.join(paths['path_data'], affix))
    logger.debug('SB%s: found %s %s files', sbid, len(fnamelist), affix)
    for fname in fnamelist:
        txt = command + ' ' + fname
        process_txt(args, txt)
        


def process_txt(args, txt):
    if args.dry_run:
        logger.warning('Dry run: skip "%s"', txt)
    else:
        logger.info('Executing "%s"', txt)
        exit_code = os.system(txt)
        logger.debug('Exit with %s: %s', exit_code, txt)



if __name__ == '__main__':
    _main()
