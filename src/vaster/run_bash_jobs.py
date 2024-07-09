#!/usr/bin/env python
"""
Copyright (C) VAST 2024
"""
from vaster.structure import DataBasic
from vaster.vtools import measure_running_time, process_txt
from vaster.prepare_data import download_selavy

import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import glob
import time
import sys
import os
import argparse
import logging

logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        prog='RunBashJobs', 
        description='Run bash jobs', 
        epilog='Example usage: run_bash_jobs 62043 -b 0', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('sbids', type=int, nargs='+', help='input sbid for processing, number only')
    parser.add_argument('-b', '--beams', type=int, nargs='+', default=None, 
                        help='input beams for processing, number only, leave it blank for all of beams')
    parser.add_argument('--dir', type=str, default='.', help='where those SBIDs folders are stored')
    parser.add_argument('-p', '--parallel', type=int, default=1, help='parallel processing num')
    parser.add_argument('--clean', action='store_true', help='Delete any relevant files before re-submit')
    parser.add_argument('--dry-run', action='store_true', help='perform a dry run')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='make it verbose')
    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)
    logger.info('Total of %s SBIDs to prepare', len(args.sbids))

    job_scripts = []

    for i, sbid in enumerate(args.sbids):
        logger.info("Processing observation SB%s (%s/%s)", sbid, i+1, len(args.sbids))
        databasic = DataBasic(sbid, args.dir)
        args.databasic = databasic 
        args.paths = databasic.paths
        logger.debug(args.paths)

        download_selavy(args, args.paths)

        if args.beams is None:
            beams = [idx for idx in range(databasic.nbeam)]
        else: 
            beams = args.beams

        for idx in beams:
            fname = extract_jobname(args, idx)
            logger.info(f'SB{sbid} beam{idx:02d}: find {fname} for submission')
            if args.clean:
                clean_data(args, sbid, affix=f"*beam{idx:02d}*", command="rm -r")

            job_scripts.append(fname)

    logger.info('Total of %s jobs', len(job_scripts))
    logger.info(job_scripts)
    

    if args.dry_run:
        logger.warning('Dry run: Skip submitting below %s scripts', len(job_scripts))
        logger.warning(job_scripts)
    else:       
        
        with ProcessPoolExecutor(max_workers=args.parallel) as executor:
            # Submit all jobs to the executor
            futures = {executor.submit(run_job, script): script for script in job_scripts}

            # As each job completes, print its result
            for future in as_completed(futures):
                script = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    logger.warning(f"{script} generated an exception: {exc}")
                else:
                    logger.warning(f"{script} completed with result:\n{result}")


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
        

def extract_jobname(args, idx, ):
    fname = os.path.join( args.paths['path_scripts'], f'bash_PROCESSING_beam{idx:02d}.sh' )
    return fname



def clean_data(args, sbid, affix, command):
    logger.debug('SB%s: clean %s', sbid, affix)
    fnamelist = []
    fnamelist += glob.glob(os.path.join(args.paths['path_data'], affix))
    fnamelist += glob.glob(os.path.join(args.paths['path_models'], affix))
    fnamelist += glob.glob(os.path.join(args.paths['path_images'], affix))
    logger.debug('SB%s: found %s %s files', sbid, len(fnamelist), affix)
    for fname in fnamelist:
        txt = command + ' ' + fname
        process_txt(args, txt)
        

def run_job(script_name):
    log_file = script_name.replace('.sh', '.log')
    log_file = log_file.replace('/scripts/', '/logfiles/')
    logger.info('Saving to logfile %s', log_file)
    try:
        with open(log_file, 'w') as log:
            process = subprocess.Popen(['bash', script_name], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                print(line, end='')  # Print each line of output immediately
                log.write(line)
            process.wait()
        return f"Success: Output saved to {log_file}"
    except subprocess.CalledProcessError as e:
        with open(log_file, 'a') as log:
            log.write(f"\nError: {e.stderr}")
        print(f"Error running {script_name}: {e.stderr}")
        return f"Error: Output saved to {log_file}"


if __name__ == '__main__':
    _main()
