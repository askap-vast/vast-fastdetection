#!/usr/bin/env python
"""
Copyright (C) VAST 2024
"""
from vaster.structure import DataBasic
from vaster.vtools import measure_running_time, process_txt

import subprocess

import glob
import time
import sys
import os
import argparse
import logging

logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>, Raghav Girgaonkar <raghav@uwm.edu>"


def _main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        prog='SubmitSlurm', 
        description='Submit slurm jobs', 
        epilog='Example usage: submit_slurm_jobs 62043 -b 0', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('sbids', type=int, nargs='+', help='input sbid for processing, number only')
    parser.add_argument('-b', '--beams', type=int, nargs='+', default=None, 
                        help=r'''input beams for processing, number only, leave it blank for all of beams;
                                using {3..10} to constract a list''')
    parser.add_argument('--dir', type=str, default='.', help='where those SBIDs folders are stored')
    parser.add_argument('--steps', type=str, nargs='+', default=['FIXDATA', 'MODELING', 'IMGFAST', 'SELCAND', 'CLNDATA'], 
                        help='tasks to process, following the order')
    parser.add_argument('--nodes', type=str, default='', help='Specify the nodes to use, e.g. "node01,node02"')
    parser.add_argument('--sleep', type=float, default=0, help='Job submission sleep time between two sbids, unit of hours')
    parser.add_argument('--clean', action='store_true', help='Delete any relevant files before re-submit')
    parser.add_argument('--dry-run', action='store_true', help='perform a dry run')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='make it verbose')
    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)
    logger.info('Total of %s SBIDs to prepare', len(args.sbids))

    for i, sbid in enumerate(args.sbids):
        logger.info("Processing observation SB%s (%s/%s)", sbid, i+1, len(args.sbids))
        databasic = DataBasic(sbid, args.dir)
        args.databasic = databasic 
        args.paths = databasic.paths
        logger.debug(args.paths)

        if args.beams is None:
            beams = [idx for idx in range(databasic.nbeam)]
        else: 
            beams = args.beams

        for idx in beams:
            fnamelist = extract_joblist(args, idx)
            logger.info(f'SB{sbid} beam{idx:02d}: find {len(fnamelist)} jobs for submission')
            if args.clean:
                clean_data(args, sbid, affix=f"*beam{idx:02d}*", command="rm -r")

            if args.dry_run:
                logger.warning(f'Dry run: SB{sbid} beam{idx:02d}: will submit below scripts in order')
                logger.warning(fnamelist)
            else:
                job_id_list = submit_joblist(fnamelist,nodes=args.nodes)
                write_scancel_scripts(args, idx, job_id_list)
        
        if i + 1 < len(args.sbids) and args.sleep > 0:
            sleep_seconds = int(args.sleep * 3600)
            logger.info('Sleep %s hours (= %s seconds) before submitting the next sbid...', args.sleep, sleep_seconds)
            time.sleep(sleep_seconds)

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
        

def extract_joblist(args, idx, ):
    fnames = args.databasic.files[idx]
    fnames = fnames[fnames['function'] == 'run']

    fnamelist = []
    for step in args.steps:
        fname = fnames[fnames['step'] == step]['fname'][0]
        if not os.path.isfile(fname):
            logger.error('file %s does not exist', fname)
            break

        fnamelist.append(fname)

    logger.debug(fnamelist)
    return fnamelist


def write_scancel_scripts(args, idx, job_id_list):
    savename = os.path.join(args.paths['path_scripts'], f'kill_beam{idx:02d}_jobs.sh')
    logger.info('Writing scancel scripts to %s', savename)
    with open(savename, 'w') as fw:
        cmd = 'scancel ' + ' '.join(job_id_list)
        logger.info(cmd)
        fw.write(cmd)



def clean_data(args, sbid, affix, command):
    logger.debug('SB%s: clean %s', sbid, affix)
    fnamelist = []
    if 'FIXDATA' in args.steps:
        fnamelist += glob.glob(os.path.join(args.paths['path_data'], affix+".corrected"))
    if 'MODELING' in args.steps:
        fnamelist += glob.glob(os.path.join(args.paths['path_models'], affix))
    if 'IMGFAST' in args.steps:
        fnamelist += glob.glob(os.path.join(args.paths['path_images'], affix))
    if 'SELCAND' in args.steps:
        fnamelist += glob.glob(os.path.join(args.paths['path_cand'], affix))
    logger.debug('SB%s: found %s %s files', sbid, len(fnamelist), affix)
    for fname in fnamelist:
        txt = command + ' ' + fname
        process_txt(args, txt)
        


def submit_joblist(fnamelist,nodes=None):
    job_id_list = []
    for i, fname in enumerate(fnamelist):
        if i == 0:
            job_id = submit_onejob(fname,nodes=nodes)
        else:
            job_id = submit_onejob(fname, nodes=nodes, dependency=True, dep_job_id=pre_job_id)

        if job_id is None:
            break 

        job_id_list.append(job_id)
        pre_job_id = job_id
        time.sleep(1)

    return job_id_list



def submit_onejob(fname, nodes=None,dependency=False, dep_job_id=None):
    if nodes and (('MODELING' in fname) or ('FIXDATA' in fname)):
        if dependency:
            cmd = ['sbatch', f'--nodelist={nodes}', '-d', f'afterok:{dep_job_id}', fname]
        else:
            cmd = ['sbatch', '--nodelist=', f'{nodes}', fname]
    else:
        if dependency:
            cmd = ['sbatch', '-d', f'afterok:{dep_job_id}', fname]
        else:
            cmd = ['sbatch', fname]

    logger.info('Executing "%s"', cmd)
    job = subprocess.run(cmd, capture_output=True, text=True)
    if job.returncode == 0:
        job_id = job.stdout.strip().split()[-1]
        logger.info('Job %s submitted successfully with ID %s', fname, job_id)
        return job_id
    
    else:
        logger.warning('Failed to submit job %s', fname)
        logger.warning(f'Error: {job.stderr.strip()}')
        return None



if __name__ == '__main__':
    _main()
