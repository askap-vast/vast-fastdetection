#!/usr/bin/env python
"""
Copyright (C) Swinburne 2024
"""
import os
import sys
import time
import argparse
import glob
import subprocess

import logging
logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        prog='UploadCands', 
        description='Upload realtime candidates', 
        epilog='Example usage: python upload_realtime_candidates.py -h', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('--basedir', type=str, default='.', help='Base directory to save output candidates')
    parser.add_argument('--untar', action='store_true', help='untar files')
    parser.add_argument('--rename', action='store_true', help='Rename candidate files to expected format')
    parser.add_argument('--rename-pattern', type=str, default='output_', )
    parser.add_argument('--mvfiles', action='store_true', help='move files')
    parser.add_argument('--upload', action='store_true', help='Upload cadidates files to web app')
    parser.add_argument('--clean', action='store_true', help='clean folder')
    parser.add_argument('--mvfolder', action='store_true', help='move folder')
    parser.add_argument('--skip-sbids', type=str, nargs='+', default=[], help='Skip these sbids, example 54028')
    parser.add_argument('--saved-sbids-txt', type=str, default='uploaded_sbids.txt', help='Save uploaded SBIDs')
    parser.add_argument("--sleep", type=int, default=600, help="Sleep time in seconds between checks")
    parser.add_argument("--token", type=str, default='', help='token to upload candidate files')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run')
    parser.add_argument('-v', '--verbose', action='store_true',help='make it verbose')

    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)


    # ===============
    # main program 
    # ===============

    while True:
        logger.info("")
        logger.info("=======================================================")
        logger.info("================ %s ==================", time.strftime("%Y-%m-%d %H:%M:%S"))
        sbids = gather_sbids(args)
        clean_sbids_list = clean_sbids(args, sbids)
        
        for sbid in clean_sbids_list:
            if not has_all_tarballs(args, sbid):
                continue

            if args.untar:
                untar_files(args, sbid)

            if args.rename:
                rename_files(args, sbid)

            if args.mvfiles:
                move_files(args, sbid)

            if args.upload:
                upload_files(args, sbid)

            if args.clean:
                remove_output_folders(args, sbid)

            if args.mvfolder:
                move_folder(args, sbid)
                    
        logger.info(f"Sleeping for {args.sleep} seconds...")
        time.sleep(args.sleep)

    end_time = time.time()
    measure_running_time(start_time, end_time)


def gather_sbids(args):
    '''
    Gather SBIDs in base directory
    '''
    sbids = [d for d in os.listdir(args.basedir) if os.path.isdir(os.path.join(args.basedir, d)) and d.isdigit()]
    logger.info(f"Found SBIDs: {sbids}")
    return sbids


def clean_sbids(args, sbids):
    if os.path.exists(args.saved_sbids_txt):
        with open(args.saved_sbids_txt, "r") as f:
            uploaded_sbids = set(f.read().splitlines())
    else:
        uploaded_sbids = set()
    cleaned_sbids = [sbid for sbid in sbids if sbid not in args.skip_sbids and sbid not in uploaded_sbids]

    logger.info(f"Skiping defined SBIDs: {args.skip_sbids}")
    logger.info(f"Skiping uploaded SBIDs: {uploaded_sbids}")
    logger.info(f"Cleaned SBIDs: {cleaned_sbids}")
    return cleaned_sbids


def has_all_tarballs(args, sbid, expected_count=36):
    sbid_dir = os.path.join(args.basedir, sbid)
    tar_files = [f for f in os.listdir(sbid_dir) if f.endswith(".tar.gz")]
    if len(tar_files) != expected_count:
        logger.warning("********************************")
        logger.warning("********************************")
        logger.warning(f"******************************** Skipping SBID {sbid}: found {len(tar_files)} tar.gz files, expected {expected_count}.")
        logger.warning("********************************")
        logger.warning("********************************")
        return False
    return True


def untar_files(args, sbid):
    sbid_dir = os.path.join(args.basedir, sbid)
    for file in os.listdir(sbid_dir):
        if file.endswith(".tar.gz"):
            file_path = os.path.join(sbid_dir, file)
            cmd = ["tar", "xvf", file_path, "-C", sbid_dir]
            logger.info(f"Executing: {' '.join(cmd)}")
            if not args.dry_run:
                subprocess.run(cmd, check=True)


def rename_files(args, sbid):
    sbid_dir = os.path.join(args.basedir, sbid)
    for root, _, files in os.walk(sbid_dir):
        for file in files:
            if args.rename_pattern in file:
                old_path = os.path.join(root, file)
                new_path = os.path.join(root, file.replace(args.rename_pattern, ""))
                logger.info(f"Renaming {old_path} -> {new_path}")
                if not args.dry_run:
                    os.rename(old_path, new_path)


def move_files(args, sbid):
    sbid_dir = os.path.join(args.basedir, sbid)
    candidates_dir = os.path.join(sbid_dir, "candidates")
    os.makedirs(candidates_dir, exist_ok=True)
    for subdir in os.listdir(sbid_dir):
        subdir_path = os.path.join(sbid_dir, subdir)
        if os.path.isdir(subdir_path) and subdir != "candidates":
            for file in os.listdir(subdir_path):
                file_path = os.path.join(subdir_path, file)
                cmd = ["mv", file_path, candidates_dir]
                logger.info(f"Executing: {' '.join(cmd)}")
                if not args.dry_run:
                    subprocess.run(cmd, check=True)


def upload_files(args, sbid):
    candidates_dir = os.path.join(args.basedir, sbid, "candidates")
    if not os.path.exists(candidates_dir):
        return
    
    cmd = [
        "python", "/home/yuwang/scripts/vaster_webapp/ywangvaster_webapp/upload_cand.py", 
        "--base_url", "http://vaster.duckdns.org:80", 
        "--token", args.token, 
        "--project_id", "realtime", 
        "--observation_id", f'SB{sbid}', 
        "--data_directory", candidates_dir
    ]
    logger.info(f"Executing: {' '.join(cmd)}")
    if not args.dry_run:
        subprocess.run(cmd, check=True)
        with open(args.saved_sbids_txt, "a") as f:
            logger.info(f"Write {sbid} into {args.saved_sbids_txt}")
            f.write(f"{sbid}\n")


def remove_output_folders(args, sbid):
    sbid_dir = os.path.join(args.basedir, sbid)
    for item in os.listdir(sbid_dir):
        item_path = os.path.join(sbid_dir, item)
        if os.path.isdir(item_path) and item.startswith("output"):
            cmd = ["rm", "-r", item_path]
            logger.info(f"Executing: {' '.join(cmd)}")
            if not args.dry_run:
                subprocess.run(cmd, check=True)


def move_folder(args, sbid):
    sbid_dir = os.path.join(args.basedir, sbid)
    sbid_dir_new = os.path.join(args.basedir, f"SB{sbid}")
    cmd = ["mv",  sbid_dir, sbid_dir_new]
    logger.info(f"Executing: {' '.join(cmd)}")
    if not args.dry_run:
        subprocess.run(cmd, check=True)



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


def create_dir(dir, ):
    if not os.path.exists(dir):
        os.makedirs(dir, exist_ok=True)
        logger.info('Create new directory %s', dir)


def measure_running_time(start_time, end_time, nround=2):
    total_time = end_time - start_time
    if total_time <= 60:
        logger.info('Running time %s seconds', round(total_time, nround))
    elif total_time <= 60*60:
        logger.info('Running time %s minutes', round(total_time/60, nround))
    elif total_time <= 60*60*24:
        logger.info('Running time %s hours', round(total_time/60/60, nround))
    else:
        logger.info('Running time %s days', round(total_time/60/60/24, nround))




if __name__ == '__main__':
    _main()
