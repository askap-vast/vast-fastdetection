#!/usr/bin/env python
"""
Copyright (C) Swinburne 2024
"""
import os
import time
import argparse
import glob
import subprocess
import requests
import pandas as pd

import random
import gspread
from oauth2client.service_account import ServiceAccountCredentials

from astropy.io import fits

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
    parser.add_argument('--notify-slack', action='store_true', help='Send a slack message when finish uploading')
    parser.add_argument('--clean', action='store_true', help='clean output tmp folder')
    parser.add_argument('--deep-clean', action='store_true', help='clean candidates folder, start with a fresh run')
    parser.add_argument('--mvfolder', action='store_true', help='move folder')
    parser.add_argument('--skip-sbids', type=str, nargs='+', default=[], help='Skip these sbids, example 54028')
    parser.add_argument('--saved-sbids-txt', type=str, default='uploaded_sbids.txt', help='Save uploaded SBIDs')
    parser.add_argument('--skipped-sbids-txt', type=str, default='skipped_sbids.txt', help='Skipped SBIDs')
    parser.add_argument("--sleep", type=int, default=600, help="Sleep time in seconds between checks")
    parser.add_argument("--token", type=str, default='', help='token to upload candidate files')
    parser.add_argument("--webhook", type=str, default='', help='WEBHOOK URL to send message in slack')
    parser.add_argument('--force', action='store_true', help='force to process SBIDs even it doesnt have 36 beams')
    parser.add_argument('--assign', action='store_true', help='automatic assign people to classify candidates')
    parser.add_argument('--select', action='store_true', help='select promissing candidates with defined filtering metrics')
    parser.add_argument('--metadata', action='store_true', help='read fits metadata')
    parser.add_argument('--update-sheet', action='store_true', help='automatic update info to google spreadsheet')
    parser.add_argument('--skip-projects', type=str, nargs='+', default=[], help='Skip project id example AS116; only works with --metadata')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run')
    parser.add_argument('-v', '--verbose', action='store_true',help='make it verbose')

    args = parser.parse_args()

    make_verbose(args)
    logger.info(args)

    filtering = (
        'md_deep > 0.1'                # modulation index > 0.1
        ' & '
        '('
            'deep_sep_arcsec < 2'      # with a crossmatch (remove sidelobes)
            ' | '
            'deep_sep_arcsec > 20'     # or isolated (no crossmatch)
            ' | '
            'deep_peak_flux < 0.002'   # or a faint source < 2mJy nearby 
        ')'
    )

    # ===============
    # main program 
    # ===============

    loop_count = 0

    while True:
        logger.info("")
        logger.info("===================== Loop %s =========================", loop_count)
        logger.info("================ %s ==================", time.strftime("%Y-%m-%d %H:%M:%S"))
        sbids = gather_sbids(args)
        clean_sbids_list = clean_sbids(args, sbids)

        sheet_id = '1xd1h4k9GtlAEH4TkUEDBQ6uYGvGNaDhEeB4Xf8WAiIw'
        creds_path='/fred/oz330/realtime/vaster-471803-f266a065caa2.json'

        if loop_count % 50 == 0:
            heartbeat_slack(args)

        loop_count += 1
        
        for sbid in clean_sbids_list:
            args.sbid_dir = os.path.join(args.basedir, sbid)
            args.candidates_dir = os.path.join(args.sbid_dir, "candidates")
            args.candidates_all_dir = os.path.join(args.sbid_dir, "candidates_all")
            args.sbid_dir_new = os.path.join(args.basedir, f"SB{sbid}")

            if args.deep_clean:
                deep_clean(args, sbid)

            if args.assign:
                uids, last_sbids, unames = get_available_slack_user_ids(sheet_id, creds_path, sub_sheet_idx=0)
                logger.info("%s: Available Users: %s", sbid, unames)
                logger.info("%s: Available User IDs: %s", sbid, uids)
                logger.info("%s: Last SBIDs: %s", sbid, last_sbids)
            else:
                uids = None

            if not has_all_tarballs(args, sbid) and not args.force:
                logger.info(f"{sbid}: Sleeping for 60 seconds...")
                time.sleep(60)
                continue

            if args.untar:
                untar_files(args, sbid)

            if args.rename:
                rename_files(args, sbid)

            if args.mvfiles:
                move_files(args, sbid)

            ncands = len(glob.glob(os.path.join(args.candidates_dir, "*slices*gif")))
            logger.info('Total number of candidates (before filtering): %s', ncands)

            if args.metadata:
                fits_path = glob.glob(os.path.join(args.candidates_dir, f"SB{sbid}_beam*_std.fits"))[0]
                logger.info('%s: Read metadata from fits %s', sbid, fits_path)
                metadata = read_fits_metadata(fits_path)
                if metadata['project_id'] in args.skip_projects:
                    logger.warning("********************************")
                    logger.warning("********************************")
                    logger.warning("******************************** %s: Skip sbid %s in project %s", sbid, sbid, metadata['project_id'])
                    logger.warning("********************************")
                    logger.warning("********************************")
                    with open(args.skipped_sbids_txt, "a") as f:
                        logger.info(f"Write {sbid} into {args.skipped_sbids_txt}")
                        f.write(f"{sbid}\n")
                    continue

            if args.select:
                ncands = select_cands(args, sbid, filtering)
                logger.info('Filtering metrics: %s', filtering)
                logger.info('Total candidates (after filtering): %s', ncands)

            if args.upload:
                upload_files(args, sbid)

            if uids is not None:
                uid = select_next_user(uids, last_sbids)
                uname = get_uname_from_uid(uid, uids, unames)
                logger.info('%s: Assign %s to user %s with ID %s', sbid, sbid, uname, uid)
            else:
                uid = None
                uname = None

            if args.notify_slack:
                notify_slack(args, sbid, ncands, userid=uid)

            if uid is not None and args.update_sheet:
                update_user_assignment(sheet_id, creds_path, uid, sbid, sub_sheet_idx=0)

            if args.update_sheet:
                write_fits_metadata_to_google_sheet(args, sbid, metadata, sheet_id, creds_path, uname, ncands, sub_sheet_idx=1)

            if args.clean:
                remove_output_folders(args, sbid)

            if args.mvfolder:
                move_folder(args, sbid)

            logger.info(f"{sbid}: Sleeping for {args.sleep} seconds...")
            time.sleep(args.sleep)
                    
        logger.info(f"Sleeping for {args.sleep} seconds...")
        time.sleep(args.sleep)

    end_time = time.time()
    measure_running_time(start_time, end_time)


def gather_sbids(args):
    '''
    Gather SBIDs in base directory
    '''
    sbids = [d for d in os.listdir(args.basedir) if os.path.isdir(os.path.join(args.basedir, d)) and d.isdigit()]
    sbids.sort()
    logger.info(f"Found {len(sbids)} SBIDs: {sbids}")
    return sbids


def clean_sbids(args, sbids):
    if os.path.exists(args.saved_sbids_txt):
        with open(args.saved_sbids_txt, "r") as f:
            uploaded_sbids = list(set(f.read().splitlines()))
            uploaded_sbids.sort()

        args.uploaded_sbids = uploaded_sbids

    else:
        uploaded_sbids = list(set())

    if os.path.exists(args.skipped_sbids_txt):
        with open(args.skipped_sbids_txt, "r") as f:
            skipped_sbids = list(set(f.read().splitlines()))
            skipped_sbids.sort()
    else:
        skipped_sbids = list(set())

    cleaned_sbids = [sbid for sbid in sbids if sbid not in args.skip_sbids and sbid not in uploaded_sbids and sbid not in skipped_sbids]
    cleaned_sbids.sort()

    logger.info(f"Skiping {len(args.skip_sbids)} args defined SBIDs: {args.skip_sbids}")
    logger.info(f"Skiping {len(skipped_sbids)} txt defined SBIDs: {skipped_sbids}")
    logger.info(f"Skiping {len(uploaded_sbids)} txt uploaded SBIDs: {uploaded_sbids}")
    logger.info(f"Cleaned SBIDs: {cleaned_sbids}")
    return cleaned_sbids


def has_all_tarballs(args, sbid, expected_count=36):
    sbid_dir = args.sbid_dir
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
    sbid_dir = args.sbid_dir
    for file in os.listdir(sbid_dir):
        if file.endswith(".tar.gz"):
            file_path = os.path.join(sbid_dir, file)
            cmd = ["tar", "xvf", file_path, "-C", sbid_dir]
            logger.info(f"Executing: {' '.join(cmd)}")
            if not args.dry_run:
                subprocess.run(cmd, check=True)


def rename_files(args, sbid):
    sbid_dir = args.sbid_dir
    for root, _, files in os.walk(sbid_dir):
        for file in files:
            if args.rename_pattern in file:
                old_path = os.path.join(root, file)
                new_path = os.path.join(root, file.replace(args.rename_pattern, ""))
                logger.info(f"Renaming {old_path} -> {new_path}")
                if not args.dry_run:
                    os.rename(old_path, new_path)


def move_files(args, sbid):
    sbid_dir = args.sbid_dir
    candidates_dir = args.candidates_dir
    os.makedirs(candidates_dir, exist_ok=True)

    for subdir in os.listdir(sbid_dir):
        subdir_path = os.path.join(sbid_dir, subdir)
        if os.path.isdir(subdir_path) and subdir.startswith("output"):
            for file in os.listdir(subdir_path):
                file_path = os.path.join(subdir_path, file)
                cmd = ["mv", file_path, candidates_dir]
                logger.info(f"Executing: {' '.join(cmd)}")
                if not args.dry_run:
                    subprocess.run(cmd, check=True)


def select_cands(args, sbid, filtering):
    """Copy *final.csv files from candidates_dir to candidates_all_dir, then filter in-place.

    filtering : str
        A pandas.DataFrame.query string, e.g., "SNR > 10 and chi2 < 3".
    """
    candidates_dir = args.candidates_dir
    if not os.path.exists(candidates_dir):
        logger.warning(f"select_cands: candidates_dir not found: {candidates_dir}")
        return

    candidates_all_dir = args.candidates_all_dir
    os.makedirs(candidates_all_dir, exist_ok=True)

    ncands = 0

    final_csv_fnames = glob.glob(os.path.join(candidates_dir, "*final.csv"))
    for src in final_csv_fnames:
        dst = os.path.join(candidates_all_dir, os.path.basename(src))
        if os.path.exists(dst):
            logger.info(f"Skipping copy, already exists: {dst}")
            continue
        cmd = ["cp", src, dst]
        logger.info(f"Executing: {' '.join(cmd)}")
        if not args.dry_run:
            subprocess.run(cmd, check=True)

    for fpath in final_csv_fnames:
        try:
            df = pd.read_csv(fpath)
        except Exception as e:
            logger.error(f"Failed to read {fpath}: {e}")
            continue

        try:
            filtered = df.query(filtering) 
        except Exception as e:
            logger.error(f"Bad filtering expression '{filtering}': {e}")
            continue

        logger.info(f"Filtering {fpath}: {len(df)} -> {len(filtered)} cands")
        ncands += len(filtered)

        if not args.dry_run:
            try:
                filtered.to_csv(fpath, index=False)
            except Exception as e:
                logger.error(f"Failed to write filtered CSV to {fpath}: {e}")

    return ncands


def upload_files(args, sbid):
    candidates_dir = args.candidates_dir
    if not os.path.exists(candidates_dir):
        logger.warning(f"upload_files: candidates_dir not found: {candidates_dir}")
        return
    
    cmd = [
        "python", "/fred/oz330/software/vaster_webapp/ywangvaster_webapp/upload_cand.py", 
        "--base_url", "http://vaster.duckdns.org:80", 
        "--token", args.token, 
        "--project_id", f'SB{sbid}', 
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
    sbid_dir = args.sbid_dir
    for item in os.listdir(sbid_dir):
        item_path = os.path.join(sbid_dir, item)
        if os.path.isdir(item_path) and item.startswith("output"):
            cmd = ["rm", "-r", item_path]
            logger.info(f"Executing: {' '.join(cmd)}")
            if not args.dry_run:
                subprocess.run(cmd, check=True)

def deep_clean(args, sbid):
    cmd = ["rm", "-r", args.candidates_dir]
    logger.info(f"Executing: {' '.join(cmd)}")
    if not args.dry_run:
        subprocess.run(cmd, check=True)

    cmd = ["rm", "-r", args.candidates_all_dir]
    logger.info(f"Executing: {' '.join(cmd)}")
    if not args.dry_run:
        subprocess.run(cmd, check=True)


def move_folder(args, sbid):
    sbid_dir = args.sbid_dir
    sbid_dir_new = args.sbid_dir_new
    cmd = ["mv",  sbid_dir, sbid_dir_new]
    logger.info(f"Executing: {' '.join(cmd)}")
    if not args.dry_run:
        subprocess.run(cmd, check=True)


def notify_slack(args, sbid, ncands, userid=None):
    candidates_dir = args.candidates_dir
    nbeam = len(glob.glob(os.path.join(candidates_dir, "*peak.fits")))

    if userid is not None:
        mention = f"<@{userid}>"  # ← put actual Slack user ID here
    else:
        mention = ""

    message = (
        f"----------------\n"
        f"*** SB{sbid} FINISHED:    beams={nbeam}    cands={ncands} ***    {mention} \n"
        f"----------------"
    )

    payload = {"text": message}
    if not args.dry_run:
        response = requests.post(args.webhook, json=payload)
        if response.status_code != 200:
            logger.warning("Error posting to Slack:", response.status_code, response.text)
        else:
            logger.info("Message sent to Slack: %s", message)
    else:
        logger.info("Dry run: skip posting message %s", message)


def heartbeat_slack(args):
    """
    Send a Slack heartbeat message at specified times to confirm the system is running.
    """
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    msg = f"{timestamp} - System is running normally. Total of {len(args.uploaded_sbids)} observations have been processed. "

    payload = {"text": msg}

    if not args.dry_run:
        response = requests.post(args.webhook, json=payload)
        if response.status_code != 200:
            logger.warning(f"Error posting heartbeat to Slack: {response.status_code} {response.text}")
        else:
            logger.info(f"Heartbeat sent to Slack: {msg}")
    else:
        logger.info(f"Dry run: skip posting heartbeat {msg}")


def get_available_slack_user_ids(sheet_id, creds_path, sub_sheet_idx=0):
    """
    Fetch available Slack User IDs and their Last SBID from a Google Sheet.

    Parameters
    ----------
    sheet_id : str
        The Google Sheet ID.
    creds_path : str
        Path to the service account JSON credentials.

    Returns
    -------
    uids : List[str]
        Slack User IDs with Availability == 'Y'.
    last_sbids : List[str]
        Last SBID values (can be empty or numeric) for each user.
    """
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name(creds_path, scope)
    client = gspread.authorize(creds)

    sheet = client.open_by_key(sheet_id).get_worksheet(sub_sheet_idx)
    data = sheet.get_all_records()

    uids = []
    last_sbids = []
    unames = []

    for row in data:
        if str(row.get("Availability", "")).strip().upper() == "Y":
            uids.append(str(row.get("Slack User ID")).strip())
            last_sbids.append(str(row.get("Last SBID", "")).strip())
            unames.append(str(row.get("Name")).strip())

    return uids, last_sbids, unames


def get_uname_from_uid(uid, uids, unames):
    try:
        index = uids.index(uid)
        return unames[index]
    except ValueError:
        return uid  # fallback if uid not found


def select_next_user(uids, last_sbids):
    """
    Select a user for assignment:
    - If any user has empty Last SBID, randomly choose one.
    - Otherwise, pick user with smallest Last SBID (numerically).

    Parameters
    ----------
    uids : List[str]
    last_sbids : List[str]

    Returns
    -------
    str
        Selected Slack User ID.
    """
    candidates_empty = [uid for uid, sbid in zip(uids, last_sbids) if not sbid]
    logger.info("Users with empty last sbids %s", candidates_empty)
    if candidates_empty:
        return random.choice(candidates_empty)

    # Convert SBIDs to integers and find the user with the smallest one
    valid_sbid_users = [(uid, int(sbid)) for uid, sbid in zip(uids, last_sbids) if sbid.isdigit()]
    if not valid_sbid_users:
        return random.choice(uids)  # fallback in case of bad data

    return min(valid_sbid_users, key=lambda x: x[1])[0]


def update_user_assignment(sheet_id, creds_path, user_id, sbid, sub_sheet_idx=0):
    """
    Update the sheet for a given user: increment Assigned, append SBID to SBIDs, and set Last SBID.

    Parameters
    ----------
    sheet_id : str
        The Google Sheet ID.
    creds_path : str
        Path to service account credentials.
    user_id : str
        Slack User ID to update.
    sbid : str or int
        SBID to assign.
    """
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name(creds_path, scope)
    client = gspread.authorize(creds)

    sheet = client.open_by_key(sheet_id).get_worksheet(sub_sheet_idx)
    data = sheet.get_all_records()

    header = sheet.row_values(1)  # column names
    col_idx = lambda name: header.index(name) + 1  # 1-based indexing for gspread

    for idx, row in enumerate(data):
        if str(row.get("Slack User ID", "")).strip() == user_id:
            row_num = idx + 2  # +2 for header row and 0-based index

            # Update Assigned
            current_assigned = int(row.get("Assigned", 0))
            sheet.update_cell(row_num, col_idx("Assigned"), current_assigned + 1)

            # Update SBIDs
            existing_sbids = str(row.get("SBIDs", "")).strip()
            new_sbids = f"{existing_sbids} {sbid}".strip()
            sheet.update_cell(row_num, col_idx("SBIDs"), new_sbids)

            # Update Last SBID if empty or sbid is newer
            current_last_sbid = str(row.get("Last SBID", "")).strip()
            if not current_last_sbid:
                sheet.update_cell(row_num, col_idx("Last SBID"), str(sbid))
            else:
                try:
                    if int(sbid) > int(current_last_sbid):
                        sheet.update_cell(row_num, col_idx("Last SBID"), str(sbid))
                except ValueError:
                    # In case non-integer content exists in the cell
                    logger.warning("Non-integer Last SBID found for %s, skipping update.", user_id)

            break

    logger.info('Updating spreadsheet %s finished. ', sub_sheet_idx)


def read_fits_metadata(fits_path):
    # --- Read FITS Header ---
    with fits.open(fits_path) as hdul:
        hdr = hdul[0].header
    
    metadata = {
        'field_ra': hdr.get('FIELDRA'), 
        'field_dec': hdr.get('FIELDDEC'), 
        'field_name': hdr.get('FIELD'), 
        'project_id': hdr.get("PROJECT"), 
        'obs_time': hdr.get("DATE-OBS"), 
        'cent_freq': hdr.get("CRVAL3") / 1e6,     # convert from Hz to MHz
        'bmaj': hdr.get("BMAJ") * 3600,           # convert from degree to arcsec 
        'bmin': hdr.get('BMIN') * 3600,           # convert from degree to arcsec 
        'bpa': hdr.get("BPA"),
    }

    logger.info('------ Metadata -------')
    logger.info(metadata)

    return metadata


def write_fits_metadata_to_google_sheet(args, sbid, meta_dict, sheet_id, creds_path, uname, ncands, sub_sheet_idx=1):
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name(creds_path, scope)
    client = gspread.authorize(creds)

    sheet = client.open_by_key(sheet_id).get_worksheet(sub_sheet_idx)

    candidates_dir = args.candidates_dir

    # Format values to write
    row = [
        sbid, 
        meta_dict.get("project_id", ""), 
        "'" + meta_dict.get("field_name", ""),
        "'" + meta_dict.get("field_ra", ""),
        "'" + meta_dict.get("field_dec", ""),
        meta_dict.get("obs_time", ""),
        meta_dict.get("cent_freq", ""),              # MHz
        meta_dict.get("bmaj", ""),                   # arcsec
        meta_dict.get("bmin", ""),
        meta_dict.get("bpa", ""),
        ncands, 
        uname, 
    ]

    # Append row
    sheet.append_row(row, value_input_option='USER_ENTERED')
    logger.info('Updating spreadsheet %s finished. ', sub_sheet_idx)


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


if __name__ == "__main__":
    max_restarts = 10
    restart_count = 0

    while restart_count < max_restarts:
        try:
            _main()
            break  # if main() exits cleanly, break the loop
        except Exception as e:
            restart_count += 1
            logger.error(f"Program crashed with error: {e}")
            logger.info(f"Restarting in 10 seconds... (attempt {restart_count}/{max_restarts})")
            time.sleep(10)

    if restart_count >= max_restarts:
        logger.error("Maximum number of restarts reached. Exiting.")
