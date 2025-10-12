#!/usr/bin/env python
"""
Author: Yuanming Wang
Copyright (C) VAST 2024
"""
import time
import sys
import os
import argparse
import logging
import numpy as np

import smtplib
from email.message import EmailMessage

from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u

from vaster.vtools import measure_running_time
from vaster.structure import DataBasic
from vaster.check_complete import check_sbid_compelte


logger = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        prog='SummariseCandidates', 
        description='''
        Summarise candidates from one SBID, combine them into one csv. 
        With option to upload it via rclone and/or send an email notification to an address. 
        The program will check if all 36 beams are exists, and it will until it passes checks. 
        You can get rid of check and force run the scripts by specify --force. 
        ''', 
        epilog='Example Usage: summarise_candidates 60314 -e <test@gmail.com> -rclone google', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
        )
    parser.add_argument('sbids', type=int, nargs='+', help='SBID in format of number')
    parser.add_argument('--dir', type=str, default='.', help='path where SBxxxx were stored')
    parser.add_argument('-p', '--port', type=int, default=8053, help='local host forward port')
    parser.add_argument('--host', type=str, default='localhost', help='launched base URL host for website. Could be ada.physics.usyd.edu.au for ada')
    parser.add_argument('-r', '--crossmatch-radius', type=float, default=20, help='source crossmatch radius, unit of arcsec')
    parser.add_argument('--query-radius', type=float, default=3, help='query region radius, unit of degree')
    parser.add_argument('-e', '--email', default=None, nargs='+', 
                        help='send email notification to a list of email addresses')
    parser.add_argument('-rclone', default=None, 
                        help='rclone final csv to other cloud disk, e.g. google')
    parser.add_argument('--force', action='store_true', help='force run summary even if the sbid is not completed')
    parser.add_argument('-v', '--verbose', action='store_true', help='make it verbose')
    args = parser.parse_args()

    make_verbose(args)
    logger.debug(args)

    for i, sbid in enumerate(args.sbids):
        logger.info("Processing observation SB%s (%s/%s)", sbid, i+1, len(args.sbids))
        args.sbid = sbid

        databasic = DataBasic(args.sbid, args.dir)
        args.paths = databasic.paths
        logger.debug(args.paths)

        while True:
            complete = check_sbid_compelte(args, args.sbid, args.paths, databasic.nbeam)
            args.complete = complete
            if args.force:
                logger.info('Force running SB%s', args.sbid)
            if complete or args.force:
                logger.info('SB%s complete %s', args.sbid, complete)
                run(args, databasic)
                break
            else:
                time.sleep(600)


    end_time = time.time()
    measure_running_time(start_time, end_time)


def run(args, databasic):
    base_url = f"http://{args.host}:{args.port}/SB{args.sbid}/candidates/"
    dyspec_url = f"http://{args.host}:{args.port}/dyspec/"

    catfname = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 
                            'collections', 
                            'atnf_psrcat_good_astrometry.csv') # code location
    logger.info('Read catalogue from %s', catfname)
    catname, catsrc = read_catalog(args, catfname)

    cand_list = []
    for beamid in range(databasic.nbeam):
        fname = os.path.join(args.paths['path_cand'], f"SB{args.sbid}_beam{beamid:02d}_final.csv")
        if not os.path.exists(fname):
            print(fname, "doesn't exist. ")
            continue
        
        cands = Table.read(fname)
        cands = crossmatch_onebeam(args, cands, beamid, catname, catsrc, base_url, dyspec_url)
        cand_list.append(cands)

    if len(cand_list) == 0:
        logger.error(f'no final csv for SB{args.sbid}')
        sys.exit()

    final_cands = vstack(cand_list)
    logger.debug(final_cands)
    logger.info('High priority: %s', sum(final_cands['priority'] == 'high'))
    logger.info('Mid priority: %s', sum(final_cands['priority'] == 'mid'))
    logger.info('Low priority: %s', sum(final_cands['priority'] == 'low'))

    savename = os.path.join( args.paths['path_cand'], f'SB{args.sbid}.csv' )
    final_cands.write(savename, overwrite=True)
    logger.info(f'Saved final results to %s', savename)

    if args.email is not None:
        logger.warning('Will send email use variables set in $SMTP_USER and $SMTP_PWD')
        sender_email = os.getenv('SMTP_USER')
        sender_pwd = os.getenv('SMTP_PWD')
        body = write_email_body(args, final_cands, args.complete)

        for receiver_email in args.email:
            email_alert(body, sender_email, sender_pwd, receiver_email, args)
            logger.info('Yayy!!!')

    if args.rclone is not None:
        rclone_copy(args, savename)



def read_catalog(args, fname):
    cat = Table.read(fname)
    src = SkyCoord(cat['ra_deg'], cat['dec_deg'], unit=u.degree)
    return cat['PSRJ'], src


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
    


def crossmatch_onebeam(args, cands, beamid, catname, catsrc, base_url, dyspec_url):
    sbid = args.sbid
    beam = f'beam{beamid:02d}'

    beamsrc = SkyCoord(cands['beam_ra'], cands['beam_dec'], unit=u.degree)[0]

    # select pulsars with the primary beam 
    idx = beamsrc.separation(catsrc) < args.query_radius * u.degree
    sel_catname = catname[idx]
    sel_catsrc = catsrc[idx]

    new_cols = []
    
    for cand in cands:
        '''
        # check the priority 
        if cand['bright_sep_arcmin'] <= 5 or cand['deep_num'] >= 2:
            priority = 'low'
        elif cand['deep_sep_arcsec'] > 5 and cand['deep_sep_arcsec'] < 30:
            priority = 'low'
        elif cand['peak_map_sigma'] >= 5 and cand['deep_peak_flux'] < 0.01:
            priority = 'high'
        else:
            priority = 'mid'

        if cand['deep_sep_arcsec'] <= 2:
            priority = 'high'        
        '''

        # new priority metrics
        if cand['md_deep'] <= 0.1:
             priority = 'low'
        elif cand['bright_sep_arcmin'] <= 5:
            priority = 'mid'
        elif cand['deep_sep_arcsec'] < 2:
            priority = 'high'
        elif cand['deep_sep_arcsec'] > 20:
            priority = 'high'
        elif cand['deep_peak_flux'] < 0.002:
            priority = 'high'
        else:
            priority = 'mid'
            
        
        if priority == 'high':
            logger.info('High priority candidates: SB%s %s %s', sbid, beam, cand['name'])
        
        # get the location of plots
        lc = os.path.join(base_url, "SB{}_{}_lightcurve_{}.png".format(sbid, beam, cand['name']))
        dc = os.path.join(base_url, "SB{}_{}_deepcutout_{}.png".format(sbid, beam, cand['name']))
        sl = os.path.join(base_url, "SB{}_{}_slices_{}.gif".format(sbid, beam, cand['name']))
        map1 = os.path.join(base_url, "SB{}_{}_chisquare_map2.png".format(sbid, beam))
        map2 = os.path.join(base_url, "SB{}_{}_peak_map2.png".format(sbid, beam))
        dyspec = os.path.join(dyspec_url, "SB{}_{}_{}/".format(sbid, cand['name'], beam))

        candsrc = SkyCoord(cand['ra'], cand['dec'], unit=u.degree)
        ind = candsrc.separation(sel_catsrc) < args.crossmatch_radius * u.arcsec

        if sum(ind) == 0:
            match = ''
            sep = ''
        else:
            match = sel_catname[ind][0]
            sep = candsrc.separation(sel_catsrc)[ind].arcsec[0]
        
        new_cols.append([priority, lc, dc, sl, map1, map2, dyspec, beam, sbid, match, sep])
        
    logger.info('SB%s %s: total of %s candidates, with %s has a crossmatch', sbid, beam, len(cands), sum(ind))

    cands.add_columns(np.array(new_cols).T.tolist(), 
                      names=['priority', 'lightcurve', 'deepcutout', 'slices', 'chisq_map2', 'peak_map2', 'dyspec', 'beam', 'sbid', 'KNOWN_name', 'KNOWN_sep'])
    
    return cands


def write_email_body(args, cands, complete):

    if 'PSR_name' in cands.colnames:
        key = 'PSR_name'
    elif 'KNOWN_name' in cands.colnames:
        key = 'KNOWN_name'
    else:
        key = None

    if key is None:
        known = ''
    else:
        try:
            known = sum(~cands[key].mask)
        except:
            known = len(cands)

    body = 'Hi, ' + '\n\n' \
            f'SB{args.sbid} finished processing. Complete status **{complete}**. ' + '\n' + \
            f'CANDS={len(cands)} ' + f'KNOWN={known}' +  '\n\n' \
            'The csv should be uploded to google drive!'

    logger.debug('Email body: %s', body)

    return body


def email_alert(body, sender_email, sender_pwd, receiver_email, args):

    msg = EmailMessage()
    msg.set_content(body)
    msg['subject'] = f'VASTER SB{args.sbid} finished processing'
    msg['from'] = sender_email
    msg['to'] = receiver_email

    # creates SMTP session
    s = smtplib.SMTP('smtp.gmail.com', 587)
    # start TLS for security
    s.starttls()
    # Authentication
    s.login(sender_email, sender_pwd)
    # sending the mail
    # s.sendmail(sender_email, receiver_email, msgs)
    s.send_message(msg)
    # terminating the session
    s.quit()


def rclone_copy(args, fname):
    cmd = f'rclone copy {fname} {args.rclone}: -P'
    logger.info('Executing: %s', cmd)
    try:
        os.system(cmd)
    except:
        logger.error('Cannot executing %s', cmd)


if __name__ == '__main__':
    _main()
