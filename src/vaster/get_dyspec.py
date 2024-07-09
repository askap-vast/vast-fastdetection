'''
Get dynamic spectrum for known PSR sources in each SBID

Author: Yuanming Wang
Date: 2024-06-16
'''

import pandas as pd
import argparse
import os

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(
            format='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')


def _main():
    parser = argparse.ArgumentParser(description='Create dynamic spectrum for known PSR sources for each SBID'\
        'Example Usage: python ~/scripts/vast_fastdetection/tools/run_dyspec_sbid.py 62648' )
    parser.add_argument('sbid', type=str, help='SBID, in format of pure number, e.g. 62648')
    parser.add_argument('--path', type=str, default='.', help='path to store SBID outputs')
    parser.add_argument('--dry-run', action='store_true', help='dry run')
    args = parser.parse_args()

    logger.info(args)

    fname = f'{args.path}/SB{args.sbid}/candidates/SB{args.sbid}.csv'

    if os.path.isfile(fname):

        tb = read_table(fname)
        keyname = check_keyname(tb)
        logger.info('SB%s crossmatch source keyname: %s', args.sbid, keyname)

        real = select_real(tb, key=keyname)
        cand = remove_duplication(real, key=keyname)
        run_dyspec_command(cand, args)

        logger.info('Program finished and exit as 0. ')

    else:
        logger.warning('%s does not exists!!!', fname)


def check_keyname(tb):
    if 'PSR_name' in tb.columns:
        keyname = 'PSR_name'
    elif 'KNOWN_name' in tb.columns:
        keyname = 'KNOWN_name'
    else:
        keyname = None
        logger.error('Cannot find proper known source column name')
        sys.exit()

    return keyname


def read_table(filename):
    logger.info('Reading table %s', filename)
    tb = pd.read_csv(filename)
    logger.info('Table length %s', len(tb))
    logger.debug(tb)

    return tb


def select_real(tb, key='PSR_name'):
    real = tb[ ~tb[key].isna() ]
    logger.info('Real sources %s', len(real))
    logger.debug(real)
    return real

    
def remove_duplication(tb, key='PSR_name'):

    names = tb[key].unique()
    logger.info('Real sources unique %s', len(names))
    logger.info(names)

    cand = pd.DataFrame(columns=tb.columns)
    logger.debug(cand)

    for name in names:
        rows = tb[ tb[key] == name ]
        sep = rows['beam_sep_deg']
        logger.info('source %s separation to beam centre is', name)
        logger.info(list(sep))
        ind = sep.argmin()
        logger.info('cloest index %s', ind)

        cand.loc[len(cand)] = rows.iloc[ind]
        
    logger.debug(cand)

    return cand


def run_dyspec_command(cand, args):
    '''
    bash ~/dyspectra/run_dyspec_bash.sh J083720.17-460248.44 SB45577 beam26
    '''
    namelist = cand['name']
    sbidlist = cand['sbid']
    beamlist = cand['beam']

    for name, sbid, beam in zip(namelist, sbidlist, beamlist):
        cmdtxt = f'bash ~/scripts/dyspectra/run_dyspec_bash.sh {name} SB{sbid} {beam}'
        logger.info(cmdtxt)

        if args.dry_run:
            logger.info('Dry run! Skip running... ')
        else:
            os.system(cmdtxt)
            os.system('rm *.last')

    logger.info('run_dyspec_command finished!')



if __name__ == "__main__":
    _main()
    








