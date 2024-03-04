#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 20:41:50 2023
@author: ywan3191
Get ready for various scripts. 
For each beam, we want: 
1. Query CASDA to get visibility download url + selavy catalogue download url. 
Generate **bash_GETDATA_beam??.sh** and run the script to download visibilities 
and catalogues. (Double check the completion of visibilities downloading and 
untar visibilities) 
2. Run slurm_FIXDATA_beam??.sh (including rescale and fix dir)
3. Run slurm_MODELING_beam??.sh (deep modeling)
4. Run slurm_IMGFAST_beam??.sh (fast imaging, in either time interval)
5. Run slurm_SELCAND_beam??.sh (select candidates)
"""


from astroquery.utils.tap.core import TapPlus

import os
import sys
import getpass
import requests
import xmltodict

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
sh = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
sh.setFormatter(formatter)
logger.handlers.clear()
logger.addHandler(sh)


sbid = sys.argv[-2]  # number only
path = sys.argv[-1] # output parent location 

loc = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) # code location

logger.info("Processing observation SB%s", sbid)
logger.info("Saving outptus to %s", path)
logger.info("Using code in %s", loc)

############################
# Build file saving system structure 
############################
path_data = os.path.join(path, 'data') # saving visibilities, selavy catalogues 
path_models = os.path.join(path, 'models') 
path_images = os.path.join(path, 'images')
path_cand = os.path.join(path, 'candidates')
path_fits = os.path.join(path, 'fitsfiles')

path_scripts = os.path.join(path, 'scripts')
path_logs = os.path.join(path, 'logfiles')


def create_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
        logger.info('Create new directory %s', dir)
    else:
        logger.warning('Directory %s exists.', dir)


create_dir(path)
create_dir(path_data)
create_dir(path_models)
create_dir(path_images)
create_dir(path_cand)
create_dir(path_scripts)
create_dir(path_logs)
create_dir(path_fits)


############################
# Find visibilities, selavy catalogues and fits images from CASDA
############################
tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
job = tap.launch_job_async("SELECT * FROM ivoa.obscore WHERE obs_id='{}' AND (dataproduct_type='{}' OR dataproduct_subtype='{}' OR dataproduct_subtype='{}') ".format(
    'ASKAP-'+str(sbid), 'visibility', 'catalogue.continuum.component', 'cont.restored.t0'))

r = job.get_results()

vis = r[r['dataproduct_type'] == 'visibility']
cat = r[r['dataproduct_subtype'] == 'catalogue.continuum.component']
img = r[r['dataproduct_subtype'] == 'cont.restored.t0']

vis.sort('filename')

logger.info('Found {} visibilities'.format(len(vis)))
logger.info(vis['filename'])
print('')
logger.info('Found {} selavy catalogues'.format(len(cat)))
logger.info(cat['filename'])
print('')
logger.info('Found {} images'.format(len(img)))
logger.warning(img['filename'])
print('')

if len(vis) != 36:
    print("")
    logger.warning('Number of visibilities is not 36. Exit with error. ')
    logger.warning('You can download the visibility manually. ')
    sys.exit()


############################
# Get download urls
############################
# username = 'wym20131028@gmail.com'
username = os.getenv('OPAL_USER')
logger.info('OPAL username: %s', username)
password = os.getenv('OPAL_PWD')
# username = input("Enter OPAL username: ")
# password = getpass.getpass(str("Enter OPAL password: "))


def get_url(access_url):
    session = requests.Session()
    session.auth = (username, password)

    s = session.get(access_url)
    data_dict = xmltodict.parse(s.content)

    return data_dict['VOTABLE']['RESOURCE'][0]['TABLE']['DATA']['TABLEDATA']['TR'][0]['TD'][1]



############################
# Generate download_selavy.sh
############################
savename = os.path.join(path_scripts, 'download_selavy.sh')

with open(savename, 'w') as fw:
    fw.write("#!/bin/bash" + '\n')
    fw.write('\n')
    fw.write('# Generate automatically from a python script' + '\n')
    fw.write('# Download selavy components catalogue for SB{}'.format(sbid) + '\n')
    fw.write('# You can run this in terminal directly, simply use "bash {}" '.format(
        savename) + '\n')
    fw.write('\n\n\n')

    for i, access_url in enumerate(cat['access_url']):
        url = get_url(access_url)
        filename = cat[i]['filename']
        path_file = os.path.join(path_data, filename)

        text = 'wget -O {} {} -t 0 -c'.format(path_file, url)

        fw.write("echo " + '\n')
        fw.write("echo progress {}/{}".format(i+1, len(cat)) + '\n')
        fw.write("echo " + text + '\n')
        fw.write(text + '\n')
        fw.write("sleep 1s" + '\n')
        fw.write('\n')

logger.info('Writing {}'.format(savename))


############################
# Generate download_mosaic_images.sh
############################
savename = os.path.join(path_scripts, 'download_mosaic_images.sh')

with open(savename, 'w') as fw:
    fw.write("#!/bin/bash" + '\n')
    fw.write('\n')
    fw.write('# Generate automatically from a python script' + '\n')
    fw.write('# Download mosaiced fits images for SB{}'.format(sbid) + '\n')
    fw.write('# You can run this in terminal directly, simply use "bash {}" '.format(
        savename) + '\n')
    fw.write('\n\n\n')

    for i, access_url in enumerate(img['access_url']):
        url = get_url(access_url)
        filename = img[i]['filename']
        path_file = os.path.join(path_data, filename)

        text = 'wget -O {} {} -t 0 -c'.format(path_file, url)

        fw.write("echo " + '\n')
        fw.write("echo progress {}/{}".format(i+1, len(img)) + '\n')
        fw.write("echo " + text + '\n')
        fw.write(text + '\n')
        fw.write("sleep 1s" + '\n')
        fw.write('\n')

logger.info('Writing {}'.format(savename))


# ###########################
# # Generate bash_CHECKDATA.sh
# ###########################

savename = os.path.join(path_scripts, 'bash_CHECKDATA.sh')

with open(savename, 'w') as fw:
    fw.write("#!/bin/bash" + '\n')
    fw.write('\n')
    fw.write('# Generate automatically from a python script' + '\n')
    fw.write('# Download visibility for SB{} '.format(sbid) + '\n')
    fw.write('# You can run this in terminal directly, simply use "bash {}" '.format(
        savename) + '\n')
    fw.write('\n\n\n')

    for idx in range(36):

        logger.info("now is running beam{:02d}".format(idx))

        if 'beam{:02d}'.format(idx) not in vis[idx]['filename']:
            logger.warning('no. {} -- beam number/order might be wrong. Continue running...'.format(idx))


        url = get_url(vis[idx]['access_url'])
        filename = vis[idx]['filename']
        path_file = os.path.join(path_data, filename)

        text = 'wget -O {} {} -t 0 -c'.format(path_file, url)
        fw.write("echo " + '\n')
        fw.write("echo download data from {} progress {}/{}".format(sbid, idx+1, len(vis)) + '\n')
        fw.write(text + '\n')
        fw.write('\n')

print('Writing {}'.format(savename))



# ############################
# # Generate scripts for one beam data
# ############################

for idx in range(36):
    savename = os.path.join(path_scripts, 'bash_PROCESSING_beam{:02d}.sh'.format(idx))

    with open(savename, 'w') as fw:
        fw.write("#!/bin/bash" + '\n')
        fw.write('\n')
        fw.write('# Generate automatically from a python script' + '\n')
        fw.write('# Prcessing data for SB{} beam{:02d}'.format(sbid, idx) + '\n')
        fw.write("# It will download calibrated visibilities from CASDA" + '\n')
        fw.write('# You can run this in terminal directly, simply use "bash {}" '.format(
            savename) + '\n')
        fw.write('\n\n\n')

        if 'beam{:02d}'.format(idx) not in vis[idx]['filename']:
            logger.warning('No. {} -- beam number/order might be wrong. Continue running...'.format(idx))


        affix = 'SB{}_beam{:02d}'.format(sbid, idx)

        # ======== Download calibrated visibilities data =============

        url = get_url(vis[idx]['access_url'])
        filename = vis[idx]['filename']
        path_file = os.path.join(path_data, filename)

        text = 'wget -O {} {} -t 0 -c'.format(path_file, url)
        fw.write("echo beam{:02d}: Download calibrated visibilities data".format(idx) + '\n')
        fw.write(text + '\n')
        fw.write('\n')

        # ======== Untar the data get measurement sets ============

        text = 'tar xvf {} -C {}'.format(path_file, path_data)
        fw.write("echo beam{:02d}: Untar the data to measurement sets".format(idx) + '\n')
        fw.write(text + '\n')
        fw.write('\n')

        # ======== Fix ms scaling and beam position ============

        filename = vis[idx]['filename'][:-4]
        path_file = os.path.join(path_data, filename)

        text = 'python {} {} {}'.format(
            os.path.join(loc, 'tools', 'askapsoft_rescale.py'), 
            path_file, path_file+'.corrected'
            )        
        fw.write("echo beam{:02d}: Fix the measurement sets flux scaling".format(idx) + '\n')
        fw.write(text + '\n')
        fw.write('\n')

        text = 'python {} {}'.format(
            os.path.join(loc, 'tools', 'fix_dir.py'), 
            path_file+'.corrected'
            )    
        fw.write("echo beam{:02d}: Fix the measurement sets pointing".format(idx) + '\n')
        fw.write(text + '\n')
        fw.write('\n')

        # ======== Running CASA for sky model creation and subtraction ===========

        filename = vis[idx]['filename'][:-4] + '.corrected'
        path_file = os.path.join(path_data, filename)

        text = 'casa --log2term --logfile {} --nogui -c {} {} {}'.format(
            os.path.join(path_logs, 'casa_MODELING_{}.log'.format(affix)), 
            os.path.join(loc, 'imaging', 'model_making.py'), 
            path_file, 
            affix
        )
        fw.write("echo beam{:02d}: Create sky model and subtract...".format(idx) + '\n')
        fw.write('cd ' + path_models + ' \n')
        fw.write(text + '\n')
        fw.write('\n')

        # ======= Running CASA for model-subtracted short images creation =========

        filename = vis[idx]['filename'][:-4] + '.corrected'
        path_file = os.path.join(path_data, filename)

        text = 'casa --log2term --logfile {} --nogui -c {} {} {} {}'.format(
            os.path.join(path_logs, 'casa_IMGFAST_{}.log'.format(affix)), 
            os.path.join(loc, 'imaging', 'short_imaging.py'), 
            path_file, 
            affix, 
            10
        )
        fw.write("echo beam{:02d}: Create model-subtracted short images...".format(idx) + '\n')
        fw.write('cd ' + path_images + ' \n')
        fw.write(text + '\n')
        fw.write('\n')

        # ====== Candidates selection ===========

        text = 'python {} {} {} {} {} {} {} > {} 2>&1'.format(
            os.path.join(loc, 'select_candidates.py'), # scripts
            os.path.join(path_models, affix+'.image.tt0.fits'), # deep image
            os.path.join(path_data, cat[0]['filename']), # selavy catalogue
            path_images, # short images location
            'beam{:02d}'.format(idx), # beam number
            path_cand, # output directory 
            affix, # affix
            os.path.join(path_logs, 'bash_SELCAND_{}.log'.format(affix)), # output log files 
        )
        fw.write("echo beam{:02d}: Select candidates...".format(idx) + '\n')
        fw.write(text + '\n')
        fw.write('\n')

        # ======= Clean work tree ============

        filename = vis[idx]['filename'][:-4]
        path_file = os.path.join(path_data, filename)

        fw.write("echo beam{:02d}: Clean intermediate products...".format(idx) + '\n')
        fw.write('rm {}.tar'.format(path_file) + '\n')
        fw.write('rm -r {}'.format(path_file) + '\n')
        fw.write('mv {} {}'.format(os.path.join(path_models, '*beam{:02d}*.fits'.format(idx)), path_fits) + '\n')
        fw.write('mv {} {}'.format(os.path.join(path_images, '*beam{:02d}*.fits'.format(idx)), path_fits) + '\n')
        fw.write('rm -r {}/*beam{:02d}*'.format(path_models, idx) + '\n')
        fw.write('rm -r {}/*beam{:02d}*'.format(path_images, idx) + '\n')
        fw.write('mv {} {}'.format(os.path.join(path_fits, '*beam{:02d}*image*fits'.format(idx)), path_models) + '\n')
        fw.write('mv {} {}'.format(os.path.join(path_fits, '*beam{:02d}*fits'.format(idx)), path_images) + '\n')

        fw.write('\n')

        fw.write('echo beam{:02d}: Finished!'.format(idx))


    logger.info('Writing {}'.format(savename))
