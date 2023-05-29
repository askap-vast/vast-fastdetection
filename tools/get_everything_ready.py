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
4. Run slurm_IMAFAST_beam??.sh (fast imaging, in either time interval)
5. Run slurm_SELCAND_beam??.sh (select candidates)
"""


from astroquery.utils.tap.core import TapPlus

import os
import sys
import getpass
import requests
import xmltodict


sbid = sys.argv[-2]  # number only
loc = sys.argv[-1] # output location 


############################
# Find visibilities, selavy catalogues and fits images from CASDA
############################
tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
job = tap.launch_job_async("SELECT * FROM ivoa.obscore WHERE obs_id='{}' AND (dataproduct_type='{}' OR dataproduct_subtype='{}' OR dataproduct_subtype='{}') ".format(
    sbid, 'visibility', 'catalogue.continuum.component', 'cont.restored.t0'))

r = job.get_results()

vis = r[r['dataproduct_type'] == 'visibility']
cat = r[r['dataproduct_subtype'] == 'catalogue.continuum.component']
img = r[r['dataproduct_subtype'] == 'cont.restored.t0']

print('Found {} visibilities'.format(len(vis)))
print(vis['filename'])
print('')
print('Found {} selavy catalogues'.format(len(cat)))
print(cat['filename'])
print('')
print('Found {} images'.format(len(img)))
print(img['filename'])
print('')

if len(vis) != 36:
    print("")
    print('Number of visibilities is not 36. Exit with error. ')
    print('You can download the visibility manually. ')
    sys.exit()


############################
# Get download urls
############################
username = 'wym20131028@gmail.com'
print('OPAL username:', username)
password = getpass.getpass(str("Enter OPAL password: "))


def get_url(access_url):
    session = requests.Session()
    session.auth = (username, password)

    s = session.get(access_url)
    data_dict = xmltodict.parse(s.content)

    return data_dict['VOTABLE']['RESOURCE'][0]['TABLE']['DATA']['TABLEDATA']['TR'][0]['TD'][1]


############################
# Generate bash_GETDATA_beam??.sh
############################

for idx in range(36):
    savename = os.path.join(loc, 'bash_GETDATA_beam{:02d}.sh'.format(idx))

    with open(savename, 'w') as fw:
        fw.write("#!/bin/bash" + '\n')
        fw.write('\n')
        fw.write('# Generate automatically from a python script' + '\n')
        fw.write('# Download and untar visibility for SB{} beam{:02d}'.format(sbid, idx) + '\n')
        fw.write('# You can run this in terminal directly, simply use "bash {}" '.format(
            savename) + '\n')
        fw.write('\n\n\n')

        if 'beam{:02d}'.format(idx) not in vis[idx]['filename']:
            print('WARNING: no. {} -- beam number/order might be wrong. Continue running...'.format(idx))

        url = get_url(vis[idx]['access_url'])
        filename = vis[idx]['filename']

        text = 'wget -O {} {} -t 0'.format(filename, url)
        fw.write("echo " + '\n')
        fw.write(text + '\n')
        fw.write(text + ' -c' + '\n')
        fw.write('\n')

        text = 'tar xvf {}'.format(filename)
        fw.write("echo " + '\n')
        fw.write(text + '\n')
        fw.write('\n')

    print('Generate {} finished. '.format(savename))


############################
# Generate download_selavy.sh
############################
savename = os.path.join(loc, 'download_selavy.sh')

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

        text = 'wget -O {} {} -t 0'.format(filename, url)

        fw.write("echo " + '\n')
        fw.write("echo progress {}/{}".format(i+1, len(cat)) + '\n')
        fw.write("echo " + text + '\n')
        fw.write(text + '\n')
        fw.write("sleep 1s" + '\n')
        fw.write('\n')

print('Generate {} finished. '.format(savename))


############################
# Generate download_mosaic_images.sh
############################
savename = os.path.join(loc, 'download_mosaic_images.sh')

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

        text = 'wget -O {} {} -t 0'.format(filename, url)

        fw.write("echo " + '\n')
        fw.write("echo progress {}/{}".format(i+1, len(img)) + '\n')
        fw.write("echo " + text + '\n')
        fw.write(text + '\n')
        fw.write("sleep 1s" + '\n')
        fw.write('\n')

print('Generate {} finished. '.format(savename))



