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
3. Run slurm_MODDEEP_beam??.sh (deep modeling)
4. Run slurm_IMAFAST_beam??.sh (fast imaging, in either time interval)
5. Run slurm_SELCAND_beam??.sh (select candidates)
"""


from astroquery.utils.tap.core import TapPlus

import sys
import getpass
import requests
import xmltodict


sbid = sys.argv[-2]
savename = sys.argv[-1]


parallel = 10


# find visibilities from CASDA
tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
job = tap.launch_job_async("SELECT * FROM ivoa.obscore WHERE obs_id='{}' AND (dataproduct_type='{}' OR dataproduct_subtype='{}') ".format(sbid, 'visibility', 'catalogue.continuum.component'))

r = job.get_results()

vis = r[r['dataproduct_type'] == 'visibility']
cat = r[r['dataproduct_subtype'] == 'catalogue.continuum.component']

print('Found {} visibilities'.format(len(vis)))
print(r['filename'])


if len(vis) != 36:
    print("")
    print('Number of visibilities is not 36. Exit with error. ')
    print('You can download the visibility manually. ')
    sys.exit()
    

# get download url
username = 'wym20131028@gmail.com'
print('OPAL username:', username)
password = getpass.getpass(str("Enter OPAL password: "))


def get_url(access_url):
    session = requests.Session()
    session.auth = (username, password)
    
    s = session.get(access_url)
    data_dict = xmltodict.parse(s.content)
    
    return data_dict['VOTABLE']['RESOURCE'][0]['TABLE']['DATA']['TABLEDATA']['TR'][0]['TD'][1]
    
    
    
# and save it to a bash script 

with open(savename, 'w') as fw:
    fw.write("#!/bin/bash" + '\n')
    fw.write('\n')
    fw.write('# Generate automatically from a python script' + '\n')
    fw.write('# Test for all visibilities for SB{}'.format(sbid) + '\n')
    fw.write('# You can run this in the terminal directly, simply use "bash {}" '.format(savename) + '\n')
    fw.write('\n\n\n')
    

    for i, access_url in enumerate(r['access_url']):
        url = get_url(access_url)
        filename = r[i]['filename']
        
        print("Progress {}/{}...".format(i+1, len(r)))
        
        fw.write("echo " + '\n')
        fw.write("echo progress {}/{}".format(i+1, len(r)) + '\n')
        
        text = 'wget -O {} {} -t 0 -c'.format(filename, url)
        
        fw.write("echo " + text + '\n')
        
        if (i+1)%parallel == 0 or i == 35:
            fw.write(text + '\n')
        else:
            fw.write(text + ' &' + '\n')
            
        fw.write("sleep 1s" + '\n')
        fw.write('\n')


    
    
    
print('')
print('Generate {} finished. '.format(savename))








