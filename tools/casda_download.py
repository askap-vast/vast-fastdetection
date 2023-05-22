# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 18:09:10 2022

@author: ywan3191
"""

from astroquery.utils.tap.core import TapPlus

import sys
import getpass
import requests
import xmltodict


sbid = sys.argv[-2]
savename = sys.argv[-1]

#savename= 'visibility_download.sh'
datatype = 'visibility'
parallel = 10


# find visibilities from CASDA
tap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")
job = tap.launch_job_async("SELECT * FROM ivoa.obscore WHERE (obs_id='{}' AND dataproduct_type='{}')".format(sbid, datatype))

r = job.get_results()

print('Found {} visibilities'.format(len(r)))
print(r['filename'])


if len(r) != 36:
    print("")
    print('Number of visibilities is not 36. Exit with error. ')
    print('You can download the visibility manually. ')
    sys.exit()
    

# get download url
username = 'wym20131028@gmail.com'
password = getpass.getpass(str("Enter your OPAL password: "))


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







