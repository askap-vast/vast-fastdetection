#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 19:13:32 2022

@author: ywan3191
"""

# combine csv and give priority


from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u

import numpy as np
import sys
import glob
import os


sbid = sys.argv[-1]

# base_folder = "/import/ada2/ywan3191/fast_pipeline/results/SB{}/".format(sbid)
# base_url = 'ada.physics.usyd.edu.au:1028/view/fast_pipeline/results/SB{}/'.format(sbid)

base_folder = "/import/ada2/ywan3191/fast_survey/SB{}/".format(sbid)
# base_folder = 'c:/Users/wym19/OneDrive/Melbourne/06 Fast Transients GP Survey/Results/SB{}/'.format(sbid)
base_url = 'ada.physics.usyd.edu.au:1029/fast_survey/SB{}/'.format(sbid)


cand_list = []


for beam in ['beam{:02d}'.format(i) for i in range(36)]:
    
    name = os.path.join(base_folder, "SB{}_{}_final.csv".format(sbid, beam))
    if not os.path.exists(name):
        print(name, "doesn't exist. ")
        continue
    
    
    cands = Table.read(name)
    
    new_cols = []
    
    for cand in cands:
        
        # check the priority 
        if cand['bright_sep_arcmin'] <= 5 or cand['deep_num'] >= 2:
            priority = 'low'
        elif cand['deep_sep_arcsec'] > 5 and cand['deep_sep_arcsec'] < 30:
            priority = 'low'
        elif cand['peak_map_sigma'] >= 5 and cand['deep_peak_flux'] < 0.01:
            priority = 'high'
        else:
            priority = 'mid'
        
        
        if priority == 'high':
            print(sbid, beam, cand['name'], priority)
        
        # get the location of plots
        lc = os.path.join(base_url, "SB{}_{}_lightcurve_{}.png".format(sbid, beam, cand['name']))
        dc = os.path.join(base_url, "SB{}_{}_deepcutout_{}.png".format(sbid, beam, cand['name']))
        sl = os.path.join(base_url, "SB{}_{}_slices_{}.gif".format(sbid, beam, cand['name']))
        
        new_cols.append([priority, lc, dc, sl, beam, sbid])
        
        
    cands.add_columns(np.array(new_cols).T.tolist(), 
                      names=['priority', 'lightcurve', 'deepcutout', 'slices', 'beam', 'sbid'])
    
    cand_list.append(cands)
    
    
    
new_csv = vstack(cand_list)
print('High priority:', sum(new_csv['priority'] == 'high'))
print('Mid priority:', sum(new_csv['priority'] == 'mid'))
print('Low priority', sum(new_csv['priority'] == 'low'))


new_csv.write(os.path.join(base_folder, "SB{}.csv").format(sbid), overwrite=True)


