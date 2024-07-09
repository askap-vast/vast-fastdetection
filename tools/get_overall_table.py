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

# base_folder = "/import/ada2/ywan3191/fast_survey/SB{}/".format(sbid)
# base_folder = 'c:/Users/wym19/OneDrive/Melbourne/06 Fast Transients GP Survey/Results/SB{}/'.format(sbid)
# base_url = 'ada.physics.usyd.edu.au:1029/fast_survey/SB{}/'.format(sbid)

# base_folder = "/o9000/ASKAP/VAST/fast_survey/SB{}/candidates/".format(sbid)


base_folder = os.path.join(os.getcwd(), f'SB{sbid}', 'candidates')
base_url = "http://localhost:8053/SB{}/candidates/".format(sbid)
dyspec_url = "http://localhost:8053/dyspec/"
#SB62423_J180513.64-212011.36_beam31/

print(base_folder)

cand_list = []
psrcat_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'collections', 'atnf_psrcat_good_astrometry.csv') # code location
goodcat = Table.read(psrcat_path)
psrsrc = SkyCoord(goodcat['ra_deg'], goodcat['dec_deg'], unit=u.degree)


for beam in ['beam{:02d}'.format(i) for i in range(36)]:
    
    name = os.path.join(base_folder, "SB{}_{}_final.csv".format(sbid, beam))
    if not os.path.exists(name):
        print(name, "doesn't exist. ")
        continue
    
    cands = Table.read(name)
    beamsrc = SkyCoord(cands['beam_ra'], cands['beam_dec'], unit=u.degree)[0]

    # select pulsars with the primary beam 
    ind = beamsrc.separation(psrsrc) < 3*u.degree
    selpsrcat = goodcat[ind]
    selpsrsrc = psrsrc[ind]

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

        if cand['deep_sep_arcsec'] <= 2:
            priority = 'high'        
        
        if priority == 'high':
            print(sbid, beam, cand['name'], priority)
        
        # get the location of plots
        lc = os.path.join(base_url, "SB{}_{}_lightcurve_{}.png".format(sbid, beam, cand['name']))
        dc = os.path.join(base_url, "SB{}_{}_deepcutout_{}.png".format(sbid, beam, cand['name']))
        sl = os.path.join(base_url, "SB{}_{}_slices_{}.gif".format(sbid, beam, cand['name']))
        map1 = os.path.join(base_url, "SB{}_{}_chisquare_map2.png".format(sbid, beam))
        map2 = os.path.join(base_url, "SB{}_{}_peak_map2.png".format(sbid, beam))
        dyspec = os.path.join(dyspec_url, "SB{}_{}_{}/".format(sbid, cand['name'], beam))

        candsrc = SkyCoord(cand['ra'], cand['dec'], unit=u.degree)
        ind = candsrc.separation(selpsrsrc) < 20*u.arcsec

        if sum(ind) == 0:
            atnfmatch = ''
            atnfsep = ''
        else:
            atnfmatch = selpsrcat[ind]['PSRJ'][0]
            atnfsep = candsrc.separation(selpsrsrc)[ind].arcsec[0]
        
        new_cols.append([priority, lc, dc, sl, map1, map2, dyspec, beam, sbid, atnfmatch, atnfsep])
        
    print(beam+':', 'Total', len(cands), 'Pulsars', sum(ind))

    cands.add_columns(np.array(new_cols).T.tolist(), 
                      names=['priority', 'lightcurve', 'deepcutout', 'slices', 'chisq_map2', 'peak_map2', 'dyspec', 'beam', 'sbid', 'PSR_name', 'PSR_sep'])
    
    cand_list.append(cands)
    
if len(cand_list) == 0:
    print(f'no final csv for SB{sbid}')
    sys.exit()

    
    
new_csv = vstack(cand_list)
print('High priority:', sum(new_csv['priority'] == 'high'))
print('Mid priority:', sum(new_csv['priority'] == 'mid'))
print('Low priority', sum(new_csv['priority'] == 'low'))


new_csv.write(os.path.join(base_folder, "SB{}.csv").format(sbid), overwrite=True)
print(f'Saved final results to {os.path.join(base_folder, "SB{}.csv").format(sbid)}')

