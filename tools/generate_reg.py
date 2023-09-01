#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 19:07:11 2022

@author: ywan3191
"""

'''
Generate beam footprint .reg file from 36 model images
'''


import os
import sys
import glob
import math

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u



folder = sys.argv[-2]
savename = sys.argv[-1]
num = 36


    
def get_beam_fwhm(imagename):
    '''
    Get beam coordinate and fwhm size for a certain fits (beam) file
    
    imagename: str 
        fits image of certain beam
    
    Return:
        ra: str in hourangle format, e.g. 20:50:34.450
        dec: str in degree format, e.g. -30:20:45.32
        fwhm: full width at half maximum, unit of degree (diameter)
    '''
    
    hdu = fits.open(imagename)[0]
    
    beam_center = SkyCoord(hdu.header['CRVAL1'], 
                           hdu.header['CRVAL2'], unit=u.degree)
    beam_ra = beam_center.ra.to_string(
        unit=u.hourangle, sep=':', precision=3, pad=True
        )
    beam_dec = beam_center.dec.to_string(
        unit=u.degree, sep=':', precision=2, alwayssign=True, pad=True
        )
    
    fwhm = 3e8 / hdu.header['CRVAL3']/12 * 180/math.pi # unit of degree
    
    return beam_ra, beam_dec, fwhm
    
    
    
    
# =====================
# Write in a reg file
# =====================


with open(savename+'.reg', 'w') as fw:
    fw.write('global color=white' + '\n')
    fw.write('fk5' + '\n')

    # for each line 
    for i in range(num):
        
        imagename = glob.glob(os.path.join(folder, "*beam{:02d}*.fits".format(i)))
        
        if len(imagename) == 0:
            print("ERROR: No fits file for beam{:02d}".format(i))
            # sys.exit()
            continue
            
        elif len(imagename) != 1:
            print('')
            print('Number of matched images for beam{:02d} is'.format(i), len(imagename))
            print(imagename[:3])
            
        imagename = imagename[0]
        print('Select', imagename)
        
        beam_ra, beam_dec, fwhm = get_beam_fwhm(imagename)
        
        ellipse = 'ellipse({},{},{:.2f}d,{:.2f}d,0d)'.format(beam_ra, beam_dec, fwhm/2, fwhm/2)
        text = 'text({},{}) # text='.format(beam_ra, beam_dec) + '{' + str(i) + '}'
        
        fw.write(ellipse + '\n')
        fw.write(text + '\n')
    
    
    
print('Finished')
    

    





