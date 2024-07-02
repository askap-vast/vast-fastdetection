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
import aplpy

from vaster.vastfast import plot

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def _main():
    parser = ArgumentParser(description='Generate beam footprint', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--plot', action='store_true', help='Generate footprint plot')
    parser.add_argument('-o', '--output', type=str, default='output_footprint', help='Recommend format is SBID_FIELDNAME')
    parser.add_argument('-i', '--mosaic', type=str, help='mosaic image')
    parser.add_argument(dest='folder', type=str, default='.', help='Give location of models fits files')
    args = parser.parse_args()

    write_reg(args.output, args.folder)

    if args.plot:
        beam_position, radius = get_beam_radius(args.output+'.reg')
        plot_footprint(args, beam_position, radius)
        

    print('Finished')



    
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
# Writing in a reg file
# =====================

def write_reg(savename, folder, num=36):

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



def get_beam_radius(beam_reg):
    '''Get beam centre coordinates and fwhm radius from ds9 reg file
    
    Return:
        beam_position: a list of SkyCoords
        radius: int, fwhm in unit of degree
    '''

    with open(beam_reg) as f:
        beam_position = []
        radius = []

        for line in f:
            if line[:7] != 'ellipse':
                continue
            ra, dec, width, height, _ = line[8:].strip().split(',')

            beam_position.append(ra + ' ' + dec)


    beam_position = SkyCoord(beam_position, unit=(u.hourangle, u.degree))
    radius = float(width[:-1])
    
    return beam_position, radius



def plot_footprint(args, beam_position, radius):

    fitsname, savename = args.mosaic, args.output, 
    
    # read the image
    f = aplpy.FITSFigure(fitsname, figsize=(8, 8), dpi=100)
    
    # fix the wcs dimension issue
    plot.fix_aplpy_fits(f)
    
    # choose a color map (and scale if you want) to plot
    f.show_colorscale(cmap='plasma')

    f.show_circles(xw=beam_position.ra.degree, yw=beam_position.dec.degree, 
                    radius=radius, coords_frame='world', color='white', alpha=1)

    for num, coord in enumerate(beam_position):
        f.add_label(x=coord.ra.degree, y=coord.dec.degree, text=num, color="white", size="x-large")

    f.set_title(savename)

    f.savefig(filename=savename+'_footprint.png', dpi=100)
    
    
    
        
if __name__ == '__main__':
    _main()    

    





