#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

"""
Generate the beam, field, SBID file for one target 
"""

from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
from astropy import units as u

import numpy as np

import glob
import os
import sys


reglist = glob.glob("/o9000/ASKAP/VAST/fast_test/detection_results/footprint/*.reg")

name = sys.argv[-1]
outputname = "output"


# beam position

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
    
    return beam_position, 1.2*radius



# check single coordinates

SOURCE, SBID, FIELD, BEAM, SEP = [], [], [], [], []
name = [name]
unique_src = SkyCoord(name, unit=(u.hourangle, u.degree))

for beam_reg in reglist:

#    print(beam_reg)
    
    sbid = beam_reg.split('/')[-1].split('_')[0]
    field = beam_reg.split('/')[-1][len(sbid)+1:-4]
    
#    print(sbid, field)

    beam_position, radius = get_beam_radius(beam_reg)
    idxc, idxcatalog, d2d, d3d = unique_src.search_around_sky(beam_position, radius*u.deg)
    
    source = name * len(idxcatalog)
    beam = ['beam{:02d}'.format(i) for i in idxc]
    sep = list(d2d.degree)
    
    SOURCE += source
    SBID += [sbid] * len(source)
    FIELD += [field] * len(source)
    BEAM += beam
    SEP += sep
    
candidates = Table(data=[SOURCE, SBID, FIELD, BEAM, SEP], 
                   names=('SOURCE', 'SBID', 'FIELD', 'BEAM', 'SEP'))
print(len(candidates))

candidates.write(outputname+".csv")
