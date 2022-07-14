#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 20:57:50 2022

@author: ywan3191
"""

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs.utils import skycoord_to_pixel
from astropy.nddata.utils import Cutout2D
from astropy.table import Table
import matplotlib.animation as animation

import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)



def plot_slices(src_name, imagelist, radius=8, vsigma=1, name='animation'):
    """Generate a gif contains a series of images cutout in given position. 
    
    src: 
        astropy Skycoord object. 
    imagelist: 
        A list of short FITS images, in correct order. 
    radius:
        unit of arcmin, the cutout size
    vsigma:
        float, the color range, in the unit of sigma
    """
    
    # get the source position 
    src = SkyCoord(src_name, unit=(u.hourangle, u.degree))
    logger.info("Get source position...")
    logger.info(src)
    
    # generate a plot, including a bunch of single images
    fig = plt.figure()
    ims = []
    
    # generate each frame
    for i, image in enumerate(imagelist):
        logger.info("Processing image number %s/%s" % (i, len(imagelist)))
        
        # open the fits
        hdu = fits.open(image)
        # get the fits data (drop-out extra dimensions)
        data = hdu[0].data.squeeze()
        # get wcs frame
        wcs = WCS(header=hdu[0].header, naxis=2)
        # generate cutout
        cutout = Cutout2D(data, 
                          position=src, 
                          size=radius*u.arcmin, 
                          wcs=wcs)
        
        # get the image rms 
        rms = np.nanstd(cutout.data) * 1e3
        # (for vmax and vmin)
        vmax = vsigma * rms
        vmin = -vsigma * rms
        
        # get the src in pixel (cutout image)
        xw, yw = cutout.input_position_cutout
        xw, yw = round(xw), round(yw)
        # get the image flux at src 
        flux = cutout.data[yw, xw] * 1e3
        
        # set the image frame
        if i == 0:
            fig.gca(projection=cutout.wcs)
            fig.gca().coords[0].set_major_formatter('hh:mm:ss')
            fig.gca().coords[1].set_major_formatter('dd:mm:ss')
            
        im = plt.imshow(cutout.data*1e3, 
                        cmap='seismic', 
                        origin='lower', 
                        vmin=vmin, 
                        vmax=vmax)
        
        # set the marker
        if cutout.data[yw, xw] > 0: # half of the range
            maker_color = 'black'
        else:
            maker_color = 'lightgray'

        sc = plt.scatter(src.ra.deg, 
                         src.dec.deg, 
                         marker='+', 
                         c=maker_color,
                         label=r'{:.2f} $\pm$ {:.2f} mJy'.format(flux, rms),
                         transform=fig.gca().get_transform('fk5')
                         )
        te = plt.text(x=0.60, y=0.95,
                      s=r'{:.2f} $\pm$ {:.2f} mJy {}'.format(flux, rms, i),
                      backgroundcolor='w',
                      transform=fig.gca().transAxes)

        ims.append([im, sc, te])

    fig.gca().set_title(name)
    fig.gca().set_xlabel('RA (J2000)')
    fig.gca().set_ylabel('DEC (J2000)')
    cbar = fig.colorbar(im)
    cbar.set_label('Flux (mJy/beam)')

    ani = animation.ArtistAnimation(fig, ims, interval=200, 
                                    blit=True, repeat_delay=1e6)
    ani.save('{}.gif'.format(name), dpi=80, writer='imagemagick')
    
    logger.info("Save image {}.gif".format(name))












