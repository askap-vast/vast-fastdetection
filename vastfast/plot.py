#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 20:57:50 2022

@author: ywan3191
"""


import numpy as np
import matplotlib.pyplot as plt

import aplpy
from skimage.feature import peak_local_max

from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from astropy.nddata.utils import Cutout2D
from astropy.table import Table
import matplotlib.animation as animation

import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class ArgumentError(Exception):
    pass



def plot_slices(src_name, imagelist, radius=5, vsigma=5, name='animation'):
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
        # logger.info("Processing image number %s/%s" % (i, len(imagelist)))
        
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
        # if cutout.data[yw, xw] > 0: # half of the range
        #     marker_color = 'black'
        # else:
        #     marker_color = 'lightgray'
        marker_color = 'gray'

        sc = plt.scatter(src.ra.deg, 
                         src.dec.deg, 
                         marker='+', 
                         c=marker_color,
                         label=r'{:.2f} $\pm$ {:.2f} mJy'.format(flux, rms),
                         transform=fig.gca().get_transform('fk5')
                         )
        te = plt.text(x=0.60, y=0.95,
                      s=r'{:.2f} $\pm$ {:.2f} mJy ({})'.format(flux, rms, i),
                      backgroundcolor='w',
                      transform=fig.gca().transAxes)

        ims.append([im, sc, te])

    fig.gca().set_title(name)
    fig.gca().set_xlabel('RA (J2000)')
    fig.gca().set_ylabel('DEC (J2000)')
    cbar = fig.colorbar(im)
    cbar.set_label('Flux density (mJy/beam)')

    ani = animation.ArtistAnimation(fig, ims, interval=200, 
                                    blit=True, repeat_delay=1e6)
    ani.save('{}.gif'.format(name), dpi=80, writer='imagemagick')
    
    logger.info("Save image {}.gif".format(name))
    
    
    
    
def fix_aplpy_fits(aplpy_obj, dropaxis=2):
    """This removes the degenerated dimensions in APLpy 2.X...
    The input must be the object returned by aplpy.FITSFigure().
    `dropaxis` is the index where to start dropping the axis (by default it assumes the 3rd,4th place).
    """
    temp_wcs = aplpy_obj._wcs.dropaxis(dropaxis)
    temp_wcs = temp_wcs.dropaxis(dropaxis)
    aplpy_obj._wcs = temp_wcs


    # Then you can just do the following every time you load a FITS figure:
    # fig = aplpy.FITSFigure('fitsfilename')
    # fix_aplpy_fits(fig)


    
    
def plot_fits(fitsname, src=None, imagename='plot_fits'):
    """Plot the fits and marker with given src
    """

    # read the image
    f = aplpy.FITSFigure(fitsname, figsize=(10, 10), dpi=100)
    
    # fix the wcs dimension issue
    fix_aplpy_fits(f)
    
    # choose a color map (and scale if you want) to plot
    f.show_colorscale(cmap='plasma')
    
    
    if src != None:
        # show selected local maximum
        f.show_markers(xw=src.ra.deg, yw=src.dec.deg, coords_frame='world', 
                       marker='o', ec='cyan')
    
    # show their statistics
    # for i in range(sum(final_idx)):
    #     f.add_label(x=cand_src[final_idx][i].ra.degree, 
    #                 y=cand_src[final_idx][i].dec.degree, 
    #                 text='{:.1%}, {:.1f} mJy'.format(md[final_idx][i], 
    #                                                  np.array(catalogue[idx]['peak_flux'])[final_idx][i]*1e3), 
    #                 horizontalalignment='left', verticalalignment='bottom', color='white')
    
    # show the beam 
    fwhm = 3e8/f._header['CRVAL3']/12 * 180/np.pi
    f.show_circles(xw=f._header['CRVAL1'], yw=f._header['CRVAL2'], 
                   radius=fwhm/2, coords_frame='world', color='white')
    
    # show the selected circle
    f.show_circles(xw=f._header['CRVAL1'], yw=f._header['CRVAL2'], 
                   radius=1.2*fwhm/2, coords_frame='world', color='white', ls='--')
    
    
    # save image
    f.savefig(filename=imagename+'.png', dpi=300)
    
    
    
    
def get_sigma_logspace(data, flux):
    # get log statistics
    logmean = np.nanmean(np.log10(data))
    logstd = np.nanstd(np.log10(data))
    
    # sigma level
    sigma = (np.log10(flux) - logmean) / logstd
    
    return sigma


def get_threshold_logspace(data, sigma=3):
    
    # get log statistics
    logmean = np.nanmean(np.log10(data))
    logstd = np.nanstd(np.log10(data))
    # get threshold using given sigma
    threshold = 10 ** (logmean+sigma*logstd)
    
    logger.info("Threshold log space rms = {}, mean = {}".format(logstd, logmean))
    logger.info('Threshold log space is {} sigma = {}'.format(sigma, threshold))

    return threshold




class Candidates:
    """Generate the vot table for final candidates
    """
    
    def __init__(self, chisq_map, peak_map, std_map):
        """chisq_map: str
            chisquare map location, should be FITS file
        """
        
        if isinstance(chisq_map, str):
            try: 
                self.fi = fits.open(chisq_map)[0]
            except FileNotFoundError:
                logger.exception("Unable to open image %s" % chisq_map)
        elif isinstance(chisq_map, fits.HDUList):
            self.fi= chisq_map
        else:
            raise ArgumentError("Do not understand input image")
            
        # read chisquare map 
        self.chisq_map = self.fi.data.squeeze()
        # read basic information from chisquare map
        self.wcs = WCS(self.fi.header, naxis=2)
        # beam center
        self.beam_center = SkyCoord(self.fi.header['CRVAL1'], 
                                    self.fi.header['CRVAL2'], unit=u.degree)
        # FWHM of primary beam
        self.fwhm = 3e8/self.fi.header['CRVAL3']/12 * 180/np.pi * u.degree
        
        
        # peak map
        self.peak_map = fits.open(peak_map)[0].data.squeeze()
        # std map
        self.std_map = fits.open(std_map)[0].data.squeeze()

        
        
        
    def local_max(self, min_distance=30, sigma=5, data=None):
        '''Find the local maxium of an image
        
        sigma: identify blobs above a specfic sigma threshold
        min_distance: pixel number of the minimal distance of two neighbours blobs
        '''
        if data == None:
            data = self.chisq_map
        
        # get threshold in log space
        threshold = get_threshold_logspace(data, sigma=sigma)
        
        # find local maximum 
        xy = peak_local_max(data, min_distance=min_distance, 
                            threshold_abs=threshold)
        
        # get coordiantes, in pixel and world 
        yp, xp = xy[:, 0], xy[:, 1]
        
        self.xp = xp
        self.yp = yp
        
        self.cand_src = pixel_to_skycoord(xp, yp, wcs=self.wcs)
        logger.info("Selected {} candidates...".format(self.cand_src.shape[0]))
        
        
        
        
    def select_candidates(self, deepcatalogue, deepimage=None, 
                          sep=30, mdlim=0.05, extlim=1.5, beamlim=1.2):
        """Select high priority candidates using deep image information
        
        deepimage: str
            deep image name/location, should be FITS format
        deepcatalogue: str
            deep catalogue location, should be aegean format (vot)
        sep: float
            radius search for deep countparts, unit of arcsec 
        mdlim: float
            lower limit of modulation index to select candidates
        extlim: float
            upper limit to select compact sources
        """
        
        # read deep image 
        # self.read_deepfits(imagename=deepimage)
        # read deep catalogue
        self.read_catalogue(catalogue=deepcatalogue)
        
        # rule out candidates outside primary beam size
        self.beamlim = beamlim
        beamidx = self.cand_src.separation(self.beam_center) < beamlim*self.fwhm/2
        logger.info("Candidates inside {} primary beam: {}".format(beamlim, 
                                                                   sum(beamidx)))
        
        # find the deep catalpgue counterparts for each candidates
        self.deepidx, self.d2d, d3d = self.cand_src.match_to_catalog_sky(self.deep_src)
        
        # calculate the modulation index 
        self.md = self.std_map[self.yp, self.xp] / self.deep_peak_flux[self.deepidx]
        logger.info("Candidates with modulation index > {:.1%}: {}".format(
            mdlim, sum(self.md > mdlim)
            ))
        
        # calculate the extend feature
        ext = (self.deep_int_flux / self.deep_peak_flux)[self.deepidx]
        logger.info("Candidates with compactness < {:.2f}: {}".format(
            extlim, sum(ext < extlim)
            ))
        
        # only select candidates that
        # 1. within the primary beam size
        # 2. have no countparts in deep image
        # 3. or have a countpart with md > 0.05 and ext < 1.5
        self.final_idx = beamidx & ((self.d2d.arcsec > sep) | ((self.md > mdlim) & (ext < extlim)))
        logger.info("Final candidates: {}".format(sum(self.final_idx)))
        
        
    
    def save_csvtable(self, tablename="cand_catalogue", savevot=False):
        """Save selected candidates to a csv table
        
        tablename: str
            saved name
        savevot: bool
            if True, will also save a vot format table 
        """
        t = Table()
        
        # source id 
        t['source_id'] = np.arange(sum(self.final_idx))
        
        # name in J005800.91-235449.00 format
        self.cand_name = ['J' + \
             self.cand_src[self.final_idx][i].ra.to_string(unit=u.hourangle, 
                                                           sep="", 
                                                           precision=2, 
                                                           pad=True) + \
             self.cand_src[self.final_idx][i].dec.to_string(sep="", 
                                                            precision=2, 
                                                            alwayssign=True, 
                                                            pad=True)
             for i in range(sum(self.final_idx))
            ]
            
        t['name'] = self.cand_name
        
        
        # ra_str in format of 00:58:00.91 
        t['ra_str'] = self.cand_src[self.final_idx].ra.to_string(
            unit=u.hourangle, sep=':', precision=2, pad=True
            )
        
        # dec_str in format of -58:21:04.3
        t['dec_str'] = self.cand_src[self.final_idx].dec.to_string(
            unit=u.degree, sep=':', precision=2, 
            alwayssign=True, pad=True
            )
        
        # ra and dec in deg
        t['ra'] = self.cand_src[self.final_idx].ra.degree
        t['dec'] = self.cand_src[self.final_idx].dec.degree
        
        # read the chisq value at each pixel
        t['chi_square'] = self.chisq_map[self.yp, self.xp][self.final_idx]
        # calculate the sigma
        t['chi_square_sigma'] = get_sigma_logspace(self.chisq_map, 
                                                   np.array(t['chi_square']))
        
        # read the peak value at each pixel 
        t['peak_map'] = self.peak_map[self.yp, self.xp][self.final_idx]
        # calculate the sigma
        t['peak_map_sigma'] = get_sigma_logspace(self.peak_map, 
                                                 np.array(t['peak_map']))
        
        # read the peak value at each pixel 
        t['std_map'] = self.std_map[self.yp, self.xp][self.final_idx]
        
        # modulation index using deep source flux
        t['md_deep'] = self.md[self.final_idx]
        
        
        # save csv table
        t.write("{}.csv".format(tablename))
        logger.info("Save csv {}".format(tablename))
        
        # save vot table
        if savevot:
            t.write("{}.vot".format(tablename), table_id="candidates", format="votable")
            logger.info("Save vot {}".format(tablename))
        
        
        
        
        
    # def read_fits(self, imagename):
        
    #     if isinstance(imagename, str):
    #         try: 
    #             self.fi = fits.open(imagename)[0]
    #         except FileNotFoundError:
    #             logger.exception("Unable to open image %s" % imagename)
    #     elif isinstance(imagename, fits.HDUList):
    #         self.fi = imagename
    #     else:
    #         raise ArgumentError("Do not understand input image")
             
    #     # beam center
    #     self.beam_center = SkyCoord(self.fi.header['CRVAL1'], 
    #                                 self.fi.header['CRVAL2'], unit=u.degree)
             
    #     # FWHM of primary beam
    #     self.fwhm = 3e8/self.fi.header['CRVAL3']/12 * 180/np.pi * u.degree
        
        
        
        
    def read_catalogue(self, catalogue, tabletype="aegean"):
        """Read the deep image catalogue
        
        
        tabletype: "aegean" or "selavy"
            current only support aegean
        """
        
        logger.info("Read deep catalogue {}...".format(catalogue))
        self.catalogue = Table.read(catalogue)
        logger.info(self.catalogue.info)
        
        
        # for Aegean convention
        # deep catalogue sources coordinates
        self.deep_src = SkyCoord(self.catalogue['ra'], 
                                 self.catalogue['dec'], 
                                 unit=u.degree)
        
        # peak flux
        self.deep_peak_flux = np.array(self.catalogue['peak_flux'])
        
        # integrated flux
        self.deep_int_flux = np.array(self.catalogue['int_flux'])
        
        
        
    def plot_fits(self, fitsname, imagename='plot_fits'):
        """Plot the fits and marker with given src
        """
        
        # read the image
        f = aplpy.FITSFigure(fitsname, figsize=(10, 10), dpi=100)
        
        # fix the wcs dimension issue
        fix_aplpy_fits(f)
        
        # choose a color map (and scale if you want) to plot
        f.show_colorscale(cmap='plasma')
        
        
        if hasattr(self, 'final_idx'):
            # show selected local maximum
            f.show_markers(xw=self.cand_src[self.final_idx].ra.deg, 
                           yw=self.cand_src[self.final_idx].dec.deg, 
                           coords_frame='world', 
                           marker='o', 
                           ec='cyan')
            # show info 
            for i in range(sum(self.final_idx)):
                f.add_label(x=self.cand_src[self.final_idx][i].ra.degree, 
                            y=self.cand_src[self.final_idx][i].dec.degree, 
                            text='{:.1%}, {:.1f} mJy'.format(
                                self.md[self.final_idx][i], 
                                self.deep_peak_flux[self.deepidx][self.final_idx][i]*1e3), 
                                horizontalalignment='left', 
                                verticalalignment='bottom', 
                                color='white')
        
        elif hasattr(self, 'cand_src'):
                # show selected local maximum
                f.show_markers(xw=self.cand_src.ra.deg, 
                               yw=self.cand_src.dec.deg, 
                               coords_frame='world', 
                               marker='o', 
                               ec='cyan')
        
        
        # show the beam 
        fwhm = 3e8/f._header['CRVAL3']/12 * 180/np.pi
        f.show_circles(xw=f._header['CRVAL1'], yw=f._header['CRVAL2'], 
                       radius=fwhm/2, coords_frame='world', color='white')
        
        if hasattr(self, 'beamlim'):
        # show the selected circle
            f.show_circles(xw=f._header['CRVAL1'], 
                           yw=f._header['CRVAL2'], 
                           radius=self.beamlim*fwhm/2, 
                           coords_frame='world', 
                           color='white', 
                           ls='--')
        
        
        # save image
        f.savefig(filename=imagename+'.png', dpi=300)









