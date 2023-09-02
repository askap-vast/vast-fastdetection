#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 20:57:50 2022

@author: ywan3191
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator

import aplpy
from skimage.feature import peak_local_max
from scipy.stats import chi2, norm

from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, search_around_sky
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from astropy.nddata.utils import Cutout2D
from astropy.table import Table, vstack
from astropy.time import Time

import warnings
warnings.filterwarnings("ignore")

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
    logger.info("Get source position ({}, {})...".format(src.ra.degree, src.dec.degree))
    
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

    fig.gca().set_title(src_name)
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
    f = aplpy.FITSFigure(fitsname, figsize=(8, 8))
    
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
    f.savefig(filename=imagename+'.png', dpi=100)
    
    
    
    
def plot_cutout(src_name, fitsname, radius=5, name='cutout'):
    """Plot cutout png from deep image
    src_name: str
        in format of "Jxxxx-xxxx"
    fitsname: str
        the location of the deep fits image
    radius: float
        cutout size, in unit of arcmin
    vsigma: float
        plot color range in unit of sigma
    """
    
    # get the source position 
    logger.info("Plotting deep cutout source {}...".format(src_name))
    src = SkyCoord(src_name, unit=(u.hourangle, u.degree))
    logger.info("Get source position...")
    logger.info(src)
    
    # get deep image information 
    f = aplpy.FITSFigure(fitsname, figsize=(5, 5))
    
    # fix the wcs dimension issue
    fix_aplpy_fits(f)
    
    # change unit from Jy to mJy
    f._data = f._data * 1e3
    
    f.recenter(src.ra, src.dec, radius=radius/60)
    f.show_grayscale()
    
    f.add_colorbar()
    f.colorbar.set_axis_label_text("Flux Density (mJy/beam)")
    f.colorbar.set_axis_label_pad(1)
    
    f.show_circles(src.ra, src.dec, radius=30/60/60, ec='orange')
    f.show_markers(src.ra, src.dec, c='red', marker='+', alpha=0.2)

    f.set_title(src_name)
    
    f.savefig("{}.png".format(name))
    logger.info("Save image {}.png".format(name))
    
    
    
    
def extract_lightcurve(src_name, imagelist, deep_imagename, residual_imagename=None, 
                    ksize=99):
    """plot lightcurve for each source candidate 
        ksize: kernel size for local rms calculation (square radius)
    """
    logger.info("Plotting lightcurve source {}...".format(src_name))
    # get the source position
    src = SkyCoord(src_name, unit=(u.hourangle, u.degree))
    
    # get deep image pixel position xw, yw
    deep_image = fits.open(deep_imagename)
    deep_wcs = WCS(header=deep_image[0].header, naxis=2)
    
    deep_xw, deep_yw = skycoord_to_pixel(coords=src, wcs=deep_wcs) 
    deep_xw, deep_yw = int(np.round(deep_xw)), int(np.round(deep_yw))
    
    # get deep flux density 
    deep_flux = deep_image[0].data.squeeze()[deep_yw, deep_xw]
    
    peak_flux = []
    local_rms = []
    timestamp = []
    
    # for each short image 
    for i, image in enumerate(imagelist):
        
        # open the fits
        hdu = fits.open(image)
        # get the fits data (drop-out extra dimensions)
        data = hdu[0].data.squeeze()
        # get wcs frame
        wcs = WCS(header=hdu[0].header, naxis=2)
        
        xw, yw = skycoord_to_pixel(coords=src, wcs=wcs)
        xw, yw = int(np.round(xw)), int(np.round(yw))
        
        peak_flux.append(data[yw, xw])
        
        # get observing time
        timestamp.append(hdu[0].header['DATE-OBS'])
        
        # calculate local rms 
        try:
            rms = np.nanstd(data[yw-ksize: yw+ksize, xw-ksize: xw+ksize])
            local_rms.append(rms)
        except:
            local_rms.append(np.nan)
            
    # covert them to np.array
    peak_flux = np.array(peak_flux)
    local_rms = np.array(local_rms)
    logger.info("Get peak flux shape {}, local rms shape {}".format(
        peak_flux.shape, local_rms.shape))
    
    
    # modify the peak flux (+deep image - residual image)
    
    # get residual image (if there is)
    try:
        residual_image = fits.open(residual_imagename)
        wcs = WCS(residual_image[0].header, naxis=2)
        
        xw, yw = skycoord_to_pixel(coords=src, wcs=wcs)
        xw, yw = int(np.round(xw)), int(np.round(yw))
        
        residual_flux = residual_image[0].data.squeeze()[yw, xw]
        logger.info("Use residual image information {}".format(residual_imagename))
        
    except:
        residual_flux = np.nanmean(peak_flux)
        logger.warning("No available residual image, using mean(peak_flux) instead. ")
        
    logger.info("Residual flux {} Jy/beam".format(residual_flux))
        
    # final peak flux density
    peak_flux = peak_flux + deep_flux - residual_flux 
    
    return peak_flux, local_rms, timestamp
    
    
    
    

# plot
def plot_lightcurve(flux, times, rms, title='', name='lightcurve'):
    
    fig, ax = plt.subplots()

    times = Time(times)
    times.format = 'datetime64'
    
    flux, rms = flux*1e3, rms*1e3
    
    ax.errorbar(x=times.value, y=flux, yerr=rms, color='black', marker='.', alpha=0.6)
    
    date_form = mdates.DateFormatter("%Y-%b-%d/%H:%M")
    ax.xaxis.set_major_formatter(date_form)
    
    ax.set_xlabel('Time (UTC)')
    ax.set_ylabel('Peak Flux Density (mJy/beam)')
    ax.set_title(title)
    
    # ax.set_ylim(bottom=-0.1, top=1.2*np.nanmax(np.array(flux))) 
    
    # set a reasonable time interval
    time_length = ((times[-1] - times[0]).to_value('sec') + 900) / 3600
    
    if int(time_length/4) != 0:
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=int(time_length/4)))

    fig.autofmt_xdate(rotation=15)
    
    fig.savefig("{}.png".format(name), bbox_inches='tight')
    
    
    
    
    
    
    
    
def get_sigma_logspace(data, flux):
    # get log statistics
    logmean = np.nanmean(np.log10(data))
    logstd = np.nanstd(np.log10(data))
    
    # sigma level
    sigma = (np.log10(flux) - logmean) / logstd
    
    return sigma


def get_threshold_logspace(data, sigma=5):
    
    # get log statistics
    logmean = np.nanmean(np.log10(data))
    logstd = np.nanstd(np.log10(data))
    # get threshold using given sigma
    threshold = 10 ** (logmean+sigma*logstd)
    
    logger.info("Threshold log space rms = {}, mean = {}".format(logstd, logmean))
    logger.info('Threshold log space is {} sigma = {}'.format(sigma, threshold))

    return threshold


def get_threshold_peak(value=None, sigma=None, num=70):

    # from sigma to calculate the theoritical value/threshold
    if value is None:
        value = norm.ppf(np.power(norm.cdf(sigma), 1/num))
        logger.info("Calculating peak threshold for {} images...".format(num))
        logger.info("{} sigma threshold (peak) is {:.2f}".format(sigma, value))
        return value
    
    # from threshold/given value to calculate theoritical sigma 
    elif sigma is None:
        sigma = norm.ppf(norm.cdf(value)**num)
        return sigma
    

def get_threshold_chisquare(value=None, sigma=None, num=70):

    df = num - 1 # degree of freedom

    if value is None:
        value = chi2.ppf(norm.cdf(sigma), df) / df
        logger.info("Calculating chi2 threshold for {} images...".format(num))
        logger.info("{} sigma threshold (chi2) is {:.2f}".format(sigma, value))
        return value
    
    elif sigma is None:
        sigma = norm.ppf(chi2.cdf(value*df, df))
        return sigma





## combine two/three csv to one csv/vot
def combine_csv(namelist, radius=10, 
                tablename='final_cand', savevot=True):
    """
    tables: list
        a list of location of the chisquare/Gaussian/peak csv catalogue
    radius: float
        the crossmatch radius between two catalogues, unit of arcsec
    """
    
    tables = []
    
    for name in namelist:
        
        try:
            table = Table.read(name)
            logger.info('Successfully read {}'.format(name))
        except:
            logger.warning('Cannot read {}'.format(name))
            continue    
        
        if len(table) == 0:
            logger.warning('Empty Table {}'.format(name))
            continue
        else:
            tables.append(table)
            
            
    if len(tables) == 0:
        logger.warning('Empty final table. ')
        return None
    
    stacked_table = vstack(tables, join_type='exact')
    logger.info('Stacked table length {}'.format(len(stacked_table)))
    
    src = SkyCoord(stacked_table['ra'], stacked_table['dec'], unit=u.degree)
    
    # combine unique sources
    idx1, idx2, sep2d, _ = search_around_sky(src, src, 
                                             seplimit=radius*u.arcsec)
    idx = np.full((len(stacked_table), ), True)
    idx[np.unique(idx2[idx1 < idx2])] = False
    
    final_csv = stacked_table[idx]

    # save csv table
    final_csv.write("{}.csv".format(tablename), overwrite=True)
    logger.info("Save final csv {}".format(tablename))
    logger.info("Final table length {}".format(len(final_csv)))
    
    # save vot table
    if savevot:
        final_csv.write("{}.vot".format(tablename), 
                        table_id="candidates", 
                        format="votable", 
                        overwrite=True)
        logger.info("Save final vot {}".format(tablename))
    






class Candidates:
    """Generate the vot table for final candidates
    """
    
    def __init__(self, chisq_map, peak_map, std_map, gaussian_map='', num=70):
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
        # number of images
        self.num = num
        
        
        # peak map
        self.peak_map = fits.open(peak_map)[0].data.squeeze()
        # std map
        self.std_map = fits.open(std_map)[0].data.squeeze()
        # gaussian map
        if gaussian_map != '':
            logger.info("Open Gaussian map %s" % gaussian_map)
            self.gaussian_map = fits.open(gaussian_map)[0].data.squeeze()
        else:
            logger.info("No gaussian map is input")

        
        
        
    def local_max(self, min_distance=30, sigma=5, data=None):
        '''Find the local maxium of an image
        
        sigma: identify blobs above a specfic sigma threshold
        min_distance: pixel number of the minimal distance of two neighbours blobs
        num: number of short images 
        '''
        if data == None or data == 'chisquare':
            data = self.chisq_map
            threshold = get_threshold_chisquare(sigma=sigma, num=self.num)
            
        elif data == 'peak':
            data = self.peak_map
            threshold = get_threshold_peak(sigma=sigma, num=self.num)
            
        elif data == "gaussian":
            data = self.gaussian_map
            threshold = get_threshold_logspace(data, sigma=sigma)
        
        # # get threshold in log space
        # threshold = get_threshold_logspace(data, sigma=sigma)
        
        # find local maximum 
        xy = peak_local_max(data, min_distance=min_distance, 
                            threshold_abs=threshold)
        
        # get coordiantes, in pixel and world 
        yp, xp = xy[:, 0], xy[:, 1]
        
        self.xp = xp
        self.yp = yp
        
        self.cand_src = pixel_to_skycoord(xp, yp, wcs=self.wcs)
        logger.info("Selected {} candidates...".format(self.cand_src.shape[0]))


    # Depends if we preformed primary beam correction on short images or not 
    # if no primary beam correction previously, we should apply for below factor 
    def primary_correction(self, x0:float):
        '''x0 is the distance to the beam centre, unit of degree
        '''
        sigma = self.fwhm.to_value(u.degree) / (2*np.sqrt(2*np.log(2)))
        return norm.pdf(x0, scale=sigma) / norm.pdf(0, scale=sigma)
        
        
        
        
    def select_candidates(self, deepcatalogue, deepimage=None, 
                          sep=30, mdlim=0.05, extlim=1.5, beamlim=1.2, 
                          bright=0.05):
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
        bright: float
            flux density threshold of bright sources, unit of Jy 
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
        fc = self.primary_correction(x0=self.cand_src.separation(self.beam_center).degree)
        self.md = self.std_map[self.yp, self.xp] / self.deep_peak_flux[self.deepidx] / fc # primary beam correction factor
        # self.md = self.std_map[self.yp, self.xp] / self.deep_peak_flux[self.deepidx] 
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
        
        # check number of close deep conterparts within 30 arcsec 
        idx1, _, _, _ = search_around_sky(coords1=self.cand_src, 
                                 coords2=self.deep_src, 
                                 seplimit=sep*u.arcsec)
        unique_idx, unique_counts = np.unique(idx1, return_counts=True)
        num_deep_close = np.zeros_like(self.cand_src, dtype=np.int)
        num_deep_close[unique_idx] = unique_counts
        self.deep_num = num_deep_close
        
        # check separation with bright deep sources
        deep_bright = self.deep_src[self.deep_int_flux > bright]
        if sum(self.deep_int_flux > bright) == 0:
            self.bright_sep_arcmin = np.full_like(self.cand_src, 999)
        else:
            _, d2d, _ = self.cand_src.match_to_catalog_sky(deep_bright)
            self.bright_sep_arcmin = d2d.arcmin
        
        
    
    def save_csvtable(self, tablename="cand_catalogue", savevot=False):
        """Save selected candidates to a csv table
        
        tablename: str
            saved name
        savevot: bool
            if True, will also save a vot format table 
        """
        
        # if no final candidates - don't need to save csv/vot
        
        
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
        # t['chi_square_sigma'] = get_sigma_logspace(self.chisq_map, 
        #                                            np.array(t['chi_square']))
        t['chi_square_sigma'] = get_threshold_chisquare(value=np.array(t['chi_square']), 
                                                        num=self.num)
        
        # read the peak value at each pixel 
        t['peak_map'] = self.peak_map[self.yp, self.xp][self.final_idx]
        # calculate the sigma
        # t['peak_map_sigma'] = get_sigma_logspace(self.peak_map, 
        #                                          np.array(t['peak_map']))
        t['peak_map_sigma'] = get_threshold_peak(value=np.array(t['peak_map']), 
                                                 num=self.num)
        
        # if there's Gaussian map 
        if hasattr(self, 'gaussian_map'):
            t['gaussian_map'] = self.gaussian_map[self.yp, self.xp][self.final_idx]
            t['gaussian_map_sigma'] = get_sigma_logspace(self.gaussian_map, 
                                                     np.array(t['gaussian_map']))
        else:
            t['gaussian_map']  = [np.nan] * sum(self.final_idx)
            t['gaussian_map_sigma'] = [np.nan] * sum(self.final_idx)
            
        
        # read the peak value at each pixel 
        t['std_map'] = self.std_map[self.yp, self.xp][self.final_idx]
        
        # modulation index using deep source flux
        t['md_deep'] = self.md[self.final_idx]
        
        # separaion to nearest deep counterpart
        t['deep_sep_arcsec'] = self.d2d.arcsec[self.final_idx]
        
        # number of close deep sources
        t['deep_num'] = self.deep_num[self.final_idx]
        
        # separation to bright deep source
        t['bright_sep_arcmin'] = self.bright_sep_arcmin[self.final_idx]
        
        # separation to beam center
        t['beam_sep_deg'] = self.cand_src.separation(self.beam_center).degree[self.final_idx]
        
        # beam center coordinates
        t['beam_ra']= [self.beam_center.ra.deg] * sum(self.final_idx)
        t['beam_dec']= [self.beam_center.dec.deg] * sum(self.final_idx)
        
        # name of the nearest deep counterpart
        t['deep_name'] = np.array(self.deep_name)[self.deepidx][self.final_idx]
        
        # coordinates for nearest deep counterpart 
        t['deep_ra_deg'] = self.deep_src.ra.degree[self.deepidx][self.final_idx]
        t['deep_dec_deg'] = self.deep_src.dec.degree[self.deepidx][self.final_idx]
        
        # flux density of the nearest deep counterpart
        t['deep_peak_flux'] = self.deep_peak_flux[self.deepidx][self.final_idx]
        t['deep_int_flux'] = self.deep_int_flux[self.deepidx][self.final_idx]
        
        
        
        # save csv table
        t.write("{}.csv".format(tablename), overwrite=True)
        logger.info("Save csv {}".format(tablename))
        
        # save vot table
        if savevot and len(t) != 0:
            t.write("{}.vot".format(tablename), 
                    table_id="candidates", 
                    format="votable", 
                    overwrite=True)
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
        # logger.info(self.catalogue.info)
        
        ####################
        # for Aegean convention
        ####################
        # # deep catalogue sources coordinates
        # self.deep_src = SkyCoord(self.catalogue['ra'], 
        #                          self.catalogue['dec'], 
        #                          unit=u.degree)
        
        # # peak flux
        # self.deep_peak_flux = np.array(self.catalogue['peak_flux']) # unit of Jy
        
        # # integrated flux
        # self.deep_int_flux = np.array(self.catalogue['int_flux']) # unit of Jy
        
        # # get name
        # # for selavy just read the column 'col_component_name'
        # # for aegean you want to use following code to read from scratch 
        # self.deep_name = ['J' + \
        #      src.ra.to_string(unit=u.hourangle, sep="", precision=0, pad=True) + \
        #      src.dec.to_string(sep="", precision=0, alwayssign=True, pad=True)
        #      for src in self.deep_src
        #     ]

        ####################
        # for selavy convention
        ####################
        self.deep_src = SkyCoord(self.catalogue['col_ra_deg_cont'], 
                                 self.catalogue['col_dec_deg_cont'], 
                                 unit=u.degree)
        
        self.deep_peak_flux = np.array(self.catalogue['col_flux_peak']) / 1e3 # unit of Jy
        self.deep_int_flux = np.array(self.catalogue['col_flux_int']) / 1e3 # unit of Jy
        self.deep_name = self.catalogue['col_component_name']

        
        
    def plot_fits(self, fitsname, imagename='plot_fits'):
        """Plot the fits and marker with given src
        """
        
        # read the image
        f = aplpy.FITSFigure(fitsname, figsize=(8, 8))
        
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
        f.savefig(filename=imagename+'.png', dpi=100)






class Products:
    """Generate the vot table for final candidates
    """
    
    def __init__(self, final_csv, limit=-1):
        """
        final_csv: str
            location of the csv catalogue for (final) candidates
        limit: int
            only plot first limit number of candidates. -1 means plot all candidates
        """
        cat = Table.read(final_csv)
        if limit==-1: 
            self.final_csv = cat
        else:
            cat.sort(keys=['peak_map', 'chi_square'], reverse=True)
            self.final_csv = cat[:limit]


        self.cand_name = self.final_csv['name']
        self.cand_src = SkyCoord(self.final_csv['ra'], self.final_csv['dec'], 
                                  unit=u.degree)
        
        
        
    def generate_cutout(self, fitsname, radius=5, savename='output'):
        
        # run the cutout plot one by one 
        for src_name in self.cand_name:
            plot_cutout(src_name, fitsname, 
                        radius=5, 
                        name='{}_{}'.format(savename, src_name))
            
            
            
    def generate_slices(self, imagelist, radius=5, vsigma=5, savename='output'):
        
        # run the slices plot one by one
        for src_name in self.cand_name:
            plot_slices(src_name, imagelist, 
                        radius=5, vsigma=5, 
                        name='{}_{}'.format(savename, src_name))




    def generate_lightcurve(self, imagelist, deepname, savename='output', 
                            savecsv=True):
        
        peak_flux_table = Table()
        local_rms_table = Table()
        
        # run the lightucrve plot one by one
        for i, src_name in enumerate(self.cand_name):
            peak_flux, local_rms, timestamp = extract_lightcurve(src_name=src_name, 
                                                      imagelist=imagelist, 
                                                      deep_imagename=deepname)
            
            if i == 0:
                peak_flux_table['Time'] = timestamp
                local_rms_table['Time'] = timestamp
                
            
            peak_flux_table[src_name] = peak_flux
            local_rms_table[src_name] = local_rms
            
            # plot
            plot_lightcurve(flux=peak_flux, 
                            times=timestamp, 
                            rms=local_rms, 
                            title=src_name, 
                            name='{}_{}'.format(savename, src_name))
            
        # save catalgoue 
        if savecsv:
            # save csv table
            peak_flux_table.write("{}_peak_flux.csv".format(savename))
            local_rms_table.write("{}_local_rms.csv".format(savename))

            logger.info("Save csv {}".format(savename))
            
            
            
            
            
            
            
            

