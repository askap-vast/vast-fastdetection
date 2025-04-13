#!/usr/bin/env python
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

import logging
logger = logging.getLogger(__name__)


class ArgumentError(Exception):
    pass


def plot_slices(src_name, imagelist, radius=5, vsigma=5, name='animation'):
    """
    Generate an animated GIF showing cutouts of a source from a list of FITS images.

    Parameters
    ----------
    src_name : str
        Source name, interpreted by SkyCoord (e.g. 'Jhhmmss+ddmmss').
    imagelist : list of str
        Ordered list of FITS filenames.
    radius : float
        Cutout radius in arcminutes.
    vsigma : float
        Controls color scale as ±vsigma * rms (displayed using $\pm$ in plot labels).
    name : str
        Output filename (without extension).
    """
    src = SkyCoord(src_name, unit=(u.hourangle, u.degree))
    logger.info(f"Source position: RA={src.ra.hms}, DEC={src.dec.dms}")

    fig = plt.figure()
    ims = []

    for i, image_path in enumerate(imagelist):
        with fits.open(image_path) as hdul:
            data = hdul[0].data.squeeze().copy()
            wcs = WCS(hdul[0].header, naxis=2)

        cutout = Cutout2D(data, position=src, size=radius * u.arcmin, wcs=wcs)

        rms = np.nanstd(cutout.data) * 1e3
        vmin, vmax = -vsigma * rms, vsigma * rms
        flux = cutout.data[int(round(cutout.input_position_cutout[1])),
                           int(round(cutout.input_position_cutout[0]))] * 1e3

        if i == 0:
            ax = fig.gca(projection=cutout.wcs)
            ax.coords[0].set_major_formatter('hh:mm:ss')
            ax.coords[1].set_major_formatter('dd:mm:ss')
        else:
            ax = fig.gca()

        im = ax.imshow(cutout.data * 1e3, cmap='seismic', origin='lower',
                       vmin=vmin, vmax=vmax)

        sc = ax.scatter(src.ra.deg, src.dec.deg, marker='+', c='gray',
                        label=rf'${flux:.2f} \pm {rms:.2f}$ mJy',
                        transform=ax.get_transform('fk5'))

        te = ax.text(0.60, 0.95,
                     rf'${flux:.2f} \pm {rms:.2f}$ mJy ({i})',
                     backgroundcolor='w',
                     transform=ax.transAxes)

        ims.append([im, sc, te])

    ax.set_title(src_name)
    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('DEC (J2000)')

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Flux density (mJy/beam)')

    ani = animation.ArtistAnimation(fig, ims, interval=200, blit=True, repeat_delay=1e6)
    ani.save(f'{name}.gif', dpi=80, writer='imagemagick')

    plt.close('all')
    logger.info(f"Saved animation to {name}.gif")
    
    
def save_fits_cube(src_name, imagelist, radius=5, name='cube'):
    """
    Create a FITS cube from cutouts centered at a given source position.

    Parameters
    ----------
    src_name : str
        Source position in SkyCoord-compatible format.
    imagelist : list of str
        List of FITS image paths.
    radius : float
        Size of each cutout (in arcminutes).
    name : str
        Output cube filename prefix ('.fits' added automatically).
    """
    src = SkyCoord(src_name, unit=(u.hourangle, u.degree))
    logger.info(f"Creating FITS cube for source: RA={src.ra.hms}, DEC={src.dec.dms}")

    # Determine shape from the first image
    with fits.open(imagelist[0], memmap=True) as hdul:
        data = hdul[0].data.squeeze()
        wcs = WCS(hdul[0].header, naxis=2)
        cutout = Cutout2D(data, position=src, size=radius * u.arcmin, wcs=wcs)
        ny, nx = cutout.data.shape
        header = cutout.wcs.to_header()

    # Allocate output cube array
    n_frames = len(imagelist)
    cube_data = np.empty((n_frames, ny, nx), dtype=np.float32)

    for i, image_path in enumerate(imagelist):
        with fits.open(image_path, memmap=True) as hdul:
            data = hdul[0].data.squeeze()
            wcs = WCS(hdul[0].header, naxis=2)
            cutout = Cutout2D(data, position=src, size=radius * u.arcmin, wcs=wcs)
            cube_data[i] = cutout.data

    # Write FITS cube
    hdu = fits.PrimaryHDU(data=cube_data, header=header)
    output_path = f"{name}.fits"
    hdu.writeto(output_path, overwrite=True)
    logger.info(f"Saved FITS cube to {output_path}")

    
def save_fits_cutout(src_name, image_path, radius=5, name='cutout'):
    """
    Save a FITS cutout around a source from a single image.

    Parameters
    ----------
    src_name : str
        Source position string in SkyCoord-compatible format.
    image_path : str
        Input FITS filename.
    radius : float
        Cutout radius in arcminutes.
    name : str
        Output filename prefix ('.fits' will be appended).
    """
    src = SkyCoord(src_name, unit=(u.hourangle, u.degree))
    logger.info(f"Creating FITS cutout for source: RA={src.ra.hms}, DEC={src.dec.dms}")
    with fits.open(image_path) as hdul:
        data = hdul[0].data.squeeze().copy()
        wcs = WCS(hdul[0].header, naxis=2)
        cutout = Cutout2D(data, position=src, size=radius * u.arcmin, wcs=wcs)

        hdul[0].data = cutout.data
        hdul[0].header.update(cutout.wcs.to_header())

        output_path = f'{name}.fits'
        hdul.writeto(output_path, overwrite=True)
        logger.info(f"Saved FITS cutout to {output_path}")

        
def fix_aplpy_fits(aplpy_obj, dropaxis=2):
    """
    Fix dimensionality issues in APLpy 2.X by removing extra axes.

    Parameters
    ----------
    aplpy_obj : aplpy.FITSFigure
        The APLpy object to modify in-place.
    dropaxis : int
        Starting axis index to drop (default 2).
    """
    temp_wcs = aplpy_obj._wcs.dropaxis(dropaxis)
    temp_wcs = temp_wcs.dropaxis(dropaxis)
    aplpy_obj._wcs = temp_wcs

    # Then you can just do the following every time you load a FITS figure:
    # fig = aplpy.FITSFigure('fitsfilename')
    # fix_aplpy_fits(fig)
    
    
def plot_cutout(src_name, fitsname, radius=5, name='cutout'):
    """
    Create and save a PNG plot of a cutout centered at the specified source position.

    Parameters
    ----------
    src_name : str
        Source name (SkyCoord-readable format).
    fitsname : str
        Path to the deep image FITS file.
    radius : float
        Cutout radius in arcminutes.
    name : str
        Output filename prefix for the PNG image.
    """ 
    logger.info("Plotting deep cutout source {}...".format(src_name))
    src = SkyCoord(src_name, unit=(u.hourangle, u.degree))
    logger.info("Get source position...")
    logger.info(src)
    
    f = aplpy.FITSFigure(fitsname, figsize=(5, 5))
    fix_aplpy_fits(f)

    f._data *= 1e3  # Convert flux units from Jy to mJy
    f.recenter(src.ra, src.dec, radius=radius / 60)
    f.show_grayscale()

    f.add_colorbar()
    f.colorbar.set_axis_label_text("Flux Density (mJy/beam)")
    f.colorbar.set_axis_label_pad(1)

    f.show_circles(src.ra, src.dec, radius=30 / 3600, ec='orange')
    f.show_markers(src.ra, src.dec, marker='+', c='red', alpha=0.2)
    f.set_title(src_name)

    f.savefig(f"{name}.png")
    plt.close('all')
    logger.info(f"Saved cutout image to {name}.png")
    
    
def extract_lightcurve(src_name, imagelist, deep_imagename, residual_imagename=None, ksize=99):
    """
    Extract the lightcurve of a source from a series of short FITS images.

    Parameters
    ----------
    src_name : str
        Source name or coordinates (SkyCoord-compatible).
    imagelist : list of str
        List of FITS images in temporal order.
    deep_imagename : str
        Path to the deep image FITS file.
    residual_imagename : str, optional
        Path to the residual FITS image (for flux correction).
    ksize : int
        Half-size of the square region used for local RMS estimation.

    Returns
    -------
    tuple
        Arrays of peak_flux, local_rms, and timestamps.
    """
    logger.info(f"Extracting lightcurve for source {src_name}...")
    src = SkyCoord(src_name, unit=(u.hourangle, u.degree))

    with fits.open(deep_imagename) as hdul:
        deep_data = hdul[0].data.squeeze().copy()
        deep_wcs = WCS(header=hdul[0].header, naxis=2)
        xw, yw = map(int, np.round(skycoord_to_pixel(src, deep_wcs)))
        deep_flux = deep_data[yw, xw]

    peak_flux = []
    local_rms = []
    timestamp = []

    for image_path in imagelist:
        with fits.open(image_path) as hdul:
            data = hdul[0].data.squeeze().copy()
            wcs = WCS(hdul[0].header, naxis=2)
            xw, yw = map(int, np.round(skycoord_to_pixel(src, wcs)))

            peak_flux.append(data[yw, xw])
            timestamp.append(hdul[0].header.get('DATE-OBS', 'NaT'))

            try:
                window = data[max(0, yw - ksize): yw + ksize, max(0, xw - ksize): xw + ksize]
                local_rms.append(np.nanstd(window))
            except Exception:
                local_rms.append(np.nan)

    peak_flux = np.array(peak_flux)
    local_rms = np.array(local_rms)

    # Estimate residual flux correction
    if residual_imagename:
        try:
            with fits.open(residual_imagename) as hdul:
                residual_data = hdul[0].data.squeeze().copy()
                wcs = WCS(hdul[0].header, naxis=2)
                xw, yw = map(int, np.round(skycoord_to_pixel(src, wcs)))
                residual_flux = residual_data[yw, xw]
            logger.info(f"Using residual image: {residual_imagename}")
        except Exception:
            residual_flux = np.nanmean(peak_flux)
            logger.warning("Residual image not found or error encountered. Using mean peak flux as fallback.")
    else:
        residual_flux = np.nanmean(peak_flux)
        logger.warning("No residual image provided. Using mean peak flux as fallback.")

    logger.info(f"Residual flux: {residual_flux:.3e} Jy/beam")

    corrected_flux = peak_flux + deep_flux - residual_flux
    return corrected_flux, local_rms, timestamp
    
    
def plot_lightcurve(flux, times, rms, title='', name='lightcurve'):
    """
    Plot a lightcurve with error bars and save it to a PNG file.

    Parameters
    ----------
    flux : array-like
        Peak flux values (in Jy).
    times : array-like
        Timestamps (ISO format or Time-compatible).
    rms : array-like
        Local RMS values corresponding to the flux measurements (in Jy).
    title : str
        Plot title.
    name : str
        Output filename prefix for the PNG image.
    """
    fig, ax = plt.subplots()

    times = Time(times)
    times.format = 'datetime64'

    flux_mJy = np.array(flux) * 1e3
    rms_mJy = np.array(rms) * 1e3

    ax.errorbar(x=times.value, y=flux_mJy, yerr=rms_mJy,
                color='black', marker='.', alpha=0.6)

    ax.set_xlabel('Time (UTC)')
    ax.set_ylabel('Peak Flux Density (mJy/beam)')
    ax.set_title(title)

    date_form = mdates.DateFormatter("%Y-%b-%d/%H:%M")
    ax.xaxis.set_major_formatter(date_form)

    # Choose tick spacing based on time span
    if len(times) > 1:
        span_sec = (times[-1] - times[0]).to_value('sec') + 900
        span_hr = span_sec / 3600
        interval = int(span_hr / 4)
        if interval > 0:
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=interval))

    fig.autofmt_xdate(rotation=15)
    fig.savefig(f"{name}.png", bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved lightcurve to {name}.png")


def combine_csv(namelist, radius=10, tablename='final_cand', savevot=True):
    """
    Combine multiple source candidate catalogs into one table with crossmatching.

    Parameters
    ----------
    namelist : list of str
        List of paths to input CSV/VOT tables.
    radius : float
        Maximum angular separation (in arcsec) for crossmatch merging.
    tablename : str
        Output table filename prefix (CSV and optionally VOT).
    savevot : bool
        Whether to also save the combined result as a VOTable.

    Returns
    -------
    astropy.table.Table or None
        Combined table of unique candidates, or None if input tables are empty.
    """
    tables = []
    
    for name in namelist:
        try:
            table = Table.read(name)
            if len(table) == 0:
                logger.warning(f"Empty table: {name}")
                continue
            tables.append(table)
            logger.info(f"Read: {name}")
        except Exception:
            logger.warning(f"Could not read: {name}")
            
    if len(tables) == 0:
        logger.warning('Empty final table. ')
        return None
    
    stacked_table = vstack(tables, join_type='exact')
    logger.info('Stacked table length {}'.format(len(stacked_table)))
    
    src = SkyCoord(stacked_table['ra'], stacked_table['dec'], unit=u.degree)
    
    # Crossmatch to identify duplicates
    idx1, idx2, _, _ = search_around_sky(src, src, seplimit=radius * u.arcsec)
    idx = np.full((len(stacked_table), ), True)
    idx[np.unique(idx2[idx1 < idx2])] = False
    
    final_table = stacked_table[idx]

    # save csv table
    final_table.write("{}.csv".format(tablename), overwrite=True)
    logger.info("Save final csv {}.csv".format(tablename))
    logger.info("Final table length {}".format(len(final_table)))
    
    # save vot table
    if savevot:
        final_table.write("{}.vot".format(tablename), table_id="candidates", 
                        format="votable", overwrite=True)
        logger.info("Save final vot {}.vot".format(tablename))

    return final_table
    

class Candidates:
    """
    Class to handle source candidate selection, evaluation, and table generation
    using chi-square, peak flux, Gaussian, and standard deviation maps.
    """
    def __init__(self, chisq_map, peak_map, std_map, gaussian_map='', num=70):
        """
        Initialize a Candidates instance with required input maps.

        Parameters
        ----------
        chisq_map : str or fits.HDUList
            Path to chi-square FITS image or already opened FITS HDU.
        peak_map : str
            Path to peak flux FITS image.
        std_map : str
            Path to standard deviation map (RMS).
        gaussian_map : str, optional
            Path to Gaussian map (if available).
        num : int
            Number of short images used in variability analysis.
        """
        if isinstance(chisq_map, str):
            try:
                with fits.open(chisq_map) as hdul:
                    self.chisq_map = hdul[0].data.squeeze().copy()
                    header = hdul[0].header
            except FileNotFoundError:
                logger.exception(f"Unable to open image {chisq_map}")
                raise
        elif isinstance(chisq_map, fits.HDUList):
            self.chisq_map = chisq_map[0].data.squeeze()
            header = chisq_map[0].header
        else:
            raise ArgumentError("Do not understand input image")

        self.wcs = WCS(header, naxis=2)
        self.beam_center = SkyCoord(header['CRVAL1'], header['CRVAL2'], unit=u.degree)
        self.fwhm = 3e8 / header['CRVAL3'] / 12 * 180 / np.pi * u.degree
        self.num = num
        
        self.peak_map = fits.getdata(peak_map).squeeze()
        self.std_map = fits.getdata(std_map).squeeze()

        if gaussian_map:
            logger.info(f"Opening Gaussian map: {gaussian_map}")
            self.gaussian_map = fits.getdata(gaussian_map).squeeze()
        else:
            logger.info("No Gaussian map provided.")


    def get_sigma_logspace(self, data, flux):
        logger.debug(data.shape)
        logmean = np.nanmean(np.log10(data))
        logstd = np.nanstd(np.log10(data))
        return (np.log10(flux) - logmean) / logstd


    def get_threshold_logspace(self, data, sigma=5):
        logmean = np.nanmean(np.log10(data))
        logstd = np.nanstd(np.log10(data))
        threshold = 10 ** (logmean + sigma * logstd)
        logger.info("Threshold log space rms = {}, mean = {}".format(logstd, logmean))
        logger.info('Threshold log space is {} sigma = {}'.format(sigma, threshold))
        return threshold


    def get_threshold_peak(self, value=None, sigma=None, num=70):
        if value is None:
            value = norm.isf(-norm.logcdf(sigma) / num)
            logger.info(f"{sigma} sigma peak threshold for {num} images = {value:.2f}")
            return value
        elif sigma is None:
            return norm.isf(-norm.logcdf(value) * num)


    def get_threshold_chisquare(self, value=None, sigma=None, num=70):
        df = num - 1
        if value is None:
            value = chi2.isf(norm.sf(sigma), df) / df
            logger.info(f"{sigma} sigma chi-square threshold (df={df}) = {value:.2f}")
            return value
        elif sigma is None:
            return norm.isf(chi2.sf(value * df, df))

    
    def local_max(self, min_distance=30, sigma=5, threshold=10, data=None):
        '''Identify local maxima in a given map above a sigma threshold.
        
        sigma: identify blobs above a specfic sigma threshold
        min_distance: pixel number of the minimal distance of two neighbours blobs
        num: number of short images 
        '''
        if data == None or data == 'chisquare':
            data = self.chisq_map
        elif data == 'peak':
            data = self.peak_map            
        elif data == "gaussian":
            data = self.gaussian_map
        
        # get threshold in log space
        threshold_logspace = self.get_threshold_logspace(data, sigma=sigma)
        threshold_abs = max(threshold, threshold_logspace)
        logger.info('Theortical threshold %s, logspace threshold %s, final abs threshold %s', 
                    threshold, threshold_logspace, threshold_abs)

        # find local maximum 
        coords = peak_local_max(data, min_distance=min_distance,
                                threshold_abs=threshold_abs)

        self.yp, self.xp = coords[:, 0], coords[:, 1]
        self.cand_src = pixel_to_skycoord(self.xp, self.yp, wcs=self.wcs)
        logger.info(f"Identified {len(self.cand_src)} candidate peaks.")


    # Depends if we preformed primary beam correction on short images or not 
    # if no primary beam correction previously, we should apply for below factor 
    def primary_correction(self, x0):
        """
        Compute the primary beam correction factor.

        Parameters
        ----------
        x0 : float
            Angular distance from beam center in degrees.

        Returns
        -------
        float
            Correction factor to normalize to beam center.
        """
        sigma = self.fwhm.to_value(u.degree) / (2 * np.sqrt(2 * np.log(2)))
        return norm.pdf(x0, scale=sigma) / norm.pdf(0, scale=sigma)
        
        
    def select_candidates(self, deepcatalogue, tabletype='selavy', sep=30,
                          mdlim=0.05, extlim=1.5, beamlim=1.2, bright=0.05, bright_sep=1):
        """
        Select final source candidates by comparing with deep image catalog.

        Parameters
        ----------
        deepcatalogue : str
            Path to the deep catalog (Aegean or Selavy format).
        tabletype : str
            'aegean' or 'selavy'. Affects how catalog is interpreted.
        sep : float
            Max separation from deep counterpart to count as match (arcsec).
        mdlim : float
            Minimum modulation index threshold.
        extlim : float
            Maximum allowed extent ratio (integrated/peak flux).
        beamlim : float
            Candidates must lie within this × FWHM from beam center.
        bright : float
            Threshold flux for a deep source to be considered bright (Jy).
        bright_sep : float
            Threshold separation to a bright source (arcmin).
        """
        # read deep catalogue
        self.read_catalogue(catalogue=deepcatalogue, tabletype=tabletype)
        
        # rule out candidates outside primary beam size
        self.beamlim = beamlim
        beamidx = self.cand_src.separation(self.beam_center) < beamlim*self.fwhm/2
        logger.info("Candidates inside {} primary beam: {}".format(beamlim, sum(beamidx)))
        
        # find the deep catalogue counterparts for each candidates
        self.deepidx, self.d2d, d3d = self.cand_src.match_to_catalog_sky(self.deep_src)
        # no counterparts in deep image
        nodeepidx = self.d2d.arcsec > sep
        logger.info("Candidates without deep counterparts <= %s arcsec: %s", sep, sum(nodeepidx))
        
        # calculate the modulation index 
        if tabletype == 'aegean':
            logger.info('Tabletype %s, No anti-primary beam correction...', tabletype)
            self.md = self.std_map[self.yp, self.xp] / self.deep_peak_flux[self.deepidx] 
        else:
            logger.info('Tabletype %s, Anti-primary beam correction...', tabletype)
            fc = self.primary_correction(x0=self.cand_src.separation(self.beam_center).degree)
            self.md = self.std_map[self.yp, self.xp] / self.deep_peak_flux[self.deepidx] / fc # primary beam correction factor
        
        mdidx = self.md > mdlim
        logger.info("Candidates with modulation index > {:.1%}: {}".format(mdlim, sum(mdidx)))
        
        # calculate the extend feature
        ext = (self.deep_int_flux / self.deep_peak_flux)[self.deepidx]
        extidx = ext < extlim
        logger.info("Candidates with compactness < {:.2f}: {}".format(extlim, sum(extidx)))

        # check number of close deep conterparts within 30 arcsec 
        idx1, _, _, _ = search_around_sky(coords1=self.cand_src, 
                                 coords2=self.deep_src, 
                                 seplimit=sep*u.arcsec)
        unique_idx, unique_counts = np.unique(idx1, return_counts=True)
        num_deep_close = np.zeros_like(self.cand_src, dtype=int)
        num_deep_close[unique_idx] = unique_counts
        self.deep_num = num_deep_close
        
        # check separation with bright deep sources
        deep_bright = self.deep_src[self.deep_int_flux > bright]
        if sum(self.deep_int_flux > bright) == 0:
            self.bright_sep_arcmin = np.full_like(self.cand_src, 999)
        else:
            _, d2d, _ = self.cand_src.match_to_catalog_sky(deep_bright)
            self.bright_sep_arcmin = d2d.arcmin
        
        nobrightidx = self.bright_sep_arcmin > bright_sep
        logger.info(f"Candidates without bright sources (> {bright} Jy) nearby (<= {bright_sep} arcmin): {sum(nobrightidx)}")

        # only select candidates that
        # 1. within the primary beam size
        # 2. have no countparts in deep image
        # 3. or have a countpart with md > 0.05 and ext < 1.5
        # 4. no bright sources nearby (< 1arcmin)
        self.final_idx = beamidx & (nodeepidx | (mdidx & extidx)) & nobrightidx
        logger.info("Final candidates: {}".format(sum(self.final_idx)))
        
        
    def save_csvtable(self, tablename="cand_catalogue", savevot=False):
        """Save selected candidates to a csv table
        
        tablename: str
            saved name
        savevot: bool
            if True, will also save a vot format table 
        """
        if not hasattr(self, 'final_idx') or self.final_idx.sum() == 0:
            logger.warning("No final candidates - saving empty table. ")
        
        t = Table()
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
        t['chi_square_log_sigma'] = self.get_sigma_logspace(self.chisq_map, 
                                                   np.array(t['chi_square']))
        t['chi_square_sigma'] = self.get_threshold_chisquare(value=np.array(t['chi_square']), 
                                                        num=self.num)
        
        # read the peak value at each pixel 
        t['peak_map'] = self.peak_map[self.yp, self.xp][self.final_idx]
        # calculate the sigma
        t['peak_map_log_sigma'] = self.get_sigma_logspace(self.peak_map, 
                                                 np.array(t['peak_map']))
        t['peak_map_sigma'] = self.get_threshold_peak(value=np.array(t['peak_map']), 
                                                 num=self.num)
        
        # if there's Gaussian map 
        if hasattr(self, 'gaussian_map'):
            t['gaussian_map'] = self.gaussian_map[self.yp, self.xp][self.final_idx]
            t['gaussian_map_sigma'] = self.get_sigma_logspace(self.gaussian_map, 
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
        logger.info("Save csv {}.csv".format(tablename))
        
        # save vot table
        if savevot and len(t) != 0:
            t.write("{}.vot".format(tablename), 
                    table_id="candidates", 
                    format="votable", 
                    overwrite=True)
            logger.info("Save vot {}.vot".format(tablename))
        
        
    def read_catalogue(self, catalogue, tabletype="selavy"):
        """
        Read and parse deep image catalog for source crossmatching.

        Parameters
        ----------
        catalogue : str
            Path to the catalog file (VOTable or FITS).
        tabletype : str
            Type of catalog ('aegean' or 'selavy').
        """
        logger.info(f"Reading deep catalogue: {catalogue} [format: {tabletype}]")

        if tabletype == 'aegean':
            with fits.open(catalogue) as hdul:
                self.catalogue = hdul[1].data.copy()

            self.deep_src = SkyCoord(self.catalogue['ra'], self.catalogue['dec'], unit=u.degree)
            self.deep_peak_flux = np.array(self.catalogue['peak_flux'])
            self.deep_int_flux = np.array(self.catalogue['int_flux'])

            self.deep_name = [
                'J' + src.ra.to_string(unit=u.hourangle, sep='', precision=0, pad=True) +
                src.dec.to_string(sep='', precision=0, alwayssign=True, pad=True)
                for src in self.deep_src
            ]
        
        else:
            if tabletype != 'selavy':
                logger.warning(f"Unrecognized table type '{tabletype}', defaulting to 'selavy'.")

            self.catalogue = Table.read(catalogue)
            self.deep_src = SkyCoord(self.catalogue['col_ra_deg_cont'],
                                     self.catalogue['col_dec_deg_cont'], unit=u.degree)
            self.deep_peak_flux = np.array(self.catalogue['col_flux_peak']) / 1e3  # Jy
            self.deep_int_flux = np.array(self.catalogue['col_flux_int']) / 1e3     # Jy
            self.deep_name = self.catalogue['col_component_name']
        
        
    def plot_fits(self, fitsname, imagename='plot_fits'):
        """
        Plot selected candidates on a FITS map using APLpy.

        Parameters
        ----------
        fitsname : str
            Input FITS filename.
        imagename : str
            Output PNG filename prefix.
        """
        f = aplpy.FITSFigure(fitsname, figsize=(8, 8))
        fix_aplpy_fits(f)
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
        plt.close('all')
        logger.info(f"Saved candidate overlay to {imagename}.png")


class Products:
    """
    Generate plots and cutouts for a list of selected candidate sources.
    """
    def __init__(self, final_csv, limit=-1, cand_name=None):
        """
        Parameters
        ----------
        final_csv : str
            Path to the CSV file with final candidate table.
        limit : int
            Number of top candidates to include (based on SNR/variability). -1 means all.
        cand_name : list of str, optional
            If provided, overrides default name list from table.
        """
        cat = Table.read(final_csv)
        if limit == -1:
            self.final_csv = cat
        else:
            cat.sort(keys=['peak_map', 'chi_square'], reverse=True)
            self.final_csv = cat[:limit]

        self.cand_name = cand_name if cand_name else self.final_csv['name']
        self.cand_src = SkyCoord(self.final_csv['ra'], self.final_csv['dec'], unit=u.degree)
        

    def generate_cutout(self, fitsname, radius=5, savename='output'):
        for src_name in self.cand_name:
            plot_cutout(src_name, fitsname, radius=radius,
                        name=f'{savename}_{src_name}')


    def generate_fits_cutout(self, fitsname, radius=5, savename='output'):
        for src_name in self.cand_name:
            save_fits_cutout(src_name, fitsname, radius=radius,
                             name=f'{savename}_{src_name}')


    def generate_fits_cube(self, imagelist, radius=5, savename='output'):
        for src_name in self.cand_name:
            save_fits_cube(src_name, imagelist, radius=radius,
                           name=f'{savename}_{src_name}')


    def generate_slices(self, imagelist, radius=5, vsigma=5, savename='output'):
        for src_name in self.cand_name:
            plot_slices(src_name, imagelist, radius=radius, vsigma=vsigma,
                        name=f'{savename}_{src_name}')
            

    def generate_lightcurve(self, imagelist, deepname, savename='output', savecsv=True):
        peak_flux_table = Table()
        local_rms_table = Table()

        for i, src_name in enumerate(self.cand_name):
            peak_flux, local_rms, timestamp = extract_lightcurve(
                src_name=src_name,
                imagelist=imagelist,
                deep_imagename=deepname
            )

            if i == 0:
                peak_flux_table['Time'] = timestamp
                local_rms_table['Time'] = timestamp

            peak_flux_table[src_name] = peak_flux
            local_rms_table[src_name] = local_rms

            plot_lightcurve(
                flux=peak_flux,
                times=timestamp,
                rms=local_rms,
                title=src_name,
                name=f'{savename}_{src_name}'
            )

        if savecsv:
            peak_flux_table.write(f"{savename}_peak_flux.csv", overwrite=True)
            local_rms_table.write(f"{savename}_local_rms.csv", overwrite=True)
            logger.info(f"Saved lightcurve CSVs to {savename}_peak_flux.csv")
            logger.info(f"Saved lightcurve CSVs to {savename}_local_rms.csv")
            
