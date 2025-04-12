#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate a significance cube for transients detection.
Created on Wed Jun  1 14:54:54 2022

This module includes classes for generating significance maps from FITS images and stacking them
into a cube for transient detection. The `Map` class is used to generate a single significance
map, and the `Cube` class constructs a 3D significance cube over time.

@author: ywan3191
"""

import logging
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.convolution import convolve, Gaussian1DKernel
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.interpolate import interp2d
from vaster.vastfast.fastFunc import G2D

logger = logging.getLogger(__name__)


class ArgumentError(Exception):
    """Custom exception for invalid arguments."""
    pass


class Map:
    """
    Generate a 2D significance map from a single FITS image.

    Attributes
    ----------
    data : numpy.ndarray
        Image data array.
    header : fits.Header
        FITS header.
    BMAJ, BMIN, BPA : astropy.units.Quantity
        Beam size and position angle.
    w : astropy.wcs.WCS
        World coordinate system.
    pixelscale : astropy.units.Quantity
        Pixel scale in arcseconds.
    map : numpy.ndarray
        Significance map (after smoothing).
    kernel : numpy.ndarray
        Kernel used in smoothing.
    """
    
    def __init__(self, image, idx=0):
        """
        Initialize the Map from a FITS file or HDUList.

        Parameters
        ----------
        image : str or fits.HDUList
            Input FITS file path or already opened HDU list.
        idx : int, optional
            Index of the image extension to read (default is 0).
        """
        if isinstance(image, str):
            try:
                with fits.open(image) as hdul:
                    self.header = hdul[idx].header.copy()
                    self.data = hdul[idx].data.squeeze().copy()
            except FileNotFoundError:
                logger.exception(f"Unable to open image {image}")
                raise
        elif isinstance(image, fits.HDUList):
            self.header = image[idx].header.copy()
            self.data = image[idx].data.squeeze().copy()
        else:
            raise ArgumentError("Invalid image input. Must be a file path or HDUList.")

        for key in ("BMAJ", "BMIN", "BPA"):
            if key not in self.header.keys():
                raise KeyError(f"Missing required FITS header key: {key}")

        self.BMAJ = self.header["BMAJ"] * u.deg
        self.BMIN = self.header["BMIN"] * u.deg
        self.BPA = self.header["BPA"] * u.deg
        self.w = WCS(self.header).celestial
        self.pixelscale = (proj_plane_pixel_scales(self.w)[1] * u.deg).to(u.arcsec)
        
        
    def imap(self, ktype, nx, ny, psf=None):
        """
        Generate the significance map using a chosen kernel.

        Parameters
        ----------
        ktype : str
            Type of kernel to use ('gaussian', 'psf', or 'combine').
        nx : int
            Width of the kernel (must be odd).
        ny : int
            Height of the kernel (must be odd).
        psf : str, optional
            Path to PSF image (required for 'psf' and 'combine').
        """
        if nx % 2 == 0 or ny % 2 == 0:
            raise ArgumentError("Kernel dimensions nx and ny must be odd.")

        if ktype == 'gaussian':
            kernel = self._gaussian(nx, ny)
            self.map = self._smooth(self.data, kernel)

        elif ktype == 'psf':
            kernel = self._psf(psf, nx, ny)
            self.map = self._smooth(self.data, kernel)

        elif ktype == 'combine':
            k1 = self._psf(psf, nx, ny)
            k2 = self._gaussian(nx, ny)
            self.map = self._smooth(self._smooth(self.data, k1), k2)
            kernel = self._smooth(k1, k2)

        else:
            raise ArgumentError(f"Unknown kernel type: {ktype}")

        self.kernel = kernel
    
    
    def tofits(self, fitsname):
        """
        Save the generated significance map to a FITS file.

        Parameters
        ----------
        fitsname : str
            Output filename.
        """
        hdu = fits.PrimaryHDU(data=self.map.reshape(self.data.shape), header=self.header)
        hdu.writeto(fitsname, overwrite=True)

        
    def _smooth(self, image, kernel):
        """
        Convolve image with kernel.

        Parameters
        ----------
        image : numpy.ndarray
            2D input image. 
        kernel : numpy.ndarray
            2D kernel. Ideally the (dirty) psf of FITS image

        Returns
        -------
        numpy.ndarray
            Smoothed image.
        """
        return convolve(image, kernel)
            
            
    def _gaussian(self, nx, ny):
        """
        Create a 2D Gaussian kernel.

        Parameters
        ----------
        nx : int
            Kernel width.
        ny : int
            Kernel height.

        Returns
        -------
        numpy.ndarray
            Gaussian kernel.
        """
        xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
        g = G2D(
            nx // 2,
            ny // 2,
            (self.BMAJ / self.pixelscale).to(u.dimensionless_unscaled),
            (self.BMIN / self.pixelscale).to(u.dimensionless_unscaled),
            self.BPA,
        )
        return g(xx, yy)
    
    
    def _psf(self, psf, nx, ny):
        """
        Extract a PSF kernel from a FITS file.

        Parameters
        ----------
        psf : str
            Path to PSF FITS file.
        nx : int
            Desired kernel width.
        ny : int
            Desired kernel height.

        Returns
        -------
        numpy.ndarray
            Cutout of PSF kernel.
        """
        if not isinstance(psf, str):
            raise ArgumentError("PSF must be a file path string.")

        try:
            with fits.open(psf) as hdul:
                return self._kcutout(hdul[0].data.squeeze().copy(), nx, ny)
        except FileNotFoundError:
            logger.exception(f"Unable to open PSF image {psf}")
            raise

    
    def _kcutout(self, kernel, nx, ny):
        """
        Cut a centered region from a PSF kernel.

        Parameters
        ----------
        kernel : numpy.ndarray
            Full 2D PSF image.
        nx : int
            Size of kernel width.
        ny : int
            Size of kernel height.

        Returns
        -------
        numpy.ndarray
            Centered cutout of the kernel.
        """
        # find the central position
        x, y = np.where(kernel == 1)[0][0], np.where(kernel == 1)[1][0]
        kcutout = Cutout2D(kernel, (x, y), (ny, nx))
        return kcutout.data


class Cube:
    """
    Create a significance cube from a list of FITS images.

    This class constructs a 3D cube for transient search across multiple images
    by applying kernel filtering or using raw image data.

    Attributes
    ----------
    imagelist : list of str
        List of FITS file paths.
    idx : int
        Index to access FITS extension.
    sigcube : numpy.ndarray
        Generated 3D cube of significance maps.
    oricube : numpy.ndarray
        Cube of original image data.
    """
    
    def __init__(self, imagelist, idx=0):
        """
        Initialize Cube with list of image paths.

        Parameters
        ----------
        imagelist : list of str
            List of FITS file paths.
        idx : int
            Index of FITS extension to use.
        """
        if not isinstance(imagelist, list):
            raise ArgumentError("imagelist must be a list of FITS file paths.")

        self.imagelist = imagelist
        self.idx = idx

        
    def icube(self, ktype=None, nx=None, ny=None, psflist=None, fitsfolder=None):
        """
        Build a 3D significance cube.

        Parameters
        ----------
        ktype : str, optional
            Type of kernel ('gaussian', 'psf', 'combine', or None).
        nx : int, optional
            Kernel width (must be odd).
        ny : int, optional
            Kernel height (must be odd).
        psflist : list of str, optional
            List of PSF image paths (if required).
        fitsfolder : str, optional
            Directory to save output FITS maps.
        """
        sigcube = []
        
        if ktype is None:
            logger.info("Using raw data for significance cube.")
            for image in self.imagelist:
                m = Map(image, self.idx)
                sigcube.append(m.data)
            self.sigcube = np.array(sigcube)
            return
            
        # ============
        # if using smoothed function 
        if ktype != 'gaussian':
            self._psf_init(psflist)
        else:
            psflist = [None] * len(self.imagelist)
            
        # build the cube using for loop
        for i, image in enumerate(self.imagelist):
            logger.info(f"Processing image {i}: {image}")
            m = Map(image, self.idx)
            m.imap(ktype, nx, ny, psflist[i])
            sigcube.append(m.map)
            
            if fitsfolder != None:
                m.tofits(fitsname="{}/map_{}.fits".format(fitsfolder, i))
            
        self.kernel = m.kernel
        self.sigcube = np.array(sigcube)
        
        
    def savecube(self, savename='sigcube.npy'):
        """
        Save the significance cube to a .npy file.

        Parameters
        ----------
        savename : str
            Filename for saving the cube.
        """
        np.save(savename, self.sigcube)
        
        
    def save_oricube(self, savename=None):
        """
        Save the original (unsmoothed) cube.

        Parameters
        ----------
        savename : str, optional
            Output filename for .npy array.
        """
        oricube = []
        for image in self.imagelist:
            m = Map(image, self.idx)
            oricube.append(m.data)
        self.oricube = np.array(oricube)
        
        if savename != None:
            np.save(savename, self.oricube)
        
            
    def _psf_init(self, psflist):
        """
        Validate the PSF list.

        Parameters
        ----------
        psflist : list of str
            List of PSF file paths.

        Raises
        ------
        ArgumentError
            If the list is not valid or lengths mismatch.
        """
        if not isinstance(psflist, list) or len(psflist) != len(self.imagelist):
            raise ArgumentError("PSF list must match image list in length.")

        self.psflist = psflist
            
            
            
            
    def remove_bad_images(self, sigma=2):
        """
        Remove images with excessively high RMS values.

        Parameters
        ----------
        sigma : float
            Threshold multiplier for median RMS.
        """
        rmslist = np.nanstd(self.sigcube, axis=(1, 2))
        threshold = sigma * np.nanmedian(rmslist)
        
        logger.info("Median rms level {:.0f} uJy/beam".format(np.nanmedian(rmslist)*1e6))
        logger.info("Remove > {} sigma outliers, "
                    "i.e. rms above {} Jy/beam".format(sigma, threshold))
        logger.info("Bad image rms level: ")
        logger.info(rmslist[rmslist>threshold])
        logger.info([self.imagelist[i] for i in np.where(rmslist > threshold)[0]])
        
        # remove image with rms >= threshold and remove empty images
        ind = (rmslist < threshold) & (rmslist > 0)
        self.sigcube = self.sigcube[ind]
        logger.info("Remove {} of {} images".format(sum(~ind), len(self.imagelist)))

        
class Filter:
    """
    Create a temporal kernel filter to identify transient candidates from a significance cube.

    Attributes
    ----------
    sigcube : numpy.ndarray
        The input 3D data cube with shape (time, ny, nx).
    rmscube : numpy.ndarray
        3D cube of local RMS values.
    map : numpy.ndarray
        2D filtered significance map.
    """

    def __init__(self, sigcube):
        """
        Initialize the Filter with a 3D significance cube.

        Parameters
        ----------
        sigcube : numpy.ndarray
            3D array with shape (time, ny, nx).
        """
        if not isinstance(sigcube, np.ndarray) or sigcube.ndim != 3:
            raise ArgumentError("Input cube must be a 3D numpy array.")

        self.sigcube = sigcube
        self.cube_local_rms()
            
        
    def fmap(self, ktype, width=4):
        """
        Apply a filtering method across the time axis.

        Parameters
        ----------
        ktype : str
            Type of filter to apply: 'chisquare', 'peak', 'std', or 'gaussian'.
        width : int, optional
            Width of the temporal kernel for convolution (default is 4).
        """
        self.width = width
        
        if ktype == "chisquare":
            self.map = self._chimap()
            
        elif ktype == "peak":
            self.map = self._peakmap()
            
        elif ktype == "std":
            self.map = self._stdmap()
        
        else:
            if ktype == 'gaussian':
                kernel = self._gaussian()
                
            self.map = self._filter(kernel)
        
        
    def tofits(self, fitsname, imagename):
        """
        Save the filtered 2D significance map to a FITS file.

        Parameters
        ----------
        fitsname : str
            Output FITS file name.
        imagename : str
            Input FITS image to copy header info from.
        """
        with fits.open(imagename) as hdul:
            header = hdul[0].header.copy()
            shape = hdul[0].data.shape

        hdu = fits.PrimaryHDU(data=self.map.reshape(shape), header=header)
        hdu.writeto(fitsname, overwrite=True)

        
    def _filter(self, kernel, axis=0):
        """
        Apply a 1D kernel filter across the specified axis.

        Parameters
        ----------
        kernel : numpy.ndarray
            1D convolution kernel.
        axis : int, optional
            Axis along which to apply filter (default is 0).

        Returns
        -------
        numpy.ndarray
            2D filtered map.
        """
        smoothed = np.apply_along_axis(lambda m: convolve(m, kernel), axis, self.sigcube)
        return np.nanmax(smoothed, axis=0) - np.nanmean(smoothed, axis=0)
    

    def _chimap(self):
        """
        Compute a chi-square map using local RMS.

        Returns
        -------
        numpy.ndarray
            2D chi-square map.
        """
        nu = self.sigcube.shape[0] - 1
        mean = np.nanmean(self.sigcube, axis=0)
        data = (self.sigcube - mean) / self.rmscube
        return np.sum(data**2, axis=0) / nu
        
    
    def _gaussian(self):
        """
        Create a Gaussian kernel for temporal convolution.

        Returns
        -------
        Gaussian1DKernel
            1D Gaussian kernel.
        """
        return Gaussian1DKernel(stddev=self.width)
    
    
    def _chisquare(self, peak_flux, local_rms, m=1):
        mask = np.isnan(peak_flux) + np.isnan(local_rms)
        peak_flux = np.ma.masked_array(peak_flux, mask, dtype=float)
        local_rms = np.ma.masked_array(local_rms, mask, dtype=float)
        # freedom
        nu = np.sum(~mask) - m
        mean_flux = np.average(peak_flux, weights=np.power(local_rms, -2))
        return np.sum(np.power((peak_flux - mean_flux)/local_rms, 2)) / nu


    def _peakmap(self):
        """
        Compute peak signal-to-noise ratio map.

        Returns
        -------
        numpy.ndarray
            2D peak map.
        """
        snr = self.sigcube / self.rmscube
        logger.debug('snr')
        logger.debug(snr)
        return (np.nanmax(snr, axis=0) - np.nanmedian(snr, axis=0)) 
    
    
    def _stdmap(self):
        """
        Compute standard deviation across the time axis.

        Returns
        -------
        numpy.ndarray
            2D standard deviation map.
        """   
        return np.nanstd(self.sigcube, axis=0)
    
    
    def cube_local_rms(self):
        """
        Estimate local RMS for each time slice in the cube.

        Sets
        ----
        rmscube : numpy.ndarray
            Cube of local RMS estimates.
        """
        logger.info("Calculate local rms...")
        
        rmscube = []
        for i in range(self.sigcube.shape[0]):
            data_image = self.sigcube[i]
            rms_image = self.get_local_rms(data=data_image)
            rmscube.append(rms_image)
            
        logger.info("Calculate local rms finished...")
        self.rmscube = np.array(rmscube, dtype='float32')
    
    
    def get_local_rms(self, data, window_size=299, step=69):
        """
        Estimate local RMS using a sliding window and 2D interpolation.

        Parameters
        ----------
        data : numpy.ndarray
            2D image data.
        window_size : int, optional
            Size of the sliding window (default is 299).
        step : int, optional
            Step size between window centers (default is 69).

        Returns
        -------
        numpy.ndarray
            2D interpolated RMS map.
        """        
        ridx = int(window_size/2)
        xp = np.arange(ridx+1, data.shape[0]-ridx, step=step)
        xx, yy = np.meshgrid(xp, xp)
        
        view_data = np.lib.stride_tricks.sliding_window_view(data, (window_size, window_size))
        rms_i = np.nanstd(view_data[yy-ridx, xx-ridx], axis=(-2, -1))

        # do interpolate
        f = interp2d(x=xp, y=xp, z=rms_i)
        xnew = np.arange(data.shape[0])
        rms = f(x=xnew, y=xnew)
        rms = np.array(rms, dtype='float32')
        
        return rms
        
        