#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 14:54:54 2022

@author: ywan3191
"""

"""
Generate a significance cube for transients detection
"""

import logging
import warnings
import numpy as np
from functools import partial
import dask.array as da
import os

from astropy.io import fits
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.convolution import convolve, Gaussian1DKernel
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales, pixel_to_skycoord
from astropy.utils.exceptions import AstropyWarning, AstropyDeprecationWarning

from skimage.feature import peak_local_max
from scipy.interpolate import interp2d

from .fastFunc import G2D
from .setting import CHUNK_0, CHUNK_1



warnings.filterwarnings('ignore', category=AstropyWarning, append=True)
warnings.filterwarnings('ignore',
                        category=AstropyDeprecationWarning, append=True)


# logging.basicConfig()
# logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)
logger = logging.getLogger(__name__)

class ArgumentError(Exception):
    pass


class Map:
    """Create a significance map
    
    """
    
    def __init__(self, image, idx=0):
        if isinstance(image, str):
            try: 
                self.fi = fits.open(image)
            except FileNotFoundError:
                logger.exception("Unable to open image %s" % image)
        elif isinstance(image, fits.HDUList):
            self.fi = image
        else:
            raise ArgumentError("Do not understand input image")
            
        if not (
            ("BMAJ" in self.fi[idx].header.keys())
            and ("BMIN" in self.fi[idx].header.keys())
            and ("BPA" in self.fi[idx].header.keys())
        ):

            raise KeyError("Image header does not have BMAJ, BMIN, BPA keywords")
            
        self.idx = idx
        self.BMAJ = self.fi[idx].header["BMAJ"] * u.deg
        self.BMIN = self.fi[idx].header["BMIN"] * u.deg
        self.BPA = self.fi[idx].header["BPA"] * u.deg
        self.w = WCS(self.fi[idx].header).celestial
        self.pixelscale = (proj_plane_pixel_scales(self.w)[1] * u.deg).to(u.arcsec)
        
        self.data = self.fi[idx].data.squeeze()
        self.fi.close()
        
        
        
    def imap(self, ktype, nx: int, ny: int, psf=None):
        """Output kernel shape (ny, nx) - consistent with FITS data
        """
        if nx % 2 == 0 or ny % 2 == 0:
            raise ArgumentError('nx or ny must be odd number')
        
        if ktype == 'gaussian':
            kernel = self._gaussian(nx, ny)
            simage = self._smooth(self.data, kernel)
        
        elif ktype == 'psf':
            kernel = self._psf(psf, nx, ny)
            simage = self._smooth(self.data, kernel)
            
        elif ktype == 'combine':
            k1= self._psf(psf, nx, ny)
            k2 = self._gaussian(nx, ny)
            simage = self._smooth(self._smooth(self.data, k1), k2)
            kernel = self._smooth(k1, k2)
            
        else:
            raise ArgumentError('Do not understand kernel type %s' % ktype)
            
        self.map = simage
        self.kernel = kernel
    
    
    
    def tofits(self, fitsname: str):
        """Save the significance map to FITS file 
        """
        hdu = self.fi
        data = hdu[self.idx].data 
        
        hdu[self.idx].data = self.map.reshape(data.shape)
        hdu.writeto(fitsname)
        
        
        
        
    def _smooth(self, image, kernel):
        """
        kernal: np.array 2d, ideally the (dirty) psf of FITS image
        image:  np.array 2d, read from FITS image 
        
        return
        simage: np.array 2d, significance map of the FITS image, 
                can get from the convolution of the kernal and the image
        """
        
        simage = convolve(image, kernel)
        
        return simage
            
            
            
    def _gaussian(self, nx, ny):
    # def tgaussian(self, nx, ny):
        """Build a Gaussian kernel with size of nx, ny
        """
        xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
        g = G2D(
            nx//2, ny//2, 
            (self.BMAJ/self.pixelscale).to(u.dimensionless_unscaled), 
            (self.BMIN/self.pixelscale).to(u.dimensionless_unscaled), 
            self.BPA
            )
        return g(xx, yy)
    
    
    
    def _psf(self, psf, nx, ny):
    # def tpsf(self, psf, nx, ny):
        """Build kernel with dirty beam 
        """
        if isinstance(psf, str):
            try: 
                hdu = fits.open(psf)
            except FileNotFoundError:
                logger.exception("Unable to open image %s" % psf)    
        else:
            raise ArgumentError("Do not understand input image")
            
        kernel = self._kcutout(
                hdu[0].data.squeeze(), 
                nx, 
                ny, 
            )
        return kernel
    
    
                
    def _kcutout(self, kernel, nx, ny):
        """
        kernel size can only be odd - reduce the size to a (nx, ny) array (from centre)
        kernel: np.array 2d
        nx:     int
        ny:     int
        
        return: kcutout: np.array 2d 
        """
        
        # find the central position
        x, y = np.where(kernel == 1)[0][0], np.where(kernel == 1)[1][0]
        kcutout = Cutout2D(kernel, (x, y), (ny, nx))
        
        return kcutout.data





class Cube:
    """Create a significance cube for transients detection
    
    Example usage:
        TBD
        
    Args:
        imagelist: a list of names of short FITS image. 
        idx: hdu index, default is 0 
    """
    
    def __init__(self, imagelist, idx=0):
        
        # if isinstance(imagelist, list):
        #     try:
        #         self.hdulist = [fits.open(image) for image in imagelist]
        #     except FileNotFoundError:
        #         logger.exception("Unable to open this image list")
                
        # else:
        #     raise ArgumentError("Do not understand input image list")
        
        if not isinstance(imagelist, list):
            raise ArgumentError("Do not understand input image list")
            
        self.idx = idx
        self.imagelist = imagelist
        
        
            
    def icube(self, ktype=None, nx=None, ny=None, psflist=None, fitsfolder=None):
        """To generate a cube for further processing, 
            Can be smoothed by different ktype (Gaussian, dirty beam)
            Or simply generate a cube using original dataset
            
        """
        
        # set the original cube size 
        sigcube = []
        
        if ktype == None:
            logger.info('Save original dataset to a cube...')
            
            for i, image in enumerate(self.imagelist):
                m = Map(image, self.idx)
                sigcube.append(m.data)
                
            self.sigcube = np.array(sigcube)
            
            return None
            
        # ============
        # if using smoothed function 
        if ktype != 'gaussian':
            self._psf_init(psflist)
        else:
            psflist = [None] * len(self.imagelist)
            
        # build the cube using for loop
        for i, image in enumerate(self.imagelist):
            logger.info('Processing image number %s: %s' % (i, self.imagelist[i]))
            m = Map(image, self.idx)
            m.imap(ktype, nx, ny, psflist[i])
            sigcube.append(m.map)
            
            # save the fits
            if fitsfolder != None:
                m.tofits(fitsname="{}/map_{}.fits".format(fitsfolder, i))
            
        self.kernel = m.kernel
        self.sigcube = np.array(sigcube)
       
        
        
        
    # save the significance cube
    def savecube(self, savename='sigcube.npy'):
        """Save cube to npy format. 
        """
        np.save(savename, self.sigcube)
        
        
        
    # save the original fits to a cube
    def save_oricube(self, savename=None):
        """Save the original fits to a cube
        """
        oricube = []
        
        for i, image in enumerate(self.imagelist):
            m = Map(image)
            oricube.append(m.data)
            
        self.oricube = np.array(oricube)
        
        if savename != None:
            np.save(savename, self.oricube)
        
            
        
                    
    def _psf_init(self, psflist):
        """Check if the input psf list is in correct format
        """
        if isinstance(psflist, list):
            if not len(self.imagelist) == len(psflist):
                raise ArgumentError('The length of image list is not '
                                    'consistent with the length of psf list. ')
            self.psflist = psflist
        else:
            raise ArgumentError("Do not understand input psf list")
            
            
            
            
    def remove_bad_images(self, sigma=2):
        """Get the rms threshold, to decide which image should be ruled out 
        
        sigma: float, the outlier image rms threshod. 
            images with higher rms will be ruled out
        """
        # get the rms of 
        rmslist = np.nanstd(self.sigcube, axis=(1, 2))
        # get threshold
        # threshold = np.nanmean(rmslist) + 2*np.nanstd(rmslist)
        threshold = 2 * np.nanmedian(rmslist)
        
        logger.info("Median rms level {:.0f} uJy/beam".format(np.nanmedian(rmslist)*1e6))
        logger.info("Remove > {} sigma outliers, "
                    "i.e. rms above {} Jy/beam".format(sigma, threshold))
        logger.info("Bad image rms level: ")
        logger.info(rmslist[rmslist>threshold])
        logger.info([self.imagelist[i] for i in np.where(rmslist > threshold)[0]])
        
        # remove image with rms > threshold
        self.sigcube = self.sigcube[~(rmslist>threshold)]
        logger.info("Remove {} of {} images".format(sum(rmslist>threshold), 
                                                 len(self.imagelist)))

            
        
class Filter:
    """Create a kernel filter in time axis, to select transients candidates
    """
    
    def __init__(self, sigcube):
        """Cube shape in (time, ny, nx), mp.3darray
        """
        
        if isinstance(sigcube, np.ndarray):
            if len(sigcube.shape) == 3:
                self.sigcube = sigcube
            else:
                raise ArgumentError('The dimension of the cube should be 3d. ')
        else:
            raise ArgumentError('The input cube should be np.array')
            
        # get local rms
        self.cube_local_rms()
        
    
    def fmap(self, ktype, width=4, nprocess=4):
        self.width = width

        if ktype == "chisquare":
            self.map = self._chimap()
            
        elif ktype == "peak":
            self.map = self._peakmap()
            
        elif ktype == "std":
            self.map = self._stdmap()

        elif ktype == "gaussian":
            self.map = self._gmap(nprocess)
        
        
        
    def tofits(self, fitsname: str, imagename=None):
        """Save the significance map to FITS file 
        """
        # check if there's imagename 
        if not hasattr(self, 'imagename'):
            self._readfits(imagename)
        
        hdu = self.fi
        data = hdu.data 
       
        hdu.data = self.map.reshape(data.shape)
        hdu.writeto(fitsname, overwrite=True)
        
        
        
    # def local_max(self, min_distance=30, sigma=5, imagename=None):
    #     '''Find the local maxium of an image
        
    #     sigma: identify blobs above a specfic sigma threshold
    #     min_distance: pixel number of the minimal distance of two neighbours blobs
    #     '''
    #     rms, mean = np.nanstd(self.map), np.nanmean(self.map)
    #     threshold = sigma*rms + mean
        
    #     logger.info("Threshold rms = {}, mean = {}".format(rms, mean))
    #     logger.info('Threshold is {} sigma = {}'.format(sigma, threshold))
        
    #     xy = peak_local_max(self.map, min_distance=min_distance, 
    #                         threshold_abs=threshold)
        
    #     self.xy_peak = xy
        
    #     # convert pixel to skycoord
    #     if not hasattr(self, 'imagename'):
    #         self._readfits(imagename)
        
    #     yp, xp = xy[:, 0], xy[:, 1]
        
    #     self.coord = pixel_to_skycoord(xp, yp, wcs=self.wcs)
        


    def _readfits(self, imagename):
        
        self.fi = fits.open(imagename)[0]
        self.wcs = WCS(self.fi.header)
        self.imagename = imagename
        
        
        
    def _filter(self, kernel, axis=0):
        """Convolve, get the maximum value
        """
        self.smocube =  np.apply_along_axis(lambda m: convolve(m, kernel), 
                                   axis=axis, arr=self.sigcube)

        return np.nanmax(self.smocube, axis=0) - np.nanmean(self.smocube, axis=0)
    


    def _chimap(self):
        """Chi-square map
        """
        # local_rms = np.std(self.sigcube, axis=(1, 2))
        # return np.apply_along_axis(lambda m: self._chisquare(m, local_rms), axis=0, arr=self.sigcube)
        
        # covert cubes into dask array
        da_sig = da.from_array(self.sigcube)
        da_rms = da.from_array(self.rmscube)
        logger.info("chisquare map calculation using dask -- sigcube chunksize: {}".format(da_sig.chunksize))
        logger.info("chisquare map calculation using dask -- rmscube chunksize: {}".format(da_rms.chunksize))

        # freedom
        nu = da_sig.shape[0] - 1
        # mean, rms
        mean = da.nanmean(da_sig, axis=0)
        # rms = np.nanstd(self.sigcube, axis=(1, 2))
        # rms = self.cube_local_rms()
        
        # for each data point
        # data = (self.sigcube - mean) / rms
        data = (da_sig - mean) / da_rms
        
        res = da.sum(da.square(data), axis=0) / nu
        logger.info("chi map res chunksize: {}".format(res.chunksize))
        # nres = res.compute(num_workers=1)
        nres = res.compute()
        # np.save("chimap", nres)
        return nres
        
    


    def _gaussian(self):
        return Gaussian1DKernel(stddev=self.width)

    def _gmap(self, nprocess):
        sigcube_t = self.sigcube.transpose(1,2,0).copy(order="C")
        logger.info("Transposed sigcube shape: {}".format(sigcube_t.shape))
        time_dim = sigcube_t.shape[2]
        da_sigcube_t = da.from_array(sigcube_t, chunks=(CHUNK_0,CHUNK_1,time_dim))
        logger.info("gaussian map calculation using dask -- transposed sigcube chunksize: {}".format(da_sigcube_t.chunksize))
        
        da_gmap = _get_gmap(da_sigcube_t)
        gmap = da_gmap.compute(scheduler="processes", num_workers=nprocess)
        # np.save("gaussianmap", gmap)
        return gmap
        

    def _chisquare(self, peak_flux, local_rms, m=1):
        mask = np.isnan(peak_flux) + np.isnan(local_rms)
        peak_flux = np.ma.masked_array(peak_flux, mask, dtype=float)
        local_rms = np.ma.masked_array(local_rms, mask, dtype=float)
        # freedom
        nu = np.sum(~mask) - m
        mean_flux = np.average(peak_flux, weights=np.power(local_rms, -2))
    
        return np.sum(np.power((peak_flux - mean_flux)/local_rms, 2)) / nu


    def _peakmap(self):
        """Get peak map, sensitive to single flare event
        """

        # convert cubes to dask array
        da_sig = da.from_array(self.sigcube)
        da_rms = da.from_array(self.rmscube)
        logger.info("peak map calculation using dask -- sigcube chunksize: {}".format(da_sig.chunksize))
        logger.info("peak map calculation using dask -- rmscube chunksize: {}".format(da_rms.chunksize))

        snr = da_sig / da_rms
        
        # return (np.nanmax(self.sigcube, axis=0) - np.nanmedian(self.sigcube, axis=0)) 
        res = (da.nanmax(snr, axis=0) - da.nanmedian(snr, axis=0)) 
        logger.info("peak map calculation using dask -- result chunksize: {}".format(res.chunksize))
        # nres = res.compute(num_workers=1)
        nres = res.compute()
        # np.save("peakmap", nres)
        return nres
    
    
    def _stdmap(self):
        """Get the standard deviation map, useful for modulation index 
        """        
        # convert cubes to dask array
        da_sig = da.from_array(self.sigcube)
        logger.info("std map calculation using dask -- sigcube chunksize: {}".format(da_sig.chunksize))

        res = da.nanstd(da_sig, axis=0)
        logger.info("std map calculation using dask -- result chunksize: {}".format(res.chunksize))
        # nres = res.compute(num_workers=1)
        nres = res.compute()
        # np.save("stdmap", nres)
        return nres
    
    
    
    def cube_local_rms(self):
        """input: sigcube
            output: rms cube
        """
        logger.info("Calculating local rms...")
        
        rmscube = np.zeros(self.sigcube.shape, dtype=np.float32)
        rmscube[:] = np.nan
        for i in range(self.sigcube.shape[0]):
            data_image = self.sigcube[i]
            rms_image = self.get_local_rms(data=data_image)
            rmscube[i] = rms_image
            
        logger.info("Finish calculating local rms.")
            
        self.rmscube = rmscube
    
    
    
    def get_local_rms(self, data, window_size=299, step=69):
        """Calculate the local rms
            get the rms in several key pixel, and do interpolate 2d
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
        

### extra functions needed for gaussian map 
def _conv_1d(arr, kernel):
        arr_res = convolve(arr, kernel)
        return arr_res

def _process_block(arr1, block_info=None):
    kernel = Gaussian1DKernel(stddev=4)
    res = np.apply_along_axis(_conv_1d, axis=2, arr=arr1, kernel=kernel)
    print(block_info)
    print("process: {}".format(os.getpid()))
    return res


def _get_gmap(sigcube):
    # tt = da.map_blocks(_process_block, sigcube, chunks=(100,100,40))
    tt = da.map_blocks(_process_block, sigcube)
    logger.info("gaussian map calculation using dask -- blocked sigcube chunksize: {}".format(tt.chunksize))
    res = da.nanmax(tt, axis=2) - da.nanmean(tt, axis=2)  
    logger.info("gaussian map calculation using dask -- result chunksize: {}".format(res.chunksize))
    return res  