#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 28 16:51:48 2022

@author: ywan3191
"""

'''
Build a 2D kernel; 
Convolve with the FITS images; 
Get a significance cube. 
'''

import numpy as np
from astropy import units as u


pa_offset = 90 * u.deg
# pa_offset = 0 * u.deg


class G2D:
    """2D Gaussian for use as a kernel.
    Example usage:
        create the kernel:
        g = G2D(x0, y0, fwhm_x, fwhm_y, PA)
        and return the kernel:
        g(x, y)
    Args:
        x0 (float): the mean x coordinate (pixels)
        y0 (float): the mean y coordinate (pixels)
        fwhm_x (float): the FWHM in the x coordinate (pixels)
        fwhm_y (float): the FWHM in the y coordinate (pixels)
        pa (float): the position angle of the Gaussian (E of N) as a Quantity or in
            radians.
    """

    def __init__(self, x0: float, y0: float, fwhm_x: float, fwhm_y: float, pa: float):
        self.x0 = x0
        self.y0 = y0
        self.fwhm_x = fwhm_x
        self.fwhm_y = fwhm_y
        # adjust the PA to agree with the selavy convention
        # E of N
        self.pa = pa - pa_offset
        self.sigma_x = self.fwhm_x / 2 / np.sqrt(2 * np.log(2))
        self.sigma_y = self.fwhm_y / 2 / np.sqrt(2 * np.log(2))

        self.a = (
            np.cos(self.pa) ** 2 / 2 / self.sigma_x ** 2
            + np.sin(self.pa) ** 2 / 2 / self.sigma_y ** 2
        )
        self.b = (
            np.sin(2 * self.pa) / 2 / self.sigma_x ** 2
            - np.sin(2 * self.pa) / 2 / self.sigma_y ** 2
        )
        self.c = (
            np.sin(self.pa) ** 2 / 2 / self.sigma_x ** 2
            + np.cos(self.pa) ** 2 / 2 / self.sigma_y ** 2
        )

    def __call__(self, x: float, y: float) -> np.ndarray:
        """Return the kernel evaluated at given pixel coordinates.
        Args:
            x (float): x coordinate for evaluation
            y (float): y coordinate for evaluation
        Returns:
            np.ndarray: the kernel evaluated at the given coordinates
        """
        return np.exp(
            -self.a * (x - self.x0) ** 2
            - self.b * (x - self.x0) * (y - self.y0)
            - self.c * (y - self.y0) ** 2
        )
    
    
    
    