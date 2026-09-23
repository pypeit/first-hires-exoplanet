#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:12:20 2020

Additional code for the pyodine package by René Tronsgaard Rasmussen,
to find cosmics and bad pixels in the spectra and return a bad pixel mask.

@author: Paul Heeren
"""

import numpy as np
import matplotlib.pyplot as plt
from pyodine.components import Observation #, SummedObservation

from . import correct_spec


class BadPixelMask:
    """A class of a bad pixel mask
    
    :param spec: The spectrum for which to compute the mask.
    :type spec: :class:`Observation`
    :param cutoff: Flux changes larger than the cutoff might be due to a 
        bad pixel (default: 0.18).
    :type cutoff: float
    :param plotting: If True, show diagnostic plots (default: False).
    :type plotting: bool
    """
    def __init__(self, spec, cutoff=0.18, plotting=False):
        self.cutoff = cutoff
        self.plotting = plotting
        
        if isinstance(spec, Observation):
            self.nord = spec.nord
            self.npix = spec.npix
        else:
            self.nord = spec.flux.shape[0]
            self.npix = spec.flux.shape[1]
        
        self.compute_mask(spec)
        
    
    def compute_mask(self, spec):
        """Check for cosmics and bad pixels.
        Based on 'find_cosmic.pro' (written by S. Reffert, LSW Heidelberg).
        The cutoff controls the behaviour of the bad_pixel finding algorithm: 
        values of 0.15--0.2 seem appropriate (lower values make it more 
        sensitive to little wobbles in the spectrum, and with values higher 
        than 0.2 it may fail to identify real spikes).
        
        :param spec: The spectrum for which to compute the mask.
        :type spec: :class:`Observation`
        """
        self.mask = np.zeros((self.nord,self.npix))
        for no in range(self.nord):
            flux = spec[no].flux
            cont = spec[no].cont
            
            norm_flux = flux / cont
            
            alt = 0.
            oldflag = 0
            
            for i in range(3,self.npix-3):
                neu = norm_flux[i] - norm_flux[i-1]
                
                if np.abs(neu - alt) >= self.cutoff:
                    newflag = 2
                    if oldflag == 0:
                        start = i-2
                    
                    if i == self.npix-3:
                        if self.plotting:
                            plt.plot(norm_flux[start:i])
                            plt.title('Order {}, pixel {}-{}'.format(
                                no, start, i))
                            plt.show()
                        for j in range(start,i):
                            self.mask[no,j] = 1
                
                else:
                    newflag = 0
                    if oldflag == 2:
                        newflag = 1
                        
                    if oldflag == 1 or (
                            oldflag ==2 and i == self.npix-3):
                        if self.plotting:
                            plt.plot(norm_flux[start:i])
                            plt.title('Order {}, pixel {}-{}'.format(
                                no, start, i))
                            plt.show()
                        for j in range(start,i):
                            self.mask[no,j] = 1
            
                alt = neu
                oldflag = newflag
            
        self.bad_pixels = np.where(self.mask==1)
    
    def __len__(self):
        """The dedicated length-method
        
        :return: The number of bad pixels in the mask.
        :rtype: int
        """
        return len(self.bad_pixels[0])