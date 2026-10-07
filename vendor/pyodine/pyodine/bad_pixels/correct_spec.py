#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:12:20 2020

@author: paul
"""

import numpy as np
from scipy.interpolate import splrep, splev
from ..components import Observation, NormalizedObservation
import logging

def correct_spectrum(spec, weights, orders):
    """Correct bad pixel regions in spectra
    Used to correct regions of spectra with weights=0 due to bad pixels etc.
    These regions are extrapolated from their boundaries.
    
    :param spec: The spectrum to correct.
    :type spec: :class:`Observation`
    :param weights: The mask which marks the bad pixels (being 0 there).
    :type weights: ndarray[nr_ord,nr_pix]
    :param orders: This array indicates which orders to correct.
    :type orders: ndarray[nr_ord_correct]
    
    :return: The corrected spectrum.
    :rtype: :class:`Observation`
    """
    for i, o in enumerate(orders):
        spec_order = spec[o]
        ind = np.where(weights[i]==0.)[0]
        
        if len(ind) > 0:
            start_pix = [ind[0]-1]
            end_pix = []
            for j in range(len(ind)-1):
                if ind[j+1] - ind[j] > 1:
                    start_pix.append(ind[j+1]-1)
                    end_pix.append(ind[j]+1)
            end_pix.append(ind[-1]+1)
            
            for j in range(len(start_pix)):
                spec_order.flux = interpolate_spec_region(spec_order.wave, spec_order.flux, start_pix[j], end_pix[j])
            
            if isinstance(spec, NormalizedObservation):
                spec[o] = spec_order.flux
            elif isinstance(spec, Observation):
                spec[o].flux = spec_order.flux
        else:
            logging.info('No zero weights: Order {}!'.format(o))
    
    return spec


def interpolate_spec_region(wave, flux, start_pix, end_pix, ext_region=5):
    """Function to interpolate bad spectral regions from their left and right 
    edges. If there are sufficiently many pixels to the order edges on both 
    sides, it is smoothly interpolated - otherwise just a linear
    interpolation is performed.
    
    :param wave: The order wavelength vector.
    :type wave: ndarray[nr_pix] or list
    :param flux: The order flux vector.
    :type flux: ndarray[nr_pix] or list
    :param start_pix: The last pixel before the bad region begins.
    :type start_pix: int
    :param end_pix: The first pixel behind the bad region.
    :type end_pix: int
    :param ext_region: How many pixels on either side to use for the 
        interpolation (optional, but really 2 are enough).
    :type ext_region: int
    
    :return: The corrected order flux vector.
    :rtype: ndarray[nr_pix]
    """
    # First check whether the extension region left or right extends past order edges
    if start_pix-ext_region+1 < 0:
        ext_region = start_pix + 1
    if end_pix+ext_region > len(wave):
        ext_region = len(wave) - end_pix
    print(ext_region)
    # If extension region now 0: Do simple linear interpolation
    if ext_region <= 1:
        flux[start_pix:end_pix+1] = np.linspace(
                flux[start_pix], flux[end_pix], end_pix-start_pix+1)
    # Otherwise interpolate more smoothly
    else:
        int_wave = wave[start_pix+1:end_pix]
        
        ind_left = np.arange(start_pix-ext_region+1, start_pix+1)
        ind_right = np.arange(end_pix, end_pix+ext_region)
        ext_wave = wave[np.concatenate((ind_left, ind_right), axis=0)]
        ext_flux = flux[np.concatenate((ind_left, ind_right), axis=0)]
        
        tck = splrep(ext_wave, ext_flux, s=0)
        int_flux = splev(int_wave, tck, der=0)
        flux[start_pix+1:end_pix] = int_flux
    
    return flux