#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 15:40:23 2022

@author: paul
"""

from astropy.io import fits
import numpy as np
import barycorrpy
import os
from astropy.coordinates import SkyCoord
from astropy import units as u


def filesplit(specname, wavename):
    
    # First get the flux
    header = fits.getheader(specname)
    data = fits.getdata(specname)
    
    flux = np.array([data[i]['flux'] for i in range(len(data))])
    
    # Then get the wavelengths
    data_wave = fits.getdata(wavename)
    
    wave = np.array([data_wave[i]['wavelength'] for i in range(len(data_wave))])
    
    
    # Now split the orders of the two fibers
    
    wave2 = wave[2::2]
    wave1 = wave[3::2]
    
    if len(flux) == 136:
        fiber2 = flux[2::2]
        fiber1 = flux[3::2]
    elif len(flux) == 135:
        fiber2 = flux[0::2]
        fiber1 = flux[1::2]
        #wave2 = wave[0::2]
        #wave1 = wave[1::2]
    else:
        print('Length of order vector was {}'.format(len(flux)))
    
    coordinates = SkyCoord(
            header['OBJ-RA'].strip() + ' ' + header['OBJ-DEC'].strip(),
            unit=(u.hourangle, u.deg)
        )
    
    bcvel, warn0, stat0 = barycorrpy.barycorrpy.get_BC_vel(
                JDUTC = header['JD-MID'],
                lat = -27.79793,
                longi = 151.85555,
                alt = 800.,
                ra = coordinates.ra.deg,
                dec = coordinates.dec.deg,
                ephemeris = 'de430'
                )
    
    header['BVC'] = bcvel[0] / 1000.
    
    header1 = header
    header2 = header
    
    header1['FIBER'] = 1
    header1['TELESCOP'] = 'T1'
    header2['FIBER'] = 2
    header2['TELESCOP'] = 'T2'
    header1['WAV_FILE'] = wavename
    header2['WAV_FILE'] = wavename
    
    cont1 = np.ones((fiber1.shape))
    cont2 = np.ones((fiber2.shape))
    
    hdu = fits.PrimaryHDU(np.stack((fiber1, cont1, wave1)), header=header1)
    hdu.writeto(specname[:-12]+'fib1.fits')
    
    hdu = fits.PrimaryHDU(np.stack((fiber2, cont2, wave2)), header=header2)
    hdu.writeto(specname[:-12]+'fib2.fits')