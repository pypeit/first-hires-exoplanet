#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 11:29:25 2020

Additional code for the pyodine package by René Tronsgaard Rasmussen,
to take the influence of tellurics into account.

@author: Paul Heeren
"""

from os.path import dirname, join
import glob
import numpy as np

from .. import components


class SimpleTellurics(components.Spectrum):
    """A simple class of a telluric mask, using the existing spectrum class
    
    :param tell_type: Telluric mask type - you can choose one of 'carmenes',
        'hitran', or 'uves'.
    :type tell_type: str, or None
    :param wave_min: The minimum wavelength of the telluric mask. If None, 
        the minimum wavelength of the chosen mask is used.
    :type wave_min: float, or None
    :param wave_max: The maximum wavelength of the telluric mask. If None, 
        the maximum wavelength of the chosen mask is used.
    :type wave_max: float, or None
    :param disp: The desired dispersion of the mask. Defaults to 0.002 A.
    :type disp: float, or None
    """
    
    def __init__(self, tell_type='carmenes', wave_min=None, wave_max=None, disp=0.002):
        
        # Which mask is desired? Load the data into a dictionary
        if tell_type=='carmenes':
            self.dict_tellurics = load_tellurics_carm()
        elif tell_type=='hitran':
            self.dict_tellurics = load_HITRAN()
        elif tell_type=='uves':
            self.dict_tellurics = load_UVES()
        else:
            raise Exception("You need to choose tell_type='carmenes', 'hitran', or 'uves'," + 
                            "but you chose {}.".format(tell_type))
        
        # Setup the wavelength and flux vectors
        self.disp = disp
        if wave_min == None:
            self.wave_min = np.min(self.dict_tellurics['wave_start'])
        if wave_max == None:
            self.wave_max = np.max(self.dict_tellurics['wave_stop'])
        else:
            self.wave_max = wave_max
        
        wave = np.arange(self.wave_min, self.wave_max, self.disp)
        flux = np.zeros(wave.shape)
        
        # Fill the wavelength and flux vectors
        for wave_start, wave_stop in zip(
                self.dict_tellurics['wave_start'],
                self.dict_tellurics['wave_stop']):
            ind = np.where((wave >= wave_start) & (wave <= wave_stop))
            flux[ind] = 1.
        super().__init__(flux, wave=wave, cont=np.zeros(flux.shape))
    
    def is_affected(self, wavelength):
        """Check if a certain wavelength is affected by a telluric feature
        of this mask
        
        :param wavelength: The wavelength(s) to check.
        :type wavelength: list, ndarray, float, or int
        
        :return: A mask with 1s where the wavelengths are affected.
        :rtype: list, or int
        """
        
        if type(wavelength) is int or type(wavelength) is float:
                ind = np.where((wavelength >= self.dict_tellurics['wave_start']) &
                               (wavelength <= self.dict_tellurics['wave_stop']))
                return len(ind[0])
        elif type(wavelength) is list or type(wavelength) is np.ndarray:
            return np.array([self.is_affected(float(w)) for w in wavelength])
        else:
            raise IndexError(type(wavelength))


def load_tellurics_carm():
    """Function to load an atlas of telluric lines from the Carmenes mask into 
    a dictionary and return that.
    
    :return: The Carmenes telluric data.
    :rtype: dict
    """
    filenames = join(dirname(__file__), 
                     'CARMENES/Wallace11_mask_0125_ext_carm_short_waveair.dat')
    tellurics = np.genfromtxt(filenames, skip_header=1, skip_footer=0, unpack=True)
    
    ind = np.where(tellurics[2,:]==1.0)
    
    dict_tellurics = {'wave_start': tellurics[1,ind[0][0::2]],
                      'wave_stop': tellurics[1,ind[0][1::2]],
                      'flux': np.array([1.0] * len(tellurics[1,ind[0][0::2]]))
                     }
    return dict_tellurics
                

def load_HITRAN():
    """Function to load an atlas of telluric lines from the HITRAN file into a 
    dictionary and return that.
    
    :return: The HITRAN telluric data.
    :rtype: dict
    """
    
    filenames = join(dirname(__file__), 'HITRAN/HITRAN_mask.dat')
    tellurics = np.genfromtxt(filenames, skip_header=1, skip_footer=0, unpack=True)
    dict_tellurics = {'wave_start': tellurics[2,:],
                      'wave_stop': tellurics[3,:],
                      'flux': tellurics[4,:]
                     }
    return dict_tellurics


def load_UVES():
    """Function to load an atlas of telluric lines from the UVES data into a 
    dictionary and return that.
    
    :return: The UVES telluric data.
    :rtype: dict
    """
    _ref_path = join(dirname(__file__), 'UVES/UVES_ident')
    #_ref_file = 'telluric_lines.txt'
    _ref_file = 'gident_*.dat'
    
    skip_header = 3
    skip_footer = 1
    
    filenames = glob.glob(join(_ref_path, _ref_file))
    tellurics = np.array([[],[],[],[],[]])
    for filename in filenames:
        tellurics = np.concatenate((tellurics, np.genfromtxt(filename, skip_header=skip_header, 
                                       skip_footer=skip_footer, unpack=True)), axis=1)
    
    dict_tellurics = {'wave_start': tellurics[1]-0.5*tellurics[3],
                      'wave_stop': tellurics[1]+0.5*tellurics[3],
                      'flux': tellurics[4]
                     }
    return dict_tellurics