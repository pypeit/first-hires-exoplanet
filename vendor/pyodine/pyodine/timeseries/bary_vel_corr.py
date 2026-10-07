"""
Created on Fri Oct  2 15:52:56 2020

Wrapper for the package barycorrpy by Shubham Kanodia and Jason Wright,
a Python version of Jason Eastman and Jason Wright’s IDL code [BaryCorr]
(http://astroutils.astronomy.ohio-state.edu/exofast/pro/exofast/bary/zbarycorr.pro) 
based on [Wright and Eastman (2014)](https://arxiv.org/pdf/1409.4774.pdfby).

This wrapper makes some functions easily accessible for the pyodine package.

@author: Paul Heeren
"""

import logging
import sys
import numpy as np

from barycorrpy.barycorrpy import get_BC_vel
from barycorrpy.barycorrpy import utc_tdb


def bvc_wrapper(bvc_dict, timeseries_dict, use_hip=True, z_meas=None, 
                solar=False):
    """A simple function to get the barycentric velocities for given 
    observation times, for a star and instrument, all defined in the 
    dictionaries of the CombinedResults object, as well as the correct time 
    (barycentric julian date in the barycentric dynamical time standard).
    
    :param bvc_dict: A dictionary with star and observatory info needed for the 
        barycentric correction: 'star_ra', 'star_dec', 'star_pmra', 
        'star_pmdec', 'star_rv0', 'star_name', 'instrument_lat', 
        'instrument_long', 'instrument_alt'.
    :type bvc_dict: dict
    :param timeseries_dict: A dictionary with timeseries info from the 
        CombinedResults object, containing e.g. Julian dates.
    :type timeseries_dict: dict
    :param usehip: If True, try and use the built-in hip catalogue of 
        barycorrpy to find star's coordinates (only works if star_name in
        bvc_dict is of format 'HIPxxxx'). Default is True.
    :type usehip: bool
    :param z_meas: If an array of absolute measured redshifts is handed to 
        z_meas, then the precise (multiplicative) algorithm of barycorrpy is 
        used (non-predictive). Otherwise, the less precise simple sum is 
        performed, without taking measured redshifts into account.
    :type z_meas: np.ndarray, list, tuple, or None
    :param solar: If True, return the barycentric correction for the Sun as 
        target. Defaults to False. (But even then, if the 'star_name' in 
            the bvc_dict is 'Sun', the solar BC correction will be done.)
    :type solar: bool
    
    :return: An array with the BC velocities, or the already corrected RVs (if
        z_meas was given).
    :rtype: np.ndarray
    :return: An array with corrected Julian dates (BJDs).
    :rtype: np.ndarray
    """
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    if not isinstance(z_meas, (list, np.ndarray, tuple, float, int)):
        z_meas = 0.0
    
    if solar is True or 'sun' in bvc_dict['star_name'].lower():
        logging.info('BVC for the Sun')
        
        # Calculate barycentric correction
        bcvel, warn0, stat0 = get_BC_vel(
                JDUTC = timeseries_dict['bary_date'],
                lat = bvc_dict['instrument_lat'],
                longi = bvc_dict['instrument_long'],
                alt = bvc_dict['instrument_alt'],
                ephemeris = 'de430',
                SolSystemTarget = 'Sun',
                zmeas = z_meas
                )
        
        # For solar observations, we do not compute any BJDs as barycorrpy
        # does not offer that (so just return the input times)
        bjdtdb = np.array(timeseries_dict['bary_date'])
    
    else:
        if use_hip is True and ('star_name' in bvc_dict and 'hip' in 
                                bvc_dict['star_name'].lower()):
            
            hip_nr = hip_from_name(bvc_dict['star_name'])
            
            logging.info('BVC through HIP number: {}'.format(hip_nr))
            
            # Calculate barycentric correction
            bcvel, warn0, stat0 = get_BC_vel(
                    JDUTC = timeseries_dict['bary_date'],
                    hip_id = hip_nr,
                    lat = bvc_dict['instrument_lat'],
                    longi = bvc_dict['instrument_long'],
                    alt = bvc_dict['instrument_alt'],
                    ephemeris = 'de430',
                    zmeas = z_meas
                    )
            
            # JDUTC to BJDTDB time converter
            bjdtdb, warn1, stat1 = utc_tdb.JDUTC_to_BJDTDB(
                    JDUTC = timeseries_dict['bary_date'],
                    hip_id = hip_nr,
                    lat = bvc_dict['instrument_lat'],
                    longi = bvc_dict['instrument_long'],
                    alt = bvc_dict['instrument_alt']
                    )
        
        else:
            
            ra    = bvc_dict['star_ra']
            dec   = bvc_dict['star_dec']
            pmra  = bvc_dict['star_pmra']
            pmdec = bvc_dict['star_pmdec']
            rv0   = bvc_dict['star_rv0']
            if not np.isfinite(pmra):
                pmra = 0.
            if not np.isfinite(pmdec):
                pmdec = 0.
            
            logging.info('BVC through coordinates:')
            logging.info('RA, DEC:       {}, {} (deg)'.format(ra, dec))
            logging.info('PM_RA, PM_DEC: {}, {} (mas/yr)'.format(pmra, pmdec))
            
            # Calculate barycentric correction
            bcvel, warn0, stat0 = get_BC_vel(
                    JDUTC = timeseries_dict['bary_date'],
                    ra = ra,
                    dec = dec,
                    pmra = pmra,
                    pmdec = pmdec,
                    rv = rv0,
                    lat = bvc_dict['instrument_lat'],
                    longi = bvc_dict['instrument_long'],
                    alt = bvc_dict['instrument_alt'],
                    ephemeris = 'de430',
                    zmeas = z_meas
                    )
            
            # JDUTC to BJDTDB time converter
            bjdtdb, warn1, stat1 = utc_tdb.JDUTC_to_BJDTDB(
                    JDUTC = timeseries_dict['bary_date'],
                    ra = ra,
                    dec = dec,
                    pmra = pmra,
                    pmdec = pmdec,
                    rv = rv0,
                    lat = bvc_dict['instrument_lat'],
                    longi = bvc_dict['instrument_long'],
                    alt = bvc_dict['instrument_alt']
                    )
        
    # JDUTC to JDTDB time converter
    #jdtdb = utc_tdb.JDUTC_to_JDTDB(utctime = Obstime,
    #                               leap_update=True)
    
    return bcvel, bjdtdb


def hip_from_name(star_name):
    """Small convenience function to return a HIP identifier from a star name
    
    :param star_name: The name of the star.
    :type star_name: str
    
    :return: The HIP identifier.
    :rtype: int
    """
    if 'hip' in star_name.lower():
        hip = star_name.lower().replace('hip', '').strip()
        if hip.isdecimal():
            hip = int(hip)
        else:
            hip = None
    else:
        hip = None
    
    return hip
