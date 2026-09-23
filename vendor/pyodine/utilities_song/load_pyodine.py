from os.path import splitext, abspath

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits as pyfits
from astropy.time import Time, TimeDelta
import h5py
from pyodine import components

from utilities_song import conf

__all__ = ["IodineTemplate", "ObservationWrapper"]


class IodineTemplate(components.IodineAtlas):
    """The iodine template class to be used in the modelling
    
    :param iodine_cell_id: The iodine cell ID to identify the I2 template
        spectrum by in the :ref:`overview_utilities_conf`, or the direct pathname to the I2
        template spectrum.
    :type iodine_cell_id: int or str
    """
    def __init__(self, iodine_cell):
        if not isinstance(iodine_cell, (int,str)):
            raise KeyError('Argument "iodine_cell" must be either int or string!')
        elif isinstance(iodine_cell, int):
            if iodine_cell in conf.my_iodine_atlases.keys():
                self.orig_filename = conf.my_iodine_atlases[iodine_cell]
            else:
                raise ValueError('Unknown iodine_cell ID!')
        elif isinstance(iodine_cell, str):
            self.orig_filename = iodine_cell
        
        with h5py.File(self.orig_filename, 'r') as h:
            flux = h['flux_normalized'][()]
            wave = h['wavelength_air'][()]    # originally: wavelength_air
        super().__init__(flux, wave)


class ObservationWrapper(components.Observation):
    """A wrapper for the representation of SONG observation spectra, based
    on the parent class :class:`pyodine.components.Observation`
    
    :param filename: The filename of the observation to load.
    :type filename: str
    :param instrument: The instrument used to obtain the observation. If None,
        the information is drawn from the Fits-header (default).
    :type instrument: :class:`components.Instrument`
    :param star: The star of the observation. If None, the information is 
        drawn from the Fits-header (default).
    :type star: :class:`components.Star`
    """

    # Custom properties
    _spec = None    # Internal storage of spectral flux
    _wave = None    # Internal storage of wavelength solution
    _cont = None    # Internal storage of extracted continuum

    def __init__(self, filename, instrument=None, star=None):
        flux, wave, cont, header = load_file(filename)

        self._flux = flux
        self._wave = wave
        self._cont = cont
        
        """ Weights added. Using this formula from dop code for now
            (the value of 0.008 is the flatfield noise - should be changed)
        if weight is None:# or len(weight) is not len(self.flux):
            #self._weight = (1./self._flux) / (1. + self._flux * 0.008**2)
            self._weight = np.ones(self._flux.shape)
        else:
            self._weight = weight"""

        self.nord = flux.shape[0]
        self.npix = flux.shape[1]

        self.orig_header = header
        self.orig_filename = abspath(filename)

        self.instrument = instrument or get_instrument(header)
        self.star = star or get_star(header)
        self.iodine_in_spectrum, self.iodine_cell_id = check_iodine_cell(header)

        # Camera details
        self.exp_time = get_exposuretime(header, self.instrument)  # or_none(header, 'EXPOSURE')
        self.flux_level = None      # FIXME: Define a flux measure
        self.gain = None            # FIXME: Not in header
        self.readout_noise = None   # FIXME: Not in header
        self.dark_current = None    # FIXME: Not in header

        # Timing
        self.time_start = Time(header['DATE-OBS'].strip(), format='isot', scale='utc')
        self.time_weighted = None

        self.bary_date = or_none(header, 'JD-MID')
        self.bary_vel_corr = or_none(header, 'BVC') * 1000.     # km/s in SONG header
        #self.topo_bary_factor = or_none(header, 'BVCFACT')
        #self.mjd_corr = or_none(header, 'MID-JD')#'MBJD')
        #self.moon_vel = or_none(header, 'MOONVEL') * 1000.  # convert to m/s
        # TODO: Implement flux check
        # TODO: Re-calculate BVC

    def __getitem__(self, order) -> components.Spectrum:
        """Return one or more spectral orders
        
        :param order: The order(s) of the spectrum to return.
        :type order: int, list, ndarray, slice
        
        :return: The desired order(s).
        :rtype: :class:`Spectrum` or list[:class:`Spectrum`]
        """
        # Return one order
        if type(order) is int or hasattr(order, '__int__'):
            flux = self._flux[order]
            wave = self._wave[order]
            cont = self._cont[order]
            #weight = self._weight[order]
            return components.Spectrum(flux, wave=wave, cont=cont)#, weight=weight)
        elif isinstance(order, (list, np.ndarray)):
            return [self.__getitem__(int(i)) for i in order]  # Return MultiOrderSpectrum instead?
        elif type(order) is slice:
            return self.__getitem__([int(i) for i in np.arange(self.nord)[order]])
        else:
            raise IndexError(type(order))


def load_file(filename) -> components.Observation:
    """A convenience function to load observation data from file
    
    :param filename: The filename of the observation to load.
    :type filename: str
    
    :return: The flux of the observation spectrum.
    :rtype: ndarray
    :return: The wavelengths of the observation spectrum.
    :rtype: ndarray
    :return: The continuum flux of the observation spectrum.
    :rtype: ndarray
    :return: The Fits-header.
    :rtype: :class:`fits.Header`
    """
    try:
        ext = splitext(filename)[1]
        if ext == '.fits':
            # Load the file
            h = pyfits.open(filename)
            header = h[0].header
            # Prepare data
            if 'OPT_DONE' in header.keys():
                flux = h[0].data[0] if header['OPT_DONE'] == 'TRUE' else h[0].data[1]
            else:
                flux = h[0].data[0] if sum(h[0].data[0].flatten()) > 0 else h[0].data[1]
            cont = h[0].data[2]
            wave = h[0].data[3]
            #weight = None

            h.close()
            # TODO: Check for `songwriter` signature
            return flux, wave, cont, header
        else:
            # Unsupported file format
            raise TypeError('Unsupported file format (%s)' % ext)
    except IOError:
        print('Could not open file %s' % filename)
    except TypeError as e:
        print(e.args[0])


def get_star(header) -> components.Star:
    """Create a star object based on header data
    
    :param header: The Fits-header.
    :type header: :class:`fits.Header`
    
    :return: The star object.
    :rtype: :class:`Star`
    """
    # TODO: Load stars from some kind of catalog based on name instead?

    name = or_none(header, 'OBJECT')
    try:
        coordinates = SkyCoord(
            header['OBJ-RA'].strip() + ' ' + header['OBJ-DEC'].strip(),
            unit=(u.hourangle, u.deg)
        )
    except Exception as e:
        # TODO: Log this event
        coordinates = None
    # Get the proper motion vector
    try:
        proper_motion = (header['S-PM-RA'], header['S-PM-DEC'])
    except Exception as e:
        # TODO: Log this event
        proper_motion = (None, None)

    return components.Star(name, coordinates=coordinates, proper_motion=proper_motion)


def get_instrument(header) -> components.Instrument:
    """Determine the instrument from the header and return Instrument object
    
    :param header: The Fits-header.
    :type header: :class:`fits.Header`
    
    :return: The instrument object.
    :rtype: :class:`Instrument`
    """
    if 'TELESCOP' in header:
        if 'Node 1' in header['TELESCOP'] and 'Spectrograph' in header['INSTRUM']:
            return conf.my_instruments['song_1']
        elif 'Node 2' in header['TELESCOP'] and 'Spectrograph' in header['INSTRUM']:
            return conf.my_instruments['song_2']
        elif 'Waltz' in header['TELESCOP']:
            return conf.my_instruments['waltz']
        elif 'Hamilton' in header['INSTRUME'] or 'HAMILTON' in header['PROGRAM'].upper() or \
        '3M-COUDE' in header['TELESCOP'].upper() or '3M-CAT' in header['PROGRAM'].upper():
            return conf.my_instruments['lick']
    else:
        if 'NEWCAM' in header['PROGRAM'] and 'hamcat' in header['VERSION']:
            return conf.my_instruments['lick']
        # TODO: Log this event
        raise TypeError('Could not determine instrument')


def check_iodine_cell(header):
    """Check the position and state of the I2 cell during the observation
    
    :param header: The Fits-header.
    :type header: :class:`fits.Header`
    
    :return: Whether or not the I2 cell was in the light path.
    :rtype: bool
    :return: The ID of the used I2 cell.
    :rtype: int, or None
    """
    # If the IODID keyword is set, we should be safe
    if 'IODID' in header.keys() and header['I2POS'] != 2:
        iodine_in_spectrum = True
        iodine_cell_id = header['IODID']
    # Otherwise, let's make a qualified guess based on the I2POS keyword
    else:
        # TODO: Log this event
        # Position 3 corresponds to id=1
        if header['I2POS'] == 3:
            iodine_in_spectrum = True
            iodine_cell_id = 1
        elif header['I2POS'] == 1:
            iodine_in_spectrum = True
            iodine_cell_id = 2
        else:
            # Position 2 lets the light pass through
            iodine_in_spectrum = False
            iodine_cell_id = None
    return iodine_in_spectrum, iodine_cell_id


def or_none(header, key, fallback_value=None):
    """A convenience function to prevent non-existent Fits-header cards from
    throwing up errors
    
    :param header: The Fits-header.
    :type header: :class:`fits.Header`
    :param key: The keyword of the header card of interest.
    :type key: str
    :param fallback_value: What to return if the header card does not exist
        (default: None).
    :type fallback_value: str, int, float, or None
    
    :return: The header card or the 'fallback_value'.
    :rtype: str, int, float, or None
    """
    try:
        return header[key]
    except KeyError:
        # TODO: Log this event
        return fallback_value


def get_exposuretime(header, instrument):
    """Get the exposure time from the fits header (this extra function is 
    neccessary to make old Lick spectra work smoothly)
    
    """
    if 'SONG' in instrument.name:
        return or_none(header, 'EXPTIME')
    elif 'EXPOSURE' in header and 'Lick' in instrument.name:
        # sometimes it's in millisecs - let's try and catch most of these times
        if header['EXPOSURE'] > 3600.:
            return header['EXPOSURE'] / 1000.
        else:
            return header['EXPOSURE']
    elif 'EXPTIME' in header and 'Lick' in instrument.name:
        return header['EXPTIME']
    else:
        return None


def get_barytime(header, instrument):
    """Get the date and time of the weighted midpoint from the fits header
    (this extra function is neccessary to make old Lick spectra work smoothly)\
    
    """
    if 'SONG' in instrument.name:
        return or_none(header, 'BJD-MID')
    elif 'Lick' in instrument.name:
        # in Lick the MID-time is only given in hrs, mins, secs
        # so we create the MID-JD manually here
        date  = header['DATE-OBS'].strip()[:10]
        stime = header['MP-START'].strip()
        mtime = header['MP-MID'].strip()
        time_start = Time(date+'T'+stime, format='isot', scale='utc')
        bary_date  = Time(date+'T'+mtime, format='isot', scale='utc')
        # check whether time_weighted is on the next full day
        if time_start > bary_date:
            bary_date += TimeDelta(1., format='jd')
        return bary_date.jd

