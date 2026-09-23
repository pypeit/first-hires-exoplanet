import numpy as np
from astropy.time import TimeDelta
import logging
import sys

import pyodine.comp_io as comp_io

__all__ = ["Spectrum", "MultiOrderSpectrum", "Observation", "Instrument", "Star",
           "NormalizedObservation", "SummedObservation", "Chunk", "ChunkArray",
           "TemplateChunk"]


class NoDataError(BaseException):
    """Use this Exception class to indicate missing data
    """
    pass  # FIXME: Is this it?


class DataMismatchError(BaseException):
    """Use this Exception class to indicate that the data is incompatible
    (e.g. trying to add a non-iodine spectrum with an iodine spectrum)
    """
    pass  # FIXME: Is this it?


class Spectrum:
    """A 1D spectrum, i.e. a set of flux values corresponding to pixels or 
    wavelengths
    
    This class serves as most basic parent class to all spectrum objects.
    
    :param flux: Flux values of the spectrum.
    :type flux: ndarray[nr_pix]
    :param wave: Wavelength values of the spectrum.
    :type wave: ndarray[nr_pix], or None
    :param cont: Continuum values of the spectrum.
    :type cont: ndarray[nr_pix], or None
    
    """
    def __init__(self, flux, wave=None, cont=None):
        if not any(flux):
            raise NoDataError('Invalid flux vector!')
        self.flux = flux
        self.wave = wave
        self.cont = cont

    def __len__(self):
        """The dedicated length-method
        
        :return: Length of the flux vector.
        :rtype: int
        """
        return len(self.flux)

    def __getitem__(self, pixels):
        """The dedicated get-method
        
        :param pixels: The pixel indices to return from the spectrum.
        :type pixels: int, list, ndarray, slice
        
        :return: A spectrum of the desired pixel indices.
        :rtype: :class:`Spectrum`
            
        """
        flux = self.flux[pixels]
        wave = self.wave[pixels] if self.wave is not None else None
        cont = self.cont[pixels] if self.cont is not None else None
        return Spectrum(flux, wave=wave, cont=cont)

    def check_wavelength_range(self, wave_start, wave_stop):
        """Check the fraction of wavelength range as supplied by the input
        arguments covered by the data.
        
        :param wave_start: Starting wavelength.
        :type wave_start: float
        :param wave_stop: Stopping wavelength.
        :type wave_stop: float
        
        :return: A value between 0.0 and 1.0, telling how big a fraction of 
            the wavelength range [wave_start:wave_stop] is covered by data.
        :rtype: float
        
        """
        # Check that wavelength data is present
        if self.wave is None:
            raise NoDataError('No wavelength data')
        # Check that parameters are valid
        if wave_start >= wave_stop:
            raise ValueError('Bad input! (wave_start >= wave_stop)')
        # Case 1: Requested range completely outside data range
        if wave_start >= self.wave[-1] or wave_stop <= self.wave[0]:
            return 0.0
        # Case 2: Requested range completely covered by data range
        elif wave_start >= self.wave[0] and wave_stop <= self.wave[-1]:
            return 1.0
        # Case 3: Requested range partially within data range
        else:
            wave_max = np.min([wave_stop, self.wave[-1]])
            wave_min = np.max([wave_start, self.wave[0]])
            return (wave_max - wave_min) / (wave_stop - wave_start)

    def get_wavelength_range(self, wave_start, wave_stop, require=None):
        """Return a wavelength inverval of the spectrum, defined by wave_start
        and wave_stop. If require='full' ensure that the full interval is
        covered (one pixel outside in each end).
        
        :param wave_start: Starting wavelength.
        :type wave_start: float
        :param wave_stop: Stopping wavelength.
        :type wave_stop: float
        :param require: If set to 'full', make sure that the whole wavelength 
            range is covered by the data (error otherwise).
        :type require: str, or None
        
        :return: The spectrum in the wavelength range.
        :rtype: :class:`Spectrum`
        
        """
        # Measure the overlap between requested range and data range
        # (and throw error if no wavelength data present)
        coverage = self.check_wavelength_range(wave_start, wave_stop)
        # Decide if coverage is acceptable
        if require == 'full' and not np.isclose(coverage, 1.0):
            raise NoDataError('Not enough data to cover wavelength range!')
        elif np.isclose(coverage, 0.0):
            raise NoDataError('No data in requested wavelength range!')
        # Determine the first pixel of the range
        if wave_start >= self.wave[0]:
            first = np.searchsorted(self.wave, wave_start, side='right') - 1
        else:
            first = 0
        # Determine the last pixel of the range
        if wave_stop <= self.wave[-1]:
            last = np.searchsorted(self.wave, wave_stop, side='left')
        else:
            last = len(self.wave)
        # Return wavelength range
        return self[first:last + 1]

    def __str__(self):
        """The dedicated string-method
        
        :return: A string with information about the data.
        :rtype: str
        
        """
        if self.wave is not None:
            return '<Spectrum ({} pixels, {:.4f}-{:.4f} Å>'.format(
                len(self),
                self.wave[0], self.wave[-1]
            )
        else:
            return '<Spectrum ({} pixels)>'.format(len(self))
    
    def compute_weight(self, weight_type='flat', rel_noise=0.008):
        """Compute and return pixel weights for the spectrum
        
        If weight_type='inverse' is used, the pixel weights are estimated from
        the flux values, with lower flux (-> absorption lines) receiving higher 
        weights. This has been included in analogy to the dop-code by D. Fisher,
        but it is not well-tested here!
        
        :param weight_type: The type of weights to compute. Either 'flat' for 
            flat weights (all ones, default), or 'inverse' for inversely 
            weighted by flux (as in dop-code, D. Fisher, Yale University).
        :type weight_type: str
        :param rel_noise: The relative noise as measured in a flatfield. Only
            required if using weight_type='inverse'.
        :type rel_noise: float
        
        :return: The computed weights array.
        :rtype: ndarray[nr_pix]
        
        """
        if weight_type == 'flat':
            return np.ones(self.flux.shape)
        elif weight_type == 'inverse':
            return 1./(self.flux * (1. + self.flux * rel_noise**2))
        else:
            raise NotImplementedError('Choose one of: "flat", "inverse"')


class MultiOrderSpectrum:
    """A spectrum with multiple orders, represented as a list of 1D 
    :class:'Spectrum' objects
    
    Base class for Observation and StellarTemplate. Final subclasses must 
    implement the __getitem__() method!
    """
    nord = None
    npix = None
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')

    @property
    def orders(self):
        return np.arange(self.nord, dtype='int')

    def __getitem__(self, order) -> Spectrum:
        """Return one spectral order"""
        raise NotImplementedError

    def __len__(self):
        """The dedicated length-method
        
        Return:
            int: The number of orders.
        """
        return self.nord

    def check_wavelength_range(self, wave_start, wave_stop):
        """Find the order with the best coverage of the wavelength range defined
        by wave_start and wave_stop. If no such order is found, raise a 
        NoDataError. Otherwise return the
        order index and coverage as a tuple.
        
        :param wave_start: Starting wavelength.
        :type wave_start: float
        :param wave_stop: Stopping wavelength.
        :type wave_stop: float
        
        :return: The index of the order best covering the wavelength range.
        :rtype: int
        :return: A value between 0.0 and 1.0, telling how big a fraction of 
            the wavelength range [wave_start:wave_stop] is covered by data.
        :rtype: float
        
        """
        selected_order = None
        best_coverage = 0.
        # Loop through orders and look for requested wavelengths
        for i in self.orders:
            try:
                coverage = self[i].check_wavelength_range(wave_start, wave_stop)
                # Select order if coverage better than previous orders
                if coverage > best_coverage:
                    best_coverage = coverage
                    selected_order = i
            except NoDataError:
                pass
        if selected_order is None:
            logging.error(NoDataError(
                'Could not find wavelength range {}-{} Å'.format(
                    wave_start, wave_stop
                )))
        return selected_order, best_coverage

    def get_wavelength_range(self, wave_start, wave_stop, require=None):
        """Loop through orders and search for a given wavelength range.
        Return a wavelength inverval of the spectrum, defined by wave_start
        and wave_stop.
        If require='full' ensure that the full interval is
        covered (one pixel outside in each end).
        If multiple results, return the one with best coverage.
        
        :param wave_start: Starting wavelength.
        :type wave_start: float
        :param wave_stop: Stopping wavelength.
        :type wave_stop: float
        :param require: If set to 'full', make sure that the whole wavelength 
            range is covered by the data (error otherwise).
        :type require: str, or None
        
        :return: The spectrum in the wavelength range.
        :rtype: :class:`Spectrum`
        
        """
        selected_order, best_coverage = \
            self.check_wavelength_range(wave_start, wave_stop)
        return self[selected_order].get_wavelength_range(
            wave_start, wave_stop, require=require
        )
    
    def compute_weight(self, weight_type='flat', rel_noise=0.008):
        """Loop through orders and compute pixel weights for the spectrum
        """
        weight = []
        for i in self.orders:
            weight.append(self[i].compute_weight(weight_type, rel_noise))
        return np.array(weight)


class IodineAtlas(Spectrum):
    """A high-resolution FTS spectrum of the iodine absorption lines
    """
    pass


class Observation(MultiOrderSpectrum):
    """A cross-dispersed spectrum of a specific star at specific time, as seen
    from a specific instrument
    """

    instrument = None
    star = None

    orig_filename = None    # Absolute path to original file
    orig_header = None      # Original FITS header

    iodine_cell_id = None
    exp_time = None         # Exposure time in seconds
    flux_level = None       # An instrument-specific measure of the flux level
    # FIXME: Flux level defined like this would also depend on the star and the BVC?
    gain = 1.0              # TODO
    readout_noise = None    # TODO
    dark_current = None     # TODO

    # Timing
    time_start = None       # Start time of the observation <astropy.time.Time>
    time_weighted = None    # Weighted mid-time of the obs. <astropy.time.Time>

    bary_date = None  # Mid-time as Barycentric Reduced Julian Date (BJD - 2400000.0)
    bary_vel_corr = None    # Barycentric velocity correction (km/s)

    @property
    def time_mid(self):
        return self.time_start + TimeDelta(0.5 * self.exp_time, format='sec')

    def save(self, filename, data, header=None):
        """Save observation in fits format
        
        :param filename: The filename to save the observation to.
        :type filename: str
        :param data: The data array to save. By leaving it as a required input
            argument here, this needs to be defined downstream in child classes.
        :type data: ndarray or list
        :param header: An instance of the original fits header when the data
            was loaded from file, or a dictionary. If None, try using the
            property orig_header.
        :type header: :class:`fits.header`, dict, or None
        """
        if not isinstance(filename, str):
            raise ValueError('No output filename as type string given!')
        
        if not isinstance(data, (np.ndarray, list)):
            raise ValueError('No data as type ndarray or list given!')
        
        if header == None:
            if hasattr(self, 'orig_header'):
                header = self.orig_header
        
        comp_io.save_fits(filename, data, add_header=header)


class Instrument:
    """A generic class to represent an instrument
    
    :param name: The name of the instrument.
    :type name: str
    :param longitude: Longitude in degrees.
    :type longitude: float, or None
    :param latitude: Latitude in degrees.
    :type latitude: float, or None
    :param altitude: Altitude in meters.
    :type altitude: float, or None
    """
    def __init__(self, name, longitude=None, latitude=None, altitude=None):
        self.name = name
        self.longitude = longitude
        self.latitude = latitude
        self.altitude = altitude


class Star:
    """Generic representation of a stellar target
    
    :param name: The name of the star.
    :type name: str
    :param coordinates: The sky coordinates of the star.
    :type coordinates: :class:`SkyCoord`, or None
    :param proper_motion: Proper motion in (RA, DEC) in mas/year.
    :type proper_motion: tuple(float,float), or None
    
    """
    def __init__(self, name, coordinates=None, proper_motion=(None, None)):
        self.name = name  # Name of the target (e.g. `Sigma Draconis`)
        self.coordinates = coordinates      # As a SkyCoord object (astropy.coordinates)
        self.proper_motion = proper_motion  # (pm_ra, pm_dec) in mas/year


class NormalizedObservation(Observation):
    """Wraps an observation object and overrides the flux with its own
    normalized value. Everything else is loaded from the provided
    observation object.
    
    Args:
        observation ()
    """
    # List of attributes that cannot be overwritten by orig_obs
    __attrs = ['orders', 'orig_obs', '_flux', '_normalized_orders', '__dir__',
               'save_norm']

    def __init__(self, observation):
        """Initialize the class by providing the original observation"""
        self.orig_obs = observation
        self._flux = [None] * len(observation)
        self._normalized_orders = []

    @property
    def orders(self):
        """Return a list of orders that have been normalized"""
        return sorted(self._normalized_orders)

    def __getattribute__(self, item):
        """
            Except for normalized flux and list of normalized orders, get all
            attributes from the the wrapped observation.
        """
        if item in NormalizedObservation.__attrs:
            return object.__getattribute__(self, item)
        return self.orig_obs.__getattribute__(item)

    def __dir__(self):
        """Set the list of attributes"""
        return NormalizedObservation.__attrs + dir(self.orig_obs)

    def __getitem__(self, order):
        """
            Fetch one or more orders from the original observation and replace
            the flux with the normalized one.
        """
        # Return one order
        if type(order) is int or hasattr(order, '__int__'):
            spec = self.orig_obs[order]
            return Spectrum(self._flux[order], spec.wave, spec.cont)
        elif type(order) is list:
            return [self.__getitem__(int(i)) for i in order]  # Return MultiOrderSpectrum instead?
        elif type(order) is slice:
            return self.__getitem__([int(i) for i in np.arange(self.nord)[order]])
        else:
            raise IndexError(type(order))

    def __setitem__(self, order, flux):
        """
            Set the flux of order i using this syntax:

            norm_obs[i] = flux
        """
        if len(flux) == self.orig_obs.npix:
            self._flux[order] = flux
            self._normalized_orders.append(order)
        else:
            raise ValueError('Flux vector length does not match observation')
    
    def save_norm(self, filename):
        """Save normalized observation (with original flux, wave and cont) in 
        fits format
        
        :param filename: The filename to save the observation to.
        :type filename: str
        """
        
        # First get the normalized flux, and the flux, wave and cont data from 
        # the original observation (all only for normalized orders)
        norm_flux = np.array([self._flux[i] for i in self._normalized_orders])
        
        orig_flux = np.array([self.orig_obs[i].flux for i in self._normalized_orders])
        
        if isinstance(self.orig_obs[self._normalized_orders[0]].wave, (np.ndarray, list)):
            wave = np.array([self.orig_obs[i].wave for i in self._normalized_orders])
        else:
            wave = np.zeros(orig_flux.shape)
        
        if isinstance(self.orig_obs[self._normalized_orders[0]].cont, (np.ndarray, list)):
            cont = np.array([self.orig_obs[i].cont for i in self._normalized_orders])
        else:
            cont = np.zeros(orig_flux.shape)
        
        # Now create the data array which should be saved to fits
        data = np.array([norm_flux, orig_flux, wave, cont])
        
        # Add the original header (if existent)
        if hasattr(self.orig_obs, 'orig_header'):
            header = self.orig_obs.orig_header
        else:
            header = None
        
        # And save, using the method from the parent class
        super().save(filename, data, header=header)


class SummedObservation(Observation):

    def __init__(self, *observations):
        """
            Initialize with one or more observations - read all properties
            except flux, exptime, (???) from this one.
        """
        self.npix = observations[0].npix
        self.nord = observations[0].nord
        self.star = observations[0].star
        self.instrument = observations[0].instrument
        self.orig_header = observations[0].orig_header
        self.bary_date = observations[0].bary_date
        self.bary_vel_corr = observations[0].bary_vel_corr
        self.orig_filename = observations[0].orig_filename
        self.all_filenames = [obs.orig_filename for obs in observations]

        self.observations = []
        self._flux = {}  # Dict with the summed flux in each order
        self.exp_time = 0.0  # Sum of exposure times
        self.time_start = observations[0].time_start
        self.iodine_in_spectrum = observations[0].iodine_in_spectrum
        self.iodine_cell_id = observations[0].iodine_cell_id
        # TODO: Handle read-out noise etc.

        for i in observations[0].orders:
            self._flux[i] = np.zeros(self.npix)

        self.add(*observations)

    def __getitem__(self, order):
        """
            Fetch one or more orders from the original observation and replace
            the flux with the normalized one.
        """
        # Return one order
        if type(order) is int or hasattr(order, '__int__'):
            flux = self._flux[order]
            wave = self.observations[0][order].wave
            cont = self.observations[0][order].cont
            return Spectrum(flux, wave=wave, cont=cont)
        elif type(order) is list:
            return [self.__getitem__(int(i)) for i in order]  # Return MultiOrderSpectrum instead?
        elif type(order) is slice:
            return self.__getitem__([int(i) for i in np.arange(self.nord)[order]])
        else:
            raise IndexError(type(order))

    @property
    def nobs(self):
        """Number of observations"""
        return len(self.observations)

    @property
    def orders(self):
        """List order numbers"""
        return list(self._flux.keys())

    def add(self, *observations):
        """Add one or more observation to the sum of observations"""
        for obs in observations:
            self.observations += [obs]
            # Add exposure time
            self.exp_time += obs.exp_time
            # Set start time
            if obs.time_start < self.time_start:
                self.time_start = obs.time_start
            # Check that iodine cell configuration is the same
            if obs.iodine_in_spectrum != self.iodine_in_spectrum or \
                    obs.iodine_cell_id != self.iodine_cell_id:
                raise DataMismatchError
            # Add the flux from each order
            for i in obs.orders:
                self._flux[i] += obs[i].flux
    
    def save(self, filename):
        """Save summed observation in fits format
        
        :param filename: The filename to save the observation to.
        :type filename: str
        """
        
        # First get the flux, wave and cont data
        orders = list(self._flux.keys())
        
        flux = np.array([self[i].flux for i in orders])
        
        if isinstance(self[orders[0]].wave, (np.ndarray, list)):
            wave = np.array([self[i].wave for i in self.orders])
        else:
            wave = np.zeros(flux.shape)
        
        if isinstance(self[orders[0]].cont, (np.ndarray, list)):
            cont = np.array([self[i].cont for i in self.orders])
        else:
            cont = np.zeros(flux.shape)
        
        # Now create the data array which should be saved to fits
        data = np.array([flux, wave, cont])
        
        # And save, using the method from the parent class
        super().save(filename, data, header=self.orig_header)


class Chunk(Spectrum):
    """A subsection of an observation, defined by order and pixels
    
    This data object will be used in the fitting procedure.
    
    :param observation: The observation from which the chunk created.
    :type observation: :class:`Observation`
    :param order: The order that the chunk sits in.
    :type order: int
    :param pixels: The pixels covered by the chunk.
    :type pixels: ndarray[chunk_size]
    :param padding: The number of pixels used to extend the chunk with 
        property 'chunk.padded' (necessary in deconvolution etc.).
    :type padding: int
    
    """

    def __init__(self, observation, order, pixels, padding=0):
        spec = observation[order][pixels]
        super().__init__(spec.flux, wave=spec.wave, cont=spec.cont)
        self.observation = observation
        self.order = order
        self.abspix = pixels
        self.padding = padding

    @property
    def pix(self):
        """A pixel vector for the chunk, centered around zero
        
        :return: The pixel vector.
        :rtype: ndarray[chunk_size]
        
        """
        n = len(self)
        return np.arange(-(n // 2), n - n // 2)

    @property
    def padded(self):
        """The chunk spectrum with padding on either side included
        
        :return: Chunk including padding.
        :rtype: :class:`Chunk`
        
        """
        if self.padding == 0:
            return self
        first = self.abspix[0] - self.padding
        last = self.abspix[-1] + self.padding
        pixels = np.arange(first, last + 1, dtype='int')
        return Chunk(self.observation, self.order, pixels)

    def __str__(self):
        """The dedicated string-method
        
        :return: Information about the chunk.
        :rtype: str
        
        """
        return '<Chunk (order:{} ; pixels:{}-{})>'.format(self.order, *self.abspix[[0, -1]])


class ChunkArray(list):
    """Behaves like a list of chunks, with the added ability to filter
    specific orders

        Example:
            $ chunkarr = ChunkArray(chunklist)
            $ chunkarr[7]  # Returns chunk number 7
            $ chunkarr.orders  # Returns a list of unique orders
            $ chunkarr.get_order(22)  # Returns all chunks from order 22
    """

    @property
    def orders(self):
        """Return the order numbers contained in the chunk array as ndarray
        
        :return: The unique order numbers.
        :rtype: ndarray
        """
        return np.unique([chunk.order for chunk in self])

    def get_order(self, order) -> list:
        """Return chunks within order
        
        :param order: The order to use.
        :type order: int
        
        :return: The list of chunks within the supplied order.
        :rtype: list[:class:`Chunk`]
        
        """
        return [chunk for chunk in self if chunk.order == order]

    def get_order_indices(self, order) -> list:
        """Return indices of chunks within order
        
        :param order: The order to use.
        :type order: int
        
        :return: The list of chunk indices within the supplied order.
        :rtype: list[int] 
        """
        return [i for i in range(len(self)) if self[i].order == order]


class TemplateChunk(Spectrum):
    """A chunk of a deconvolved template
    
    This is used in the :class:`StellarTemplate_Chunked`. Keyword 'padding' 
    defines the number of pixels used to extend the chunk with property 
    `chunk.padded` (necessary in deconvolution etc.).
    
    :param flux: The flux values of the template.
    :type flux: ndarray[nr_pix]
    :param wave: The wavelength values of the template.
    :type wave: ndarray[nr_pix]
    :param pixel:  The pixel vector centered around 0.
    :type pixel: ndarray[nr_pix]
    :param w0: The zero point of the wavelength solution used to create
        this template chunk.
    :type w0: float
    :param w1: The dispersion of the wavelength solution used to create this 
        template chunk.
    :type w1: float
    :param order: The order of the chunk.
    :type order: int
    :param pix0: The starting pixel of the chunk within the original template 
        observation order.
    :type pix0: int
    :param weight: The weight of the chunk.
    :type weight: float
    
    """

    def __init__(self, flux, wave, pixel, w0, w1, order, pix0, weight):
        super().__init__(flux, wave=wave, cont=None)
        self.pixel = pixel
        self.w0 = w0
        self.w1 = w1
        self.order = order
        self.pix0 = pix0
        self.weight = weight

    @property
    def pix(self):
        """A pixel vector for the chunk, centered around zero
        
        :return: The pixel vector.
        :rtype: ndarray[nr_pix]
        
        """
        n = len(self)
        return np.arange(-(n // 2), n - n // 2)
