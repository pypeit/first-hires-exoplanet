import numpy as np
import sys
import logging

from ..components import Spectrum, MultiOrderSpectrum, NormalizedObservation
from ..lib.misc import rebin
from ..reference import load_reference

_c = 299792458  # Speed of light in m/s


class Normalizer:
    """Abstract base class
    """
    def normalize_obs(self, observation, velocity, orders=None):
        """Normalize an Observation (multiorder) or a list of orders within
        """
        raise NotImplementedError

    def normalize_single(self, spectrum, velocity, order):
        """Normalize a single spectrum
        """
        raise NotImplementedError


class SimpleNormalizer(Normalizer):
    """A simple Normalizer
    
    :param velocity: ?
    :type velocity: float, or None
    :param reference: The reference star, used to locate the continuum. 
        Options are 'arcturus' or 'sun' (default).
    :type reference: str
    """
    # This path should point to the location of ardata.fits from
    # ftp://ftp.noao.edu/catalogs/arcturusatlas/visual/

    def __init__(self, velocity=None, reference='sun'):
        
        self.velocity = velocity
        self.reference = load_reference(reference)
        
        # Setup the logging if not existent yet
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                                format='%(message)s')
    

    def normalize_obs(self, observation, velocity, orders=None):
        """Normalize an Observation (multiorder) or a list of orders within
        
        :param observation: This observation will be normalized.
        :type observation: :class:`Observation`
        :param velocity: Velocity between the observation and reference
            spectrum.
        :type velocity: float
        :param orders: Which orders of the observation to normalize. Default 
            is None, which means all.
        :type orders: list, or None
        
        :return: The normalized observation.
        :rtype: :class:`NormalizedObservation`
        """
        
        norm_obs = NormalizedObservation(observation)
        if orders is None:
            orders = observation.orders
        for i in orders:
            norm_obs[i] = self.normalize_single(
                observation[i],
                velocity,
                return_type='flux'
            )
        return norm_obs


    def normalize_single(self, spectrum: Spectrum, velocity, return_type='spectrum'):
        """Normalize a single spectrum by comparison with the solar atlas, if
        possible. Fallback to envelope fitting.
        
        Returns either a new spectrum (default) or just the flux value, if
        keyword return_type='flux'.
        
        :param spectrum: The input spectrum to normalize.
        :type spectrum: :class:`Spectrum`
        :param velocity: Velocity between the observation and reference
            spectrum.
        :type velocity: float
        :param return_type: Whether to return as a :class:`Spectrum` 
            ('spectrum', default), or a `ndarray` with the flux values ('flux').
        :type return_type: str
        
        :return: The normalized spectrum.
        :rtype: :class:`Spectrum` or ndarray
        """
        # Divide by extracted blaze function
        if spectrum.cont is not None:
            flux = spectrum.flux / spectrum.cont
        else:
            flux = spectrum.flux
        # Try to use the solar spectrum
        try:
            flux /= self.compare_sun(spectrum.wave, flux, velocity)
        # Crude normalization - fit the upper envelope
        except ValueError:
            flux /= top(flux, 2, eps=0.00001)  # Better to do nothing?
        # Package result in a Spectrum object
        if return_type == 'flux':
            return flux
        if return_type == 'spectrum':
            return Spectrum(flux, spectrum.wave, spectrum.cont)


    def compare_sun(self, wave, flux, velocity, threshold=0.01) -> np.ndarray:
        """Compare a spectrum to the normalized reference spectrum and fit the 
        continuum
        
        :param wave: The wavelength vector of the spectrum.
        :type wave: ndarray
        :param flux: The flux vector of the spectrum.
        :type flux: ndarray
        :param velocity: Velocity between the observation and reference
            spectrum.
        :type velocity: float
        :param threshold: Only pixels this close to the continuum are used. 
            Defaults to 0.01.
        :type threshold: float
        
        :return: The fitted continuum flux.
        :rtype: ndarray
        """
        pix = np.arange(len(wave), dtype='float')

        # Resample the solar spectrum
        z = (1.0 + velocity / _c)
        work_sun = rebin(z * self.reference.wave, self.reference.flux, wave)

        # Detect the pixels where the sun is close to the continuum
        jj = np.where(np.abs(work_sun - 1.0) < threshold)
        if len(jj[0]) > 100 \
                and np.std(jj[0]) > 0.2 * len(wave) \
                and np.abs(np.std(flux[jj]) / np.mean(flux[jj])) < 3 * threshold:
            logging.debug('Fit the continuum...')
            
            yfit = np.polyval(np.polyfit(pix[jj], flux[jj], 4), pix)
            return yfit
        else:
            logging.debug('Dropping...')

            raise ValueError('Solar comparison failed - could not fit.')

    def guess_velocity(self, spec, normalize=True, delta_v=1000., maxlag=500):
        """Re-bin a spectrum to constant velocity spacing and cross-correlate
        with the solar spectrum to find the relative velocity offset.
        If given an observation or a list of spectra, vguess will return a
        median of the individual velocities.
        By default, the input spectrum is normalized by dividing with the
        supplied blaze or a fit to the upper envelope of the spectrum.
        This can be disabled using keyword 'normalize=False'.
        
        :param spec: The input spectrum.
        :type spec: :class:`Observation` or :class:`Spectrum`
        :param normalize: Whether the input spectrum should be normalized. 
            Default is True.
        :type normalize: bool
        :param delta_v: The velocity step size (in m/s).
        :type delta_v: float
        :param maxlag: The number of steps to each side in the cross-correlation.
        :type maxlag: int
        
        :return: The velocity offset.
        :rtype: float
        """
        return get_velocity_offset(spec, self.reference, normalize=normalize, 
                                   delta_v=delta_v, maxlag=maxlag)


# Shared tools:  # TODO: Move somewhere else..

def top(flux, deg=2, max_iter=40, eps=0.001) -> np.ndarray:
    """Fits a polynomial to the upper envelope of a spectrum - useful for
    normalization/blaze fitting.

    Inspired by IDL function TOP from REDUCE, Piskunov & Valenti (2002)
    
    :param flux: A flux vector.
    :type flux: ndarray
    :param deg: The polynomial degree. Defaults to 2.
    :type deg: int
    :param max_iter: The maximum number of iterations for the loop. Defaults 
        to 40.
    :type max_iter: int
    :param eps: If the step-to-step changes become smaller than this, stop the 
        loop. Defaults to 0.001.
    :type eps: float
    
    :return: The fitted continuum.
    :rtype: ndarray
    """
    n = len(flux)

    fmin = np.min(flux) - 1.0
    fmax = np.max(flux) + 1.0
    xx = np.linspace(-1.0, 1.0, n)
    ff = (flux - fmin) / (fmax - fmin)

    # Loop until converging or reaching max_iter
    i = 0
    if max_iter < 1:
        max_iter = 1  # Make sure at least one iteration
    while i < max_iter:
        ff_old = ff.copy()

        # Fit a polynomial to the data
        t = np.polyval(np.polyfit(xx, ff, deg), xx)

        # Fit another polynomial to the squared residuals
        yfit = np.polyval(np.polyfit(xx, (ff - t) ** 2, deg), xx)

        # Take the square root of the fitted polynomial wherever it is
        # above 0
        dev = np.zeros(n)
        jj = np.where(yfit > 0.0)
        dev[jj] = np.sqrt(yfit[jj])

        # Assign new values to ff
        jj = np.where((t - eps) > ff)
        ff[jj] = (t - eps)[jj]

        jj = np.where(ff < (t + dev * 3))
        ff[jj] = ff[jj]

        # Stop loop if converged
        if np.max(np.abs(ff - ff_old)) <= eps:
            break
        else:
            i += 1

    # Scale back to original axis and return
    return t * (fmax - fmin) + fmin


def get_velocity_offset(spectrum, reference, normalize=True, delta_v=1000., 
                        maxlag=500):
    """Find the velocity offset (m/s) between spectrum and reference
    
    If spectrum is an observation or a list of spectra, return value will
    be a median of the individual velocities.
    Keyword 'normalize' determines whether to normalize 'spectrum'.
    Reference is assumed to be normalized.

    Sign convention: spectrum = reference * (1 + v/c)
    
    :param spectrum: The input spectrum.
    :type spectrum: :class:`Spectrum` or :class:`MultiOrderSpectrum`
    :param reference: The reference spectrum.
    :type reference: :class:`Spectrum`
    :param normalize: Whether to normalize the input spectrum. Defaults to True.
    :type normalize: bool
    :param delta_v: The velocity step size (in m/s).
    :type delta_v: float
    :param maxlag: The number of steps to each side in the cross-correlation.
    :type maxlag: int
    
    :return: The velocity offset.
    :rtype: float
    """

    # Settings  # TODO: Add keywords for these
    #delta_v = 1000.0  # m/s #1000
    # How large does this have to be?
    #maxlag = 500  # Used to be 150 # Number of steps to each side in the cross-correlation
    
    # If first input argument is a list of spectra, do recursion
    if isinstance(spectrum, list) or isinstance(spectrum, MultiOrderSpectrum):
        results = []
        for s in spectrum:
            results.append(get_velocity_offset(s, reference, normalize=normalize,
                                               delta_v=delta_v, maxlag=maxlag))
        return np.median(results)

    # Unpack input spectrum, normalize unless keyword is set to False
    wave = spectrum.wave
    if normalize:
        if spectrum.cont is not None:
            ind0 = np.where(spectrum.cont == 0.)
            if len(ind0[0]) > 0:
                logging.info('0-values in continuum. Estimating continuum for normalization...')
                flux = spectrum.flux / top(spectrum.flux)
            else:
                try:
                    flux = spectrum.flux / spectrum.cont
                except Exception as e:
                    logging.error(e)
                    logging.error('Estimating continuum for normalization...')
                    flux = spectrum.flux / top(spectrum.flux)
        else:
            flux = spectrum.flux / top(spectrum.flux)
    else:
        flux = spectrum.flux

    # In order to use numpy.correlate() with mode 'valid', we prepare a narrow
    # and a wide vector with a difference in length corresponding to 2*maxlag

    # Number of points in the narrow vector
    n = int(np.log(wave[-1] / wave[0]) / np.log(1.0 + delta_v / _c))
    # The actual velocity step size
    delta_v = (np.power(10.0, np.log10(wave[-1] / wave[0]) / n) - 1.0) * _c
    # Create wavelength vectors with logarithmic steps
    logsteps = np.arange(-maxlag, n + maxlag, dtype='float')
    logwave_wide = wave[0] * np.power(1.0 + delta_v / _c, logsteps)
    logwave_narrow = logwave_wide[maxlag:-maxlag]

    # Rebin fluxes to new wavelength vectors
    flux = rebin(wave, flux, logwave_narrow)
    flux_ref = rebin(reference.wave, reference.flux, logwave_wide)

    # Create a velocity vector corresponding to the lag vector
    lags = np.arange(-maxlag, maxlag + 1, dtype='float')
    lagscale = lags * delta_v

    # Cross correlate the rebinned fluxes
    cc = np.correlate(flux, flux_ref, 'valid')

    # Pick the 5 points around the maximum of the cross-correlation
    ii = np.argmax(cc) + np.arange(-2, 3)
    # Fit a 2nd degree polynomial  # TODO: Use Gaussian instead?
    # (Offset subtracted in order to avoid high x-values)
    offset = np.mean(lagscale[ii])
    p = np.polyfit(lagscale[ii] - offset, cc[ii], 2)
    # Find the analytical maximum (and add back the offset)
    velocity = offset - 0.5 * p[1] / p[0]

    return velocity
