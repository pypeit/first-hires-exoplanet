from os.path import dirname
from astropy.io import fits

from ..components import Spectrum

_ref_path = dirname(__file__) + '/arcturus/ardata.fits'


# (memmap=False in order to retain data access after closing handle)
def load_reference(name='arcturus'):
    """Load the solar/arcturus reference spectrum. Required argument 'name'
    can be either 'arcturus' (default) or 'sun'.
    
    :param name: The name of the reference spectrum: 'arcturus' (default), or
        'sun'.
    :type name: str
    
    :return: The loaded reference spectrum.
    :rtype: :class:`Spectrum`
    """
    with fits.open(_ref_path, memmap=False) as h:
        ref_wave = h[1].data['wavelength']
        if name.lower() == 'sun':
            ref_flux = h[1].data['solarflux']
            return Spectrum(ref_flux, ref_wave)
        elif name.lower() == 'arcturus':
            ref_flux = h[1].data['arcturus']
            return Spectrum(ref_flux, ref_wave)
        else:
            raise ValueError('Unknown reference "{}"'.format(name))
