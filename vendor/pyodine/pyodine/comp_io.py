from astropy.io import fits as pyfits
from os.path import isfile
from os import remove


def save_fits(filename, data, add_header=None):
    """Save image data (and header) to a fits file
    
    :param filename: The filename to save the observation to. If file exists,
        it is overwritten!
    :type filename: str
    :param data: The data array to save.
    :type data: ndarray or list
    :param header: An instance of a fits header or a dictionary. If None, only
        a minimal default header is saved.
    :type header: :class:`fits.header`, dict, or None
    """
    if isfile(filename):
        remove(filename)
    hdu = pyfits.PrimaryHDU()
    hdu.data = data
    
    if isinstance(add_header, dict):
        for key, value in add_header.items():
            hdu.header[key] = value
    elif isinstance(add_header, pyfits.header.Header):
        hdu.header = add_header
    elif add_header is not None:
        raise ValueError('"add_header" must be of type dict or tuple!')
    
    hdu.writeto(filename, overwrite=True)