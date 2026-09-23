from astropy.io import fits

_ref_path = 'ardata.fits'

# Open the fits file
# (memmap=False in order to retain data access after closing handle)
with fits.open(_ref_path, memmap=False) as h:
    
    # Get wavelength
    ref_wave = h[1].data['wavelength']
    
    # Get arcturus flux
    arcturus_flux = h[1].data['arcturus']
    
    # Get solar flux
    #solar_flux = h[1].data['solarflux']
    
    # Get tellurics
    #tellurics = h[1].data['telluric']

print('Length of vector: ', len(ref_wave))
