""" The adapter gives every pixel a wavelength (phase 3, prompt 4).

PypeIt writes masked pixels as exact zeros in the *wavelength* column, not
only the flux: the first and last pixel of every order of every spec1d, and
long runs where an order runs off the chip.  `pyodine`'s velocity guess
(`template.normalize.get_velocity_offset`) divides by ``wave[0]`` and
``wave[-1]`` of each order, so an epoch loaded as phase 2 left it fails at the
first step of `model_single_observation`.  Phase 3 prompt 3 found this in the
template; the Q&A after it decided the adapter should be fixed here.

Also pinned: ``IodineTemplate`` accepts the integer atlas index `pyodine`'s
drivers pass (``IodineTemplate(Pars.i2_to_use)``), resolved through
``conf.my_iodine_atlases``, exactly as upstream's Lick adapter does.
"""

# Standard imports
import os
import glob

import numpy as np
import pytest
from astropy.io import fits

from first_hires_exoplanet.utilities_hires import load_pyodine


def _spec1d():
    files = sorted(glob.glob(os.path.join(
        load_pyodine.DEFAULT_REDUX, 'reduce_19980826', 'Science', 'spec1d_*.fits')))
    if not files:
        pytest.skip('no 1998-08-26 spec1d on disk')
    return files[0]


def test_every_wavelength_is_positive_and_monotonic():
    obs = load_pyodine.ObservationWrapper(_spec1d())
    assert np.all(obs._wave > 0), 'zero wavelengths reach pyodine'
    assert np.all(np.diff(obs._wave, axis=1) > 0)


def test_filled_pixels_carry_no_weight():
    """ The wavelengths PypeIt zeroed are filled, never fitted. """
    path = _spec1d()
    obs = load_pyodine.ObservationWrapper(path)
    with fits.open(path) as hdul:
        for hdu in hdul[1:]:
            if not hasattr(hdu, 'columns') or 'OPT_WAVE' not in hdu.columns.names:
                continue
            zero = np.asarray(hdu.data['OPT_WAVE']) <= 0
            i = obs.index_of(int(hdu.header['ECH_ORDER']))
            assert np.all(obs._weight[i][zero] == 0.)


def test_velocity_guess_runs_on_an_epoch():
    """ The call that divides by wave[0] and wave[-1] returns a number. """
    from pyodine.template.normalize import SimpleNormalizer
    obs = load_pyodine.ObservationWrapper(_spec1d())
    blue = [i for i in range(obs.nord) if np.median(obs._wave[i]) < 5000.][:4]
    v = SimpleNormalizer(reference='arcturus').guess_velocity(obs[blue])
    assert np.isfinite(v)


def test_no_zero_flux_and_zeroed_pixels_carry_no_weight():
    """ PypeIt zero-fluxes masked pixels too.  1998-08-26 19326 has a run of
    113 in order 70, and a 40-pixel chunk inside it makes `Spectrum` raise
    `NoDataError` -- for the whole epoch, because every chunk is built before
    any is fitted.  Fill the flux, keep the weight at zero. """
    path = _spec1d()
    obs = load_pyodine.ObservationWrapper(path)
    assert np.all(obs._flux != 0.)
    with fits.open(path) as hdul:
        for hdu in hdul[1:]:
            if not hasattr(hdu, 'columns') or 'OPT_COUNTS' not in hdu.columns.names:
                continue
            zero = np.asarray(hdu.data['OPT_COUNTS']) == 0
            i = obs.index_of(int(hdu.header['ECH_ORDER']))
            assert np.all(obs._weight[i][zero] == 0.)


def test_every_chunk_of_an_epoch_can_be_built():
    """ `auto_wave_comoving` against the prompt-3 template, end to end. """
    import pyodine
    from pyodine.template.base import StellarTemplate_Chunked
    h5 = os.path.join(load_pyodine.DEFAULT_REDUX, 'template',
                      'hd187123_template_19980826_deconv_os10.h5')
    if not os.path.exists(h5):
        pytest.skip('prompt-3 template not on disk')
    files = sorted(glob.glob(os.path.join(load_pyodine.DEFAULT_REDUX, 'reduce_19980826',
                                          'Science', 'spec1d_HI.19980826.19326*.fits')))
    if not files:
        pytest.skip('1998-08-26 19326 not on disk')
    obs = load_pyodine.ObservationWrapper(files[0])
    template = StellarTemplate_Chunked(h5)
    chunks = pyodine.chunks.auto_wave_comoving(obs, template,
                                               orders=template.orders_unique, padding=6)
    assert len(chunks) == len(template)


def test_iodine_template_accepts_atlas_index():
    from first_hires_exoplanet.utilities_hires import conf
    iod = load_pyodine.IodineTemplate(1)
    assert iod.orig_filename == conf.my_iodine_atlases[1]
    assert iod.wave_frame == 'vacuum'
