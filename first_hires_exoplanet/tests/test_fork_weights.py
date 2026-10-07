""" `compute_weight` must accept a propagated inverse variance.

`pyodine` offered two weightings: 'flat' (all ones) and 'inverse'
(1/(f(1+f*rel_noise^2)), an estimate of the variance *from the flux itself*,
inherited from the dop code).  Neither knows that a modern pipeline has
already propagated a real variance through bias, flat, sky and extraction.
PypeIt has: `OPT_COUNTS_IVAR` is the inverse variance of every extracted pixel,
and phase 2's adapter carries it.

This is not only about weighting quality.  `fitters/lmfit_wrapper.py:96`
builds its residual as

    (model - flux) * sqrt(abs(weight))

so the weight *is* the inverse variance as far as the fit is concerned: with
the right array there, the residual is a proper chi and the reduced chi-square
means what it says.  With 'flat' weights it does not, which is why phase 3
prompt 1 measured a median reduced chi-square of 150516 on a fit that had
visibly worked.

Zero weight is also how this project excludes data.  11% of HIRES
order-spectra are unusable, and `Spectrum.__init__` raises on an all-zero flux
vector, so a bad order must be zero-*weighted*, never zero-*fluxed*.

    conda run -n pypeit14 pytest first_hires_exoplanet/tests/test_fork_weights.py

"""

# Standard imports
import numpy as np
import pytest


@pytest.fixture
def spectrum_with_ivar():
    """ A spectrum carrying a propagated inverse variance.

    Returns
    -------
    :obj:`pyodine.components.Spectrum`
    """
    from pyodine.components import Spectrum

    flux = np.linspace(100., 1000., 10)
    return Spectrum(flux, wave=np.linspace(5000., 5001., 10), ivar=1. / flux)


def test_ivar_is_stored_and_returned(spectrum_with_ivar):
    """ The weighting the extraction already computed, used as it stands. """
    weight = spectrum_with_ivar.compute_weight(weight_type='ivar')
    assert np.allclose(weight, 1. / spectrum_with_ivar.flux)


def test_ivar_survives_slicing(spectrum_with_ivar):
    """ Chunks are slices of orders; the variance has to come along.

    `Chunk.__init__` builds itself from `observation[order][pixels]`, so an
    inverse variance that is dropped by `__getitem__` would silently become
    unavailable exactly where the fit needs it.
    """
    piece = spectrum_with_ivar[2:5]
    assert piece.ivar is not None
    assert np.allclose(piece.ivar, spectrum_with_ivar.ivar[2:5])
    assert np.allclose(piece.compute_weight(weight_type='ivar'),
                       spectrum_with_ivar.ivar[2:5])


def test_bad_pixels_get_zero_weight_not_a_nan():
    """ Masked pixels arrive as exact zeros, and must end up weightless.

    A masked pixel has zero inverse variance; an unphysical negative one, or a
    NaN from a division somewhere upstream, must not be allowed to propagate
    into `sqrt(abs(weight))` and quietly become a large weight.
    """
    from pyodine.components import Spectrum

    flux = np.array([100., 100., 100., 100., 100.])
    ivar = np.array([0.01, 0.0, -0.5, np.nan, np.inf])
    weight = Spectrum(flux, ivar=ivar).compute_weight(weight_type='ivar')

    assert np.all(np.isfinite(weight))
    assert weight[0] == pytest.approx(0.01)
    assert np.all(weight[1:] == 0.0)


def test_asking_for_ivar_without_one_is_an_error():
    """ Not silently answered with ones, which would look like it worked. """
    from pyodine.components import NoDataError, Spectrum

    spectrum = Spectrum(np.ones(5))
    with pytest.raises(NoDataError):
        spectrum.compute_weight(weight_type='ivar')


def test_the_old_weightings_are_untouched(spectrum_with_ivar):
    """ Existing instruments keep exactly the behaviour they had. """
    flat = spectrum_with_ivar.compute_weight(weight_type='flat')
    assert np.allclose(flat, np.ones(len(spectrum_with_ivar)))

    rel_noise = 0.008
    inverse = spectrum_with_ivar.compute_weight(weight_type='inverse',
                                                rel_noise=rel_noise)
    flux = spectrum_with_ivar.flux
    assert np.allclose(inverse, 1. / (flux * (1. + flux * rel_noise ** 2)))

    with pytest.raises(NotImplementedError):
        spectrum_with_ivar.compute_weight(weight_type='nonsense')


def test_an_order_loop_carries_the_ivar_through():
    """ `MultiOrderSpectrum.compute_weight` must reach the same array.

    This is the call the drivers actually make: one (nord, npix) weight array
    per observation, sliced per chunk afterwards.
    """
    from pyodine.components import MultiOrderSpectrum, Spectrum

    class _TwoOrders(MultiOrderSpectrum):
        nord, npix = 2, 4

        def __getitem__(self, order):
            flux = np.full(self.npix, 100.)
            ivar = np.full(self.npix, 0.01) * (order + 1)
            # The second order is rejected, as 11% of HIRES orders are.
            if order == 1:
                ivar = np.zeros(self.npix)
            return Spectrum(flux, ivar=ivar)

    weight = _TwoOrders().compute_weight(weight_type='ivar')
    assert weight.shape == (2, 4)
    assert np.allclose(weight[0], 0.01)
    assert np.all(weight[1] == 0.0)


def test_the_hires_adapter_speaks_the_forks_name():
    """ 'ivar' must mean the same thing on both sides of the boundary. """
    from first_hires_exoplanet.utilities_hires.load_pyodine import (
        ObservationWrapper)

    assert 'ivar' in ObservationWrapper.compute_weight.__doc__
