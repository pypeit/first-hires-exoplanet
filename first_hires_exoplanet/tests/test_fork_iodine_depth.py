""" The iodine depth parameter must act as an optical depth, not a multiplier.

`pyodine` scaled the iodine atlas with

    flux = iod_depth * (flux - 1) + 1

which is a linear scaling of line *depth*.  The physics of a cell with a
different column density is Beer-Lambert: transmission goes as ``T**alpha``.
The two agree to first order for weak lines and diverge badly for strong ones,
and phase 2 measured the HIRES cell at ``alpha = 2.59`` against the Fischer
atlas the fork ships.  At that exponent the linear form is not merely
inaccurate, it is unphysical: the atlas reaches ``flux_normalized = 0.000``,
and a saturated line scaled linearly by 2.59 returns ``-1.59`` — negative
transmission, over 5.8% of the atlas.

    conda run -n pypeit14 pytest first_hires_exoplanet/tests/test_fork_iodine_depth.py

"""

# Standard imports
import numpy as np
import pytest


#: The exponent phase 2 measured between the Fischer atlas and the HIRES cell
ALPHA_HIRES = 2.59


def test_unit_depth_is_the_identity():
    """ `iod_depth = 1` must leave the atlas exactly alone. """
    from pyodine.models.spectrum import scale_iodine_depth

    flux = np.linspace(0., 1., 101)
    assert np.allclose(scale_iodine_depth(flux, 1.0), flux)


def test_saturated_lines_stay_physical():
    """ The failure the linear form produces, at the exponent HIRES needs. """
    from pyodine.models.spectrum import scale_iodine_depth

    flux = np.array([0.0, 0.01, 0.1, 0.5, 0.9, 1.0])
    scaled = scale_iodine_depth(flux, ALPHA_HIRES)

    assert np.all(scaled >= 0.0), 'transmission went negative'
    assert np.all(scaled <= 1.0), 'transmission exceeded unity'
    # A saturated line is still saturated, not driven to -1.59 as the linear
    # form would have it.
    assert scaled[0] == 0.0
    # Deeper lines get deeper, shallow ones stay shallow: a column density,
    # not a multiplier.
    assert scaled[3] < flux[3]
    assert scaled[-1] == 1.0


def test_depths_compose_multiplicatively():
    """ Two column densities in series are one column density.

    This is what lets the HIRES adapter pre-scale the atlas by `alpha` and
    still let the fit vary `iod_depth` around 1: the two compose exactly,
    which the linear form does not do.
    """
    from pyodine.models.spectrum import scale_iodine_depth

    flux = np.linspace(0.001, 1., 51)
    twice = scale_iodine_depth(scale_iodine_depth(flux, ALPHA_HIRES), 1.1)
    once = scale_iodine_depth(flux, ALPHA_HIRES * 1.1)
    assert np.allclose(twice, once)


def test_agrees_with_the_linear_form_for_weak_lines():
    """ Backwards compatibility where the old form was defensible.

    For shallow lines the two differ at second order in the line depth, so
    nothing that was fitting sensibly before should move appreciably.
    """
    from pyodine.models.spectrum import scale_iodine_depth

    flux = 1. - np.array([1e-4, 1e-3, 5e-3])
    depth = 1.05
    linear = depth * (flux - 1.) + 1.
    assert np.allclose(scale_iodine_depth(flux, depth), linear, rtol=1e-4)


def test_the_model_never_emits_negative_iodine_transmission():
    """ End to end, through `SimpleModel.eval`, on a saturated atlas.

    The helper above could be right while the model still applied the old
    formula, so this drives the actual evaluation path with a synthetic
    observation whose atlas contains a black line core.
    """
    from pyodine.components import IodineAtlas, Observation, Chunk, Spectrum
    from pyodine.models.spectrum import SimpleModel
    from pyodine.models.lsf import SingleGaussian
    from pyodine.models.wave import LinearWaveModel
    from pyodine.models.cont import LinearContinuumModel
    from pyodine.models.base import ParameterSet

    # The chunk is a slice out of the middle of a longer order, because
    # `chunk.padded` reaches beyond the chunk and must stay inside the data.
    n_pix, first, last = 200, 70, 130
    wave = np.linspace(4999., 5002., n_pix)

    class _FakeObservation(Observation):
        """ The smallest thing `Chunk` will accept. """
        nord = 1
        npix = n_pix

        def __getitem__(self, order):
            return Spectrum(np.ones(n_pix), wave=wave, cont=np.ones(n_pix))

    # An atlas with a saturated core, on a grid wide enough to cover the
    # padded chunk after the wavelength model has stretched it.
    atlas_wave = np.linspace(4995., 5005., 20001)
    atlas_flux = 1. - np.exp(-0.5 * ((atlas_wave - 5000.5) / 0.02) ** 2)
    atlas = IodineAtlas(atlas_flux, wave=atlas_wave)

    # The sub-models are used as classes, not instances: they subclass
    # `StaticModel`, which refuses to be instantiated.
    model = SimpleModel(SingleGaussian, LinearWaveModel,
                        LinearContinuumModel, atlas, osample_factor=4,
                        conv_width=6.)
    chunk = Chunk(_FakeObservation(), 0, np.arange(first, last), padding=5)

    params = ParameterSet(velocity=0., tem_depth=1., iod_depth=ALPHA_HIRES)
    params.add(model.lsf_model.guess_params(chunk), prefix='lsf')
    params.add(ParameterSet(intercept=wave[(first + last) // 2],
                            slope=wave[1] - wave[0]), prefix='wave')
    params.add(ParameterSet(intercept=1., slope=0.), prefix='cont')

    spectrum = model.eval(chunk, params)

    assert np.all(np.isfinite(spectrum))
    assert spectrum.min() >= 0.0, \
        f'model emitted negative flux ({spectrum.min():.3f})'
