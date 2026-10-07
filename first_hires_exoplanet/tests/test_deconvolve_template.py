""" Tests of the machinery `deconvolve_template.py` measures with (phase 3,
prompt 3).

Three properties the prompt-3 numbers depend on and that would fail silently:

* masked pixels arrive with a wavelength of exactly zero, and `pyodine`
  divides by the first and last wavelength of every order;
* the velocity estimator's sign, which every "dv" in the report inherits;
* the reference Jansson with an explicit adjoint agrees with the fork's to
  machine precision for a symmetric kernel, so any difference measured with an
  asymmetric one is the adjoint and not the harness.
"""

# Standard imports
import numpy as np
import pytest

dt = pytest.importorskip('first_hires_exoplanet.deconvolve_template')


def test_fill_wavelengths_replaces_zeros_monotonically():
    """ Zeros at both ends and in a run are replaced; good pixels untouched. """
    pix = np.arange(2048, dtype=float)
    true = 5000. + 0.03 * pix + 1e-6 * pix ** 2
    wave = true.copy()
    wave[[0, -1]] = 0.
    wave[1500:1700] = 0.
    out = dt.fill_wavelengths(wave)
    good = wave > 0
    assert np.array_equal(out[good], wave[good])
    assert np.all(out > 0)
    assert np.all(np.diff(out) > 0)
    assert np.max(np.abs(out - true)) < 1e-6


def test_velocity_shift_sign():
    """ A spectrum redshifted by +v comes back as +v. """
    dt._check_velocity_sign()


def test_true_adjoint_matches_fork_for_symmetric_kernel():
    """ Identical for a symmetric kernel; different for a skewed one. """
    from pyodine.template.deconvolve import jansson
    rng = np.random.default_rng(0)
    x = np.linspace(-20, 20, 401)
    truth = 1. - 0.6 * np.exp(-0.5 * (x / 0.4) ** 2) - 0.3 * np.exp(-0.5 * ((x - 5) / 0.3) ** 2)
    sym = dt.gaussian_lsf(2.2)
    skew = dt.skewed_lsf(2.2)
    kw = dict(a=0., b=1.02, delta=0.2, chi_change=1e-6)
    for kernel, same in ((sym, True), (skew, False)):
        obs = np.convolve(truth, kernel, 'same') + rng.normal(size=len(x)) * 1e-3
        a = jansson(obs.copy(), kernel, 300, **kw)
        b = dt.jansson_true_adjoint(obs.copy(), kernel, 300, **kw)
        if same:
            assert np.max(np.abs(a - b)) < 1e-12
        else:
            assert np.max(np.abs(a - b)) > 1e-6
