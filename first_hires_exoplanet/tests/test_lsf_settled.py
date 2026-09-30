""" The instrumental profile prompt 6 settled on, and the start rule it needed.

`compare_lsf_models.py` found the four workable `pyodine` LSF models equal in
per-chunk scatter and chose the super-Gaussian for run 1; it also found that
run 1 cannot start a shape parameter at the ~1e-12 `fit_lsfs` hands it, or
`lmfit`'s relative finite-difference step leaves the whole fit at its start.
"""

# Standard imports
import pytest

from first_hires_exoplanet.utilities_hires import pyodine_parameters as hp


def test_run1_is_the_super_gaussian():
    for Pars in (hp.Parameters(), hp.Template_Parameters()):
        assert Pars.model_runs[0]['lsf_model'].name() == 'SingleGaussian'
        assert Pars.model_runs[1]['lsf_model'].name() == 'SuperGaussian'


def test_template_is_unoversampled():
    assert hp.Template_Parameters().deconvolution_pars['osample_temp'] == 1


@pytest.mark.parametrize('value, lo, hi, expect', [
    (6.2e-12, 0.0, 1.0, 0.02),      # a satellite fit_lsfs sent to zero: margin
    (-4e-13, -0.5, 0.5, -1e-3),     # a Hermite weight: the floor, sign kept
    (2.0, 1.0, 4.0, 2.0),           # an exponent well inside: untouched
    (4.0, 1.0, 4.0, 3.94),          # on a bound: pulled inside
])
def test_start_inside(value, lo, hi, expect):
    assert hp.start_inside(value, lo, hi) == pytest.approx(expect)


def test_super_bounds_are_physical():
    lo, hi = hp.lsf_bounds('SuperGaussian', 'exponent', -5.)
    assert lo >= 1.0 and hi <= 4.0
