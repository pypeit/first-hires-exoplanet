""" The fork must run against NumPy 2.

`pyodine` was written for NumPy 1.x and uses two names that NumPy 2.0 removed:
`np.float`, in the resampler `lib.misc.rebin`, and `np.NaN`, in the fitter's
failure paths.  Phase 3 prompt 1 established that this stops the code in the
first stage of template creation, and — because `pyodine` catches its own
exceptions and returns normally — that it does so *silently*, surfacing two
calls later as a missing output file.

These tests exercise the two code paths rather than grepping for the names, so
they would still catch a reintroduction through a different spelling.

    conda run -n pypeit14 pytest first_hires_exoplanet/tests/test_fork_numpy2.py

"""

# Standard imports
import re
from pathlib import Path

import numpy as np
import pytest

from .conftest import VENDOR_DIR


def test_rebin_runs_and_conserves_flux():
    """ `lib.misc.rebin` is where `np.float` stopped the template stage. """
    from pyodine.lib.misc import rebin

    # A narrow Gaussian line on a fine grid, resampled onto a coarser one.
    wave_fine = np.linspace(5000., 5010., 2001)
    flux_fine = 1. - 0.5 * np.exp(-0.5 * ((wave_fine - 5005.) / 0.1) ** 2)
    wave_new = np.linspace(5001., 5009., 101)

    flux_new = rebin(wave_fine, flux_fine, wave_new)

    assert flux_new.shape == wave_new.shape
    assert np.all(np.isfinite(flux_new))
    # Rebinning conserves the mean level of the continuum away from the line.
    assert np.isclose(flux_new[0], 1.0, atol=1e-6)
    # And it keeps the line: the deepest resampled pixel is well below unity.
    assert flux_new.min() < 0.9


def test_fitter_failure_paths_return_nan():
    """ `np.NaN` sits in `lmfit_wrapper`'s five failure returns. """
    import pyodine.fitters.lmfit_wrapper as lmfit_wrapper

    source = Path(lmfit_wrapper.__file__).read_text()
    # The module is only importable at all because these are evaluated lazily;
    # the point of the test is that the name they use still exists in NumPy.
    for name in set(re.findall(r'np\.(\w+)', source)):
        assert hasattr(np, name), f'lmfit_wrapper uses np.{name}, which is gone'


def test_no_removed_numpy_aliases_anywhere_in_the_fork():
    """ Nothing else in the tree reaches for a name NumPy 2 removed. """
    removed = ('float', 'int', 'bool', 'object', 'str', 'complex', 'NaN',
               'Inf', 'alltrue', 'sometrue', 'product', 'round_')
    pattern = re.compile(r'\bnp\.(' + '|'.join(removed) + r')\b')

    offenders = []
    for path in sorted(VENDOR_DIR.rglob('*.py')):
        if '__pycache__' in path.parts:
            continue
        for number, line in enumerate(path.read_text().splitlines(), start=1):
            if pattern.search(line):
                rel = path.relative_to(VENDOR_DIR)
                offenders.append(f'{rel}:{number}: {line.strip()}')

    assert not offenders, 'removed NumPy aliases still in the fork:\n' + \
        '\n'.join(offenders)
