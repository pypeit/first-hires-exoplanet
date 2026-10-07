""" The adapter's objects survive `pyodine`'s results writer (phase 3, prompt 5).

`results_io.create_results_dict` reads ``obs.star.proper_motion``, which the
adapter's ``Star`` did not have.  The first full-epoch run fitted all 700
chunks and then died saving them -- and reported it through `logging.info`,
so its error log stayed empty.  This fits one real chunk and saves and loads
it, which is the whole path that broke.
"""

# Standard imports
import os

import numpy as np
import pytest

from first_hires_exoplanet.utilities_hires import load_pyodine


def test_star_has_proper_motion():
    star = load_pyodine.Star('HD 187123')
    assert star.proper_motion == (None, None)
    star = load_pyodine.Star('HD 187123', pmra=1.5, pmdec=-2.5)
    assert star.proper_motion == (1.5, -2.5)


def test_one_chunk_saves_and_loads(tmp_path):
    from pyodine.fitters import results_io
    from first_hires_exoplanet import fit_one_chunk as f1
    from first_hires_exoplanet import deconvolve_template as dt
    from first_hires_exoplanet.utilities_hires import pyodine_parameters as hp
    if not os.path.exists(dt.TEMPLATE_H5[10]):
        pytest.skip('prompt-3 template not on disk')
    Pars = hp.Parameters()
    prep = f1.prepare(Pars, f1.epoch_file(), dt.TEMPLATE_H5[10])
    model, fitter = f1.build_run0(Pars, prep)
    params = f1.run0_params(Pars, prep, model, fitter)
    i = f1.pick_chunk(prep)[0]
    result = f1.fit_chunk(fitter, prep['chunks'][i], params[i], prep['weight'], i)[0]
    out = str(tmp_path / 'one.h5')
    results_io.save_results(out, [result], filetype='h5py')
    back = results_io.load_results(out)
    assert np.isclose(back['params']['velocity'][0], result.params['velocity'])
