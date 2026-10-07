""" A chunk is modelled against its own template chunk (phase 3, after prompt 6).

`pipe_lib.model_all_chunks` fits observation chunk *i* with ``chunk_ind=i``,
and `SimpleModel.eval` then takes ``stellar_template[chunk_ind]``.  That is
right only when the observation's chunks are the template's, one for one.
`chunks.auto_wave_comoving` builds chunks only for the orders it is given, so
on any subset of orders every chunk after the gap was modelled against a
template chunk from another order -- silently, because the template spline
extrapolates.  11 of the 31 cell-in HIRES epochs have a rejected iodine order.

The fix: `auto_wave_comoving` records each chunk's template index on the chunk
(``chunk.template_index``), and the model resolves the template from that.
The test builds the chunks of one epoch twice -- for all the template's orders
and for a subset -- and requires the same chunk to give the same model from
either list, each evaluated at its own list index as the driver does.
"""

# Standard imports
import os
import glob

import numpy as np
import pytest

from first_hires_exoplanet.utilities_hires import load_pyodine
from first_hires_exoplanet.utilities_hires import pyodine_parameters as hp


def _setup():
    import pyodine
    from pyodine.template.base import StellarTemplate_Chunked
    h5 = os.path.join(load_pyodine.DEFAULT_REDUX, 'template',
                      'hd187123_template_19980826_deconv_os1.h5')
    files = sorted(glob.glob(os.path.join(load_pyodine.DEFAULT_REDUX, 'reduce_19980825',
                                          'Science', 'spec1d_HI.19980825.19425*.fits')))
    if not os.path.exists(h5) or not files:
        pytest.skip('template or epoch not on disk')
    obs = load_pyodine.ObservationWrapper(files[0])
    template = StellarTemplate_Chunked(h5)
    return pyodine, obs, template


def test_subset_of_orders_uses_the_right_template_chunk():
    pyodine, obs, template = _setup()
    orders = template.orders_unique
    full = pyodine.chunks.auto_wave_comoving(obs, template, orders=orders, padding=6)
    sub = pyodine.chunks.auto_wave_comoving(obs, template, orders=orders[2:], padding=6)
    Pars = hp.Parameters()
    run = Pars.model_runs[0]
    model = pyodine.models.spectrum.SimpleModel(
        run['lsf_model'], run['wave_model'], run['cont_model'],
        load_pyodine.IodineTemplate(Pars.i2_to_use), stellar_template=template,
        osample_factor=Pars.osample_obs, conv_width=Pars.lsf_conv_width)
    # The 10th chunk of the subset is chunk 110 of the full list (two orders
    # of 50 chunks are missing in front of it)
    j_sub = 10
    j_full = j_sub + (len(full) - len(sub))
    a, b = full[j_full], sub[j_sub]
    assert a.order == b.order and np.array_equal(a.abspix, b.abspix)
    params = model.guess_params(a)
    m_full = model.eval(a, params, chunk_ind=j_full)
    m_sub = model.eval(b, params, chunk_ind=j_sub)
    assert np.allclose(m_full, m_sub), 'subset chunk modelled against the wrong template chunk'


def test_chunks_record_their_template_index():
    pyodine, obs, template = _setup()
    orders = template.orders_unique
    sub = pyodine.chunks.auto_wave_comoving(obs, template, orders=orders[2:], padding=6)
    for ch in sub:
        t = template[ch.template_index]
        assert t.order == ch.order
