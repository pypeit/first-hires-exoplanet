""" Epoch orders are matched to the template by echelle number (phase 3, prompt 7).

`pyodine` maps template order ``o`` to observation position
``o + order_correction``, one constant shift.  1998-09-13 was extracted with
34 orders and no order 71, so that shift sends template order 71 onto echelle
order 72: 50 chunks fitted against the wrong wavelengths, velocities of
millions of m/s, and no error.  `fit_all_epochs.usable_orders` drops any
template order whose position does not hold the same echelle order.  The
mapping of result chunks back to template chunks is tested with it.
"""

import numpy as np
import pytest

from first_hires_exoplanet import fit_all_epochs as fa
from first_hires_exoplanet import fit_one_chunk as f1
from first_hires_exoplanet import deconvolve_template as dt
from first_hires_exoplanet.utilities_hires import load_pyodine as lp


@pytest.fixture(scope='module')
def template():
    from pyodine.template.base import StellarTemplate_Chunked
    try:
        return StellarTemplate_Chunked(dt.TEMPLATE_H5[fa.OSAMPLE])
    except (OSError, FileNotFoundError):
        pytest.skip('deconvolved template not on disk')


def _obs(koaid):
    try:
        return lp.ObservationWrapper(f1.epoch_file(koaid))
    except FileNotFoundError:
        pytest.skip('{:s} not reduced on disk'.format(koaid))


def test_order_gap_is_not_bridged(template):
    obs = _obs('HI.19980913.29698')
    assert 71 not in list(obs.ech_orders)
    keep, rejected, corr = fa.usable_orders('HI.19980913.29698', obs, template)
    tech = fa.template_echelle_orders()
    assert 71 in rejected
    # Every kept template order sits on its own echelle order
    for to in keep:
        assert int(obs.ech_orders[int(to) + corr]) == tech[int(to)]
    # The constant shift alone would have bridged the gap
    assert int(obs.ech_orders[14 + corr]) == 72 and tech[14] == 71


def test_full_format_keeps_every_good_order(template):
    obs = _obs('HI.19980825.19425')
    keep, rejected, corr = fa.usable_orders('HI.19980825.19425', obs, template)
    assert corr == 0 and rejected == [] and len(keep) == 14


def test_template_columns_follow_the_order(template):
    from astropy.table import Table
    # Two orders, 3 and 7, chunks listed out of pixel order
    rows = []
    for pos in (7, 3):
        for k in (2, 0, 1):
            rows.append(dict(pos=pos, pix0=100 * k))
    t = Table(rows)
    with pytest.raises(RuntimeError):
        fa.template_columns(t, template, 0)            # 3 chunks where the template has 50
    n = len(template.get_order_indices(3))
    t = Table(dict(pos=np.r_[np.full(n, 7), np.full(n, 3)],
                   pix0=np.r_[np.arange(n)[::-1], np.arange(n)] * 40))
    col = fa.template_columns(t, template, 0)
    assert list(col[n:]) == list(template.get_order_indices(3))
    assert list(col[:n]) == list(template.get_order_indices(7))[::-1]
