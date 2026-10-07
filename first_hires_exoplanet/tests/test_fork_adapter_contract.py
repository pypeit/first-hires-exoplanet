""" The HIRES adapter must obey the fork's contract, now that it subclasses it.

Phase 2 wrote `utilities_hires` against stand-in base classes, because
`pyodine` was not yet vendored.  Phase 3 prompt 2 put the real fork on the
path, and one collision surfaced immediately: `components.MultiOrderSpectrum`
defines `orders` as a read-only property giving *positions* — which is what
`pyodine`'s own loops index with — while the adapter had been assigning the
echelle order *numbers* (57-93) to the same name.  Assignment raised, so this
one at least was loud; the tests below pin the resolution so it stays fixed.

They need the reduced data, and skip without it.

    conda run -n pypeit14 pytest first_hires_exoplanet/tests/test_fork_adapter_contract.py

"""

# Standard imports
import glob
import os

import numpy as np
import pytest

from first_hires_exoplanet.utilities_hires import load_pyodine


@pytest.fixture(scope='module')
def observation():
    """ One reduced HIRES epoch, as `pyodine` will see it.

    Returns
    -------
    :obj:`ObservationWrapper`
    """
    files = sorted(glob.glob(os.path.join(
        load_pyodine.DEFAULT_REDUX, 'reduce_19*', 'Science', 'spec1d_*.fits')))
    if not files:
        pytest.skip('reduced data not present on disk')
    return load_pyodine.ObservationWrapper(files[0])


def test_the_adapter_subclasses_the_real_fork():
    """ Not the phase-2 stand-ins. """
    assert load_pyodine.HAVE_PYODINE, 'vendored pyodine is not importable'


def test_orders_means_positions_as_pyodine_expects(observation):
    """ `for i in obs.orders: obs[i]` has to work, and it indexes positions. """
    assert np.array_equal(observation.orders, np.arange(observation.nord))
    for i in observation.orders:
        assert len(observation[i]) == observation.npix


def test_echelle_numbers_live_under_their_own_name(observation):
    """ And they are the physical HIRES orders, not positions. """
    assert observation.ech_orders.min() >= 30
    assert observation.ech_orders.max() <= 120
    assert len(observation.ech_orders) == observation.nord
    # `index_of` bridges the two.
    first = int(observation.ech_orders[0])
    assert observation.index_of(first) == 0


def test_each_order_carries_its_inverse_variance(observation):
    """ So a chunk sliced out of an order is weightable on its own. """
    order = observation[0]
    assert order.ivar is not None
    assert np.allclose(order.ivar, observation._weight[0])
    assert np.allclose(order.compute_weight(weight_type='ivar'), order.ivar)


def test_rejected_orders_are_weightless_but_still_loadable(observation):
    """ Phase 2's rule: zero-weighted, never zero-fluxed.

    `Spectrum.__init__` raises `NoDataError` on an all-zero flux vector, so an
    order excluded by zeroing its flux would become unloadable rather than
    merely unused.
    """
    weight = observation.compute_weight()
    rejected = [i for i in observation.orders if not np.any(weight[i] > 0)]
    if not rejected:
        pytest.skip('this epoch has no rejected orders')
    for i in rejected:
        assert np.any(observation[i].flux), 'flux was zeroed, not the weight'
        assert np.all(observation[i].compute_weight(weight_type='ivar') == 0.)


def test_the_documented_barycentric_units_hold(observation):
    """ The fork's `check_bary_units` must pass on what the adapter builds. """
    observation.check_bary_units()
    assert observation.bary_date > 2.4e6
    assert abs(observation.bary_vel_corr) < 4.0e4
