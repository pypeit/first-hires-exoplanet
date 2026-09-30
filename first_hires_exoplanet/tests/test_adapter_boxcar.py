""" The adapter can read PypeIt's boxcar extraction (phase 3, after prompt 6).

The B stars over-fill the slit; with `force_center_obj` PypeIt extracts them
with a boxcar spanning each order, and its optimal extraction's mask collapses
(8-38% of iodine-order pixels kept in most frames, against 99% for the
boxcar).  The adapter read only ``OPT_*``.  ``extraction='BOX'`` reads
``BOX_*``; the default stays ``OPT``.
"""

# Standard imports
import os
import glob

import numpy as np
import pytest
from astropy.io import fits

from first_hires_exoplanet.utilities_hires import load_pyodine


def _bstar():
    files = sorted(glob.glob(os.path.join(load_pyodine.DEFAULT_REDUX, 'reduce_19980826_lsf',
                                          'Science', 'spec1d_HI.19980826.37213*.fits')))
    if not files:
        pytest.skip('B-star reduction not on disk')
    return files[0]


def test_boxcar_columns_are_read():
    path = _bstar()
    box = load_pyodine.ObservationWrapper(path, extraction='BOX')
    opt = load_pyodine.ObservationWrapper(path)
    with fits.open(path) as hdul:
        hdu = [h for h in hdul[1:] if hasattr(h, 'columns') and h.header.get('ECH_ORDER') == 63][0]
        i = box.index_of(63)
        good = hdu.data['BOX_COUNTS'] != 0
        assert np.allclose(box._flux[i][good], hdu.data['BOX_COUNTS'][good])
        good = hdu.data['OPT_COUNTS'] != 0
        assert np.allclose(opt._flux[i][good], hdu.data['OPT_COUNTS'][good])
    assert np.mean(box._weight[1:15] > 0) > 0.9


def test_unknown_extraction_refused():
    with pytest.raises(ValueError):
        load_pyodine.ObservationWrapper(_bstar(), extraction='FOO')
