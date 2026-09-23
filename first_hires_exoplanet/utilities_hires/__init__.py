""" A `pyodine` instrument adapter for Keck/HIRES reduced with PypeIt.

Laid out like `pyodine`'s own `utilities_lick/`, so it can be dropped into a
vendored `pyodine` tree unchanged.  Phase 2 prompt 8 writes only the input
adapter; `conf.py`, `pyodine_parameters.py` and `timeseries_parameters.py`
belong to phase 3.
"""

from .load_pyodine import (ObservationWrapper, IodineTemplate, load_file,
                           get_star, get_instrument, HIRES)

__all__ = ['ObservationWrapper', 'IodineTemplate', 'load_file', 'get_star',
           'get_instrument', 'HIRES']
