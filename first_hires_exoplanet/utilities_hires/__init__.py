""" A `pyodine` instrument adapter for Keck/HIRES reduced with PypeIt.

Laid out like `pyodine`'s own `utilities_lick/`, so it can be dropped into a
vendored `pyodine` tree unchanged.  Phase 2 prompt 8 wrote the input
adapter; phase 3 prompt 4 added `conf.py`, `pyodine_parameters.py`,
`timeseries_parameters.py` and `logging.json`, so the package now exposes
what `pyodine`'s drivers read from `utilities_lick`: ``load_pyodine``,
``conf``, ``pyodine_parameters`` and ``timeseries_parameters``.

`logging.json` names no file.  `pyodine.lib.misc.setup_logging` writes the
info and error logs only where a driver passes their paths, and deletes both
file handlers otherwise; Lick's copy hard-codes ``/home/paul/`` paths that
would fail the moment a driver asked for a log.
"""

from .load_pyodine import (ObservationWrapper, IodineTemplate, load_file,
                           get_star, get_instrument, HIRES)
from . import load_pyodine, conf, pyodine_parameters, timeseries_parameters

__all__ = ['ObservationWrapper', 'IodineTemplate', 'load_file', 'get_star',
           'get_instrument', 'HIRES', 'load_pyodine', 'conf',
           'pyodine_parameters', 'timeseries_parameters']
