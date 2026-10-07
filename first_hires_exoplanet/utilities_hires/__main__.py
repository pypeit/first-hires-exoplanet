""" Run the adapter's self-check.

    conda run -n pypeit14 python -m first_hires_exoplanet.utilities_hires

A separate entry point so the check can be run without importing the module
twice, which `python -m ...load_pyodine` does once the package `__init__`
has already imported it.
"""

from .load_pyodine import _self_check

if __name__ == '__main__':
    raise SystemExit(0 if _self_check() else 1)
