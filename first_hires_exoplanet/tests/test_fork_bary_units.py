""" `bary_date` and `bary_vel_corr` must be documented as they are used.

`components.Observation` declared

    bary_date = None       # Mid-time as Barycentric Reduced Julian Date (BJD - 2400000.0)
    bary_vel_corr = None   # Barycentric velocity correction (km/s)

and both lines are wrong about the code underneath them.
`timeseries/bary_vel_corr.py` hands `bary_date` straight to
`barycorrpy.get_BC_vel(JDUTC=...)`, which wants a **full Julian Date in UTC** —
not a barycentric date, and not a reduced one.  `chunks.py` divides
`bary_vel_corr` by `astropy.constants.c`, whose value is in **m/s**.

Neither mistake announces itself.  Feeding a reduced JD produces a
barycentric correction for the wrong century, and feeding km/s produces a
velocity shift a thousand times too small; both come back as plausible
numbers.  So beyond correcting the comments this adds the one check that can
actually be made mechanically: a reduced Julian Date is three orders of
magnitude from a full one and can be recognised on sight.  A km/s value
cannot be distinguished from a small m/s one, which is exactly why the
documentation has to be right.

    conda run -n pypeit14 pytest first_hires_exoplanet/tests/test_fork_bary_units.py

"""

# Standard imports
import re
from pathlib import Path

import pytest


def _declaration(name):
    """ The source line declaring a class attribute of `Observation`.

    Parameters
    ----------
    name : :obj:`str`
        Attribute name, e.g. 'bary_date'.

    Returns
    -------
    :obj:`str`
        The full source line, stripped.
    """
    import pyodine.components as components

    source = Path(components.__file__).read_text()
    match = re.search(rf'^\s*{name}\s*=\s*None.*$', source, re.M)
    assert match is not None, f'no declaration of {name} found'
    return match.group(0).strip()


def test_bary_date_is_documented_as_a_full_jd_in_utc():
    """ Not the reduced barycentric date the comment claimed. """
    line = _declaration('bary_date')
    assert 'UTC' in line
    assert 'Reduced' not in line
    assert '2400000' not in line


def test_bary_vel_corr_is_documented_in_metres_per_second():
    """ `chunks.py` divides it by c in m/s. """
    line = _declaration('bary_vel_corr')
    assert 'm/s' in line
    assert 'km/s' not in line


def test_a_reduced_julian_date_is_rejected():
    """ The one unit error here that a machine can catch. """
    from pyodine.components import Observation

    obs = Observation()
    obs.bary_date = 51051.5          # MJD-like: a reduced date
    obs.bary_vel_corr = -12345.0
    with pytest.raises(ValueError):
        obs.check_bary_units()


def test_a_full_julian_date_passes():
    """ 1998-08-26, the night the HD 187123 template was taken. """
    from pyodine.components import Observation

    obs = Observation()
    obs.bary_date = 2451051.5
    obs.bary_vel_corr = -12345.0     # m/s; well within the +-30 km/s possible
    obs.check_bary_units()


def test_an_impossible_barycentric_velocity_is_rejected():
    """ The Earth does not move at 300 km/s. """
    from pyodine.components import Observation

    obs = Observation()
    obs.bary_date = 2451051.5
    obs.bary_vel_corr = 3.0e5
    with pytest.raises(ValueError):
        obs.check_bary_units()


def test_the_hires_adapter_satisfies_the_documented_contract():
    """ What phase 2 built must agree with what the fork now documents. """
    from first_hires_exoplanet.utilities_hires import load_pyodine

    assert 'JD(UTC)' in load_pyodine.__doc__
    assert '**m/s**' in load_pyodine.__doc__
