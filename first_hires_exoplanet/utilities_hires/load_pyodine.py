""" Present a PypeIt `spec1d` to `pyodine` as an Observation (phase 2, prompt 8).

Modelled on `pyodine`'s `utilities_lick/load_pyodine.py` at commit
`4488b0914fe5b272b787982647691045bff2604a`, and laid out as
`utilities_hires/` so it drops into a vendored `pyodine` tree unchanged.

This module builds the input only.  No forward model is attempted.

WHAT `pyodine` EXPECTS, and where each answer comes from:

  ``flux``           (nord, npix)   PypeIt `OPT_COUNTS`
  ``wave``           (nord, npix)   `OPT_WAVE`, divided by `VEL_CORR`
  ``cont``           (nord, npix)   running upper percentile of the flux
  ``weight``         (nord, npix)   inverse variance, zeroed where unusable
                                   (also per order as ``Spectrum.ivar``)
  ``bary_date``      scalar         full JD(UTC) at the exposure midpoint
  ``bary_vel_corr``  scalar         barycentric correction in **m/s**
  ``nord``, ``npix``, ``instrument``, ``star``, ``exp_time``

FOUR THINGS THIS ADAPTER HAS TO GET RIGHT, each established earlier in phase 2
and each capable of destroying the result silently:

1.  **Reference frame (prompt 1).**  PypeIt applied a *heliocentric* correction
    by default -- and its heliocentric calculation carries a sign error on the
    solar term worth +13 m/s, drifting 5 m/s over this baseline.  `pyodine`
    wants the **observed** frame: `chunks.py` shifts the template by the
    difference of the two barycentric corrections, which only works if neither
    spectrum has been corrected already.  So every wavelength here is divided
    by `VEL_CORR`, unconditionally, treating an absent value as 1.0.  That is
    exact -- `VEL_CORR` is one scalar per epoch and the round trip is lossless
    -- and it is correct whether the reduction used `refframe = heliocentric`
    or `observed`, which is the point of doing it unconditionally.

2.  **Units and epoch (prompt 1).**  `components.py` documents `bary_date` as a
    "Barycentric Reduced Julian Date (BJD - 2400000.0)" and `bary_vel_corr` as
    km/s.  Both docstrings are wrong.  `timeseries/bary_vel_corr.py` passes
    `bary_date` straight to `barycorrpy.get_BC_vel(JDUTC=...)`, which wants a
    **full JD in UTC**; and `chunks.py` divides `bary_vel_corr` by
    `astropy.constants.c` in **m/s**.  The header MJD is the exposure *start*
    (proved to 1.4 s from the recorded hour angle), so the midpoint has to be
    constructed here.

3.  **Quality (prompt 2).**  11% of order-spectra are unusable, almost all
    because `OPT_MASK` collapsed when the object profile spilled past the
    extraction aperture in poor seeing.  Those orders are given zero weight.
    Their flux is left alone: `components.Spectrum.__init__` raises
    `NoDataError` on an all-zero flux vector, so zeroing the flux would make
    the order unloadable rather than merely unused.

4.  **Vacuum, not air (prompt 7).**  PypeIt reports vacuum wavelengths.
    `pyodine`'s Lick `IodineTemplate` reads `wavelength_air` from the FTS
    atlas.  The two grids differ by 1.56 A at 5615 A, which is 83 km/s, and
    using the wrong one drops the correlation between the atlas and the
    measured HIRES cell from +0.90 to +0.08.  The `IodineTemplate` here reads
    `wavelength`.

`pyodine` is not installed in `pypeit14`, and prompt 8 is explicitly not the
place to vendor it.  The classes below therefore subclass `pyodine.components`
when it is importable and fall back to equivalent stand-ins when it is not, so
the adapter can be built and checked today and will subclass the real thing the
moment `pyodine` is vendored.  `HAVE_PYODINE` says which happened.

Self-check:

    conda run -n pypeit14 python -m first_hires_exoplanet.utilities_hires.load_pyodine

"""

# Standard imports
import os
import sys
import glob
import warnings

import numpy as np

from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation

from scipy.ndimage import percentile_filter

# The fork is source, not an installed package (see `vendor/README.md`), so it
# is imported by path as pyodine's own tutorial recommends.  Phase 2 wrote this
# module before the fork existed and fell back to stand-ins; since phase 3
# prompt 1 vendored it, the real classes are what these subclass.
_VENDOR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'vendor', 'pyodine')
if os.path.isdir(_VENDOR) and _VENDOR not in sys.path:
    sys.path.insert(0, _VENDOR)

try:                                                # pragma: no cover
    from pyodine import components
    HAVE_PYODINE = True
except ImportError:
    components = None
    HAVE_PYODINE = False


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
_PKG = os.path.dirname(_HERE)
DATA_ROOT = os.path.join(os.path.dirname(os.path.dirname(_PKG)),
                         'first-hires-exoplanet-data')
DEFAULT_REDUX = os.path.join(DATA_ROOT, 'redux')
DEFAULT_QUALITY = os.path.join(_PKG, 'data', 'order_quality.csv')
DEFAULT_ATLAS = os.path.join(DATA_ROOT, 'atlas',
                             'Fischer_Cell_May2022_downsampled3.h5')

#: Keck, from PypeIt's telescope parameters.  Also present per-file as
#: LON-OBS / LAT-OBS / ALT-OBS, which is what the code actually reads.
KECK_FALLBACK = dict(lon=-155.47833333333335, lat=19.82833333333333,
                     elv=4159.99999999977)

#: Window, in pixels, of the running continuum
CONT_WIN = 301
CONT_PCT = 95

#: Detector constants.  The raw headers give CCDGN01 = 4.8 e-/ADU and
#: CCDRN01 = 6.0 e-, while `keck_hires.py` hard-codes 1.9 and 2.8 and phase 1
#: deliberately did not change them.  The reduction that produced these spec1d
#: files therefore used PypeIt's numbers, so the inverse variance already in
#: the file is on that scale.  Both are recorded so phase 3 can decide; nothing
#: here rescales anything.
PYPEIT_GAIN, PYPEIT_RON = 1.9, 2.8

#: The Beer-Lambert exponent that turns the Fischer atlas into the HIRES cell
#: (prompt 7).  Applied by `IodineTemplate` when `scale_depth=True`, because
#: pyodine's own `iod_depth` is a LINEAR scaling and cannot express it: at this
#: exponent a saturated atlas line would go to negative transmission.
ATLAS_ALPHA = 2.59


class NoDataError(Exception):
    """ Raised for an order with no usable flux, as pyodine's does. """


if not HAVE_PYODINE:                                # pragma: no cover

    class _Spectrum:
        """ Stand-in for `pyodine.components.Spectrum`. """

        def __init__(self, flux, wave=None, cont=None):
            if not np.any(flux):
                raise NoDataError('Invalid flux vector!')
            self.flux, self.wave, self.cont = flux, wave, cont

        def __len__(self):
            return len(self.flux)

        def __getitem__(self, pixels):
            return _Spectrum(
                self.flux[pixels],
                None if self.wave is None else self.wave[pixels],
                None if self.cont is None else self.cont[pixels])

    class _Base:
        """ Stand-in for `components.Observation` / `components.IodineAtlas`. """

    _SpectrumCls, _ObsBase, _AtlasBase = _Spectrum, _Base, _Base
else:                                               # pragma: no cover
    _SpectrumCls = components.Spectrum
    _ObsBase = components.Observation
    _AtlasBase = components.IodineAtlas


# ---------------------------------------------------------------------------
# Small helpers, named as in utilities_lick
# ---------------------------------------------------------------------------

def or_none(header, key):
    """ Header value, or None when the card is absent. """
    return header[key] if key in header else None


class HIRES:
    """ The instrument, in the shape `pyodine` expects of one. """

    name = 'Keck/HIRES'
    longitude = KECK_FALLBACK['lon']
    latitude = KECK_FALLBACK['lat']
    altitude = KECK_FALLBACK['elv']

    def __init__(self, header=None):
        if header is not None:
            self.longitude = float(header.get('LON-OBS', self.longitude))
            self.latitude = float(header.get('LAT-OBS', self.latitude))
            self.altitude = float(header.get('ALT-OBS', self.altitude))
            self.detector = str(header.get('DETECTOR', '')).strip()
            self.decker = str(header.get('DECKER', '')).strip()
            self.binning = str(header.get('BINNING', '')).strip()
        # The pre-2004 single-Tektronix HIRES, as PypeIt's keck_hires_orig
        self.name = 'Keck/HIRES (orig)'

    def __repr__(self):
        return '<{:s}>'.format(self.name)


class Star:
    """ The star, in the shape `pyodine` expects of one. """

    def __init__(self, name, coordinates=None, pmra=np.nan, pmdec=np.nan,
                 rv0=np.nan):
        self.name = name
        self.coordinates = coordinates
        self.pmra, self.pmdec, self.rv0 = pmra, pmdec, rv0

    def __repr__(self):
        return '<Star {:s}>'.format(self.name)


def get_instrument(header):
    """ Build the instrument from a spec1d primary header. """
    return HIRES(header)


def get_star(header):
    """ Build the star from a spec1d primary header.

    PypeIt writes RA and DEC as decimal degrees, unlike the Lick loader's
    sexagesimal strings.  `TARGET` carries the observer's own label, which on
    these nights is variously `Star+Iodine`, `k24`, `k25`, `K17`, `K21` and
    `HD 187123`; the raw `TARGNAME` is `187123` throughout, with a typo of
    `H187123` on 1998-08-12.  Neither is a usable catalogue name, so the name
    is normalised here.
    """
    ra, dec = or_none(header, 'RA'), or_none(header, 'DEC')
    coord = None
    if ra is not None and dec is not None:
        coord = SkyCoord(ra=float(ra) * u.deg, dec=float(dec) * u.deg)
    return Star('HD 187123', coordinates=coord)


def continuum(flux, window=CONT_WIN, pct=CONT_PCT):
    """ Running upper-percentile continuum for one order. """
    win = min(window, (len(flux) // 2) * 2 - 1)
    return np.maximum(percentile_filter(flux, pct, size=win), 1e-9)


def bary_velocity(mjd_start, exptime, coord, lon, lat, elv):
    """ Barycentric correction at the exposure midpoint, in **m/s**.

    The header MJD is the exposure start (prompt 1), so the midpoint is
    constructed.  1998 HIRES had no exposure meter, so the geometric midpoint
    is the best available and leaves roughly 1 m/s of irreducible systematic.

    `barycorrpy` -- which `pyodine`'s own timeseries stage uses -- is not
    installed in `pypeit14`.  `astropy` is used instead; prompt 1 established
    that astropy's *barycentric* correction is the sound one and that PypeIt's
    heliocentric carries a sign error.  Phase 3 should install `barycorrpy` so
    that the initial guess and the final correction come from one source.
    """
    loc = EarthLocation.from_geodetic(lon=lon * u.deg, lat=lat * u.deg,
                                      height=elv * u.m)
    t_mid = Time(mjd_start + 0.5 * exptime / 86400., format='mjd', location=loc)
    return float(coord.radial_velocity_correction(
        kind='barycentric', obstime=t_mid).to(u.m / u.s).value)


def bary_julian_date(mjd_start, exptime):
    """ Full JD(UTC) at the exposure midpoint, which is what `bary_date` is.

    Not a BJD and not reduced, whatever `components.py` line 315 says: the
    value is passed to `barycorrpy.get_BC_vel(JDUTC=...)`, which converts it
    itself.  Handing it a reduced JD would place the observation in 1858.
    """
    return float(mjd_start + 0.5 * exptime / 86400. + 2400000.5)


# ---------------------------------------------------------------------------
# Quality flags
# ---------------------------------------------------------------------------

def load_quality(path=DEFAULT_QUALITY):
    """ Prompt 2's per-(epoch, order) quality table.

    Returns:
        dict: (koaid, order) -> bool usable.  Empty when the file is absent.
    """
    if not os.path.isfile(path):
        warnings.warn('No quality table at {:s}; every order will be '
                      'weighted.'.format(path))
        return {}
    from astropy.table import Table
    tbl = Table.read(path)
    # Written as 0/1 rather than True/False, because ascii.csv round-trips
    # booleans as the strings 'True'/'False', both of which are truthy.
    return {(str(r['koaid']), int(r['order'])): bool(int(r['use']))
            for r in tbl}


# ---------------------------------------------------------------------------
# Reading a spec1d
# ---------------------------------------------------------------------------

def load_file(filename, quality=None):
    """ Read a PypeIt spec1d into the rectangular arrays `pyodine` wants.

    Read with plain `astropy.io.fits` rather than `pypeit.specobjs`, so the
    adapter has no PypeIt dependency and can live inside a vendored `pyodine`.

    Args:
        filename (str): a PypeIt spec1d file.
        quality (dict, optional): output of :func:`load_quality`.

    Returns:
        tuple: (flux, wave, cont, weight, orders, header), the first four
        shaped (nord, npix).
    """
    quality = {} if quality is None else quality
    with fits.open(filename) as hdul:
        header = hdul[0].header
        koaid = str(header.get('FILENAME', '')).replace('.fits', '')

        rows = []
        for hdu in hdul[1:]:
            if not hasattr(hdu, 'columns') or 'OPT_WAVE' not in hdu.columns.names:
                continue
            order = int(hdu.header['ECH_ORDER'])
            # VEL_CORR is per-object and identical across orders, but read it
            # per extension so a file that somehow disagrees is handled, not
            # assumed away.  Absent means the reduction ran refframe=observed.
            vel_corr = float(hdu.header.get('VEL_CORR', 1.0) or 1.0)
            d = hdu.data
            rows.append((order, vel_corr,
                         np.asarray(d['OPT_WAVE'], dtype=float),
                         np.asarray(d['OPT_COUNTS'], dtype=float),
                         np.asarray(d['OPT_COUNTS_IVAR'], dtype=float),
                         np.asarray(d['OPT_MASK'], dtype=bool)))

    if not rows:
        raise NoDataError('No echelle orders in {:s}'.format(filename))
    rows.sort(key=lambda r: r[0])

    npix = len(rows[0][2])
    nord = len(rows)
    flux = np.zeros((nord, npix))
    wave = np.zeros((nord, npix))
    cont = np.ones((nord, npix))
    weight = np.zeros((nord, npix))
    orders = np.array([r[0] for r in rows], dtype=int)

    for i, (order, vel_corr, w, f, iv, m) in enumerate(rows):
        # (1) back to the observed frame -- see the module docstring
        wave[i] = w / vel_corr
        flux[i] = f
        # Masked pixels arrive as exact zeros, not flagged gaps (phase 1)
        good = m & (iv > 0) & (w > 0) & np.isfinite(f) & (f != 0)
        if good.sum() > 50:
            cont[i] = continuum(np.where(good, f, np.median(f[good])))
        # (3) inverse variance, zeroed on bad pixels and on orders prompt 2
        # rejects.  The flux is deliberately left intact: an all-zero flux
        # vector makes components.Spectrum raise instead of merely being
        # ignored.
        usable = quality.get((koaid, order), True)
        weight[i] = np.where(good & usable, iv, 0.)

    return flux, wave, cont, weight, orders, header


# ---------------------------------------------------------------------------
# The wrapper
# ---------------------------------------------------------------------------

class ObservationWrapper(_ObsBase):
    """ A PypeIt spec1d of HD 187123, in `pyodine`'s Observation shape.

    Args:
        filename (str): the spec1d file.
        instrument: override the instrument built from the header.
        star: override the star built from the header.
        quality (dict, optional): prompt 2's flags; loaded by default.
    """

    _flux = None
    _wave = None
    _cont = None
    _weight = None

    def __init__(self, filename, instrument=None, star=None, quality=None):
        if quality is None:
            quality = load_quality()
        flux, wave, cont, weight, orders, header = load_file(filename, quality)

        self._flux, self._wave, self._cont, self._weight = flux, wave, cont, weight
        # NOT `self.orders`: `components.MultiOrderSpectrum` defines that as a
        # read-only property giving *positions* (0..nord-1), which is what
        # pyodine's loops index with.  Phase 2 wrote this against a stand-in
        # base class and reused the name for the echelle order NUMBERS (57-93);
        # the collision only surfaced when phase 3 prompt 2 put the real fork
        # on the path.  The physical numbers live under their own name.
        self.ech_orders = orders
        self.nord, self.npix = flux.shape

        self.orig_header = header
        self.orig_filename = os.path.abspath(filename)
        self.koaid = str(header.get('FILENAME', '')).replace('.fits', '')

        self.instrument = instrument or get_instrument(header)
        self.star = star or get_star(header)

        # The iodine state is not in the spec1d; it is in the raw frame.
        self.iodine_in_spectrum = self._iodine_state()
        self.iodine_cell_id = 'keck_hires' if self.iodine_in_spectrum else None

        self.exp_time = float(header['EXPTIME'])
        self.flux_level = float(np.median(flux[weight > 0])) \
            if np.any(weight > 0) else np.nan
        # PypeIt's values, which are what produced the inverse variance in this
        # file, not the raw header's 4.8 / 6.0.  See the note on the constants.
        self.gain = PYPEIT_GAIN
        self.readout_noise = PYPEIT_RON
        self.dark_current = None

        mjd = float(header['MJD'])
        self.time_start = Time(mjd, format='mjd', scale='utc')
        # No exposure meter on 1998 HIRES, so the flux-weighted midpoint is not
        # recoverable and the geometric one is the best available.
        self.time_weighted = Time(mjd + 0.5 * self.exp_time / 86400.,
                                  format='mjd', scale='utc')

        # (2) full JD(UTC) at the midpoint, and m/s -- not BJD, not km/s
        self.bary_date = bary_julian_date(mjd, self.exp_time)
        self.bary_vel_corr = bary_velocity(
            mjd, self.exp_time, self.star.coordinates,
            self.instrument.longitude, self.instrument.latitude,
            self.instrument.altitude)

    def _iodine_state(self):
        """ Was the cell in the beam?  Read from the raw frame if we have it. """
        pattern = os.path.join(DATA_ROOT, 'raw', '*', self.koaid + '.fits')
        hits = sorted(glob.glob(pattern))
        if not hits:
            warnings.warn('No raw frame for {:s}; iodine state unknown, '
                          'assuming cell IN.'.format(self.koaid))
            return True
        return bool(fits.getheader(hits[0]).get('IODIN'))

    def __len__(self):
        return self.nord

    def __getitem__(self, order):
        """ One order, several, or a slice -- indexed by POSITION, not by
        echelle order number, which is what `pyodine` iterates over.  Use
        :meth:`index_of` to go from an echelle order number to a position.
        """
        if isinstance(order, (int, np.integer)) or hasattr(order, '__int__'):
            i = int(order)
            if HAVE_PYODINE:                        # pragma: no cover
                return _SpectrumCls(self._flux[i], wave=self._wave[i],
                                    cont=self._cont[i], ivar=self._weight[i])
            return _SpectrumCls(self._flux[i], wave=self._wave[i],
                                cont=self._cont[i])
        if isinstance(order, (list, np.ndarray)):
            return [self.__getitem__(int(i)) for i in order]
        if isinstance(order, slice):
            return self.__getitem__([int(i) for i in np.arange(self.nord)[order]])
        raise IndexError(type(order))

    def index_of(self, ech_order):
        """ Position in the arrays of a given echelle order number. """
        hits = np.where(self.ech_orders == int(ech_order))[0]
        if len(hits) != 1:
            raise IndexError('Order {:d} is not in {:s}'.format(
                int(ech_order), self.koaid))
        return int(hits[0])

    def compute_weight(self, weight_type='inverse-variance', rel_noise=0.008):
        """ Pixel weights.

        `pyodine` originally offered 'flat' (ones) and 'inverse'
        (1/(f(1+f*rel_noise^2)), from the dop code).  Neither knows that PypeIt
        has already propagated a full inverse variance through the extraction,
        nor that prompt 2 rejected 11% of order-spectra.  The default here is
        that inverse variance, with rejected orders and masked pixels at zero;
        the two `pyodine` options remain available for comparison.

        Phase 3 prompt 2 taught the fork the same thing, where it is spelled
        ``ivar``.  That spelling is accepted here too, so the same string
        means the same array on both sides of the boundary; the longer name
        this module was written with in phase 2 still works.
        """
        if weight_type in ('inverse-variance', 'ivar'):
            return self._weight.copy()
        if weight_type == 'flat':
            return np.ones_like(self._flux)
        if weight_type == 'inverse':
            with np.errstate(divide='ignore', invalid='ignore'):
                w = 1. / (self._flux * (1. + self._flux * rel_noise ** 2))
            return np.where(np.isfinite(w), w, 0.)
        raise NotImplementedError(
            'Choose one of: "ivar" (= "inverse-variance"), "flat", "inverse"')

    @property
    def usable_orders(self):
        """ Echelle order numbers with any non-zero weight. """
        return self.ech_orders[np.any(self._weight > 0, axis=1)]

    def __repr__(self):
        return ('<HIRES {:s}  {:d} orders x {:d} pix  {:.0f} s  '
                'I2 {:s}  BVC {:+.3f} km/s>'.format(
                    self.koaid, self.nord, self.npix, self.exp_time,
                    'in' if self.iodine_in_spectrum else 'out',
                    self.bary_vel_corr / 1e3))


# ---------------------------------------------------------------------------
# The atlas
# ---------------------------------------------------------------------------

class IodineTemplate(_AtlasBase):
    """ The FTS iodine atlas, on the VACUUM grid and scaled to the HIRES cell.

    Two deliberate departures from `utilities_lick.IodineTemplate`, both
    established in prompt 7:

    * it reads ``wavelength``, not ``wavelength_air``.  PypeIt reports vacuum,
      the two grids differ by 83 km/s, and using the air one drops the
      correlation against the measured HIRES cell from +0.90 to +0.08;
    * it optionally raises the transmission to ``ATLAS_ALPHA``.  The HIRES cell
      is 2.6x optically thicker than the cell this atlas scanned, and the
      difference is a single Beer-Lambert exponent.  `pyodine`'s `iod_depth`
      applies depth *linearly*, which at this exponent sends 5.8% of the atlas
      to negative transmission, so it cannot absorb the difference itself.

    Args:
        iodine_cell (str): path to the atlas HDF5 file.
        scale_depth (bool): apply the Beer-Lambert rescaling.
    """

    def __init__(self, iodine_cell=DEFAULT_ATLAS, scale_depth=True):
        import h5py

        self.orig_filename = iodine_cell
        with h5py.File(iodine_cell, 'r') as h:
            flux = np.array(h['flux_normalized'], dtype=float)
            # 'wavelength' is the vacuum grid; 'wavelength_air' is the one
            # upstream's adapters read.  Phase 3 prompt 2 gave the fork's
            # `IodineAtlas` a recorded `wave_frame` so that the choice is
            # visible to anything downstream instead of being implicit here.
            wave = np.array(h[components.IodineAtlas.WAVE_FRAMES['vacuum']]
                            if HAVE_PYODINE else h['wavelength'], dtype=float)

        self.alpha = ATLAS_ALPHA if scale_depth else 1.0
        if scale_depth:
            flux = np.clip(flux, 0., None) ** self.alpha

        if HAVE_PYODINE:                            # pragma: no cover
            super().__init__(flux, wave, wave_frame='vacuum')
        else:
            self.flux, self.wave, self.cont = flux, wave, None
            self.wave_frame = 'vacuum'

    def __repr__(self):
        return '<IodineTemplate {:.1f}-{:.1f} A vacuum, alpha={:.2f}>'.format(
            self.wave.min(), self.wave.max(), self.alpha)


# ---------------------------------------------------------------------------
# Self-check
# ---------------------------------------------------------------------------

def _self_check(redux_dir=DEFAULT_REDUX):           # pragma: no cover
    """ Load every reduced epoch and verify the contract. """
    files = sorted(glob.glob(os.path.join(redux_dir, 'reduce_19*', 'Science',
                                          'spec1d_*.fits')))
    print('\n=== pyodine input adapter, self-check ===')
    print('  pyodine importable : {!s}'.format(HAVE_PYODINE))
    print('  spec1d files found : {:d}'.format(len(files)))
    quality = load_quality()
    print('  quality flags      : {:d} (epoch, order) pairs'.format(len(quality)))

    obs = [ObservationWrapper(f, quality=quality) for f in files]
    print('\n  {:18s} {:>5s} {:>5s} {:>7s} {:>5s} {:>7s} {:>14s} {:>9s}'.format(
        'koaid', 'nord', 'npix', 'exp', 'I2', 'usable', 'bary_date', 'BVC m/s'))
    for o in obs[:6] + obs[-3:]:
        print('  {:18s} {:5d} {:5d} {:7.0f} {:>5s} {:7d} {:14.5f} {:+9.1f}'.format(
            o.koaid, o.nord, o.npix, o.exp_time,
            'in' if o.iodine_in_spectrum else 'out',
            len(o.usable_orders), o.bary_date, o.bary_vel_corr))
    print('  ... {:d} epochs total'.format(len(obs)))

    print('\n  Contract checks')
    ok = True

    shapes = {(o.nord, o.npix) for o in obs}
    print('    rectangular (nord, npix)      : {!s}'.format(shapes))
    ok &= all(o._flux.shape == o._wave.shape == o._weight.shape for o in obs)

    # The observed frame: undoing VEL_CORR must move the wavelengths
    f0 = files[0]
    with fits.open(f0) as h:
        vc = float(h[1].header.get('VEL_CORR', 1.0) or 1.0)
        raw = np.asarray(h[1].data['OPT_WAVE'], dtype=float)
    i = obs[0].index_of(int(fits.getheader(f0, 1)['ECH_ORDER']))
    shift = (obs[0]._wave[i][raw > 0] / raw[raw > 0] - 1.).mean() * 299792.458
    print('    VEL_CORR removed              : {:+.3f} km/s (expect '
          '{:+.3f})'.format(shift, (1. / vc - 1.) * 299792.458))
    ok &= abs(shift - (1. / vc - 1.) * 299792.458) < 1e-6

    print('    bary_date is a full JD        : {!s}'.format(
        all(2400000. < o.bary_date < 2500000. for o in obs)))
    ok &= all(2400000. < o.bary_date < 2500000. for o in obs)

    bvc = np.array([o.bary_vel_corr for o in obs])
    print('    bary_vel_corr in m/s          : {:+.0f} to {:+.0f}'.format(
        bvc.min(), bvc.max()))
    ok &= bool(np.all(np.abs(bvc) < 4.0e4) and np.any(np.abs(bvc) > 1.0e3))

    print('    midpoint, not start           : {:.1f} s after MJD'.format(
        (obs[0].bary_date - (float(obs[0].orig_header['MJD']) + 2400000.5))
        * 86400.))

    zero = sum(int(np.sum(~np.any(o._weight > 0, axis=1))) for o in obs)
    print('    orders zero-weighted          : {:d} of {:d}'.format(
        zero, sum(o.nord for o in obs)))

    # No order may have been made unloadable
    bad = []
    for o in obs:
        for k in range(o.nord):
            try:
                o[k]
            except Exception as exc:
                bad.append((o.koaid, int(o.ech_orders[k]), str(exc)))
    print('    every order still loadable    : {!s}{:s}'.format(
        not bad, '' if not bad else '  <- {:d} failed'.format(len(bad))))
    ok &= not bad

    w = obs[0].compute_weight()
    print('    compute_weight default        : inverse variance, '
          '{:.1%} of pixels non-zero'.format(float(np.mean(w > 0))))
    ok &= w.shape == obs[0]._flux.shape

    if os.path.isfile(DEFAULT_ATLAS):
        try:
            atlas = IodineTemplate()
            print('    atlas (vacuum, scaled)        : {!r}'.format(atlas))
        except ImportError:
            print('    atlas                         : h5py not in this env; '
                  'run in `astro` to check')
    else:
        print('    atlas                         : not on disk, skipped')

    print('\n  {:s}'.format('ALL CHECKS PASS' if ok else 'SOME CHECKS FAILED'))
    return ok


if __name__ == '__main__':                          # pragma: no cover
    _self_check()
