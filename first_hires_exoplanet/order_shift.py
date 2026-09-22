""" How far do the echelle orders move between nights? (phase 2, prompt 3)

1997-12-23 and 1997-12-24 have no iodine-free B1 flat of their own, and their
echelle angles (+0.0149, +0.0160) put them more than PypeIt's 0.01 tolerance
away from every clean B1 flat in the discovery era.  PypeIt therefore refuses
to group any flat with them and the reduction dies with "No frames of
type=trace provided".

The question that decides whether a flat from another night can be forced in is
not what PypeIt's tolerance says -- that is a grouping heuristic, not a physical
law -- but how far the orders actually move on the detector between those
configurations.  A trace frame is only useful if it puts the order edges where
the science frame's orders really are.

The cross-disperser angle sets the spatial position of the orders and barely
moves across the era (-0.538 to -0.562).  The echelle angle sets which
wavelengths land in each order, and because the cross-disperser's dispersion
depends on wavelength, changing it moves the orders spatially as well.  How
much is measurable, and this measures it.

Method, following `redux/checks/order_shift_jul14_vs_jul16.py` from phase 1:
collapse a frame along the dispersion axis to get the spatial order profile,
then cross-correlate one night's profile against another's.  ThAr arcs are used
because every discovery-era night has one, including the two December nights,
and they illuminate every order.

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.order_shift

"""

# Standard imports
import os
import glob
import argparse

import numpy as np

from astropy.io import fits
from astropy.table import Table


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RAW = os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                           'first-hires-exoplanet-data', 'raw')

#: Columns to collapse over.  A central strip, away from the blaze roll-off at
#: either end of each order (phase 1, prompt 5).
STRIP = (800, 1300)

#: Largest shift searched, in binned spatial pixels.  The B1 orders are only
#: ~10.3 binned pixels wide, so a shift beyond a few pixels already means the
#: trace lands on the wrong part of the order.
MAXLAG = 120


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

def spatial_profile(path):
    """ Collapse a frame along dispersion into its spatial order profile.

    Args:
        path (str): a raw FITS frame.

    Returns:
        `numpy.ndarray`_: median-subtracted profile along the spatial axis.
    """
    img = fits.getdata(path, memmap=False).astype(float)
    prof = np.median(img[:, STRIP[0]:STRIP[1]], axis=1)
    return prof - np.median(prof)


def xcorr_shift(ref, other, maxlag=MAXLAG):
    """ Shift, in binned spatial pixels, that best aligns `other` onto `ref`.

    Args:
        ref, other (`numpy.ndarray`_): spatial profiles.
        maxlag (int): largest shift searched.

    Returns:
        tuple: (shift, normalised peak correlation)
    """
    lags = np.arange(-maxlag, maxlag + 1)
    cc = np.array([
        np.sum(ref[max(0, l):len(ref) + min(0, l)]
               * other[max(0, -l):len(other) + min(0, -l)])
        for l in lags])
    norm = np.sqrt(np.sum(ref * ref) * np.sum(other * other))
    return int(lags[np.argmax(cc)]), float(cc.max() / norm)


def find_arcs(raw_root, decker='B1'):
    """ One ThAr arc per night, through the given decker.

    Args:
        raw_root (str): the raw data tree.
        decker (str): decker to require.

    Returns:
        `astropy.table.Table`_: one row per night.
    """
    rows = []
    for night_dir in sorted(glob.glob(os.path.join(raw_root, '199*'))):
        best = None
        for path in sorted(glob.glob(os.path.join(night_dir, 'HI.*.fits'))):
            hdr = fits.getheader(path)
            if not (hdr.get('LAMPCAT1') or hdr.get('LAMPCAT2')):
                continue
            if str(hdr.get('DECKNAME', '')).strip() != decker:
                continue
            # Prefer the longest exposure: more lines, a cleaner correlation.
            if best is None or float(hdr.get('ELAPTIME', 0)) > best[1]:
                best = (path, float(hdr.get('ELAPTIME', 0)), hdr)
        if best is None:
            continue
        path, elap, hdr = best
        rows.append(dict(night=os.path.basename(night_dir), frame=os.path.basename(path),
                         path=path, exptime=elap,
                         ech=float(hdr.get('ECHANGL', np.nan)),
                         xd=float(hdr.get('XDANGL', np.nan))))
    return Table(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(options=None):
    parser = argparse.ArgumentParser(
        description='Measure the spatial order shift between nights.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--raw', type=str, default=DEFAULT_RAW,
                        help='Raw data tree')
    parser.add_argument('--reference', type=str, default='1998jul14',
                        help='Night whose flats are the phase-1 trace source')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):

    arcs = find_arcs(pargs.raw)
    if pargs.reference not in list(arcs['night']):
        raise RuntimeError('No {:s} arc found'.format(pargs.reference))
    ref_row = arcs[arcs['night'] == pargs.reference][0]
    ref = spatial_profile(ref_row['path'])

    print('\n=== Spatial order shift against {:s} ==='.format(pargs.reference))
    print('  reference: {:s}  echangle {:+.5f}  xdangle {:+.4f}'.format(
        ref_row['frame'], ref_row['ech'], ref_row['xd']))
    print('  Both are B1 ThAr arcs, collapsed over columns '
          '{:d}-{:d} and cross-correlated.'.format(*STRIP))
    print()
    print('  {:12s} {:>10s} {:>9s} {:>7s} {:>7s}  {:s}'.format(
        'night', 'd_echangle', 'd_xdangle', 'shift', 'corr', 'within PypeIt tol?'))

    out = []
    for row in arcs:
        prof = spatial_profile(row['path'])
        shift, corr = xcorr_shift(ref, prof)
        dech = row['ech'] - ref_row['ech']
        dxd = row['xd'] - ref_row['xd']
        legal = abs(dech) <= 1e-2 + 1e-3 * abs(ref_row['ech'])
        print('  {:12s} {:+10.5f} {:+9.4f} {:7d} {:7.3f}  {:s}'.format(
            row['night'], dech, dxd, shift, corr, 'yes' if legal else 'NO'))
        out.append(dict(night=row['night'], dech=dech, dxd=dxd, shift=shift,
                        corr=corr, legal=legal))

    tbl = Table(out)
    inside = tbl[np.asarray(tbl['legal'], dtype=bool)]
    outside = tbl[~np.asarray(tbl['legal'], dtype=bool)]
    print('\n  Nights PypeIt would accept : shift '
          '{:+d} to {:+d} px'.format(int(np.min(inside['shift'])),
                                     int(np.max(inside['shift']))))
    if len(outside):
        print('  Nights PypeIt rejects      : shift '
              '{:+d} to {:+d} px'.format(int(np.min(outside['shift'])),
                                         int(np.max(outside['shift']))))
    print('\n  The B1 orders are ~10.3 binned pixels wide, so a shift of more')
    print('  than ~2 px puts a borrowed trace on the wrong part of the order.')
    return tbl


if __name__ == '__main__':
    main(parse_args())
