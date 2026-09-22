""" What calibrations exist for each discovery-era night (phase 2, prompt 3).

Phase 1 reduced five nights of one July run, and every one of them borrowed its
flats from 1998-07-14 because those were the only clean B1 flats within PypeIt's
angle tolerances of that run.  Prompt 3 extends the reduction over nine months,
so the question has to be asked again for each night rather than assumed: does
this night have its own arcs?  Its own flats?  If not, which night's will
PypeIt actually accept?

PypeIt groups frames into configurations by binning, decker, dispname,
echangle, xdangle and filter1 (`keck_hires.py`).  The angle comparisons carry
tolerances -- echangle `rtol=1e-3, atol=1e-2`, xdangle `rtol=1e-2, atol=1e-1` --
so a flat taken on another night is usable only if its echelle angle is within
about 0.01 of the science frame's.  Across the discovery era the echelle angle
runs from +0.0149 (1997-12) to -0.0051 (1998-08), a span of 0.02, so this is a
real constraint and not a formality: the December 1997 frames cannot share
calibrations with the August 1998 ones.

This module reads the archive survey already in the repository
(`data/koa_hd187123_nights_all.csv`, 3124 frames over 21 nights) and reports,
per night, the frames that are in the planet-search configuration and what role
each can play.  It queries nothing and downloads nothing.

It also checks the frames on disk against what the archive says about them,
because on this dataset the two disagree.  Five of 1998-08-12's eighteen
"Narrowflat" frames are recorded by KOA as `flatlamp`, iodine out, hatch
closed -- a perfect clean flat by every column in the survey -- and contain
nothing at all: the cross-disperser cover (`XCOVOPEN`) was shut, so they hold
bias and read noise, mean 768 counts with a standard deviation of 1.4 against
17,400 and 12,000 for a real one.  PypeIt declines to frametype them and the
reduction fails several minutes later with "No frames of type=trace provided",
which says nothing about the cause.  Calibrations cannot be chosen from
metadata alone.

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.calib_inventory

"""

# Standard imports
import os
import argparse

import numpy as np

from astropy.table import Table


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SURVEY = os.path.join(_HERE, 'data', 'koa_hd187123_nights_all.csv')

#: The Marcy planet-search configuration, as established in phase 1 prompt 3.
#: Anything outside it belongs to another program sharing the night.
CONFIG = dict(binning='1,2', xdispers='RED', fil1name='kv370',
              numamps='1', ccdgain='F')

#: PypeIt's configuration tolerances for `keck_hires`
#: (`pypeit/spectrographs/keck_hires.py`, the `meta` definitions)
ECH_RTOL, ECH_ATOL = 1e-3, 1e-2
XD_RTOL, XD_ATOL = 1e-2, 1e-1


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def in_config(tbl):
    """ Boolean mask: frames in the planet-search configuration.

    Args:
        tbl (`astropy.table.Table`_): the archive survey.

    Returns:
        `numpy.ndarray`_: boolean mask.
    """
    keep = np.ones(len(tbl), dtype=bool)
    for key, val in CONFIG.items():
        keep &= np.array([str(x).strip() == val for x in tbl[key]])
    return keep


def classify(row):
    """ What role can this frame play in a reduction?

    KOA's `koaimtyp` is the archive's own classification and is not always what
    PypeIt will decide, so the lamp and hatch keywords are used directly.

    Args:
        row (`astropy.table.Row`_): one frame.

    Returns:
        str: 'science', 'arc', 'flat', 'bias', or 'other'.
    """
    imtyp = str(row['koaimtyp']).strip()
    if imtyp == 'object' and str(row['targname']).strip() == '187123':
        return 'science'
    if str(row['lampcat1']).strip() == 'T' or str(row['lampcat2']).strip() == 'T':
        return 'arc'
    if imtyp in ('flatlamp', 'trace') or str(row['lampqtz2']).strip() == 'T':
        # A flat with the hatch OPEN sees the sky as well as the lamp, and a
        # flat with the iodine cell IN carries the I2 forest into the flat
        # field.  Phase 1 rejected both; the distinction is kept here.
        clean = (str(row['hatopen']).strip() == 'F'
                 and str(row['iodin']).strip() != 'T')
        return 'flat' if clean else 'flat_dirty'
    if imtyp in ('bias', 'dark'):
        return 'bias'
    return 'other'


def angles_compatible(ech_a, ech_b, xd_a, xd_b):
    """ Would PypeIt place two frames in the same configuration, on angles?

    Args:
        ech_a, ech_b (float): echelle angles.
        xd_a, xd_b (float): cross-disperser angles.

    Returns:
        bool
    """
    ech_ok = np.abs(ech_a - ech_b) <= ECH_ATOL + ECH_RTOL * np.abs(ech_b)
    xd_ok = np.abs(xd_a - xd_b) <= XD_ATOL + XD_RTOL * np.abs(xd_b)
    return bool(ech_ok and xd_ok)


# ---------------------------------------------------------------------------
# Inventory
# ---------------------------------------------------------------------------

def build(survey_file):
    """ Per-night inventory of usable frames.

    Args:
        survey_file (str): the archive survey CSV.

    Returns:
        tuple: (per-frame table restricted to the configuration, per-night table)
    """
    tbl = Table.read(survey_file)
    tbl['night'] = [str(d)[:10] for d in tbl['date_obs']]
    tbl = tbl[in_config(tbl)]
    tbl['role'] = [classify(r) for r in tbl]

    rows = []
    for night in sorted(set(tbl['night'])):
        sel = tbl[tbl['night'] == night]
        sci = sel[sel['role'] == 'science']
        if len(sci) == 0:
            continue
        deck = str(sci['deckname'][0]).strip()
        same_deck = np.array([str(d).strip() == deck for d in sel['deckname']])
        counts = {r: int(np.sum(sel['role'] == r))
                  for r in ['science', 'arc', 'flat', 'flat_dirty', 'bias']}
        # A flat is only useful if it was taken through the science decker.
        # Most nights of this program flat-fielded through B2 while observing
        # through B1, so the undifferentiated count is badly misleading.
        counts['flat'] = int(np.sum((sel['role'] == 'flat') & same_deck))
        counts['arc'] = int(np.sum((sel['role'] == 'arc') & same_deck))
        rows.append(dict(
            night=night,
            decker='/'.join(sorted(set(str(d).strip() for d in sci['deckname']))),
            n_sci=counts['science'],
            n_arc=counts['arc'],
            n_flat=counts['flat'],
            n_flat_otherdeck=int(np.sum(sel['role'] == 'flat')) - counts['flat'],
            n_flat_dirty=counts['flat_dirty'],
            n_bias=counts['bias'],
            ech_sci=float(np.median(sci['echangl'])),
            xd_sci=float(np.median(sci['xdangl'])),
        ))
    return tbl, Table(rows)


def flat_donors(tbl, nights):
    """ For each night, which nights could legally supply clean flats?

    Args:
        tbl (`astropy.table.Table`_): per-frame table from :func:`build`.
        nights (`astropy.table.Table`_): per-night table from :func:`build`.

    Returns:
        dict: night -> list of (donor night, n_flats, delta echangle)
    """
    flats = tbl[tbl['role'] == 'flat']
    out = {}
    for row in nights:
        donors = []
        for donor in sorted(set(flats['night'])):
            cand = flats[flats['night'] == donor]
            # A donor night's flats must match on decker too
            cand = cand[[str(d).strip() == row['decker'] for d in cand['deckname']]]
            if len(cand) == 0:
                continue
            ok = [angles_compatible(float(c['echangl']), row['ech_sci'],
                                    float(c['xdangl']), row['xd_sci'])
                  for c in cand]
            if any(ok):
                donors.append((donor, int(np.sum(ok)),
                               float(np.median(cand['echangl'])) - row['ech_sci']))
        out[row['night']] = donors
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def verify_flats(raw_root):
    """ Check the clean flats on disk against what the archive claims.

    Args:
        raw_root (str): the raw data tree.

    Returns:
        `astropy.table.Table`_: one row per candidate flat found on disk.
    """
    import glob
    from astropy.io import fits

    rows = []
    for path in sorted(glob.glob(os.path.join(raw_root, '*', 'HI.*.fits'))):
        hdr = fits.getheader(path)
        if str(hdr.get('TARGNAME', '')).strip() == '187123':
            continue
        if hdr.get('LAMPCAT1') or hdr.get('LAMPCAT2'):
            continue
        if not (hdr.get('LAMPQTZ2')
                or str(hdr.get('LAMPNAME', '')).strip() == 'quartz1'):
            continue
        if hdr.get('IODIN') or hdr.get('HATOPEN'):
            continue
        rows.append(dict(
            frame=os.path.basename(path),
            night=str(hdr.get('DATE-OBS', '')).strip(),
            decker=str(hdr.get('DECKNAME', '')).strip(),
            ech=float(hdr.get('ECHANGL', np.nan)),
            xcovopen=bool(hdr.get('XCOVOPEN')),
            mean=float(hdr.get('IM01MN01', np.nan)),
            sigma=float(hdr.get('IM01SD01', np.nan)),
        ))
    return Table(rows) if rows else None


def parse_args(options=None):
    parser = argparse.ArgumentParser(
        description='Calibration inventory for the discovery-era nights.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--survey', type=str, default=DEFAULT_SURVEY,
                        help='Archive survey CSV')
    parser.add_argument('--raw', type=str, default=os.path.join(
                            os.path.dirname(os.path.dirname(_HERE)),
                            'first-hires-exoplanet-data', 'raw'),
                        help='Raw data tree, for the on-disk flat check')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):

    tbl, nights = build(pargs.survey)

    print('\n=== 1. Frames in the planet-search configuration, per night ===')
    print('  arc, flat = taken through the SAME decker as that night\'s science')
    print('  otherdk   = clean flats taken through a different decker, which')
    print('              PypeIt will not use for this configuration')
    print('  flat      = quartz, hatch closed, iodine cell out (usable)')
    print('  flat_dty  = quartz, but hatch open or iodine cell in (rejected')
    print('              in phase 1 prompt 4)')
    print()
    print(f'  {"night":12s} {"deck":5s} {"sci":>4s} {"arc":>4s} {"flat":>5s} '
          f'{"otherdk":>7s} {"f_dty":>6s} {"echangle":>9s}')
    for row in nights:
        print(f'  {row["night"]:12s} {row["decker"]:5s} {row["n_sci"]:4d} '
              f'{row["n_arc"]:4d} {row["n_flat"]:5d} {row["n_flat_otherdeck"]:7d} '
              f'{row["n_flat_dirty"]:6d} {row["ech_sci"]:9.5f}')

    print('\n=== 2. Nights with no arc of their own ===')
    for row in nights:
        if row['n_arc'] == 0:
            print(f'  {row["night"]}  (echangle {row["ech_sci"]:+.5f})')

    print('\n=== 3. Where clean flats of each decker exist at all ===')
    flats = tbl[tbl['role'] == 'flat']
    for deck in ['B1', 'B2']:
        sel = flats[[str(d).strip() == deck for d in flats['deckname']]]
        by_night = [(n, int(np.sum(sel['night'] == n)))
                    for n in sorted(set(sel['night']))]
        print(f'  {deck}: {len(by_night)} nights -- '
              + ', '.join(f'{n} ({c})' for n, c in by_night))

    print('\n=== 4. Which nights can legally supply each night\'s flats ===')
    print('  PypeIt accepts a flat only if the deckers match and')
    print(f'  |d echangle| <= {ECH_ATOL} + {ECH_RTOL}*|ech|.')
    donors = flat_donors(tbl, nights)
    for row in nights:
        got = donors[row['night']]
        own = ' [own]' if row['n_flat'] >= 4 else ''
        if len(got) == 0:
            print(f'  {row["night"]} ({row["decker"]}): NO legal donor')
            continue
        best = ', '.join(f'{d} ({n}, {de:+.4f})' for d, n, de in got)
        print(f'  {row["night"]} ({row["decker"]}): {best}{own}')

    print('\n=== 5. Echelle-angle span of the discovery era ===')
    print(f'  {np.min(nights["ech_sci"]):+.5f} to {np.max(nights["ech_sci"]):+.5f}'
          f'  (span {np.ptp(nights["ech_sci"]):.5f}, tolerance {ECH_ATOL})')
    print('  => nights more than 0.01 apart in echelle angle cannot share')
    print('     calibrations, so there is no single flat set for the era.')

    print('\n=== 6. Do the flats on disk contain anything? ===')
    flats = verify_flats(pargs.raw)
    if flats is None:
        print(f'  No raw frames under {pargs.raw}; skipping.')
        return
    dead = flats[~np.asarray(flats['xcovopen'], dtype=bool)]
    live = flats[np.asarray(flats['xcovopen'], dtype=bool)]
    print(f'  {len(flats)} clean flats on disk by every archive column')
    print(f'    cross-disperser cover OPEN : {len(live):3d}  '
          f'sigma {np.min(live["sigma"]):.0f}-{np.max(live["sigma"]):.0f}')
    if len(dead):
        print(f'    cross-disperser cover SHUT : {len(dead):3d}  '
              f'sigma {np.min(dead["sigma"]):.1f}-{np.max(dead["sigma"]):.1f}'
              f'  <- EMPTY, unusable')
        for night in sorted(set(dead['night'])):
            n = int(np.sum(dead['night'] == night))
            tot = int(np.sum(flats['night'] == night))
            print(f'      {night}: {n} of {tot}')
    else:
        print('    cross-disperser cover SHUT :   0')


if __name__ == '__main__':
    main(parse_args())
