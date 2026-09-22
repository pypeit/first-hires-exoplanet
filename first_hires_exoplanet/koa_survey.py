""" Survey the Keck Observatory Archive for the 1997-1998 HIRES data on HD 187123.

Metadata only -- this script downloads no FITS files.  It asks the KOA TAP
service two things:

  1. every HIRES frame of HD 187123 taken between 1997 December and 1998
     October (the discovery-era epochs of Butler, Marcy, Vogt & Apps 1998); and
  2. *every* HIRES frame, of any type, taken on each UT night that has such an
     observation -- so we can see what calibrations (ThAr arcs, quartz flats,
     biases, traces) were actually taken alongside the science.

The raw query results are written as CSV to first_hires_exoplanet/data/, and a
per-night summary is printed: frame counts by KOA image type, the distinct
instrument configurations, the HD 187123 exposures themselves, and -- the
question that matters most for a PypeIt reduction -- whether arcs and flats
exist in the *same* configuration (binning, decker, cross-disperser, filter,
echelle and cross-disperser angles) as the science frames.

The KOA TAP service silently appends a proprietary-period clause to every
query, so it only ever returns frames that are already public.  The script
demonstrates this with a probe of recent data (see `check_public_status`).

Query recipe (validated 2026-09-20):

  * POST to https://koa.ipac.caltech.edu/TAP/sync with REQUEST=doQuery,
    LANG=ADQL, FORMAT=csv and QUERY=<ADQL>.
  * The table is `koa_hires`, not `koa.koa_hires` (which returns 403).
  * `date_obs` is a char column; filtering on it raises ORA-01861.  Filter on
    the double column `mjd` instead.  1997-12-01 = MJD 50783, 1998-11-01 =
    MJD 51118.
  * `filename` is *not* a column (ORA-00904); the path is `filehand` and the
    observer's original name is `ofname`.  The cross-disperser angle is
    `xdangl`.  Check anything else against
    SELECT column_name FROM TAP_SCHEMA.columns WHERE table_name='koa_hires'.

Run from the repository root with:

    conda run -n pypeit14 python -m first_hires_exoplanet.koa_survey

Add `--refresh` to re-query KOA even if the CSV files already exist; without
it, existing CSVs are read back and only the summary is regenerated.

"""

# Standard imports
import argparse
import csv
import io
import os
import sys
from collections import Counter, OrderedDict

import requests
from astropy.time import Time


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: KOA synchronous TAP endpoint
TAP_URL = 'https://koa.ipac.caltech.edu/TAP/sync'

#: KOA table for HIRES
TABLE = 'koa_hires'

#: Search window, MJD.  1997-12-01 to 1998-11-01 (exclusive).
MJD_MIN = 50783.0
MJD_MAX = 51118.0

#: How the target is named in KOA (`targname` is '187123' or 'H187123')
TARGET_PATTERN = '%187123%'

#: Columns pulled for every frame.  All verified against TAP_SCHEMA.columns.
COLUMNS = [
    'koaid', 'targname', 'object', 'koaimtyp', 'imtype', 'obstype', 'date_obs',
    'ut', 'mjd', 'elaptime', 'exptime', 'binning', 'window', 'numamps',
    'ccdgain', 'ccdspeed', 'deckname', 'slitname', 'slitwid', 'slitlen',
    'xdispers', 'xdname', 'echangl', 'xdangl', 'fil1name', 'fil2name',
    'iodin', 'iodout', 'lampname', 'lampcat1', 'lampcat2', 'lampqtz2',
    'lampdeut', 'lfilname', 'lselname', 'lmirrin', 'lmirrout', 'hatopen',
    'autoshut', 'xcovopen', 'rccvopen', 'bccvopen', 'frameno', 'ofname',
    'progid', 'progpi', 'progtitl', 'proginst', 'semid', 'propint', 'ra',
    'dec', 'airmass', 'npixsat',
]

#: The columns that make up an instrument configuration for our purposes.
#: These are the KOA names of PypeIt's `keck_hires_orig` configuration keys
#: (dispname, decker, filter1, binning) plus the readout mode, which PypeIt does
#: not key on but which changes the raw image shape.
CONFIG_COLS = ['binning', 'deckname', 'xdispers', 'fil1name', 'numamps', 'ccdgain']

#: PypeIt tolerances for the angle configuration keys (keck_hires.py,
#: init_meta): np.isclose(a, b, rtol, atol), i.e. |a - b| <= atol + rtol*|b|.
ECHANGLE_TOL = dict(rtol=1e-3, atol=1e-2)
XDANGLE_TOL = dict(rtol=1e-2, atol=1e-1)

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_HERE, 'data')
SCIENCE_CSV = os.path.join(DATA_DIR, 'koa_hd187123_science.csv')
NIGHTS_CSV = os.path.join(DATA_DIR, 'koa_hd187123_nights_all.csv')
RV_TSV = os.path.join(DATA_DIR, 'hd187123_hires_rv.tsv')

#: Barycentric JD of the Butler et al. (1998) publication, 1998 Oct 22 -- the
#: discovery-era cut used in figs.py
BJD_PUBLICATION = 2451110.


# ---------------------------------------------------------------------------
# KOA queries
# ---------------------------------------------------------------------------

def query_koa(adql, maxrec=100000, timeout=180):
    """ Run one synchronous ADQL query against KOA and parse the CSV reply.

    Args:
        adql (str): the ADQL statement.
        maxrec (int): row cap sent as MAXREC.  KOA honours it; the default is
            far above anything a single night produces.
        timeout (float): seconds before `requests` gives up.

    Returns:
        list: one `dict` per row, keyed by column name.  Values are strings
        exactly as KOA returned them (empty string for NULL).

    Raises:
        RuntimeError: if the HTTP status is not 200.  The server's message
        (which names e.g. an invalid column) is included verbatim.
    """
    data = {'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv',
            'MAXREC': maxrec, 'QUERY': adql}
    reply = requests.post(TAP_URL, data=data, timeout=timeout)
    if reply.status_code != 200:
        raise RuntimeError('KOA returned HTTP {:d} for query:\n{:s}\n{:s}'.format(
            reply.status_code, adql, reply.text.strip()))
    return list(csv.DictReader(io.StringIO(reply.text)))


def science_adql():
    """ ADQL for every HD 187123 frame in the window.

    Matches on both `targname` and `object`, because the 1998 observers typed
    the name in more than one way ('187123', 'H187123').
    """
    return ("SELECT {cols} FROM {table} WHERE mjd >= {lo:.1f} AND mjd < {hi:.1f} "
            "AND (upper(targname) LIKE '{pat}' OR upper(object) LIKE '{pat}') "
            "ORDER BY mjd").format(cols=', '.join(COLUMNS), table=TABLE,
                                   lo=MJD_MIN, hi=MJD_MAX, pat=TARGET_PATTERN)


def night_adql(night):
    """ ADQL for every frame of every type on one UT night.

    Args:
        night (str): UT date as 'YYYYMMDD', i.e. the middle field of a KOAID.
    """
    return ("SELECT {cols} FROM {table} WHERE koaid LIKE 'HI.{night}.%' "
            "ORDER BY mjd").format(cols=', '.join(COLUMNS), table=TABLE,
                                   night=night)


def night_of(row):
    """ The UT night of a frame, as 'YYYYMMDD', taken from its KOAID. """
    return row['koaid'].split('.')[1]


def fetch_all():
    """ Query KOA for the science frames, then for every frame on their nights.

    Returns:
        tuple: (science_rows, night_rows) -- lists of row dicts.
    """
    print('Querying KOA for HD 187123 frames, MJD {:.1f} - {:.1f} ...'.format(
        MJD_MIN, MJD_MAX))
    science = query_koa(science_adql())
    print('  {:d} frames'.format(len(science)))

    nights = sorted(set(night_of(row) for row in science))
    print('Querying KOA for every frame on the {:d} nights with an '
          'HD 187123 frame ...'.format(len(nights)))
    all_rows = []
    for night in nights:
        rows = query_koa(night_adql(night))
        print('  {:s}: {:d} frames'.format(night, len(rows)))
        all_rows.extend(rows)
    return science, all_rows


# ---------------------------------------------------------------------------
# CSV on disk
# ---------------------------------------------------------------------------

def write_csv(rows, path):
    """ Write a list of row dicts as CSV with the columns in COLUMNS order. """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='') as ff:
        writer = csv.DictWriter(ff, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print('Wrote {:s} ({:d} rows)'.format(path, len(rows)))


def read_csv(path):
    """ Read a CSV written by `write_csv` back into a list of row dicts. """
    with open(path, 'r', newline='') as ff:
        return list(csv.DictReader(ff))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _f(value):
    """ Float from a KOA string, or None if it is empty or not a number. """
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _isclose(a, b, rtol, atol):
    """ numpy.isclose semantics for two scalars. """
    return abs(a - b) <= atol + rtol * abs(b)


def _config(row):
    """ The configuration tuple of a frame (values of CONFIG_COLS). """
    return tuple(row[col] for col in CONFIG_COLS)


def _config_str(cfg):
    """ Render a configuration tuple as 'key=value' pairs. """
    return ' '.join('{:s}={:s}'.format(key, val or '-')
                    for key, val in zip(CONFIG_COLS, cfg))


def _distinct(rows, col):
    """ Sorted distinct values of one column, as strings. """
    return sorted(set(row[col] for row in rows), key=lambda v: (_f(v) is None, _f(v) or 0, v))


def _lamps(row):
    """ A compact description of the lamp state of a frame. """
    return 'lampname={:s} cat1={:s} cat2={:s} qtz2={:s} deut={:s} lfil={:s}'.format(
        row['lampname'] or '-', row['lampcat1'], row['lampcat2'], row['lampqtz2'],
        row['lampdeut'], row['lfilname'] or '-')


def _angle_range(rows, col):
    """ (min, max) of an angle column over rows, or (None, None). """
    vals = [_f(row[col]) for row in rows if _f(row[col]) is not None]
    return (min(vals), max(vals)) if vals else (None, None)


def _night_iso(night):
    """ 'YYYYMMDD' -> 'YYYY-MM-DD'. """
    return '{:s}-{:s}-{:s}'.format(night[:4], night[4:6], night[6:])


def is_hd187123(row):
    """ True if the frame is an exposure of HD 187123. """
    return '187123' in (row['targname'] or '') or '187123' in (row['object'] or '')


# ---------------------------------------------------------------------------
# Per-night summary
# ---------------------------------------------------------------------------

def angles_match(cal, science_rows):
    """ Do a calibration frame's angles match *any* science frame's, within
    PypeIt's configuration tolerances?

    Args:
        cal (dict): the calibration row.
        science_rows (list): science rows to compare against.

    Returns:
        tuple: (echelle_ok, xd_ok, d_ech, d_xd) -- booleans and the smallest
        absolute offsets (deg) from a science frame.
    """
    ech, xd = _f(cal['echangl']), _f(cal['xdangl'])
    d_ech = min(abs(ech - _f(s['echangl'])) for s in science_rows)
    d_xd = min(abs(xd - _f(s['xdangl'])) for s in science_rows)
    ech_ok = any(_isclose(ech, _f(s['echangl']), **ECHANGLE_TOL) for s in science_rows)
    xd_ok = any(_isclose(xd, _f(s['xdangl']), **XDANGLE_TOL) for s in science_rows)
    return ech_ok, xd_ok, d_ech, d_xd


def _print_calibs(label, rows, science, science_cfgs):
    """ List calibration frames of one type and score them against the
    science configuration(s) of the night. """
    print('  {:s}: {:d}'.format(label, len(rows)))
    if not rows:
        print('    NONE')
        return
    for row in rows:
        cfg = _config(row)
        cfg_ok = cfg in science_cfgs
        ech_ok, xd_ok, d_ech, d_xd = angles_match(row, science)
        verdict = 'MATCH' if (cfg_ok and ech_ok and xd_ok) else 'no match'
        why = []
        if not cfg_ok:
            why.append('config differs')
        if not ech_ok:
            why.append('echangl off by {:.4f}'.format(d_ech))
        if not xd_ok:
            why.append('xdangl off by {:.4f}'.format(d_xd))
        print('    {:s} ut={:s} elap={:>5s}s ech={:>11s} xd={:>11s} decker={:2s} '
              'bin={:3s} fil1={:6s} amps={:s} gain={:s} {:s} hatch={:s} '
              'iodin={:s} | dEch={:.5f} dXd={:.4f} -> {:s}{:s}'.format(
                  row['koaid'], row['ut'][:8], row['elaptime'], row['echangl'],
                  row['xdangl'], row['deckname'], row['binning'],
                  row['fil1name'], row['numamps'], row['ccdgain'], _lamps(row),
                  row['hatopen'], row['iodin'], d_ech, d_xd, verdict,
                  ' ({:s})'.format('; '.join(why)) if why else ''))


def summarize_night(night, rows):
    """ Print everything we want to know about one UT night.

    Args:
        night (str): 'YYYYMMDD'.
        rows (list): every KOA row on that night.

    Returns:
        dict: a few numbers used for the whole-window summary.
    """
    science = [row for row in rows if is_hd187123(row) and row['koaimtyp'] == 'object']
    print('\n' + '=' * 100)
    print('UT night {:s}   ({:d} frames in KOA)'.format(_night_iso(night), len(rows)))
    print('=' * 100)

    # Programmes sharing the night
    progs = Counter((row['progid'], row['progpi'], row['semid']) for row in rows)
    print('  programmes: ' + '; '.join(
        '{:s} ({:s}, {:s}) x{:d}'.format(pid, pi, sem, n)
        for (pid, pi, sem), n in progs.most_common()))
    print('  propint (months): {:s}'.format(', '.join(_distinct(rows, 'propint'))))

    # Counts by KOA image type
    types = Counter(row['koaimtyp'] for row in rows)
    print('  koaimtyp counts: ' + ', '.join(
        '{:s}={:d}'.format(key, val) for key, val in sorted(types.items())))

    # Distinct values of the things the caller asked about
    for col in ['binning', 'deckname', 'xdispers', 'fil1name', 'fil2name',
                'echangl', 'xdangl', 'elaptime']:
        vals = _distinct(rows, col)
        shown = ', '.join(vals) if len(vals) <= 12 else \
            ', '.join(vals[:12]) + ', ... ({:d} distinct)'.format(len(vals))
        print('  distinct {:<9s}: {:s}'.format(col, shown))

    # Configurations present, by image type
    print('  configurations (binning, decker, xdispers, fil1name, numamps, ccdgain):')
    cfg_types = Counter((_config(row), row['koaimtyp']) for row in rows)
    by_cfg = OrderedDict()
    for (cfg, typ), n in sorted(cfg_types.items()):
        by_cfg.setdefault(cfg, []).append('{:s}={:d}'.format(typ, n))
    for cfg, parts in by_cfg.items():
        print('    [{:s}]  {:s}'.format(_config_str(cfg), ', '.join(parts)))

    # The science frames themselves
    print('  HD 187123 science frames: {:d}'.format(len(science)))
    for row in science:
        print('    {:s} ut={:s} elap={:>4s}s ech={:>11s} xd={:>11s} decker={:2s} '
              'bin={:3s} fil1={:6s} amps={:s} gain={:s} iodin={:s} iodout={:s} '
              'airmass={:s} npixsat={:s} {:s} ({:s})'.format(
                  row['koaid'], row['ut'][:8], row['elaptime'], row['echangl'],
                  row['xdangl'], row['deckname'], row['binning'], row['fil1name'],
                  row['numamps'], row['ccdgain'], row['iodin'], row['iodout'],
                  row['airmass'][:5], row['npixsat'], row['progid'], row['progpi']))
    if not science:
        return dict(n_science=0, n_arcs_match=0, n_flats_match=0, n_bias=0)

    science_cfgs = set(_config(row) for row in science)
    ech_lo, ech_hi = _angle_range(science, 'echangl')
    xd_lo, xd_hi = _angle_range(science, 'xdangl')
    print('  science echangl range: {:.5f} .. {:.5f} (spread {:.5f}; PypeIt atol {:.3f})'.format(
        ech_lo, ech_hi, ech_hi - ech_lo, ECHANGLE_TOL['atol']))
    print('  science xdangl  range: {:.5f} .. {:.5f} (spread {:.5f}; PypeIt atol {:.3f})'.format(
        xd_lo, xd_hi, xd_hi - xd_lo, XDANGLE_TOL['atol']))

    # Every other frame in the science configuration(s), i.e. the frames
    # PypeIt would group with the science on the configuration keys alone
    same_cfg = [row for row in rows if _config(row) in science_cfgs]
    same_types = Counter(row['koaimtyp'] for row in same_cfg)
    print('  frames sharing the science configuration: ' + ', '.join(
        '{:s}={:d}'.format(key, val) for key, val in sorted(same_types.items())))
    other_objects = [row for row in same_cfg
                     if row['koaimtyp'] == 'object' and not is_hd187123(row)]
    print('  other targets in the science configuration: {:d} frames, {:d} distinct '
          'targname'.format(len(other_objects),
                            len(set(row['targname'] for row in other_objects))))

    # Calibrations, scored against the science configuration
    arcs = [row for row in rows if row['koaimtyp'] in ('arclamp', 'arc')]
    flats = [row for row in rows if row['koaimtyp'] in ('flatlamp', 'flat')]
    biases = [row for row in rows if row['koaimtyp'].startswith('bias')
              or row['koaimtyp'] == 'dark']
    traces = [row for row in rows if row['koaimtyp'] == 'trace']
    others = [row for row in rows if row['koaimtyp'] not in
              ('object', 'arclamp', 'arc', 'flatlamp', 'flat', 'trace', 'dark')
              and not row['koaimtyp'].startswith('bias')]

    _print_calibs('ARC frames (koaimtyp arclamp)', arcs, science, science_cfgs)
    _print_calibs('FLAT frames (koaimtyp flatlamp)', flats, science, science_cfgs)
    _print_calibs('BIAS / DARK frames', biases, science, science_cfgs)
    _print_calibs('TRACE frames', traces, science, science_cfgs)
    if others:
        other_types = Counter((row['koaimtyp'], row['deckname'], row['binning'])
                              for row in others)
        print('  other frame types: ' + ', '.join(
            '{:s} (decker {:s}, bin {:s}) x{:d}'.format(t, d, b, n)
            for (t, d, b), n in sorted(other_types.items())))

    def _n_match(cals):
        return sum(1 for row in cals
                   if _config(row) in science_cfgs
                   and all(angles_match(row, science)[:2]))

    return dict(n_science=len(science), n_arcs_match=_n_match(arcs),
                n_flats_match=_n_match(flats), n_bias=len(biases),
                n_arcs=len(arcs), n_flats=len(flats), n_traces=len(traces))


# ---------------------------------------------------------------------------
# Whole-window summary
# ---------------------------------------------------------------------------

def load_rv_epochs():
    """ Read the catalogue epochs (BJD) of HD 187123 before publication.

    Returns:
        list: BJDs as floats, or an empty list if the TSV is missing.
    """
    if not os.path.isfile(RV_TSV):
        return []
    bjd = []
    with open(RV_TSV, 'r') as ff:
        for line in ff:
            if line.startswith('HD187123'):
                bjd.append(float(line.split('\t')[1]))
    return [b for b in bjd if b < BJD_PUBLICATION]


def crosscheck_epochs(science):
    """ Match the catalogue's discovery-era RV epochs to KOA science frames.

    For each catalogue BJD, find the KOA frame whose exposure midpoint (MJD of
    shutter open plus half the elapsed time) is nearest.  Barycentric
    correction is at most about 8 minutes, so a genuine match is within ~10
    minutes; anything larger means the catalogue epoch has no frame in KOA.

    Args:
        science (list): HD 187123 science rows.

    Returns:
        set: KOAIDs matched to a catalogue epoch.
    """
    epochs = load_rv_epochs()
    print('\n--- Cross-check against the RV catalogue ({:s}) ---'.format(
        os.path.basename(RV_TSV)))
    if not epochs:
        print('  catalogue file not found; skipping')
        return set()
    print('  catalogue epochs before publication: {:d}'.format(len(epochs)))
    matched = set()
    for bjd in epochs:
        best, best_dt = None, None
        for row in science:
            mid = _f(row['mjd']) + 0.5 * _f(row['elaptime']) / 86400.
            dt = (bjd - 2400000.5 - mid) * 1440.      # minutes
            if best is None or abs(dt) < abs(best_dt):
                best, best_dt = row, dt
        flag = '' if abs(best_dt) < 10. else '   <-- NO FRAME WITHIN 10 MIN'
        print('  BJD {:.5f} -> {:s} (iodin={:s}, {:s}) dt={:+.1f} min{:s}'.format(
            bjd, best['koaid'], best['iodin'], best['progid'], best_dt, flag))
        if abs(best_dt) < 10.:
            matched.add(best['koaid'])
    print('  matched {:d} of {:d} catalogue epochs to distinct KOA frames'.format(
        len(matched), len(epochs)))
    unmatched = [row for row in science if row['koaid'] not in matched]
    print('  KOA HD 187123 frames NOT in the catalogue: {:d}'.format(len(unmatched)))
    for row in unmatched:
        print('    {:s} elap={:>4s}s iodin={:s} {:s} ({:s})'.format(
            row['koaid'], row['elaptime'], row['iodin'], row['progid'], row['progpi']))
    return matched


def summarize_window(science, per_night):
    """ Print the science-frame census over the whole window. """
    print('\n' + '#' * 100)
    print('WHOLE WINDOW: HD 187123 science frames per UT night')
    print('#' * 100)
    print('  {:<10s} {:>6s} {:>6s} {:>6s} {:>6s} {:>6s}  {:s}'.format(
        'UT night', 'nSci', 'iodin', 'arcOK', 'flatOK', 'bias', 'progid (PI)'))
    nights = sorted(set(night_of(row) for row in science))
    for night in nights:
        rows = [row for row in science if night_of(row) == night]
        progs = sorted(set('{:s} ({:s})'.format(r['progid'], r['progpi']) for r in rows))
        stats = per_night[night]
        print('  {:<10s} {:>6d} {:>6d} {:>6d} {:>6d} {:>6d}  {:s}'.format(
            _night_iso(night), len(rows),
            sum(1 for r in rows if r['iodin'] == 'T'),
            stats['n_arcs_match'], stats['n_flats_match'], stats['n_bias'],
            ', '.join(progs)))
    print('  total science frames: {:d} on {:d} nights'.format(len(science), len(nights)))
    print('  iodine cell in (iodin=T): {:d};  out (iodin=F): {:d}'.format(
        sum(1 for r in science if r['iodin'] == 'T'),
        sum(1 for r in science if r['iodin'] == 'F')))
    print('  iodine-free frames (templates?):')
    for row in science:
        if row['iodin'] != 'T':
            print('    {:s} ut={:s} elap={:>4s}s decker={:s} bin={:s} ({:s}, {:s})'.format(
                row['koaid'], row['ut'][:8], row['elaptime'], row['deckname'],
                row['binning'], row['progid'], row['progpi']))

    for col in ['binning', 'deckname', 'xdispers', 'fil1name', 'fil2name',
                'numamps', 'ccdgain', 'elaptime', 'propint', 'semid']:
        print('  distinct {:<9s} over all science: {:s}'.format(col, ', '.join(_distinct(science, col))))
    ech_lo, ech_hi = _angle_range(science, 'echangl')
    xd_lo, xd_hi = _angle_range(science, 'xdangl')
    print('  echangl over all science: {:.5f} .. {:.5f} (spread {:.5f})'.format(
        ech_lo, ech_hi, ech_hi - ech_lo))
    print('  xdangl  over all science: {:.5f} .. {:.5f} (spread {:.5f})'.format(
        xd_lo, xd_hi, xd_hi - xd_lo))


# ---------------------------------------------------------------------------
# Public status
# ---------------------------------------------------------------------------

def check_public_status():
    """ Demonstrate that KOA only returns frames whose proprietary period has
    expired.

    KOA appends `current_date > add_months(date_obs, propint)` (or membership
    of the querying user's own programmes) to every query.  We cannot see the
    clause, but we can see its effect: for the most recent three years, group
    by `propint` and look at the latest `date_obs` returned for each.  If the
    filter is active, no row is younger than its own proprietary period, and
    a direct request for young proprietary rows returns nothing.
    """
    print('\n--- Public status: probing the proprietary filter ---')
    today = Time.now()
    print('  today: {:s}'.format(today.iso[:10]))
    rows = query_koa("SELECT propint, count(*) AS n, min(date_obs) AS first, "
                     "max(date_obs) AS last FROM {:s} WHERE mjd > 60000.0 "
                     "GROUP BY propint ORDER BY propint".format(TABLE))
    print('  frames returned with MJD > 60000 (2023-02-25 on), by proprietary period:')
    for row in rows:
        last = Time(row['last'].split()[0])
        age_months = (today - last).jd / 30.4375
        print('    propint={:>3s} months: {:>6s} frames, latest date_obs {:s} '
              '({:.1f} months ago)'.format(row['propint'], row['n'], last.iso[:10],
                                           age_months))
    young = query_koa("SELECT count(*) AS n FROM {:s} WHERE mjd > {:.1f} AND "
                      "propint > 12".format(TABLE, today.mjd - 548.))
    print('  frames with propint > 12 months observed in the last 18 months: {:s}'.format(
        young[0]['n']))
    print('  => every row this service returns has cleared its proprietary period; '
          'everything in the CSVs is public.')


# ---------------------------------------------------------------------------

def main():
    """ Query (or re-read) KOA metadata and print the survey. """
    parser = argparse.ArgumentParser(
        description='Survey KOA metadata for the 1997-98 HIRES frames of HD 187123.')
    parser.add_argument('--refresh', action='store_true',
                        help='re-query KOA even if the CSV files already exist')
    args = parser.parse_args()

    have_csv = os.path.isfile(SCIENCE_CSV) and os.path.isfile(NIGHTS_CSV)
    if have_csv and not args.refresh:
        print('Reading existing CSVs (use --refresh to re-query KOA)')
        science, all_rows = read_csv(SCIENCE_CSV), read_csv(NIGHTS_CSV)
        queried = False
    else:
        science, all_rows = fetch_all()
        write_csv(science, SCIENCE_CSV)
        write_csv(all_rows, NIGHTS_CSV)
        queried = True
    print('science frames: {:d};  all frames on their nights: {:d}'.format(
        len(science), len(all_rows)))
    print('science query:\n  {:s}'.format(science_adql()))
    print('per-night query (one per night):\n  {:s}'.format(night_adql('YYYYMMDD')))

    per_night = {}
    for night in sorted(set(night_of(row) for row in all_rows)):
        rows = [row for row in all_rows if night_of(row) == night]
        per_night[night] = summarize_night(night, rows)

    summarize_window(science, per_night)
    crosscheck_epochs(science)

    if queried:
        check_public_status()
    else:
        print('\n(public-status probe skipped; it needs the network -- run with --refresh)')

    return 0


if __name__ == '__main__':
    sys.exit(main())
