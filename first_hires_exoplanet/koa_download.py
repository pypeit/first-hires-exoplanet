""" Download one night of the 1998 July HIRES run on HD 187123 from KOA.

Fetches the raw (level-0) FITS frames for a PypeIt reduction of the 1998
July 16 (UT) observation of HD 187123, together with the calibration frames
that the prompt-3 survey and prompt-4 source audit (claude_prompts/
data_phase1_prompt.md) identified as the ones we need:

  1998-07-16  the single HD 187123 science frame (400 s, iodine cell in);
              both B1 ThAr arcs bracketing the night; all 16 B2 2 s quartz
              flats -- 8 with the hatch open and 8 with it closed, so the two
              can be compared; and both B1 3 s iodine flats.
  1998-07-14  the four iodine-free, hatch-closed, decker-B1 3 s "Narrowflat"
              frames that are the only clean flats in the science slit within
              PypeIt's angle tolerances of the whole July 15-19 run; the two
              hatch-closed B1 ThAr arcs nearest in time to those flats; and
              three 0 s biases, for inspection.

Every frame in these selections is in the Marcy planet-search configuration
(binning 1,2; one amplifier; ccdgain F; cross-disperser RED; filter kv370).
The daytime Mercury frames of T. Bida's program N01H, which share these
nights in a different setup (binning 1,1, two amplifiers, filter bg24a), are
excluded.  The one exception is 1998-07-16's only "bias" (`koaimtyp =
bias_lamp_on`), which is a Bida frame: it is fetched, for inspection only,
into a separate `inspect_only/` directory that no reduction should read.

The frame selection is expressed as explicit ADQL fragments in `SELECTIONS`
below, so the list is derived from the archive rather than hand-transcribed.
Where the selection is "the N nearest in time to the flats", the query
returns the candidates and the choice is made here, in `nearest_in_time`.

The raw frames are large and must **never** land inside the git repository.
The default output root is a sibling directory of the repository, and the
script refuses to write anywhere under the repository.

Retrieval recipe (validated 2026-09-20):

  * Metadata: POST to https://koa.ipac.caltech.edu/TAP/sync with
    REQUEST=doQuery, LANG=ADQL, FORMAT=csv, QUERY=<ADQL>, table `koa_hires`.
    The archive path of a frame is the `filehand` column (not `filename`,
    which does not exist).  `filesize_mb` exists but is NULL for these 1998
    frames, so sizes are checked from the download itself.
  * Data: GET https://koa.ipac.caltech.edu/cgi-bin/getKOA/nph-getKOA?filehand=
    <filehand>.  A good reply has Content-Type `image/x-fits`; a *wrong*
    filehand returns an HTML error page with HTTP status 200, so every
    download is verified: Content-Type, size against Content-Length, the
    file starts with `SIMPLE`, `astropy.io.fits` opens it, and the size on
    disk equals the size implied by the FITS header.  A failed download is
    deleted, never left in place.

Idempotent: a file that is already present and passes the FITS verification
is skipped, so a re-run costs only the metadata queries.

Run from the repository root with:

    conda run -n pypeit14 python -m first_hires_exoplanet.koa_download

Options: `--night 19980716` (repeatable) restricts to one night;
`--outroot DIR` changes the output root; `--dry-run` queries KOA and prints
the manifest without downloading anything.

"""

# Standard imports
import argparse
import csv
import io
import os
import sys

import requests
from astropy.io import fits


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: KOA synchronous TAP endpoint (metadata)
TAP_URL = 'https://koa.ipac.caltech.edu/TAP/sync'

#: KOA file retrieval endpoint (data)
GET_URL = 'https://koa.ipac.caltech.edu/cgi-bin/getKOA/nph-getKOA'

#: KOA table for HIRES
TABLE = 'koa_hires'

#: Columns pulled for every selected frame.  `filehand` is the archive path
#: that the retrieval endpoint needs; the rest describe the frame for the
#: manifest.
COLUMNS = ['koaid', 'filehand', 'koaimtyp', 'ut', 'mjd', 'elaptime', 'deckname',
           'binning', 'numamps', 'hatopen', 'iodin', 'lampname', 'echangl',
           'xdangl', 'object', 'targname', 'progid', 'ofname']

#: The repository root (parent of this package).  Nothing is written below it.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: Default output root: a sibling of the repository, outside git.
DEFAULT_OUTROOT = os.path.join(os.path.dirname(REPO_ROOT),
                               'first-hires-exoplanet-data', 'raw')

#: UT night ('YYYYMMDD', the middle field of a KOAID) -> output sub-directory.
#: 07-16 was reduced first (prompts 5-9); 07-15/17/18/19 are the rest of the
#: five-night run, added for prompt 10.  07-14 is not part of the run: it is
#: D. Latham's night, and supplies the only clean B1 flats in the science
#: configuration (see the prompt-4 Report).
NIGHT_DIRS = {'19980714': '1998jul14', '19980715': '1998jul15',
              '19980716': '1998jul16', '19980717': '1998jul17',
              '19980718': '1998jul18', '19980719': '1998jul19'}

#: The remaining nights of the run.  Each needs only its own science frames and
#: whatever B1 ThAr arcs it has; the flats come from 07-14, which is the only
#: source of clean B1 flats.  1998-07-19 has no calibrations of any kind and
#: borrows the 07-18 arc (see the prompt-3 Report).
REST_OF_RUN = ('19980715', '19980717', '19980718', '19980719')

#: Sub-directory, under the output root, for frames we only want to look at
INSPECT_DIR = 'inspect_only'

#: The Marcy planet-search configuration, shared by every frame we reduce
MARCY_CONFIG = "binning = '1,2' AND numamps = 1 AND xdispers = 'RED' AND fil1name = 'kv370'"

#: Frame selections.  Each entry is (label, night, ADQL WHERE fragment,
#: keep) where `keep` is None (take every row) or (n, 'nearest_to', label):
#: keep the n rows nearest in time to the rows of that other selection on the
#: same night.  The night restriction and MARCY_CONFIG are added by
#: `selection_adql`, except where `config` is given explicitly.
SELECTIONS = [
    # ---- 1998-07-16: the night we reduce -------------------------------
    dict(label='science', night='19980716',
         where="koaimtyp = 'object' AND (targname LIKE '%187123%' OR object LIKE '%187123%')",
         keep=None),
    dict(label='arc_B1', night='19980716',
         where="koaimtyp = 'arclamp' AND deckname = 'B1' AND lampname = 'ThAr1'",
         keep=None),
    dict(label='flat_B2_quartz', night='19980716',
         where="koaimtyp = 'flatlamp' AND deckname = 'B2' AND iodin = 'F' AND lampqtz2 = 'T'",
         keep=None),
    dict(label='flat_B1_iodine', night='19980716',
         where="koaimtyp = 'flatlamp' AND deckname = 'B1' AND iodin = 'T'",
         keep=None),
    # The only "bias" on 07-16 is T. Bida's 1,1 two-amplifier frame.  Inspect
    # only; it goes to INSPECT_DIR and is not in the Marcy configuration.
    dict(label='bias_lamp_on_BIDA', night='19980716',
         where="koaimtyp = 'bias_lamp_on'", config='1 = 1',
         subdir=INSPECT_DIR, keep=None),
    # ---- 1998-07-14: the clean B1 flats and their neighbours ---------------
    dict(label='flat_B1_clean', night='19980714',
         where="koaimtyp = 'flatlamp' AND deckname = 'B1' AND iodin = 'F' AND hatopen = 'F'",
         keep=None),
    dict(label='arc_B1', night='19980714',
         where="koaimtyp = 'arclamp' AND deckname = 'B1' AND lampname = 'ThAr1' AND hatopen = 'F'",
         keep=(2, 'nearest_to', 'flat_B1_clean')),
    dict(label='bias', night='19980714',
         where="koaimtyp = 'bias' AND elaptime = 0",
         keep=(3, 'nearest_to', 'flat_B1_clean')),
]

# ---- the remaining nights of the run: science plus any B1 ThAr arc ---------
for _night in REST_OF_RUN:
    SELECTIONS.append(dict(
        label='science', night=_night,
        where="koaimtyp = 'object' AND (targname LIKE '%187123%' "
              "OR object LIKE '%187123%')",
        keep=None))
    SELECTIONS.append(dict(
        label='arc_B1', night=_night,
        where="koaimtyp = 'arclamp' AND deckname = 'B1' AND lampname = 'ThAr1'",
        keep=None))

#: Chunk size for streaming downloads, bytes
CHUNK = 1 << 20

#: Seconds before `requests` gives up on one call
TIMEOUT = 300


# ---------------------------------------------------------------------------
# KOA metadata
# ---------------------------------------------------------------------------

def query_koa(adql, timeout=180):
    """ Run one synchronous ADQL query against KOA and parse the CSV reply.

    Args:
        adql (str): the ADQL statement.
        timeout (float): seconds before `requests` gives up.

    Returns:
        list: one `dict` per row, keyed by column name; values are strings.

    Raises:
        RuntimeError: if the HTTP status is not 200, with the server's
            message verbatim.
    """
    data = {'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': adql}
    reply = requests.post(TAP_URL, data=data, timeout=timeout)
    if reply.status_code != 200:
        raise RuntimeError('KOA returned HTTP {:d} for query:\n{:s}\n{:s}'.format(
            reply.status_code, adql, reply.text.strip()))
    return list(csv.DictReader(io.StringIO(reply.text)))


def selection_adql(sel):
    """ The full ADQL statement for one entry of SELECTIONS. """
    config = sel.get('config', MARCY_CONFIG)
    return ("SELECT {cols} FROM {table} WHERE koaid LIKE 'HI.{night}.%' "
            "AND {config} AND {where} ORDER BY mjd").format(
                cols=', '.join(COLUMNS), table=TABLE, night=sel['night'],
                config=config, where=sel['where'])


def nearest_in_time(candidates, anchors, n):
    """ The n candidate rows nearest in time (MJD) to any anchor row.

    Args:
        candidates (list): rows to choose from.
        anchors (list): rows to measure distance to.
        n (int): how many to keep.

    Returns:
        list: the chosen rows, in time order.
    """
    def _dist(row):
        return min(abs(float(row['mjd']) - float(a['mjd'])) for a in anchors)
    chosen = sorted(candidates, key=_dist)[:n]
    return sorted(chosen, key=lambda r: float(r['mjd']))


def select_frames(nights):
    """ Query KOA for every selection on the requested nights.

    Args:
        nights (list): UT nights as 'YYYYMMDD'.

    Returns:
        list: one dict per frame to fetch -- the KOA row plus 'label' and
        'outdir' (relative to the output root).
    """
    chosen_by_label = {}
    frames = []
    for sel in SELECTIONS:
        if sel['night'] not in nights:
            continue
        adql = selection_adql(sel)
        rows = query_koa(adql)
        key = (sel['night'], sel['label'])
        if sel['keep'] is not None:
            n, _, anchor_label = sel['keep']
            anchors = chosen_by_label[(sel['night'], anchor_label)]
            kept = nearest_in_time(rows, anchors, n)
        else:
            kept = rows
        chosen_by_label[key] = kept
        print('{:s} {:<18s} {:>3d} of {:>3d} rows kept  <- {:s}'.format(
            sel['night'], sel['label'], len(kept), len(rows), adql))
        subdir = sel.get('subdir', NIGHT_DIRS[sel['night']])
        for row in kept:
            row = dict(row)
            row['label'] = sel['label']
            row['outdir'] = subdir
            frames.append(row)
    return frames


# ---------------------------------------------------------------------------
# Download and verification
# ---------------------------------------------------------------------------

def fits_expected_size(path):
    """ The size a FITS file should have, from its headers.

    The end of the last HDU's data, rounded up to the 2880 byte FITS record.
    A truncated or padded file will not match.
    """
    with fits.open(path, memmap=False) as hdul:
        info = hdul[-1].fileinfo()
        end = info['datLoc'] + info['datSpan']
    return -(-end // 2880) * 2880


def verify_fits(path):
    """ Check that `path` is a complete, readable FITS file.

    Returns:
        tuple: (ok, message).  `ok` is True only if the file starts with
        'SIMPLE', opens with astropy, has a data array in its primary HDU,
        and has exactly the size its headers imply.
    """
    if not os.path.isfile(path):
        return False, 'missing'
    size = os.path.getsize(path)
    with open(path, 'rb') as ff:
        head = ff.read(6)
    if head != b'SIMPLE':
        return False, 'does not start with SIMPLE (first bytes {!r})'.format(head)
    try:
        # memmap=False: these 1998 frames are 16-bit with BZERO/BSCALE, which
        # astropy refuses to memory-map
        with fits.open(path, memmap=False) as hdul:
            shape = hdul[0].data.shape if hdul[0].data is not None else None
        expected = fits_expected_size(path)
    except Exception as err:  # astropy raises a variety of things for junk
        return False, 'astropy cannot open it: {!r}'.format(err)
    if shape is None:
        return False, 'primary HDU has no data'
    if size != expected:
        return False, 'size {:d} != {:d} implied by FITS headers'.format(size, expected)
    return True, 'FITS ok, shape {} , {:d} bytes'.format(shape, size)


def download(filehand, path):
    """ Stream one KOA frame to disk and verify it.

    Args:
        filehand (str): the KOA archive path (`filehand` column).
        path (str): destination file.

    Returns:
        tuple: (ok, message).  On failure nothing is left at `path`.
    """
    tmp = path + '.part'
    try:
        with requests.get(GET_URL, params={'filehand': filehand}, stream=True,
                          timeout=TIMEOUT) as reply:
            ctype = reply.headers.get('Content-Type', '')
            clen = reply.headers.get('Content-Length')
            if reply.status_code != 200:
                return False, 'HTTP {:d}'.format(reply.status_code)
            if 'fits' not in ctype.lower():
                snippet = reply.raw.read(300).decode('ascii', 'replace').strip()
                return False, 'Content-Type {!r}, not FITS; body starts: {!r}'.format(
                    ctype, snippet)
            with open(tmp, 'wb') as ff:
                for chunk in reply.iter_content(chunk_size=CHUNK):
                    ff.write(chunk)
        size = os.path.getsize(tmp)
        if clen is not None and int(clen) != size:
            return False, 'got {:d} bytes, Content-Length said {:s}'.format(size, clen)
        ok, msg = verify_fits(tmp)
        if not ok:
            return False, msg
        os.replace(tmp, path)
        return True, msg + ' (Content-Type {:s})'.format(ctype)
    except (requests.RequestException, OSError) as err:
        return False, repr(err)
    finally:
        if os.path.isfile(tmp):
            os.remove(tmp)


def check_outroot(outroot):
    """ Refuse an output root inside the git repository. """
    real_out = os.path.realpath(outroot)
    real_repo = os.path.realpath(REPO_ROOT)
    if real_out == real_repo or real_out.startswith(real_repo + os.sep):
        raise SystemExit('Refusing to download into the repository: {:s}'.format(outroot))


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def print_manifest(frames, outroot):
    """ Print one line per frame with its status and the total bytes. """
    print('\n' + '=' * 118)
    print('MANIFEST  (output root {:s})'.format(outroot))
    print('=' * 118)
    print('{:<22s} {:<18s} {:<12s} {:<8s} {:>6s} {:<3s} {:<5s} {:<5s} {:>9s}  {:s}'.format(
        'koaid', 'label', 'koaimtyp', 'UT', 'elap', 'dck', 'hatch', 'iodin', 'bytes', 'status'))
    total = 0
    n_ok = 0
    for row in frames:
        size = os.path.getsize(row['path']) if os.path.isfile(row['path']) else 0
        total += size
        n_ok += row['ok']
        print('{:<22s} {:<18s} {:<12s} {:<8s} {:>6s} {:<3s} {:<5s} {:<5s} {:>9d}  {:s}'.format(
            row['koaid'], row['label'], row['koaimtyp'], row['ut'][:8], row['elaptime'],
            row['deckname'], row['hatopen'], row['iodin'], size, row['status']))
    print('-' * 118)
    print('{:d} of {:d} frames verified on disk; {:d} bytes ({:.1f} MB) total'.format(
        n_ok, len(frames), total, total / 1e6))
    by_dir = {}
    for row in frames:
        by_dir.setdefault(os.path.dirname(row['path']), []).append(row)
    for outdir, rows in sorted(by_dir.items()):
        print('  {:s}: {:d} files'.format(outdir, len(rows)))


# ---------------------------------------------------------------------------

def main():
    """ Select frames from KOA, download the missing ones, print a manifest. """
    parser = argparse.ArgumentParser(
        description='Download the 1998-07-16 HIRES frames of HD 187123 and their '
                    'calibrations (plus the 1998-07-14 clean B1 flats) from KOA.')
    parser.add_argument('--night', action='append', choices=sorted(NIGHT_DIRS),
                        help='UT night to fetch (repeatable; default: all)')
    parser.add_argument('--outroot', default=DEFAULT_OUTROOT,
                        help='output root; frames go in <outroot>/<night dir>/ '
                             '(default: {:s})'.format(DEFAULT_OUTROOT))
    parser.add_argument('--dry-run', action='store_true',
                        help='query KOA and print the manifest; download nothing')
    args = parser.parse_args()

    nights = args.night or sorted(NIGHT_DIRS)
    check_outroot(args.outroot)

    print('Selecting frames from KOA ...')
    frames = select_frames(nights)
    print('{:d} frames selected'.format(len(frames)))

    n_fail = 0
    for row in frames:
        outdir = os.path.join(args.outroot, row['outdir'])
        row['path'] = os.path.join(outdir, row['koaid'])
        ok, msg = verify_fits(row['path'])
        if ok:
            row['ok'], row['status'] = True, 'present, ' + msg
            continue
        if args.dry_run:
            row['ok'], row['status'] = False, 'would download ' + row['filehand']
            continue
        os.makedirs(outdir, exist_ok=True)
        print('  fetching {:s} ({:s}) ...'.format(row['koaid'], row['label']), flush=True)
        ok, msg = download(row['filehand'], row['path'])
        if not ok:
            # One retry, then report loudly
            print('    FAILED ({:s}); retrying once'.format(msg), flush=True)
            ok, msg = download(row['filehand'], row['path'])
        row['ok'] = ok
        row['status'] = ('downloaded, ' if ok else 'FAILED: ') + msg
        if not ok:
            n_fail += 1
            print('    FAILED again: {:s}'.format(msg), flush=True)

    print_manifest(frames, args.outroot)
    if n_fail:
        print('\n*** {:d} download(s) FAILED -- see the manifest ***'.format(n_fail))
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
