""" Reduce the five nights of the 1998 July HIRES run on HD 187123.

Each night is reduced in isolation, the same way 1998-07-16 was in prompt 8, so
that night-to-night differences reflect the data rather than the recipe:

  * the night's own science frames and B1 ThAr arcs;
  * the four clean B1 `Narrowflat` frames from 1998-07-14, which are the only
    iodine-free, hatch-closed flats in the science configuration anywhere in the
    run (see the prompt-4 Report);
  * for 1998-07-19, which has no calibrations at all, the 07-18 arc.

Each night is staged as a directory of symlinks so `pypeit_setup` sees exactly
the frames we intend and nothing else.  The three reduction parameters
established in prompt 8 are then injected into the generated `.pypeit` file:

    find_trim_edge = 1,1   the B1 orders are only ~10 binned px wide, and the
                           inherited 3,3 masks six of them
    skip_skysub = True     HD 187123 fills the slit; the global sky fit has no
                           object-free pixels
    no_local_sky = True    ... and neither does the local one

Nothing here is written inside the git repository: the raw frames, the staging
directories and the reductions all live in the sibling data tree.

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.reduce_run
    conda run -n pypeit14 python -m first_hires_exoplanet.reduce_run --night 19980719
    conda run -n pypeit14 python -m first_hires_exoplanet.reduce_run --summarise-only

"""

# Standard imports
import os
import sys
import glob
import json
import shutil
import argparse
import subprocess

import numpy as np


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                         'first-hires-exoplanet-data')
RAW_ROOT = os.path.join(DATA_ROOT, 'raw')
REDUX_ROOT = os.path.join(DATA_ROOT, 'redux')

#: The night that supplies the clean B1 flats for every other night
FLAT_NIGHT = '1998jul14'

#: Basenames of the four clean B1 `Narrowflat` frames on FLAT_NIGHT.  The
#: other 07-14 frames (arcs, biases) are deliberately left out: each night
#: uses its own arcs, so the wavelength solution is same-night wherever one
#: exists.
FLAT_FRAMES = ('HI.19980714.07711.fits', 'HI.19980714.07810.fits',
               'HI.19980714.07910.fits', 'HI.19980714.55408.fits')

#: UT night -> (raw sub-directory, extra frames borrowed from elsewhere).
#: 1998-07-19 took no calibrations whatsoever; it borrows the 07-18 arc, which
#: is within PypeIt's echangle/xdangle tolerances (prompt 3 Report).
NIGHTS = {
    '19980715': ('1998jul15', ()),
    '19980716': ('1998jul16', ()),
    '19980717': ('1998jul17', ()),
    '19980718': ('1998jul18', ()),
    '19980719': ('1998jul19', (os.path.join('1998jul18', 'HI.19980718.54587.fits'),)),
}

#: The prompt-8 parameter block, injected before the `# Setup` section
PARAM_BLOCK = """[reduce]
    [[findobj]]
        # The B1 (3.5") orders are only ~10.3 binned pixels wide; the inherited
        # find_trim_edge = 3,3 masks six of them and object finding fails.
        find_trim_edge = 1,1
        # HD 187123 (V=7.9) fills the slit, so the global sky fit has no
        # object-free pixels to work with.
        skip_skysub = True
    [[skysub]]
        # ... and neither does the local sky model during extraction.
        no_local_sky = True
"""

#: Frames that are science rather than calibration, by KOAID prefix
SCIENCE_TARGET = '187123'


# ---------------------------------------------------------------------------
# Staging
# ---------------------------------------------------------------------------

def stage_night(night, overwrite=False):
    """ Build a symlink directory holding exactly the frames for one night.

    Args:
        night (str): UT night, 'YYYYMMDD'.
        overwrite (bool, optional): rebuild an existing staging directory.

    Returns:
        str: the staging directory.
    """
    subdir, borrowed = NIGHTS[night]
    stage = os.path.join(REDUX_ROOT, 'stage_{:s}'.format(night))
    if os.path.isdir(stage):
        if not overwrite:
            return stage
        shutil.rmtree(stage)
    os.makedirs(stage)

    sources = sorted(glob.glob(os.path.join(RAW_ROOT, subdir, 'HI.*.fits')))
    sources += [os.path.join(RAW_ROOT, f) for f in borrowed]
    sources += [os.path.join(RAW_ROOT, FLAT_NIGHT, f) for f in FLAT_FRAMES]

    for src in sources:
        if not os.path.isfile(src):
            raise RuntimeError('Missing raw frame: {:s}'.format(src))
        os.symlink(src, os.path.join(stage, os.path.basename(src)))
    print('  staged {:d} frames in {:s}'.format(len(sources), stage))
    return stage


# ---------------------------------------------------------------------------
# Setup and reduction
# ---------------------------------------------------------------------------

def _run(cmd, cwd=None, log=None):
    """ Run a command, tee-ing its output to `log`, and raise if it fails. """
    with open(log, 'w') as fh:
        proc = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError('{:s} failed ({:d}); see {:s}'.format(
            cmd[0], proc.returncode, log))


def science_setup(setup_dir):
    """ Return the generated `.pypeit` file whose data block holds science frames.

    `pypeit_setup` splits the B1 science configuration from the B2 flats, and
    the setup letters are not stable between runs, so the science setup has to
    be identified by content rather than by name.

    Args:
        setup_dir (str): directory `pypeit_setup` wrote into.

    Returns:
        str: path to the science `.pypeit` file.
    """
    hits = []
    for path in sorted(glob.glob(os.path.join(setup_dir, '*', '*.pypeit'))):
        with open(path) as fh:
            body = fh.read()
        n_sci = 0
        for line in body.splitlines():
            stripped = line.strip()
            if not stripped.startswith('HI.') or '|' not in stripped:
                continue
            # The columns are '|'-separated and whitespace-padded; the
            # frametype is the second field and may be a comma-separated list.
            ftypes = [t.strip() for t in stripped.split('|')[1].split(',')]
            if 'science' in ftypes:
                n_sci += 1
        if n_sci:
            hits.append((path, n_sci))
    if len(hits) != 1:
        raise RuntimeError('Expected one setup with science frames in {:s}, '
                           'found {:d}'.format(setup_dir, len(hits)))
    return hits[0][0]


def inject_params(pypeit_file):
    """ Add the prompt-8 reduction parameters to a generated `.pypeit` file. """
    with open(pypeit_file) as fh:
        text = fh.read()
    if 'find_trim_edge' in text:
        return
    marker = '\n# Setup\n'
    if marker not in text:
        raise RuntimeError('No "# Setup" section in {:s}'.format(pypeit_file))
    with open(pypeit_file, 'w') as fh:
        fh.write(text.replace(marker, '\n' + PARAM_BLOCK + marker, 1))


def reduce_night(night, overwrite=False):
    """ Stage, set up and reduce one night.

    Returns:
        str: the reduction directory.
    """
    print('=' * 72)
    print('NIGHT {:s}'.format(night))
    stage = stage_night(night, overwrite=overwrite)

    setup_dir = os.path.join(REDUX_ROOT, 'setup_{:s}'.format(night))
    redux_dir = os.path.join(REDUX_ROOT, 'reduce_{:s}'.format(night))
    for path in (setup_dir, redux_dir):
        if overwrite and os.path.isdir(path):
            shutil.rmtree(path)
    os.makedirs(setup_dir, exist_ok=True)
    os.makedirs(redux_dir, exist_ok=True)

    print('  pypeit_setup ...')
    _run(['pypeit_setup', '-s', 'keck_hires_orig', '-r', stage, '-c', 'all',
          '-d', setup_dir], log=os.path.join(setup_dir, 'pypeit_setup.log'))

    src = science_setup(setup_dir)
    dst = os.path.join(redux_dir, os.path.basename(src))
    shutil.copy(src, dst)
    inject_params(dst)

    # A reduction must describe exactly one pass: run_pypeit -o appends to an
    # existing spec1d instead of replacing it (see the prompt-9 log).
    science = os.path.join(redux_dir, 'Science')
    if os.path.isdir(science):
        shutil.rmtree(science)

    print('  run_pypeit ...')
    _run(['run_pypeit', os.path.basename(dst), '-o'], cwd=redux_dir,
         log=os.path.join(redux_dir, 'run_pypeit.log'))
    return redux_dir


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def summarise(night):
    """ Collect per-night reduction statistics.

    Returns:
        dict or None: summary, or None if the night has no reduction.
    """
    from pypeit import specobjs
    from pypeit.wavecalib import WaveCalib

    redux_dir = os.path.join(REDUX_ROOT, 'reduce_{:s}'.format(night))
    spec1ds = sorted(glob.glob(os.path.join(redux_dir, 'Science', 'spec1d_*.fits')))
    if not spec1ds:
        return None

    wcs = sorted(glob.glob(os.path.join(redux_dir, 'Calibrations', 'WaveCalib_*.fits')))
    rms = np.array([])
    if wcs:
        wc = WaveCalib.from_file(wcs[0])
        rms = np.array([w.rms for w in wc.wv_fits if w.rms is not None])

    frames = []
    for path in spec1ds:
        sobjs = specobjs.SpecObjs.from_fitsfile(path)
        orders, snr = [], []
        for obj in sobjs:
            flux, sig = obj.OPT_COUNTS, obj.OPT_COUNTS_SIG
            good = np.isfinite(flux) & np.isfinite(sig) & (sig > 0)
            if good.sum() < 200:
                continue
            orders.append(int(obj.ECH_ORDER))
            snr.append(float(np.median(flux[good] / sig[good])))
        if not orders:
            continue
        frames.append(dict(
            spec1d=os.path.basename(path),
            n_orders=len(orders),
            order_min=min(orders), order_max=max(orders),
            snr_min=min(snr), snr_med=float(np.median(snr)), snr_max=max(snr),
            duplicated=len(orders) != len(set(orders))))

    return dict(night=night, n_frames=len(frames), frames=frames,
                rms_min=float(rms.min()) if rms.size else None,
                rms_med=float(np.median(rms)) if rms.size else None,
                rms_max=float(rms.max()) if rms.size else None,
                n_solutions=int(rms.size))


def print_summary(summaries):
    """ Print the night-to-night comparison table. """
    print()
    print('=' * 96)
    print('NIGHT-TO-NIGHT STABILITY, 1998 July run')
    print('=' * 96)
    print('{:<10s} {:>6s} {:>7s} {:>9s} {:>8s} {:>8s} {:>8s}   {:s}'.format(
        'night', 'frames', 'orders', 'orders', 'S/N min', 'S/N med', 'S/N max',
        'wavecal RMS (px)'))
    for s in summaries:
        if s is None:
            continue
        for i, fr in enumerate(s['frames']):
            rms = ''
            if i == 0 and s['rms_med'] is not None:
                rms = '{:.3f} / {:.3f} / {:.3f}  (n={:d})'.format(
                    s['rms_min'], s['rms_med'], s['rms_max'], s['n_solutions'])
            print('{:<10s} {:>6s} {:>7d} {:>9s} {:>8.1f} {:>8.1f} {:>8.1f}   {:s}'.format(
                s['night'] if i == 0 else '', '{:d}/{:d}'.format(i + 1, s['n_frames']),
                fr['n_orders'], '{:d}-{:d}'.format(fr['order_min'], fr['order_max']),
                fr['snr_min'], fr['snr_med'], fr['snr_max'], rms))
    print('=' * 96)


# ---------------------------------------------------------------------------

def main():
    """ Reduce every night of the run and report stability. """
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--night', action='append', choices=sorted(NIGHTS),
                        help='UT night to reduce (repeatable; default: all)')
    parser.add_argument('--overwrite', action='store_true',
                        help='rebuild staging, setup and reduction directories')
    parser.add_argument('--summarise-only', action='store_true',
                        help='skip the reductions; just report what is on disk')
    args = parser.parse_args()

    nights = args.night or sorted(NIGHTS)

    if not args.summarise_only:
        for night in nights:
            reduce_night(night, overwrite=args.overwrite)

    summaries = [summarise(n) for n in sorted(NIGHTS)]
    print_summary(summaries)

    out = os.path.join(REDUX_ROOT, 'run_summary.json')
    with open(out, 'w') as fh:
        json.dump([s for s in summaries if s], fh, indent=2)
    print('Wrote {:s}'.format(out))


if __name__ == '__main__':
    main()
