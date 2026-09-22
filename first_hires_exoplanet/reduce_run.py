""" Reduce the discovery-era HIRES nights on HD 187123.

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

Phase 2 prompt 3 extends this over the whole nine-month discovery era.  The
flats are the reason that is not a simple loop.  Clean flats -- quartz, hatch
closed, iodine cell out -- exist on most nights, but almost always through the
**B2** decker: this program flat-fielded through the wider slit as a matter of
course.  Clean flats in the **B1** science decker exist on only five nights of
the entire era: 1998-07-14 (4), 1998-08-12 (18), 1998-08-17 (2), 1998-09-17
(10) and 1998-09-18 (8).  PypeIt will only accept a flat whose echelle angle is
within about 0.01 of the science frame's, and the era spans 0.021, so no single
flat set covers it.

The rule used here, applied uniformly:

  * a night uses **its own** clean flats of its own decker when it has at least
    four of them;
  * otherwise it borrows **1998-08-12**'s B1 flats, which are legal for every
    B1 night of 1998 and sit closest to the middle of the era;
  * a **borrowed** flat set is capped at four, the four nearest in time to the
    night's science frames, so a borrowing night's recipe does not depend on
    how many flats the donor happened to take.  A night using its own flats
    stages all of them: its whole raw directory is staged, and the download
    selection already restricts that directory to science frames, arcs and
    clean flats;
  * **1997-12-23 and 1997-12-24 have no legal donor at all.**  Their echelle
    angles (+0.0149, +0.0160) are more than 0.01 from every clean B1 flat in
    the era.  They are attempted anyway, and what happens is the answer to
    prompt 3's "which need intervention".

Unlike 1998-07-19, every other discovery-era night took ThAr arcs of its own,
so no night after July borrows a wavelength solution.

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.reduce_run
    conda run -n pypeit14 python -m first_hires_exoplanet.reduce_run --night 19980826
    conda run -n pypeit14 python -m first_hires_exoplanet.reduce_run --era
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
    # ---- the rest of the discovery era, phase 2 prompt 3 ------------------
    '19971223': ('1997dec23', ()),
    '19971224': ('1997dec24', ()),
    '19980618': ('1998jun18', ()),
    '19980812': ('1998aug12', ()),
    '19980817': ('1998aug17', ()),
    '19980818': ('1998aug18', ()),
    '19980825': ('1998aug25', ()),
    '19980826': ('1998aug26', ()),
    '19980912': ('1998sep12', ()),
    '19980913': ('1998sep13', ()),
    '19980914': ('1998sep14', ()),
    '19980915': ('1998sep15', ()),
    '19980916': ('1998sep16', ()),
    # 1998-09-17 is the only B2 night and took no B2 arc; it borrows the
    # 1998-09-15 B2 arc, which is within tolerance (see the prompt-3 Report).
    '19980917': ('1998sep17', (os.path.join('1998sep15', 'HI.19980915.60867.fits'),)),
    '19980918': ('1998sep18', ()),
}

#: The five nights of the 1998 July run, reduced in phase 1
JULY_RUN = ('19980715', '19980716', '19980717', '19980718', '19980719')

#: The nights added in phase 2 prompt 3
ERA_NIGHTS = tuple(n for n in sorted(NIGHTS) if n not in JULY_RUN)

#: The flat donor for every 1998 B1 night without four clean B1 flats of its own
DONOR_NIGHT = '19980812'

#: How many flats every reduction uses, so the recipe does not vary
N_FLATS = 4

#: Nights whose flats must be forced into the generated `.pypeit` file, with
#: the frametype to force.  Only December 1997 needs this.
#:
#: Two independent things stop PypeIt calibrating those nights on its own.  The
#: echelle angle (+0.0149, +0.0160) is outside the 0.01 grouping tolerance of
#: every clean B1 flat in the era; and, decisively, `dispname` is a hard config
#: key with no tolerance, and `keck_hires.py` overrides the header's
#: `XDISPERS = 'RED'` to `RED97` for everything before 1997-12-31, on the MAKEE
#: DRP's claim that a different cross-disperser was fitted.  No 1998 frame can
#: ever join a 1997 configuration.
#:
#: So December is calibrated entirely from its own night.  The only B1 quartz
#: flat it has is hatch-closed, cover-open and **iodine-in**, at exactly the
#: science echelle angle.  That is the right frame to trace with -- the cell
#: modulates the spectrum along dispersion without moving the orders -- and
#: `order_shift.py` confirms the geometry independently: the December arcs sit
#: 1 binned pixel from 1998-07-14's, at correlation 0.986 and 0.994, closer
#: than several nights PypeIt accepts without complaint.
#:
#: What an iodine-in flat must NOT do is serve as a pixel flat: its I2 forest
#: would be divided out of the science frames between 5000 and 6200 A, which is
#: exactly the signal phase 3 depends on.  `DECEMBER_PARAMS` therefore switches
#: pixel-flatting off for these two nights.  The illumination flat is a smooth
#: fit along the slit and the narrow I2 lines do not survive it, so that stays.
FORCED_CALIBS = {
    '19971223': {'HI.19971223.15114.fits': 'pixelflat,illumflat,trace'},
    '19971224': {'HI.19971224.14230.fits': 'pixelflat,illumflat,trace'},
}

#: Extra parameters for the December nights (see FORCED_CALIBS)
DECEMBER_PARAMS = """[scienceframe]
    [[process]]
        # The only B1 flat these nights have was taken with the iodine cell in.
        # It traces the orders correctly, but its I2 forest must never reach the
        # science frames: dividing an iodine-in flat into iodine-in science
        # would partly cancel the 5000-6200 A absorption that phase 3 measures
        # velocities from.  So the flat is used for tracing only.
        #
        # Both corrections have to go together.  PypeIt validates that a
        # slit-illumination or spectral flat-field correction is only applied
        # alongside the pixel flat (pypeitpar.py:575), so use_illumflat = True
        # with use_pixelflat = False is rejected outright.
        #
        # The cost is the pixel-to-pixel correction, which phase 1 measured at
        # +/-3%, and the slit-illumination correction.  Both are smooth and
        # multiplicative; neither moves a line centre, so neither biases a
        # velocity.  These two nights are nonetheless flat-fielded differently
        # from the other eighteen and must be treated as such.
        use_pixelflat = False
        use_illumflat = False
"""

#: A lit quartz flat is far above the bias level; an unlit one is not.  The
#: clean flats of this era sit at IM01MN01 = 6,700-20,400 and the unlit ones at
#: 763-768, which is the bias.
#:
#: Test the LEVEL, not the spatial scatter.  Scatter looks like the obvious
#: discriminator -- a real B1 flat has bright orders against dark gaps and
#: reaches IM01SD01 = 4,300-12,400 where an unlit frame gives 1.4 -- but it is
#: decker-dependent: the wider B2 slit fills the frame far more evenly and its
#: perfectly good flats sit at 149-171.  A scatter threshold set on B1 throws
#: every B2 flat away.
MIN_FLAT_LEVEL = 3000.

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

#: Frames that are science rather than calibration.  Matched as a SUBSTRING of
#: the raw TARGNAME card, not for equality: on 1998-08-12 the observer typed
#: `TARGNAME = 'H187123'`.  KOA's own `targname` column normalises that to
#: `187123`, so the archive survey and the downloader both find the frames and
#: only a reader of the raw header sees the typo.  With an equality test that
#: night -- the flat donor for the whole August-September block, and the one
#: carrying the cell-in/cell-out pair -- silently has no science frames.
SCIENCE_TARGET = '187123'


# ---------------------------------------------------------------------------
# Staging
# ---------------------------------------------------------------------------

#: PypeIt's echelle-angle tolerance for grouping frames into one configuration
#: (`keck_hires.py`: rtol 1e-3, atol 1e-2).  A flat outside this of the science
#: frame lands in a different configuration and is simply not seen.
ECH_RTOL, ECH_ATOL = 1e-3, 1e-2


def frame_roles(raw_dir):
    """ Classify every raw frame in a directory from its own FITS header.

    The archive survey is used to decide what to *download*; what to *reduce*
    is decided from the headers of the files actually on disk, so a staging
    directory can never disagree with the frames it contains.

    Args:
        raw_dir (str): a night's raw directory.

    Returns:
        list: dicts with 'path', 'role', 'decker', 'mjd'.
    """
    from astropy.io import fits

    out = []
    for path in sorted(glob.glob(os.path.join(raw_dir, 'HI.*.fits'))):
        hdr = fits.getheader(path)
        decker = str(hdr.get('DECKNAME', '')).strip()
        mjd = float(hdr.get('MJD', 0.))
        ech = float(hdr.get('ECHANGL', np.nan))
        if SCIENCE_TARGET in str(hdr.get('TARGNAME', '')):
            role = 'science'
        elif hdr.get('LAMPCAT1') or hdr.get('LAMPCAT2'):
            role = 'arc'
        elif (hdr.get('LAMPQTZ2') or str(hdr.get('LAMPNAME', '')).strip() == 'quartz1'):
            # Clean means hatch closed and the iodine cell out of the beam: a
            # hatch-open flat sees the sky, and an iodine-in flat would imprint
            # the I2 forest on the flat field and partly divide it out of the
            # science frames, which is fatal for phase 3.
            #
            # XCOVOPEN matters just as much and is not in any archive
            # classification.  Five of 1998-08-12's eighteen "Narrowflat"
            # frames were taken with the cross-disperser cover shut: the
            # archive calls them flatlamp, iodine out, hatch closed, and they
            # contain nothing but bias -- mean 768 counts with a standard
            # deviation of 1.4, against 17,400 and 12,000 for a real one.
            # PypeIt silently declines to frametype them, and the reduction
            # dies much later with "No frames of type=trace provided".
            lit = (hdr.get('XCOVOPEN')
                   and float(hdr.get('IM01MN01', 0.)) > MIN_FLAT_LEVEL)
            clean = ((not hdr.get('HATOPEN')) and (not hdr.get('IODIN'))
                     and lit)
            role = 'flat' if clean else 'flat_dirty'
        else:
            role = 'other'
        out.append(dict(path=path, role=role, decker=decker, mjd=mjd,
                        ech=ech))
    return out


def choose_flats(night):
    """ The four clean flats this night's reduction will use.

    Args:
        night (str): UT night, 'YYYYMMDD'.

    Returns:
        tuple: (list of paths, provenance string)
    """
    subdir = NIGHTS[night][0]
    own = frame_roles(os.path.join(RAW_ROOT, subdir))
    science = [f for f in own if f['role'] == 'science']
    if not science:
        raise RuntimeError('No science frame for {:s}'.format(night))
    deck = science[0]['decker']
    anchors = [f['mjd'] for f in science]

    def nearest(cands):
        return sorted(cands, key=lambda f: min(abs(f['mjd'] - a) for a in anchors))

    if night in FORCED_CALIBS:
        forced = [os.path.join(RAW_ROOT, subdir, b) for b in FORCED_CALIBS[night]]
        missing = [f for f in forced if not os.path.isfile(f)]
        if missing:
            raise RuntimeError('Missing forced calibration: {:s}'.format(
                ', '.join(missing)))
        return forced, '{:s} own {:s}, iodine-in, forced as trace'.format(night, deck)

    mine = [f for f in own if f['role'] == 'flat' and f['decker'] == deck]
    if len(mine) >= N_FLATS:
        chosen = nearest(mine)[:N_FLATS]
        return [f['path'] for f in chosen], '{:s} own {:s}'.format(night, deck)


    donor = frame_roles(os.path.join(RAW_ROOT, NIGHTS[DONOR_NIGHT][0]))
    theirs = [f for f in donor if f['role'] == 'flat' and f['decker'] == deck]
    if len(theirs) < N_FLATS:
        return [], 'NONE ({:s} has {:d} own {:s} flats, donor {:s} has {:d})'.format(
            night, len(mine), deck, DONOR_NIGHT, len(theirs))
    chosen = nearest(theirs)[:N_FLATS]

    # Would PypeIt even put these in the same configuration?  If not, say so
    # here rather than letting it surface as an opaque "No frames of
    # type=trace" three minutes into the reduction.
    ech_sci = float(np.median([f['ech'] for f in science]))
    dech = max(abs(f['ech'] - ech_sci) for f in chosen)
    if dech > ECH_ATOL + ECH_RTOL * abs(ech_sci):
        return ([f['path'] for f in chosen],
                'ILLEGAL: nearest donor {:s} is d_echangle {:.4f} away, '
                'outside PypeIt\'s {:.3f} tolerance -- expect no trace frames'
                .format(DONOR_NIGHT, dech, ECH_ATOL))
    return ([f['path'] for f in chosen],
            '{:s} {:s} (d_echangle {:.4f})'.format(DONOR_NIGHT, deck, dech))


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
    if night in JULY_RUN:
        # Phase 1's five nights keep the exact flat set they were reduced with,
        # so their results stay comparable with what is already reported.
        sources += [os.path.join(RAW_ROOT, FLAT_NIGHT, f) for f in FLAT_FRAMES]
    else:
        flats, provenance = choose_flats(night)
        print('  flats: {:s}'.format(provenance))
        sources += [f for f in flats if f not in sources]

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


def force_frametypes(pypeit_file, forced):
    """ Set an explicit frametype on frames `pypeit_setup` left untyped.

    `pypeit_setup` writes every staged frame into the data block but comments
    out, and gives `frametype = None` to, anything it could not classify.  An
    iodine-in quartz flat is one of those.  The frametype column is
    authoritative at `run_pypeit` time, so writing it here is enough.

    Args:
        pypeit_file (str): the generated file, modified in place.
        forced (dict): basename -> frametype string.

    Returns:
        list: basenames that were found and rewritten.
    """
    with open(pypeit_file) as fh:
        lines = fh.read().split('\n')

    done = []
    for i, line in enumerate(lines):
        bare = line.lstrip('#').strip()
        for base, ftype in forced.items():
            if not bare.startswith(base):
                continue
            fields = bare.split('|')
            fields[0] = base
            fields[1] = ftype
            fields[-1] = '0'
            lines[i] = ' ' + ' | '.join(f.strip() for f in fields)
            done.append(base)
    if len(done) != len(forced):
        raise RuntimeError('Could not find {:s} in {:s}'.format(
            ', '.join(set(forced) - set(done)), pypeit_file))

    with open(pypeit_file, 'w') as fh:
        fh.write('\n'.join(lines))
    return done


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
    if night in FORCED_CALIBS:
        forced = force_frametypes(dst, FORCED_CALIBS[night])
        print('  forced frametypes: {:s}'.format(', '.join(forced)))
        with open(dst) as fh:
            text = fh.read()
        with open(dst, 'w') as fh:
            fh.write(text.replace('\n# Setup\n',
                                  '\n' + DECEMBER_PARAMS + '# Setup\n', 1))

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
            decker=str(sobjs.header.get('DECKER', '')).strip(),
            mjd=float(sobjs.header.get('MJD', 0.)),
            exptime=float(sobjs.header.get('EXPTIME', 0.)),
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
    print('NIGHT-TO-NIGHT STABILITY, HD 187123 discovery era')
    print('=' * 96)
    print('{:<10s} {:>6s} {:>4s} {:>7s} {:>9s} {:>8s} {:>8s} {:>8s}   {:s}'.format(
        'night', 'frames', 'deck', 'orders', 'orders', 'S/N min', 'S/N med',
        'S/N max', 'wavecal RMS (px)'))
    for s in summaries:
        if s is None:
            continue
        for i, fr in enumerate(s['frames']):
            rms = ''
            if i == 0 and s['rms_med'] is not None:
                rms = '{:.3f} / {:.3f} / {:.3f}  (n={:d})'.format(
                    s['rms_min'], s['rms_med'], s['rms_max'], s['n_solutions'])
            print('{:<10s} {:>6s} {:>4s} {:>7d} {:>9s} {:>8.1f} {:>8.1f} {:>8.1f}   {:s}'.format(
                s['night'] if i == 0 else '', '{:d}/{:d}'.format(i + 1, s['n_frames']),
                fr.get('decker', ''), fr['n_orders'],
                '{:d}-{:d}'.format(fr['order_min'], fr['order_max']),
                fr['snr_min'], fr['snr_med'], fr['snr_max'], rms))
    missing = [s for s in summaries if s is None]
    print('=' * 96)


# ---------------------------------------------------------------------------

def main():
    """ Reduce every night of the run and report stability. """
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--night', action='append', choices=sorted(NIGHTS),
                        help='UT night to reduce (repeatable; default: all)')
    parser.add_argument('--era', action='store_true',
                        help='reduce only the nights added in phase 2 prompt 3')
    parser.add_argument('--overwrite', action='store_true',
                        help='rebuild staging, setup and reduction directories')
    parser.add_argument('--summarise-only', action='store_true',
                        help='skip the reductions; just report what is on disk')
    args = parser.parse_args()

    nights = args.night or (list(ERA_NIGHTS) if args.era
                            else sorted(NIGHTS))

    failures = {}
    if not args.summarise_only:
        for night in nights:
            try:
                reduce_night(night, overwrite=args.overwrite)
            except Exception as exc:
                # One night failing must not abandon the other eighteen; the
                # point of prompt 3 is to find out which nights need help.
                failures[night] = str(exc)
                print('  FAILED: {:s}'.format(str(exc)))
    if failures:
        print('\n{:d} night(s) failed:'.format(len(failures)))
        for night, msg in sorted(failures.items()):
            print('  {:s}: {:s}'.format(night, msg))

    summaries = {n: summarise(n) for n in sorted(NIGHTS)}
    print_summary([summaries[n] for n in sorted(NIGHTS)])

    absent = [n for n in sorted(NIGHTS) if summaries[n] is None]
    if absent:
        print('\nNo reduction on disk for {:d} night(s): {:s}'.format(
            len(absent), ', '.join(absent)))
    print('{:d} of {:d} nights reduced; {:d} science frames extracted.'.format(
        len(NIGHTS) - len(absent), len(NIGHTS),
        sum(s['n_frames'] for s in summaries.values() if s)))

    out = os.path.join(REDUX_ROOT, 'run_summary.json')
    with open(out, 'w') as fh:
        json.dump([s for s in summaries.values() if s], fh, indent=2)
    print('Wrote {:s}'.format(out))


if __name__ == '__main__':
    main()
