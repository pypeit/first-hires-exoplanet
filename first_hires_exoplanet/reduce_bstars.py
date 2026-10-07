""" Reduce the 1998-08-26 B stars observed through the iodine cell (phase 3,
after prompt 6).

The Q&A after prompt 5 decided the B stars are reduced before the iodine
model is examined on them.  They are the frames prompt 4 fetched into
``raw/1998aug26_lsf/``: HR 7236 (08:21 UT), HR 8634 (10:20, 11:23, 13:16 UT)
and HR 838 (15:25 UT), 5-20 s each, cell in, in exactly the science
configuration (B1, echelle -0.0051163, cross-disperser -0.544).

They are reduced as 1998-08-26's HD 187123 frames were -- same ThAr arc, same
donor flats from 1998-08-12, same parameters (`reduce_run.PARAM_BLOCK`), plus
one: ``force_center_obj``, because these stars over-fill the slit (below) --
and the night's calibrations are copied in rather than rebuilt, since they
come from exactly those frames.  Everything is staged and reduced beside the
night's own directories, never inside them, so the HD 187123 reduction is
untouched.

The iodine-in quartz flat taken that night is *not* reduced here: a lamp fills
the slit uniformly, which this extraction (built for a star) would treat as
sky.  It is not needed for what the B stars are for.

Run with:

    conda run --no-capture-output -n pypeit14 \\
        python -m first_hires_exoplanet.reduce_bstars

"""

# Standard imports
import os
import glob
import shutil
import argparse

from first_hires_exoplanet import reduce_run as rr


NIGHT = '19980826'
LSF_RAW = os.path.join(rr.RAW_ROOT, '1998aug26_lsf')
STAGE = os.path.join(rr.REDUX_ROOT, 'stage_{:s}_lsf'.format(NIGHT))
SETUP = os.path.join(rr.REDUX_ROOT, 'setup_{:s}_lsf'.format(NIGHT))
REDUX = os.path.join(rr.REDUX_ROOT, 'reduce_{:s}_lsf'.format(NIGHT))
NIGHT_REDUX = os.path.join(rr.REDUX_ROOT, 'reduce_{:s}'.format(NIGHT))
NIGHT_STAGE = os.path.join(rr.REDUX_ROOT, 'stage_{:s}'.format(NIGHT))


def calibration_frames():
    """ The arc and flats the night's own reduction used, read from its file. """
    pyp = glob.glob(os.path.join(NIGHT_REDUX, '*.pypeit'))[0]
    frames = []
    with open(pyp) as fh:
        for line in fh:
            parts = [p.strip() for p in line.split('|')]
            if len(parts) > 2 and parts[0].startswith('HI.') and 'science' not in parts[1]:
                frames.append(os.path.realpath(os.path.join(NIGHT_STAGE, parts[0])))
    return frames


def bstar_frames():
    """ The B-star exposures; not the iodine quartz flat. """
    from astropy.io import fits
    out = []
    for f in sorted(glob.glob(os.path.join(LSF_RAW, 'HI.*.fits'))):
        if str(fits.getheader(f).get('IMAGETYP', '')).strip() == 'object':
            out.append(f)
    return out


def stage(overwrite=False):
    if os.path.isdir(STAGE):
        if not overwrite:
            return STAGE
        shutil.rmtree(STAGE)
    os.makedirs(STAGE)
    frames = calibration_frames() + bstar_frames()
    for src in frames:
        if not os.path.isfile(src):
            raise RuntimeError('Missing raw frame: {:s}'.format(src))
        os.symlink(src, os.path.join(STAGE, os.path.basename(src)))
    print('  staged {:d} frames ({:d} calibration, {:d} B-star) in {:s}'.format(
        len(frames), len(calibration_frames()), len(bstar_frames()), STAGE))
    return STAGE


#: A B star at V ~ 3-4 in 5-20 s over-fills this 3.5-arcsec slit, and
#: PypeIt's peak-finding then finds it in almost no order: the first attempt
#: found HR 8634 (37323) in one order of 37 and the echelle extraction
#: stopped.  The remedy is the one built for Shane/Hamspec's slit-filling
#: stars (PypeIt `FindObjPar.force_center_obj`, present at 017bece06): skip
#: peak-finding and extract one object at each order's centre, boxcar across
#: the order.  Sky subtraction is already off (`reduce_run.PARAM_BLOCK`).
#: Applied to the B stars only; the HD 187123 reductions are unchanged.
FORCE_CENTER = '        force_center_obj = True\n'


def force_center(pypeit_file):
    """ Add ``force_center_obj = True`` to the ``[[findobj]]`` block. """
    with open(pypeit_file) as fh:
        text = fh.read()
    if 'force_center_obj' in text:
        return
    marker = '    [[findobj]]\n'
    if marker not in text:
        raise RuntimeError('No [[findobj]] block in {:s}'.format(pypeit_file))
    with open(pypeit_file, 'w') as fh:
        fh.write(text.replace(marker, marker + FORCE_CENTER, 1))


def reduce(overwrite=False):
    stage(overwrite=overwrite)
    for path in (SETUP, REDUX):
        if overwrite and os.path.isdir(path):
            shutil.rmtree(path)
    os.makedirs(SETUP, exist_ok=True)
    os.makedirs(REDUX, exist_ok=True)

    print('  pypeit_setup ...')
    rr._run(['pypeit_setup', '-s', 'keck_hires_orig', '-r', STAGE, '-c', 'all',
             '-d', SETUP], log=os.path.join(SETUP, 'pypeit_setup.log'))
    src = rr.science_setup(SETUP)
    dst = os.path.join(REDUX, os.path.basename(src))
    shutil.copy(src, dst)
    rr.inject_params(dst)
    force_center(dst)

    # Reuse the night's calibrations: built from exactly these arc and flats
    calib = os.path.join(REDUX, 'Calibrations')
    if not os.path.isdir(calib):
        shutil.copytree(os.path.join(NIGHT_REDUX, 'Calibrations'), calib)

    science = os.path.join(REDUX, 'Science')
    if os.path.isdir(science):
        shutil.rmtree(science)

    print('  run_pypeit ...')
    rr._run(['run_pypeit', os.path.basename(dst), '-o'], cwd=REDUX,
            log=os.path.join(REDUX, 'run_pypeit.log'))
    spec1d = sorted(glob.glob(os.path.join(science, 'spec1d_*.fits')))
    print('  {:d} spec1d files in {:s}'.format(len(spec1d), science))
    return spec1d


def main():
    p = argparse.ArgumentParser(description='Reduce the 1998-08-26 B stars')
    p.add_argument('--overwrite', action='store_true')
    args = p.parse_args()
    reduce(overwrite=args.overwrite)


if __name__ == '__main__':
    main()
