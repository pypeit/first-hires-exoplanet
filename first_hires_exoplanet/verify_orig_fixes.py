#!/usr/bin/env python
"""
Verification driver for the ``keck_hires_orig`` fixes (prompt 6).

Runs PypeIt's ``keck_hires_orig`` spectrograph class over the real 1998
frames downloaded in prompt 5 and prints everything the fixes touch:
metadata, raw-image sections, detector parameters, the bad-pixel mask,
``check_spectrograph`` and automatic frame typing.  It also prints a
regression check on the post-2004 ``keck_hires`` class.

Run it before and after editing PypeIt and diff the outputs:

    conda run -n pypeit14 python first_hires_exoplanet/verify_orig_fixes.py > before.txt
    ... edit PypeIt ...
    conda run -n pypeit14 python first_hires_exoplanet/verify_orig_fixes.py > after.txt

Nothing is written to disk; ``--par-dump FILE`` writes the fully resolved
``keck_hires`` parameter tree for a before/after diff of the mosaic class.
"""
import argparse
import glob
import io
import os
import sys
import warnings

import numpy as np
from astropy.io import fits

warnings.filterwarnings('ignore')

RAW = '/Users/xavier/Projects/PypeIt/first-hires-exoplanet-data/raw'
SCIENCE = f'{RAW}/1998jul16/HI.19980716.52985.fits'
ARC = f'{RAW}/1998jul16/HI.19980716.54321.fits'
FLAT_B1 = f'{RAW}/1998jul14/HI.19980714.07711.fits'
FLAT_B2 = f'{RAW}/1998jul16/HI.19980716.54482.fits'
TWOAMP = f'{RAW}/inspect_only/HI.19980716.15422.fits'
TEST_FRAMES = [('science 400s B1', SCIENCE), ('ThAr arc 10s B1', ARC),
               ('clean B1 flat 3s (07-14)', FLAT_B1), ('B2 flat 2s hatch closed', FLAT_B2),
               ('2-amp N01H frame', TWOAMP)]
META_KEYS = ['idname', 'binning', 'decker', 'dispname', 'echangle', 'xdangle', 'exptime', 'mjd']


def col_ranges(sec_img):
    """Return {amp: (first_col, last_col)} for a section image."""
    out = {}
    for amp in np.unique(sec_img[sec_img > 0]):
        cols = np.where(np.any(sec_img == amp, axis=0))[0]
        out[int(amp)] = (int(cols.min()), int(cols.max()), int(cols.size))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--par-dump', default=None, help='write the resolved keck_hires par tree here')
    args = ap.parse_args()

    from pypeit import io as pio
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit.core import parse
    from pypeit.metadata import PypeItMetaData

    spec = load_spectrograph('keck_hires_orig')
    print(f'PypeIt from {os.path.dirname(sys.modules["pypeit"].__file__)}')

    # 1. Metadata for every test frame
    print('\n=== get_meta_value ===')
    for label, f in TEST_FRAMES:
        vals = {k: spec.get_meta_value(f, k) for k in META_KEYS}
        print(f'{label:28s} {os.path.basename(f)}')
        print('    ' + '  '.join(f'{k}={vals[k]!r}' for k in META_KEYS))

    # 2. get_rawimage sections on the science frame and the 2-amp frame
    print('\n=== get_rawimage ===')
    for label, f in [('science', SCIENCE), ('2-amp', TWOAMP)]:
        det, img, hdu, exptime, rdsec, ossec = spec.get_rawimage(f, 1)
        hdr = hdu[0].header
        print(f'{label}: NAXIS1={hdr["NAXIS1"]} PREPIX={hdr["PREPIX"]} POSTPIX={hdr["POSTPIX"]} '
              f'NUMAMPS={hdr["NUMAMPS"]} image {img.shape} exptime {exptime}')
        print(f'    datasec cols (amp: first..last, n) : {col_ranges(rdsec)}')
        print(f'    oscansec cols (amp: first..last, n): {col_ranges(ossec)}')
        print(f'    gain {det["gain"]} ronoise {det["ronoise"]} binning {det["binning"]!r} '
              f'specaxis {det["specaxis"]} platescale {det["platescale"]}')
        # Column medians at the datasec boundaries, to prove the sections land on the data
        med = np.median(img, axis=0)
        for amp, (c0, c1, n) in col_ranges(rdsec).items():
            print(f'    amp{amp} column medians just outside/inside datasec: '
                  f'[{c0-1}]={med[c0-1]:.0f} [{c0}]={med[c0]:.0f} ... [{c1}]={med[c1]:.0f} [{c1+1}]={med[c1+1]:.0f}')
        hdu.close()

    # 3. Detector par
    print('\n=== get_detector_par(1, hdu) on science ===')
    hdu = pio.fits_open(SCIENCE)
    det = spec.get_detector_par(1, hdu=hdu)
    print(f'    gain {det["gain"]} ronoise {det["ronoise"]} numamplifiers {det["numamplifiers"]} '
          f'binning {det["binning"]!r} platescale {det["platescale"]}')
    hdu.close()

    # 4. order_platescale with the frame's binning string
    binning = spec.get_meta_value(SCIENCE, 'binning')
    print(f'\n=== order_platescale(binning={binning!r}) === '
          f'{spec.order_platescale(np.arange(3), binning)[0]:.3f} arcsec/binned spatial pixel')

    # 5. BPM
    print('\n=== bpm(science, 1) ===')
    log_stream = io.StringIO()
    import logging
    handler = logging.StreamHandler(log_stream)
    logging.getLogger('pypeit').addHandler(handler)
    bpm = spec.bpm(SCIENCE, 1)
    logging.getLogger('pypeit').removeHandler(handler)
    print(f'    shape {bpm.shape}  masked fraction {bpm.mean():.4%}')
    full_cols = np.where(np.all(bpm == 1, axis=1))[0]
    print(f'    fully masked spectral columns (trimmed, 0-indexed): {full_cols.tolist()}')
    part = np.where((bpm.sum(axis=1) > 0) & (bpm.sum(axis=1) < bpm.shape[1]))[0]
    print(f'    partially masked spectral columns: n={part.size}, '
          f'ranges {compress_ranges(part)}')
    for c in (1010, 1030):   # ink-spot centre column after / before the prescan-origin fix
        rows = np.where(bpm[c] == 1)[0]
        print(f'    spatial rows flagged in trimmed spectral column {c}: {compress_ranges(rows)}')
    warns = [l for l in log_stream.getvalue().splitlines() if 'WARN' in l.upper() or 'warn' in l]
    print(f'    log warnings: {warns if warns else "none"}')
    # BPM with no example file (shape only) still works
    bpm0 = spec.bpm(None, 1, shape=(2048, 2048))
    print(f'    bpm(None, 1, shape=(2048,2048)) -> shape {bpm0.shape} masked {bpm0.mean():.4%} '
          f'fully masked cols {np.where(np.all(bpm0 == 1, axis=1))[0].tolist()}')

    # 6. check_spectrograph
    print('\n=== check_spectrograph ===')
    try:
        spec.check_spectrograph(SCIENCE)
        print('    science frame: passes (no exception)')
    except Exception as e:
        print(f'    science frame: RAISED {type(e).__name__}: {e}')
    # Fake a post-2004 date to exercise the other branch
    hdr = fits.getheader(SCIENCE)
    hdr['DATE-OBS'] = '2010-01-01'
    try:
        spec.check_spectrograph(hdr)
        print('    post-2004 header on keck_hires_orig: passes (NO exception)')
    except Exception as e:
        print(f'    post-2004 header on keck_hires_orig: RAISED {type(e).__name__}: {e}')

    # 7. Frame typing over every frame in the reduction directories
    print('\n=== automatic frame typing (PypeItMetaData.get_frame_types) ===')
    files = sorted(glob.glob(f'{RAW}/1998jul1[46]/*.fits'))
    par = spec.default_pypeit_par()
    print(f'    scienceframe exprng {par["scienceframe"]["exprng"]}  '
          f'standardframe exprng {par["calibrations"]["standardframe"]["exprng"]}')
    fitstbl = PypeItMetaData(spec, par, files=files, strict=True)
    fitstbl.get_frame_types(flag_unknown=True)
    for row in fitstbl.table:
        print(f'    {row["filename"]:26s} {str(row["exptime"]):>6s}s {row["decker"]:3s} hatch={row["hatch"]!s:5s} '
              f'idname={row["idname"]!s:8s} -> {row["frametype"]}')
    from pypeit.core import standard
    sel = np.array([os.path.basename(SCIENCE) in f for f in fitstbl.table['filename']])
    ra, dec = fitstbl.table['ra'][sel][0], fitstbl.table['dec'][sel][0]
    print(f'    is HD 187123 (ra={ra:.4f}, dec={dec:.4f}) within 10 arcmin of an archive standard? '
          f'{standard.get_archive_standard(ra, dec, tol=10., check=True)}')

    # 8. Regression on the post-2004 class
    print('\n=== keck_hires (post-2004) regression ===')
    hires = load_spectrograph('keck_hires')
    fake = fits.Header()
    fake['BINNING'] = '2,1'
    print(f'    keck_hires.compound_meta(BINNING="2,1", binning) -> {hires.compound_meta([fake], "binning")!r}')
    fake['BINNING'] = '1,2'
    print(f'    keck_hires.compound_meta(BINNING="1,2", binning) -> {hires.compound_meta([fake], "binning")!r}')
    fake_orig = fits.Header()
    fake_orig['BINNING'] = '1,2'
    print(f'    keck_hires_orig.compound_meta(BINNING="1,2", binning) -> {spec.compound_meta([fake_orig], "binning")!r}')
    hpar = hires.default_pypeit_par()
    print(f'    keck_hires scienceframe exprng {hpar["scienceframe"]["exprng"]}, '
          f'ndet {hires.ndet}, MRO {[c.__name__ for c in type(hires).__mro__[:3]]}')
    if args.par_dump:
        with open(args.par_dump, 'w') as fh:
            fh.write('\n'.join(hpar.to_config()))
        print(f'    resolved keck_hires par tree written to {args.par_dump}')


def compress_ranges(idx):
    idx = np.asarray(idx)
    if idx.size == 0:
        return '[]'
    breaks = np.where(np.diff(idx) > 1)[0]
    starts = np.r_[idx[0], idx[breaks + 1]]
    ends = np.r_[idx[breaks], idx[-1]]
    return ', '.join(f'{s}-{e}' if s != e else f'{s}' for s, e in zip(starts, ends))


if __name__ == '__main__':
    main()
