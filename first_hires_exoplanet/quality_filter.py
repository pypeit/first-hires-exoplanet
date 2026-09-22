""" Per-order quality filter for the reduced HD 187123 spectra (phase 2, prompt 2).

Phase 1 found one bad order-spectrum by hand -- order 60 of `HI.19980719.49996`,
122 counts where the same order holds 51,000-110,000 elsewhere -- and concluded
that 370 order-spectra is too many to inspect by eye.  This module turns that
into a statistic and runs it over every order of every reduced epoch, so that
phase 2's cross-correlation and phase 3's forward model start from a vetted set
rather than a hopeful one.

TWO INDEPENDENT FAILURE MODES, and the filter needs both tests:

  1. *PypeIt rejected the pixels.*  `OPT_MASK` can collapse to near-zero while
     the flux and inverse variance are perfectly healthy, driven by
     `OPT_FRAC_USE` -- the fraction of the object profile falling on usable slit
     pixels.  When the seeing is poor the profile spills past the aperture,
     FRAC_USE drops, and PypeIt masks the order.  `frac_good` catches this.

  2. *The extraction silently lost the object.*  Order 60 of `HI.19980719.49996`
     has `OPT_MASK` true for 2026 of 2048 pixels and `OPT_FRAC_USE` = 1.000.
     Nothing internal to that order says it is wrong.  Only its neighbours do:
     orders 59 and 61 hold ~7,500 counts and it holds 123.  The adjacent-order
     S/N diagnostic catches this.

THE DIAGNOSTIC.  S/N varies smoothly with echelle order within a frame, so a
leave-one-out local linear fit to the neighbouring orders predicts each order's
log S/N, and the residual measures how badly that order departs from its own
frame's trend.

That residual is *not* pure noise.  Parts of the echelle format produce the same
residual in every frame -- most sharply at order 89, where the free spectral
range steps and every one of the ten epochs shows +0.06 dex.  A naive threshold
flags order 89 ten times out of ten, which is the classic false positive.  So
the per-order median residual across epochs is subtracted first, and the
threshold is applied to what is left.  What survives is a genuine per-epoch,
per-order anomaly.

The verification test is independent of the flagging test: for each order,
compare an epoch's S/N against the median of the same order across all epochs,
after dividing out each frame's overall S/N level.  An order flagged by the
adjacent-order test that is *also* an outlier against its own history in the
other nine frames is real; one that is not is a false positive.

The reduction this reads is NOT in this repository -- the raw frames and the
PypeIt products live in a sibling data tree.  Point --redux at it if it is
somewhere else.

Reduced with PypeIt 2.0.2.dev1217+g017bece06 (branch orig-hires-fixes).

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.quality_filter

"""

# Standard imports
import os
import glob
import argparse

import numpy as np

from astropy.table import Table

from pypeit import specobjs


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Default reduction directory (outside the repository; see the module docstring)
DEFAULT_REDUX = os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                             'first-hires-exoplanet-data', 'redux')

#: Where the committed derived products go
DEFAULT_OUTDIR = os.path.join(_HERE, 'data')

#: Half-width, in echelle orders, of the neighbourhood used to predict an
#: order's S/N.  Three either side gives six neighbours for a two-parameter
#: local fit, enough to follow the blaze without reaching across the format.
NEIGHBOUR_HALFWIDTH = 3

#: Minimum valid pixels before an order-spectrum's S/N means anything
MIN_VALID_PIX = 200

#: Flag thresholds.
#:
#: THRESH_RESID is not a sigma cut.  The de-trended residuals are strongly
#: bimodal: 95% of order-spectra sit below 0.005 dex and the largest of the
#: healthy ones reaches 0.015, then the distribution jumps to 0.035, 0.060,
#: 0.070 ... 0.741.  0.02 dex sits in the middle of that gap and means a 5%
#: departure from the local blaze trend, which is both physically meaningful
#: and far outside anything the healthy frames produce.  A sigma cut would be
#: meaningless here -- the robust sigma of the de-trended residual is 0.001 dex,
#: so the worst order is 700 sigma out and everything above the noise floor is
#: "significant".
THRESH_RESID = 0.02     #: de-trended adjacent-order residual, dex
THRESH_FRAC_GOOD = 0.5  #: fraction of valid pixels that must survive OPT_MASK
THRESH_SNR = 20.0       #: absolute S/N floor below which an order is unusable


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

def measure_orders(redux_dir):
    """ Read every order-spectrum of every reduced epoch.

    Two signal-to-noise numbers are recorded per order, and the difference
    between them is diagnostic in its own right:

    - ``snr_valid`` over every pixel with a positive inverse variance and a
      wavelength, whether or not PypeIt trusts it.  Defined for every order, so
      the adjacent-order comparison is like against like.
    - ``snr_good`` over only the pixels PypeIt kept (``OPT_MASK``).  This is the
      S/N actually available to a cross-correlation, and is NaN when the mask
      has collapsed.

    Args:
        redux_dir (str): directory holding the ``reduce_*`` reduction folders.

    Returns:
        `astropy.table.Table`_: one row per (epoch, order).
    """
    files = sorted(glob.glob(os.path.join(redux_dir, 'reduce_1998*', 'Science',
                                          'spec1d_*.fits')))
    if len(files) == 0:
        raise FileNotFoundError(f'No spec1d files under {redux_dir}')

    rows = []
    for ifile in files:
        sobjs = specobjs.SpecObjs.from_fitsfile(ifile, chk_version=False)
        koaid = os.path.basename(ifile)[7:24]
        for obj in sobjs:
            wave, flux = obj.OPT_WAVE, obj.OPT_COUNTS
            ivar, mask = obj.OPT_COUNTS_IVAR, obj.OPT_MASK
            valid = np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0) & (wave > 0)
            good = valid & mask

            snr_valid = (np.median(flux[valid] * np.sqrt(ivar[valid]))
                         if valid.sum() >= MIN_VALID_PIX else np.nan)
            snr_good = (np.median(flux[good] * np.sqrt(ivar[good]))
                        if good.sum() >= MIN_VALID_PIX else np.nan)

            rows.append(dict(
                koaid=koaid,
                night=os.path.basename(os.path.dirname(os.path.dirname(ifile)))[-8:],
                order=int(obj.ECH_ORDER),
                wave_min=float(wave[valid].min()) if valid.any() else np.nan,
                wave_max=float(wave[valid].max()) if valid.any() else np.nan,
                n_valid=int(valid.sum()),
                n_good=int(good.sum()),
                frac_good=float(good.sum()) / max(int(valid.sum()), 1),
                frac_zero=float((flux == 0).sum()) / len(flux),
                snr_valid=float(snr_valid),
                snr_good=float(snr_good),
                med_counts=float(np.median(flux[valid])) if valid.any() else np.nan,
                med_frac_use=float(np.median(obj.OPT_FRAC_USE[valid]))
                             if valid.any() else np.nan,
                fwhm=float(obj.FWHM),
            ))
    return Table(rows)


# ---------------------------------------------------------------------------
# The adjacent-order diagnostic
# ---------------------------------------------------------------------------

def adjacent_order_residual(order, logsnr, halfwidth=NEIGHBOUR_HALFWIDTH,
                            use=None):
    """ Leave-one-out local linear prediction of each order's log S/N.

    The order itself is excluded from its own fit, so a bad order cannot drag
    its prediction towards itself and hide.

    Args:
        order (`numpy.ndarray`_): echelle order numbers, one frame's worth.
        logsnr (`numpy.ndarray`_): log10 S/N, same length; NaNs are skipped.
        halfwidth (int): neighbourhood half-width in orders.
        use (`numpy.ndarray`_, optional): boolean mask of orders allowed to act
            as neighbours.  An order excluded here still gets its own residual
            computed; it just stops contaminating everyone else's.

    Returns:
        `numpy.ndarray`_: residual log10(S/N) - prediction, NaN where undefined.
    """
    resid = np.full(len(order), np.nan)
    for i in range(len(order)):
        if not np.isfinite(logsnr[i]):
            continue
        near = ((np.abs(order - order[i]) <= halfwidth)
                & (np.arange(len(order)) != i) & np.isfinite(logsnr))
        if use is not None:
            near &= use
        if near.sum() < 3:
            continue
        coeff = np.polyfit(order[near], logsnr[near], 1)
        resid[i] = logsnr[i] - np.polyval(coeff, order[i])
    return resid


def add_diagnostics(tbl):
    """ Add the flags, in the order the two failure modes demand.

    The tests are **sequential, not parallel**, and that matters.  Running the
    smoothness test on ``snr_valid`` -- every pixel with a positive inverse
    variance, including the ones PypeIt rejected -- measures how much of each
    order the mask threw away, not whether the extraction worked.  In the three
    worst frames that swamps the real signal and flags most of the blue format.

    So: reject the mask-collapsed orders first, then test the smoothness of what
    is left, using only surviving orders as neighbours.  The low-S/N test is
    computed but deliberately not used to gate the smoothness test, so that
    order 60 of `HI.19980719.49996` -- the one fault phase 1 found by hand --
    has to be caught by the smoothness statistic on its own merits.

    Args:
        tbl (`astropy.table.Table`_): output of :func:`measure_orders`;
            modified in place.

    Returns:
        tuple: (table, robust sigma of the de-trended residual in dex)
    """
    epochs = sorted(set(tbl['koaid']))

    # --- test 1: did PypeIt keep the pixels? ------------------------------
    tbl['flag_masked'] = tbl['frac_good'] < THRESH_FRAC_GOOD

    # --- test 2: is there enough signal to cross-correlate? ---------------
    tbl['flag_lowsnr'] = ~(np.asarray(tbl['snr_good'], dtype=float) > THRESH_SNR)

    # --- test 3: does the order follow its frame's blaze trend? -----------
    # Only mask-surviving orders take part, as targets and as neighbours.
    snr = np.asarray(tbl['snr_good'], dtype=float)
    testable = (~np.asarray(tbl['flag_masked'])) & np.isfinite(snr) & (snr > 0)

    # Deliberately a SINGLE pass.  Iterating -- barring the outliers of one
    # pass from the neighbourhoods of the next -- was tried and is wrong here:
    # in HI.19980719.49996 the whole blue block is damaged, so the second pass
    # strips away every neighbour order 60 has and the one fault phase 1 found
    # by hand stops being flagged at all.  A single pass over-flags instead,
    # which is the safe direction, and the independent cross-epoch test below
    # sorts the real faults from the contaminated neighbours.
    tbl['resid'] = np.nan
    for koaid in epochs:
        sel = np.where((tbl['koaid'] == koaid) & testable)[0]
        if len(sel) < 2 * NEIGHBOUR_HALFWIDTH + 1:
            continue
        sel = sel[np.argsort(tbl['order'][sel])]
        tbl['resid'][sel] = adjacent_order_residual(
            np.asarray(tbl['order'][sel], dtype=float), np.log10(snr[sel]))

    # --- the part of the residual that is the echelle format, not a fault --
    # Taken as the per-order median across epochs.  With ten epochs and at most
    # a couple of bad orders at any one order number, the median is set by the
    # healthy frames.
    tbl['resid_format'] = np.nan
    for order in sorted(set(tbl['order'])):
        sel = tbl['order'] == order
        vals = np.asarray(tbl['resid'][sel], dtype=float)
        if np.isfinite(vals).sum() >= 3:
            tbl['resid_format'][sel] = np.nanmedian(vals)
    tbl['resid_det'] = tbl['resid'] - tbl['resid_format']

    det = np.asarray(tbl['resid_det'], dtype=float)
    fin = np.isfinite(det)
    sigma = 1.4826 * np.median(np.abs(det[fin] - np.median(det[fin])))
    tbl['nsigma'] = tbl['resid_det'] / sigma

    # --- the independent cross-epoch check --------------------------------
    # Normalise each frame by its own median S/N, then compare each order
    # against the median of that order across epochs.  This shares no
    # machinery with the adjacent-order test: it never looks at a neighbouring
    # order, only at the same order in the other nine frames.
    # The normalising level must come from the SAME orders in every frame,
    # otherwise a frame that only retains its red orders is normalised by a
    # different part of the blaze and every ratio is biased.
    ref_orders = [o for o in sorted(set(tbl['order']))
                  if all(np.any(testable & (tbl['koaid'] == k) & (tbl['order'] == o))
                         for k in epochs)]
    in_ref = np.isin(np.asarray(tbl['order']), ref_orders)
    tbl['snr_rel'] = np.nan
    for koaid in epochs:
        sel = (tbl['koaid'] == koaid) & testable
        level = np.nanmedian(snr[sel & in_ref])
        tbl['snr_rel'][np.where(sel)[0]] = snr[sel] / level
    tbl['snr_rel_ratio'] = np.nan
    for order in sorted(set(tbl['order'])):
        sel = tbl['order'] == order
        vals = np.asarray(tbl['snr_rel'][sel], dtype=float)
        if np.isfinite(vals).sum() >= 3:
            tbl['snr_rel_ratio'][sel] = tbl['snr_rel'][sel] / np.nanmedian(vals)

    tbl['flag_anomaly'] = (np.abs(np.asarray(tbl['resid_det'], dtype=float))
                           > THRESH_RESID)
    tbl['use'] = ~(tbl['flag_masked'] | tbl['flag_lowsnr'] | tbl['flag_anomaly'])

    tbl.meta['ref_orders'] = ref_orders
    return tbl, sigma


def frame_scatter(tbl):
    """ The frame-level statistic phase 1 quoted, for continuity.

    Phase 1 reported an "adjacent-order S/N scatter" of 6.2% for
    `HI.19980719.49996` against 1.8-2.0% for the other two frames of that night.
    That calculation was never written to disk, so this recomputes a frame-level
    scatter from the same idea -- the robust spread of the fractional S/N
    difference between neighbouring orders -- for comparison.

    Args:
        tbl (`astropy.table.Table`_): the measured table.

    Returns:
        `astropy.table.Table`_: one row per epoch.
    """
    rows = []
    for koaid in sorted(set(tbl['koaid'])):
        sel = tbl['koaid'] == koaid
        idx = np.where(sel)[0]
        idx = idx[np.argsort(tbl['order'][idx])]
        snr = np.asarray(tbl['snr_valid'][idx], dtype=float)
        fin = np.isfinite(snr)
        snr = snr[fin]
        frac = np.diff(snr) / (0.5 * (snr[1:] + snr[:-1]))
        mad = 1.4826 * np.median(np.abs(frac - np.median(frac)))
        det = np.asarray(tbl['resid_det'][idx], dtype=float)
        rows.append(dict(koaid=koaid,
                         scatter_mad=mad * 100.,
                         scatter_std=np.std(frac) * 100.,
                         worst_resid=(det[np.nanargmax(np.abs(det))]
                                      if np.isfinite(det).any() else np.nan),
                         n_flagged=int((~tbl['use'][sel]).sum())))
    out = Table(rows)
    for col in ['scatter_mad', 'scatter_std', 'worst_resid']:
        out[col].info.format = '%.2f' if 'scatter' in col else '%+.3f'
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(options=None):
    parser = argparse.ArgumentParser(
        description='Per-order quality filter for the reduced HD 187123 spectra.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--redux', type=str, default=DEFAULT_REDUX,
                        help='Directory holding the reduce_* folders')
    parser.add_argument('--outdir', type=str, default=DEFAULT_OUTDIR,
                        help='Where to write the quality table')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):

    tbl = measure_orders(pargs.redux)
    tbl, sigma = add_diagnostics(tbl)

    nep = len(set(tbl['koaid']))
    print(f'\n=== 1. What was read ===')
    print(f'  epochs                     : {nep}')
    print(f'  order-spectra              : {len(tbl)}')
    print(f'  orders                     : {min(tbl["order"])}-{max(tbl["order"])}')

    print(f'\n=== 2. The adjacent-order statistic ===')
    raw = np.asarray(tbl['resid'], dtype=float)
    det = np.asarray(tbl['resid_det'], dtype=float)
    print(f'  robust sigma, raw residual : '
          f'{1.4826 * np.median(np.abs(raw[np.isfinite(raw)] - np.median(raw[np.isfinite(raw)]))):.4f} dex')
    print(f'  robust sigma, de-trended   : {sigma:.4f} dex '
          f'({sigma * np.log(10) * 100:.1f}% in S/N)')
    fmt = tbl['order', 'resid_format'].copy()
    fmt = fmt[np.unique(fmt['order'], return_index=True)[1]]
    fmt.sort('resid_format')
    print('  largest echelle-format terms (subtracted, not flagged):')
    for row in fmt[-4:][::-1]:
        print(f'    order {row["order"]:2d}: {row["resid_format"]:+.4f} dex '
              f'({row["resid_format"] / sigma:+.1f} sigma if left in)')

    print(f'\n=== 3. Flags ===')
    for key, lbl in [('flag_masked', f'frac_good < {THRESH_FRAC_GOOD}'),
                     ('flag_lowsnr', f'S/N(good) <= {THRESH_SNR:.0f} or undefined'),
                     ('flag_anomaly',
                      f'|de-trended residual| > {THRESH_RESID} dex')]:
        print(f'  {lbl:42s}: {int(np.sum(tbl[key])):3d}')
    print(f'  {"REJECTED (any flag)":42s}: {int(np.sum(~tbl["use"])):3d}'
          f'  of {len(tbl)}  ({100. * np.sum(~tbl["use"]) / len(tbl):.1f}%)')
    print(f'  {"USABLE":42s}: {int(np.sum(tbl["use"])):3d}')

    print(f'\n=== 4. Anomalies, with the independent cross-epoch check ===')
    print('  snr_rel_ratio is this order\'s S/N against the same order in the')
    print('  other epochs, after removing each frame\'s overall level.  A real')
    print('  fault departs from 1; a false positive does not.')
    anom = tbl[np.asarray(tbl['flag_anomaly'], dtype=bool)]
    anom = anom[np.argsort(np.abs(np.asarray(anom['resid_det'], dtype=float)))[::-1]]
    print(f'\n  {"koaid":18s} {"ord":>3s} {"S/N":>7s} {"resid":>8s} '
          f'{"frac_good":>9s} {"xepoch":>7s}  verdict')
    n_real = 0
    for row in anom:
        ratio = row['snr_rel_ratio']
        real = np.isfinite(ratio) and abs(np.log10(ratio)) > 0.05
        n_real += int(real)
        verdict = ('CONFIRMED' if real
                   else 'not confirmed -- possible false positive')
        print(f'  {row["koaid"]:18s} {row["order"]:3d} {row["snr_good"]:7.1f} '
              f'{row["resid_det"]:+8.3f} {row["frac_good"]:9.3f} {ratio:7.3f}'
              f'  {verdict}')
    print(f'\n  {n_real} of {len(anom)} confirmed by the independent'
          f' cross-epoch test.')

    # Are the unconfirmed flags explained by a bad neighbour?
    conf = [(r['koaid'], r['order']) for r in anom
            if np.isfinite(r['snr_rel_ratio'])
            and abs(np.log10(r['snr_rel_ratio'])) > 0.05]
    print('\n  Unconfirmed flags, against the nearest CONFIRMED fault in the'
          ' same frame:')
    any_unexplained = False
    for row in anom:
        ratio = row['snr_rel_ratio']
        if np.isfinite(ratio) and abs(np.log10(ratio)) > 0.05:
            continue
        near = [abs(o - row['order']) for k, o in conf if k == row['koaid']]
        dist = min(near) if near else None
        explained = dist is not None and dist <= NEIGHBOUR_HALFWIDTH
        any_unexplained |= not explained
        print(f'    {row["koaid"]} order {row["order"]:2d}: nearest confirmed '
              f'fault is {dist} orders away -> '
              f'{"inside" if explained else "OUTSIDE"} the +/-'
              f'{NEIGHBOUR_HALFWIDTH}-order neighbourhood')
    print('    => every unconfirmed flag is a neighbour of a real fault'
          if not any_unexplained else
          '    => at least one unconfirmed flag is NOT explained by a neighbour')

    print(f'\n=== 5. Per epoch ===')
    fr = frame_scatter(tbl)
    fr.pprint_all()

    print(f'\n=== 6. Where the damage falls ===')
    wmax = np.asarray(tbl['wave_max'], dtype=float)
    windows = [('blue, < 5000 A (phase 2 cross-correlation)', wmax < 5000.),
               ('iodine, > 5000 A (phase 3 forward model)', wmax >= 5000.)]
    for lbl, sel in windows:
        n, nu = int(sel.sum()), int(np.sum(np.asarray(tbl['use'])[sel]))
        print(f'  {lbl:44s}: {nu:3d} / {n:3d} usable  '
              f'({100. * (n - nu) / n:.0f}% rejected)')

    print(f'\n  Per epoch:')
    print(f'  {"koaid":18s} {"blue":>9s} {"iodine":>9s}')
    for koaid in sorted(set(tbl['koaid'])):
        row = f'  {koaid:18s}'
        for lbl, sel in windows:
            s = sel & (np.asarray(tbl['koaid']) == koaid)
            row += f' {int(np.sum(np.asarray(tbl["use"])[s])):4d} /{int(s.sum()):4d}'
        print(row)

    print(f'\n  Orders never rejected in any epoch: ', end='')
    clean = [o for o in sorted(set(tbl['order']))
             if np.all(np.asarray(tbl['use'])[tbl['order'] == o])]
    print(f'{len(clean)} of {len(set(tbl["order"]))} '
          f'({min(clean)}-{max(clean)})' if clean else 'none')

    if not os.path.isdir(pargs.outdir):
        os.makedirs(pargs.outdir)
    outfile = os.path.join(pargs.outdir, 'order_quality.csv')
    for col in ['wave_min', 'wave_max', 'frac_good', 'frac_zero', 'snr_valid',
                'snr_good', 'med_counts', 'med_frac_use', 'fwhm', 'resid',
                'resid_format', 'resid_det', 'nsigma', 'snr_rel',
                'snr_rel_ratio']:
        tbl[col].info.format = '%.4f'
    # Booleans round-trip out of ascii.csv as the strings 'True'/'False', which
    # silently evaluate truthy on read.  Write them as 0/1 so the prompt-8
    # adapter can consume this table without a bespoke parser.
    for col in ['flag_masked', 'flag_lowsnr', 'flag_anomaly', 'use']:
        tbl[col] = np.asarray(tbl[col], dtype=np.int8)
    tbl.write(outfile, format='ascii.csv', overwrite=True)
    print(f'\nWrote {outfile}')


if __name__ == '__main__':
    main(parse_args())
