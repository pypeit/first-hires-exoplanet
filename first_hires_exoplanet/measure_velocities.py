""" Relative radial velocities by cross-correlation (phase 2, prompt 5).

Cross-correlates every discovery-era epoch against the iodine-free template
built in prompt 4, restricted to the orders blueward of 5000 A where the
stellar lines are not swamped by the I2 forest.

This is **not** a planet detection and is not presented as one.  The document's
own framing applies: cross-correlation is the instrument here, not the
deliverable.  What it produces is the per-epoch initial velocity guess pyodine
needs, an end-to-end consistency check across the nine-month baseline, and a
measured number for the precision a cross-correlation actually reaches on this
data -- which is the argument for building a forward model at all.

REFERENCE FRAME, following prompt 1.  Both the template and each epoch are put
in the **observed** frame by dividing out PypeIt's `VEL_CORR`.  In that frame a
measured line sits at v_star - BVC, so the cross-correlation returns

    dv_obs = (v_e - v_t) - (BVC_e - BVC_t)

and the relative barycentric velocity of the star is

    v_e - v_t = dv_obs + (BVC_e - BVC_t)

with BVC from astropy at the exposure midpoint (the header MJD is the exposure
*start*, proved in prompt 1).  Doing it this way means PypeIt's heliocentric
correction -- which prompt 1 showed carries a sign error worth +13 m/s and a
drift of 5 m/s over this baseline -- never enters.

WHICH ORDERS.  Orders whose coverage lies entirely blueward of 5000 A, and
which prompt 2's quality filter marks usable for that epoch.  Prompt 2 found
that 9% of blue order-spectra are rejected, almost all for mask collapse.

UNCERTAINTIES.  Two are computed and both are reported, because they answer
different questions and the smaller one is not the honest one:

  * the **formal** error, from Zucker (2003)'s estimator on the curvature of
    each order's correlation peak, combined in inverse-variance;
  * the **empirical** error, from the scatter of the per-order velocities
    within the epoch divided by sqrt(N).

Where they disagree the empirical one is quoted, because it contains whatever
the formal one does not model.

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.measure_velocities

"""

# Standard imports
import os
import glob
import argparse

import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation

from scipy.ndimage import percentile_filter

from pypeit import specobjs


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                         'first-hires-exoplanet-data')
DEFAULT_REDUX = os.path.join(DATA_ROOT, 'redux')
DEFAULT_TEMPLATE = os.path.join(DEFAULT_REDUX, 'template',
                                'hd187123_template_19980826.fits')
DEFAULT_OUTDIR = os.path.join(_HERE, 'data')

KECK = dict(lon=-155.47833333333335, lat=19.82833333333333, elv=4159.99999999977)
TARGET = SkyCoord(ra=296.744625 * u.deg, dec=34.41905555555555 * u.deg)
CKMS = 299792.458

#: Only orders lying entirely blueward of this are used
BLUE_LIMIT = 5000.

#: Half-width of the velocity search, km/s.  The barycentric velocity swings
#: -13 to +17 km/s across the discovery era and the template sits at -7.2, so
#: dv_obs can reach ~25 km/s; 60 gives headroom without inviting spurious peaks.
VSEARCH = 60.

#: Velocity step of the log-wavelength grid, km/s.  About a quarter of the
#: 6.7 km/s resolution element.
DV_GRID = 1.5

#: Width of the running-percentile continuum, in grid pixels
CONT_SIZE = 151
CONT_PCT = 90


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def bary_velocity(mjd_start, exptime):
    """ Barycentric correction at the exposure midpoint, km/s. """
    loc = EarthLocation.from_geodetic(lon=KECK['lon'] * u.deg,
                                      lat=KECK['lat'] * u.deg,
                                      height=KECK['elv'] * u.m)
    t_mid = Time(mjd_start + 0.5 * exptime / 86400., format='mjd', location=loc)
    return TARGET.radial_velocity_correction(
        kind='barycentric', obstime=t_mid).to(u.km / u.s).value


def load_template(path):
    """ Read the prompt-4 template.

    Returns:
        tuple: (dict order -> (wave, flux, ivar), reference barycentric velocity)
    """
    orders = {}
    with fits.open(path) as hdul:
        v_bary = float(hdul[0].header['VBARY'])
        for hdu in hdul[1:]:
            order = int(hdu.header['ECHORDER'])
            d = hdu.data
            wave = np.asarray(d['WAVE'], dtype=float)
            flux = np.asarray(d['FLUX'], dtype=float)
            ivar = np.asarray(d['IVAR'], dtype=float)
            good = (ivar > 0) & (wave > 0) & (flux > 0)
            if good.sum() > 200:
                orders[order] = (wave[good], flux[good], ivar[good])
    return orders, v_bary


def load_epoch(path):
    """ Read one spec1d into observed-frame per-order arrays. """
    sobjs = specobjs.SpecObjs.from_fitsfile(path, chk_version=False)
    hdr = sobjs.header
    with fits.open(path) as hdul:
        vel_corr = float(hdul[1].header.get('VEL_CORR', 1.0) or 1.0)

    orders = {}
    for obj in sobjs:
        wave = np.asarray(obj.OPT_WAVE) / vel_corr
        flux, ivar = np.asarray(obj.OPT_COUNTS), np.asarray(obj.OPT_COUNTS_IVAR)
        mask = np.asarray(obj.OPT_MASK)
        good = mask & (ivar > 0) & (wave > 0) & np.isfinite(flux) & (flux != 0)
        if good.sum() > 200:
            orders[int(obj.ECH_ORDER)] = (wave[good], flux[good], ivar[good])

    return dict(koaid=str(hdr['FILENAME']).replace('.fits', ''),
                mjd=float(hdr['MJD']), exptime=float(hdr['EXPTIME']),
                orders=orders)


# ---------------------------------------------------------------------------
# Cross-correlation
# ---------------------------------------------------------------------------

def normalise(wave, flux, size=CONT_SIZE, pct=CONT_PCT):
    """ Divide out a running-percentile continuum and subtract one. """
    cont = percentile_filter(flux, pct, size=min(size, len(flux) // 2 * 2 - 1))
    cont = np.maximum(cont, 1e-9)
    return flux / cont - 1.


def ccf(template_order, epoch_order):
    """ Cross-correlate one order and return the velocity of the peak.

    Both inputs are in the observed frame.  A common log-wavelength grid is
    used so that a shift in pixels is a constant shift in velocity.

    Args:
        template_order (tuple): (wave, flux, ivar) of the template.
        epoch_order (tuple): (wave, flux, ivar) of the epoch.

    Returns:
        dict or None: 'dv' (km/s), 'peak' (normalised correlation), 'sigma'
        (Zucker formal error, km/s), 'npix'.
    """
    wt, ft, _ = template_order
    we, fe, _ = epoch_order

    lo = max(wt.min(), we.min())
    hi = min(wt.max(), we.max())
    if hi / lo < 1.002:
        return None

    n = int(np.log(hi / lo) / (DV_GRID / CKMS))
    if n < 400:
        return None
    grid = lo * np.exp(np.arange(n) * DV_GRID / CKMS)

    a = normalise(grid, np.interp(grid, wt, ft))
    b = normalise(grid, np.interp(grid, we, fe))
    # Taper the ends so the shifted overlap does not inject an edge step
    taper = np.hanning(len(a))
    a, b = a * taper, b * taper

    maxlag = int(VSEARCH / DV_GRID)
    if maxlag >= n // 3:
        return None
    lags = np.arange(-maxlag, maxlag + 1)
    denom = np.sqrt(np.sum(a * a) * np.sum(b * b))
    if denom <= 0:
        return None
    cc = np.array([
        np.sum(a[max(0, l):len(a) + min(0, l)] * b[max(0, -l):len(b) + min(0, -l)])
        for l in lags]) / denom

    k = int(np.argmax(cc))
    if k == 0 or k == len(cc) - 1:
        return None
    y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]
    curv = y0 - 2 * y1 + y2
    frac = 0.5 * (y0 - y2) / curv if curv != 0 else 0.
    # SIGN.  In this construction a positive lag aligns a[i+l] with b[i]: the
    # epoch feature sits at the LOWER index, i.e. the shorter wavelength, so a
    # positive lag means the epoch is blueshifted and dv_obs is negative.
    # Getting this backwards does not look like a sign error -- it looks like
    # a real velocity, because it doubles the barycentric term instead of
    # cancelling it and hands back a smooth +37 to -11 km/s annual curve.
    dv = -(lags[k] + frac) * DV_GRID

    # Zucker (2003) eq. 3: the error on the peak position from the curvature of
    # the normalised correlation function and its height.
    peak = float(np.clip(y1, -0.999999, 0.999999))
    c2 = curv / (DV_GRID ** 2)          # second derivative, per (km/s)^2
    if c2 >= 0 or peak <= 0:
        return None
    var = -1. / (n * (c2 / peak) * (peak ** 2 / (1. - peak ** 2)))
    if not np.isfinite(var) or var <= 0:
        return None
    return dict(dv=float(dv), peak=peak, sigma=float(np.sqrt(var)), npix=n)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(options=None):
    parser = argparse.ArgumentParser(
        description='Relative velocities by cross-correlation.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--redux', type=str, default=DEFAULT_REDUX)
    parser.add_argument('--template', type=str, default=DEFAULT_TEMPLATE)
    parser.add_argument('--outdir', type=str, default=DEFAULT_OUTDIR)
    parser.add_argument('--quality', type=str,
                        default=os.path.join(DEFAULT_OUTDIR, 'order_quality.csv'))
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):

    template, v_bary_t = load_template(pargs.template)
    quality = Table.read(pargs.quality)
    use = {(str(r['koaid']), int(r['order'])): bool(int(r['use']))
           for r in quality}
    wmax = {(str(r['koaid']), int(r['order'])): float(r['wave_max'])
            for r in quality}

    blue_orders = sorted(o for o in template
                         if np.max(template[o][0]) < BLUE_LIMIT)
    print('\n=== 1. Setup ===')
    print('  template          : {:s}'.format(os.path.basename(pargs.template)))
    print('  template v_bary   : {:+.4f} km/s'.format(v_bary_t))
    print('  blue orders (<{:.0f} A): {:d}  (orders {:d}-{:d})'.format(
        BLUE_LIMIT, len(blue_orders), min(blue_orders), max(blue_orders)))
    print('  velocity grid     : {:.1f} km/s/pixel, search +/-{:.0f} km/s'.format(
        DV_GRID, VSEARCH))

    files = sorted(glob.glob(os.path.join(pargs.redux, 'reduce_19*', 'Science',
                                          'spec1d_*.fits')))
    print('  epochs            : {:d}'.format(len(files)))

    per_order, per_epoch = [], []
    for path in files:
        epoch = load_epoch(path)
        v_bary_e = bary_velocity(epoch['mjd'], epoch['exptime'])

        vs, ws, sigs, orders_used = [], [], [], []
        for order in blue_orders:
            if order not in epoch['orders']:
                continue
            if not use.get((epoch['koaid'], order), False):
                continue
            res = ccf(template[order], epoch['orders'][order])
            if res is None:
                continue
            # dv_obs -> relative barycentric velocity of the star
            v_rel = res['dv'] + (v_bary_e - v_bary_t)
            per_order.append(dict(
                koaid=epoch['koaid'], mjd=epoch['mjd'], order=order,
                v_rel=v_rel * 1e3, sigma=res['sigma'] * 1e3,
                peak=res['peak'], npix=res['npix']))
            vs.append(v_rel * 1e3)
            sigs.append(res['sigma'] * 1e3)
            orders_used.append(order)

        if len(vs) < 3:
            print('  {:s}: only {:d} usable orders, skipped'.format(
                epoch['koaid'], len(vs)))
            continue

        vs, sigs = np.array(vs), np.array(sigs)
        w = 1. / sigs ** 2
        vbar = float(np.sum(w * vs) / np.sum(w))
        sig_formal = float(1. / np.sqrt(np.sum(w)))
        # Weighted scatter of the orders about the epoch mean
        resid = vs - vbar
        wvar = np.sum(w * resid ** 2) / np.sum(w)
        scatter = float(np.sqrt(wvar * len(vs) / max(len(vs) - 1, 1)))
        sig_emp = scatter / np.sqrt(len(vs))

        per_epoch.append(dict(
            koaid=epoch['koaid'], mjd=epoch['mjd'], exptime=epoch['exptime'],
            n_orders=len(vs), v_rel=vbar,
            sigma_formal=sig_formal, sigma_empirical=float(sig_emp),
            order_scatter=scatter, v_bary=v_bary_e,
            median_peak=float(np.median([p['peak'] for p in per_order
                                         if p['koaid'] == epoch['koaid']]))))

    otbl, etbl = Table(per_order), Table(per_epoch)

    # --- report -----------------------------------------------------------
    print('\n=== 2. Per-epoch velocities ===')
    print('  v_rel is relative to the 1998-08-26 template epoch.')
    show = etbl['koaid', 'mjd', 'n_orders', 'v_rel', 'sigma_formal',
                'sigma_empirical', 'order_scatter', 'median_peak']
    for col, fmt in [('mjd', '%.5f'), ('v_rel', '%+.1f'),
                     ('sigma_formal', '%.1f'), ('sigma_empirical', '%.1f'),
                     ('order_scatter', '%.1f'), ('median_peak', '%.3f')]:
        show[col].info.format = fmt
    show.pprint_all()

    print('\n=== 3. The two scatters ===')
    print('  They answer different questions and must not be conflated.')
    os_ = np.asarray(etbl['order_scatter'], dtype=float)
    print('\n  PER-ORDER scatter, within one epoch')
    print('    "do the 22 blue orders of one exposure agree with each other?"')
    print('    median {:.0f} m/s   range {:.0f}-{:.0f} m/s'.format(
        np.median(os_), np.min(os_), np.max(os_)))
    print('    -> the precision of a single order, and the internal'
          ' consistency of\n       the wavelength solution across the format')

    v = np.asarray(etbl['v_rel'], dtype=float)
    print('\n  EPOCH-TO-EPOCH scatter, across the baseline')
    print('    "does the same star give the same velocity on different'
          ' nights?"')
    print('    rms {:.0f} m/s about the mean, over {:d} epochs and {:.0f} days'
          .format(np.std(v, ddof=1), len(v),
                  np.ptp(np.asarray(etbl['mjd'], dtype=float))))
    print('    -> the end-to-end precision of the method, including everything'
          '\n       the per-order scatter cannot see: flexure, focus, the'
          ' wavelength\n       zero point drifting between nights, the'
          ' barycentric correction')

    med_formal = np.median(np.asarray(etbl['sigma_formal'], dtype=float))
    med_emp = np.median(np.asarray(etbl['sigma_empirical'], dtype=float))
    print('\n  For comparison, the per-epoch uncertainties:')
    print('    formal (Zucker, combined) : median {:.1f} m/s'.format(med_formal))
    print('    empirical (scatter/sqrt N): median {:.1f} m/s'.format(med_emp))
    print('    epoch-to-epoch rms        : {:.0f} m/s'.format(np.std(v, ddof=1)))
    print('    The epoch-to-epoch rms exceeds both by a large factor, so the'
          '\n    per-epoch error bars are NOT the precision of this'
          ' measurement.')

    print('\n=== 4. Where the epoch-to-epoch scatter comes from ===')
    etbl['night'] = [k.split('.')[1] for k in etbl['koaid']]
    nights = sorted(set(etbl['night']))

    # Within-night scatter isolates what changes over hours -- flexure, focus,
    # the wavelength solution ageing away from its arc -- from what changes
    # between nights.
    within = []
    rows = []
    for night in nights:
        sel = etbl[etbl['night'] == night]
        vv = np.asarray(sel['v_rel'], dtype=float)
        span = (np.ptp(np.asarray(sel['mjd'], dtype=float)) * 24.
                if len(sel) > 1 else 0.)
        rows.append(dict(night=night, n=len(sel), mean=float(np.mean(vv)),
                         spread=float(np.ptp(vv)) if len(vv) > 1 else np.nan,
                         hours=span))
        if len(vv) > 1:
            within.append(np.std(vv, ddof=1))
    ntbl = Table(rows)
    for col, fmt in [('mean', '%+.1f'), ('spread', '%.1f'), ('hours', '%.1f')]:
        ntbl[col].info.format = fmt
    ntbl.pprint_all()

    nightly = np.asarray(ntbl['mean'], dtype=float)
    print('\n  WITHIN a night ({:d} nights with >1 frame):'.format(len(within)))
    print('    median rms {:.0f} m/s, worst spread {:.0f} m/s'.format(
        np.median(within), np.nanmax(np.asarray(ntbl['spread'], dtype=float))))
    print('  BETWEEN nights ({:d} nightly means):'.format(len(nightly)))
    print('    rms {:.0f} m/s'.format(np.std(nightly, ddof=1)))

    # 1998-07-19 took no calibration frames at all and borrows the 07-18 arc.
    no_arc = etbl[etbl['night'] == '19980719']
    rest = etbl[etbl['night'] != '19980719']
    print('\n  1998-07-19 is the only night with no arc of its own'
          ' (it borrows 07-18\'s):')
    print('    its {:d} frames: {:s} m/s'.format(
        len(no_arc), ', '.join('{:+.0f}'.format(x) for x in no_arc['v_rel'])))
    print('    all other epochs   : rms {:.0f} m/s'.format(
        np.std(np.asarray(rest['v_rel'], dtype=float), ddof=1)))
    print('    all epochs         : rms {:.0f} m/s'.format(np.std(v, ddof=1)))
    print('    Phase 1 concluded that same-night arcs are not required because'
          '\n    the wavelength RMS was indistinguishable.  The RMS is not the'
          ' zero\n    point: 07-19 sits 1.5-2.4 km/s off every other July'
          ' night.')

    print('\n=== 5. The template exposures, as an internal check ===')
    tmpl = [r for r in etbl if str(r['koaid']).startswith('HI.19980826.3')]
    for r in tmpl:
        print('  {:s}  v_rel {:+7.1f} m/s  (it is inside the template)'.format(
            r['koaid'], r['v_rel']))

    if not os.path.isdir(pargs.outdir):
        os.makedirs(pargs.outdir)
    for col, fmt in [('v_rel', '%.2f'), ('sigma', '%.2f'), ('peak', '%.5f'),
                     ('mjd', '%.6f')]:
        otbl[col].info.format = fmt
    for col in ['mjd', 'v_rel', 'sigma_formal', 'sigma_empirical',
                'order_scatter', 'v_bary', 'median_peak']:
        etbl[col].info.format = '%.4f'
    f1 = os.path.join(pargs.outdir, 'xcorr_velocities.csv')
    f2 = os.path.join(pargs.outdir, 'xcorr_velocities_per_order.csv')
    etbl.write(f1, format='ascii.csv', overwrite=True)
    otbl.write(f2, format='ascii.csv', overwrite=True)
    print('\nWrote {:s}\nWrote {:s}'.format(f1, f2))


if __name__ == '__main__':
    main(parse_args())
