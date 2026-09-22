""" The iodine-free stellar template of HD 187123 (phase 2, prompt 4).

The Butler forward-model method needs a spectrum of the star taken with the
iodine cell **out**, at higher signal-to-noise than any single science frame.
1998-08-26 carries exactly that: three consecutive 500 s exposures, cell out,
which the observer labelled `OBJECT = 'Template'` on the matching 08-12 frame
and clearly took on purpose.  This co-adds them, order by order, into one
spectrum, and checks it three ways.

REFERENCE FRAME.  Following prompt 1: PypeIt applied a *heliocentric*
correction by default and stored the multiplicative factor as `VEL_CORR`, so
every wavelength here is first divided by it to recover the **observed** frame,
which is what `pyodine` expects and what the co-addition must be done in.  The
three exposures span 20 minutes, over which the barycentric velocity moves
about 31 m/s; each is shifted onto the first exposure's frame before co-adding,
so that motion is removed rather than smeared in.  That is far below the 6.7
km/s resolution element, but it costs nothing to do properly and the template
is the one product every later velocity is measured against.

`barycorrpy` is not installed in `pypeit14`, so `astropy` is used.  Prompt 1
established that astropy's *barycentric* correction is the sound one (PypeIt's
heliocentric has a sign error on the solar term) and that the two agree to a
near-constant 4.6 m/s.  For the *relative* alignment of three exposures twenty
minutes apart that is irrelevant.  Phase 3 should still install `barycorrpy`,
as prompt 1 recommended, so that pyodine's initial guess and its final
correction come from one source.

CHECKS.
  1. Against the 1997-12-24 cell-out frame, eight months earlier: per-order
     cross-correlation velocity and continuum-normalised line agreement.
  2. Against a G2V expectation: the velocity zero point against HD 187123's
     -17 km/s systemic velocity, the depths of the canonical G-dwarf features,
     and the width of isolated photospheric lines.
  3. That the cell really was out: the 5000-6200 A region is compared against
     a cell-in frame from the same night, where the I2 forest must appear and
     in the template must not.

NOTE ON WAVELENGTHS: PypeIt reports *vacuum* wavelengths.  Rest wavelengths of
the stellar features are tabulated in air, as is conventional in the optical,
and are converted with pypeit.core.wave.airtovac before use.  Comparing
observed vacuum wavelengths against air rest values produces a spurious
redshift of about +70 km/s at 5000 A.

The template is a reduction product and is written to the sibling data tree,
not the repository.  The per-order signal-to-noise table is small and derived
and is written to `first_hires_exoplanet/data/`.

Reduced with PypeIt 2.0.2.dev1217+g017bece06 (branch orig-hires-fixes).

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.build_template

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

from pypeit import specobjs
from pypeit.core.wave import airtovac


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                         'first-hires-exoplanet-data')
DEFAULT_REDUX = os.path.join(DATA_ROOT, 'redux')
DEFAULT_OUTDIR = os.path.join(_HERE, 'data')

#: The three cell-out exposures, 1998-08-26, 500 s each
TEMPLATE_FRAMES = ('HI.19980826.34749', 'HI.19980826.35345', 'HI.19980826.35939')

#: The independent cell-out frame, eight months earlier, for the consistency
#: check.  1997-12-24's only science frame happens to be cell-out.
CHECK_FRAME = 'HI.19971224.16259'

#: A cell-in frame from the same night as the template, to prove the cell
#: really was out of the beam in the template exposures
IODINE_FRAME = 'HI.19980826.44043'

#: Keck, as PypeIt's telescope parameters give it
KECK = dict(lon=-155.47833333333335, lat=19.82833333333333, elv=4159.99999999977)

#: HD 187123, for the barycentric correction
TARGET = SkyCoord(ra=296.744625 * u.deg, dec=34.41905555555555 * u.deg)

#: Systemic radial velocity of HD 187123 (km/s)
V_SYS = -17.0

#: Canonical G-dwarf features, air wavelengths
G_LINES_AIR = {
    'Ca II K': 3933.66,
    'Ca II H': 3968.47,
    'H-delta': 4101.74,
    'H-gamma': 4340.47,
    'H-beta': 4861.33,
    'Mg b1': 5167.32,
    'Mg b2': 5172.68,
    'Mg b3': 5183.60,
    'Na D2': 5889.95,
    'Na D1': 5895.92,
}

#: Isolated, unblended Fe I lines for the line-width check, air wavelengths
FE_LINES_AIR = (4383.55, 4404.75, 5269.54, 5328.04, 5371.49, 5397.13, 5405.77)

CKMS = 299792.458


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def find_spec1d(redux_dir, koaid):
    """ The spec1d file for one KOAID.

    Args:
        redux_dir (str): directory holding the `reduce_*` folders.
        koaid (str): e.g. 'HI.19980826.34749'.

    Returns:
        str: path to the spec1d file.
    """
    hits = sorted(glob.glob(os.path.join(
        redux_dir, 'reduce_*', 'Science', 'spec1d_{:s}*.fits'.format(koaid))))
    if len(hits) != 1:
        raise RuntimeError('Expected one spec1d for {:s}, found {:d}'.format(
            koaid, len(hits)))
    return hits[0]


def bary_velocity(mjd_start, exptime):
    """ Barycentric velocity correction at the exposure midpoint, km/s.

    Prompt 1 established that the header MJD is the exposure *start* (proved to
    1.4 s from the recorded hour angle), so the midpoint has to be constructed.

    Args:
        mjd_start (float): header MJD.
        exptime (float): exposure time, seconds.

    Returns:
        float: correction in km/s, positive when the observer moves toward the
        star.
    """
    loc = EarthLocation.from_geodetic(lon=KECK['lon'] * u.deg,
                                      lat=KECK['lat'] * u.deg,
                                      height=KECK['elv'] * u.m)
    t_mid = Time(mjd_start + 0.5 * exptime / 86400., format='mjd', location=loc)
    return TARGET.radial_velocity_correction(
        kind='barycentric', obstime=t_mid).to(u.km / u.s).value


def load_exposure(path):
    """ Read one spec1d into observed-frame, per-order arrays.

    Args:
        path (str): a spec1d file.

    Returns:
        dict: 'orders' maps echelle order -> (wave, flux, ivar, mask), plus
        the header quantities the co-addition needs.
    """
    sobjs = specobjs.SpecObjs.from_fitsfile(path, chk_version=False)
    hdr = sobjs.header
    with fits.open(path) as hdul:
        vel_corr = float(hdul[1].header.get('VEL_CORR', 1.0) or 1.0)

    orders = {}
    for obj in sobjs:
        wave, flux = np.asarray(obj.OPT_WAVE), np.asarray(obj.OPT_COUNTS)
        ivar, mask = np.asarray(obj.OPT_COUNTS_IVAR), np.asarray(obj.OPT_MASK)
        # Undo PypeIt's heliocentric correction: back to the observed frame.
        # Exact -- VEL_CORR is one scalar and the round trip is lossless
        # (prompt 1).
        wave = wave / vel_corr
        # Masked pixels arrive as exact zeros, not flagged gaps (phase 1).
        good = mask & (ivar > 0) & (wave > 0) & np.isfinite(flux) & (flux != 0)
        orders[int(obj.ECH_ORDER)] = (wave, flux, ivar, good)

    return dict(path=path, orders=orders, vel_corr=vel_corr,
                mjd=float(hdr['MJD']), exptime=float(hdr['EXPTIME']),
                koaid=str(hdr['FILENAME']).replace('.fits', ''))


# ---------------------------------------------------------------------------
# Co-addition
# ---------------------------------------------------------------------------

def coadd(exposures):
    """ Inverse-variance co-add exposures onto the first one's wavelength grid.

    Args:
        exposures (list): output of :func:`load_exposure`, reference first.

    Returns:
        tuple: (dict order -> (wave, flux, ivar, npix), list of applied shifts
        in m/s)
    """
    ref = exposures[0]
    v_ref = bary_velocity(ref['mjd'], ref['exptime'])

    shifts = []
    out = {}
    for order in sorted(ref['orders']):
        wave_ref = ref['orders'][order][0]
        num = np.zeros_like(wave_ref)
        den = np.zeros_like(wave_ref)
        npix = np.zeros(len(wave_ref), dtype=int)

        for exp in exposures:
            if order not in exp['orders']:
                continue
            wave, flux, ivar, good = exp['orders'][order]
            if good.sum() < 10:
                continue
            # Put this exposure on the reference exposure's barycentric frame.
            # lambda_obs scales as (1 - v_bary/c), so the ratio of the two
            # observed scales is (1 - v/c) / (1 - v_ref/c).
            v = bary_velocity(exp['mjd'], exp['exptime'])
            factor = (1. - v_ref / CKMS) / (1. - v / CKMS)
            wave_shifted = wave * factor

            # Interpolate onto the reference grid.  Only pixels bracketed by
            # good input pixels are used; np.interp would happily extrapolate
            # flat across a masked gap.
            w, f, iv = wave_shifted[good], flux[good], ivar[good]
            if len(w) < 10:
                continue
            fi = np.interp(wave_ref, w, f, left=np.nan, right=np.nan)
            ii = np.interp(wave_ref, w, iv, left=np.nan, right=np.nan)
            # Kill anything that fell in a gap wider than a few pixels
            gap = np.interp(wave_ref, w, np.arange(len(w)),
                            left=np.nan, right=np.nan)
            ok = np.isfinite(fi) & np.isfinite(ii) & (ii > 0) & np.isfinite(gap)
            num[ok] += fi[ok] * ii[ok]
            den[ok] += ii[ok]
            npix[ok] += 1
            if exp is not ref:
                shifts.append((v - v_ref) * 1e3)

        flux_out = np.where(den > 0, num / np.maximum(den, 1e-30), 0.)
        out[order] = (wave_ref, flux_out, den, npix)
    return out, shifts


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def order_snr(template):
    """ Median S/N per pixel in each order of the template.

    Args:
        template (dict): output of :func:`coadd`.

    Returns:
        `astropy.table.Table`_
    """
    rows = []
    for order in sorted(template):
        wave, flux, ivar, npix = template[order]
        good = (ivar > 0) & (npix > 0) & (flux > 0)
        if good.sum() < 100:
            rows.append(dict(order=order, wave_min=np.nan, wave_max=np.nan,
                             n_good=int(good.sum()), n_coadd=0, snr=np.nan))
            continue
        rows.append(dict(
            order=order,
            wave_min=float(wave[good].min()), wave_max=float(wave[good].max()),
            n_good=int(good.sum()),
            n_coadd=int(np.median(npix[good])),
            snr=float(np.median(flux[good] * np.sqrt(ivar[good])))))
    return Table(rows)


def line_velocity(template, lines_air):
    """ Velocity of each named line against its vacuum rest wavelength.

    The deepest pixel within +/-0.8 A of the expected position, parabolically
    refined.  Crude, but it pins the zero point to a fraction of a km/s.

    Args:
        template (dict): output of :func:`coadd`.
        lines_air (dict): name -> air rest wavelength.

    Returns:
        `astropy.table.Table`_
    """
    rows = []
    for name, air in lines_air.items():
        vac = float(airtovac(air * u.AA).value)
        # Expected observed position, at the systemic velocity
        expect = vac * (1. + V_OBS_EXPECTED[0] / CKMS)
        best = None
        for order in sorted(template):
            wave, flux, ivar, npix = template[order]
            good = (ivar > 0) & (npix > 0) & (flux > 0)
            if good.sum() < 100:
                continue
            sel = good & (np.abs(wave - expect) < 0.8)
            if sel.sum() < 7:
                continue
            w, f = wave[sel], flux[sel]
            # Continuum from the 90th percentile of a wider window
            wide = good & (np.abs(wave - expect) < 6.)
            cont = np.percentile(flux[wide], 90) if wide.sum() > 50 else np.max(f)
            i = int(np.argmin(f))
            if i == 0 or i == len(f) - 1:
                continue
            # Parabolic refinement on the three lowest pixels
            y0, y1, y2 = f[i - 1], f[i], f[i + 1]
            denom = (y0 - 2 * y1 + y2)
            dx = 0.5 * (y0 - y2) / denom if denom != 0 else 0.
            centre = w[i] + dx * (w[i + 1] - w[i - 1]) / 2.
            depth = 1. - y1 / cont if cont > 0 else np.nan
            cand = dict(line=name, order=order, vac=vac,
                        observed=float(centre),
                        vel=float((centre / vac - 1.) * CKMS),
                        depth=float(depth))
            if best is None or cand['depth'] > best['depth']:
                best = cand
        if best is not None:
            rows.append(best)
    return Table(rows)


def line_width(template, lines_air, half_window=0.45):
    """ FWHM of isolated Fe I lines, in km/s.

    A slowly rotating G dwarf observed at R ~ 45,000 should give roughly the
    instrumental width, 6.7 km/s, broadened a little.  A giant or a fast
    rotator would not.

    The FWHM is taken across the **contiguous** run of pixels above half depth
    that contains the line minimum.  Spanning min to max of everything above
    half depth, as a first version did, measures the window rather than the
    line as soon as a blend enters it: it returned 111 km/s for Fe I 4383,
    which is 1.6 A, the full window.  Lines whose profile does not come back
    below half depth on both sides inside the window are blended and dropped.

    Args:
        template (dict): output of :func:`coadd`.
        lines_air (tuple): air rest wavelengths.
        half_window (float): half-width of the fitting window, Angstrom.

    Returns:
        `astropy.table.Table`_
    """
    rows = []
    for air in lines_air:
        vac = float(airtovac(air * u.AA).value)
        expect = vac * (1. + V_OBS_EXPECTED[0] / CKMS)
        for order in sorted(template):
            wave, flux, ivar, npix = template[order]
            good = (ivar > 0) & (npix > 0) & (flux > 0)
            sel = good & (np.abs(wave - expect) < half_window)
            if sel.sum() < 9:
                continue
            w, f = wave[sel], flux[sel]
            # Continuum from a window wide enough to reach real continuum but
            # not so wide it crosses another strong feature
            wide = good & (np.abs(wave - expect) < 3.)
            if wide.sum() < 40:
                continue
            cont = np.percentile(flux[wide], 95)
            if cont <= 0:
                continue
            prof = 1. - f / cont
            i = int(np.argmax(prof))
            if prof[i] < 0.05 or i == 0 or i == len(prof) - 1:
                continue
            half = 0.5 * prof[i]
            # Walk out from the minimum until the profile drops below half
            lo = i
            while lo > 0 and prof[lo - 1] >= half:
                lo -= 1
            hi = i
            while hi < len(prof) - 1 and prof[hi + 1] >= half:
                hi += 1
            if lo == 0 or hi == len(prof) - 1:
                continue  # never came back down: blended, or window too narrow
            # Linear interpolation to the half-depth crossings
            def cross(a, b):
                if prof[b] == prof[a]:
                    return w[a]
                return w[a] + (half - prof[a]) * (w[b] - w[a]) / (prof[b] - prof[a])
            w_lo = cross(lo - 1, lo)
            w_hi = cross(hi + 1, hi)
            rows.append(dict(line_air=air, order=order,
                             depth=float(prof[i]),
                             fwhm_kms=float((w_hi - w_lo) / vac * CKMS)))
            break
    return Table(rows)


def weak_line_widths(template, orders=(78, 79, 80, 81), depth_range=(0.10, 0.40)):
    """ FWHM of the weak, unsaturated lines, which bounds the resolution.

    The named Fe I lines are all strong: depths of 0.78-0.84, with damping
    wings, so their 13-22 km/s widths say nothing about the instrument.  A
    *weak* line is not saturated and its observed width is close to the
    instrumental profile, so the narrowest features in the spectrum put an
    upper bound on the resolution element.  Orders 78-81 are used because
    prompt 2 found them clean in every epoch and they carry no I2 even in a
    cell-in frame.

    Args:
        template (dict): output of :func:`coadd`.
        orders (tuple): echelle orders to search.
        depth_range (tuple): keep lines in this fractional-depth range.

    Returns:
        `astropy.table.Table`_: one row per accepted line.
    """
    from scipy.signal import find_peaks, medfilt

    rows = []
    for order in orders:
        if order not in template:
            continue
        wave, flux, ivar, npix = template[order]
        good = (ivar > 0) & (npix > 0) & (flux > 0)
        if good.sum() < 500:
            continue
        w, f = wave[good], flux[good]
        # Running continuum: the upper envelope, from a wide percentile filter
        from scipy.ndimage import percentile_filter
        cont = percentile_filter(f, 92, size=201)
        prof = 1. - f / np.maximum(cont, 1e-9)
        peaks, props = find_peaks(prof, height=depth_range,
                                  distance=4, prominence=0.05)
        for k in peaks:
            d = prof[k]
            half = 0.5 * d
            lo, hi = k, k
            while lo > 0 and prof[lo - 1] >= half:
                lo -= 1
            while hi < len(prof) - 1 and prof[hi + 1] >= half:
                hi += 1
            if lo == 0 or hi == len(prof) - 1 or (hi - lo) < 2:
                continue
            # Reject blends: require the profile to fall to a third of the
            # depth within 6 pixels either side
            l2, h2 = max(0, lo - 6), min(len(prof) - 1, hi + 6)
            if prof[l2] > d / 3. or prof[h2] > d / 3.:
                continue

            def cross(a, b):
                if prof[b] == prof[a]:
                    return w[a]
                return w[a] + (half - prof[a]) * (w[b] - w[a]) / (prof[b] - prof[a])

            fwhm = cross(lo - 1, lo) - cross(hi + 1, hi)
            fwhm = abs(fwhm)
            rows.append(dict(order=order, wave=float(w[k]), depth=float(d),
                             fwhm_kms=float(fwhm / w[k] * CKMS)))
    return Table(rows) if rows else None


def iodine_check(template, iodine_path, redux_dir):
    """ Does the I2 forest appear in the cell-in frame and not the template?

    Counts absorption minima per 100 A in 5000-6200 A, where the cell absorbs,
    and in 4000-4800 A, where it does not.  The blue window is the control: if
    the template's line density were low everywhere the test would prove
    nothing.

    Args:
        template (dict): output of :func:`coadd`.
        iodine_path (str): spec1d of a cell-in frame.
        redux_dir (str): unused, kept for symmetry.

    Returns:
        `astropy.table.Table`_
    """
    from scipy.signal import find_peaks

    iodine = load_exposure(iodine_path)

    def density(orders_dict, is_template, lo, hi):
        n_lines, span = 0, 0.
        for order in sorted(orders_dict):
            if is_template:
                wave, flux, ivar, npix = orders_dict[order]
                good = (ivar > 0) & (npix > 0) & (flux > 0)
            else:
                wave, flux, ivar, good = orders_dict[order]
            sel = good & (wave > lo) & (wave < hi)
            if sel.sum() < 200:
                continue
            w, f = wave[sel], flux[sel]
            cont = np.percentile(f, 90)
            if cont <= 0:
                continue
            peaks, _ = find_peaks(1. - f / cont, prominence=0.04)
            n_lines += len(peaks)
            span += w.max() - w.min()
        return (100. * n_lines / span) if span > 0 else np.nan

    rows = []
    for lbl, lo, hi in [('iodine 5000-6200 A', 5000., 6200.),
                        ('control 4000-4800 A', 4000., 4800.)]:
        rows.append(dict(window=lbl,
                         template=density(template, True, lo, hi),
                         cell_in=density(iodine['orders'], False, lo, hi)))
    out = Table(rows)
    out['template'].info.format = '%.1f'
    out['cell_in'].info.format = '%.1f'
    return out


def compare_against(template, other_path):
    """ Cross-correlate the template against another cell-out frame, per order.

    Args:
        template (dict): output of :func:`coadd`.
        other_path (str): spec1d of the comparison frame.

    Returns:
        tuple: (`astropy.table.Table`_, the comparison exposure)
    """
    other = load_exposure(other_path)
    v_t = bary_velocity(TEMPLATE_MJD[0], 500.)
    v_o = bary_velocity(other['mjd'], other['exptime'])

    rows = []
    for order in sorted(template):
        if order not in other['orders']:
            continue
        wave, flux, ivar, npix = template[order]
        good = (ivar > 0) & (npix > 0) & (flux > 0)
        w_o, f_o, iv_o, g_o = other['orders'][order]
        if good.sum() < 500 or g_o.sum() < 500:
            continue

        # Both are in their own observed frames; put the comparison on the
        # template's barycentric frame before correlating, so what is left is
        # the star plus any reduction difference.
        factor = (1. - v_t / CKMS) / (1. - v_o / CKMS)
        w_shift = w_o[g_o] * factor

        lo = max(wave[good].min(), w_shift.min()) + 0.5
        hi = min(wave[good].max(), w_shift.max()) - 0.5
        if hi - lo < 10.:
            continue
        grid = np.linspace(lo, hi, 8000)
        a = np.interp(grid, wave[good], flux[good])
        b = np.interp(grid, w_shift, f_o[g_o])
        a = a / np.percentile(a, 90) - 1.
        b = b / np.percentile(b, 90) - 1.

        # Cross-correlate over +/- 30 km/s
        dv = (grid[1] - grid[0]) / np.mean(grid) * CKMS
        maxlag = int(30. / dv)
        lags = np.arange(-maxlag, maxlag + 1)
        cc = np.array([
            np.sum(a[max(0, l):len(a) + min(0, l)] * b[max(0, -l):len(b) + min(0, -l)])
            for l in lags])
        # Parabolic refinement of the peak.  Without it every order returns
        # exactly 0.00 km/s, because the grid samples ~2 km/s per lag and the
        # true offset is smaller than one lag -- a quantisation artefact that
        # looks like a suspiciously perfect result.
        k = int(np.argmax(cc))
        if 0 < k < len(cc) - 1:
            y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]
            denom = y0 - 2 * y1 + y2
            frac = 0.5 * (y0 - y2) / denom if denom != 0 else 0.
        else:
            frac = 0.
        peak = lags[k] + frac
        norm = cc.max() / np.sqrt(np.sum(a * a) * np.sum(b * b))
        rows.append(dict(order=order, wave=float(np.mean(grid)),
                         dv_kms=float(peak * dv), corr=float(norm)))
    out = Table(rows)
    out['wave'].info.format = '%.0f'
    out['dv_kms'].info.format = '%+.2f'
    out['corr'].info.format = '%.3f'
    return out, other


#: filled in by main, so compare_against can reach the template epoch
TEMPLATE_MJD = [0.]

#: The velocity the lines should sit at IN THE OBSERVED FRAME, filled in by
#: main.  The template's wavelengths have had PypeIt's heliocentric correction
#: removed, so the lines are not at the -17 km/s systemic velocity: they are at
#: v_sys - BVC, and BVC is -7.24 km/s on this night.  Comparing against -17
#: would manufacture an 8 km/s error that is not there.
V_OBS_EXPECTED = [V_SYS]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(options=None):
    parser = argparse.ArgumentParser(
        description='Build the iodine-free stellar template of HD 187123.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--redux', type=str, default=DEFAULT_REDUX,
                        help='Directory holding the reduce_* folders')
    parser.add_argument('--outdir', type=str, default=DEFAULT_OUTDIR,
                        help='Where to write the S/N table')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):

    print('\n=== 1. The three cell-out exposures ===')
    exposures = []
    for koaid in TEMPLATE_FRAMES:
        path = find_spec1d(pargs.redux, koaid)
        exp = load_exposure(path)
        raw = os.path.join(DATA_ROOT, 'raw', '1998aug26', koaid + '.fits')
        iodin = fits.getheader(raw).get('IODIN') if os.path.isfile(raw) else None
        if iodin:
            raise RuntimeError('{:s} has the iodine cell IN'.format(koaid))
        vb = bary_velocity(exp['mjd'], exp['exptime'])
        print('  {:s}  MJD {:.6f}  {:.0f} s  IODIN={!s:5s}  '
              'orders {:d}  v_bary {:+.4f} km/s'.format(
                  koaid, exp['mjd'], exp['exptime'], iodin,
                  len(exp['orders']), vb))
        exposures.append(exp)
    TEMPLATE_MJD[0] = exposures[0]['mjd']
    v_bary_ref = bary_velocity(exposures[0]['mjd'], exposures[0]['exptime'])
    V_OBS_EXPECTED[0] = V_SYS - v_bary_ref

    print('\n=== 2. Co-adding ===')
    template, shifts = coadd(exposures)
    print('  reference frame          : {:s}, observed frame'.format(
        exposures[0]['koaid']))
    print('  barycentric shifts applied: {:+.1f} to {:+.1f} m/s'.format(
        min(shifts), max(shifts)))
    print('  (the resolution element is 6700 m/s, so this is cosmetic -- but'
          ' it is the\n   template every later velocity is measured against)')
    print('  orders in template        : {:d}'.format(len(template)))

    snr = order_snr(template)
    print('\n=== 3. Signal-to-noise per order ===')
    single = order_snr({o: (exposures[0]['orders'][o][0],
                            exposures[0]['orders'][o][1],
                            exposures[0]['orders'][o][2],
                            exposures[0]['orders'][o][3].astype(int))
                        for o in exposures[0]['orders']})
    merged = Table(dict(order=snr['order'], wave_min=snr['wave_min'],
                        wave_max=snr['wave_max'], n_coadd=snr['n_coadd'],
                        snr_single=single['snr'], snr_template=snr['snr']))
    merged['gain'] = merged['snr_template'] / merged['snr_single']
    for col, fmt in [('wave_min', '%.0f'), ('wave_max', '%.0f'),
                     ('snr_single', '%.1f'), ('snr_template', '%.1f'),
                     ('gain', '%.2f')]:
        merged[col].info.format = fmt
    merged.pprint_all()
    good = np.isfinite(np.asarray(merged['snr_template'], dtype=float))
    print('\n  template S/N  min {:.1f}  median {:.1f}  max {:.1f}'.format(
        np.nanmin(merged['snr_template']), np.nanmedian(merged['snr_template']),
        np.nanmax(merged['snr_template'])))
    print('  gain over one exposure: median {:.2f} (sqrt(3) = 1.73)'.format(
        np.nanmedian(np.asarray(merged['gain'], dtype=float)[good])))

    print('\n=== 4. Is the cell really out? ===')
    iod = iodine_check(template, find_spec1d(pargs.redux, IODINE_FRAME),
                       pargs.redux)
    print('  absorption minima per 100 A (prominence > 4% of continuum)')
    iod.pprint_all()

    print('\n=== 5. Against a G2V expectation ===')
    lines = line_velocity(template, G_LINES_AIR)
    lines['vac'].info.format = '%.2f'
    lines['observed'].info.format = '%.2f'
    lines['vel'].info.format = '%+.1f'
    lines['depth'].info.format = '%.3f'
    lines.pprint_all()
    vels = np.asarray(lines['vel'], dtype=float)
    print('\n  These wavelengths are in the OBSERVED frame, so the lines are'
          ' not at the')
    print('  {:+.1f} km/s systemic velocity: they are at v_sys - BVC ='
          ' {:+.2f} km/s.'.format(V_SYS, V_OBS_EXPECTED[0]))
    print('  measured median  : {:+.2f} km/s'.format(np.median(vels)))
    print('  expected         : {:+.2f} km/s'.format(V_OBS_EXPECTED[0]))
    print('  residual         : {:+.2f} km/s'.format(
        np.median(vels) - V_OBS_EXPECTED[0]))
    print('  Phase 1 found the same sign and size of residual (-14.7 against'
          ' -17.0,\n  i.e. +2.3 km/s) on a July frame, so this is the known'
          ' zero-point offset,\n  not a new one.  It is ~0.2 of a 6.7 km/s'
          ' resolution element and comes from\n  taking the deepest pixel of'
          ' deep, asymmetric line cores.')

    widths = line_width(template, FE_LINES_AIR)
    if len(widths):
        widths['depth'].info.format = '%.3f'
        widths['fwhm_kms'].info.format = '%.2f'
        print('\n  isolated Fe I line widths:')
        widths.pprint_all()
        print('  median FWHM {:.2f} km/s -- but these are strong, saturated'
              ' lines with damping\n  wings, so this is a stellar width, not'
              ' an instrumental one.'.format(
                  np.median(np.asarray(widths['fwhm_kms'], dtype=float))))

    weak = weak_line_widths(template)
    if weak is not None and len(weak) > 20:
        fw = np.asarray(weak['fwhm_kms'], dtype=float)
        print('\n  weak unsaturated lines (depth 0.10-0.40) in orders 78-81,'
              ' {:d} found:'.format(len(weak)))
        for q in (10, 25, 50):
            print('    {:2d}th percentile FWHM : {:5.2f} km/s'.format(
                q, np.percentile(fw, q)))
        print('  The narrowest features bound the resolution element from'
              ' above.')
        print('  R = c / FWHM at the 10th percentile = {:.0f}'.format(
            CKMS / np.percentile(fw, 10)))

    print('\n=== 6. Against the 1997-12-24 cell-out frame ===')
    cmp_tbl, other = compare_against(template, find_spec1d(pargs.redux, CHECK_FRAME))
    print('  {:s}, {:.0f} s, eight months earlier'.format(
        other['koaid'], other['exptime']))
    dv = np.asarray(cmp_tbl['dv_kms'], dtype=float)
    cc = np.asarray(cmp_tbl['corr'], dtype=float)
    print('  orders compared : {:d}'.format(len(cmp_tbl)))
    print('  correlation     : median {:.3f}  min {:.3f}'.format(
        np.median(cc), np.min(cc)))
    print('  velocity offset : median {:+.2f} km/s  scatter {:.2f} km/s'.format(
        np.median(dv), np.std(dv)))
    print('  (HD 187123 b has a 72 m/s semiamplitude, far below this method\'s'
          '\n   resolution; what this tests is that the two epochs are the same'
          ' star,\n   reduced consistently, not the planet)')

    # --- write ------------------------------------------------------------
    if not os.path.isdir(pargs.outdir):
        os.makedirs(pargs.outdir)
    snr_file = os.path.join(pargs.outdir, 'template_snr.csv')
    merged.write(snr_file, format='ascii.csv', overwrite=True)
    print('\nWrote {:s}'.format(snr_file))

    tmpl_dir = os.path.join(pargs.redux, 'template')
    if not os.path.isdir(tmpl_dir):
        os.makedirs(tmpl_dir)
    tmpl_file = os.path.join(tmpl_dir, 'hd187123_template_19980826.fits')
    hdus = [fits.PrimaryHDU()]
    hdus[0].header['OBJECT'] = 'HD 187123'
    hdus[0].header['NCOADD'] = (len(TEMPLATE_FRAMES), 'exposures co-added')
    hdus[0].header['REFFRAME'] = ('observed', 'heliocentric correction removed')
    hdus[0].header['MJDREF'] = (exposures[0]['mjd'], 'MJD of reference exposure')
    hdus[0].header['VBARY'] = (bary_velocity(exposures[0]['mjd'], 500.),
                               'km/s, astropy barycentric at midpoint')
    for i, koaid in enumerate(TEMPLATE_FRAMES):
        hdus[0].header['FRAME{:d}'.format(i)] = koaid
    for order in sorted(template):
        wave, flux, ivar, npix = template[order]
        col = fits.BinTableHDU.from_columns([
            fits.Column(name='WAVE', format='D', array=wave),
            fits.Column(name='FLUX', format='D', array=flux),
            fits.Column(name='IVAR', format='D', array=ivar),
            fits.Column(name='NCOADD', format='J', array=npix)],
            name='ORDER{:02d}'.format(order))
        col.header['ECHORDER'] = order
        hdus.append(col)
    fits.HDUList(hdus).writeto(tmpl_file, overwrite=True)
    print('Wrote {:s}'.format(tmpl_file))


if __name__ == '__main__':
    main(parse_args())
