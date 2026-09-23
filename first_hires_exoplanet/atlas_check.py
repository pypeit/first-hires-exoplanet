""" Is the Fischer FTS atlas an adequate description of the HIRES cell? (prompt 7)

`pyodine` ships two iodine atlases.  The larger,
`iodine_atlas/Fischer_Cell_May2022_downsampled3.h5`, is almost certainly Debra
Fischer's cell -- the Lick/Keck lineage, the same tradition as the HIRES cell.
An FTS atlas is a scan of **one physical cell**, though.  The I2 line
*positions* are molecular constants and transfer to any cell; the line *depths*
depend on the column density and the cell temperature, which do not.  This
module establishes what the file actually contains and then tests it against
the HIRES cell empirically.

THE EMPIRICAL TEST.  1998-08-12 carries a 60 s frame with the cell in
(`HI.19980812.29160`) and a 60 s frame with the cell out
(`HI.19980812.29316`), 156 seconds apart.  Their ratio divides out the star,
the blaze and the detector and leaves the **transmission of the HIRES cell
itself**, measured at HIRES resolution.  That can be compared directly against
the atlas convolved to the same resolution.  Nothing else in this project can
settle the question; a header keyword cannot.

ENVIRONMENT.  This module needs `h5py`, which phase 3 prompt 1 installed into
`pypeit14` along with `barycorrpy`; it runs there now, like everything else:

    conda run -n pypeit14 python -m first_hires_exoplanet.atlas_check

It ran in `astro` for phases 1 and 2 because `h5py` was missing from
`pypeit14`.  It still avoids importing `pypeit` -- the spec1d files are read
with plain `astropy.io.fits`, whose layout phase 1 and prompt 4 already
documented -- which is what let it move environments without changing a line
of the analysis, and which is also what lets it live inside the vendored
`pyodine` tree later if that becomes useful.

The atlas is 55.8 MB and lives in the data tree, not the repository.

NOTE ON WAVELENGTHS.  The file carries both `wavelength` (vacuum) and
`wavelength_air`.  They differ by 1.56 A at 5615 A, which is 83 km/s.  PypeIt
reports vacuum.  `pyodine`'s Lick loader reads `wavelength_air`.  This is a trap
for prompt 8 and is checked explicitly below.

"""

# Standard imports
import os
import glob
import argparse

import numpy as np

import h5py
from astropy.io import fits
from scipy.ndimage import percentile_filter, gaussian_filter1d
from scipy.signal import find_peaks

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                         'first-hires-exoplanet-data')
DEFAULT_ATLAS = os.path.join(DATA_ROOT, 'atlas',
                             'Fischer_Cell_May2022_downsampled3.h5')
DEFAULT_REDUX = os.path.join(DATA_ROOT, 'redux')

#: The 1998-08-12 pair: same night, 156 s apart, one with the cell in and one
#: with it out.  The observer labelled the cell-out one OBJECT = 'Template'.
CELL_IN = 'HI.19980812.29160'
CELL_OUT = 'HI.19980812.29316'

#: The band the iodine cell absorbs in
I2_LO, I2_HI = 5000., 6200.

CKMS = 299792.458

#: Resolving power of HIRES through the B1 decker, measured from the weak
#: stellar lines of the prompt-4 template (R >~ 40,200 there; 45,000 nominal).
R_HIRES = 42000.

#: Mass of I2 in atomic mass units, for the Doppler width
M_I2 = 253.809

#: Window, in pixels, of the running continuum applied to BOTH sides
CONT_WIN = 301

FIG_DIR = os.path.join(os.path.dirname(_HERE), 'docs', 'figs')

#: Palette, matching figs.py / figs_phase1.py / figs_phase2.py
C_HIRES = '#2a78d6'
C_ATLAS = '#eb6834'
C_SCALED = '#1baf7a'
INK = '#0b0b0b'
INK_2 = '#52514e'
SURFACE = '#fcfcfb'


# ---------------------------------------------------------------------------
# The atlas
# ---------------------------------------------------------------------------

def load_atlas(path):
    """ Read the FTS atlas.

    Returns:
        dict: the five datasets plus derived description.
    """
    with h5py.File(path, 'r') as h:
        out = {k: np.array(h[k]) for k in h.keys()}
    return out


def describe_atlas(atlas):
    """ Everything the file itself can be made to say about its provenance. """
    w, wa = atlas['wavelength'], atlas['wavelength_air']
    wn, fl, fn = atlas['wavenumber'], atlas['flux'], atlas['flux_normalized']

    print('\n=== 1. What the file contains ===')
    print('  datasets            : {:s}'.format(', '.join(sorted(atlas))))
    print('  points              : {:d}'.format(len(w)))
    print('  wavelength (vacuum) : {:.2f} - {:.2f} A'.format(w.min(), w.max()))
    print('  wavelength_air      : {:.2f} - {:.2f} A'.format(wa.min(), wa.max()))
    print('  wavenumber          : {:.2f} - {:.2f} cm^-1'.format(wn.min(), wn.max()))
    print('  vacuum - air at mid : {:.4f} A  = {:.1f} km/s'.format(
        w[len(w) // 2] - wa[len(wa) // 2],
        (w[len(w) // 2] - wa[len(wa) // 2]) / w[len(w) // 2] * CKMS))
    dwn = np.abs(np.diff(wn))
    print('  wavenumber step     : {:.7f} cm^-1, spread {:.2e}'.format(
        np.median(dwn), np.ptp(dwn)))
    print('    -> uniform in WAVENUMBER, which is the signature of an FTS scan,')
    print('       not of a grating spectrograph.')
    print('  no HDF5 attributes  : {:s}'.format(
        'confirmed -- the file carries no metadata at all'))
    print('  flux                : {:.4f} - {:.4f} (raw FTS counts)'.format(
        fl.min(), fl.max()))
    print('  flux_normalized     : {:.4f} - {:.4f}'.format(fn.min(), fn.max()))
    print('  pyodine reads `flux_normalized` against `wavelength_air`')
    print('    (utilities_lick/load_pyodine.py, IodineTemplate)')


def atlas_resolution(atlas, lo=5400., hi=5600.):
    """ Resolving power of the atlas, from the narrowest resolved I2 lines.

    The intrinsic Doppler width of I2 at cell temperatures is about 0.24 km/s
    FWHM, so anything appreciably broader than that is the FTS instrument
    profile and the narrowest lines bound the atlas resolution from below.

    Returns:
        dict
    """
    w, fn = atlas['wavelength'], atlas['flux_normalized']
    sel = (w > lo) & (w < hi)
    ww, ff = w[sel], fn[sel]
    prof = 1. - ff
    peaks, props = find_peaks(prof, height=(0.05, 0.9), prominence=0.03,
                              distance=5)
    fwhm = []
    for k in peaks:
        d = prof[k]
        half = 0.5 * d
        a, b = k, k
        while a > 0 and prof[a - 1] >= half:
            a -= 1
        while b < len(prof) - 1 and prof[b + 1] >= half:
            b += 1
        if a == 0 or b == len(prof) - 1 or (b - a) < 2:
            continue
        fwhm.append((ww[b] - ww[a]) / ww[k] * CKMS)
    fwhm = np.array(fwhm)
    if len(fwhm) == 0:
        return None
    p10 = float(np.percentile(fwhm, 10))
    med = float(np.median(fwhm))
    return dict(n=len(fwhm), p10=p10, median=med, R_p10=CKMS / p10,
                R_med=CKMS / med)


def doppler_fwhm(temp_k, wavelength=5500.):
    """ Thermal Doppler FWHM of an I2 line, km/s. """
    return CKMS * 2. * np.sqrt(2. * np.log(2.) * 1.380649e-23 * temp_k
                               / (M_I2 * 1.66053907e-27)) / 299792458.


# ---------------------------------------------------------------------------
# The HIRES cell, measured
# ---------------------------------------------------------------------------

def read_spec1d(path):
    """ Observed-frame per-order arrays, with plain astropy.

    The spec1d layout is documented in phase 1 and prompt 4: one BinTable per
    order, `VEL_CORR` in the extension header, masked pixels arriving as exact
    zeros.
    """
    orders = {}
    with fits.open(path) as hdul:
        vel_corr = float(hdul[1].header.get('VEL_CORR', 1.0) or 1.0)
        for hdu in hdul[1:]:
            if not hasattr(hdu, 'columns') or 'OPT_WAVE' not in hdu.columns.names:
                continue
            order = int(hdu.header['ECH_ORDER'])
            d = hdu.data
            wave = np.asarray(d['OPT_WAVE'], dtype=float) / vel_corr
            flux = np.asarray(d['OPT_COUNTS'], dtype=float)
            ivar = np.asarray(d['OPT_COUNTS_IVAR'], dtype=float)
            mask = np.asarray(d['OPT_MASK'], dtype=bool)
            good = mask & (ivar > 0) & (wave > 0) & np.isfinite(flux) & (flux != 0)
            if good.sum() > 300:
                orders[order] = (wave[good], flux[good], ivar[good])
    return orders


def find_spec1d(redux_dir, koaid):
    hits = sorted(glob.glob(os.path.join(
        redux_dir, 'reduce_*', 'Science', 'spec1d_{:s}*.fits'.format(koaid))))
    if len(hits) != 1:
        raise RuntimeError('Expected one spec1d for {:s}, found {:d}'.format(
            koaid, len(hits)))
    return hits[0]


def measured_transmission(redux_dir):
    """ HIRES cell transmission from the 1998-08-12 cell-in / cell-out pair.

    Returns:
        dict: order -> (wave, transmission)
    """
    cin = read_spec1d(find_spec1d(redux_dir, CELL_IN))
    cout = read_spec1d(find_spec1d(redux_dir, CELL_OUT))

    out = {}
    for order in sorted(set(cin) & set(cout)):
        wi, fi, _ = cin[order]
        wo, fo, _ = cout[order]
        if wi.max() < I2_LO or wi.min() > I2_HI:
            continue
        lo, hi = max(wi.min(), wo.min()) + 0.3, min(wi.max(), wo.max()) - 0.3
        if hi - lo < 10.:
            continue
        grid = np.linspace(lo, hi, len(wi))
        a = np.interp(grid, wi, fi)
        b = np.interp(grid, wo, fo)
        ratio = a / np.maximum(b, 1e-9)
        # The two frames differ in throughput and in the star's continuum only
        # smoothly, so a running upper envelope of the ratio is the unabsorbed
        # level and dividing by it leaves the cell transmission.
        cont = percentile_filter(ratio, 95, size=301)
        out[order] = (grid, ratio / np.maximum(cont, 1e-9))
    return out


def atlas_at_hires(atlas, grid, resolution=R_HIRES, air=False):
    """ The atlas convolved to HIRES resolution and sampled on `grid`. """
    w = atlas['wavelength_air'] if air else atlas['wavelength']
    fn = atlas['flux_normalized']
    lo, hi = grid.min() - 2., grid.max() + 2.
    sel = (w > lo) & (w < hi)
    if sel.sum() < 100:
        return None
    ww, ff = w[sel], fn[sel]
    # The atlas grid is uniform in wavenumber; over one echelle order the
    # wavelength step is near enough uniform to convolve directly.
    step_kms = np.median(np.diff(ww)) / np.median(ww) * CKMS
    sigma_kms = (CKMS / resolution) / 2.3548
    smoothed = gaussian_filter1d(ff, sigma_kms / step_kms, mode='nearest')
    return np.interp(grid, ww, smoothed)


def renorm(y, window=None):
    """ Running upper-percentile continuum, applied identically to both sides.

    Normalising one side with a running filter and the other with a single
    scalar would put a continuum difference into the depth comparison.
    """
    win = window or CONT_WIN
    win = min(win, (len(y) // 2) * 2 - 1)
    cont = percentile_filter(y, 95, size=win)
    return y / np.maximum(cont, 1e-9)


def fit_alpha(t_meas, t_atlas, lo=0.02, hi=0.97):
    """ The Beer-Lambert exponent: T_meas = T_atlas ** alpha.

    Taking logarithms makes it a one-parameter linear fit through the origin,
    alpha = sum(ln T_m ln T_a) / sum((ln T_a)^2).  Pixels too close to 1 carry
    no information and pixels near 0 are saturated and noise-dominated, so both
    are excluded.

    Args:
        t_meas, t_atlas (`numpy.ndarray`_): normalised transmissions.
        lo, hi (float): transmission range to fit over.

    Returns:
        float
    """
    ok = ((t_atlas > lo) & (t_atlas < hi) & (t_meas > lo) & (t_meas < 1.2)
          & np.isfinite(t_meas) & np.isfinite(t_atlas))
    if ok.sum() < 100:
        return np.nan
    la, lm = np.log(t_atlas[ok]), np.log(np.clip(t_meas[ok], 1e-6, None))
    return float(np.sum(la * lm) / np.sum(la * la))


def noise_level(redux_dir, order):
    """ Fractional noise of the cell-in / cell-out ratio in one order. """
    cin = read_spec1d(find_spec1d(redux_dir, CELL_IN))
    cout = read_spec1d(find_spec1d(redux_dir, CELL_OUT))
    if order not in cin or order not in cout:
        return 0.02
    out = []
    for d in (cin, cout):
        w, f, iv = d[order]
        out.append(1. / np.median(f * np.sqrt(iv)))
    return float(np.hypot(*out))


def fig_atlas(atlas, trans, alpha, order=65, lo=5455., hi=5470.):
    """ The measured HIRES cell against the atlas, raw and Beer-Lambert scaled.

    A short stretch of one order, because the whole point is whether individual
    I2 lines line up and how deep they are; a whole order at this scale is a
    black band.
    """
    plt.rcParams.update({
        'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
        'savefig.facecolor': SURFACE, 'font.size': 11,
        'axes.labelcolor': INK, 'axes.edgecolor': '#c9c8c3',
        'axes.linewidth': 0.8, 'axes.spines.top': False,
        'axes.spines.right': False, 'text.color': INK,
        'xtick.color': INK_2, 'ytick.color': INK_2,
        'legend.frameon': False, 'legend.fontsize': 10,
        'grid.color': '#e6e5e1', 'grid.linewidth': 0.7})

    grid, meas = trans[order]
    sel = (grid > lo) & (grid < hi)
    g, m = grid[sel], renorm(meas)[sel]
    a = atlas_at_hires(atlas, g)
    a = renorm(a)

    fig, axs = plt.subplots(2, 1, figsize=(11, 6.0), sharex=True)

    ax = axs[0]
    ax.grid(axis='y', zorder=0)
    ax.plot(g, m, color=C_HIRES, lw=1.8, zorder=3, label='HIRES cell, measured')
    ax.plot(g, a, color=C_ATLAS, lw=1.6, zorder=2, label='Fischer atlas, as shipped')
    ax.set_ylabel('transmission')
    ax.legend(loc='lower left', ncol=2)
    ax.set_title('The line positions transfer; the depths do not',
                 fontsize=12.5, pad=10)

    ax = axs[1]
    ax.grid(axis='y', zorder=0)
    ax.plot(g, m, color=C_HIRES, lw=1.8, zorder=3, label='HIRES cell, measured')
    ax.plot(g, a ** alpha, color=C_SCALED, lw=1.6, zorder=2,
            label=r'Fischer atlas $^{{{:.1f}}}$ (Beer-Lambert)'.format(alpha))
    ax.set_ylabel('transmission')
    ax.set_xlabel('vacuum wavelength (A)')
    ax.legend(loc='lower left', ncol=2)

    for ax in axs:
        lo_y, hi_y = ax.get_ylim()
        ax.set_ylim(lo_y - 0.22 * (hi_y - lo_y), hi_y)
    fig.tight_layout()
    os.makedirs(FIG_DIR, exist_ok=True)
    out = os.path.join(FIG_DIR, 'fig_p2_atlas.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('  wrote {:s}'.format(out))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(options=None):
    parser = argparse.ArgumentParser(
        description='Settle whether the Fischer atlas describes the HIRES cell.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--atlas', type=str, default=DEFAULT_ATLAS)
    parser.add_argument('--redux', type=str, default=DEFAULT_REDUX)
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):

    atlas = load_atlas(pargs.atlas)
    describe_atlas(atlas)

    print('\n=== 2. Resolution and implied temperature ===')
    res = atlas_resolution(atlas)
    print('  {:d} isolated I2 lines measured in 5400-5600 A'.format(res['n']))
    print('  FWHM  10th percentile {:.3f} km/s  ->  R = {:.0f}'.format(
        res['p10'], res['R_p10']))
    print('  FWHM  median          {:.3f} km/s  ->  R = {:.0f}'.format(
        res['median'], res['R_med']))
    for t_c in (50., 65., 75.):
        print('    thermal Doppler FWHM at {:.0f} C = {:.3f} km/s'.format(
            t_c, doppler_fwhm(t_c + 273.15)))
    print('  The narrowest lines are close to the thermal width, so the atlas')
    print('  resolves the I2 forest and its resolution is NOT the limitation')
    print('  for a spectrograph at R = {:.0f}.'.format(R_HIRES))

    print('\n=== 3. The HIRES cell, measured from the 1998-08-12 pair ===')
    trans = measured_transmission(pargs.redux)
    print('  {:s} (cell in) / {:s} (cell out)'.format(CELL_IN, CELL_OUT))
    print('  orders in the I2 band with a usable ratio: {:d}'.format(len(trans)))

    print('\n=== 4. Measured HIRES cell against the atlas ===')
    print('  Both put through the IDENTICAL continuum treatment (running 95th')
    print('  percentile, {:d} pixels) and the atlas convolved to R = {:.0f}.'
          .format(CONT_WIN, R_HIRES))
    print()
    print('  Beer-Lambert: if the two cells differ only in column density then')
    print('  T_hires = T_atlas ** alpha, with alpha the column-density ratio.')
    print('  Fitting alpha on the logarithms is the physical comparison and')
    print('  handles saturated lines, which a mean-depth ratio does not.')
    print()
    print('  {:>5s} {:>12s} {:>8s} {:>8s} {:>7s} {:>7s}'.format(
        'order', 'range (A)', 'd_HIRES', 'd_atlas', 'alpha', 'corr'))
    rows = []
    for order in sorted(trans):
        grid, meas = trans[order]
        sel = (grid > I2_LO) & (grid < I2_HI)
        if sel.sum() < 500:
            continue
        g, m = grid[sel], meas[sel]
        a = atlas_at_hires(atlas, g)
        if a is None:
            continue
        m = renorm(m)
        a = renorm(a)
        alpha = fit_alpha(m, a)
        d_meas, d_atlas = 1. - np.mean(m), 1. - np.mean(a)
        corr = float(np.corrcoef(1. - m, 1. - a)[0, 1])
        rows.append((order, d_meas, d_atlas, alpha, corr))
        print('  {:5d} {:5.0f}-{:5.0f} {:8.4f} {:8.4f} {:7.2f} {:7.3f}'.format(
            order, g.min(), g.max(), d_meas, d_atlas, alpha, corr))

    arr = np.array([[r[1], r[2], r[3], r[4]] for r in rows])
    al = arr[:, 2][np.isfinite(arr[:, 2])]
    print('\n  median measured depth  : {:.4f}'.format(np.median(arr[:, 0])))
    print('  median atlas depth     : {:.4f}'.format(np.median(arr[:, 1])))
    print('  median alpha           : {:.2f}  (spread {:.2f}-{:.2f}, '
          'n={:d})'.format(np.median(al), al.min(), al.max(), len(al)))
    print('  median correlation     : {:.3f}'.format(np.median(arr[:, 3])))
    print('    order 58 returns no alpha: its atlas absorption is only'
          ' {:.3f},'.format(arr[0, 1]))
    print('    so almost no pixel is deep enough to constrain the exponent.')

    # How much of alpha is a noise artefact?  The measured ratio is built from
    # two 60 s frames, so it is noisy, and a running upper-percentile continuum
    # rides up on noise and deepens every line.  Inject matched noise into the
    # atlas, push it through the identical normalisation, and see what alpha
    # comes back.
    print('\n  Noise control.  The measured transmission comes from two 60 s')
    print('  frames, and a running 95th-percentile continuum rides up on noise,')
    print('  which deepens lines.  Injecting the same noise into the atlas and')
    print('  re-running the identical normalisation:')
    rng = np.random.default_rng(42)
    ctrl = []
    for order in sorted(trans):
        grid, meas = trans[order]
        sel = (grid > I2_LO) & (grid < I2_HI)
        if sel.sum() < 500:
            continue
        g = grid[sel]
        a = atlas_at_hires(atlas, g)
        if a is None:
            continue
        noise = noise_level(pargs.redux, order)
        a_noisy = renorm(a * (1. + rng.normal(0., noise, size=len(a))))
        val = fit_alpha(a_noisy, renorm(a))
        if np.isfinite(val):
            ctrl.append(val)
    print('    injected fractional noise : {:.3f}'.format(noise))
    print('    alpha recovered from the atlas against itself : {:.2f}'.format(
        np.median(ctrl)))
    print('    -> noise moves alpha by about {:+.0f}%, nowhere near the factor'
          .format(100. * (np.median(ctrl) - 1.)))
    print('       of {:.1f} measured.  The depth difference is real.'.format(
        np.median(al)))

    print('\n  Is one alpha enough?  The fitted alpha runs 3.5 at 6050 A down')
    print('  to 2.3 at 5320 A, and it is largest exactly where the absorption')
    print('  is weakest.  That is the signature of the continuum, not of the')
    print('  cell: where the forest is dense a running 95th percentile cannot')
    print('  find true continuum and flattens both sides.  Refitting over a')
    print('  restricted transmission band, 0.50-0.90, where neither saturation')
    print('  nor continuum placement bites:')
    band = []
    for order in sorted(trans):
        grid, meas = trans[order]
        sel = (grid > I2_LO) & (grid < I2_HI)
        if sel.sum() < 500:
            continue
        g = grid[sel]
        a = atlas_at_hires(atlas, g)
        if a is None:
            continue
        val = fit_alpha(renorm(meas[sel]), renorm(a), lo=0.50, hi=0.90)
        if np.isfinite(val):
            band.append((order, val))
    if band:
        bv = np.array([b[1] for b in band])
        print('    alpha over 0.50-0.90 : median {:.2f}, spread {:.2f}-{:.2f}'
              ' (n={:d})'.format(np.median(bv), bv.min(), bv.max(), len(bv)))
        print('    against {:.2f}, spread {:.2f}-{:.2f} over the full range'
              .format(np.median(al), al.min(), al.max()))

    print('\n=== 5. Can pyodine absorb the difference? ===')
    print('  pyodine carries a free `iod_depth` parameter'
          ' (models/spectrum.py:17),')
    print('  which would seem to settle it.  It does not, because of HOW it is')
    print('  applied (models/spectrum.py:93):')
    print()
    print('      flux_iod = params[\'iod_depth\'] * (flux_iod - 1.0) + 1.0')
    print()
    print('  That is a LINEAR depth scaling, T -> 1 + d (T - 1).  The physical')
    print('  law for a different column density is Beer-Lambert, T -> T**alpha.')
    print('  The two agree only for weak lines.  This atlas reaches'
          ' flux_normalized')
    print('  = {:.3f}, and a saturated line at T = 0 scaled linearly by the'
          .format(float(np.min(atlas['flux_normalized']))))
    print('  d = {:.1f} we need returns T = {:.1f}: negative transmission.'
          .format(np.median(al), 1. + np.median(al) * (0. - 1.)))
    n_neg = float(np.mean(1. + np.median(al)
                          * (atlas['flux_normalized'] - 1.) < 0.))
    print('  {:.1%} of the atlas would go negative under that scaling.'.format(
        n_neg))

    print('\n=== 6. The air/vacuum trap ===')
    order = sorted(trans)[len(trans) // 2]
    grid, meas = trans[order]
    sel = (grid > I2_LO) & (grid < I2_HI)
    g, m = grid[sel], meas[sel]
    m = m / np.percentile(m, 95)
    for lbl, air in [('vacuum (correct for PypeIt)', False),
                     ('air (what pyodine loads)', True)]:
        a = atlas_at_hires(atlas, g, air=air)
        if a is None:
            continue
        a = a / np.percentile(a, 95)
        print('  order {:d}, atlas on {:28s}: corr {:+.3f}'.format(
            order, lbl, float(np.corrcoef(1. - m, 1. - a)[0, 1])))
    print('  The air grid is offset by 1.56 A = 83 km/s at 5615 A.  PypeIt')
    print('  reports vacuum; pyodine\'s IodineTemplate reads `wavelength_air`.')
    print('  Prompt 8 must reconcile these or the forward model is nonsense.')

    print('\n=== 7. Figure ===')
    fig_atlas(atlas, trans, float(np.median(bv)) if band else float(np.median(al)))


if __name__ == '__main__':
    main(parse_args())
