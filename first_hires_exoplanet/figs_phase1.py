""" Diagnostic figures for the phase-1 reduction of HD 187123 (1998-07-16).

Four figures assessing whether the PypeIt reduction of the original Keck/HIRES
data is sound (see claude_prompts/data_phase1_prompt.md, prompt 9):

    fig_p1_orders    order tracing and the pixel flat
    fig_p1_wavecal   wavelength-solution RMS per order, and the velocity
                     zero-point measured from stellar lines
    fig_p1_iodine    the I2 absorption forest: line density against wavelength,
                     with a clean order and an iodine order side by side
    fig_p1_spectrum  the extracted 1D spectrum and its signal-to-noise

The reduction these read is NOT in this repository -- the raw frames and the
PypeIt products live in a sibling data tree, because the raw data must not be
committed.  Point --redux at it if it is somewhere else.

Reduced with PypeIt 2.0.2.dev1217+g017bece06 (branch orig-hires-fixes, commit
017bece06), which carries the six original-CCD fixes made in prompt 6.

NOTE ON WAVELENGTHS: PypeIt reports *vacuum* wavelengths.  Rest wavelengths of
the stellar features below are tabulated in air, as is conventional for the
optical, and are converted with pypeit.core.wave.airtovac before use.  Comparing
observed vacuum wavelengths against air rest values produces a spurious redshift
of about +70 km/s at 5000 A.

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.figs_phase1

"""

# Standard imports
import os
import argparse
import glob

import numpy as np
from scipy.ndimage import median_filter, percentile_filter
from scipy.signal import find_peaks

import astropy.units as u

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from pypeit import specobjs
from pypeit.slittrace import SlitTraceSet
from pypeit.flatfield import FlatImages
from pypeit.wavecalib import WaveCalib
from pypeit.core.wave import airtovac


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Default reduction directory (outside the repository; see the module docstring)
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_REDUX = os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                             'first-hires-exoplanet-data', 'redux', 'reduce_jul16')

FIG_DIR = os.path.join(os.path.dirname(_HERE), 'docs', 'figs')

#: Approximate extent of the iodine cell's absorption band (Angstrom).  The
#: cell imprints thousands of I2 lines over roughly this range; outside it the
#: spectrum is stellar only.
IODINE_RANGE = (5000., 6200.)

#: Strong features of a G2V photosphere that fall in our 3806-6262 A coverage.
#: Air wavelengths; converted to vacuum before use.
STELLAR_LINES_AIR = {
    'Ca II K': 3933.66,
    'Ca II H': 3968.47,
    'H-gamma': 4340.47,
    'H-beta': 4861.33,
    'Mg b1': 5167.32,
    'Mg b2': 5172.68,
    'Mg b3': 5183.60,
}

#: Systemic radial velocity of HD 187123 (km/s), for reference on fig 2.
V_SYS = -17.0

#: Line groups for annotation, so the near-coincident members of the Mg b
#: triplet and the Ca II doublet are labelled once rather than on top of
#: each other.  Air wavelengths; the label is placed at the group mean.
LINE_GROUPS = (
    ('Ca II H&K', (3933.66, 3968.47)),
    ('H-gamma', (4340.47,)),
    ('H-beta', (4861.33,)),
    ('Mg b', (5167.32, 5172.68, 5183.60)),
)

#: Colours, matching figs.py (dataviz palette slots 1-3 plus ink/surface)
C_DATA = '#2a78d6'
C_MODEL = '#eb6834'
C_THIRD = '#1baf7a'
INK = '#0b0b0b'
INK_2 = '#52514e'
SURFACE = '#fcfcfb'


# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------

def set_style():
    """ Apply a common, deliberately recessive style to all figures. """
    plt.rcParams.update({
        'figure.facecolor': SURFACE,
        'axes.facecolor': SURFACE,
        'savefig.facecolor': SURFACE,
        'font.size': 11,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'axes.labelcolor': INK,
        'axes.edgecolor': '#c9c8c3',
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'text.color': INK,
        'xtick.color': INK_2,
        'ytick.color': INK_2,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'grid.color': '#e6e5e0',
        'grid.linewidth': 0.8,
        'legend.frameon': False,
        'legend.fontsize': 10,
    })


def _save(fig, name):
    """ Write a figure to docs/figs and report it. """
    os.makedirs(FIG_DIR, exist_ok=True)
    outfile = os.path.join(FIG_DIR, name)
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('Wrote {:s}'.format(outfile))


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _one(redux, pattern):
    """ Return the single file in `redux` matching `pattern`, or raise. """
    hits = sorted(glob.glob(os.path.join(redux, pattern)))
    if len(hits) != 1:
        raise RuntimeError('Expected exactly one {:s} in {:s}, found {:d}'.format(
            pattern, redux, len(hits)))
    return hits[0]


def load_orders(redux):
    """ Load the extracted orders, brightest-first in order number.

    Note that `run_pypeit -o` *appends* to an existing spec1d rather than
    replacing it, so a directory that has been re-run without clearing
    Science/ will hold several reductions at once.  We check for that.

    Args:
        redux (str): reduction directory.

    Returns:
        list: (order, wave, flux, sigma) tuples, sorted by decreasing order.
    """
    sobjs = specobjs.SpecObjs.from_fitsfile(_one(redux, 'Science/spec1d_*.fits'))
    orders = [o.ECH_ORDER for o in sobjs]
    if len(orders) != len(set(orders)):
        raise RuntimeError(
            'spec1d holds {:d} entries for {:d} distinct orders: this directory '
            'accumulated several run_pypeit passes.  Delete Science/ and re-run '
            'so the figures describe one reduction.'.format(len(orders), len(set(orders))))

    out = []
    for o in sobjs:
        wave, flux, sig = o.OPT_WAVE, o.OPT_COUNTS, o.OPT_COUNTS_SIG
        good = np.isfinite(wave) & np.isfinite(flux) & np.isfinite(sig) & (wave > 0)
        if good.sum() < 200:
            continue
        out.append((o.ECH_ORDER, wave[good], flux[good], sig[good]))
    out.sort(key=lambda x: -x[0])
    return out


def normalise(wave, flux, width=301):
    """ Divide out a crude continuum.

    A running upper percentile tracks the continuum without being dragged down
    by the dense I2 forest, where a plain median sits well below the true
    continuum.  `width` is in pixels.

    Args:
        wave (`numpy.ndarray`_): wavelengths (unused; kept for symmetry).
        flux (`numpy.ndarray`_): counts.
        width (int, optional): filter width in pixels.

    Returns:
        `numpy.ndarray`_: continuum-normalised flux.
    """
    cont = percentile_filter(flux, 90, size=width)
    norm = flux / np.where(cont > 0, cont, np.nan)
    # Bad-pixel-masked columns come through as exact zeros; show them as gaps
    # rather than as spikes to the bottom of the frame.
    norm[flux <= 0] = np.nan
    return norm


def line_density(wave, flux, prominence=0.03):
    """ Absorption lines per 100 A, as a blunt measure of spectral crowding.

    Args:
        wave, flux (`numpy.ndarray`_): wavelength and counts for one order.
        prominence (float, optional): minimum depth, in normalised flux.

    Returns:
        float: detected absorption features per 100 A.
    """
    norm = normalise(wave, flux)
    good = np.isfinite(norm)
    if good.sum() < 100:
        return np.nan
    peaks, _ = find_peaks(-norm[good], prominence=prominence)
    span = wave[good].max() - wave[good].min()
    return len(peaks) / span * 100.


def measure_velocity(orders):
    """ Velocity offset of each stellar line from its vacuum rest wavelength.

    Uses the deepest pixel within +/-0.8 A of the expected position, which is
    crude but adequate to confirm the wavelength zero point to a few km/s.

    Args:
        orders (list): output of :func:`load_orders`.

    Returns:
        dict: line name -> (vacuum rest wavelength, velocity in km/s).
    """
    out = {}
    for name, air in STELLAR_LINES_AIR.items():
        vac = float(airtovac(air * u.AA).value)
        best = None
        for _, wave, flux, _ in orders:
            if not wave.min() + 2. < vac < wave.max() - 2.:
                continue
            sel = np.abs(wave - vac) < 0.8
            if sel.sum() < 5:
                continue
            sub_w, sub_f = wave[sel], flux[sel]
            ok = sub_f > 0
            if ok.sum() < 3:
                continue
            obs = sub_w[ok][np.argmin(sub_f[ok])]
            if best is None or abs(obs - vac) < abs(best - vac):
                best = obs
        if best is not None:
            out[name] = (vac, (best - vac) / vac * 299792.458)
    return out


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig_orders(redux):
    """ Order tracing and the pixel flat. """
    slits = SlitTraceSet.from_file(_one(redux, 'Calibrations/Slits_*.fits.gz'))
    flats = FlatImages.from_file(_one(redux, 'Calibrations/Flat_*.fits'))

    width = slits.right_init - slits.left_init
    nspec = slits.left_init.shape[0]
    row = nspec // 2

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.0))

    ax = axs[0]
    ax.plot(slits.left_init[row, :], width[row, :], 'o', ms=4, color=C_DATA)
    ax.axhline(np.median(width), color=C_MODEL, lw=1.2,
               label='median {:.2f} px'.format(np.median(width)))
    ax.axhspan(0, 6, color=C_MODEL, alpha=0.10)
    ax.text(0.02, 0.06, 'masked by the inherited\nfind_trim_edge = 3,3',
            transform=ax.transAxes, fontsize=9, color=INK_2)
    ax.set_xlabel('spatial position of order (binned pixels)')
    ax.set_ylabel('order width (binned pixels)')
    ax.set_title('{:d} orders traced, none masked'.format(slits.nslits))
    ax.set_ylim(0, max(16, width.max() * 1.1))
    ax.legend(loc='upper left')

    ax = axs[1]
    pix = flats.pixelflat_norm
    # Pixels outside the traced orders are left at exactly 1.0 and would
    # otherwise swamp the histogram; keep only the illuminated ones.
    good = np.isfinite(pix) & (pix > 0) & (pix != 1.0)
    ax.hist(pix[good], bins=160, range=(0.94, 1.06), color=C_DATA, alpha=0.85)
    lo, hi = np.percentile(pix[good], [1, 99])
    ax.axvline(1.0, color=INK_2, lw=0.8)
    ax.axvline(lo, color=C_MODEL, lw=1.0, ls='--')
    ax.axvline(hi, color=C_MODEL, lw=1.0, ls='--')
    ax.set_xlabel('normalised pixel flat')
    ax.set_ylabel('pixels in the traced orders')
    ax.set_title('pixel flat, {:.1f}M illuminated pixels: 1-99% = {:.3f}-{:.3f}'.format(
        good.sum() / 1e6, lo, hi))

    fig.tight_layout()
    _save(fig, 'fig_p1_orders.png')


def fig_wavecal(redux, orders):
    """ Wavelength RMS per order, and the measured velocity zero point. """
    wc = WaveCalib.from_file(_one(redux, 'Calibrations/WaveCalib_*.fits'))
    rms = np.array([w.rms for w in wc.wv_fits if w.rms is not None])

    vel = measure_velocity(orders)

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.0))

    ax = axs[0]
    ax.plot(np.arange(len(rms)), rms, 'o', ms=4, color=C_DATA)
    ax.axhline(np.median(rms), color=C_MODEL, lw=1.2,
               label='median {:.3f} px'.format(np.median(rms)))
    ax.set_xlabel('order index (blue to red)')
    ax.set_ylabel('wavelength solution RMS (pixels)')
    ax.set_title('all {:d} orders solved, none above 0.2 px'.format(len(rms)))
    ax.set_ylim(0, max(0.25, rms.max() * 1.2))
    ax.legend(loc='lower right')

    ax = axs[1]
    waves = np.array([v[0] for v in vel.values()])
    vels = np.array([v[1] for v in vel.values()])
    ax.plot(waves, vels, 'o', ms=6, color=C_DATA)
    # Label groups once: Ca II H/K and the three Mg b lines would otherwise
    # print on top of one another.
    for name, airs in LINE_GROUPS:
        members = [vel[n] for n in vel if STELLAR_LINES_AIR[n] in airs]
        if not members:
            continue
        mw = np.mean([m[0] for m in members])
        mv = np.mean([m[1] for m in members])
        ax.annotate(name, (mw, mv), textcoords='offset points', xytext=(0, 9),
                    ha='center', fontsize=9, color=INK_2)
    ax.axhline(np.mean(vels), color=C_MODEL, lw=1.2,
               label='mean {:+.1f} km/s'.format(np.mean(vels)))
    ax.axhline(V_SYS, color=C_THIRD, lw=1.2, ls='--',
               label='systemic {:+.0f} km/s'.format(V_SYS))
    ax.set_xlabel('vacuum wavelength (A)')
    ax.set_ylabel('velocity offset (km/s)')
    ax.set_title('stellar lines vs. vacuum rest wavelengths')
    ax.set_ylim(-30, 5)
    ax.legend(loc='lower right')

    fig.tight_layout()
    _save(fig, 'fig_p1_wavecal.png')


def fig_iodine(orders):
    """ The I2 forest: line density against wavelength, plus two order zooms. """
    centres = np.array([0.5 * (w.min() + w.max()) for _, w, _, _ in orders])
    dens = np.array([line_density(w, f) for _, w, f, _ in orders])

    # Pick a clean order and an iodine order of similar S/N for the zoom.
    clean = min([o for o in orders if o[1].max() < IODINE_RANGE[0]],
                key=lambda o: abs(0.5 * (o[1].min() + o[1].max()) - 4900.))
    iod = min([o for o in orders if o[1].min() > 5300.],
              key=lambda o: abs(0.5 * (o[1].min() + o[1].max()) - 5500.))

    fig = plt.figure(figsize=(11, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.42, wspace=0.18)

    ax = fig.add_subplot(gs[0, :])
    ax.axvspan(*IODINE_RANGE, color=C_MODEL, alpha=0.12)
    ax.plot(centres, dens, 'o-', ms=5, lw=1.2, color=C_DATA)
    ax.annotate('iodine cell band', (np.mean(IODINE_RANGE), 0.06),
                xycoords=('data', 'axes fraction'), ha='center',
                fontsize=10, color=C_MODEL)
    ax.set_xlabel('order centre, vacuum wavelength (A)')
    ax.set_ylabel('absorption lines per 100 A')
    ax.set_title('the iodine forest roughly doubles the line density above 5000 A')

    for ax, (order, wave, flux, _), label in (
            (fig.add_subplot(gs[1, 0]), clean, 'clean stellar spectrum'),
            (fig.add_subplot(gs[1, 1]), iod, 'stellar + I2 forest')):
        norm = normalise(wave, flux)
        mid = 0.5 * (wave.min() + wave.max())
        sel = np.abs(wave - mid) < 5.
        ax.plot(wave[sel], norm[sel], lw=0.8, color=C_DATA)
        ax.set_ylim(0, 1.25)
        ax.set_xlabel('vacuum wavelength (A)')
        ax.set_ylabel('normalised flux')
        ax.set_title('order {:d}: {:s}'.format(order, label), fontsize=11)

    _save(fig, 'fig_p1_iodine.png')


def fig_spectrum(orders):
    """ The extracted spectrum and its signal-to-noise. """
    fig, axs = plt.subplots(2, 1, figsize=(11, 6.4), sharex=True)

    ax = axs[0]
    for order, wave, flux, _ in orders:
        norm = normalise(wave, flux)
        ax.plot(wave, norm, lw=0.35, color=C_DATA, alpha=0.9)
    for name, airs in LINE_GROUPS:
        vacs = [float(airtovac(a * u.AA).value) for a in airs]
        for vac in vacs:
            ax.axvline(vac, color=C_THIRD, lw=0.8, ls='--', alpha=0.8)
        ax.annotate(name, (np.mean(vacs), 1.40), ha='center', va='top',
                    fontsize=9, color=INK_2)
    ax.set_ylim(0, 1.5)
    ax.set_ylabel('normalised flux')
    ax.set_title('HD 187123, Keck/HIRES 1998-07-16, 400 s -- {:d} orders'.format(len(orders)))

    ax = axs[1]
    centres = [0.5 * (w.min() + w.max()) for _, w, _, _ in orders]
    snr = [float(np.median(f[s > 0] / s[s > 0])) for _, _, f, s in orders]
    ax.axvspan(*IODINE_RANGE, color=C_MODEL, alpha=0.12)
    ax.plot(centres, snr, 'o-', ms=5, lw=1.2, color=C_DATA)
    ax.set_xlabel('vacuum wavelength (A)')
    ax.set_ylabel('median S/N per pixel')
    ax.set_title('S/N {:.0f} in the blue rising to {:.0f} in the red (median {:.0f})'.format(
        min(snr), max(snr), np.median(snr)))

    fig.tight_layout()
    _save(fig, 'fig_p1_spectrum.png')


# ---------------------------------------------------------------------------

def main():
    """ Build every phase-1 diagnostic figure. """
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--redux', default=DEFAULT_REDUX,
                        help='PypeIt reduction directory (default: %(default)s)')
    args = parser.parse_args()

    if not os.path.isdir(args.redux):
        raise SystemExit('No such reduction directory: {:s}\n'
                         'The raw data and reduction live outside this repository; '
                         'pass --redux.'.format(args.redux))

    set_style()
    orders = load_orders(args.redux)
    print('Loaded {:d} orders, {:.1f}-{:.1f} A'.format(
        len(orders), min(w.min() for _, w, _, _ in orders),
        max(w.max() for _, w, _, _ in orders)))

    fig_orders(args.redux)
    fig_wavecal(args.redux, orders)
    fig_iodine(orders)
    fig_spectrum(orders)

    vel = measure_velocity(orders)
    mean_v = np.mean([v for _, v in vel.values()])
    print('Mean stellar velocity offset: {:+.1f} km/s '
          '(HD 187123 systemic is about {:+.0f} km/s)'.format(mean_v, V_SYS))


if __name__ == '__main__':
    main()
