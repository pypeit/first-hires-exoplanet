""" Assessment and figures for the phase-2 cross-correlation velocities.

Prompt 6 of claude_prompts/data_phase2_prompt.md: assess what prompt 5
produced, compare it against the published 3.097-day Keplerian and against the
modern catalogue, and say plainly what the comparison does and does not
demonstrate.

    fig_p2_timeseries  our velocities and the modern catalogue's over the same
                       nine months, on their own scales
    fig_p2_phasefold   both folded on the 3.0966-day period: the catalogue
                       traces the Keplerian, ours does not
    fig_p2_compare     ours against the catalogue epoch by epoch, with the 1:1
                       line the comparison would follow if we could see it
    fig_p2_precision   what each stage of the measurement achieves, against
                       what the planet requires

The modern velocities are Teklu et al. 2025, A&A, 702, A68, VizieR
J/A+A/702/A68 table A.1, stored at data/hd187123_hires_rv.tsv.

NOTE ON UNITS: VizieR declares the RV columns as km/s.  That is wrong -- they
are relative velocities in m/s, as figs.py already records.

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.figs_phase2

"""

# Standard imports
import os
import argparse

import numpy as np

from astropy.table import Table

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_HERE, 'data')
FIG_DIR = os.path.join(os.path.dirname(_HERE), 'docs', 'figs')

VEL_FILE = os.path.join(DATA_DIR, 'xcorr_velocities.csv')
ORD_FILE = os.path.join(DATA_DIR, 'xcorr_velocities_per_order.csv')
CAT_FILE = os.path.join(DATA_DIR, 'hd187123_hires_rv.tsv')

#: Orbital period of HD 187123 b (d), the refined modern value
PERIOD = 3.0965828

#: Published semiamplitude (m/s), Butler et al. (1998)
K_PUBLISHED = 72.

#: Largest separation, in days, at which one of our epochs is taken to be the
#: same observation as a catalogue row.  Our times are JD(UTC) at mid-exposure
#: and the catalogue's are BJD(TDB); the two differ by the light-travel term,
#: at most about 0.006 d for this target.
MATCH_TOL = 0.02

#: Colours: slots 1-3 of the dataviz skill's validated categorical palette, in
#: fixed order, plus its light-mode ink and surface tokens.  Same assignment as
#: figs.py and figs_phase1.py, so a reader moving between them keeps the
#: mapping: blue is always this project's own measurement, orange the published
#: or modern comparison.
C_OURS = '#2a78d6'
C_CAT = '#eb6834'
C_THIRD = '#1baf7a'
INK = '#0b0b0b'
INK_2 = '#52514e'
SURFACE = '#fcfcfb'
GRID = '#e6e5e1'


def set_style():
    """ The recessive house style shared with figs.py and figs_phase1.py. """
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
        'legend.frameon': False,
        'legend.fontsize': 10,
        'grid.color': GRID,
        'grid.linewidth': 0.7,
    })


def save(fig, name):
    """ Write a figure into docs/figs. """
    os.makedirs(FIG_DIR, exist_ok=True)
    out = os.path.join(FIG_DIR, name)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('  wrote {:s}'.format(out))
    return out


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load_catalogue(path=CAT_FILE):
    """ The modern HIRES velocities, discovery era only.

    Returns:
        `astropy.table.Table`_: columns bjd, rv, e_rv (m/s).
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith('HD187123'):
                continue
            parts = [p.strip() for p in line.split('\t')]
            rows.append(dict(bjd=float(parts[1]), rv=float(parts[2]),
                             e_rv=float(parts[3])))
    return Table(rows)


def load_ours(path=VEL_FILE):
    """ Our cross-correlation velocities, with a JD(UTC) mid-exposure time. """
    tbl = Table.read(path)
    tbl['jd_mid'] = (np.asarray(tbl['mjd'], dtype=float)
                     + 0.5 * np.asarray(tbl['exptime'], dtype=float) / 86400.
                     + 2400000.5)
    return tbl


def match(ours, cat, tol=MATCH_TOL):
    """ Pair our epochs with catalogue rows by time.

    Returns:
        `astropy.table.Table`_: one row per matched pair.
    """
    rows = []
    for r in ours:
        d = np.abs(np.asarray(cat['bjd'], dtype=float) - r['jd_mid'])
        i = int(np.argmin(d))
        if d[i] > tol:
            continue
        rows.append(dict(koaid=r['koaid'], jd=r['jd_mid'], dt_days=float(d[i]),
                         ours=float(r['v_rel']),
                         sigma=float(r['sigma_empirical']),
                         cat=float(cat['rv'][i]), e_cat=float(cat['e_rv'][i])))
    return Table(rows)


def fit_circular(t, v, period=PERIOD, e=None):
    """ Least-squares circular Keplerian at a fixed period.

    The model is v = gamma + A cos(2 pi t / P) + B sin(2 pi t / P), which is
    linear in its three parameters, so there is no minimiser to get stuck.

    Args:
        t (`numpy.ndarray`_): times, days.
        v (`numpy.ndarray`_): velocities, m/s.
        period (float): period to force, days.
        e (`numpy.ndarray`_, optional): uncertainties for weighting.

    Returns:
        dict: 'K', 'sigma_K', 'phase', 'gamma', 'rms', 'model' (callable).
    """
    ph = 2. * np.pi * t / period
    M = np.column_stack([np.ones_like(t), np.cos(ph), np.sin(ph)])
    w = np.ones_like(t) if e is None else 1. / np.asarray(e, dtype=float) ** 2
    W = np.diag(w)
    cov = np.linalg.inv(M.T @ W @ M)
    p = cov @ (M.T @ W @ v)
    gamma, A, B = p
    K = float(np.hypot(A, B))
    # Propagate to K = sqrt(A^2+B^2)
    dK = np.array([0., A / K, B / K]) if K > 0 else np.zeros(3)
    sigma_K = float(np.sqrt(dK @ cov @ dK))
    resid = v - M @ p
    # Scale the formal error by the fit quality, so a poor fit does not report
    # a confident K
    dof = max(len(t) - 3, 1)
    chi2_red = float(np.sum(w * resid ** 2) / dof)
    if e is None:
        sigma_K *= np.sqrt(chi2_red)
    return dict(K=K, sigma_K=sigma_K, gamma=float(gamma),
                phase=float(np.arctan2(B, A)),
                rms=float(np.std(resid, ddof=1)), chi2_red=chi2_red,
                model=lambda tt: (gamma + A * np.cos(2 * np.pi * tt / period)
                                  + B * np.sin(2 * np.pi * tt / period)))


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig_timeseries(pairs, cat, ours):
    """ The two velocity sets over the same nine months, on their own scales. """
    fig, axs = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True)
    t0 = 2450800.

    ax = axs[0]
    ax.grid(axis='y', zorder=0)
    ax.errorbar(np.asarray(ours['jd_mid']) - t0, np.asarray(ours['v_rel']),
                yerr=np.asarray(ours['sigma_empirical']), fmt='o', ms=6,
                color=C_OURS, ecolor=C_OURS, elinewidth=1.4, capsize=0,
                mfc=SURFACE, mew=1.8, zorder=3)
    ax.axhline(0., color=INK_2, lw=0.8, ls=':', zorder=1)
    ax.set_ylabel('relative velocity (m/s)')
    # One series per panel, so it is named in place rather than in a legend box
    # that would have to sit on top of the data.
    ax.text(0.012, 0.95, 'this reduction, cross-correlation',
            transform=ax.transAxes, va='top', color=C_OURS, fontsize=11)
    ax.text(0.012, 0.845, 'epoch-to-epoch rms {:.0f} m/s'.format(
        np.std(np.asarray(ours['v_rel'], dtype=float), ddof=1)),
        transform=ax.transAxes, va='top', color=INK_2, fontsize=10)

    ax = axs[1]
    ax.grid(axis='y', zorder=0)
    sel = np.asarray(cat['bjd']) < 2451110.
    ax.errorbar(np.asarray(cat['bjd'])[sel] - t0, np.asarray(cat['rv'])[sel],
                yerr=np.asarray(cat['e_rv'])[sel], fmt='o', ms=6,
                color=C_CAT, ecolor=C_CAT, elinewidth=1.4, capsize=0,
                mfc=SURFACE, mew=1.8, zorder=3)
    ax.axhline(0., color=INK_2, lw=0.8, ls=':', zorder=1)
    ax.set_ylabel('relative velocity (m/s)')
    ax.set_xlabel('BJD - 2450800 (days)')
    ax.text(0.012, 0.95, 'Teklu et al. 2025, the same frames',
            transform=ax.transAxes, va='top', color=C_CAT, fontsize=11)
    ax.text(0.012, 0.845, 'median error bar {:.1f} m/s'.format(
        np.median(np.asarray(cat['e_rv'])[sel])),
        transform=ax.transAxes, va='top', color=INK_2, fontsize=10)

    for ax in axs:
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.30 * (hi - lo))
    fig.suptitle('The same photons, measured two ways', y=0.985, fontsize=12.5)
    fig.tight_layout()
    return save(fig, 'fig_p2_timeseries.png')


def fig_phasefold(pairs, cat):
    """ Both sets folded on the 3.0966-day period. """
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
    sel = np.asarray(cat['bjd']) < 2451110.
    tc, vc, ec = (np.asarray(cat['bjd'])[sel], np.asarray(cat['rv'])[sel],
                  np.asarray(cat['e_rv'])[sel])
    fit = fit_circular(tc, vc, e=ec)

    grid = np.linspace(0., 1., 400)
    curve = fit['model'](grid * PERIOD)

    for ax, (t, v, e, colour, label) in zip(axs, [
            (tc, vc, ec, C_CAT, 'Teklu et al. 2025'),
            (np.asarray(pairs['jd']), np.asarray(pairs['ours']),
             np.asarray(pairs['sigma']), C_OURS, 'this reduction')]):
        ax.grid(axis='y', zorder=0)
        ax.plot(grid, curve - fit['gamma'], color=INK_2, lw=1.6, zorder=2,
                label='K = {:.0f} m/s Keplerian'.format(fit['K']))
        ax.errorbar((t / PERIOD) % 1., v - np.mean(v), yerr=e, fmt='o', ms=6,
                    color=colour, ecolor=colour, elinewidth=1.4, capsize=0,
                    mfc=SURFACE, mew=1.8, zorder=3, label=label)
        ax.set_xlabel('orbital phase (P = {:.4f} d)'.format(PERIOD))
        ax.legend(loc='upper right')
    axs[0].set_ylabel('velocity (m/s)')
    fig.suptitle('Folded on the published period', y=0.99, fontsize=12.5)
    fig.tight_layout()
    return save(fig, 'fig_p2_phasefold.png')


def fig_compare(pairs):
    """ Ours against the catalogue, epoch by epoch.

    Equal aspect is essential and is the whole message: on a square where one
    axis is the other's units, our points collapse into a vertical stripe
    because our scatter is fifteen times the catalogue's.
    """
    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    x = np.asarray(pairs['cat'], dtype=float)
    y = np.asarray(pairs['ours'], dtype=float)
    x = x - np.mean(x)
    y = y - np.mean(y)
    lim = 1.12 * max(np.abs(np.concatenate([x, y])).max(), 1.)
    span = 3. * np.std(x, ddof=1)

    ax.grid(zorder=0)
    # The band the catalogue occupies, so the reader can see how little of the
    # frame the real signal fills
    ax.axvspan(-span, span, color=C_CAT, alpha=0.10, zorder=1, lw=0)
    ax.plot([-lim, lim], [-lim, lim], color=INK_2, lw=1.4, zorder=2)
    ax.errorbar(x, y, yerr=np.asarray(pairs['sigma']),
                xerr=np.asarray(pairs['e_cat']), fmt='o', ms=7,
                color=C_OURS, ecolor=C_OURS, elinewidth=1.3, capsize=0,
                mfc=SURFACE, mew=1.8, zorder=3)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xlabel('Teklu et al. 2025 (m/s)')
    ax.set_ylabel('this reduction (m/s)')

    ax.annotate('1:1 — what agreement\nwould look like',
                xy=(0.62 * lim, 0.62 * lim), xytext=(0.20 * lim, 0.82 * lim),
                color=INK_2, fontsize=10, ha='left',
                arrowprops=dict(arrowstyle='-', color=INK_2, lw=0.9))
    ax.text(0.0, -0.93 * lim,
            'the catalogue spans only ±{:.0f} m/s\n(shaded)'.format(span),
            color=C_CAT, fontsize=10, ha='center', va='center')
    r = np.corrcoef(x, y)[0, 1]
    ax.text(-0.96 * lim, 0.94 * lim,
            '{:d} matched epochs\nr = {:+.2f}'.format(len(x), r),
            color=INK_2, fontsize=10, va='top')
    ax.set_title('Epoch by epoch, on one scale', fontsize=12.5, pad=12)
    fig.tight_layout()
    return save(fig, 'fig_p2_compare.png')


def fig_precision(ours, orders, cat):
    """ What each stage achieves, against what the planet requires. """
    fig, ax = plt.subplots(figsize=(8.6, 4.6))

    per_order = np.median(np.asarray(ours['order_scatter'], dtype=float))
    epoch = np.std(np.asarray(ours['v_rel'], dtype=float), ddof=1)
    formal = np.median(np.asarray(ours['sigma_formal'], dtype=float))
    sel = np.asarray(cat['bjd']) < 2451110.
    modern = np.median(np.asarray(cat['e_rv'])[sel])

    labels = ['formal error\non one epoch',
              'scatter between\norders, one exposure',
              'scatter between\nepochs (achieved)',
              'modern pipeline\n(Teklu 2025)']
    vals = [formal, per_order, epoch, modern]
    colours = [C_THIRD, C_THIRD, C_OURS, C_CAT]

    ax.grid(axis='x', zorder=0)
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, height=0.52, color=colours, zorder=3)
    ax.axvline(K_PUBLISHED, color=INK, lw=1.6, ls='--', zorder=4)
    for yi, v in zip(y, vals):
        # A surface-coloured box so a value label crossing the K line stays
        # readable instead of being struck through by it.
        ax.text(v * 1.12, yi, '{:.0f} m/s'.format(v), va='center',
                color=INK_2, fontsize=10, zorder=5,
                bbox=dict(facecolor=SURFACE, edgecolor='none', pad=1.5))
    ax.text(K_PUBLISHED * 1.15, 2.62,
            'K = {:.0f} m/s\nthe planet'.format(K_PUBLISHED),
            color=INK, fontsize=10, va='center', zorder=5,
            bbox=dict(facecolor=SURFACE, edgecolor='none', pad=2.0))
    ax.set_xscale('log')
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel('velocity precision (m/s, log scale)')
    ax.set_xlim(0.7, 3000.)
    ax.set_title('What cross-correlation reaches, and what the planet needs',
                 fontsize=12.5, pad=14)
    fig.tight_layout()
    return save(fig, 'fig_p2_precision.png')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(options=None):
    parser = argparse.ArgumentParser(
        description='Assess the phase-2 velocities and make the figures.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outdir', type=str, default=FIG_DIR)
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):

    set_style()
    ours = load_ours()
    orders = Table.read(ORD_FILE)
    cat = load_catalogue()
    pairs = match(ours, cat)

    print('\n=== 1. Matching our epochs to the modern catalogue ===')
    era = np.asarray(cat['bjd']) < 2451110.
    print('  our epochs                 : {:d}'.format(len(ours)))
    print('  catalogue rows, all        : {:d}'.format(len(cat)))
    print('  catalogue rows, discovery  : {:d}'.format(int(era.sum())))
    print('  matched within {:.3f} d     : {:d}'.format(MATCH_TOL, len(pairs)))
    print('  worst time difference      : {:.4f} d ({:.1f} min)'.format(
        np.max(pairs['dt_days']), np.max(pairs['dt_days']) * 1440.))
    unmatched = [r['koaid'] for r in ours
                 if r['koaid'] not in set(pairs['koaid'])]
    print('  our epochs with no match   : {:d}'.format(len(unmatched)))
    for k in unmatched:
        print('    {:s}'.format(k))

    print('\n=== 2. Against the modern catalogue, epoch by epoch ===')
    x = np.asarray(pairs['cat'], dtype=float)
    y = np.asarray(pairs['ours'], dtype=float)
    x, y = x - np.mean(x), y - np.mean(y)
    diff = y - x
    print('  catalogue spread (rms)     : {:7.1f} m/s'.format(np.std(x, ddof=1)))
    print('  ours (rms)                 : {:7.1f} m/s'.format(np.std(y, ddof=1)))
    print('  difference (rms)           : {:7.1f} m/s'.format(np.std(diff, ddof=1)))
    r = float(np.corrcoef(x, y)[0, 1])
    # If we were measuring the catalogue's signal plus independent noise of our
    # own size, the correlation would be signal/sqrt(signal^2 + noise^2).  That
    # is the number r has to beat before it means anything.
    r_expect = np.std(x, ddof=1) / np.hypot(np.std(x, ddof=1), np.std(y, ddof=1))
    print('  correlation coefficient    : {:+7.2f}'.format(r))
    print('    expected if we measured the signal + our own noise: {:+.2f}'
          .format(r_expect))
    print('    standard error on r with {:d} points              : {:.2f}'
          .format(len(x), 1. / np.sqrt(len(x) - 3)))
    print('    -> weakly positive, about 1.5 standard errors above the')
    print('       signal-plus-noise expectation.  Not evidence of anything;'
          ' with 30\n       points this is the kind of r noise produces'
          ' routinely.')
    print('  Our scatter is {:.0f}x the catalogue\'s, so the difference is'
          ' essentially\n  our own noise: the catalogue signal is buried in'
          ' it.'.format(np.std(y, ddof=1) / np.std(x, ddof=1)))

    print('\n=== 3. Against the published 3.097-day Keplerian ===')
    tc, vc, ec = (np.asarray(cat['bjd'])[era], np.asarray(cat['rv'])[era],
                  np.asarray(cat['e_rv'])[era])
    fc = fit_circular(tc, vc, e=ec)
    print('  Catalogue, period forced to {:.4f} d:'.format(PERIOD))
    print('    K = {:6.1f} +/- {:.1f} m/s   (published {:.0f} m/s)'.format(
        fc['K'], fc['sigma_K'], K_PUBLISHED))
    print('    residual rms = {:.1f} m/s'.format(fc['rms']))
    print('    detection significance = {:.1f} sigma'.format(
        fc['K'] / fc['sigma_K']))

    fo = fit_circular(np.asarray(pairs['jd'], dtype=float),
                      np.asarray(pairs['ours'], dtype=float))
    print('\n  Ours, same period forced:')
    print('    K = {:6.1f} +/- {:.1f} m/s'.format(fo['K'], fo['sigma_K']))
    print('    residual rms = {:.1f} m/s'.format(fo['rms']))
    print('    detection significance = {:.1f} sigma'.format(
        fo['K'] / fo['sigma_K']))
    print('    95% upper limit on K   = {:.0f} m/s'.format(
        fo['K'] + 2. * fo['sigma_K']))
    print('\n  Read that fitted K carefully.  {:.0f} m/s is {:.1f}x the'
          ' published value'.format(fo['K'], fo['K'] / K_PUBLISHED))
    print('  and differs from it by {:.1f} sigma, so it is not a measurement'
          ' of the'.format(abs(fo['K'] - K_PUBLISHED) / fo['sigma_K']))
    print('  planet -- it is the amplitude our noise happens to put at this'
          ' period.')
    print('  Forcing the fit removes almost nothing: the residual rms falls'
          ' only from')
    print('  {:.0f} to {:.0f} m/s.  The honest statement is an upper limit,'
          ' and it sits\n  {:.0f}x above the signal.'.format(
              np.std(np.asarray(pairs['ours'], dtype=float), ddof=1),
              fo['rms'], (fo['K'] + 2. * fo['sigma_K']) / K_PUBLISHED))

    print('\n=== 4. What this does and does not demonstrate ===')
    print('  DOES: the reduction is sound end to end.  {:d} epochs over {:.0f}'
          ' days,'.format(len(ours),
                          np.ptp(np.asarray(ours['jd_mid'], dtype=float))))
    print('        every one yielding a velocity, with the three template'
          ' exposures')
    print('        returning zero to 17 m/s and the orders of a single'
          ' exposure')
    print('        agreeing to {:.0f} m/s.'.format(
        np.median(np.asarray(ours['order_scatter'], dtype=float))))
    print('  DOES NOT: detect the planet, or confirm the published orbit.'
          '  At {:.0f} m/s'.format(np.std(y, ddof=1)))
    print('        against a {:.0f} m/s semiamplitude, a non-detection was'
          ' the expected'.format(K_PUBLISHED))
    print('        outcome and is what we got.  Nothing here would have'
          ' changed if the')
    print('        planet were not there.')

    print('\n=== 5. What it implies for phase 3 ===')
    need = np.std(y, ddof=1) / (K_PUBLISHED / 3.)
    print('  To detect K = {:.0f} m/s at 3 sigma per-epoch precision must'
          ' improve by'.format(K_PUBLISHED))
    print('  about {:.0f}x, from {:.0f} m/s to ~{:.0f} m/s.'.format(
        need, np.std(y, ddof=1), K_PUBLISHED / 3.))
    print('  The per-order scatter is already {:.0f} m/s, so the photons are'
          ' not the'.format(np.median(np.asarray(
              ours['order_scatter'], dtype=float))))
    print('  limit -- the wavelength zero point is (prompt 5).  That is'
          ' exactly what')
    print('  the iodine cell fixes, which is the case for phase 3.')
    print('  Teklu et al. reach {:.1f} m/s on these same frames, so the data'
          ' themselves'.format(np.median(ec)))
    print('  carry the signal; only our wavelength calibration does not.')

    print('\n=== 6. Figures ===')
    fig_timeseries(pairs, cat, ours)
    fig_phasefold(pairs, cat)
    fig_compare(pairs)
    fig_precision(ours, orders, cat)


if __name__ == '__main__':
    main(parse_args())
