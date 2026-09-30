""" Assessment and figures for the phase-3 iodine velocities.

Prompt 8 of claude_prompts/data_phase3_prompt.md: assess what prompt 7
produced (``data/iodine_velocities.csv``) against the published 3.097-day
Keplerian and the modern catalogue, *exactly as phase 2 prompt 6 did* so the
two are comparable, and say what the residuals contain.

Comparability is kept by reusing `figs_phase2` itself: the same catalogue
loader, the same discovery-era cut, the same time matching (``MATCH_TOL``),
the same forced-period circular fit (`fit_circular`, linear, period
3.0965828 d), the same unweighted-and-rescaled K error, the same correlation
test and the same house style.  What phase 3 adds, each labelled as such:

  * the same fit weighted with prompt 7's ``sig_epoch``;
  * a period search that is not told the period (Lomb-Scargle);
  * the amplitude with the orbital phase fixed to the catalogue's;
  * a leave-one-night-out check of K;
  * what the residuals correlate with.

    fig_ph3_timeseries   ours and the catalogue over the nine months, on the
                         same scale this time
    fig_ph3_phasefold    both folded on the period, with residuals
    fig_ph3_compare      epoch by epoch, equal aspect
    fig_ph3_precision    each stage against what the planet requires, with
                         phase 2 beside it
    fig_ph3_residuals    what the residuals contain, and the period search

The figure names avoid ``fig_p3_*``, which prompt 3 already uses.

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.figs_phase3

"""

# Standard imports
import os
import argparse

import numpy as np
from scipy import stats

from astropy.table import Table

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from first_hires_exoplanet import figs_phase2 as f2
from first_hires_exoplanet.figs_phase2 import (C_OURS, C_CAT, C_THIRD, INK, INK_2,
                                               SURFACE, PERIOD, K_PUBLISHED)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR = f2.DATA_DIR
FIG_DIR = f2.FIG_DIR
VEL_FILE = os.path.join(DATA_DIR, 'iodine_velocities.csv')
ORD_FILE = os.path.join(DATA_DIR, 'iodine_velocities_per_order.csv')
SCI_FILE = os.path.join(DATA_DIR, 'koa_hd187123_science.csv')
C_PHASE2 = '#8a8a8a'

#: Discovery-era cut, as phase 2
ERA_END = 2451110.

#: Period-search grid (days)
PMIN, PMAX = 0.5, 30.


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load_ours(path=VEL_FILE):
    """ Prompt 7's velocities, with the columns `figs_phase2.match` reads. """
    t = Table.read(path)
    t['jd_mid'] = np.asarray(t['jd'], dtype=float)      # JD(UTC), mid-exposure
    t['v_rel'] = np.asarray(t['v_bary'], dtype=float)
    t['night'] = [str(n) for n in t['night']]
    t['sigma_empirical'] = np.asarray(t['sig_epoch'], dtype=float)
    sci = Table.read(SCI_FILE)
    key = {str(k).replace('.fits', ''): i for i, k in enumerate(sci['koaid'])}
    idx = [key[str(k)] for k in t['koaid']]
    t['airmass'] = np.asarray(sci['airmass'][idx], dtype=float)
    t['exptime'] = np.asarray(sci['exptime'][idx], dtype=float)
    return t


def phase2_velocities():
    """ Phase 2's cross-correlation velocities with the same columns. """
    return f2.load_ours()


# ---------------------------------------------------------------------------
# Statistics beyond phase 2
# ---------------------------------------------------------------------------

def fit_fixed_phase(t, v, phase, e=None, period=PERIOD):
    """ Amplitude and offset with the orbital phase held at the catalogue's.

    v = gamma + K cos(2 pi t / P - phase): linear in gamma and K.
    """
    x = np.cos(2. * np.pi * t / period - phase)
    M = np.column_stack([np.ones_like(t), x])
    w = np.ones_like(t) if e is None else 1. / np.asarray(e) ** 2
    cov = np.linalg.inv(M.T @ (w[:, None] * M))
    p = cov @ (M.T @ (w * v))
    resid = v - M @ p
    chi2_red = float(np.sum(w * resid ** 2) / max(len(t) - 2, 1))
    sig = float(np.sqrt(cov[1, 1]))
    if e is None:
        sig *= np.sqrt(chi2_red)
    return dict(K=float(p[1]), sigma_K=sig, gamma=float(p[0]), chi2_red=chi2_red,
                rms=float(np.std(resid, ddof=1)))


def weighted_fit(t, v, e):
    """ `figs_phase2.fit_circular` weighted, K error rescaled only if chi2 > 1. """
    f = f2.fit_circular(t, v, e=e)
    f['sigma_K_scaled'] = f['sigma_K'] * np.sqrt(max(f['chi2_red'], 1.))
    return f


def delta_chi2(t, v, e, fit):
    """ Constant model against the circular Keplerian: F-test, 2 extra parameters. """
    w = 1. / np.asarray(e) ** 2
    c = np.sum(w * v) / np.sum(w)
    chi2_0 = float(np.sum(w * (v - c) ** 2))
    chi2_1 = float(np.sum(w * (v - fit['model'](t)) ** 2))
    n = len(t)
    F = ((chi2_0 - chi2_1) / 2.) / (chi2_1 / (n - 3))
    p = float(stats.f.sf(F, 2, n - 3))
    return dict(chi2_0=chi2_0, chi2_1=chi2_1, F=float(F), p=p,
                nsig=float(stats.norm.isf(p / 2.)) if p > 0 else np.inf)


def period_search(t, v, e):
    """ Lomb-Scargle over PMIN-PMAX days, not told the period. """
    from astropy.timeseries import LombScargle
    ls = LombScargle(t, v, e)
    freq = np.linspace(1. / PMAX, 1. / PMIN, 60000)
    power = ls.power(freq)
    k = int(np.argmax(power))
    p_true = ls.power(np.array([1. / PERIOD]))[0]
    fap = float(ls.false_alarm_probability(power[k], minimum_frequency=freq[0],
                                           maximum_frequency=freq[-1]))
    # The next-highest peak that is not the best one or its immediate wings
    far = np.abs(freq - freq[k]) > 0.02
    k2 = int(np.argmax(np.where(far, power, 0.)))
    return dict(freq=freq, power=power, best_period=float(1. / freq[k]),
                best_power=float(power[k]), power_at_true=float(p_true), fap=fap,
                second_period=float(1. / freq[k2]), second_power=float(power[k2]))


def leave_one_night_out(t, v, e, night):
    rows = []
    for n in np.unique(night):
        s = night != n
        f = f2.fit_circular(t[s], v[s], e=e[s])
        rows.append(dict(night=n, K=f['K'], nfr=int(np.sum(~s))))
    return Table(rows)


def residual_drivers(ours, resid):
    """ Spearman correlation of the residuals (and |residual|) with each epoch quantity. """
    cols = [('dv_wave', 'iodine - PypeIt wavelength zero point'),
            ('n_orders', 'orders fitted'),
            ('sig_photon', 'photon-noise error'),
            ('redchi_meas', 'model misfit, chi2/dof'),
            ('bvc', 'barycentric correction'),
            ('airmass', 'airmass'),
            ('exptime', 'exposure time'),
            ('jd', 'time')]
    rows = []
    for c, label in cols:
        x = np.asarray(ours[c], dtype=float)
        ok = np.isfinite(x) & np.isfinite(resid)
        r, p = stats.spearmanr(x[ok], resid[ok])
        ra, pa = stats.spearmanr(x[ok], np.abs(resid[ok]))
        rows.append(dict(quantity=c, label=label, rho=float(r), p=float(p),
                         rho_abs=float(ra), p_abs=float(pa)))
    return Table(rows)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig_timeseries(ours, cat, fcat):
    """ The two velocity sets over the same nine months, on the same scale. """
    fig, axs = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True, sharey=True)
    t0 = 2450800.
    sel = np.asarray(cat['bjd']) < ERA_END
    tt = np.linspace(2450800., ERA_END, 20000)
    for ax, (t, v, e, col, name, note) in zip(axs, [
            (np.asarray(ours['jd_mid']), np.asarray(ours['v_rel']), np.asarray(ours['sig_epoch']),
             C_OURS, 'this reduction, iodine forward model (phase 3)',
             'epoch-to-epoch rms {:.0f} m/s; median error {:.0f} m/s'.format(
                 np.std(ours['v_rel'], ddof=1), np.median(ours['sig_epoch']))),
            (np.asarray(cat['bjd'])[sel], np.asarray(cat['rv'])[sel],
             np.asarray(cat['e_rv'])[sel], C_CAT, 'Teklu et al. 2025, the same frames',
             'median error bar {:.1f} m/s'.format(np.median(np.asarray(cat['e_rv'])[sel])))]):
        ax.grid(axis='y', zorder=0)
        vv = v - np.mean(v)
        ax.errorbar(t - t0, vv, yerr=e, fmt='o', ms=6, color=col, ecolor=col,
                    elinewidth=1.4, capsize=0, mfc=SURFACE, mew=1.8, zorder=3)
        ax.axhline(0., color=INK_2, lw=0.8, ls=':', zorder=1)
        ax.set_ylabel('relative velocity (m/s)')
        ax.text(0.012, 0.95, name, transform=ax.transAxes, va='top', color=col, fontsize=11)
        ax.text(0.012, 0.845, note, transform=ax.transAxes, va='top', color=INK_2, fontsize=10)
    axs[1].set_xlabel('JD - 2450800 (days)')
    axs[0].set_ylim(-190, 230)
    fig.suptitle('The same photons, measured two ways, on one scale', y=0.985, fontsize=12.5)
    fig.tight_layout()
    return f2.save(fig, 'fig_ph3_timeseries.png')


def fig_phasefold(ours, cat, fcat):
    fig, axs = plt.subplots(2, 2, figsize=(11, 6.2), sharex=True,
                            gridspec_kw=dict(height_ratios=[3, 1.3]))
    sel = np.asarray(cat['bjd']) < ERA_END
    grid = np.linspace(0., 1., 400)
    curve = fcat['model'](grid * PERIOD) - fcat['gamma']
    for j, (t, v, e, col, label) in enumerate([
            (np.asarray(cat['bjd'])[sel], np.asarray(cat['rv'])[sel],
             np.asarray(cat['e_rv'])[sel], C_CAT, 'Teklu et al. 2025'),
            (np.asarray(ours['jd_mid']), np.asarray(ours['v_rel']),
             np.asarray(ours['sig_epoch']), C_OURS, 'this reduction (phase 3)')]):
        f = f2.fit_circular(t, v, e=e)
        ph = (t / PERIOD) % 1.
        vv = v - f['gamma']
        ax = axs[0, j]
        ax.grid(axis='y', zorder=0)
        ax.plot(grid, curve, color=INK_2, lw=1.6, zorder=2,
                label='catalogue Keplerian, K = {:.0f} m/s'.format(fcat['K']))
        ax.errorbar(ph, vv, yerr=e, fmt='o', ms=6, color=col, ecolor=col, elinewidth=1.4,
                    capsize=0, mfc=SURFACE, mew=1.8, zorder=3,
                    label='{:s}: K = {:.0f} ± {:.0f} m/s'.format(label, f['K'], f['sigma_K']))
        ax.set_ylim(-170, 200)
        ax.legend(loc='upper right', fontsize=9)
        ax = axs[1, j]
        ax.grid(axis='y', zorder=0)
        r = v - fcat['model'](t)
        r -= np.average(r, weights=1. / e ** 2)
        ax.errorbar(ph, r, yerr=e, fmt='o', ms=4, color=col, ecolor=col, elinewidth=1.2,
                    capsize=0, mfc=SURFACE, mew=1.4, zorder=3)
        ax.axhline(0., color=INK_2, lw=0.8)
        ax.set_ylim(-110, 110)
        ax.text(0.02, 0.92, 'minus the catalogue Keplerian: rms {:.1f} m/s'.format(
            np.std(r, ddof=1)), transform=ax.transAxes, va='top', color=INK_2, fontsize=9)
        ax.set_xlabel('orbital phase (P = {:.4f} d)'.format(PERIOD))
    axs[0, 0].set_ylabel('velocity (m/s)')
    axs[1, 0].set_ylabel('residual (m/s)')
    fig.suptitle('Folded on the published period', y=0.99, fontsize=12.5)
    fig.tight_layout()
    return f2.save(fig, 'fig_ph3_phasefold.png')


def fig_compare(pairs, pairs2):
    """ Ours against the catalogue, epoch by epoch, equal aspect; phase 2 in grey. """
    fig, ax = plt.subplots(figsize=(6.6, 6.0))
    x = np.asarray(pairs['cat'], dtype=float)
    y = np.asarray(pairs['ours'], dtype=float)
    x, y = x - np.mean(x), y - np.mean(y)
    lim = 200.
    ax.grid(zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color=INK_2, lw=1.4, zorder=2)
    # Phase 2 on the same axes: almost all of it falls outside the frame
    x2 = np.asarray(pairs2['cat'], dtype=float)
    y2 = np.asarray(pairs2['ours'], dtype=float)
    x2, y2 = x2 - np.mean(x2), y2 - np.mean(y2)
    inside = np.abs(y2) < lim
    ax.plot(x2[inside], y2[inside], 's', ms=5, color=C_PHASE2, alpha=0.6, zorder=2,
            label='phase 2, cross-correlation ({:d} of {:d} inside the frame)'.format(
                int(inside.sum()), len(y2)))
    ax.errorbar(x, y, yerr=np.asarray(pairs['sigma']), xerr=np.asarray(pairs['e_cat']),
                fmt='o', ms=7, color=C_OURS, ecolor=C_OURS, elinewidth=1.3, capsize=0,
                mfc=SURFACE, mew=1.8, zorder=3, label='phase 3, iodine forward model')
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xlabel('Teklu et al. 2025 (m/s)')
    ax.set_ylabel('this reduction (m/s)')
    r = np.corrcoef(x, y)[0, 1]
    ax.text(-0.96 * lim, 0.94 * lim, '{:d} matched epochs\nr = {:+.2f}'.format(len(x), r),
            color=INK_2, fontsize=10, va='top')
    ax.legend(loc='lower right', fontsize=8.5)
    ax.set_title('Epoch by epoch, on one scale', fontsize=12.5, pad=12)
    fig.tight_layout()
    return f2.save(fig, 'fig_ph3_compare.png')


def fig_precision(ours, resid_rms, phase2_rms, modern):
    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    labels = ['photon noise,\none epoch',
              'adopted error,\none epoch',
              'scatter between\norders, one exposure',
              'residual about the\nKeplerian (achieved)',
              'phase 2: scatter\nbetween epochs',
              'modern pipeline\n(Teklu 2025)']
    vals = [np.median(ours['sig_photon']), np.median(ours['sig_epoch']),
            np.median(ours['order_scatter']), resid_rms, phase2_rms, modern]
    colours = [C_THIRD, C_THIRD, C_THIRD, C_OURS, C_PHASE2, C_CAT]
    ax.grid(axis='x', zorder=0)
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, height=0.52, color=colours, zorder=3)
    ax.axvline(K_PUBLISHED, color=INK, lw=1.6, ls='--', zorder=4)
    ax.axvline(50., color=INK_2, lw=1.0, ls=':', zorder=4)
    for yi, v in zip(y, vals):
        ax.text(v * 1.12, yi, '{:.0f} m/s'.format(v) if v >= 2 else '{:.1f} m/s'.format(v),
                va='center', color=INK_2, fontsize=10, zorder=5,
                bbox=dict(facecolor=SURFACE, edgecolor='none', pad=1.5))
    ax.text(K_PUBLISHED * 1.12, y[0] + 0.45, 'K = {:.0f} m/s'.format(K_PUBLISHED), color=INK,
            fontsize=10, va='center', zorder=5,
            bbox=dict(facecolor=SURFACE, edgecolor='none', pad=2.0))
    ax.text(50. / 1.12, y[0] + 0.45, 'phase-3 target\n50 m/s', color=INK_2, fontsize=9,
            va='center', ha='right', zorder=5,
            bbox=dict(facecolor=SURFACE, edgecolor='none', pad=2.0))
    ax.set_xscale('log')
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel('velocity precision (m/s, log scale)')
    ax.set_xlim(0.7, 3000.)
    ax.set_title('What the iodine forward model reaches, and what the planet needs',
                 fontsize=12.5, pad=14)
    fig.tight_layout()
    return f2.save(fig, 'fig_ph3_precision.png')


def fig_residuals(ours, resid, ps, ps_res, drivers):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.2))
    e = np.asarray(ours['sig_epoch'])
    t0 = 2450800.

    ax = axs[0, 0]
    ax.grid(axis='y', zorder=0)
    ax.errorbar(np.asarray(ours['jd_mid']) - t0, resid, yerr=e, fmt='o', ms=5, color=C_OURS,
                ecolor=C_OURS, elinewidth=1.2, capsize=0, mfc=SURFACE, mew=1.6, zorder=3)
    ax.axhline(0., color=INK_2, lw=0.8)
    ax.set_xlabel('JD - 2450800 (days)')
    ax.set_ylabel('residual about the catalogue Keplerian (m/s)')
    ax.set_title('(a) residuals in time: chi2/dof {:.2f}'.format(
        np.sum((resid / e) ** 2) / (len(e) - 1)), loc='left', fontsize=10)

    ax = axs[0, 1]
    ax.grid(axis='y', zorder=0)
    x = np.asarray(ours['dv_wave'], dtype=float)
    ax.errorbar(x, resid, yerr=e, fmt='o', ms=5, color=C_OURS, ecolor=C_OURS, elinewidth=1.2,
                capsize=0, mfc=SURFACE, mew=1.6, zorder=3)
    ax.axhline(0., color=INK_2, lw=0.8)
    d = drivers[drivers['quantity'] == 'dv_wave'][0]
    ax.set_xlabel('PypeIt - iodine wavelength zero point (m/s)')
    ax.set_ylabel('residual (m/s)')
    ax.set_title('(b) against the ThAr error phase 2 could not remove: '
                 'Spearman {:+.2f}'.format(d['rho']), loc='left', fontsize=10)

    ax = axs[1, 0]
    ax.grid(axis='y', zorder=0)
    ax.errorbar(np.asarray(ours['n_orders']) + np.random.default_rng(1).uniform(
        -0.15, 0.15, len(resid)), resid, yerr=e, fmt='o', ms=5, color=C_OURS, ecolor=C_OURS,
        elinewidth=1.2, capsize=0, mfc=SURFACE, mew=1.6, zorder=3)
    ax.axhline(0., color=INK_2, lw=0.8)
    ax.set_xlabel('iodine orders fitted')
    ax.set_ylabel('residual (m/s)')
    ax.set_title('(c) against the number of usable orders', loc='left', fontsize=10)

    ax = axs[1, 1]
    ax.plot(1. / ps['freq'], ps['power'], color=C_OURS, lw=1.0,
            label='velocities: best {:.4f} d, FAP {:.1e}'.format(ps['best_period'], ps['fap']))
    ax.plot(1. / ps_res['freq'], ps_res['power'], color=C_THIRD, lw=0.9, alpha=0.9,
            label='residuals: best {:.2f} d, FAP {:.2f}'.format(ps_res['best_period'],
                                                              ps_res['fap']))
    ax.axvline(PERIOD, color=INK, lw=1.0, ls='--')
    ax.set_xscale('log')
    ax.set_xlabel('period (days)')
    ax.set_ylabel('Lomb-Scargle power')
    ax.legend(fontsize=8.5, loc='upper left')
    ax.set_title('(d) period search, not told the period (dashed: 3.0966 d)', loc='left',
                 fontsize=10)

    fig.tight_layout()
    return f2.save(fig, 'fig_ph3_residuals.png')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(options=None):
    parser = argparse.ArgumentParser(
        description='Assess the phase-3 velocities and make the figures.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outdir', type=str, default=FIG_DIR)
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):
    f2.set_style()
    ours = load_ours()
    xc = phase2_velocities()
    cat = f2.load_catalogue()
    era = np.asarray(cat['bjd']) < ERA_END
    pairs = f2.match(ours, cat)
    pairs2 = f2.match(xc, cat)

    print('\n=== 1. Matching our epochs to the modern catalogue ===')
    print('  our epochs                 : {:d}'.format(len(ours)))
    print('  catalogue rows, discovery  : {:d}'.format(int(era.sum())))
    print('  matched within {:.3f} d     : {:d}'.format(f2.MATCH_TOL, len(pairs)))
    print('  worst time difference      : {:.1f} min'.format(np.max(pairs['dt_days']) * 1440.))
    unmatched = [r['koaid'] for r in ours if r['koaid'] not in set(pairs['koaid'])]
    for k in unmatched:
        r = ours[list(ours['koaid']).index(k)]
        print('  no catalogue row: {:s} ({:.0f} s exposure, {:d} orders, v = {:+.1f} +- {:.1f})'.format(
            k, r['exptime'], r['n_orders'], r['v_rel'], r['sig_epoch']))

    print('\n=== 2. Against the modern catalogue, epoch by epoch ===')
    x = np.asarray(pairs['cat'], dtype=float)
    y = np.asarray(pairs['ours'], dtype=float)
    x, y = x - np.mean(x), y - np.mean(y)
    diff = y - x
    print('  catalogue spread (rms)     : {:7.1f} m/s'.format(np.std(x, ddof=1)))
    print('  ours (rms)                 : {:7.1f} m/s'.format(np.std(y, ddof=1)))
    print('  difference (rms)           : {:7.1f} m/s'.format(np.std(diff, ddof=1)))
    print('  difference / our error     : chi2/dof {:.2f}'.format(
        np.sum((diff / np.hypot(pairs['sigma'], pairs['e_cat'])) ** 2) / (len(diff) - 1)))
    r = float(np.corrcoef(x, y)[0, 1])
    noise = np.std(diff, ddof=1)
    r_expect = np.std(x, ddof=1) / np.hypot(np.std(x, ddof=1), noise)
    print('  correlation coefficient    : {:+7.2f}   (phase 2: {:+.2f})'.format(
        r, float(np.corrcoef(pairs2['cat'], pairs2['ours'])[0, 1])))
    print('    expected if we measured the signal + our own noise: {:+.2f}'.format(r_expect))
    print('    standard error with {:d} points: {:.2f}; r = 0 rejected at {:.1f} sigma'.format(
        len(x), 1. / np.sqrt(len(x) - 3), np.arctanh(r) * np.sqrt(len(x) - 3)))
    slope = np.polyfit(x, y, 1)[0]
    print('  slope, ours on catalogue   : {:.2f}'.format(slope))

    print('\n=== 3. Against the published 3.097-day Keplerian ===')
    tc, vc, ec = (np.asarray(cat['bjd'])[era], np.asarray(cat['rv'])[era],
                  np.asarray(cat['e_rv'])[era])
    fc = f2.fit_circular(tc, vc, e=ec)
    print('  Catalogue, period forced: K = {:.1f} +- {:.1f} m/s, rms {:.1f}, {:.0f} sigma'.format(
        fc['K'], fc['sigma_K'], fc['rms'], fc['K'] / fc['sigma_K']))
    t = np.asarray(ours['jd_mid'], dtype=float)
    v = np.asarray(ours['v_rel'], dtype=float)
    e = np.asarray(ours['sig_epoch'], dtype=float)
    fo = f2.fit_circular(t, v)
    print('\n  Ours, exactly as phase 2 (unweighted, error rescaled by the scatter):')
    print('    K = {:6.1f} +/- {:.1f} m/s, residual rms {:.1f} m/s, {:.1f} sigma'.format(
        fo['K'], fo['sigma_K'], fo['rms'], fo['K'] / fo['sigma_K']))
    fo2 = f2.fit_circular(np.asarray(pairs['jd']), np.asarray(pairs['ours']))
    print('    the {:d} matched epochs only: K = {:.1f} +/- {:.1f}'.format(
        len(pairs), fo2['K'], fo2['sigma_K']))
    fx = f2.fit_circular(np.asarray(pairs2['jd']), np.asarray(pairs2['ours']))
    print('    phase 2 for comparison: K = {:.0f} +/- {:.0f} m/s, rms {:.0f}'.format(
        fx['K'], fx['sigma_K'], fx['rms']))
    fw = weighted_fit(t, v, e)
    print('  Ours, weighted with sig_epoch (phase 3 addition):')
    print('    K = {:6.1f} +/- {:.1f} m/s (chi2/dof {:.2f}; rescaled +/- {:.1f}), {:.1f} sigma'.format(
        fw['K'], fw['sigma_K'], fw['chi2_red'], fw['sigma_K_scaled'],
        fw['K'] / fw['sigma_K_scaled']))
    dphase = np.degrees(np.angle(np.exp(1j * (fw['phase'] - fc['phase']))))
    ph_err = np.degrees(fw['sigma_K_scaled'] / fw['K'])
    print('    orbital phase against the catalogue: {:+.1f} deg (+- ~{:.0f} deg)'.format(dphase, ph_err))
    ff = fit_fixed_phase(t, v, fc['phase'], e=e)
    print('  Phase fixed to the catalogue: K = {:.1f} +/- {:.1f} m/s ({:.1f} sigma)'.format(
        ff['K'], ff['sigma_K'] * np.sqrt(max(ff['chi2_red'], 1.)),
        ff['K'] / (ff['sigma_K'] * np.sqrt(max(ff['chi2_red'], 1.)))))
    dc = delta_chi2(t, v, e, fw)
    print('  Constant vs Keplerian (weighted): chi2 {:.1f} -> {:.1f}; F = {:.1f}, p = {:.1e} '
          '({:.1f} sigma)'.format(dc['chi2_0'], dc['chi2_1'], dc['F'], dc['p'], dc['nsig']))
    print('  K against the published {:.0f} m/s: {:+.1f} sigma; against the catalogue {:.1f}: '
          '{:+.1f} sigma'.format(K_PUBLISHED, (fw['K'] - K_PUBLISHED) / fw['sigma_K_scaled'],
                                 fc['K'], (fw['K'] - fc['K']) / fw['sigma_K_scaled']))

    ps = period_search(t, v, e)
    print('\n  Period search {:.1f}-{:.0f} d, not told the period:'.format(PMIN, PMAX))
    print('    best period {:.4f} d, power {:.3f}, false-alarm probability {:.1e}'.format(
        ps['best_period'], ps['best_power'], ps['fap']))
    print('    next peak {:.4f} d, power {:.3f}; power at 3.0966 d {:.3f}'.format(
        ps['second_period'], ps['second_power'], ps['power_at_true']))

    lno = leave_one_night_out(t, v, e, np.asarray(ours['night']))
    print('  Leave one night out: K from {:.1f} to {:.1f} m/s; the night that moves it most: {:s} '
          '({:d} frames, K {:.1f})'.format(lno['K'].min(), lno['K'].max(),
                                           lno['night'][np.argmax(np.abs(lno['K'] - fw['K']))],
                                           lno['nfr'][np.argmax(np.abs(lno['K'] - fw['K']))],
                                           lno['K'][np.argmax(np.abs(lno['K'] - fw['K']))]))

    print('\n=== 4. What the residuals contain ===')
    resid = v - fc['model'](t)
    resid -= np.average(resid, weights=1. / e ** 2)
    resid_own = v - fw['model'](t)
    print('  about the catalogue Keplerian (offset only): rms {:.1f} m/s, chi2/dof {:.2f}'.format(
        np.std(resid, ddof=1), np.sum((resid / e) ** 2) / (len(e) - 1)))
    print('  about our own fit: rms {:.1f} m/s, chi2/dof {:.2f}'.format(
        np.std(resid_own, ddof=1), np.sum((resid_own / e) ** 2) / (len(e) - 3)))
    print('  robust rms {:.1f}; largest |residual|: {:s}'.format(
        f2_robust(resid), ', '.join('{:s} {:+.0f} ({:.1f} sigma)'.format(
            ours['koaid'][i], resid[i], resid[i] / e[i]) for i in np.argsort(-np.abs(resid))[:4])))
    # Within-night and between-night parts
    night = np.asarray(ours['night'])
    nights = np.unique(night)
    nm = np.array([np.average(resid[night == n], weights=1. / e[night == n] ** 2) for n in nights])
    within = np.concatenate([resid[night == n] - m for n, m in zip(nights, nm)
                             if np.sum(night == n) > 1])
    nmulti = sum(np.sum(night == n) > 1 for n in nights)
    print('  within a night ({:d} nights, {:d} frames): rms about the nightly mean {:.1f} m/s'.format(
        nmulti, len(within), np.sqrt(np.sum(within ** 2) / (len(within) - nmulti))))
    print('  between nights ({:d} nightly means): rms {:.1f} m/s'.format(len(nm), np.std(nm, ddof=1)))
    for n in nights:
        s = night == n
        if s.sum() > 2:
            print('    {:s}: {:s}'.format(n, ', '.join('{:+.0f}'.format(q) for q in resid[s])))
    drivers = residual_drivers(ours, resid)
    print('  Spearman correlation of the residual (and of |residual|) with:')
    for d in drivers:
        print('    {:<42s} {:+.2f} (p {:.2f})   |r|: {:+.2f} (p {:.2f})'.format(
            d['label'], d['rho'], d['p'], d['rho_abs'], d['p_abs']))
    # Borrowed-arc nights
    for n in ('19980719', '19980917'):
        s = night == n
        print('  borrowed-arc night {:s}: residuals {:s} m/s (phase 2 put it {:s} off)'.format(
            n, ', '.join('{:+.0f}'.format(q) for q in resid[s]),
            ', '.join('{:+.0f}'.format(q) for q in np.asarray(ours['v_xcorr'][s]))))
    ps_res = period_search(t, resid, e)
    print('  residual period search: best {:.3f} d, FAP {:.2f}'.format(ps_res['best_period'],
                                                                      ps_res['fap']))
    # Against the catalogue directly (their 1.2 m/s is effectively truth)
    dmatch = diff - np.mean(diff)
    print('  ours - catalogue, matched: rms {:.1f} m/s; that is our error, at chi2/dof {:.2f} '
          'against sig_epoch'.format(np.std(dmatch, ddof=1), np.sum(
              (dmatch / np.asarray(pairs['sigma'])) ** 2) / (len(dmatch) - 1)))

    print('\n=== 5. Achieved precision and what it means ===')
    ach = float(np.std(resid, ddof=1))
    print('  per-epoch precision achieved : {:.1f} m/s (residual about the Keplerian)'.format(ach))
    print('  predicted sigma_K at that precision, N = {:d}: {:.1f} m/s -> {:.1f} sigma for K = {:.0f}'.format(
        len(t), ach * np.sqrt(2. / len(t)), K_PUBLISHED / (ach * np.sqrt(2. / len(t))), K_PUBLISHED))
    print('  against phase 2 ({:.0f} m/s): {:.0f}x better; against the goal (50 m/s): {:.1f}x; '
          'against the catalogue ({:.1f} m/s): {:.0f}x short'.format(
              np.std(np.asarray(pairs2['ours']), ddof=1),
              np.std(np.asarray(pairs2['ours']), ddof=1) / ach, 50. / ach, np.median(ec),
              ach / np.median(ec)))

    print('\n=== 6. Figures ===')
    fig_timeseries(ours, cat, fc)
    fig_phasefold(ours, cat, fc)
    fig_compare(pairs, pairs2)
    fig_precision(ours, ach, np.std(np.asarray(pairs2['ours']), ddof=1), np.median(ec))
    fig_residuals(ours, resid, ps, ps_res, drivers)

    # The residual table, beside the velocities
    out = Table()
    out['koaid'] = ours['koaid']
    out['jd'] = ours['jd_mid']
    out['v_bary'] = ours['v_rel']
    out['sig_epoch'] = ours['sig_epoch']
    out['v_keplerian'] = fc['model'](t) + np.average(v - fc['model'](t), weights=1. / e ** 2)
    out['resid'] = resid
    cmap = {k: (c, ec_) for k, c, ec_ in zip(pairs['koaid'], pairs['cat'], pairs['e_cat'])}
    out['v_catalogue'] = [cmap[k][0] if k in cmap else np.nan for k in ours['koaid']]
    out['e_catalogue'] = [cmap[k][1] if k in cmap else np.nan for k in ours['koaid']]
    for c in out.colnames[1:]:
        out[c].info.format = '.4f' if c == 'jd' else '.2f'
    out.write(os.path.join(DATA_DIR, 'iodine_keplerian_residuals.csv'), overwrite=True)
    print('  wrote {:s}'.format(os.path.join(DATA_DIR, 'iodine_keplerian_residuals.csv')))


def f2_robust(x):
    from first_hires_exoplanet.pyodine_chunk_stats import robust_sigma
    return float(robust_sigma(np.asarray(x)))


if __name__ == '__main__':
    main(parse_args())
