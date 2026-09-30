""" Diagnose what is left after the iodine forward model (phase 3, prompt 9).

The prompt's questions, each answered from prompt 7's velocities and prompt
8's Keplerian residuals, beside phase 2's cross-correlation on the same
frames:

  1. THE BORROWED-ARC NIGHTS.  1998-07-19 (the 07-18 arc) and 1998-09-17 (the
     09-15 B2 arc).  Do they come back into line?  The forward model fits its
     own wavelength solution, so it also *measures* PypeIt's: ``dv_wave`` is
     the iodine-fitted wavelength zero point minus PypeIt's, per chunk, and
     ``zp = -dv_wave`` is PypeIt's error in velocity units.
  2. THE INTRA-NIGHT DRIFTS.  1998-08-25 (four frames over 6.7 h) and 08-26
     (three HD 187123 frames, plus ten B-star frames through the cell).  The
     slope with time of phase 2's residual, of PypeIt's zero-point error, and
     of phase 3's residual.
  3. THE DECEMBER NIGHTS.  1997-12-23 and 12-24, reduced with no pixel flat
     (trace from an iodine-in flat).  12-24's one frame is cell-out, so it has
     no iodine velocity; 12-23's is compared with every other epoch on the
     quantities a flat-field error would move.
  4. 1998-09-13.  Why 34 orders, and what the iodine fit makes of the epoch.
  5. WHAT IS LEFT.  The seasonal correlation prompt 8 flagged, tested within
     nights where BVC changes but the season does not; and what the residual
     budget becomes when the weakest orders and epochs are left out
     (re-combined with `pyodine`'s own combination, as prompt 7).

Writes ``data/diagnose_nights.csv``, ``data/diagnose_variants.csv`` and
``docs/figs/fig_p9_diagnosis.png``.

Run with:

    conda run -n pypeit14 python -m first_hires_exoplanet.diagnose_phase3

"""

# Must precede any import of matplotlib, including pyodine's own
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# Standard imports
import re
import glob
import json
import argparse

import numpy as np
from scipy import stats
from astropy.io import fits
from astropy.table import Table

from first_hires_exoplanet import figs_phase2 as f2
from first_hires_exoplanet import fit_all_epochs as fa
from first_hires_exoplanet import reduce_run as rr
from first_hires_exoplanet.pyodine_chunk_stats import robust_sigma


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR = f2.DATA_DIR
FIG_DIR = f2.FIG_DIR
BORROWED = {'19980719': '1998-07-18 arc', '19980917': '1998-09-15 B2 arc'}
NO_PIXFLAT = ('19971223', '19971224')
RED_ORDERS = (58, 59, 60, 61)
C_OURS, C_CAT, C_THIRD, INK, INK_2 = f2.C_OURS, f2.C_CAT, f2.C_THIRD, f2.INK, f2.INK_2
C_PHASE2 = '#8a8a8a'


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load():
    vt = Table.read(os.path.join(DATA_DIR, 'iodine_velocities.csv'))
    vt['night'] = [str(n) for n in vt['night']]
    kr = Table.read(os.path.join(DATA_DIR, 'iodine_keplerian_residuals.csv'))
    assert list(kr['koaid']) == list(vt['koaid'])
    vt['resid'] = kr['resid']
    cat = f2.load_catalogue()
    era = np.asarray(cat['bjd']) < 2451110.
    orb = f2.fit_circular(np.asarray(cat['bjd'])[era], np.asarray(cat['rv'])[era],
                          e=np.asarray(cat['e_rv'])[era])
    t = np.asarray(vt['jd'], dtype=float)
    r2 = np.asarray(vt['v_xcorr'], dtype=float) - orb['model'](t)
    vt['resid_p2'] = r2 - np.mean(r2)
    vt['zp'] = -np.asarray(vt['dv_wave'], dtype=float)      # PypeIt - iodine, m/s
    vt['ut_h'] = [int(str(k).split('.')[2]) / 3600. for k in vt['koaid']]
    ch = Table.read(os.path.join(DATA_DIR, 'iodine_chunks.csv'))
    return vt, ch, orb


def fit_line(x, y, e=None):
    """ Weighted straight line; slope and its error (rescaled if chi2 > 1). """
    x, y = np.asarray(x, float), np.asarray(y, float)
    w = np.ones_like(x) if e is None else 1. / np.asarray(e, float) ** 2
    M = np.column_stack([np.ones_like(x), x])
    cov = np.linalg.inv(M.T @ (w[:, None] * M))
    p = cov @ (M.T @ (w * y))
    chi2 = float(np.sum(w * (y - M @ p) ** 2) / max(len(x) - 2, 1))
    s = float(np.sqrt(cov[1, 1] * (max(chi2, 1.) if e is not None else chi2)))
    return float(p[1]), s


# ---------------------------------------------------------------------------
# Calibration quality per night, from the reductions themselves
# ---------------------------------------------------------------------------

def arc_quality(night):
    """ The night's arc: its frame, UT, the arxiv cross-correlation, the fit RMS.

    ``cc`` is PypeIt's cross-correlation of each order's arc against its
    archived template (``autoid.reidentify``), read from ``run_pypeit.log``:
    a direct measure of how much the arc looks like a HIRES ThAr at all.
    """
    red = os.path.join(rr.REDUX_ROOT, 'reduce_' + night)
    pyp = glob.glob(os.path.join(red, '*.pypeit'))[0]
    arc = None
    with open(pyp) as fh:
        for line in fh:
            parts = [p.strip() for p in line.split('|')]
            if len(parts) > 2 and parts[0].startswith('HI.') and 'arc' in parts[1]:
                arc = parts[0].replace('.fits', '')
    cc = []
    with open(os.path.join(red, 'run_pypeit.log'), errors='replace') as fh:
        for line in fh:
            m = re.search(r'shift = [-0-9.]+, stretch = [0-9.]+, cc = ([0-9.]+)', line)
            if m:
                cc.append(float(m.group(1)))
    summ = {n['night']: n for n in json.load(open(os.path.join(rr.REDUX_ROOT, 'run_summary.json')))}
    s = summ[night]
    return dict(arc=arc, arc_day=arc.split('.')[1] if arc else '',
                arc_ut_h=int(arc.split('.')[2]) / 3600. if arc else np.nan,
                arc_cc_med=float(np.median(cc)) if cc else np.nan,
                arc_cc_frac_low=float(np.mean(np.asarray(cc) < 0.8)) if cc else np.nan,
                wave_rms_med=float(s['rms_med']), wave_rms_max=float(s['rms_max']),
                n_solutions=int(s['n_solutions']))


def hours_from_arc(vt, q):
    """ Science mid-time minus arc time, hours (borrowed arcs a day or more away). """
    out = []
    for r in vt:
        a = q[r['night']]
        if not a['arc']:
            out.append(np.nan)
            continue
        d_days = (np.datetime64('{}-{}-{}'.format(a['arc_day'][:4], a['arc_day'][4:6], a['arc_day'][6:]))
                  - np.datetime64('{}-{}-{}'.format(r['night'][:4], r['night'][4:6], r['night'][6:])))
        out.append(r['ut_h'] - (a['arc_ut_h'] + 24. * d_days.astype(int)))
    return np.array(out)


# ---------------------------------------------------------------------------
# The PypeIt wavelength solution, as the iodine sees it
# ---------------------------------------------------------------------------

def pypeit_wave_quality(ch):
    """ Per epoch and order: median and within-order robust scatter of zp. """
    rows = []
    g = ch[ch['good'] == 'True'] if ch['good'].dtype.kind == 'U' else ch[ch['good']]
    for k in np.unique(g['koaid']):
        s = g[g['koaid'] == k]
        for o in np.unique(s['order']):
            q = s[s['order'] == o]
            if len(q) < 10:
                continue
            rows.append(dict(koaid=k, order=int(o), zp=float(-np.median(q['dv_wave_c'])),
                             zp_within=float(robust_sigma(np.asarray(q['dv_wave_c']))),
                             chi2=float(np.median(q['redchi_meas'])), n=len(q)))
    return Table(rows)


# ---------------------------------------------------------------------------
# Re-combination variants (prompt 7's combination)
# ---------------------------------------------------------------------------

def recombine(vt, ch, drop_orders=(), drop_epochs=(), drop_chunks=None):
    koaids = [k for k in vt['koaid'] if k not in drop_epochs]
    ok = np.isfinite(ch['velocity']) & np.isfinite(ch['e_velocity'])
    if drop_chunks is not None:
        ok &= ~drop_chunks
    for o in drop_orders:
        ok &= ch['order'] != o
    nchunk = 700
    V = np.full((len(koaids), nchunk), np.nan)
    for i, k in enumerate(koaids):
        s = ch[ok & (ch['koaid'] == k)]
        V[i, np.asarray(s['tchunk'])] = np.asarray(s['velocity'], dtype=float)
    rv, aux = fa.combine(V, 14 - len(drop_orders))
    sel = np.isin(vt['koaid'], koaids)
    dbvc = np.asarray(vt['bvc'][sel] - vt['bvc_template'][sel])
    return np.asarray(koaids), rv['rv'] + dbvc, rv['rv_err'], sel


def variant_stats(name, t, v, e_int, sig_epoch, orb):
    r = v - orb['model'](t)
    r -= np.mean(r)
    f = f2.fit_circular(t, v)
    return dict(variant=name, n=len(t), resid_rms=float(np.std(r, ddof=1)),
                resid_robust=float(robust_sigma(r)),
                chi2_dof=float(np.sum((r / sig_epoch) ** 2) / (len(r) - 1)),
                K=f['K'], sigma_K=f['sigma_K'], nsig=f['K'] / f['sigma_K'],
                med_err_pyodine=float(np.median(e_int)))


# ---------------------------------------------------------------------------
# Tellurics
# ---------------------------------------------------------------------------

#: `pyodine`'s CARMENES telluric mask (Wallace et al. 2011 atlas), whose
#: first column is the VACUUM wavelength; interval edges are flagged 1.
TELL_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'vendor',
                         'pyodine', 'pyodine', 'tellurics', 'CARMENES',
                         'Wallace11_mask_0125_ext_carm_short_waveair.dat')


def telluric_intervals(path=TELL_FILE):
    t = np.genfromtxt(path, skip_header=1, unpack=True)
    i = np.where(t[2] == 1.)[0]
    return t[0, i[0::2]], t[0, i[1::2]]


def telluric_flags(ch, lo, hi):
    """ Does a telluric line fall inside each chunk, in each epoch?

    Chunks are comoving with the star, tellurics are fixed in the observed
    frame, and PypeIt's wavelengths are observed-frame vacuum -- the frame of
    the mask.  A chunk spans 40 pixels; its half-width is 20 pixels at the
    order's local dispersion (from neighbouring chunk centres).
    """
    flag = np.zeros(len(ch), dtype=bool)
    for k in np.unique(ch['koaid']):
        for o in np.unique(ch['order'][ch['koaid'] == k]):
            idx = np.where((ch['koaid'] == k) & (ch['order'] == o))[0]
            idx = idx[np.argsort(np.asarray(ch['pixc'][idx]))]
            w = np.asarray(ch['wave'][idx], float)
            px = np.asarray(ch['pixc'][idx], float)
            disp = np.gradient(w, px) if len(idx) > 1 else np.full(len(idx), 0.04)
            hw = 20. * np.abs(disp)
            for j, ii in enumerate(idx):
                flag[ii] = np.any((hi >= w[j] - hw[j]) & (lo <= w[j] + hw[j]))
    return flag


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(vt, bstar, wq, outfile):
    import matplotlib.pyplot as plt
    f2.set_style()
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.4))

    # (a) nightly: PypeIt zero-point error against the residuals
    ax = axs[0, 0]
    ax.grid(zorder=0)
    borrowed = np.isin(vt['night'], list(BORROWED))
    for sel, mk, lab in ((~borrowed, 'o', 'own arc'), (borrowed, 'D', 'borrowed arc')):
        ax.plot(vt['zp'][sel], vt['resid_p2'][sel], mk, ms=6, color=C_PHASE2, mfc='none',
                label='phase 2 residual, ' + lab)
        ax.errorbar(vt['zp'][sel], vt['resid'][sel], yerr=vt['sig_epoch'][sel], fmt=mk, ms=6,
                    color=C_OURS, ecolor=C_OURS, elinewidth=1.1, capsize=0, mfc='white', mew=1.6,
                    label='phase 3 residual, ' + lab)
    ax.axhline(0., color=INK_2, lw=0.8)
    ax.set_xlabel('PypeIt - iodine wavelength zero point (m/s)')
    ax.set_ylabel('residual about the Keplerian (m/s)')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_title('(a) the ThAr zero point drives phase 2, not phase 3', loc='left', fontsize=10)

    # (b) intra-night: 08-25 and 08-26
    ax = axs[0, 1]
    ax.grid(zorder=0)
    for night, mk in (('19980825', 'o'), ('19980826', 's')):
        s = vt['night'] == night
        lab = '{}-{}'.format(night[4:6], night[6:])
        ax.plot(vt['ut_h'][s], vt['zp'][s] - np.mean(vt['zp'][s]), mk + '--', color=C_CAT,
                ms=6, label='PypeIt zero-point error, HD 187123, ' + lab)
        ax.plot(vt['ut_h'][s], vt['resid_p2'][s] - np.mean(vt['resid_p2'][s]), mk + ':',
                color=C_PHASE2, ms=6, mfc='none', label='phase 2 residual, ' + lab)
        ax.errorbar(vt['ut_h'][s], vt['resid'][s] - np.mean(vt['resid'][s]),
                    yerr=vt['sig_epoch'][s], fmt=mk + '-', color=C_OURS, ms=6, mfc='white',
                    mew=1.6, label='phase 3 residual, ' + lab)
    ax.plot(bstar['ut_h'], bstar['zp'] - np.mean(bstar['zp']), '^', color=C_THIRD, ms=6,
            label='PypeIt zero-point error, B stars, 08-26')
    ax.axhline(0., color=INK_2, lw=0.8)
    ax.set_xlabel('UT (hours)')
    ax.set_ylabel('m/s, each series about its nightly mean')
    ax.set_ylim(-1300, 1300)
    ax.legend(fontsize=7, loc='lower left', ncol=1)
    ax.set_title('(b) the drift within a night is PypeIt\'s, and phase 3 removes it',
                 loc='left', fontsize=10)

    # (c) the PypeIt solution within each order, as the iodine sees it
    ax = axs[1, 0]
    ax.grid(zorder=0)
    orders = np.unique(wq['order'])
    typ = np.array([np.median(wq['zp_within'][wq['order'] == o]) for o in orders])
    lo = np.array([np.percentile(wq['zp_within'][wq['order'] == o], 10) for o in orders])
    hi = np.array([np.percentile(wq['zp_within'][wq['order'] == o], 90) for o in orders])
    ax.fill_between(orders, lo, hi, color=C_PHASE2, alpha=0.25, label='all epochs, 10-90%')
    ax.plot(orders, typ, '-', color=C_PHASE2, label='all epochs, median')
    for k, col, lab in (('HI.19980913.29698', C_CAT, '1998-09-13'),
                        ('HI.19971223.17175', C_THIRD, '1997-12-23 (no pixel flat)')):
        s = wq[wq['koaid'] == k]
        ax.plot(s['order'], s['zp_within'], 'o-', color=col, label=lab)
    ax.set_yscale('log')
    ax.set_xlabel('echelle order')
    ax.set_ylabel('PypeIt zero point, chunk-to-chunk scatter in an order (m/s)')
    ax.legend(fontsize=8)
    ax.set_title('(c) PypeIt\'s wavelength solution within an order, measured by the iodine',
                 loc='left', fontsize=10)

    # (d) the seasonal term: residual against BVC, and within nights
    ax = axs[1, 1]
    ax.grid(zorder=0)
    ax.errorbar(vt['bvc'] / 1e3, vt['resid'], yerr=vt['sig_epoch'], fmt='o', color=C_OURS,
                ms=5, mfc='white', mew=1.5, ecolor=C_OURS, elinewidth=1.1, capsize=0,
                label='all epochs')
    ax.axhline(0., color=INK_2, lw=0.8)
    ax.set_xlabel('barycentric correction (km/s)')
    ax.set_ylabel('residual about the Keplerian (m/s)')
    ax.legend(fontsize=8)
    ax.set_title('(d) the seasonal correlation prompt 8 flagged', loc='left', fontsize=10)

    fig.tight_layout()
    fig.savefig(outfile, dpi=130)
    plt.close(fig)
    print('  wrote {:s}'.format(outfile))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parser(options=None):
    p = argparse.ArgumentParser(description='Diagnose what is left (phase 3, prompt 9)')
    p.add_argument('--datadir', type=str, default=DATA_DIR)
    p.add_argument('--figdir', type=str, default=FIG_DIR)
    return p.parse_args() if options is None else p.parse_args(options)


def main(pargs):
    vt, ch, orb = load()
    nights = sorted(set(vt['night']) | set(NO_PIXFLAT))
    q = {n: arc_quality(n) for n in nights}
    vt['h_from_arc'] = hours_from_arc(vt, q)
    wq = pypeit_wave_quality(ch)
    e = np.asarray(vt['sig_epoch'])

    # ---- per-night table ---------------------------------------------------
    rows = []
    for n in nights:
        s = vt['night'] == n
        a = q[n]
        row = dict(night=n, n_epochs=int(s.sum()),
                   arc_source=BORROWED.get(n, 'own'),
                   flat='iodine-in trace, no pixel flat' if n in NO_PIXFLAT else 'donor 1998-08-12',
                   arc=a['arc'], arc_cc_med=a['arc_cc_med'], arc_cc_frac_low=a['arc_cc_frac_low'],
                   wave_rms_med=a['wave_rms_med'], wave_rms_max=a['wave_rms_max'],
                   n_solutions=a['n_solutions'])
        if s.any():
            w = 1. / e[s] ** 2
            row.update(zp_mean=float(np.mean(vt['zp'][s])), zp_span=float(np.ptp(vt['zp'][s])),
                       resid_p2_mean=float(np.mean(vt['resid_p2'][s])),
                       resid_p2_span=float(np.ptp(vt['resid_p2'][s])),
                       resid_mean=float(np.sum(w * vt['resid'][s]) / np.sum(w)),
                       resid_err=float(1. / np.sqrt(np.sum(w))),
                       resid_span=float(np.ptp(vt['resid'][s])),
                       h_from_arc=float(np.mean(vt['h_from_arc'][s])))
        rows.append(row)
    nt = Table(rows)
    for c in nt.colnames:
        if nt[c].dtype.kind == 'f':
            nt[c].info.format = '.3f' if c.startswith(('arc_cc', 'wave_rms')) else '.1f'
    nt.write(os.path.join(pargs.datadir, 'diagnose_nights.csv'), overwrite=True)
    print('\n=== Nights ===')
    nt.pprint(max_width=250, max_lines=40)

    # ---- 1. borrowed arcs -------------------------------------------------
    print('\n=== 1. The borrowed-arc nights ===')
    own = ~np.isin(vt['night'], list(BORROWED))
    print('  PypeIt zero-point error |zp|: own-arc epochs median {:.0f} m/s (rms {:.0f}); '
          'borrowed: {}'.format(np.median(np.abs(vt['zp'][own])), np.std(vt['zp'][own]),
                                ', '.join('{:+.0f}'.format(z) for z in vt['zp'][~own])))
    print('  phase 3 residual: own-arc rms {:.1f} m/s; borrowed-arc epochs {}'.format(
        np.std(vt['resid'][own], ddof=1), ', '.join('{:+.0f} ({:+.1f} sigma)'.format(r, r / s)
                                                    for r, s in zip(vt['resid'][~own], e[~own]))))
    print('  phase 2 residual: own-arc rms {:.0f} m/s; borrowed-arc epochs {}'.format(
        np.std(vt['resid_p2'][own], ddof=1), ', '.join('{:+.0f}'.format(r) for r in vt['resid_p2'][~own])))
    rho, p = stats.spearmanr(vt['zp'], vt['resid'])
    rho2, p2 = stats.spearmanr(vt['zp'], vt['resid_p2'])
    sl, sle = fit_line(vt['zp'], vt['resid'], e)
    print('  residual against PypeIt zero point: phase 2 Spearman {:+.2f} (p {:.0e}); phase 3 '
          '{:+.2f} (p {:.2f}), slope {:+.4f} +- {:.4f}'.format(rho2, p2, rho, p, sl, sle))
    ok = np.isfinite(vt['h_from_arc'])
    rz, pz = stats.spearmanr(np.abs(vt['h_from_arc'][ok]), np.abs(vt['zp'][ok]))
    print('  |PypeIt zero-point error| against hours from the arc: Spearman {:+.2f} (p {:.3f})'.format(rz, pz))

    # ---- 2. intra-night ---------------------------------------------------
    print('\n=== 2. Drifts within a night ===')
    bt = Table.read(os.path.join(DATA_DIR, 'bstar_chunks.csv'))
    bt = bt[(bt['run'] == 1) & (bt['ok'] == 'True')] if bt['ok'].dtype.kind == 'U' \
        else bt[(bt['run'] == 1) & bt['ok']]
    brow = []
    for fr in np.unique(bt['frame']):
        s = bt[bt['frame'] == fr]
        brow.append(dict(frame=fr, star=s['star'][0], ut_h=int(fr.split('.')[2]) / 3600.,
                         zp=float(-np.median(s['dv_wave_c'])),
                         zp_err=float(1.2533 * robust_sigma(np.asarray(s['dv_wave_c'])) / np.sqrt(len(s)))))
    bstar = Table(brow)
    for n in [m for m in nights if np.sum(vt['night'] == m) > 1]:
        s = vt['night'] == n
        h = np.asarray(vt['ut_h'][s])
        out = []
        for col, err in (('resid_p2', None), ('zp', None), ('resid', e[s])):
            y = np.asarray(vt[col][s], float)
            if s.sum() > 2:
                b, be = fit_line(h, y, err)
                out.append('{:s} {:+6.0f} +- {:4.0f} m/s/h (span {:4.0f})'.format(col, b, be, np.ptp(y)))
            else:
                out.append('{:s} span {:4.0f}'.format(col, np.ptp(y)))
        print('  {:s} ({:d} frames, {:.1f} h): {:s}'.format(n, int(s.sum()), np.ptp(h), '; '.join(out)))
    b, be = fit_line(bstar['ut_h'], bstar['zp'], bstar['zp_err'])
    pair = [abs(bstar['zp'][i + 1] - bstar['zp'][i]) for i in range(0, len(bstar) - 1, 2)]
    print('  08-26 B stars ({:d} frames, {:.1f} h): PypeIt zero point {:s}; slope {:+.0f} +- {:.0f} '
          'm/s/h; back-to-back pairs differ by {:s} m/s'.format(
              len(bstar), np.ptp(bstar['ut_h']), ', '.join('{:+.0f}'.format(z) for z in bstar['zp']),
              b, be, ', '.join('{:.0f}'.format(p_) for p_ in pair)))
    s = vt['night'] == '19980826'
    print('  08-26 HD 187123 zero points {:s} at UT {:s}'.format(
        ', '.join('{:+.0f}'.format(z) for z in vt['zp'][s]),
        ', '.join('{:.2f}'.format(h) for h in vt['ut_h'][s])))
    wn = [(vt['resid'][vt['night'] == n] - np.average(vt['resid'][vt['night'] == n],
                                                        weights=1. / e[vt['night'] == n] ** 2)) /
          e[vt['night'] == n] for n in nights if np.sum(vt['night'] == n) > 1]
    wn = np.concatenate(wn)
    nmulti = sum(np.sum(vt['night'] == n) > 1 for n in nights)
    print('  phase 3 residuals about nightly means, all multi-frame nights: chi2/dof {:.2f} '
          '({:d} frames, {:d} nights)'.format(np.sum(wn ** 2) / (len(wn) - nmulti), len(wn), nmulti))

    # ---- 3. December ------------------------------------------------------
    print('\n=== 3. The December nights (no pixel flat) ===')
    print('  1997-12-24: HI.19971224.16259 is cell-out -- no iodine velocity can be made from it.')
    k = 'HI.19971223.17175'
    r = vt[vt['koaid'] == k][0]
    full = vt[(vt['n_orders'] == 14) & (vt['koaid'] != k)]
    for c in ('redchi_meas', 'c2c_scatter', 'order_scatter', 'sig_epoch', 'n_good'):
        pct = stats.percentileofscore(np.asarray(full[c], float), float(r[c]))
        print('  {:<14s} 12-23 {:8.1f}; other 14-order epochs median {:8.1f} (12-23 at percentile {:.0f})'.format(
            c, r[c], np.median(full[c]), pct))
    print('  residual {:+.1f} +- {:.1f} m/s ({:+.1f} sigma); PypeIt zero point {:+.0f} m/s; '
          '180 days before every other epoch'.format(r['resid'], r['sig_epoch'],
                                                     r['resid'] / r['sig_epoch'], r['zp']))
    d = wq[wq['koaid'] == k]
    others = wq[np.isin(wq['koaid'], list(full['koaid']))]
    ratio = [(o, float(d['chi2'][d['order'] == o][0] / np.median(others['chi2'][others['order'] == o])))
             for o in d['order']]
    print('  chunk chi2 by order, 12-23 / other epochs: {:s}'.format(
        ', '.join('{:d}:{:.2f}'.format(o, x) for o, x in ratio)))
    wz = [(o, float(d['zp_within'][d['order'] == o][0] / np.median(others['zp_within'][others['order'] == o])))
          for o in d['order']]
    print('  PypeIt within-order zero-point scatter, 12-23 / others: median {:.2f} (range {:.2f}-{:.2f})'.format(
        np.median([x for _, x in wz]), min(x for _, x in wz), max(x for _, x in wz)))

    # ---- 4. 1998-09-13 ----------------------------------------------------
    print('\n=== 4. 1998-09-13 ===')
    a = q['19980913']
    print('  arc {:s}: arxiv cross-correlation median cc {:.2f} ({:.0%} of orders < 0.8); '
          'other nights median cc {:.2f}'.format(a['arc'], a['arc_cc_med'], a['arc_cc_frac_low'],
                                                 np.median([q[n]['arc_cc_med'] for n in nights if n != '19980913'])))
    print('  {:d} wavelength solutions, worst RMS {:.3f} px; orders 92, 86 and 71 masked '
          'BADWVCALIB+BADTILTCALIB+BADFLATCALIB -> 34 extracted'.format(a['n_solutions'], a['wave_rms_max']))
    k = 'HI.19980913.29698'
    r = vt[vt['koaid'] == k][0]
    for c in ('redchi_meas', 'c2c_scatter', 'order_scatter', 'sig_epoch', 'n_good'):
        pct = stats.percentileofscore(np.asarray(full[c], float), float(r[c]))
        print('  {:<14s} 09-13 {:8.1f}; 14-order epochs median {:8.1f} (percentile {:.0f})'.format(
            c, r[c], np.median(full[c]), pct))
    print('  residual {:+.1f} +- {:.1f} m/s ({:+.1f} sigma); PypeIt zero point {:+.0f} m/s'.format(
        r['resid'], r['sig_epoch'], r['resid'] / r['sig_epoch'], r['zp']))
    d = wq[wq['koaid'] == k]
    wz = [(int(o), float(d['zp_within'][d['order'] == o][0] / np.median(others['zp_within'][others['order'] == o])))
          for o in d['order']]
    print('  PypeIt within-order zero-point scatter, 09-13 / others: {:s}'.format(
        ', '.join('{:d}:{:.1f}'.format(o, x) for o, x in wz)))
    ratio = [(int(o), float(d['chi2'][d['order'] == o][0] / np.median(others['chi2'][others['order'] == o])))
             for o in d['order']]
    print('  chunk chi2 by order, 09-13 / others: {:s}'.format(
        ', '.join('{:d}:{:.2f}'.format(o, x) for o, x in ratio)))
    po = Table.read(os.path.join(DATA_DIR, 'iodine_velocities_per_order.csv'))
    s = po[po['koaid'] == k]
    dev = np.asarray(s['v_bary']) - r['v_bary']
    typ = {int(o): robust_sigma(np.asarray(po['v_bary'][(po['order'] == o)]) -
                                np.asarray([vt['v_bary'][list(vt['koaid']).index(kk)]
                                            for kk in po['koaid'][po['order'] == o]]))
           for o in s['order']}
    print('  order velocity - epoch velocity, in units of that order\'s usual scatter: {:s}'.format(
        ', '.join('{:d}:{:+.1f}'.format(int(o), x / typ[int(o)]) for o, x in zip(s['order'], dev))))

    # ---- 5. what is left --------------------------------------------------
    print('\n=== 5. What is left ===')
    t = np.asarray(vt['jd'], float)
    res = np.asarray(vt['resid'], float)
    bvc = np.asarray(vt['bvc'], float) / 1e3
    sl, sle = fit_line(bvc, res, e)
    st, ste = fit_line(t, res, e)
    print('  residual vs BVC: {:+.2f} +- {:.2f} m/s per km/s; vs time: {:+.3f} +- {:.3f} m/s per day'.format(
        sl, sle, st, ste))
    print('  BVC and time over the 1998 season (excluding December): Pearson r = {:+.2f}'.format(
        np.corrcoef(bvc[t > 2450900], t[t > 2450900])[0, 1]))
    # Within nights: BVC changes by Earth's rotation, the season does not
    dres, dbvc = [], []
    for n in nights:
        s = np.where(vt['night'] == n)[0]
        for i in range(len(s)):
            for j in range(i + 1, len(s)):
                dres.append(res[s[j]] - res[s[i]])
                dbvc.append(bvc[s[j]] - bvc[s[i]])
    dres, dbvc = np.array(dres), np.array(dbvc)
    if len(dres) > 2:
        b = float(np.sum(dres * dbvc) / np.sum(dbvc ** 2))
        be = float(np.std(dres - b * dbvc, ddof=1) / np.sqrt(np.sum(dbvc ** 2)))
        print('  within nights ({:d} pairs, |dBVC| up to {:.2f} km/s): {:+.1f} +- {:.1f} m/s per km/s '
              '(seasonal fit predicts {:+.1f})'.format(len(dres), np.max(np.abs(dbvc)), b, be, sl))
    # Leave out the 1997 epoch: is the "season" just one point?
    s = t > 2450900
    sl2, sle2 = fit_line(bvc[s], res[s], e[s])
    print('  without 1997-12-23: {:+.2f} +- {:.2f} m/s per km/s'.format(sl2, sle2))

    # Tellurics: an external, Keplerian-free criterion
    lo, hi = telluric_intervals()
    tell = telluric_flags(ch, lo, hi)
    ch['telluric'] = tell
    print('  telluric lines (CARMENES mask) inside a chunk, fraction of chunk-epochs by order: {:s}'.format(
        ', '.join('{:d}:{:.0%}'.format(int(o), np.mean(tell[ch['order'] == o])) for o in np.unique(ch['order']))))
    frac = np.zeros(700)
    sig = np.full(700, np.nan)
    order_of = np.zeros(700, dtype=int)
    for tc in np.unique(ch['tchunk']):
        m = ch['tchunk'] == tc
        frac[tc] = np.mean(tell[m])
        sig[tc] = np.nanmedian(ch['sigma_chunk'][m])
        order_of[tc] = int(ch['order'][m][0])
    fin = np.isfinite(sig) & (sig < 999.)
    for lab, osel in (('all orders', np.ones(700, bool)), ('orders 60-62', np.isin(order_of, (60, 61, 62)))):
        a = fin & osel & (frac > 0)
        c = fin & osel & (frac == 0)
        u, pu = stats.mannwhitneyu(sig[a], sig[c]) if a.sum() > 3 and c.sum() > 3 else (np.nan, np.nan)
        print('  chunk time-series sigma, {:s}: telluric-crossed {:.0f} m/s ({:d} chunks) vs clean {:.0f} m/s '
              '({:d}); Mann-Whitney p {:.1e}'.format(lab, np.median(sig[a]), int(a.sum()), np.median(sig[c]),
                                                    int(c.sum()), pu))
    # Does the crossing move with the season?  Chunks crossed in some epochs
    # and not others are the ones whose velocity a telluric could drag
    part = fin & (frac > 0) & (frac < 1)
    print('  chunks crossed in some epochs only: {:d}; their sigma {:.0f} m/s'.format(
        int(part.sum()), np.median(sig[part]) if part.any() else np.nan))

    # Re-combination variants
    rows = []
    kk, v, eint, sel = recombine(vt, ch)
    print('  re-combination check: max |v - prompt 7| = {:.2e} m/s'.format(
        np.max(np.abs(v - np.asarray(vt['v_bary'])))))
    rows.append(variant_stats('prompt 7 (all)', t, v, eint, e, orb))
    few = list(vt['koaid'][vt['n_orders'] <= 8])
    kk, v, eint, sel = recombine(vt, ch, drop_epochs=few)
    rows.append(variant_stats('without the {:d} epochs of <= 8 orders'.format(len(few)), t[sel], v, eint,
                              e[sel], orb))
    kk, v, eint, sel = recombine(vt, ch, drop_orders=RED_ORDERS)
    rows.append(variant_stats('without orders 58-61', t, v, eint, e, orb))
    kk, v, eint, sel = recombine(vt, ch, drop_orders=RED_ORDERS, drop_epochs=few)
    rows.append(variant_stats('both', t[sel], v, eint, e[sel], orb))
    s = ~np.isin(vt['night'], list(BORROWED) + ['19980913', '19971223'])
    kk, v, eint, sel = recombine(vt, ch, drop_epochs=list(vt['koaid'][~s]))
    rows.append(variant_stats('without the borrowed-arc, December and 09-13 epochs', t[sel], v, eint,
                              e[sel], orb))
    kk, v, eint, sel = recombine(vt, ch, drop_chunks=tell)
    rows.append(variant_stats('telluric-crossed chunks masked', t, v, eint, e, orb))
    kk, v, eint, sel = recombine(vt, ch, drop_chunks=tell, drop_orders=(58, 59))
    rows.append(variant_stats('tellurics masked, orders 58-59 dropped', t, v, eint, e, orb))
    kk, v, eint, sel = recombine(vt, ch, drop_chunks=tell, drop_orders=(58, 59), drop_epochs=few)
    rows.append(variant_stats('tellurics masked, 58-59 and the <= 8-order epochs dropped', t[sel], v, eint,
                              e[sel], orb))
    var = Table(rows)
    for c in var.colnames:
        if var[c].dtype.kind == 'f':
            var[c].info.format = '.2f'
    var.pprint(max_width=200)
    var.write(os.path.join(pargs.datadir, 'diagnose_variants.csv'), overwrite=True)

    make_figure(vt, bstar, wq, os.path.join(pargs.figdir, 'fig_p9_diagnosis.png'))


if __name__ == '__main__':
    main(parser())
