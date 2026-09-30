""" Fit one epoch end to end (phase 3, prompt 5).

`pyodine`'s own driver, `pyodine_model_observations.model_single_observation`,
run on one HIRES epoch through ``utilities_hires`` -- run 0 (single-Gaussian
LSF) and run 1 (whatever ``pyodine_parameters`` sets: Lick's multi-Gaussian
when prompt 5 ran, the super-Gaussian since prompt 6) over every chunk of the
template -- against both prompt-3 templates, as the Q&A after prompt 3
decided.  Then an
analysis of what it produced:

  * the per-chunk velocity scatter, robust and not, beside the photon-noise
    error each chunk would have if the model fitted;
  * how many chunks fail, in every sense the fit can fail;
  * the LSF against position on the detector, beside the same night's ThAr;
  * the fitted wavelength solution against PypeIt's;
  * the fitted continuum against PypeIt's extracted blaze (``OPT_FLAT``);
  * the runtime.

WHICH EPOCH.  ``HI.19980825.19425``: cell in, 430 s, the night before the
template, and the highest-S/N cell-in frame whose 14 iodine orders all pass
phase 2's quality filter.  That last condition is not taste.  The driver fits
chunk *i* of the observation against chunk *i* of the template, but
`auto_wave_comoving` builds observation chunks only for the orders it is
given -- so passing the driver a subset of orders silently pairs every chunk
after the gap with the wrong template chunk.  Prompt 4's epoch, 19326, has
two rejected orders and cannot be run without that trap; see the Report.  Its
raw OBJECT reads 'Gl 83.1', a stale keyword: the pointing (19:46:58
+34:25:10) and TARGNAME are HD 187123's.

FAILURES ARE READ FROM THE INFO LOG.  `model_single_observation` catches every
exception, logs it with `logging.info` -- not `logging.error` -- and returns
normally.  Its error log stays empty when it has failed.  This driver checks
that every result file exists and that the info log has no traceback.

THE NOISE.  As decided after prompt 3: weights are PypeIt's inverse
variance, every chi-square is also reported divided by ``NOISE_SCALE**2``,
and a chunk's photon-noise velocity error is its `lmfit` error with the
`scale_covar` rescaling undone and the measured noise put in:
``sigma * noise_scale / sqrt(redchi)`` (validated by Monte Carlo on one chunk
in prompt 4).

Run with:

    conda run --no-capture-output -n pypeit14 \\
        python -m first_hires_exoplanet.fit_one_epoch

"""

# Must precede any import of matplotlib, including pyodine's own
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# Standard imports
import time
import argparse

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

from first_hires_exoplanet import utilities_hires
from first_hires_exoplanet.utilities_hires import load_pyodine as lp
from first_hires_exoplanet.utilities_hires import pyodine_parameters as hires_pars
from first_hires_exoplanet import deconvolve_template as dt
from first_hires_exoplanet import fit_one_chunk as f1
from first_hires_exoplanet.pyodine_chunk_stats import summarise, robust_sigma

import pyodine_model_observations as pmo
from pyodine.fitters import results_io
from pyodine.models import lsf as plsf
from pyodine.models.base import ParameterSet
from pyodine.template.base import StellarTemplate_Chunked


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_HERE, 'data')
FIG_DIR = os.path.join(os.path.dirname(_HERE), 'docs', 'figs')
RUN_ROOT = os.path.join(lp.DATA_ROOT, 'pyodine_runs')

EPOCH = 'HI.19980825.19425'
CKMS = 299792.458

#: Chunks whose chi-square exceeds this multiple of their order's median are
#: counted as misfits
MISFIT_FACTOR = 3.
#: Chunk velocities further than this many robust sigma from the median are
#: counted as outliers
OUTLIER_NSIG = 5.


# ---------------------------------------------------------------------------
# Running the driver
# ---------------------------------------------------------------------------

def run_driver(obs_file, temp_file, outdir, Pars, orders=None):
    """ `model_single_observation`, with its silent failures made loud.

    ``orders``, if given, are the observation's order positions to fit (the
    driver's own argument); safe for a subset since fork change 5.

    Returns:
        dict: result file per run, runtimes, and the info-log counts
    """
    os.makedirs(outdir, exist_ok=True)
    res = {r: os.path.join(outdir, 'run{:d}.h5'.format(r)) for r in Pars.model_runs}
    info = os.path.join(outdir, 'info.log')
    err = os.path.join(outdir, 'errors.log')
    for p in list(res.values()) + [info, err]:
        if os.path.exists(p):
            os.remove(p)
    t0 = time.time()
    pmo.model_single_observation(utilities_hires, Pars, obs_file, temp_file,
                                 orders=orders, plot_dir=outdir, res_names=[res[r] for r in sorted(res)],
                                 error_log=err, info_log=info, quiet=True)
    total = time.time() - t0
    with open(info) as fh:
        log = fh.read()
    if 'Traceback' in log or 'Something went wrong' in log:
        tail = [l for l in log.splitlines() if 'Error' in l or 'error' in l]
        raise RuntimeError('pyodine failed (see {:s}): {:s}'.format(
            info, tail[-1] if tail else 'traceback in info log'))
    missing = [p for p in res.values() if not os.path.exists(p)]
    if missing:
        raise RuntimeError('pyodine wrote no result for {}'.format(missing))
    # Per-run wall time from the result files (each includes its plots)
    t_run, prev = {}, t0
    for r in sorted(res):
        t_run[r] = os.path.getmtime(res[r]) - prev
        prev = os.path.getmtime(res[r])
    counts = {}
    for key, label in (('no_uncertainties', 'Number of chunks with no uncertainties: '),
                       ('chauvenet', 'Number of chunks with outliers: '),
                       ('nan_redchi', 'Number of chunks with nan fitted red. Chi2: ')):
        counts[key] = [int(l.split(label)[1]) for l in log.splitlines() if label in l]
    return dict(res=res, total=total, t_run=t_run, counts=counts, log=info)


# ---------------------------------------------------------------------------
# Per-chunk measurements
# ---------------------------------------------------------------------------

def blaze(obs_file, obs):
    """ PypeIt's extracted flat (``OPT_FLAT``) per order position, NaN where masked. """
    out = np.full_like(obs._flux, np.nan)
    with fits.open(obs_file) as hdul:
        for hdu in hdul[1:]:
            if not hasattr(hdu, 'columns') or 'OPT_FLAT' not in hdu.columns.names:
                continue
            i = obs.index_of(int(hdu.header['ECH_ORDER']))
            f = np.asarray(hdu.data['OPT_FLAT'], dtype=float)
            out[i] = np.where(f > 0, f, np.nan)
    return out


def lsf_shape(x, y):
    """ FWHM, centroid and half-maximum asymmetry of a sampled LSF, pixels.

    Asymmetry is (right half-width - left half-width) / FWHM, measured from
    the peak: positive means a red wing.
    """
    k = int(np.argmax(y))
    half = 0.5 * y[k]
    lo = k
    while lo > 0 and y[lo - 1] >= half:
        lo -= 1
    hi = k
    while hi < len(y) - 1 and y[hi + 1] >= half:
        hi += 1
    if lo == 0 or hi == len(y) - 1:
        return np.nan, np.nan, np.nan
    xl = np.interp(half, [y[lo - 1], y[lo]], [x[lo - 1], x[lo]])
    xr = np.interp(half, [y[hi + 1], y[hi]], [x[hi + 1], x[hi]])
    fwhm = xr - xl
    centroid = np.sum(x * y) / np.sum(y)
    return float(fwhm), float(centroid), float(((xr - x[k]) - (x[k] - xl)) / fwhm)


def chunk_table(res_file, run_id, Pars, obs, template, iod, flat, weight, k,
                tmpl_index=None):
    """ One row per chunk: everything fitted, and everything it is compared with.

    ``tmpl_index[i]`` is result chunk *i*'s template chunk; needed whenever the
    epoch was fitted over a subset of the template's orders (else chunk *i*).
    """
    r = results_io.load_results(res_file)
    ch, par, err = r['chunks'], r['params'], r['errors']
    nchunk = len(ch['order'])
    reports = [s.decode() if isinstance(s, bytes) else str(s) for s in r['reports']]
    lsf_model = Pars.model_runs[run_id]['lsf_model']
    if 'lsf_setup_dict' in Pars.model_runs[run_id]:
        lsf_model.adapt_LSF(Pars.model_runs[run_id]['lsf_setup_dict'])
    x_lsf = lsf_model.generate_x(20, Pars.lsf_conv_width)
    lsf_names = [n for n in par if n.startswith('lsf_')]

    rows = []
    for i in range(nchunk):
        pos, p0, p1 = int(ch['order'][i]), int(ch['firstpix'][i]), int(ch['lastpix'][i])
        n = p1 - p0 + 1
        pix = np.arange(-(n // 2), n - n // 2)
        pc = p0 + n // 2                                    # pixel at pix = 0
        row = dict(run=run_id, chunk=i, order=int(obs.ech_orders[pos]), pos=pos,
                   pix0=p0, pix1=p1, pixc=pc, wave=float(obs._wave[pos, pc]),
                   redchi=float(r['redchi2'][i]),
                   masked=int(np.sum(weight[pos, p0:p1 + 1] <= 0)),
                   no_unc=('uncertainties could not be estimated' in reports[i]),
                   failed=(reports[i] == 'Chunk failed...'))
        for name in par:
            row[name] = float(par[name][i])
            row['e_' + name] = float(err[name][i]) if err[name][i] is not None else np.nan
        row['redchi_meas'] = row['redchi'] / k ** 2
        row['sigv_noise'] = row['e_velocity'] * k / np.sqrt(row['redchi']) \
            if row['redchi'] > 0 else np.nan

        # Wavelength: the fitted straight line against PypeIt's solution
        w_fit = row['wave_intercept'] + row['wave_slope'] * pix
        w_pyp = obs._wave[pos, p0:p1 + 1]
        row['dv_wave_c'] = float((row['wave_intercept'] / obs._wave[pos, pc] - 1.) * CKMS * 1e3)
        row['dv_wave_mean'] = float(np.mean(w_fit / w_pyp - 1.) * CKMS * 1e3)
        row['slope_ratio'] = float(row['wave_slope'] / np.polyfit(pix, w_pyp, 1)[0])

        # Continuum: undo the model's normalisations by the chunk means, so
        # the fitted continuum is the flux the model puts at a line-free pixel
        half = n / 2. + Pars.chunk_padding
        wlo = row['wave_intercept'] - row['wave_slope'] * half
        whi = row['wave_intercept'] + row['wave_slope'] * half
        sel = (iod.wave >= wlo) & (iod.wave <= whi)
        t_mean = float(np.mean(iod.flux[sel])) if sel.any() else np.nan
        ti = i if tmpl_index is None else int(tmpl_index[i])
        s_mean = float(np.mean(template[ti].flux))
        # How much structure each component puts in this chunk: the RMS of
        # the (alpha-scaled) atlas and of the template about their means
        row['iod_contrast'] = float(np.std(iod.flux[sel]) / t_mean) if sel.any() else np.nan
        row['tem_contrast'] = float(np.std(template[ti].flux) / s_mean)
        c_true = row['cont_intercept'] * t_mean ** (-row['iod_depth']) * \
            (row['tem_depth'] * (1. / s_mean - 1.) + 1.)
        row['cont_true'] = float(c_true)
        row['cont_env'] = float(obs._cont[pos, pc])
        row['flat'] = float(flat[pos, pc])

        # LSF: the fitted profile, reduced to three numbers
        lp_ = ParameterSet({nm[4:]: row[nm] for nm in lsf_names})
        try:
            y = lsf_model.eval(x_lsf, lp_)
            row['lsf_w'], row['lsf_cen'], row['lsf_asym'] = lsf_shape(x_lsf, y)
        except Exception:                                    # noqa: BLE001
            row['lsf_w'] = row['lsf_cen'] = row['lsf_asym'] = np.nan
        rows.append(row)
    return Table(rows)


def flag_failures(t, Pars):
    """ Every way a chunk can fail, as boolean columns. """
    t['bad_velocity'] = ~np.isfinite(t['velocity']) | ~np.isfinite(t['e_velocity'])
    at_bound = np.zeros(len(t), dtype=bool)
    if 'lsf_fwhm' in t.colnames:
        at_bound |= (np.abs(t['lsf_fwhm'] - 0.5) < 1e-6) | (np.abs(t['lsf_fwhm'] - 4.0) < 1e-6)
    at_bound |= np.abs(t['iod_depth'] - 0.1) < 1e-6
    t['at_bound'] = at_bound
    med = {o: np.nanmedian(t['redchi'][t['order'] == o]) for o in np.unique(t['order'])}
    t['misfit'] = np.array([r > MISFIT_FACTOR * med[o] for r, o in zip(t['redchi'], t['order'])])
    ok = ~t['bad_velocity']
    v = np.asarray(t['velocity'], dtype=float)
    rs = robust_sigma(v[ok])
    t['v_outlier'] = ok & (np.abs(v - np.nanmedian(v[ok])) > OUTLIER_NSIG * rs)
    t['good'] = ok & ~t['at_bound'] & ~t['misfit'] & ~t['v_outlier'] & ~t['failed']
    return t


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def misfit_drivers(t):
    """ Which component's structure does the misfit follow?

    Spearman rank correlation of chi-square with the iodine contrast and with
    the stellar contrast of each good chunk, and the partial correlation of
    each with the other held fixed (on ranks).
    """
    from scipy.stats import spearmanr, rankdata
    g = t[t['good'] & np.isfinite(t['iod_contrast']) & np.isfinite(t['tem_contrast'])]
    if len(g) < 10:
        return dict(n=len(g), rho_iod=np.nan, rho_tem=np.nan, rho_iod_tem=np.nan,
                    partial_iod=np.nan, partial_tem=np.nan)
    c = np.log(np.asarray(g['redchi_meas'], dtype=float))
    ri, rs, rc = (rankdata(np.asarray(g[k], dtype=float)) for k in
                  ('iod_contrast', 'tem_contrast', 'redchi_meas'))

    def partial(y, x, z):
        # correlation of y and x after regressing both on z (ranks)
        ry = y - np.polyval(np.polyfit(z, y, 1), z)
        rx = x - np.polyval(np.polyfit(z, x, 1), z)
        return float(np.corrcoef(ry, rx)[0, 1])
    return dict(n=len(g),
                rho_iod=float(spearmanr(g['iod_contrast'], c)[0]),
                rho_tem=float(spearmanr(g['tem_contrast'], c)[0]),
                rho_iod_tem=float(spearmanr(g['iod_contrast'], g['tem_contrast'])[0]),
                partial_iod=partial(rc, ri, rs), partial_tem=partial(rc, rs, ri))


def summarise_run(t, expected):
    """ The numbers the prompt asks for, for one template and one run. """
    g = t[t['good']]
    v_all = summarise(t['velocity'])
    v_good = summarise(g['velocity'])
    per_order = []
    for o in np.unique(g['order']):
        s = g[g['order'] == o]
        per_order.append((o, np.median(s['velocity']), robust_sigma(s['velocity']), len(s)))
    per_order = np.array(per_order)
    return dict(
        n=len(t), n_good=len(g),
        failed=int(np.sum(t['failed'])), no_unc=int(np.sum(t['no_unc'])),
        bad_velocity=int(np.sum(t['bad_velocity'])), at_bound=int(np.sum(t['at_bound'])),
        misfit=int(np.sum(t['misfit'])), v_outlier=int(np.sum(t['v_outlier'])),
        v_median=v_good['median'], v_robust=v_good['robust'], v_std=v_good['std'],
        v_std_all=v_all['std'], v_robust_all=v_all['robust'],
        v_expected=expected,
        sigv_stderr=float(np.nanmedian(g['e_velocity'])),
        sigv_noise=float(np.nanmedian(g['sigv_noise'])),
        redchi_meas=float(np.nanmedian(t['redchi_meas'])),
        order_scatter=float(robust_sigma(per_order[:, 1])) if len(per_order) > 2 else np.nan,
        within_order=float(np.median(per_order[:, 2])) if len(per_order) else np.nan,
        epoch_error=float(v_good['robust'] / np.sqrt(max(len(g), 1))),
        per_order=per_order)


def wave_summary(t):
    """ The fitted wavelength solution against PypeIt's, by order. """
    g = t[t['good']]
    rows = []
    for o in np.unique(g['order']):
        s = g[g['order'] == o]
        # Linear trend of the offset along the order, m/s per 1000 pixels
        slope = np.polyfit(s['pixc'], s['dv_wave_c'], 1)[0] * 1000. if len(s) > 5 else np.nan
        rows.append(dict(order=int(o), dv_median=float(np.median(s['dv_wave_c'])),
                         dv_robust=float(robust_sigma(s['dv_wave_c'])),
                         dv_trend=float(slope),
                         slope_ratio=float(np.median(s['slope_ratio']))))
    return Table(rows)


def continuum_summary(t):
    """ The fitted continuum against PypeIt's blaze and the adapter's envelope.

    Within each order the ratio of the fitted (de-normalised) continuum to
    ``OPT_FLAT`` should be smooth -- the star's energy distribution over the
    lamp's -- so its chunk-to-chunk scatter about a quadratic is the
    continuum's own noise.  The raw ``cont_intercept`` is shown the same way,
    to show what the de-normalisation removes.
    """
    g = t[t['good'] & np.isfinite(t['flat']) & (t['flat'] > 0)]
    rows = []
    for o in np.unique(g['order']):
        s = g[g['order'] == o]
        if len(s) < 8:
            continue
        out = dict(order=int(o))
        for name, col in (('true', 'cont_true'), ('raw', 'cont_intercept')):
            ratio = np.asarray(s[col] / s['flat'], dtype=float)
            ratio /= np.median(ratio)
            fit = np.polyval(np.polyfit(s['pixc'], ratio, 2), s['pixc'])
            out['scatter_' + name] = float(robust_sigma(ratio - fit))
        env = np.asarray(s['cont_true'] / s['cont_env'], dtype=float)
        out['true_over_env'] = float(np.median(env))
        out['true_over_env_scatter'] = float(robust_sigma(env))
        rows.append(out)
    return Table(rows)


def lsf_summary(t, arc):
    """ The fitted LSF against position, beside the same night's ThAr. """
    g = t[t['good'] & np.isfinite(t['lsf_w'])]
    bins = np.linspace(0, 2048, 5)
    rows = []
    for o in np.unique(g['order']):
        s = g[g['order'] == o]
        a = arc[arc['order'] == o]
        for lo, hi in zip(bins[:-1], bins[1:]):
            ss = s[(s['pixc'] >= lo) & (s['pixc'] < hi)]
            aa = a[(a['pixel'] >= lo) & (a['pixel'] < hi)]
            if len(ss) < 3:
                continue
            rows.append(dict(order=int(o), pix_lo=int(lo), pix_hi=int(hi), n=len(ss),
                             fwhm=float(np.median(ss['lsf_w'])),
                             asym=float(np.median(ss['lsf_asym'])),
                             centroid=float(np.median(ss['lsf_cen'])),
                             arc_fwhm=float(np.median(aa['fwhm_px'])) if len(aa) >= 3 else np.nan,
                             n_arc=len(aa)))
    return Table(rows)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(tabs, sums, lsf_tabs, outfile):
    """ Six panels for the two templates' run 1, with run 0 for the LSF. """
    import matplotlib.pyplot as plt
    INK, C1, C2, C3 = '#0b0b0b', '#2a78d6', '#eb6834', '#1baf7a'
    cols = {10: C2, 1: C1}
    fig, axs = plt.subplots(2, 3, figsize=(17, 9.5))

    # (a) chunk velocities along the band
    ax = axs[0, 0]
    for os_ in (10, 1):
        t = tabs[(os_, 1)]
        g = t[t['good']]
        ax.plot(g['wave'], g['velocity'], '.', ms=3, alpha=0.5, color=cols[os_],
                label='osample {:d}: robust sigma {:.0f} m/s'.format(os_, sums[(os_, 1)]['v_robust']))
    exp = sums[(10, 1)]['v_expected']
    ax.axhline(exp, color=INK, lw=1, ls='--', label='expected {:+.0f} m/s'.format(exp))
    ax.set_ylim(exp - 2500, exp + 2500)
    ax.set_xlabel('wavelength (A)')
    ax.set_ylabel('chunk velocity, run 1 (m/s)')
    ax.legend(fontsize=8)
    ax.set_title('(a) per-chunk velocities, good chunks', loc='left', fontsize=10)

    # (b) velocity error: reported vs photon-noise-only vs scatter
    ax = axs[0, 1]
    t = tabs[(10, 1)]
    g = t[t['good']]
    bins = np.logspace(0.5, 4, 50)
    ax.hist(g['e_velocity'], bins=bins, color=C2, alpha=0.6, label='lmfit error (misfit-scaled)')
    ax.hist(g['sigv_noise'], bins=bins, color=C3, alpha=0.6, label='photon-noise error')
    ax.axvline(sums[(10, 1)]['v_robust'], color=INK, ls='--', label='observed robust scatter')
    ax.set_xscale('log')
    ax.set_xlabel('m/s')
    ax.set_ylabel('chunks')
    ax.legend(fontsize=8)
    ax.set_title('(b) what the errors say vs what the chunks do (os 10, run 1)', loc='left', fontsize=10)

    # (c) chi2 by order
    ax = axs[0, 2]
    for os_ in (10, 1):
        t = tabs[(os_, 1)]
        orders = np.unique(t['order'])
        med = [np.median(t['redchi_meas'][t['order'] == o]) for o in orders]
        ax.plot(orders, med, 'o-', color=cols[os_], label='osample {:d}'.format(os_))
    ax.set_yscale('log')
    ax.set_xlabel('echelle order')
    ax.set_ylabel('median chi2/dof, measured noise')
    ax.legend(fontsize=8)
    ax.set_title('(c) the model misfit, by order (run 1)', loc='left', fontsize=10)

    # (d) LSF width vs position (run 0, os 10) and ThAr
    ax = axs[1, 0]
    ls = lsf_tabs[(10, 0)]
    cmap = plt.get_cmap('viridis')
    orders = np.unique(ls['order'])
    for j, o in enumerate(orders):
        s = ls[ls['order'] == o]
        c = cmap(j / max(len(orders) - 1, 1))
        xc = 0.5 * (s['pix_lo'] + s['pix_hi'])
        ax.plot(xc, s['fwhm'], 'o-', color=c, ms=4, lw=1)
        ax.plot(xc, s['arc_fwhm'], 'x--', color=c, ms=4, lw=0.6)
    ax.plot([], [], 'o-', color='0.4', label='star, run 0 fit')
    ax.plot([], [], 'x--', color='0.4', label='ThAr, same night')
    ax.set_xlabel('detector pixel (dispersion)')
    ax.set_ylabel('LSF FWHM (px)')
    ax.legend(fontsize=8)
    ax.set_title('(d) LSF width along the orders (colour: order {:d}-{:d})'.format(
        orders.min(), orders.max()), loc='left', fontsize=10)

    # (e) wavelength offset vs PypeIt
    ax = axs[1, 1]
    for os_ in (10, 1):
        t = tabs[(os_, 1)]
        g = t[t['good']]
        ax.plot(g['wave'], g['dv_wave_c'], '.', ms=3, alpha=0.5, color=cols[os_],
                label='osample {:d}: median {:+.0f} m/s'.format(os_, np.median(g['dv_wave_c'])))
    ax.set_xlabel('wavelength (A)')
    ax.set_ylabel('fitted - PypeIt wavelength (m/s)')
    ax.legend(fontsize=8)
    ax.set_title('(e) wavelength solution vs PypeIt (run 1)', loc='left', fontsize=10)

    # (f) continuum vs PypeIt blaze
    ax = axs[1, 2]
    t = tabs[(10, 1)]
    g = t[t['good'] & np.isfinite(t['flat'])]
    for j, o in enumerate(np.unique(g['order'])):
        s = g[g['order'] == o]
        r = s['cont_true'] / s['flat']
        ax.plot(s['pixc'], r / np.median(r), '-', color=cmap(j / 13.), lw=1)
    ax.set_xlabel('detector pixel (dispersion)')
    ax.set_ylabel('fitted continuum / PypeIt blaze (normalised)')
    ax.set_ylim(0.8, 1.2)
    ax.set_title('(f) continuum vs PypeIt OPT_FLAT (os 10, run 1)', loc='left', fontsize=10)

    fig.tight_layout()
    fig.savefig(outfile, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def parser(options=None):
    p = argparse.ArgumentParser(description='Fit one epoch end to end (phase 3, prompt 5)',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--epoch', type=str, default=EPOCH)
    p.add_argument('--skip-run', action='store_true',
                   help='Analyse existing pyodine results without re-running')
    p.add_argument('--datadir', type=str, default=DATA_DIR)
    p.add_argument('--figdir', type=str, default=FIG_DIR)
    return p.parse_args() if options is None else p.parse_args(options)


def fmt(s):
    return ('  {n:d} chunks: {failed:d} failed, {no_unc:d} without uncertainties, '
            '{bad_velocity:d} bad velocity, {at_bound:d} at a bound, {misfit:d} misfits '
            '(chi2 > {mf:.0f}x order median), {v_outlier:d} velocity outliers '
            '(> {ns:.0f} robust sigma) -> {n_good:d} good').format(mf=MISFIT_FACTOR,
                                                                   ns=OUTLIER_NSIG, **s)


def main(pargs):
    Pars = hires_pars.Parameters()
    k = Pars.noise_scale
    obs_file = f1.epoch_file(pargs.epoch)
    obs = lp.ObservationWrapper(obs_file)
    weight = obs.compute_weight(Pars.weight_type)
    flat = blaze(obs_file, obs)
    iod = lp.IodineTemplate(Pars.i2_to_use)
    wc = os.path.join(lp.DEFAULT_REDUX, 'reduce_' + pargs.epoch[3:11], 'Calibrations',
                      'WaveCalib_A_0_DET01.fits')
    arc = dt.arc_line_widths(wc)
    print('Epoch {:s}: {:s}'.format(pargs.epoch, obs_file))
    print('ThAr, same night: {:d} clean lines, median FWHM {:.3f} px'.format(
        len(arc), np.median(arc['fwhm_px'])))

    tabs, sums, lsf_tabs, wsums, csums, runs = {}, {}, {}, {}, {}, {}
    for os_ in (10, 1):
        template = StellarTemplate_Chunked(dt.TEMPLATE_H5[os_])
        outdir = os.path.join(RUN_ROOT, '{:s}_os{:d}'.format(pargs.epoch, os_))
        print('=== template oversampling {:d} -> {:s}'.format(os_, outdir))
        if not pargs.skip_run:
            info = run_driver(obs_file, dt.TEMPLATE_H5[os_], outdir, Pars)
            runs[os_] = info
            print('  runtime {:.1f} s total; run 0 {:.1f} s, run 1 {:.1f} s (each incl. plots)'.format(
                info['total'], info['t_run'][0], info['t_run'][1]))
            print('  driver log: no uncertainties {}, Chauvenet {}, nan redchi {}'.format(
                info['counts']['no_uncertainties'], info['counts']['chauvenet'],
                info['counts']['nan_redchi']))
        exp = f1.expected_velocity(dict(obs=obs, template=template))
        for run_id in sorted(Pars.model_runs):
            res = os.path.join(outdir, 'run{:d}.h5'.format(run_id))
            t = chunk_table(res, run_id, Pars, obs, template, iod, flat, weight, k)
            t = flag_failures(t, Pars)
            t['osample'] = os_
            s = summarise_run(t, exp['expected'])
            tabs[(os_, run_id)], sums[(os_, run_id)] = t, s
            lsf_tabs[(os_, run_id)] = lsf_summary(t, arc)
            wsums[(os_, run_id)] = wave_summary(t)
            csums[(os_, run_id)] = continuum_summary(t)
            print(' -- run {:d}'.format(run_id))
            print(fmt(s))
            print('  velocity, good chunks: median {:+.1f} m/s (expected {:+.1f}); robust sigma '
                  '{:.1f}, std {:.1f}; all chunks: std {:.0f}, robust {:.1f}'.format(
                      s['v_median'], s['v_expected'], s['v_robust'], s['v_std'],
                      s['v_std_all'], s['v_robust_all']))
            print('  median per-chunk error: lmfit {:.0f} m/s, photon-noise {:.1f} m/s; '
                  'epoch error (robust/sqrt n) {:.1f} m/s'.format(
                      s['sigv_stderr'], s['sigv_noise'], s['epoch_error']))
            print('  order-to-order scatter of order medians {:.1f} m/s; median within-order '
                  'robust sigma {:.1f} m/s; median chi2/dof (measured noise) {:.1f}'.format(
                      s['order_scatter'], s['within_order'], s['redchi_meas']))
            w = wsums[(os_, run_id)]
            print('  wavelength vs PypeIt: median over orders {:+.0f} m/s, order-to-order '
                  '{:.0f} m/s, within-order robust {:.0f} m/s, trend {:+.0f} m/s/kpx, '
                  'dispersion ratio {:.5f}'.format(np.median(w['dv_median']),
                                                   robust_sigma(w['dv_median']),
                                                   np.median(w['dv_robust']),
                                                   np.median(w['dv_trend']),
                                                   np.median(w['slope_ratio'])))
            c = csums[(os_, run_id)]
            print('  continuum vs PypeIt blaze: chunk-to-chunk scatter {:.2%} de-normalised, '
                  '{:.2%} raw; de-normalised / envelope {:.3f} +- {:.3f}'.format(
                      np.median(c['scatter_true']), np.median(c['scatter_raw']),
                      np.median(c['true_over_env']), np.median(c['true_over_env_scatter'])))
            md = misfit_drivers(t)
            print('  misfit vs structure ({:d} chunks, Spearman): iodine {:+.2f}, stellar '
                  '{:+.2f} (iodine-stellar {:+.2f}); partial: iodine {:+.2f}, stellar '
                  '{:+.2f}'.format(md['n'], md['rho_iod'], md['rho_tem'], md['rho_iod_tem'],
                                   md['partial_iod'], md['partial_tem']))
            s['misfit_drivers'] = md
            if run_id == 1:
                print('  order  n   median v   robust sigma   photon sigma   chi2/dof')
                g = t[t['good']]
                for o in np.unique(g['order']):
                    q = g[g['order'] == o]
                    print('  {:3d} {:4d} {:+9.0f} {:12.0f} {:14.0f} {:10.1f}'.format(
                        o, len(q), np.median(q['velocity']), robust_sigma(q['velocity']),
                        np.median(q['sigv_noise']), np.median(q['redchi_meas'])))
            ls = lsf_tabs[(os_, run_id)]
            ok = np.isfinite(ls['arc_fwhm'])
            print('  LSF: FWHM median {:.3f} px (range over bins {:.3f}-{:.3f}), asymmetry '
                  '{:+.3f}, centroid {:+.3f} px; star/ThAr {:.3f}'.format(
                      np.median(ls['fwhm']), ls['fwhm'].min(), ls['fwhm'].max(),
                      np.median(ls['asym']), np.median(ls['centroid']),
                      np.median(ls['fwhm'][ok] / ls['arc_fwhm'][ok])))

    allt = vstack([tabs[key] for key in sorted(tabs)], metadata_conflicts='silent')
    # Six significant figures, except where the value needs more: a
    # wavelength to 1e-5 A is 0.5 m/s
    for c in allt.colnames:
        if allt[c].dtype.kind == 'f':
            allt[c].info.format = '.10g' if c in ('wave', 'wave_intercept', 'e_wave_intercept') \
                else '.6g'
    allt.write(os.path.join(pargs.datadir, 'one_epoch_chunks.csv'), overwrite=True)
    rows = []
    for key in sorted(sums):
        s = dict(sums[key])
        s.pop('per_order')
        for kk, vv in s.pop('misfit_drivers').items():
            s['misfit_' + kk] = vv
        s.update(osample=key[0], run=key[1])
        if key[0] in runs:
            s.update(runtime_total=runs[key[0]]['total'], runtime_run=runs[key[0]]['t_run'][key[1]])
        rows.append(s)
    Table(rows).write(os.path.join(pargs.datadir, 'one_epoch_summary.csv'), overwrite=True)
    lsfall = []
    for key in sorted(lsf_tabs):
        tt = lsf_tabs[key]
        tt['osample'], tt['run'] = key
        lsfall.append(tt)
    vstack(lsfall).write(os.path.join(pargs.datadir, 'one_epoch_lsf.csv'), overwrite=True)
    fig = os.path.join(pargs.figdir, 'fig_p5_one_epoch.png')
    make_figure(tabs, sums, lsf_tabs, fig)
    print('Figure: {:s}'.format(fig))


if __name__ == '__main__':
    main(parser())
