""" Fit all 31 cell-in epochs and build the velocity table (phase 3, prompt 7).

Every cell-in frame of HD 187123 goes through `pyodine`'s own driver,
`model_single_observation`, exactly as prompt 5 ran one epoch, with the
configuration prompt 6 settled: run 0 a single Gaussian, run 1 the
super-Gaussian, the oversampling-1 deconvolved template, weights from
PypeIt's inverse variance, the Fischer atlas at the adapter's alpha with a
free per-chunk ``iod_depth`` (the Q&A before prompt 7: the current model,
answer (a)).

WHICH ORDERS.  Each epoch is fitted over the iodine orders phase 2's quality
filter marks usable (``data/order_quality.csv``), passed to the driver as its
``orders`` argument.  Before fork change 5 that silently paired chunks after a
rejected order with the wrong template chunks; since, it is the intended use.
Result chunks are mapped back to template chunks by order and position in the
order (`template_columns`), never by list position.

COMBINING CHUNKS.  An epoch's velocity is `pyodine`'s own combination,
`timeseries.combine_vels.combine_chunk_velocities` -- the iSONG algorithm:
each chunk's time series is offset to the epoch means, weighted by its own
scatter over the time series, and re-weighted per epoch by how far it
deviates.  It takes an (epoch x chunk) array with NaN where a chunk was not
fitted; rejected orders are NaN.  Two simpler estimators are kept beside it
so the effect of each step can be seen: the plain median of an epoch's good
chunks (prompt 5's estimator) and the median after removing each chunk's
offset (`pyodine`'s ``mdvel``).

VELOCITY FRAME.  `pyodine`'s chunk velocity is the observed-frame shift of
the template onto the epoch.  In that frame the star sits at v - BVC, so the
velocity relative to the template in the barycentre is

    v_bary = v + BVC_epoch - BVC_template

additively, as phase 2 prompt 5 did (the multiplicative cross-term is < 2 m/s
here).  Its zero point is the template's plus the combination's chunk-offset
convention; only differences between epochs mean anything.

THE SCATTERS, separately (as phase 2 prompt 5):
  * per chunk -- the robust scatter of an epoch's offset-corrected chunk
    velocities;
  * per order -- the scatter of an epoch's order means about each other;
  * epoch to epoch -- the rms of the epoch velocities over the baseline, and
    its within-night and between-night parts.

Run with (the fits, in shards, then the analysis):

    conda run --no-capture-output -n pypeit14 \\
        python -m first_hires_exoplanet.fit_all_epochs --run --shard 0 --nshard 3
    ...
    conda run --no-capture-output -n pypeit14 \\
        python -m first_hires_exoplanet.fit_all_epochs

"""

# Must precede any import of matplotlib, including pyodine's own
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# Standard imports
import time
import argparse

import numpy as np
from astropy.table import Table, vstack

from first_hires_exoplanet import utilities_hires  # noqa: F401 (puts pyodine on the path)
from first_hires_exoplanet.utilities_hires import load_pyodine as lp
from first_hires_exoplanet.utilities_hires import pyodine_parameters as hires_pars
from first_hires_exoplanet import deconvolve_template as dt
from first_hires_exoplanet import fit_one_chunk as f1
from first_hires_exoplanet import fit_one_epoch as f5
from first_hires_exoplanet.pyodine_chunk_stats import robust_sigma

from pyodine.template.base import StellarTemplate_Chunked
from pyodine.timeseries import combine_vels


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR = f5.DATA_DIR
FIG_DIR = f5.FIG_DIR
RUN_DIR = os.path.join(f5.RUN_ROOT, 'epochs')
OSAMPLE = 1
RUN = 1                                     # the run whose velocities are used

#: `pyodine`'s default (SONG) weighting parameters, copied: the module-level
#: default dict is mutated by `combine_chunk_velocities` when it is used
#: implicitly (``del pars['good_orders']``), so a second call would fail.
WEIGHTING = dict(combine_vels._weighting_pars)

#: A pair of frames within a night closer than this (hours) is treated as
#: seeing the same stellar velocity: the planet moves at most
#: 2 pi K / P = 6 m/s per hour
PAIR_HOURS = 2.


# ---------------------------------------------------------------------------
# The epochs
# ---------------------------------------------------------------------------

def epochs():
    """ The cell-in science frames that phase 2 reduced, in time order. """
    xc = Table.read(os.path.join(DATA_DIR, 'xcorr_velocities.csv'))
    sci = Table.read(os.path.join(DATA_DIR, 'koa_hd187123_science.csv'))
    cell = {str(k).replace('.fits', ''): str(v).strip() == 'T'
            for k, v in zip(sci['koaid'], sci['iodin'])}
    out = [str(k) for k in xc['koaid'] if cell.get(str(k), False)]
    return sorted(out, key=lambda k: float(xc['mjd'][list(xc['koaid']).index(k)]))


_TEMPLATE_ECH = {}


def template_echelle_orders():
    """ Template order -> echelle order, from a frame with the full 37 orders.

    The template was built from 1998-08-26 frames with the full format, and
    its order numbers are positions in that format.
    """
    if not _TEMPLATE_ECH:
        ref = lp.ObservationWrapper(f1.epoch_file(f5.EPOCH))
        if len(ref.ech_orders) != 37:
            raise RuntimeError('reference epoch does not have the full format')
        _TEMPLATE_ECH.update({p: int(e) for p, e in enumerate(ref.ech_orders)})
    return _TEMPLATE_ECH


def usable_orders(koaid, obs, template):
    """ Template orders (the driver's ``orders``) whose echelle order is usable.

    `pyodine` maps template order ``o`` to observation position
    ``o + order_correction``: one constant shift, which assumes the
    observation's orders are contiguous.  1998-09-13 was extracted with 34
    orders and no order 71, so from that gap on the shift is wrong and
    template order 71 lands on echelle order 72.  A template order is used
    only if its position really holds the same echelle order, and phase 2's
    quality filter passes it.

    Returns:
        tuple: template orders, echelle orders rejected, order correction
    """
    q = Table.read(os.path.join(DATA_DIR, 'order_quality.csv'))
    q = q[q['koaid'] == koaid]
    use = {int(o): bool(u) for o, u in zip(q['order'], q['use'])}
    tech = template_echelle_orders()
    tind = template.get_order_indices(template[0].order)
    obs_order_min, _ = obs.check_wavelength_range(template[tind[0]].w0, template[tind[-1]].w0)
    corr = int(obs_order_min - template[tind[0]].order)
    keep, rejected = [], []
    for to in template.orders_unique:
        want = tech[int(to)]
        pos = int(to) + corr
        have = int(obs.ech_orders[pos]) if 0 <= pos < len(obs.ech_orders) else -1
        if have == want and use.get(want, False):
            keep.append(int(to))
        else:
            rejected.append(want)
    return np.array(keep), rejected, corr


def outdir(koaid):
    return os.path.join(RUN_DIR, '{:s}_os{:d}'.format(koaid, OSAMPLE))


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------

def done(koaid):
    """ Both result files written and no traceback in the info log. """
    d = outdir(koaid)
    info = os.path.join(d, 'info.log')
    if not all(os.path.exists(os.path.join(d, 'run{:d}.h5'.format(r))) for r in (0, 1)) \
            or not os.path.exists(info):
        return False
    with open(info) as fh:
        return 'Traceback' not in fh.read()


def run_epochs(koaids, resume=False):
    template = StellarTemplate_Chunked(dt.TEMPLATE_H5[OSAMPLE])
    rows = []
    for j, koaid in enumerate(koaids):
        if resume and done(koaid):
            print('[{:d}/{:d}] {:s}: done'.format(j + 1, len(koaids), koaid), flush=True)
            continue
        Pars = hires_pars.Parameters()
        obs_file = f1.epoch_file(koaid)
        obs = lp.ObservationWrapper(obs_file)
        orders, rejected, corr = usable_orders(koaid, obs, template)
        # The driver's example-chunk plots index the chunk list directly, and
        # an epoch with few usable orders has fewer chunks than they name
        nchunk = sum(len(template.get_order_indices(o)) for o in orders)
        for run in Pars.model_runs.values():
            if run.get('plot_chunks') is not None:
                run['plot_chunks'] = [c for c in run['plot_chunks'] if c < nchunk]
        print('[{:d}/{:d}] {:s}: {:d} orders (rejected echelle orders {}), order correction {:d}'.format(
            j + 1, len(koaids), koaid, len(orders), rejected, corr), flush=True)
        t0 = time.time()
        info = f5.run_driver(obs_file, dt.TEMPLATE_H5[OSAMPLE], outdir(koaid), Pars,
                             orders=orders)
        print('   {:.1f} s'.format(time.time() - t0), flush=True)
        rows.append(dict(koaid=koaid, runtime=info['total'], n_orders=len(orders)))
    return rows


# ---------------------------------------------------------------------------
# Per-epoch chunk tables
# ---------------------------------------------------------------------------

def template_columns(t, template, corr):
    """ The template chunk index of every result chunk.

    Result chunks carry their observation order position and first pixel.
    Within an order they are in pixel order, and `auto_wave_comoving` makes
    exactly as many as the template has in that order, so the k-th chunk of
    observation position p is the k-th template chunk of template order
    p - corr.
    """
    out = np.full(len(t), -1, dtype=int)
    for pos in np.unique(t['pos']):
        sel = np.where(t['pos'] == pos)[0]
        sel = sel[np.argsort(np.asarray(t['pix0'][sel]))]
        tind = template.get_order_indices(int(pos) - corr)
        if len(tind) != len(sel):
            raise RuntimeError('order position {:d}: {:d} chunks against {:d} in the template'.format(
                int(pos), len(sel), len(tind)))
        out[sel] = tind
    return out


def epoch_chunks(koaid, Pars, template, iod):
    """ Prompt 5's per-chunk table for one epoch's run 1, with template columns. """
    obs_file = f1.epoch_file(koaid)
    obs = lp.ObservationWrapper(obs_file)
    weight = obs.compute_weight(Pars.weight_type)
    flat = f5.blaze(obs_file, obs)
    _, _, corr = usable_orders(koaid, obs, template)
    res = os.path.join(outdir(koaid), 'run{:d}.h5'.format(RUN))
    # First pass for the positions only, then the table with the right template
    from pyodine.fitters import results_io
    r = results_io.load_results(res)
    pos = Table(dict(pos=np.asarray(r['chunks']['order'], dtype=int),
                     pix0=np.asarray(r['chunks']['firstpix'], dtype=int)))
    tcol = template_columns(pos, template, corr)
    t = f5.chunk_table(res, RUN, Pars, obs, template, iod, flat, weight, Pars.noise_scale,
                       tmpl_index=tcol)
    t = f5.flag_failures(t, Pars)
    t['tchunk'] = tcol
    t['koaid'] = koaid
    meta = dict(koaid=koaid, bary_date=float(obs.bary_date), bvc=float(obs.bary_vel_corr),
                night=koaid[3:11], corr=corr)
    return t, meta


# ---------------------------------------------------------------------------
# Combination
# ---------------------------------------------------------------------------

def velocity_matrix(tabs, nchunk, column='velocity', good_only=False):
    """ (epoch x template chunk) array, NaN where a chunk was not fitted or failed. """
    V = np.full((len(tabs), nchunk), np.nan)
    for i, t in enumerate(tabs):
        ok = ~t['bad_velocity'] & ~t['failed']
        if good_only:
            ok &= t['good']
        V[i, np.asarray(t['tchunk'][ok])] = np.asarray(t[column][ok], dtype=float)
    return V


#: A chunk needs this many epochs to define its offset and its time-series
#: scatter.  `pyodine`'s `robust.mean` divides by zero on a column with one
#: finite value (``sqrt(len(good) - 1)``), so thinner columns are blanked.
MIN_EPOCHS = 3


def combine(V, nr_orders):
    """ `pyodine`'s chunk combination, quietly, on chunks seen in enough epochs. """
    import logging
    logging.getLogger().setLevel(logging.WARNING)
    V = V.copy()
    thin = np.sum(np.isfinite(V), axis=0) < MIN_EPOCHS
    V[:, thin] = np.nan
    if thin.any():
        print('  combine: {:d} chunks seen in < {:d} epochs blanked'.format(int(thin.sum()), MIN_EPOCHS))
    rv, aux, _ = combine_vels.combine_chunk_velocities(V, nr_orders,
                                                        weighting_pars=dict(WEIGHTING))
    return rv, aux


def per_order(V, offsets, torder, weights=None):
    """ Each epoch's offset-corrected velocity per order (median of its chunks). """
    Vc = V - offsets[None, :]
    orders = np.unique(torder)
    out = np.full((V.shape[0], len(orders)), np.nan)
    nn = np.zeros_like(out, dtype=int)
    for j, o in enumerate(orders):
        sel = torder == o
        with np.errstate(all='ignore'):
            out[:, j] = np.nanmedian(Vc[:, sel], axis=1)
        nn[:, j] = np.sum(np.isfinite(Vc[:, sel]), axis=1)
    return orders, out, nn


def night_pairs(v, err, jd, night, hours=PAIR_HOURS):
    """ Differences of frames on one night closer than ``hours``. """
    rows = []
    for i in range(len(v)):
        for j in range(i + 1, len(v)):
            dt_h = abs(jd[j] - jd[i]) * 24.
            if night[i] == night[j] and dt_h < hours:
                rows.append(dict(i=i, j=j, dt_h=dt_h, dv=v[j] - v[i],
                                 sig=np.hypot(err[i], err[j])))
    return Table(rows) if rows else None


def excess_from_pairs(pairs):
    """ Extra per-epoch scatter the internal errors do not account for.

    sigma_extra^2 = < dv^2 / 2 - (sig_i^2 + sig_j^2) / 2 >, floored at zero.
    """
    if pairs is None or len(pairs) == 0:
        return np.nan
    val = np.mean(pairs['dv'] ** 2 / 2. - pairs['sig'] ** 2 / 2.)
    return float(np.sqrt(max(val, 0.)))


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(vt, W, aux, outfile):
    import matplotlib.pyplot as plt
    INK, C1, C2, C3, GR = '#0b0b0b', '#2a78d6', '#eb6834', '#1baf7a', '#8a8a8a'
    fig, axs = plt.subplots(2, 2, figsize=(14, 9.5))

    # (a) epoch velocities in time, both methods
    ax = axs[0, 0]
    t0 = 2450000.
    ax.errorbar(vt['jd'] - t0, vt['v_xcorr'] - np.median(vt['v_xcorr']), yerr=vt['sig_xcorr'],
                fmt='s', ms=4, color=GR, alpha=0.7, label='phase 2 cross-correlation')
    ax.errorbar(vt['jd'] - t0, vt['v_bary'] - np.median(vt['v_bary']), yerr=vt['sig_epoch'],
                fmt='o', ms=5, color=C1, label='iodine forward model (this prompt)')
    ax.set_xlabel('JD - 2450000')
    ax.set_ylabel('relative barycentric velocity (m/s)')
    ax.legend(fontsize=8)
    ax.set_title('(a) the 31 cell-in epochs', loc='left', fontsize=10)

    # (b) the scatter hierarchy
    ax = axs[0, 1]
    bins = np.logspace(0, 3.7, 40)
    ax.hist(vt['c2c_scatter'], bins=bins, color=C2, alpha=0.6, label='per chunk (robust, within epoch)')
    ax.hist(vt['order_scatter'], bins=bins, color=C3, alpha=0.6, label='per order (between order means)')
    ax.hist(vt['sig_epoch'], bins=bins, color=C1, alpha=0.6, label='adopted epoch error')
    ax.axvline(np.std(vt['v_bary']), color=INK, ls='--', label='epoch-to-epoch rms')
    ax.set_xscale('log')
    ax.set_xlabel('m/s')
    ax.set_ylabel('epochs')
    ax.legend(fontsize=8)
    ax.set_title('(b) the three scatters', loc='left', fontsize=10)

    # (c) iodine vs cross-correlation difference against the wavelength zero point
    ax = axs[1, 0]
    d = vt['v_xcorr'] - vt['v_bary']
    ax.plot(-vt['dv_wave'], d - np.median(d), 'o', color=C1)
    lim = np.array([-1500, 2500])
    ax.plot(lim - np.median(-vt['dv_wave']), lim - np.median(-vt['dv_wave']), '-', color=GR, lw=1,
            label='one to one')
    ax.set_xlabel('PypeIt - iodine wavelength zero point (m/s)')
    ax.set_ylabel('cross-correlation - iodine velocity (m/s, offset)')
    ax.legend(fontsize=8)
    ax.set_title('(c) what phase 2 got wrong is the ThAr zero point', loc='left', fontsize=10)

    # (d) chunk sigma over the time series against wavelength
    ax = axs[1, 1]
    with np.errstate(all='ignore'):
        wave = np.nanmedian(W, axis=0)
    ax.plot(wave, aux['chunk_sigma'], '.', ms=3, color=C2)
    ax.set_yscale('log')
    ax.set_xlabel('wavelength (A)')
    ax.set_ylabel('chunk time-series sigma (m/s)')
    ax.set_title('(d) how well each chunk repeats, epoch to epoch', loc='left', fontsize=10)

    fig.tight_layout()
    fig.savefig(outfile, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def analyse(koaids, pargs):
    Pars = hires_pars.Parameters()
    template = StellarTemplate_Chunked(dt.TEMPLATE_H5[OSAMPLE])
    iod = lp.IodineTemplate(Pars.i2_to_use)
    nchunk = len(template)
    torder = np.array([template[i].order for i in range(nchunk)])
    bvc_t = float(template.bary_vel_corr)

    tabs, metas = [], []
    for koaid in koaids:
        t, m = epoch_chunks(koaid, Pars, template, iod)
        tabs.append(t)
        metas.append(m)
        print('  {:s}: {:d} chunks, {:d} good'.format(koaid, len(t), int(np.sum(t['good']))), flush=True)

    V_all = velocity_matrix(tabs, nchunk)
    V_good = velocity_matrix(tabs, nchunk, good_only=True)
    W = velocity_matrix(tabs, nchunk, column='wave')
    rv, aux = combine(V_all, len(np.unique(torder)))
    rv_g, aux_g = combine(V_good, len(np.unique(torder)))

    bvc = np.array([m['bvc'] for m in metas])
    jd = np.array([m['bary_date'] for m in metas])
    night = np.array([m['night'] for m in metas])
    dbvc = bvc - bvc_t

    # Per-order velocities, and their scatter within each epoch
    orders, vo, no = per_order(V_all, aux['chunk_offsets'], torder)
    # Template order -> echelle order, from the fitted chunks
    ech = {}
    for t in tabs:
        for e, tc in zip(t['order'], t['tchunk']):
            ech[int(torder[tc])] = int(e)

    # Plain median of good chunks (prompt 5), and the photon-noise floor
    v_med = np.array([np.median(t['velocity'][t['good']]) for t in tabs])
    sig_phot = np.array([1. / np.sqrt(np.sum(1. / np.asarray(t['sigv_noise'][t['good']]) ** 2))
                         for t in tabs])
    n_good = np.array([int(np.sum(t['good'])) for t in tabs])
    n_fit = np.array([len(t) for t in tabs])
    n_ord = np.sum(no > 0, axis=1)
    order_scatter = np.array([robust_sigma(vo[i][np.isfinite(vo[i])]) for i in range(len(tabs))])
    order_std = np.array([np.nanstd(vo[i], ddof=1) for i in range(len(tabs))])
    sig_order = order_std / np.sqrt(n_ord)
    dv_wave = np.array([np.median(t['dv_wave_c'][t['good']]) for t in tabs])
    redchi = np.array([np.nanmedian(t['redchi_meas']) for t in tabs])

    xc = Table.read(os.path.join(DATA_DIR, 'xcorr_velocities.csv'))
    xmap = {str(k): i for i, k in enumerate(xc['koaid'])}
    xi = np.array([xmap[k] for k in koaids])

    v_bary = rv['rv'] + dbvc
    sig_int = np.maximum(rv['rv_err'], sig_order)
    pairs = night_pairs(v_bary, sig_int, jd, night)
    extra = excess_from_pairs(pairs)
    sig_epoch = np.hypot(sig_int, extra if np.isfinite(extra) else 0.)

    vt = Table()
    vt['koaid'] = koaids
    vt['night'] = night
    vt['jd'] = jd
    vt['mjd'] = jd - 2400000.5
    vt['v_bary'] = v_bary
    vt['sig_epoch'] = sig_epoch
    vt['sig_internal'] = sig_int
    vt['sig_pyodine'] = rv['rv_err']
    vt['sig_order'] = sig_order
    vt['sig_photon'] = sig_phot
    vt['v_obs'] = rv['rv']
    vt['v_bary_good'] = rv_g['rv'] + dbvc
    vt['v_bary_mdvel'] = rv['mdvel'] + dbvc
    vt['v_bary_median'] = v_med + dbvc
    vt['bvc'] = bvc
    vt['bvc_template'] = bvc_t
    vt['c2c_scatter'] = rv['c2c_scatter']
    vt['order_scatter'] = order_scatter
    vt['n_orders'] = n_ord
    vt['n_chunks'] = n_fit
    vt['n_good'] = n_good
    vt['redchi_meas'] = redchi
    vt['dv_wave'] = dv_wave
    vt['v_xcorr'] = np.asarray(xc['v_rel'][xi], dtype=float)
    vt['sig_xcorr'] = np.asarray(xc['sigma_empirical'][xi], dtype=float)
    for c in vt.colnames:
        if vt[c].dtype.kind == 'f':
            vt[c].info.format = '.4f' if c in ('jd', 'mjd') else '.2f'
    vt.meta['comments'] = [
        'HD 187123, 31 cell-in HIRES epochs, pyodine forward model (phase 3 prompt 7).',
        'v_bary: relative barycentric velocity, m/s, arbitrary zero point (template + chunk offsets).',
        'sig_epoch: adopted error = hypot(max(sig_pyodine, sig_order), sigma_extra from close night pairs).',
    ]
    vt.write(os.path.join(pargs.datadir, 'iodine_velocities.csv'), overwrite=True)

    ot = Table()
    ot['koaid'] = np.repeat(koaids, len(orders))
    ot['order'] = np.tile([ech.get(int(o), -1) for o in orders], len(koaids))
    ot['v_bary'] = (vo + dbvc[:, None]).ravel()
    ot['n_chunks'] = no.ravel()
    ot = ot[ot['n_chunks'] > 0]
    ot['v_bary'].info.format = '.2f'
    ot.write(os.path.join(pargs.datadir, 'iodine_velocities_per_order.csv'), overwrite=True)

    keep = ['koaid', 'order', 'tchunk', 'pixc', 'wave', 'velocity', 'e_velocity', 'sigv_noise',
            'redchi_meas', 'iod_depth', 'lsf_w', 'dv_wave_c', 'no_unc', 'misfit', 'v_outlier',
            'good']
    ct = vstack([t[keep] for t in tabs], metadata_conflicts='silent')
    ct['sigma_chunk'] = aux['chunk_sigma'][np.asarray(ct['tchunk'])]
    ct['offset_chunk'] = aux['chunk_offsets'][np.asarray(ct['tchunk'])]
    for c in ct.colnames:
        if ct[c].dtype.kind == 'f':
            ct[c].info.format = '.8g' if c == 'wave' else '.5g'
    ct.write(os.path.join(pargs.datadir, 'iodine_chunks.csv'), overwrite=True)
    # The full per-chunk tables (every fitted parameter) stay beside the runs
    vstack(tabs, metadata_conflicts='silent').write(
        os.path.join(RUN_DIR, 'all_epochs_chunks_run{:d}.fits'.format(RUN)), overwrite=True)

    report(vt, pairs, extra, aux, aux_g, tabs, V_all, torder, orders, vo, ech)
    make_figure(vt, W, aux, os.path.join(pargs.figdir, 'fig_p7_epochs.png'))
    return vt


def report(vt, pairs, extra, aux, aux_g, tabs, V, torder, orders, vo, ech):
    v = np.asarray(vt['v_bary'])
    print()
    print('=== {:d} epochs, {:d} nights'.format(len(vt), len(np.unique(vt['night']))))
    print('chunks fitted {:d}, good {:d}; median chi2/dof (measured noise) {:.1f}'.format(
        int(np.sum(vt['n_chunks'])), int(np.sum(vt['n_good'])), np.nanmedian(vt['redchi_meas'])))
    print()
    print('PER CHUNK, within an epoch (robust, offset-corrected): median {:.0f} m/s, range {:.0f}-{:.0f}'.format(
        np.median(vt['c2c_scatter']), vt['c2c_scatter'].min(), vt['c2c_scatter'].max()))
    raw = [robust_sigma(np.asarray(t['velocity'][t['good']])) for t in tabs]
    off = aux['chunk_offsets']
    cor = [robust_sigma(np.asarray(t['velocity'][t['good']]) - off[np.asarray(t['tchunk'][t['good']])])
           for t in tabs]
    print('   good chunks only: {:.0f} m/s without chunk offsets (prompt 5 statistic), {:.0f} m/s '
          'with them removed'.format(np.median(raw), np.nanmedian(cor)))
    print('   chunk time-series sigma: median {:.0f} m/s (the chunk-offset model is what removes '
          'the common misfit)'.format(np.nanmedian(aux['chunk_sigma'])))
    print('   chunk offsets: robust spread {:.0f} m/s'.format(robust_sigma(aux['chunk_offsets'][
        np.isfinite(aux['chunk_offsets'])])))
    print('PER ORDER, within an epoch (order means about each other): median robust {:.0f} m/s, '
          'range {:.0f}-{:.0f}'.format(np.median(vt['order_scatter']), vt['order_scatter'].min(),
                                        vt['order_scatter'].max()))
    print('EPOCH TO EPOCH: rms {:.1f} m/s, robust {:.1f} m/s over {:d} epochs'.format(
        np.std(v, ddof=1), robust_sigma(v), len(v)))
    nights = np.unique(vt['night'])
    nm = np.array([np.mean(v[vt['night'] == n]) for n in nights])
    within = [np.std(v[vt['night'] == n], ddof=1) for n in nights if np.sum(vt['night'] == n) > 1]
    spread = [np.ptp(v[vt['night'] == n]) for n in nights if np.sum(vt['night'] == n) > 1]
    print('   within a night ({:d} nights with >1 frame): median rms {:.1f}, worst spread {:.1f} m/s'.format(
        len(within), np.median(within), np.max(spread)))
    print('   between nights ({:d} nightly means): rms {:.1f} m/s'.format(len(nm), np.std(nm, ddof=1)))
    for est in ('v_bary_good', 'v_bary_mdvel', 'v_bary_median'):
        print('   other estimator {:s}: rms {:.1f}, robust {:.1f}'.format(
            est, np.std(vt[est], ddof=1), robust_sigma(np.asarray(vt[est]))))
    print('   phase 2 cross-correlation, same epochs: rms {:.1f} m/s'.format(np.std(vt['v_xcorr'], ddof=1)))
    print()
    print('ERRORS (median): photon {:.1f}; pyodine rv_err {:.1f}; order-based {:.1f}; '
          'internal {:.1f}; adopted {:.1f} m/s'.format(
              np.median(vt['sig_photon']), np.median(vt['sig_pyodine']), np.median(vt['sig_order']),
              np.median(vt['sig_internal']), np.median(vt['sig_epoch'])))
    if pairs is not None:
        print('   close night pairs (< {:.0f} h): {:d}; rms difference {:.1f} m/s against {:.1f} '
              'expected from internal errors -> sigma_extra {:.1f} m/s'.format(
                  PAIR_HOURS, len(pairs), np.sqrt(np.mean(pairs['dv'] ** 2)),
                  np.sqrt(np.mean(pairs['sig'] ** 2)), extra))
        for p in pairs:
            print('     {:s} - {:s}: {:5.2f} h, dv {:+7.1f} +- {:5.1f}'.format(
                vt['koaid'][p['j']], vt['koaid'][p['i']], p['dt_h'], p['dv'], p['sig']))
    # A magnitude check only (trap 8), not the assessment, which is prompt 8's:
    # the modern catalogue's circular orbit, offset free
    from first_hires_exoplanet import figs_phase2
    cat = figs_phase2.load_catalogue()
    orb = figs_phase2.fit_circular(np.asarray(cat['bjd']), np.asarray(cat['rv']),
                                   e=np.asarray(cat['e_rv']))
    model = orb['model'](np.asarray(vt['jd']))
    res = v - model
    res -= np.mean(res)
    print('MAGNITUDE CHECK: minus the catalogue orbit (K = {:.1f} m/s), offset only: rms {:.1f} m/s, '
          'chi2/dof with sig_epoch {:.2f} (prompt 8 assesses this properly)'.format(
              orb['K'], np.std(res, ddof=1), np.sum((res / vt['sig_epoch']) ** 2) / (len(res) - 1)))
    for est in ('v_bary_good', 'v_bary_mdvel', 'v_bary_median', 'v_xcorr'):
        r = np.asarray(vt[est]) - model
        print('   the same for {:s}: rms {:.1f} m/s'.format(est, np.std(r - np.mean(r), ddof=1)))
    print()
    d = np.asarray(vt['v_xcorr'] - vt['v_bary'])
    ok = np.isfinite(d)
    r = np.corrcoef(d[ok], -np.asarray(vt['dv_wave'])[ok])[0, 1]
    fit = np.polyfit(-np.asarray(vt['dv_wave'])[ok], d[ok], 1)
    print('CROSS-CHECK: cross-correlation - iodine vs PypeIt - iodine wavelength zero point: '
          'r = {:+.3f}, slope {:.3f}; residual rms {:.1f} m/s'.format(
              r, fit[0], np.std(d[ok] - np.polyval(fit, -np.asarray(vt['dv_wave'])[ok]), ddof=2)))
    print()
    print('koaid              night     n_ord n_good  v_bary  sig_ep  sig_int sig_ph  c2c  ord_sc  chi2  dv_wave  v_xcorr')
    for row in vt:
        print('{:s} {:s} {:3d} {:5d} {:+8.1f} {:6.1f} {:7.1f} {:6.1f} {:5.0f} {:6.0f} {:5.1f} {:+7.0f} {:+8.1f}'.format(
            row['koaid'], row['night'], row['n_orders'], row['n_good'], row['v_bary'], row['sig_epoch'],
            row['sig_internal'], row['sig_photon'], row['c2c_scatter'], row['order_scatter'],
            row['redchi_meas'], row['dv_wave'], row['v_xcorr']))
    print()
    print('PER ORDER over the time series (order velocity - epoch velocity):')
    for j, o in enumerate(orders):
        dd = vo[:, j] + (vt['v_bary'] - vt['v_obs']) - vt['v_bary']
        dd = dd[np.isfinite(dd)]
        print('  order {:3d}: {:2d} epochs, rms {:6.1f}, robust {:6.1f} m/s'.format(
            ech.get(int(o), -1), len(dd), np.std(dd, ddof=1), robust_sigma(dd)))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def parser(options=None):
    p = argparse.ArgumentParser(description='Fit all cell-in epochs (phase 3, prompt 7)',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--run', action='store_true', help='Fit the epochs (else analyse)')
    p.add_argument('--shard', type=int, default=0)
    p.add_argument('--nshard', type=int, default=1)
    p.add_argument('--resume', action='store_true', help='Skip epochs already fitted cleanly')
    p.add_argument('--only', type=str, nargs='*', default=None, help='KOAIDs to fit')
    p.add_argument('--datadir', type=str, default=DATA_DIR)
    p.add_argument('--figdir', type=str, default=FIG_DIR)
    return p.parse_args() if options is None else p.parse_args(options)


def main(pargs):
    koaids = epochs()
    print('{:d} cell-in epochs'.format(len(koaids)))
    if pargs.run:
        todo = pargs.only if pargs.only else koaids[pargs.shard::pargs.nshard]
        rows = run_epochs(todo, resume=pargs.resume)
        Table(rows).write(os.path.join(RUN_DIR, 'runtimes_shard{:d}.csv'.format(pargs.shard)),
                          overwrite=True)
        return
    analyse(pargs.only if pargs.only else koaids, pargs)


if __name__ == '__main__':
    main(parser())
