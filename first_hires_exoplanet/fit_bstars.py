""" The iodine cell seen through B stars (phase 3, after prompt 6).

The Q&A after prompt 5 decided (question 1, option b): before building on an
LSF chosen against a known-wrong atlas, fit the B stars with the atlas free to
rescale per chunk, then compare LSF models there.  A rapidly rotating B star
is almost featureless, so through the cell it is -- to first order -- the
iodine spectrum of *the HIRES cell*, at the stellar LSF, on the science
night.  That separates the two ways the atlas could be wrong:

* its line **depths** (one Beer-Lambert alpha for the whole band): then a
  per-chunk alpha fits the B stars to the noise, and fixing each chunk's
  alpha in the stellar fit removes prompt 5's iodine correlation;
* its line **pattern** (a different cell): then even per-chunk alpha leaves
  a misfit that follows the iodine, in the B stars as in the star.

The fitting is `pyodine`'s own hot-star path, as `create_template` uses it:
the template's chunk grid (so B-star chunk *i* is template chunk *i* is epoch
chunk *i*), `SimpleModel` with no stellar template, and
``Template_Parameters``' runs and constraints (velocity and template depth
fixed; ``iod_depth`` free per chunk).  ``iod_depth`` is an exponent on top of
the adapter's alpha = 2.59, so a chunk's effective alpha is 2.59 x
``iod_depth``.

Four parts:

  1. every B-star frame, run 0 and run 1 (the settled super-Gaussian): the
     alpha map, its frame-to-frame consistency, the B-star chi-square and
     what it follows, the LSF, and the wavelength solution against PypeIt;
  2. the noise of these frames, from two consecutive HR 8634 exposures
     differenced pixel by pixel (as prompt 3 did for the template);
  3. the five LSF models on the two HR 8634 frames nearest the template;
  4. the payoff: prompt 5's epoch refitted with each chunk's ``iod_depth``
     fixed to its B-star value, against the same fit with it free.

Run with:

    conda run --no-capture-output -n pypeit14 \\
        python -m first_hires_exoplanet.fit_bstars

"""

# Must precede any import of matplotlib, including pyodine's own
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# Standard imports
import glob
import time
import argparse
import logging

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

from first_hires_exoplanet.utilities_hires import load_pyodine as lp
from first_hires_exoplanet.utilities_hires import pyodine_parameters as hp
from first_hires_exoplanet import deconvolve_template as dt
from first_hires_exoplanet import fit_one_chunk as f1
from first_hires_exoplanet import fit_one_epoch as fe
from first_hires_exoplanet import compare_lsf_models as cl
from first_hires_exoplanet.pyodine_chunk_stats import robust_sigma

import pyodine
import pipe_lib
from pyodine.fitters.lmfit_wrapper import LmfitWrapper
from pyodine.models import lsf as plsf
from pyodine.models.base import ParameterSet
from pyodine.template.base import StellarTemplate_Chunked


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_HERE, 'data')
FIG_DIR = os.path.join(os.path.dirname(_HERE), 'docs', 'figs')
BSTAR_SCIENCE = os.path.join(lp.DEFAULT_REDUX, 'reduce_19980826_lsf', 'Science')
CKMS = 299792.458

#: The pair nearest the template exposures (09:39-09:59 UT), consecutive
PAIR = ('HI.19980826.37213', 'HI.19980826.37323')


# ---------------------------------------------------------------------------
# Fitting a B star: pyodine's hot-star path
# ---------------------------------------------------------------------------

#: PypeIt extracted the B stars with `force_center_obj` (a boxcar across each
#: order); their optimal-extraction mask collapses on a slit-filling star, the
#: boxcar's does not (99% of iodine-order pixels kept).  See `reduce_bstars`.
EXTRACTION = 'BOX'

#: Frames whose median S/N over the iodine orders is below this are not
#: fitted.  HR 838 at 20 s (55734) peaks at 773 ADU: the star missed the slit.
MIN_SNR = 50.


def bstar_files():
    return sorted(glob.glob(os.path.join(BSTAR_SCIENCE, 'spec1d_*.fits')))


def open_bstar(path):
    return lp.ObservationWrapper(path, extraction=EXTRACTION)


def targname(path):
    raw = os.path.join(lp.DATA_ROOT, 'raw', '1998aug26_lsf', koaid(path) + '.fits')
    return str(fits.getheader(raw).get('TARGNAME', '?')).strip().upper()


def frame_snr(obs, positions):
    s = []
    for i in positions:
        w = obs._weight[i] > 0
        if w.any():
            s.append(np.median(obs._flux[i][w] * np.sqrt(obs._weight[i][w])))
    return float(np.median(s)) if s else 0.


def koaid(path):
    return os.path.basename(path)[7:24]


def build_chunks(obs, Pars, positions):
    """ The template's grid: `auto_equal_width` with the template parameters. """
    return pyodine.chunks.auto_equal_width(
        obs, width=Pars.chunk_width, padding=Pars.chunk_padding, orders=positions,
        chunks_per_order=Pars.chunks_per_order, pix_offset0=Pars.pix_offset0)


class OstarTrial(hp.Template_Parameters):
    """ ``Template_Parameters`` with run 1 set to one of prompt 6's LSF models.

    Velocity and template depth stay fixed (a B star has no template);
    ``iod_depth`` stays free per chunk, which is the point.
    """

    def __init__(self, name=None):
        super().__init__()
        if name is not None:
            model, setup, _ = cl.MODELS[name]
            self.model_runs[1] = hp._run_dict(model, lsf_setup_dict=setup)
        self.trial = name or 'super'

    def constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        if run_id == 0:
            return super().constrain_parameters(lmfit_params, 0, run_results, fitter)
        median_lsf = run_results[0]['median_pars'].filter('lsf')
        lsf_start = fitter.fit_lsfs(self.model_runs[0]['lsf_model'], median_lsf)
        model_name = self.model_runs[1]['lsf_model'].name()
        active = getattr(fitter.model.lsf_model, 'pars_dict', {})
        for i in range(len(lmfit_params)):
            lmfit_params[i]['velocity'].set(vary=False)
            lmfit_params[i]['tem_depth'].set(vary=False)
            lmfit_params[i]['iod_depth'].set(value=run_results[0]['median_pars']['iod_depth'],
                                             min=0.1)
            lmfit_params[i]['wave_slope'].set(value=run_results[0]['wave_slope_fit'][i])
            lmfit_params[i]['wave_intercept'].set(value=run_results[0]['wave_intercept_fit'][i])
            lmfit_params[i]['cont_intercept'].set(
                value=run_results[0]['results'][i].params['cont_intercept'])
            lmfit_params[i]['cont_slope'].set(
                value=run_results[0]['results'][i].params['cont_slope'])
            for p, v in lsf_start.items():
                if model_name == 'HermiteGaussian':
                    if p != 'fwhm' and not active.get(p, 1):
                        lmfit_params[i]['lsf_' + p].set(value=0., vary=False)
                        continue
                    lo, hi = (0.5, 4.0) if p == 'fwhm' else (-0.5, 0.5)
                else:
                    lo, hi = hp.lsf_bounds(model_name, p, v)
                lmfit_params[i]['lsf_' + p].set(value=hp.start_inside(v, lo, hi), min=lo, max=hi)
        return lmfit_params


def fit_ostar(obs, Pars, iod, positions):
    """ Every run of ``Pars`` on one B-star observation, as the drivers do it.

    Returns:
        tuple: (chunks, {run_id: list of fit results}, runtime)
    """
    t0 = time.time()
    chunks = build_chunks(obs, Pars, positions)
    weight = obs.compute_weight(weight_type=Pars.weight_type)
    n_bad = 0
    if Pars.bad_pixel_mask:
        # `pyodine`'s own cosmic / bad-pixel finder, as ``Template_Parameters``
        # asks.  `create_template` computes it on the *template* observation
        # and applies it to the hot star's weights; here it is computed on
        # the B-star frame itself, whose own defects are what matter.  The
        # boxcar extraction does not mask bad columns the optimal one did.
        mask = pyodine.bad_pixels.BadPixelMask(obs, cutoff=Pars.bad_pixel_cutoff)
        bad = (mask.mask == 1.) & (weight > 0)
        n_bad = int(bad.sum())
        weight[bad] = 0.
    chunk_weight = [np.array(weight[c.order, c.abspix[0]:c.abspix[-1] + 1]) for c in chunks]
    run_results, out = {}, {}
    for run_id, run in Pars.model_runs.items():
        lsf_model = run['lsf_model']
        if isinstance(run.get('lsf_setup_dict'), dict):
            lsf_model.adapt_LSF(run['lsf_setup_dict'])
        model = pyodine.models.spectrum.SimpleModel(
            lsf_model, run['wave_model'], run['cont_model'], iod, stellar_template=None,
            osample_factor=Pars.osample_obs, conv_width=Pars.lsf_conv_width)
        fitter = LmfitWrapper(model)
        starting = [model.guess_params(c) for c in chunks]
        for name, key in (('wave_slope', 'pre_wave_slope_deg'),
                          ('wave_intercept', 'pre_wave_intercept_deg')):
            if run.get(key, 0) > 0:
                poly = pyodine.lib.misc.smooth_parameters_over_orders(starting, name, chunks,
                                                                      deg=run[key])
                for i in range(len(chunks)):
                    starting[i][name] = poly[i]
        lmfit_params = [fitter.convert_params(p, to_lmfit=True) for p in starting]
        run_results[run_id] = {'velocity_guess': 0.}
        lmfit_params = Pars.constrain_parameters(lmfit_params, run_id, run_results, fitter)
        results, chunk_w, failed, chauv, rchi = pipe_lib.model_all_chunks(
            chunks, chunk_weight, fitter, lmfit_params, None,
            use_chauvenet=run.get('use_chauvenet_pixels', False), compute_redchi2=True,
            use_progressbar=False)
        rr = run_results[run_id]
        rr['results'] = results
        names = list(results[0].params.keys())
        rr['median_pars'] = ParameterSet({n: np.nanmedian([r.params[n] for r in results])
                                          for n in names})
        for par in ('wave_slope', 'wave_intercept'):
            rr[par + '_fit'] = pyodine.lib.misc.smooth_fitresult_over_orders(
                results, par, deg=run.get(par + '_deg', 3))
        out[run_id] = results
    return chunks, out, time.time() - t0


def table_ostar(chunks, results, obs, iod, run_id, Pars, k):
    """ One row per chunk of a B-star fit. """
    lsf_model = Pars.model_runs[run_id]['lsf_model']
    x_lsf = lsf_model.generate_x(20, Pars.lsf_conv_width)
    rows = []
    for i, (c, r) in enumerate(zip(chunks, results)):
        lr = r.lmfit_result
        row = dict(chunk=i, order=int(obs.ech_orders[c.order]), pixc=int(c.abspix[len(c) // 2]),
                   wave=float(obs._wave[c.order, c.abspix[len(c) // 2]]),
                   ok=lr is not None and r.errors['iod_depth'] is not None)
        if lr is None:
            rows.append(row)
            continue
        p, e = r.params, r.errors
        row.update(iod_depth=float(p['iod_depth']),
                   e_iod_depth=float(e['iod_depth']) if e['iod_depth'] is not None else np.nan,
                   redchi=float(lr.redchi), redchi_meas=float(lr.redchi) / k ** 2,
                   no_unc='uncertainties could not be estimated' in r.report)
        row['alpha'] = lp.ATLAS_ALPHA * row['iod_depth']
        n = len(c)
        pix = np.arange(-(n // 2), n - n // 2)
        w_pyp = obs._wave[c.order, c.abspix[0]:c.abspix[-1] + 1]
        row['dv_wave_c'] = float((p['wave_intercept'] / w_pyp[n // 2] - 1.) * CKMS * 1e3)
        half = n / 2. + Pars.chunk_padding
        sel = (iod.wave >= p['wave_intercept'] - p['wave_slope'] * half) & \
              (iod.wave <= p['wave_intercept'] + p['wave_slope'] * half)
        row['iod_contrast'] = float(np.std(iod.flux[sel]) / np.mean(iod.flux[sel])) \
            if sel.any() else np.nan
        try:
            y = lsf_model.eval(x_lsf, p.filter('lsf'))
            row['lsf_w'], row['lsf_cen'], row['lsf_asym'] = fe.lsf_shape(x_lsf, y)
        except Exception:                                    # noqa: BLE001
            row['lsf_w'] = row['lsf_cen'] = row['lsf_asym'] = np.nan
        rows.append(row)
    t = Table(rows)
    for c in ('iod_depth', 'e_iod_depth', 'redchi', 'redchi_meas', 'alpha', 'dv_wave_c',
              'iod_contrast', 'lsf_w', 'lsf_cen', 'lsf_asym'):
        if c not in t.colnames:
            t[c] = np.nan
        t[c] = np.where(t['ok'], t[c], np.nan)
    return t


def pair_noise(file_a, file_b, positions, smooth=101):
    """ Measured / propagated noise from two consecutive exposures, pixel by
    pixel (the method of `deconvolve_template.noise_check`). """
    from scipy.ndimage import median_filter
    ea, eb = open_bstar(file_a), open_bstar(file_b)
    widths = []
    for i in positions:
        f1_, iv1 = ea._flux[i], ea._weight[i]
        f2_, iv2 = eb._flux[i], eb._weight[i]
        ok = (iv1 > 0) & (iv2 > 0)
        if ok.sum() < 500:
            continue
        ratio = np.where(ok, f1_ / np.where(f2_ != 0, f2_, 1.), np.nan)
        scale = median_filter(np.nan_to_num(ratio, nan=np.nanmedian(ratio)), smooth)
        z = ((f1_ - scale * f2_) / np.sqrt(1. / np.where(ok, iv1, np.inf)
                                            + scale ** 2 / np.where(ok, iv2, np.inf)))[ok]
        widths.append(1.4826 * np.median(np.abs(z - np.median(z))))
    return float(np.median(widths)), len(widths)


# ---------------------------------------------------------------------------
# The payoff: the epoch with the B-star alpha
# ---------------------------------------------------------------------------

class AlphaFixedParameters(hp.Parameters):
    """ The settled observation parameters, with run 1's ``iod_depth`` fixed
    per chunk to the B stars' value (by template chunk index). """

    def __init__(self, iod_by_template_chunk):
        super().__init__()
        self.iod = iod_by_template_chunk
        self.vel_analysis_plots = None
        for r in self.model_runs.values():
            r['plot_analysis'] = False
            r['plot_success'] = False

    def constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        lmfit_params = super().constrain_parameters(lmfit_params, run_id, run_results, fitter)
        if run_id == 1:
            n_fixed = 0
            for i in range(len(lmfit_params)):
                ti = getattr(run_results[0]['results'][i].chunk, 'template_index', i)
                v = self.iod.get(ti, np.nan)
                if np.isfinite(v):
                    lmfit_params[i]['iod_depth'].set(value=float(v), vary=False)
                    n_fixed += 1
            logging.info('iod_depth fixed from the B stars in {} chunks'.format(n_fixed))
        return lmfit_params


class AlphaFreeParameters(AlphaFixedParameters):
    """ The same configuration with ``iod_depth`` left free: the control. """

    def constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        return hp.Parameters.constrain_parameters(self, lmfit_params, run_id, run_results, fitter)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(bt, alpha_map, lsfcmp, payoff, outfile):
    import matplotlib.pyplot as plt
    INK, C1, C2, C3 = '#0b0b0b', '#2a78d6', '#eb6834', '#1baf7a'
    fig, axs = plt.subplots(2, 3, figsize=(17, 9.5))

    ax = axs[0, 0]
    for f in np.unique(bt['frame']):
        s = bt[(bt['frame'] == f) & np.isfinite(bt['alpha'])]
        ax.plot(s['wave'], s['alpha'], '.', ms=2, alpha=0.3, color='0.5')
    ok = np.isfinite(alpha_map['alpha'])
    ax.plot(alpha_map['wave'][ok], alpha_map['alpha'][ok], '.', ms=4, color=C1,
            label='median over the B-star frames')
    ax.axhline(lp.ATLAS_ALPHA, color=C2, ls='--', label='adapter alpha = {:.2f}'.format(lp.ATLAS_ALPHA))
    ax.set_ylim(0, 6)
    ax.set_xlabel('wavelength (A)')
    ax.set_ylabel('effective alpha (2.59 x iod_depth)')
    ax.legend(fontsize=8)
    ax.set_title('(a) the per-chunk Beer-Lambert exponent, B stars', loc='left', fontsize=10)

    ax = axs[0, 1]
    g = bt[(bt['run'] == 1) & np.isfinite(bt['redchi_meas'])]
    ax.plot(g['iod_contrast'], g['redchi_meas'], '.', ms=3, alpha=0.4, color=C1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('iodine contrast of the chunk (RMS / mean, scaled atlas)')
    ax.set_ylabel('chi2/dof, measured noise')
    ax.set_title('(b) B-star misfit against iodine structure (run 1)', loc='left', fontsize=10)

    ax = axs[0, 2]
    ok = np.isfinite(alpha_map['lsf_w'])
    sc = ax.scatter(alpha_map['pixc'][ok], alpha_map['lsf_w'][ok], c=alpha_map['order'][ok],
                    s=8, cmap='viridis')
    fig.colorbar(sc, ax=ax, label='order')
    ax.set_xlabel('detector pixel')
    ax.set_ylabel('LSF FWHM (px), B stars, run 1')
    ax.set_title('(c) the LSF from the B stars', loc='left', fontsize=10)

    ax = axs[1, 0]
    names = list(dict.fromkeys(lsfcmp['model']))
    ax.plot(range(len(names)), [np.median(lsfcmp['redchi_meas'][lsfcmp['model'] == n])
                                for n in names], 'o', color=C1)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names)
    ax.set_yscale('log')
    ax.set_ylabel('median chi2/dof (measured), HR 8634 pair')
    ax.set_title('(d) LSF models on pure iodine', loc='left', fontsize=10)

    ax = axs[1, 1]
    for lab, col in (('iod_depth free', C2), ('iod_depth fixed (B stars)', C1)):
        t = payoff[lab]
        g = t[t['good']]
        ax.plot(g['wave'], g['velocity'], '.', ms=3, alpha=0.5, color=col,
                label='{:s}: robust sigma {:.0f} m/s'.format(lab, robust_sigma(g['velocity'])))
    ax.set_xlabel('wavelength (A)')
    ax.set_ylabel('chunk velocity, run 1 (m/s)')
    ax.legend(fontsize=8)
    ax.set_title('(e) the epoch, with and without the B-star alpha', loc='left', fontsize=10)

    ax = axs[1, 2]
    for lab, col in (('iod_depth free', C2), ('iod_depth fixed (B stars)', C1)):
        t = payoff[lab]
        g = t[t['good']]
        ax.plot(g['iod_contrast'], g['redchi_meas'], '.', ms=3, alpha=0.4, color=col, label=lab)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('iodine contrast of the chunk')
    ax.set_ylabel('chi2/dof, measured noise')
    ax.legend(fontsize=8)
    ax.set_title('(f) does the iodine correlation go away?', loc='left', fontsize=10)

    fig.tight_layout()
    fig.savefig(outfile, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def parser(options=None):
    p = argparse.ArgumentParser(description='The iodine cell seen through B stars',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--datadir', type=str, default=DATA_DIR)
    p.add_argument('--figdir', type=str, default=FIG_DIR)
    return p.parse_args() if options is None else p.parse_args(options)


def main(pargs):
    files = bstar_files()
    if not files:
        raise FileNotFoundError('No B-star spec1d in {:s}; run reduce_bstars'.format(BSTAR_SCIENCE))
    Pars = OstarTrial()
    iod = lp.IodineTemplate(Pars.i2_to_use)
    positions = list(range(Pars.temp_order_range[0], Pars.temp_order_range[1] + 1))
    usable = []
    for f in files:
        snr = frame_snr(open_bstar(f), positions)
        if snr < MIN_SNR:
            print('Skipping {:s} ({:s}): median S/N {:.0f} in the iodine orders'.format(
                koaid(f), targname(f), snr))
        else:
            usable.append(f)
    files = usable

    # 2. the noise of these frames
    pair = [f for f in files if koaid(f) in PAIR]
    k, n_ord = pair_noise(pair[0], pair[1], positions)
    print('Noise, {:s} - {:s}, pixel by pixel: measured / propagated = {:.3f} over {:d} '
          'orders (HD 187123 template exposures: {:.3f})'.format(
              koaid(pair[0]), koaid(pair[1]), k, n_ord, hp.NOISE_SCALE))

    # 1. every frame, runs 0 and 1 (super-Gaussian)
    tabs = []
    for f in files:
        obs = open_bstar(f)
        star = targname(f)
        chunks, res, rt = fit_ostar(obs, Pars, iod, positions)
        for run_id in (0, 1):
            t = table_ostar(chunks, res[run_id], obs, iod, run_id, Pars, k)
            t['frame'], t['star'], t['run'] = koaid(f), star, run_id
            tabs.append(t)
        t1 = tabs[-1]
        print('  {:s} {:8s} {:4.0f} s fit: alpha median {:.2f} (robust {:.2f}); chi2 (meas) '
              'median {:.2f}; LSF {:.3f} px; no-unc {:d}; wave vs PypeIt {:+.0f} m/s'.format(
                  koaid(f), star, rt, np.nanmedian(t1['alpha']), robust_sigma(t1['alpha']),
                  np.nanmedian(t1['redchi_meas']), np.nanmedian(t1['lsf_w']),
                  int(np.nansum(t1['no_unc'].astype(float))), np.nanmedian(t1['dv_wave_c'])),
              flush=True)
    bt = vstack(tabs, metadata_conflicts='silent')

    # The alpha map: median over frames, per template chunk, run 1
    r1 = bt[bt['run'] == 1]
    rows = []
    for i in np.unique(r1['chunk']):
        s = r1[r1['chunk'] == i]
        a = np.asarray(s['alpha'], dtype=float)
        a = a[np.isfinite(a)]
        rows.append(dict(chunk=int(i), order=int(s['order'][0]), pixc=int(s['pixc'][0]),
                         wave=float(s['wave'][0]), n=len(a),
                         alpha=float(np.median(a)) if len(a) else np.nan,
                         alpha_scatter=float(robust_sigma(a)) if len(a) > 2 else np.nan,
                         iod_depth=float(np.median(a)) / lp.ATLAS_ALPHA if len(a) else np.nan,
                         e_iod_depth=float(np.nanmedian(s['e_iod_depth'])),
                         redchi_meas=float(np.nanmedian(s['redchi_meas'])),
                         iod_contrast=float(np.nanmedian(s['iod_contrast'])),
                         lsf_w=float(np.nanmedian(s['lsf_w'])),
                         lsf_asym=float(np.nanmedian(s['lsf_asym'])),
                         dv_wave_c=float(np.nanmedian(s['dv_wave_c']))))
    amap = Table(rows)

    from scipy.stats import spearmanr
    ok = np.isfinite(amap['alpha']) & np.isfinite(amap['redchi_meas'])
    print('Alpha map ({:d} chunks): median {:.3f}, robust spread between chunks {:.3f}; '
          'median frame-to-frame scatter per chunk {:.3f}; median chunk error {:.3f}'.format(
              int(ok.sum()), np.median(amap['alpha'][ok]), robust_sigma(amap['alpha'][ok]),
              np.nanmedian(amap['alpha_scatter']),
              np.nanmedian(amap['e_iod_depth']) * lp.ATLAS_ALPHA))
    print('  B-star chi2 (measured noise): median {:.2f}; Spearman with iodine contrast '
          '{:+.2f}; alpha vs wavelength {:+.2f}; alpha vs iodine contrast {:+.2f}'.format(
              np.median(amap['redchi_meas'][ok]),
              spearmanr(amap['iod_contrast'][ok], amap['redchi_meas'][ok])[0],
              spearmanr(amap['wave'][ok], amap['alpha'][ok])[0],
              spearmanr(amap['iod_contrast'][ok], amap['alpha'][ok])[0]))
    print('  LSF from B stars: FWHM median {:.3f} px; asymmetry {:+.3f}; wavelength vs PypeIt '
          '{:+.0f} m/s (order-to-order {:.0f})'.format(
              np.nanmedian(amap['lsf_w']), np.nanmedian(amap['lsf_asym']),
              np.nanmedian(amap['dv_wave_c']),
              robust_sigma([np.nanmedian(amap['dv_wave_c'][amap['order'] == o])
                            for o in np.unique(amap['order'])])))

    # 3. the five LSF models on the HR 8634 pair
    lrows = []
    for name in cl.MODELS:
        TP = OstarTrial(name)
        for f in pair:
            obs = open_bstar(f)
            chunks, res, rt = fit_ostar(obs, TP, iod, positions)
            t = table_ostar(chunks, res[1], obs, iod, 1, TP, k)
            nlsf = cl.MODELS[name][2]
            g = t[np.isfinite(t['redchi_meas'])]
            chi2 = np.sum(g['redchi_meas'] * (40 - (5 + nlsf)))
            lrows.append(dict(model=name, frame=koaid(f), runtime=rt,
                              redchi_meas=float(np.median(g['redchi_meas'])),
                              no_unc=int(np.nansum(t['no_unc'].astype(float))),
                              alpha=float(np.nanmedian(t['alpha'])),
                              alpha_robust=float(robust_sigma(t['alpha'])),
                              lsf_w=float(np.nanmedian(t['lsf_w'])),
                              lsf_asym=float(np.nanmedian(t['lsf_asym'])),
                              bic=float(chi2 + len(g) * (5 + nlsf) * np.log(40))))
            print('  LSF {:10s} {:s}: chi2 {:.2f}  alpha {:.3f} +- {:.3f}  FWHM {:.3f}  asym '
                  '{:+.3f}  no-unc {:d}  BIC {:.4g}'.format(
                      name, koaid(f), lrows[-1]['redchi_meas'], lrows[-1]['alpha'],
                      lrows[-1]['alpha_robust'], lrows[-1]['lsf_w'], lrows[-1]['lsf_asym'],
                      lrows[-1]['no_unc'], lrows[-1]['bic']), flush=True)
    lsfcmp = Table(lrows)

    # 4. the payoff: prompt 5's epoch, super-Gaussian, os 1, iod_depth fixed vs free
    iod_map = {int(r['chunk']): r['iod_depth'] for r in amap}
    obs_file = f1.epoch_file(fe.EPOCH)
    obs = lp.ObservationWrapper(obs_file)
    template = StellarTemplate_Chunked(dt.TEMPLATE_H5[1])
    exp = f1.expected_velocity(dict(obs=obs, template=template))
    weight = obs.compute_weight(hp.Parameters().weight_type)
    flat = fe.blaze(obs_file, obs)
    payoff, prow = {}, []
    for lab, P in (('iod_depth free', AlphaFreeParameters(iod_map)),
                   ('iod_depth fixed (B stars)', AlphaFixedParameters(iod_map))):
        outdir = os.path.join(fe.RUN_ROOT, 'bstar_alpha_{:s}_{:s}'.format(
            'fixed' if 'fixed' in lab else 'free', fe.EPOCH))
        info = fe.run_driver(obs_file, dt.TEMPLATE_H5[1], outdir, P)
        t = fe.chunk_table(os.path.join(outdir, 'run1.h5'), 1, P, obs, template, iod, flat,
                           weight, P.noise_scale)
        t = fe.flag_failures(t, P)
        s = fe.summarise_run(t, exp['expected'])
        md = fe.misfit_drivers(t)
        payoff[lab] = t
        prow.append(dict(config=lab, n_good=s['n_good'], v_robust=s['v_robust'],
                         order_scatter=s['order_scatter'], within_order=s['within_order'],
                         redchi_meas=s['redchi_meas'], v_median=s['v_median'],
                         v_expected=s['v_expected'], epoch_error=s['epoch_error'],
                         partial_iod=md['partial_iod'], partial_tem=md['partial_tem'],
                         no_unc=s['no_unc'], runtime=info['total']))
        print('  {:26s} good {:3d}  robust {:6.1f}  order-to-order {:6.1f}  chi2 {:5.1f}  '
              'v {:+7.1f} (exp {:+.1f} +- {:.1f})  partial iod {:+.2f} tem {:+.2f}'.format(
                  lab, s['n_good'], s['v_robust'], s['order_scatter'], s['redchi_meas'],
                  s['v_median'], s['v_expected'], s['epoch_error'], md['partial_iod'],
                  md['partial_tem']), flush=True)
    b = cl.paired_bootstrap(payoff['iod_depth free'], payoff['iod_depth fixed (B stars)'])
    print('  paired, fixed against free: dsigma {:+.1f} [{:+.1f}, {:+.1f}] m/s on {:d} chunks; '
          'dv median {:+.1f}'.format(b['dsig'], b['lo'], b['hi'], b['n'], b['dv_median']))
    ptab = Table(prow)
    ptab.meta.update({'dsig_' + kk: vv for kk, vv in b.items()})

    for tab in (bt, amap, lsfcmp, ptab):
        for c in tab.colnames:
            if tab[c].dtype.kind == 'f':
                tab[c].info.format = '.10g' if c == 'wave' else '.6g'
    bt.write(os.path.join(pargs.datadir, 'bstar_chunks.csv'), overwrite=True)
    amap.write(os.path.join(pargs.datadir, 'bstar_alpha_map.csv'), overwrite=True)
    lsfcmp.write(os.path.join(pargs.datadir, 'bstar_lsf_models.csv'), overwrite=True)
    ptab.write(os.path.join(pargs.datadir, 'bstar_alpha_payoff.csv'), overwrite=True)
    make_figure(bt, amap, lsfcmp, payoff, os.path.join(pargs.figdir, 'fig_p6b_bstars.png'))


if __name__ == '__main__':
    main(parser())
