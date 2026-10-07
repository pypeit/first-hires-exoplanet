""" Settle the instrumental profile (phase 3, prompt 6).

Every LSF model `pyodine` offers, fitted to the prompt-5 epoch through the
driver prompt 5 used, and compared on the numbers that decide the final
precision.

THE MODELS.  All share run 0 -- a single Gaussian, the first wavelength
solution -- and differ in run 1, which starts from run 0's medians exactly as
Lick's scheme does:

  single      SingleGaussian: FWHM only (1 parameter)
  super       SuperGaussian: width, exponent and two satellites (4)
  multi_lick  MultiGaussian_Lick: Lick's 11-Gaussian layout, 10 amplitudes,
              not re-centred (10)
  multi_song  MultiGaussian: SONG's layout, 10 amplitudes (10)
  hermite     HermiteGaussian: FWHM and the Hermite weights SONG enables,
              3-8 (7)

and each is carried into a run 2 with `FixedLSF`: run 1's LSFs smoothed over
neighbouring chunks (+/-160 pixels, +/-3 orders; Lick's radii) and then held
fixed -- the Butler (1996) / Lick practice that trades per-chunk freedom for
per-chunk noise.  The order separation for that smoothing is HIRES's, 71
physical pixels (36 binned; measured from the spec1d traces), not Lick's 15.

Two departures from the upstream configurations, both because the upstream
setting would not test the model it names:

* SONG's Hermite constraints bound every weight to within 2e-13 of a start of
  ~1e-9: the "Hermite" model is then a single Gaussian.  Here the enabled
  weights may range over [-0.5, 0.5].
* The super-Gaussian's exponent and width are kept physical (exponent 1-4,
  sigma 0.2-3 px); Lick's generic p +/- (|p| + 0.3) allows a negative one.

Both templates of prompt 3 are carried, as the Q&A after prompt 3 decided.
Weights are PypeIt's; every chi-square is also given against the measured
noise (Q&A after prompt 3, question 4).

WHAT DECIDES.  The per-chunk velocity scatter of the good chunks, robust, and
the scatter between order medians -- the two that feed an epoch velocity.
Differences are tested by a paired bootstrap over the chunks both models
call good.  Failures and chi-square are secondary; a single epoch's median
against the prediction (+/- ~13 m/s) is reported but cannot decide.

Run with:

    conda run --no-capture-output -n pypeit14 \\
        python -m first_hires_exoplanet.compare_lsf_models

"""

# Must precede any import of matplotlib, including pyodine's own
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# Standard imports
import time
import argparse
import logging

import numpy as np
from astropy.table import Table, vstack

from first_hires_exoplanet.utilities_hires import load_pyodine as lp
from first_hires_exoplanet.utilities_hires import pyodine_parameters as hires_pars
from first_hires_exoplanet import deconvolve_template as dt
from first_hires_exoplanet import fit_one_chunk as f1
from first_hires_exoplanet import fit_one_epoch as fe
from first_hires_exoplanet.pyodine_chunk_stats import robust_sigma

from pyodine import models
from pyodine.template.base import StellarTemplate_Chunked


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_HERE, 'data')
FIG_DIR = os.path.join(os.path.dirname(_HERE), 'docs', 'figs')

EPOCH = fe.EPOCH

#: Physical pixels between adjacent iodine orders (36 binned, binning 2)
ORDER_SEPARATION = 71

_HERMITE_SETUP = {'weight_1': 0, 'weight_2': 0, 'weight_3': 1, 'weight_4': 1,
                  'weight_5': 1, 'weight_6': 1, 'weight_7': 1, 'weight_8': 1,
                  'weight_9': 0}
_MULTI_SONG_SETUP = {
    'positions': [-2.9, -2.5, -1.9, -1.4, -1.0, 0.0, 1.0, 1.4, 1.9, 2.5, 2.9],
    'sigmas':    [0.9, 0.9, 0.9, 0.9, 0.9, 0.6, 0.9, 0.9, 0.9, 0.9, 0.9]}

#: name -> (LSF model, setup dict, number of LSF parameters)
MODELS = {
    'single':     (models.lsf.SingleGaussian, None, 1),
    'super':      (models.lsf.SuperGaussian, None, 4),
    'multi_lick': (models.lsf.MultiGaussian_Lick, hires_pars._multigauss_setup_dict, 10),
    'multi_song': (models.lsf.MultiGaussian, _MULTI_SONG_SETUP, 10),
    'hermite':    (models.lsf.HermiteGaussian, _HERMITE_SETUP, 7),
}


# ---------------------------------------------------------------------------
# Parameters for one trial
# ---------------------------------------------------------------------------

class TrialParameters(hires_pars.Parameters):
    """ ``utilities_hires`` parameters with run 1 set to one LSF model, and a
    run 2 that smooths it and holds it fixed. """

    def __init__(self, name):
        super().__init__()
        # The trial configurations are the prompt-6 comparison, independent of
        # what utilities_hires has since settled on
        model, setup, _ = MODELS[name]
        self.trial = name
        self.model_runs = {
            0: hires_pars._run_dict(models.lsf.SingleGaussian, pre_wave_deg=3),
            1: hires_pars._run_dict(model, lsf_setup_dict=setup),
            2: hires_pars._run_dict(models.lsf.FixedLSF),
        }
        self.model_runs[2].update({'smooth_lsf_run': 1, 'smooth_pixels': 160,
                                   'smooth_orders': 3,
                                   'order_separation': ORDER_SEPARATION,
                                   'smooth_manual_redchi': False, 'smooth_osample': 0})
        for r in self.model_runs.values():
            r['plot_analysis'] = False          # the comparison makes its own
            r['plot_success'] = False
        self.vel_analysis_plots = None

    def _lsf_bounds(self, name, value):
        """ Bounds for one LSF parameter of this trial's run-1 model. """
        model = self.model_runs[1]['lsf_model'].name()
        if model == 'HermiteGaussian':
            return (0.5, 4.0) if name == 'fwhm' else (-0.5, 0.5)
        return hires_pars.lsf_bounds(model, name, value)

    #: See `pyodine_parameters.start_inside`: found here, kept there
    _start = staticmethod(hires_pars.start_inside)

    def constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        if run_id == 0:
            return super().constrain_parameters(lmfit_params, 0, run_results, fitter)
        prev = run_id - 1
        logging.info('Constraining parameters for RUN {} ({})'.format(run_id, self.trial))
        if run_id == 1:
            median_lsf = run_results[0]['median_pars'].filter('lsf')
            lsf_start = fitter.fit_lsfs(self.model_runs[0]['lsf_model'], median_lsf)
            active = getattr(fitter.model.lsf_model, 'pars_dict', {})
        for i in range(len(lmfit_params)):
            lmfit_params[i]['velocity'].set(value=run_results[prev]['median_pars']['velocity'])
            lmfit_params[i]['iod_depth'].set(value=run_results[prev]['median_pars']['iod_depth'])
            lmfit_params[i]['tem_depth'].set(value=run_results[prev]['median_pars']['tem_depth'])
            lmfit_params[i]['wave_slope'].set(value=run_results[prev]['wave_slope_fit'][i])
            lmfit_params[i]['wave_intercept'].set(value=run_results[prev]['wave_intercept_fit'][i])
            lmfit_params[i]['cont_intercept'].set(
                value=run_results[prev]['results'][i].params['cont_intercept'])
            lmfit_params[i]['cont_slope'].set(
                value=run_results[prev]['results'][i].params['cont_slope'])
            if run_id == 1:
                for p, v in lsf_start.items():
                    if self.model_runs[1]['lsf_model'].name() == 'HermiteGaussian' \
                            and p != 'fwhm' and not active.get(p, 1):
                        lmfit_params[i]['lsf_' + p].set(value=0., vary=False)
                        continue
                    lo, hi = self._lsf_bounds(p, v)
                    lmfit_params[i]['lsf_' + p].set(value=self._start(v, lo, hi), min=lo, max=hi)
            else:
                for p in ('order', 'pixel0', 'amplitude'):
                    lmfit_params[i]['lsf_' + p].vary = False
        return lmfit_params


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def paired_bootstrap(ta, tb, nboot=2000, seed=6):
    """ Robust-sigma difference (b - a) on the chunks both call good, with a
    68% interval from resampling those chunks. """
    ga = {int(c): v for c, v, g in zip(ta['chunk'], ta['velocity'], ta['good']) if g}
    gb = {int(c): v for c, v, g in zip(tb['chunk'], tb['velocity'], tb['good']) if g}
    common = sorted(set(ga) & set(gb))
    va = np.array([ga[c] for c in common])
    vb = np.array([gb[c] for c in common])
    rng = np.random.default_rng(seed)
    d = []
    for _ in range(nboot):
        k = rng.integers(0, len(common), len(common))
        d.append(robust_sigma(vb[k]) - robust_sigma(va[k]))
    dv = vb - va
    return dict(n=len(common), dsig=float(robust_sigma(vb) - robust_sigma(va)),
                lo=float(np.percentile(d, 16)), hi=float(np.percentile(d, 84)),
                dv_median=float(np.median(dv)), dv_robust=float(robust_sigma(dv)))


def bic(t, n_lsf, k):
    """ Summed chi-square (measured noise) plus k ln n, over good chunks.

    With chi2/dof ~ 20 the misfit is systematic, so this rewards flexibility
    more than it measures truth; reported, not used to decide.
    """
    g = t[t['good']]
    npix = 40
    chi2 = np.sum(g['redchi_meas'] * (npix - (7 + n_lsf)))
    return float(chi2 + len(g) * (7 + n_lsf) * np.log(npix))


def make_figure(tabs, sums, outfile):
    """ Scatter and order-to-order scatter by model, run and template. """
    import matplotlib.pyplot as plt
    INK, C1, C2 = '#0b0b0b', '#2a78d6', '#eb6834'
    names = list(MODELS)
    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    for ax, key, label in ((axs[0], 'v_robust', 'per-chunk robust sigma (m/s)'),
                           (axs[1], 'order_scatter', 'scatter of order medians (m/s)'),
                           (axs[2], 'redchi_meas', 'median chi2/dof, measured noise')):
        for os_, col, dx in ((10, C2, -0.12), (1, C1, 0.12)):
            for run, mk in ((1, 'o'), (2, 's')):
                y = [sums[(n, os_)][run][key] for n in names]
                ax.plot(np.arange(len(names)) + dx + (0.05 if run == 2 else 0), y, mk,
                        color=col, mfc=col if run == 1 else 'none', ms=7,
                        label='os {:d}, run {:d}{}'.format(os_, run,
                                                          ' (fixed, smoothed)' if run == 2 else ''))
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names)
        ax.set_ylabel(label)
    axs[0].legend(fontsize=8)
    axs[0].set_title('(a) what the chunks do', loc='left', fontsize=10)
    axs[1].set_title('(b) what an epoch inherits', loc='left', fontsize=10)
    axs[2].set_title('(c) how well the model fits', loc='left', fontsize=10)
    fig.tight_layout()
    fig.savefig(outfile, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def parser(options=None):
    p = argparse.ArgumentParser(description='Compare LSF models (phase 3, prompt 6)',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--models', nargs='+', default=list(MODELS))
    p.add_argument('--skip-run', action='store_true')
    p.add_argument('--datadir', type=str, default=DATA_DIR)
    p.add_argument('--figdir', type=str, default=FIG_DIR)
    return p.parse_args() if options is None else p.parse_args(options)


def main(pargs):
    obs_file = f1.epoch_file(EPOCH)
    obs = lp.ObservationWrapper(obs_file)
    base = hires_pars.Parameters()
    k = base.noise_scale
    weight = obs.compute_weight(base.weight_type)
    flat = fe.blaze(obs_file, obs)
    iod = lp.IodineTemplate(base.i2_to_use)
    arc = dt.arc_line_widths(os.path.join(lp.DEFAULT_REDUX, 'reduce_' + EPOCH[3:11],
                                          'Calibrations', 'WaveCalib_A_0_DET01.fits'))

    tabs, sums, rows = {}, {}, []
    for os_ in (10, 1):
        template = StellarTemplate_Chunked(dt.TEMPLATE_H5[os_])
        exp = f1.expected_velocity(dict(obs=obs, template=template))
        for name in pargs.models:
            Pars = TrialParameters(name)
            outdir = os.path.join(fe.RUN_ROOT, 'lsf_{:s}_{:s}_os{:d}'.format(name, EPOCH, os_))
            if not pargs.skip_run:
                info = fe.run_driver(obs_file, dt.TEMPLATE_H5[os_], outdir, Pars)
                runtime = info['total']
                t_run = info['t_run']
            else:
                runtime, t_run = np.nan, {r: np.nan for r in Pars.model_runs}
            sums[(name, os_)] = {}
            for run in (1, 2):
                res = os.path.join(outdir, 'run{:d}.h5'.format(run))
                t = fe.chunk_table(res, run, Pars, obs, template, iod, flat, weight, k)
                t = fe.flag_failures(t, Pars)
                t['model'], t['osample'] = name, os_
                s = fe.summarise_run(t, exp['expected'])
                s['misfit_drivers'] = fe.misfit_drivers(t)
                s['bic'] = bic(t, MODELS[name][2] if run == 1 else 0, k)
                tabs[(name, os_, run)] = t
                sums[(name, os_)][run] = s
                ls = fe.lsf_summary(t, arc) if run == 1 else None
                row = dict(model=name, osample=os_, run=run, n_lsf=MODELS[name][2] if run == 1 else 0,
                           runtime=runtime, runtime_run=t_run.get(run, np.nan),
                           n_good=s['n_good'], no_unc=s['no_unc'], misfit=s['misfit'],
                           v_outlier=s['v_outlier'], at_bound=s['at_bound'],
                           v_robust=s['v_robust'], v_std=s['v_std'], v_std_all=s['v_std_all'],
                           order_scatter=s['order_scatter'], within_order=s['within_order'],
                           v_median=s['v_median'], v_expected=s['v_expected'],
                           epoch_error=s['epoch_error'], sigv_noise=s['sigv_noise'],
                           sigv_stderr=s['sigv_stderr'], redchi_meas=s['redchi_meas'],
                           bic=s['bic'], partial_iod=s['misfit_drivers']['partial_iod'],
                           partial_tem=s['misfit_drivers']['partial_tem'],
                           lsf_fwhm=float(np.median(ls['fwhm'])) if ls is not None and len(ls) else np.nan,
                           lsf_asym=float(np.median(ls['asym'])) if ls is not None and len(ls) else np.nan,
                           lsf_cen=float(np.median(ls['centroid'])) if ls is not None and len(ls) else np.nan)
                rows.append(row)
                print('{:10s} os{:<2d} run {:d}: good {:3d}  robust {:6.1f}  order-to-order {:6.1f}  '
                      'within {:6.1f}  chi2 {:5.1f}  v {:+7.1f} (exp {:+.1f}, +-{:.1f})  '
                      'no-unc {:2d}  runtime {:5.1f}'.format(
                          name, os_, run, s['n_good'], s['v_robust'], s['order_scatter'],
                          s['within_order'], s['redchi_meas'], s['v_median'], s['v_expected'],
                          s['epoch_error'], s['no_unc'], runtime), flush=True)

    summary = Table(rows)
    # Paired comparisons against the prompt-5 configuration (Lick
    # multi-Gaussian, run 1) for the same template
    pairs = []
    for row in summary:
        ref = tabs[('multi_lick', row['osample'], 1)]
        cmp_ = tabs[(row['model'], row['osample'], row['run'])]
        b = paired_bootstrap(ref, cmp_)
        pairs.append(b)
    for key in ('n', 'dsig', 'lo', 'hi', 'dv_median', 'dv_robust'):
        summary['vs_lick_' + key] = [b[key] for b in pairs]
    # An analysis-only pass has no runtimes; keep those of the last full run
    old_file = os.path.join(pargs.datadir, 'lsf_models_summary.csv')
    if pargs.skip_run and os.path.exists(old_file):
        old = Table.read(old_file)
        look = {(r['model'], int(r['osample']), int(r['run'])): r for r in old}
        for col in ('runtime', 'runtime_run'):
            summary[col] = [look[(r['model'], int(r['osample']), int(r['run']))][col]
                            if (r['model'], int(r['osample']), int(r['run'])) in look
                            else np.nan for r in summary]
    for c in summary.colnames:
        if summary[c].dtype.kind == 'f':
            summary[c].info.format = '.4g'
    summary.write(os.path.join(pargs.datadir, 'lsf_models_summary.csv'), overwrite=True)
    allt = vstack([tabs[key] for key in sorted(tabs)], metadata_conflicts='silent')
    keep = ['model', 'osample', 'run', 'chunk', 'order', 'pixc', 'wave', 'velocity',
            'e_velocity', 'sigv_noise', 'redchi_meas', 'good', 'no_unc', 'misfit',
            'v_outlier', 'lsf_w', 'lsf_cen', 'lsf_asym', 'iod_depth', 'tem_depth',
            'dv_wave_c']
    allt = allt[keep]
    for c in allt.colnames:
        if allt[c].dtype.kind == 'f':
            allt[c].info.format = '.10g' if c == 'wave' else '.6g'
    allt.write(os.path.join(pargs.datadir, 'lsf_models_chunks.csv'), overwrite=True)

    print('\nPaired against multi_lick run 1 (same template): change in robust sigma, '
          '68% bootstrap interval; shift and scatter of chunk velocities')
    for row in summary:
        print('  {:10s} os{:<2d} run {:d}: n {:3d}  dsigma {:+6.1f} [{:+6.1f}, {:+6.1f}]  '
              'dv median {:+6.1f}  dv robust {:6.1f}'.format(
                  row['model'], row['osample'], row['run'], row['vs_lick_n'],
                  row['vs_lick_dsig'], row['vs_lick_lo'], row['vs_lick_hi'],
                  row['vs_lick_dv_median'], row['vs_lick_dv_robust']))
    # The template, paired: the same model and run, oversampling 1 against 10
    trows = []
    print('\nPaired, template oversampling 1 against 10 (same model, run 1): change in '
          'robust sigma, 68% interval; shift and scatter of chunk velocities')
    for name in pargs.models:
        b = paired_bootstrap(tabs[(name, 10, 1)], tabs[(name, 1, 1)])
        b['model'] = name
        trows.append(b)
        print('  {:10s} n {:3d}  dsigma {:+6.1f} [{:+6.1f}, {:+6.1f}]  dv median {:+6.1f}  '
              'dv robust {:6.1f}'.format(name, b['n'], b['dsig'], b['lo'], b['hi'],
                                         b['dv_median'], b['dv_robust']))
    Table(trows)['model', 'n', 'dsig', 'lo', 'hi', 'dv_median', 'dv_robust'].write(
        os.path.join(pargs.datadir, 'lsf_models_template_pairs.csv'), overwrite=True)
    make_figure(tabs, sums, os.path.join(pargs.figdir, 'fig_p6_lsf_models.png'))


if __name__ == '__main__':
    main(parser())
