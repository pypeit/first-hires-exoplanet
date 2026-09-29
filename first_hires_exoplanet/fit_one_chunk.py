""" Fit one chunk of one order of one epoch (phase 3, prompt 4).

The first forward-model fit of a HIRES observation.  It follows
`pyodine_model_observations.model_single_observation` step for step -- load
the epoch, the deconvolved template and the atlas through ``utilities_hires``;
the order correction, the weights, the velocity guess, `auto_wave_comoving`
chunking, the model's own starting guesses smoothed over the order, and
``Parameters.constrain_parameters`` for run 0 -- and then fits **one** chunk,
with the driver's Chauvenet re-fit.  Nothing is fitted to any other chunk.

WHICH CHUNK.  Epoch ``HI.19980826.19326``: cell in, 400 s, the night the
template was taken, so the relative velocity is small and predictable.  Order
63, the middle of the iodine band (orders 58-71), and the middle chunk of that
order.  Run 0's model: single-Gaussian LSF, linear wavelength, linear
continuum -- run 1 needs the median of run 0 over every chunk and so cannot be
a one-chunk fit.

WHICH TEMPLATE.  Both of prompt 3's, as the Q&A after prompt 3 decided: the
same chunk is fitted once against the Lick-default template (oversampling 10)
and once against the unsplined one (oversampling 1).

EVERY PARAMETER IS CHECKED AGAINST SOMETHING INDEPENDENT:

  velocity        the relative barycentric velocity (barycorrpy for the epoch,
                  astropy for the template) plus the planet's motion between
                  the two times, from a circular fit to the modern catalogue
  lsf_fwhm        the ThAr width of order 63 (prompt 3)
  iod_depth       1, if the adapter's alpha = 2.59 turns the atlas into the
                  HIRES cell
  tem_depth       1, if the template's line depths are right
  wave_*          PypeIt's ThAr wavelength solution for the same pixels
  cont_*          the chunk's own flux level

UNCERTAINTIES.  `lmfit` scales its covariance by the reduced chi-square
(``scale_covar=True``), so the 2.35x overstatement of PypeIt's inverse
variance (prompt 3) cancels out of every standard error; the chi-square is
reported both as fitted and divided by ``NOISE_SCALE**2``, as decided after
prompt 3.  The standard errors are then tested by Monte Carlo: the same chunk
is refitted to its own best-fit model plus noise at the *measured* level.

Run with:

    conda run --no-capture-output -n pypeit14 \\
        python -m first_hires_exoplanet.fit_one_chunk

"""

# Must precede any import of matplotlib, including pyodine's own
import os
os.environ.setdefault('MPLBACKEND', 'Agg')

# Standard imports
import copy
import time
import argparse

import numpy as np
from astropy.table import Table

from first_hires_exoplanet.utilities_hires import load_pyodine as lp
from first_hires_exoplanet.utilities_hires import pyodine_parameters as hires_pars
from first_hires_exoplanet import deconvolve_template as dt
from first_hires_exoplanet import figs_phase2

import pyodine
from pyodine.fitters.lmfit_wrapper import LmfitWrapper
from pyodine.models.base import ParameterSet
from pyodine.template.base import StellarTemplate_Chunked
from pyodine.template.normalize import SimpleNormalizer


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_HERE, 'data')
FIG_DIR = os.path.join(os.path.dirname(_HERE), 'docs', 'figs')

EPOCH = 'HI.19980826.19326'
ORDER = 63
CKMS = 299792.458


def epoch_file(koaid=EPOCH):
    sci = os.path.join(lp.DEFAULT_REDUX, 'reduce_' + koaid[3:11], 'Science')
    hit = [p for p in os.listdir(sci) if p.startswith('spec1d_' + koaid) and p.endswith('.fits')]
    if len(hit) != 1:
        raise FileNotFoundError('spec1d for {:s}'.format(koaid))
    return os.path.join(sci, hit[0])


# ---------------------------------------------------------------------------
# Preparation: model_single_observation up to the fitting loop
# ---------------------------------------------------------------------------

def prepare(Pars, obs_file, temp_file):
    """ Everything `model_single_observation` builds before it fits.

    Returns:
        dict
    """
    obs = lp.ObservationWrapper(obs_file)
    template = StellarTemplate_Chunked(temp_file)
    iod = lp.IodineTemplate(Pars.i2_to_use)

    orders = template.orders_unique
    template_ind = template.get_order_indices(template[0].order)
    obs_order_min, coverage = obs.check_wavelength_range(
        template[template_ind[0]].w0, template[template_ind[-1]].w0)
    order_correction = obs_order_min - template[template_ind[0]].order

    weight = obs.compute_weight(weight_type=Pars.weight_type, rel_noise=Pars.rel_noise)

    normalizer = SimpleNormalizer(reference=Pars.ref_spectrum)
    ref_velocity = normalizer.guess_velocity(
        obs[Pars.velgues_order_range[0]:Pars.velgues_order_range[1]],
        delta_v=Pars.delta_v, maxlag=Pars.maxlag)
    obs_velocity = ref_velocity - template.velocity_offset

    chunks = pyodine.chunks.auto_wave_comoving(
        obs, template, orders=orders, padding=Pars.chunk_padding,
        order_correction=order_correction, delta_v=Pars.chunk_delta_v)

    return dict(obs=obs, template=template, iod=iod, orders=orders,
                order_correction=order_correction, coverage=coverage,
                weight=weight, ref_velocity=ref_velocity,
                obs_velocity=obs_velocity, chunks=chunks)


def run0_params(Pars, prep, model, fitter):
    """ Starting parameters for every chunk, exactly as run 0 builds them.

    No chunk is fitted here: this is the driver's guess, its polynomial
    smoothing of the wavelength guesses over each order, and
    ``constrain_parameters``.
    """
    run = Pars.model_runs[0]
    chunks = prep['chunks']
    starting = [model.guess_params(ch) for ch in chunks]
    for name, key in (('wave_slope', 'pre_wave_slope_deg'),
                      ('wave_intercept', 'pre_wave_intercept_deg')):
        if run.get(key, 0) > 0:
            poly = pyodine.lib.misc.smooth_parameters_over_orders(
                starting, name, chunks, deg=run[key])
            for i in range(len(chunks)):
                starting[i][name] = poly[i]
    lmfit_params = [fitter.convert_params(p, to_lmfit=True) for p in starting]
    run_results = {0: {'velocity_guess': prep['obs_velocity']}}
    return Pars.constrain_parameters(lmfit_params, 0, run_results, fitter)


def build_run0(Pars, prep):
    """ Run 0's model and fitter. """
    run = Pars.model_runs[0]
    model = pyodine.models.spectrum.SimpleModel(
        run['lsf_model'], run['wave_model'], run['cont_model'], prep['iod'],
        stellar_template=prep['template'], osample_factor=Pars.osample_obs,
        conv_width=Pars.lsf_conv_width)
    return model, LmfitWrapper(model)


def pick_chunk(prep, ech_order=ORDER):
    """ The chunk nearest the middle of an order with no masked pixel.

    The literal middle chunk of order 63 in this epoch has 14 of its 40
    pixels at zero weight (a PypeIt-masked run near detector column 1050),
    which would test the machinery on a third less data than it will usually
    have.  The rule is stated and applied blind to the fit.

    Returns:
        tuple: (chosen index, literal middle index, zero-weight pixels in
        the literal middle chunk)
    """
    pos = prep['obs'].index_of(ech_order)
    idx = [i for i, ch in enumerate(prep['chunks']) if ch.order == pos]
    mid = idx[len(idx) // 2]
    w = prep['weight']

    def nbad(i):
        ch = prep['chunks'][i]
        return int(np.sum(w[ch.order, ch.abspix[0]:ch.abspix[-1] + 1] <= 0))
    clean = [i for i in idx if nbad(i) == 0]
    if not clean:
        raise RuntimeError('No unmasked chunk in order {:d}'.format(ech_order))
    best = min(clean, key=lambda i: abs(i - mid))
    return best, mid, nbad(mid)


# ---------------------------------------------------------------------------
# The fit
# ---------------------------------------------------------------------------

def fit_chunk(fitter, chunk, lmfit_params, weight, i, use_chauvenet=True):
    """ `pipe_lib.model_all_chunks`, for one chunk -- faithfully.

    The driver re-fits after the Chauvenet criterion only ``if
    any(mask)==False``, i.e. only when *every* pixel fails it; the evident
    intent is "if any pixel fails".  As written, ``use_chauvenet_pixels``
    almost never does anything.  This mirrors the driver, bug included, and
    reports how many pixels the criterion would have rejected.

    Returns:
        tuple: (result, pixel weights used, pixels re-fitted without,
        pixels the criterion flags)
    """
    ch_w = np.array(weight[chunk.order, chunk.abspix[0]:chunk.abspix[-1] + 1])
    ch_w[ch_w <= 0.] = 0.
    result = fitter.fit(chunk, lmfit_params, weight=ch_w, chunk_ind=i)
    n_out, n_flag = 0, 0
    if result.lmfit_result is not None and use_chauvenet:
        mask, mask_true, mask_false = pyodine.lib.misc.chauvenet_criterion(result.residuals)
        n_flag = int(np.sum(~mask))
        if not any(mask):
            n_out = len(mask_false[0])
            ch_w[mask_false] = 0.
            result = fitter.fit(chunk, lmfit_params, weight=ch_w, chunk_ind=i)
    return result, ch_w, n_out, n_flag


def parameter_table(result, start, noise_scale):
    """ Every parameter: value, standard error, start, bounds, varied. """
    lr = result.lmfit_result
    rows = []
    for name, p in lr.params.items():
        rows.append(dict(param=name, value=float(p.value),
                         stderr=float(p.stderr) if p.stderr is not None else np.nan,
                         start=float(start[name].value),
                         min=float(p.min), max=float(p.max), vary=bool(p.vary)))
    tab = Table(rows)
    tab.meta.update(dict(redchi=float(lr.redchi),
                         redchi_measured=float(lr.redchi) / noise_scale ** 2,
                         ndata=int(lr.ndata), nvarys=int(lr.nvarys),
                         nfev=int(lr.nfev), success=bool(lr.success),
                         message=str(lr.message)))
    return tab


def monte_carlo(fitter, chunk, lmfit_params, ch_w, i, result, noise_scale, nmc,
                start='best', seed=4):
    """ Refit the same chunk to its best-fit model plus measured noise.

    The noise is drawn at ``noise_scale`` times PypeIt's propagated sigma --
    the level actually measured between exposures -- while the weights stay
    PypeIt's, exactly as in the real fit.  Two questions, two starts:

    * ``start='best'``: from the best fit.  Tests the standard errors: if
      `lmfit`'s rescaled covariance is honest, the scatter of each parameter
      matches its standard error, and the reduced chi-square is
      ``noise_scale**2``.
    * ``start='guess'``: from the driver's own starting guess.  Tests whether
      the minimum is found at all from where `pyodine` starts.  A reduced
      chi-square above ``noise_scale**2`` means fits stopping short of it.
    """
    if start == 'best':
        lmfit_params = copy.deepcopy(result.lmfit_result.params)
    rng = np.random.default_rng(seed)
    model_flux = result.fitted_spectrum.flux
    sig = np.where(ch_w > 0, 1. / np.sqrt(np.where(ch_w > 0, ch_w, 1.)), 0.) * noise_scale
    names = [n for n, p in result.lmfit_result.params.items() if p.vary]
    draws, redchi = {n: [] for n in names}, []
    for _ in range(nmc):
        fake = copy.copy(chunk)
        fake.flux = model_flux + rng.normal(size=len(model_flux)) * sig
        r = fitter.fit(fake, lmfit_params, weight=ch_w, chunk_ind=i)
        if r.lmfit_result is None:
            continue
        redchi.append(r.lmfit_result.redchi)
        for n in names:
            draws[n].append(r.lmfit_result.params[n].value)
    rows = [dict(param=n, mc_std=float(np.std(draws[n], ddof=1)),
                 mc_mean=float(np.mean(draws[n])),
                 stderr=float(result.lmfit_result.params[n].stderr),
                 ratio=float(np.std(draws[n], ddof=1) /
                             result.lmfit_result.params[n].stderr))
            for n in names]
    tab = Table(rows)
    tab.meta['n'] = len(redchi)
    tab.meta['redchi_median'] = float(np.median(redchi))
    tab.meta['frac_converged'] = float(np.mean(np.asarray(redchi) < 2 * noise_scale ** 2))
    tab['start'] = start
    return tab


# ---------------------------------------------------------------------------
# The independent checks
# ---------------------------------------------------------------------------

def template_epoch_jd(template_fits=dt.TEMPLATE_FITS):
    """ Mean mid-exposure JD(UTC) of the three co-added template exposures. """
    from astropy.io import fits
    hdr = fits.getheader(template_fits)
    mids = []
    for k in range(int(hdr['NCOADD'])):
        sci = os.path.join(lp.DEFAULT_REDUX, 'reduce_19980826', 'Science')
        f = [p for p in os.listdir(sci) if p.startswith('spec1d_' + hdr['FRAME{:d}'.format(k)])
             and p.endswith('.fits')][0]
        o = lp.ObservationWrapper(os.path.join(sci, f))
        mids.append(o.bary_date)
    return float(np.mean(mids))


def expected_velocity(prep):
    """ The velocity the fit should find, and its two parts.

    `pyodine` shifts the template by ``velocity`` onto the observation, in the
    observed frame, so it should find the change of the star's observed
    velocity between the template and the epoch: the planet's motion minus
    the change of the barycentric correction.  The planet's motion comes from
    a circular orbit fitted to the modern catalogue at the modern period
    (`figs_phase2.fit_circular`).  Catalogue times are JD(UTC) whatever the
    column says (phase 2 trap 5), as are ours.
    """
    cat = figs_phase2.load_catalogue()
    orbit = figs_phase2.fit_circular(np.asarray(cat['bjd']), np.asarray(cat['rv']),
                                     e=np.asarray(cat['e_rv']))
    t_obs = prep['obs'].bary_date
    t_tmp = template_epoch_jd()
    d_planet = float(orbit['model'](t_obs) - orbit['model'](t_tmp))
    d_bary = float(prep['template'].bary_vel_corr - prep['obs'].bary_vel_corr)
    return dict(expected=d_planet + d_bary, planet=d_planet, bary=d_bary,
                K=orbit['K'], t_obs=t_obs, t_tmp=t_tmp)


def wavelength_check(result, chunk):
    """ The fitted wavelength solution against PypeIt's, over the chunk. """
    p = result.params
    w_fit = p['wave_intercept'] + p['wave_slope'] * chunk.pix
    w_pyp = chunk.wave
    dv = (w_fit / w_pyp - 1.) * CKMS * 1e3
    slope_pyp = np.polyfit(chunk.pix, w_pyp, 1)[0]
    return dict(dv_mean=float(np.mean(dv)), dv_span=float(np.ptp(dv)),
                slope_ratio=float(p['wave_slope'] / slope_pyp))


def arc_fwhm(ech_order=ORDER):
    t = Table.read(os.path.join(DATA_DIR, 'deconv_arc_widths.csv'))
    return float(np.median(t['fwhm_px'][t['order'] == ech_order]))


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(fits_, outfile):
    """ Data, model and residuals of the one chunk, for both templates. """
    import matplotlib.pyplot as plt
    INK, C1, C2 = '#0b0b0b', '#2a78d6', '#eb6834'
    fig, axs = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                            gridspec_kw=dict(height_ratios=[2, 1]))
    first = True
    for os_, col in ((10, C2), (1, C1)):
        f = fits_[os_]
        r, ch, w, k = f['result'], f['chunk'], f['weight'], f['noise_scale']
        wave = r.fitted_spectrum.wave
        if first:
            axs[0].plot(wave, ch.flux, 'o', ms=3.5, color=INK, label='observed', zorder=3)
            first = False
        axs[0].plot(wave, r.fitted_spectrum.flux, '-', color=col, lw=1.2,
                    label='model, template oversampling {:d}'.format(os_))
        sig = np.where(w > 0, 1. / np.sqrt(np.where(w > 0, w, 1.)), np.nan) * k
        axs[1].plot(wave, r.residuals / sig, 'o-', ms=3, lw=0.8, color=col,
                    label='oversampling {:d}: chi2/dof = {:.2f} (measured noise)'.format(
                        os_, f['table'].meta['redchi_measured']))
    axs[0].set_ylabel('counts')
    axs[0].legend(fontsize=8)
    axs[0].set_title('{:s}, order {:d}, chunk {:d}: run 0 (single Gaussian)'.format(
        EPOCH, ORDER, fits_[10]['index']), loc='left', fontsize=10)
    axs[1].axhline(0, color='0.5', lw=0.8)
    axs[1].set_ylabel('residual / measured sigma')
    axs[1].set_xlabel('vacuum wavelength, observed frame, from the fitted solution (A)')
    axs[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outfile, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def parser(options=None):
    p = argparse.ArgumentParser(
        description='Fit one chunk of one order of one epoch (phase 3, prompt 4)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--nmc', type=int, default=30,
                   help='Monte Carlo refits of the same chunk')
    p.add_argument('--datadir', type=str, default=DATA_DIR)
    p.add_argument('--figdir', type=str, default=FIG_DIR)
    return p.parse_args() if options is None else p.parse_args(options)


def main(pargs):
    Pars = hires_pars.Parameters()
    k = Pars.noise_scale
    obs_file = epoch_file()
    fits_, rows, mcrows = {}, [], []
    for os_ in (10, 1):
        print('=== template oversampling {:d}: {:s}'.format(os_, dt.TEMPLATE_H5[os_]))
        t0 = time.time()
        prep = prepare(Pars, obs_file, dt.TEMPLATE_H5[os_])
        model, fitter = build_run0(Pars, prep)
        lmfit_params = run0_params(Pars, prep, model, fitter)
        i, mid, mid_bad = pick_chunk(prep)
        ch = prep['chunks'][i]
        start = copy.deepcopy(lmfit_params[i])
        t_prep = time.time() - t0

        t0 = time.time()
        result, ch_w, n_out, n_flag = fit_chunk(fitter, ch, lmfit_params[i], prep["weight"], i)
        t_fit = time.time() - t0
        if result.lmfit_result is None:
            raise RuntimeError('The fit failed')
        tab = parameter_table(result, start, k)

        print('  {:d} chunks built; fitting chunk {:d} only (order {:d}, pixels {:d}-{:d}, '
              '{:.2f}-{:.2f} A)'.format(len(prep['chunks']), i, ORDER, ch.abspix[0],
                                        ch.abspix[-1], ch.wave[0], ch.wave[-1]))
        print('  (the literal middle chunk, {:d}, has {:d} of 40 pixels masked)'.format(
            mid, mid_bad))
        print('  order correction {:d}; velocity guess vs reference {:+.3f} km/s, '
              'vs template {:+.3f} km/s'.format(prep['order_correction'],
                                                prep['ref_velocity'] / 1e3,
                                                prep['obs_velocity'] / 1e3))
        print('  prep {:.1f} s, fit {:.2f} s, nfev {:d}, success {}, Chauvenet outliers {:d}, '
              '(flagged {:d}), zero-weight pixels {:d}'.format(t_prep, t_fit, tab.meta['nfev'],
                                              tab.meta['success'], n_out, n_flag, int(np.sum(ch_w == 0))))
        print('  reduced chi2: {:.3f} as fitted (PypeIt ivar), {:.2f} against the measured '
              'noise; ndata {:d}, nvarys {:d}'.format(tab.meta['redchi'],
                                                      tab.meta['redchi_measured'],
                                                      tab.meta['ndata'], tab.meta['nvarys']))
        for r in tab:
            print('    {:15s} {:>16.8g} +- {:<12.4g} start {:>14.8g}  [{:g}, {:g}]{}'.format(
                r['param'], r['value'], r['stderr'], r['start'], r['min'], r['max'],
                '' if r['vary'] else '  FIXED'))

        exp = expected_velocity(prep)
        wv = wavelength_check(result, ch)
        p, e = result.params, result.errors
        print('  checks:')
        print('    velocity  {:+8.1f} +- {:.1f} m/s   expected {:+8.1f} (barycentric {:+.1f}, '
              'planet {:+.1f}, K = {:.1f})'.format(p['velocity'], e['velocity'],
                                                   exp['expected'], exp['bary'],
                                                   exp['planet'], exp['K']))
        print('    lsf_fwhm  {:.3f} +- {:.3f} px   ThAr, order {:d}: {:.3f} px'.format(
            p['lsf_fwhm'], e['lsf_fwhm'], ORDER, arc_fwhm()))
        print('    iod_depth {:.4f} +- {:.4f}   tem_depth {:.4f} +- {:.4f}   (both 1 if right)'.format(
            p['iod_depth'], e['iod_depth'], p['tem_depth'], e['tem_depth']))
        print('    wavelength vs PypeIt: mean {:+.1f} m/s, span over the chunk {:.1f} m/s, '
              'dispersion ratio {:.6f}'.format(wv['dv_mean'], wv['dv_span'], wv['slope_ratio']))
        print('    continuum {:.0f} at chunk centre; chunk median flux {:.0f}'.format(
            p['cont_intercept'], np.median(ch.flux)))

        for start in ('best', 'guess'):
            mc = monte_carlo(fitter, ch, lmfit_params[i], ch_w, i, result, k, pargs.nmc,
                             start=start)
            print('  Monte Carlo from the {:s}, {:d} refits of this chunk at the measured '
                  'noise: median redchi {:.3f} (expect {:.3f}); {:.0f}% reach it'.format(
                      'best fit' if start == 'best' else "driver's starting guess",
                      mc.meta['n'], mc.meta['redchi_median'], k ** 2,
                      100 * mc.meta['frac_converged']))
            for r in mc:
                print('    {:15s} MC scatter {:12.4g}   stderr {:12.4g}   ratio {:.2f}'.format(
                    r['param'], r['mc_std'], r['stderr'], r['ratio']))
            mc['osample'] = os_
            mc.meta = {}
            mcrows.append(mc)

        tab['osample'] = os_
        tab['chunk'] = i
        rows.append(tab)
        fits_[os_] = dict(result=result, chunk=ch, weight=ch_w, noise_scale=k,
                          table=tab, index=i, expected=exp, wave=wv)

    from astropy.table import vstack
    out = vstack(rows, metadata_conflicts='silent')
    for os_ in (10, 1):
        f = fits_[os_]
        out.meta['redchi_os{:d}'.format(os_)] = f['table'].meta['redchi']
        out.meta['redchi_measured_os{:d}'.format(os_)] = f['table'].meta['redchi_measured']
        out.meta['v_expected_os{:d}'.format(os_)] = f['expected']['expected']
    out.meta['epoch'] = EPOCH
    out.meta['order'] = ORDER
    out.meta['noise_scale'] = k
    out.write(os.path.join(pargs.datadir, 'one_chunk_fit.csv'), overwrite=True)
    vstack(mcrows, metadata_conflicts='silent').write(
        os.path.join(pargs.datadir, 'one_chunk_mc.csv'), overwrite=True)
    fig = os.path.join(pargs.figdir, 'fig_p4_one_chunk.png')
    make_figure(fits_, fig)
    print('Figure: {:s}'.format(fig))


if __name__ == '__main__':
    main(parser())
