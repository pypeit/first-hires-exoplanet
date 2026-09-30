""" Parameters for the HIRES iodine forward model (phase 3, prompt 4).

Modelled on `utilities_lick/pyodine_parameters.py` of `pyodine`
4488b0914fe5b272b787982647691045bff2604a, and read by the same drivers
(`pyodine_model_observations.py`, `pyodine_create_templates.py`) through the
same two classes.  Where a value is Lick's, it is Lick's because nothing
measured here says otherwise; where it differs, the comment says what said so.

What differs from Lick, and why:

* ``weight_type = 'ivar'``.  PypeIt propagates an inverse variance; prompt 2
  taught ``compute_weight`` to use it.  Lick uses flat weights.
* ``noise_scale = 0.426``.  PypeIt's inverse variance overstates the noise:
  exposures differenced pixel by pixel scatter at 0.426 of the propagated
  sigma (prompt 3, ``deconvolve_template.noise_check``).  Decided after
  prompt 3 (Q&A, question 4): the weights are left as PypeIt made them --
  relative weighting is unaffected by a uniform factor -- and every chi-square
  is divided by ``noise_scale**2`` when it is reported.  The parameter
  uncertainties need no correction: `lmfit.minimize` runs with its default
  ``scale_covar=True``, which already rescales them by sqrt(reduced chi2).
  `pyodine` itself never reads this attribute.
* ``wavelength_scale = 'vacuum'`` and ``i2_to_use = 1``.  PypeIt reports
  vacuum; see ``conf.py``.
* ``velgues_order_range = (15, 37)``: the positions, in ascending echelle
  order, of orders 72-93, which are blueward of 5000 A and free of iodine.
  Lick's (2, 15) would land inside the iodine band here.
* ``maxlag = 500`` for observations (Lick: 10000).  The guess against
  `pyodine`'s air-wavelength reference comes back near +75 km/s (prompt 3),
  and 500 steps of 1 km/s covers it with room to spare.
* ``number_cores = 4`` (Lick: 12), and ``MPLBACKEND=Agg`` must be set by any
  driver before import: plotting then forking hangs on macOS (prompt 1).
* Template creation copies prompt 3's deconvolution exactly, including its
  chunking (every chunk that fits, 700 over orders 58-71).  Prompt 3 also
  wrote a template at ``osample_temp = 1``, which alone reproduces the
  observation; the Q&A after prompt 3 carries both forward, and the value
  here stays Lick's until prompt 5 decides.

The LSF model, settled in prompt 6 (`compare_lsf_models.py`).  Run 0 is a
single Gaussian started at the ThAr width (2.2 px); run 1 is `pyodine`'s
**super-Gaussian**.  On the prompt-5 epoch the four workable models -- single,
super, Lick's and SONG's multi-Gaussians -- are indistinguishable in
per-chunk velocity scatter (paired bootstrap, all within +/-13 m/s and every
68% interval including zero); the super-Gaussian ties for the lowest scatter
and has the lowest chi-square, is centred (Lick's multi-Gaussian is not, and
carries a +0.03 px centroid), and has 4 parameters to their 10.  The
Hermite model is worse and biased; a run 2 with the smoothed LSF held fixed
does not help.  The single Gaussian is the equivalent fallback.

Run 1's shape parameters are started with :func:`start_inside`: `fit_lsfs`
puts every parameter a Gaussian does not need at ~1e-12, where `lmfit`'s
relative finite-difference step cannot move it and the whole fit returns its
start -- which is what the super-Gaussian and Hermite trials of prompt 6 did
until this was fixed.
"""

from pyodine import models

import logging
import os
import sys


utilities_dir_path = os.path.dirname(os.path.realpath(__file__))

###############################################################################
## Instrument-specific setup of the LSF models.  Lick's, pending prompt 6.   #
###############################################################################

# For the MultiGaussian model: 11 positions & sigmas for the central Gauss and
# the satellites, in pixels.  Lick's values; Butler et al. (1996) place the
# satellites 0.5 pixel apart.
_multigauss_setup_dict = {
        'positions': [-2.4, -2.1, -1.6, -1.1, -0.6, 0.0, 0.6, 1.1, 1.6, 2.1, 2.4],
        'sigmas':    [ 0.3,  0.3,  0.3,  0.3,  0.3, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3]
        }

#: Starting FWHM of the single-Gaussian LSF, pixels: the ThAr median over the
#: iodine orders of 1998-08-26 (prompt 3).  Also Lick's value, by coincidence.
LSF_FWHM_START = 2.2

#: Measured / propagated noise (prompt 3).  Divide reported chi-square by the
#: square of this; see the module docstring.
NOISE_SCALE = 0.426

#: Physical bounds of the super-Gaussian's parameters, pixels (prompt 6).
#: Lick's generic p +/- (|p| + 0.3) would allow a negative exponent.
SUPER_BOUNDS = {'sigma': (0.2, 3.0), 'exponent': (1.0, 4.0),
                'left': (0.0, 1.0), 'right': (0.0, 1.0)}


def start_inside(value, lo, hi, floor=1e-3, margin=0.02):
    """ A starting value `lmfit`'s leastsq can move (prompt 6).

    MINPACK's finite-difference step is relative to the value, so a parameter
    started at ~1e-12 -- where `fit_lsfs` puts every shape parameter a
    Gaussian does not need -- has a zero derivative, and the fit returns its
    start.  No smaller than ``floor`` in magnitude, and ``margin`` of the
    range inside either bound.
    """
    v = value if abs(value) >= floor else (floor if value >= 0 else -floor)
    span = hi - lo
    return float(min(max(v, lo + margin * span), hi - margin * span))


def lsf_bounds(model_name, name, value):
    """ Bounds for one run-1 LSF parameter: physical for the super-Gaussian,
    Lick's rule otherwise. """
    if model_name == 'SuperGaussian':
        return SUPER_BOUNDS[name]
    if model_name == 'SingleGaussian':
        return 0.5, 4.0
    return value - abs(value) - 0.3, value + abs(value) + 0.3


def _run_dict(lsf_model, lsf_setup_dict=None, pre_wave_deg=0):
    """ One entry of ``model_runs``, with Lick's bookkeeping settings. """
    run = {
        'lsf_model': lsf_model,
        'wave_model': models.wave.LinearWaveModel,
        'cont_model': models.cont.LinearContinuumModel,
        'pre_wave_slope_deg': pre_wave_deg,
        'pre_wave_intercept_deg': pre_wave_deg,
        'use_chauvenet_pixels': True,
        'save_result': True,
        'save_filetype': 'h5py',
        'wave_slope_deg': 3,
        'wave_intercept_deg': 3,
        'plot_success': True,
        'plot_analysis': True,
        'plot_chunks': [150, 250, 400],
        'plot_lsf_pars': True,
        'save_median_pars': True,
    }
    if lsf_setup_dict is not None:
        run['lsf_setup_dict'] = lsf_setup_dict
    return run


class Parameters:
    """The control commands for modelling an observation.

    See the module docstring for every departure from Lick.
    """

    def __init__(self):

        # Setup the logging if not existent yet
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                                format='%(message)s')

        # General parameters
        self.osample_obs = 4                    # Oversample factor for the observation modeling
        self.lsf_conv_width = 6.                # LSF is evaluated over this many pixels (times 2)
        self.number_cores = 4                   # Number of processor cores for multiprocessing

        self.log_config_file = os.path.join(utilities_dir_path, 'logging.json')
        self.log_level = logging.INFO

        self.use_progressbar = False

        # Tellurics: none masked.  The iodine band carries weak H2O lines
        # near 5900-6000 A; not yet assessed.
        self.telluric_mask = None
        self.tell_wave_range = (None, 6500)
        self.tell_dispersion = 0.002

        # Chunking: the observation's chunks follow the template's, shifted
        # by the relative barycentric velocity
        self.chunking_algorithm = 'auto_wave_comoving'
        self.order_range = (None, None)         # (None, None): the template's orders
        self.chunk_width = 40
        self.chunk_padding = 6
        self.chunks_per_order = None
        self.chunk_delta_v = None

        # Reference spectrum and the first velocity guess
        self.ref_spectrum = 'arcturus'
        self.velgues_order_range = (15, 37)     # orders 72-93: blue of 5000 A, no I2
        self.delta_v = 1000.
        self.maxlag = 500

        self.normalize_chunks = False

        # Weighting of pixels
        self.bad_pixel_mask = False
        self.bad_pixel_cutoff = 0.22
        self.correct_obs = False
        self.weight_type = 'ivar'               # PypeIt's propagated inverse variance
        self.rel_noise = 0.008                  # only for weight_type='inverse'
        self.noise_scale = NOISE_SCALE          # reporting only; see module docstring

        # I2 atlas
        self.i2_to_use = 1                      # conf.my_iodine_atlases
        self.wavelength_scale = 'vacuum'

        self.vel_analysis_plots = -1

        # Run 0: single Gaussian, the first wavelength solution.
        # Run 1: super-Gaussian (prompt 6), started from run 0's median LSF.
        self.model_runs = {
            0: _run_dict(models.lsf.SingleGaussian, pre_wave_deg=3),
            1: _run_dict(models.lsf.SuperGaussian),
        }

    def constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        """Constrain the lmfit parameters for each run (Lick's scheme).

        :param lmfit_params: A list of :class:`lmfit.Parameters` objects
            for each chunk.
        :type lmfit_params: list[:class:`lmfit.Parameters`]
        :param run_id: The current run_id.
        :type run_id: int
        :param run_results: Results from previous modelling runs.
        :type run_results: dict
        :param fitter: The fitter used in the modelling.
        :type fitter: :class:`LmfitWrapper`

        :return: The updated list of :class:`lmfit.Parameters` objects.
        :rtype: list[:class:`lmfit.Parameters`]
        """
        logging.info('')
        logging.info('Constraining parameters for RUN {}'.format(run_id))

        if run_id == 0:
            for i in range(len(lmfit_params)):
                lmfit_params[i]['velocity'].set(
                        value=run_results[run_id]['velocity_guess'])
                lmfit_params[i]['lsf_fwhm'].set(
                        value=LSF_FWHM_START, min=0.5, max=4.0)
                # iod_depth is a Beer-Lambert exponent applied on top of the
                # adapter's alpha = 2.59 (prompt 2, change 1), so 1 is the
                # expectation; keep it positive
                lmfit_params[i]['iod_depth'].set(min=0.1)

        elif run_id == 1:
            median_lsf_pars = run_results[0]['median_pars'].filter('lsf')
            lsf_fit_pars = fitter.fit_lsfs(self.model_runs[0]['lsf_model'], median_lsf_pars)
            logging.info('')
            logging.info('Fitted LSF parameters:')
            logging.info(lsf_fit_pars)

            for i in range(len(lmfit_params)):
                lmfit_params[i]['velocity'].set(
                        value=run_results[0]['median_pars']['velocity'])
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[0]['median_pars']['iod_depth'])
                lmfit_params[i]['tem_depth'].set(
                        value=run_results[0]['median_pars']['tem_depth'])
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[0]['wave_slope_fit'][i])
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[0]['wave_intercept_fit'][i])
                lmfit_params[i]['cont_intercept'].set(
                        value=run_results[0]['results'][i].params['cont_intercept'])
                lmfit_params[i]['cont_slope'].set(
                        value=run_results[0]['results'][i].params['cont_slope'])
                for p in lsf_fit_pars.keys():
                    lo, hi = lsf_bounds(self.model_runs[1]['lsf_model'].name(), p,
                                        lsf_fit_pars[p])
                    lmfit_params[i]['lsf_' + p].set(
                        value=start_inside(lsf_fit_pars[p], lo, hi), min=lo, max=hi)

        return lmfit_params


class Template_Parameters:
    """The control commands for template creation.

    The deconvolution settings are prompt 3's, which are Lick's; the O-star
    (B-star) modelling runs are Lick's scheme, pending prompt 6.
    """

    def __init__(self):

        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                                format='%(message)s')

        self.osample_obs = 4
        self.lsf_conv_width = 6.

        self.log_config_file = os.path.join(utilities_dir_path, 'logging.json')
        self.log_level = logging.INFO

        self.use_progressbar = True

        self.telluric_mask = None
        self.tell_wave_range = (None, 6500)
        self.tell_dispersion = 0.002

        # Chunking, as prompt 3: every 40-pixel chunk that fits in orders
        # 58-71 (positions 1-14), 12 pixels of padding -- 700 chunks
        self.chunking_algorithm = 'auto_equal_width'
        self.temp_order_range = (1, 14)
        self.chunk_width = 40
        self.chunk_padding = 12
        self.chunks_per_order = None
        self.pix_offset0 = None
        self.wavelength_dict = self.chunk_wavelengths

        self.ref_spectrum = 'arcturus'
        self.velgues_order_range = (15, 37)
        self.delta_v = 1000.
        self.maxlag = 500

        self.normalize_chunks = False

        self.bad_pixel_mask = True
        self.bad_pixel_cutoff = 0.22
        self.correct_obs = True
        self.weight_type = 'ivar'
        self.rel_noise = 0.008
        self.noise_scale = NOISE_SCALE

        self.i2_to_use = 1
        self.wavelength_scale = 'vacuum'

        self.jansson_run_model = 1
        self.chunk_weights_redchi = False
        # Prompt 3.  At osample_temp = 10 the template does not reproduce the
        # observation (chi2 5.2 against the measured noise) because the
        # deconvolver splines barely-sampled pixels; at 1 it does (chi2 1.1).
        # Both templates exist.  Prompts 5 and 6 decided for 1: lower
        # per-chunk scatter for every workable LSF model (prompt 6, paired).
        self.deconvolution_pars = {
                'osample_temp': 1,
                'jansson_niter': 1200,
                'jansson_zerolevel': 0.00,
                'jansson_contlevel': 1.02,
                'jansson_conver': 0.2,
                'jansson_chi_change': 1e-6,
                'lsf_conv_width': self.lsf_conv_width,
                }

        self.jansson_lsf_smoothing = {
                'do_smoothing': False,
                'smooth_lsf_run': 1,
                'smooth_pixels': 160,
                'smooth_orders': 3,
                'order_separation': 15,
                'smooth_manual_redchi': False,
                }

        self.model_runs = {
            0: _run_dict(models.lsf.SingleGaussian, pre_wave_deg=3),
            1: _run_dict(models.lsf.SuperGaussian),
        }

    def constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        """Constrain the lmfit parameters when modelling the hot star (Lick's).

        Velocity and template depth are fixed: a hot star has no lines to
        move and no template.
        """
        logging.info('')
        logging.info('Constraining parameters for RUN {}'.format(run_id))

        if run_id == 0:
            for i in range(len(lmfit_params)):
                lmfit_params[i]['velocity'].set(vary=False)
                lmfit_params[i]['tem_depth'].set(vary=False)
                lmfit_params[i]['lsf_fwhm'].set(
                        value=LSF_FWHM_START, min=0.5, max=4.0)
                lmfit_params[i]['iod_depth'].set(min=0.1)

        elif run_id == 1:
            median_lsf_pars = run_results[0]['median_pars'].filter('lsf')
            lsf_fit_pars = fitter.fit_lsfs(self.model_runs[0]['lsf_model'], median_lsf_pars)
            logging.info('')
            logging.info('Fitted LSF parameters:')
            logging.info(lsf_fit_pars)

            for i in range(len(lmfit_params)):
                lmfit_params[i]['velocity'].set(vary=False)
                lmfit_params[i]['tem_depth'].set(vary=False)
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[0]['median_pars']['iod_depth'])
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[0]['wave_slope_fit'][i])
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[0]['wave_intercept_fit'][i])
                lmfit_params[i]['cont_intercept'].set(
                        value=run_results[0]['results'][i].params['cont_intercept'])
                lmfit_params[i]['cont_slope'].set(
                        value=run_results[0]['results'][i].params['cont_slope'])
                for p in lsf_fit_pars.keys():
                    lo, hi = lsf_bounds(self.model_runs[1]['lsf_model'].name(), p,
                                        lsf_fit_pars[p])
                    lmfit_params[i]['lsf_' + p].set(
                        value=start_inside(lsf_fit_pars[p], lo, hi), min=lo, max=hi)

        return lmfit_params

    # Start and end wavelengths of chunks if the wavelength_defined chunking
    # algorithm is used in template creation
    chunk_wavelengths = {
            'start_wave': [],
            'end_wave': []
            }
