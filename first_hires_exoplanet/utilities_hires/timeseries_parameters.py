""" Parameters for combining chunk velocities into a timeseries (phase 3,
prompt 4).

Modelled on `utilities_lick/timeseries_parameters.py` of `pyodine`
4488b0914fe5b272b787982647691045bff2604a.  Nothing here is exercised until
prompt 7; the values are Lick's except where a known property of these data
says otherwise.

What differs from Lick, and why:

* ``use_hip_for_bvc = False``.  The barycentric correction is computed from
  coordinates, not by parsing a Hipparcos number out of the star name: the
  1997-98 headers name HD 187123 inconsistently ('187123', 'H187123', and a
  program code in OBJECT), and a name that parses wrongly would fail silently.
* The barycentric correction is still recomputed here with `barycorrpy`
  (``compute_bvc = 'precise'``), which is what the adapter already supplies
  per epoch.  PypeIt's own heliocentric correction is never used: it carries a
  sign error worth +13 m/s (phase 2, trap 6).
* ``txt_flux_chunk`` indexes chunks 251-253 of the 700 that the prompt-3
  templates carry; they exist, which is all Lick's values guarantee there.
* ``rv_err`` inherits `lmfit`'s ``scale_covar`` rescaling, so the 2.35x
  overstatement of PypeIt's inverse variance (prompt 3) does not propagate
  into it.  Any chi-square this stage reports must be divided by
  ``pyodine_parameters.NOISE_SCALE**2``.
"""

import logging
import os

utilities_dir_path = os.path.dirname(os.path.realpath(__file__))


class Timeseries_Parameters:

    def __init__(self):

        # Logging options
        self.log_config_file = os.path.join(utilities_dir_path, 'logging.json')
        self.log_level = logging.INFO

        # Rejection lists name the modelling results, not the observations
        self.reject_type = 'res_files'

        # Barycentric correction: barycorrpy, from coordinates (see docstring)
        self.compute_bvc = 'precise'            # 'precise', 'predictive', or 'no'
        self.use_hip_for_bvc = False
        self.solar_bvc = False

        # Weighting algorithm ('song' or 'lick'); Lick's choice and settings
        self.weighting_algorithm = 'song'
        self.weighting_pars_song = {
                'good_chunks': None,
                'good_orders': None,
                'sig_limit_low': 4.,
                'sig_limit_up': 1000.,
                'sig_correct': 1000.,
                'reweight_alpha': 1.8,
                'reweight_beta': 8.0,
                'reweight_sigma': 2.0,
                'weight_correct': 0.01,
                }
        self.weighting_pars_lick = {
                'percentile': 0.997,
                'maxchi': 100000000.,
                'min_counts': 1000.,
                'default_sigma': 1000.,
                'useage_percentile': 0.997,
                }

        # Chromatic index
        self.do_crx = True
        self.crx_pars = {
                'crx_sigma': 0.,
                'crx_iterative': False,
                'crx_max_iters': 10
                }

        # Text output
        self.txt_outkeys = ['bary_date', 'rv_bc', 'rv_err']
        self.txt_delimiter = '\t'
        self.txt_header = ''
        self.txt_outformat = ['%10.5f', '%6.4f', '%3.4f']
        self.txt_detailed = False
        self.txt_flux_chunk = [251, 252, 253]

        self.save_comb_res = True
        self.plot_analysis = True
        self.print_outliers = True
